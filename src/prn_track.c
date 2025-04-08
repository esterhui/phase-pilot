/* tone tracker
 *
 *
 * $Id: prn_track.c 128 2015-03-23 16:40:14Z esterhui $
 *
 */
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "st_datatypes.h"
#include "st_dsp.h"
#include "st_io.h"

#define QUARTER_CYCLE (M_PI / 2.0) // 90 degrees (but in radians)
#define TRACKING_THRESHOLD (0.0)   // 2 sigma [V/V]

#define CMC_FIT_SEC (5.0)

#define SNPRINTF_LENGTH (512)

// This is how many accumulations need to be processed before estimating
// the noise floor
#define PROMPT_WINDOW_SIZE (200)

void displayUsage(char *progname) {
  fprintf(stderr, "%s [-a -p -n -i -d] < inputfile > outputfile\n\n", progname);
  fprintf(stderr, "Reads IF data from stdin and "
                  "writes I/Q accumulations to stdout\n\n");
  fprintf(stderr, "$Revision: 128 $\n\n");
#ifdef FITTER
  fprintf(stderr, "FITTING ENABLED\n");
#else
  fprintf(stderr, "FITTING DISABLED\n");
#endif
  fprintf(stderr, "OPTIONS\n");
  fprintf(stderr, "\t-a acqfile    :   File containing initial parameters\n");
  fprintf(stderr,
          "\t-p prnfile    :   File containing PRN to correlate against\n");
  fprintf(stderr,
          "\t-n integ_periods  :   Number of code periods to integrate over\n");
  fprintf(stderr,
          "\t-i numaccum   :   Number of integrations to perform, 0=inf\n");
  fprintf(stderr, "\t-m channel    :   If multiple channels in data, use "
                  "channel m, default m=0\n");
  fprintf(stderr, "\t-d pll,dll    :   Disable feedback for dll,pll\n");
  fprintf(stderr, "\t-s lag spacing:   Earlty to Late lag spacing in chips "
                  "(default is 1.0 chips)\n");
  fprintf(stderr, "\t-f loops      :   Run FLL for LOOPS integrations, then "
                  "switch to PLL. If loops is 0, will always use FLL\n");
  fprintf(stderr, "\t-q            :   Invert I/Q\n");

  return;
}

int main(int argc, char **argv) {

  double lag_spacing_chips = 0.5; // lag spacing in chips
  double kp = 1.0 / (1.0 + 1.0);
  // FIXME, assuming no filter effects; When using matlab
  // simulations (with no filter, but bandlimited to fs/2),
  // According to Brooks 95-6, eq 6.7 (page 6-2), kp is
  // computed for 9.6MHz BPF P-code by taking 0.9/(1.28+1.15)
  // where:
  //  0.9 = Peak correlation amplitude of P-code
  //  1.28,1.15 = Slope of derivative of correlation waveform
  //              at Early and Late correlators.
  //  For C/A code and a wide-band filter (>> 2MHz), kp
  //  effectively evaluates to 1/(1+1) = 0.5;
  unsigned int samples_to_use = 0;
  double integ_periods = 1;
  unsigned int cmc_counter = 0;
  int cmc_error = 0;
  unsigned long long int sample_counter = 0;
  int signal_locked = 1;
  char *acqfile = NULL;
  FILE *acqfid = NULL;
  char *prnfile = NULL;
  FILE *fid_datahdr = NULL;
  double integration_time_seconds = 0;
  double carrier_residual_s = 0, code_residual_chips, carrier_residual_rad;
  double carrier_error_radians = 0;
  double car_freq_orig = 0, cod_freq_orig = 0;
  double nominal_phase_cycles;
  double nominal_cycles_per_sample;
  double IE, QE, IP = 0, QP = 0, IL, QL;
  double IPM1 = 0, QPM1 = 0, cross = 0, dot = 0;
  int freq_update = 0;
  double Ae, Ap, Al;
  double vsnr = 0, noise = 0;
  double chip_offset;
  unsigned int accum_count = 0, numaccum = 0, epocs_per_sec = 0;
  double carrier_phase_s = 0, carrier_phase_s_center, carrier_phase_cycles = 0;
  double carrier_observable = 0, code_observable = 0;
  double chisq = 0;
  double *cmc = NULL, *cmc_time = NULL;
  double receivertime_s_center = 0;
  double modelchips = 0, modeltime_s_center, modelchips_center = 0;
  double old_modeltime_s_center = 0;
  double nominal_chips_per_epoch;
  unsigned long long samples_skipped;
  double tau_delta_tx = 0;
  double bitsign = 1;
  int chan_to_use = 0;
  int do_pll = 1; // By default we do only pll
  double freq_error = 0;
  double coef[4];         // For saving fitting coefficients
  double cmc_rate_hz = 0; // Code-carrier rate in Hz
  double datarate_hz = 0; // Data rate in Hz
  double cmc_fix_hz = 0;  // Frequency fix for CMC offset
  double current_cmc = 0;
  double elapsed_time_s = 0;
  unsigned long long epoch_counter =
      0; // Counts how many integrations we've done

  double car_dop;
  double initial_code_dop;
  double prompt_window[PROMPT_WINDOW_SIZE];

  int c;
  int flipiq = 0;

  complex_double_t *pdata = NULL;
  data_header_t *datahdr = NULL;
  acq_parameters_t *acqhdr = NULL;

  cr_t *pcr = NULL;
  cr_t *pcr_out = NULL;
  corr_t corr;

  uint8_t *pprn = NULL;

  double car_freq_adj = 0, cod_freq_adj = 0, car_freq_adj_old = 0;
  int disable_pll_feedback = 0;
  int disable_dll_feedback = 0;
  double rss_prompt = 0;
  int read_backwards = 0;
  int dir_sign = -1;
  unsigned long long fll_loops = 0;

  loop_t *carrier = (loop_t *)malloc(sizeof(loop_t));
  loop_t *code = (loop_t *)malloc(sizeof(loop_t));
  pcr = (cr_t *)malloc(sizeof(cr_t));

  memset(carrier, 0, sizeof(loop_t));
  memset(code, 0, sizeof(loop_t));
  memset(pcr, 0, sizeof(cr_t));

  if (pcr == NULL || carrier == NULL || code == NULL) {
    perror("Can't allocate memory");
    goto alldone;
  }

  /* --- Scan input arguments --- */

  opterr = 0;

  while ((c = getopt(argc, argv, "hqi:f:p:d:a:n:m:s:")) != -1) {
    switch (c) {
    case 'f':
      sscanf(optarg, "%llu", &fll_loops);
      do_pll = 0;
      break;
    case 'a':
      acqfile = optarg;
      if ((acqfid = fopen(acqfile, "r")) == NULL) {
        perror(acqfile);
        goto alldone;
      }

      break;
    case 'p':
      prnfile = optarg;
      break;
    case 'q':
      flipiq = 1;
      break;
    case 'd':
      if (strstr(optarg, "dll") != NULL) {
        disable_dll_feedback = 1;
      }
      if (strstr(optarg, "pll") != NULL) {
        disable_pll_feedback = 1;
      }
      break;
    case 'n':
      sscanf(optarg, "%lf", &integ_periods);
      assert(integ_periods > 0);
      break;
    case 's':
      sscanf(optarg, "%lf", &lag_spacing_chips);
      /* Since user entered Early-to-Late spacing,
         divide by 2 here to get actual E,P,L spacing */
      lag_spacing_chips /= 2.0;
      break;
    case 'i':
      sscanf(optarg, "%d", &numaccum);
      break;
    case 'm':
      sscanf(optarg, "%d", &chan_to_use);
      break;
    case 'h':
      displayUsage(argv[0]);
      goto alldone;
    default:
      displayUsage(argv[0]);
      goto alldone;
    }
  }

  // Parse data header
  char hdr_file[SNPRINTF_LENGTH];
  snprintf(hdr_file, SNPRINTF_LENGTH, "iq_meta.txt");
  if ((fid_datahdr = fopen(hdr_file, "r")) == NULL) {
    perror("Can't open header file");
    goto alldone;
  }

  datahdr = parseDataHeader(fid_datahdr);
  if (datahdr == NULL)
    goto alldone;

  // Parse acquisition header
  acqhdr = parseAcquisitionHeader(acqfid);
  if (acqhdr == NULL)
    goto alldone;
  acqhdr->fs = datahdr->fs;

  if (datahdr->rf_freq < 0) {
    read_backwards = 1;
    // Pseudorange uses this sign
    dir_sign = -1;
  } else {
    read_backwards = 0;
    // If we're reading things backwards,
    // compute pseudorange differently
    dir_sign = 1;
  }

  // Parse prn file
  if ((pprn = read_prncode(prnfile, acqhdr, read_backwards)) == NULL) {
    fprintf(stderr, "Couldn't read PRN file: %s\n", prnfile);
    goto alldone;
  }
  // FIXME: Hack for old PRSR files with incorrect RF_FREQ setting
  /*
  if ((datahdr->rf_freq<1e9 || datahdr->rf_freq>3e9) && datahdr->rf_freq>0) {
      datahdr->rf_freq=2032e6;
      fprintf(stderr,"WARNING, changing rf_freq to %lf Hz\n",
          datahdr->rf_freq);
  }
  */
  // If rf frequency is negative (ie. reading back data file backwards,
  // then flip the IF
  if (datahdr->rf_freq < 0) {
    acqhdr->intermediate_frequency =
        -(acqhdr->carrier_frequency + datahdr->rf_freq);
  } else {
    acqhdr->intermediate_frequency =
        (acqhdr->carrier_frequency - datahdr->rf_freq);
  }

  fprintf(stderr, "Nominal IF is %lf Hz\n", acqhdr->intermediate_frequency);
  // Compute number of samples in integration. Here we always try to
  // have an odd number of samples to make time-tags land on an integer
  // sample instead of 1/2 way inbetween.
  // (JPL 95-6 page 4-2 eq 4.5)
  samples_to_use = round(integ_periods * acqhdr->samples_per_code_period);
  // If this is an even number, add 1 to make integration odd
  if (samples_to_use % 2 == 0) {
    pcr->num_samples = samples_to_use + 1;
  } else {
    // If odd, recompute closest even
    samples_to_use =
        round((integ_periods * acqhdr->samples_per_code_period) / 2.0) * 2;
    pcr->num_samples = samples_to_use + 1;
  }

  // Init PLL stuff
  integration_time_seconds = 1.0 / datahdr->fs * samples_to_use;

  if (integration_time_seconds > 0.2) {
    loop_init(carrier, integration_time_seconds, 0.4, 0.4);
    loop_init(code, integration_time_seconds, 0.1, 0.1);
  } else if (integration_time_seconds > 0.08) {
    loop_init(carrier, integration_time_seconds, 2.0, 2.0);
    loop_init(code, integration_time_seconds, 0.2, 0.2);
  } else if (integration_time_seconds > 15e-3) {
    // works for 4g
    // loop_init(carrier,integration_time_seconds,18.0,5.0);
    // loop_init(code,integration_time_seconds,1.0,1.0);
    // loop_init(carrier,integration_time_seconds,28.0,6.0);
    // loop_init(code,integration_time_seconds,1.0,1.0);
    // works for 1g

    // tracks through 292 sec problem
    // loop_init(carrier,integration_time_seconds,10.0,4.0);
    // loop_init(code,integration_time_seconds,0.4,0.4);
    // works for < 1g
    loop_init(carrier, integration_time_seconds, 5.0, 2.0);
    loop_init(code, integration_time_seconds, 0.2, 0.2);
  } else if (integration_time_seconds >= 5e-3) {
    loop_init(carrier, integration_time_seconds, 14.0, 4.0);
    loop_init(code, integration_time_seconds, 14.0, 4.0);
  } else {
    loop_init(carrier, integration_time_seconds, 18.0, 5.0);
    loop_init(code, integration_time_seconds, 18.0, 5.0);
  }

  // Other init stuff
  pcr->fs = datahdr->fs;
  pcr->fc = acqhdr->doppler;
  pcr->start_phase = acqhdr->carrier_phase;

  if (acqhdr->codes_per_databit > 0) {
    datarate_hz =
        acqhdr->coderate / acqhdr->codelength / acqhdr->codes_per_databit;
  }

  epocs_per_sec = (unsigned int)ceil(datahdr->fs / samples_to_use);

  cmc = (double *)malloc(sizeof(double) * epocs_per_sec * CMC_FIT_SEC);
  cmc_time = (double *)malloc(sizeof(double) * epocs_per_sec * CMC_FIT_SEC);

  if (cmc == NULL || cmc_time == NULL) {
    perror("Can't allocate memory");
    goto alldone;
  }

  // Corr init stuff

  corr.number = 3;

  // ---- NARROW SPACING ---------
  corr.start_chips = -lag_spacing_chips; // Early correlator is at -0.25 chips
  corr.start_chips_orig = corr.start_chips;
  corr.spacing_chips = lag_spacing_chips;

  // print header -------------------------------------------------
  printAccumHeader(stderr, datahdr, acqhdr, integ_periods, lag_spacing_chips);

  // Can't integrate over data bits
  if (acqhdr->codes_per_databit > 0) {
    assert(integ_periods <= acqhdr->codes_per_databit);
  }

  // This is how many chips we have in our accumulation window nominally
  // Note, this is NOT tied to samples. Eg, if we have 1023 chips
  // per code period, and we are integrating over 20 of these,
  // it is 1023*20 chips
  nominal_chips_per_epoch = integ_periods * acqhdr->codelength;

  fprintf(stdout, "time_s  IE  QE  IP  QP  IL  QL  ");
  fprintf(stdout, "doppler_hz  coderate_hz  phase_s  range_s  ");
  fprintf(stdout, "epoch_offset_ns  num_samp  vsnr  noise\n");

  // Save off the original frequency, need this for
  // PLL to work properly
  car_freq_orig = pcr->fc;

  // if (read_backwards) acqhdr->coderate*=-1;
  cod_freq_orig = acqhdr->coderate;

  /* Compute carrier doppler, scale to code doppler (Ftx/Fchip) */
  car_dop = dir_sign * (pcr->fc - acqhdr->intermediate_frequency);
  initial_code_dop = dir_sign * car_dop / (datahdr->rf_freq / cod_freq_orig);
  if (acqhdr->intermediate_frequency == 2.558e6) {
    fprintf(stderr, "WARNING!!!! Using ISS/connect HACK for L5 data\n");
    initial_code_dop = -initial_code_dop;
  }
  acqhdr->coderate += initial_code_dop;

  // This is the time in seconds that a perfect integration (over n-epocs
  // will take)
  nominal_cycles_per_sample = acqhdr->intermediate_frequency / acqhdr->fs;

  samples_skipped = acqhdr->sample_number;

  // Account for the start offset of the prompt correlator,
  // if we are off by, say 0.5 chips, the start time is not zero
  // seconds!

  // Setup start cycles to be tied to clock
  // carrier_phase_cycles = acqhdr->sample_number * nominal_cycles_per_sample;
  // carrier_phase_cycles -= floor(carrier_phase_cycles);


  fprintf(stderr,"Sky freq is %f\n",acqhdr->carrier_frequency);
  // ------- Loop through all data and run prcessing on blocks
  while (((pdata = parseDataPayload(stdin, datahdr, acqhdr->sample_number,
                                    pcr->num_samples, chan_to_use, flipiq)) !=
          NULL) &&
         signal_locked && (accum_count < numaccum || numaccum == 0)) {

    // Advance the data pointer
    pcr->data = pdata;

    // counter rotation ------------------------------------------
    pcr_out = counterRotate(pcr);

    // remove PRN code & acucm------------------------------------
    prnMultiplyAccumulate(&corr, pcr_out, pprn, acqhdr);
    accum_count++;

    IE = corr.data[0].real;
    QE = corr.data[0].imag;
    IP = corr.data[1].real;
    QP = corr.data[1].imag;
    IL = corr.data[2].real;
    QL = corr.data[2].imag;

    free(corr.data);

    // Keep a 'moving window' of 100 samples of rss(IP,QP)
    rss_prompt = sqrt(IP * IP + QP * QP);

    prompt_window[(accum_count - 1) % PROMPT_WINDOW_SIZE] = rss_prompt;

    if (accum_count >= PROMPT_WINDOW_SIZE) {
      vsnr = compute_snr(rss_prompt, prompt_window, PROMPT_WINDOW_SIZE, &noise);
    }

    // Observable extraction -----------------------------------------
    // FIXME: All this shit below should go into a
    // computeObservable(tr) function,
    // where tr is at which UTC time (or sample time?)
    // to evaluate the observable at

    // -------------------------------------------------------------------
    // ----------------------------- Code phase --------------------------
    // -------------------------------------------------------------------

    // Now take integrated phase and subtract out nominal phase,
    // divide by chip rate to convert to seconds
    carrier_error_radians = atan2(QP, IP);
    /* Compute early, prompt, late ampl */
    Ae = cos(carrier_error_radians) * IE + sin(carrier_error_radians) * QE;
    Ap = cos(carrier_error_radians) * IP + sin(carrier_error_radians) * QP;
    Al = cos(carrier_error_radians) * IL + sin(carrier_error_radians) * QL;

    /* Compute tracking error */
    /* Brooks 95-6, eq 6.7 (page 6-2)
     */
    code_residual_chips = kp * (Al - Ae) / (Ap);

    // Compute model at center of INTEGRATION window
    modelchips_center = modelchips + ((acqhdr->chips_per_sample) *
                                      (pcr->num_samples - 1) / 2.0);

    // Compute model time in seconds at center of INTEGRATION,
    // Here we also add in the residual
    modeltime_s_center =
        (modelchips_center - code_residual_chips) / cod_freq_orig;

    if (read_backwards)
      modeltime_s_center *= -1.0;

    // Keep track of how many chips have gone by in our model
    modelchips += corr.integrated_phase;

    // -------------------------------------------------------------------
    // -------------------------- Carrier phase --------------------------
    // -------------------------------------------------------------------
    nominal_phase_cycles =
        nominal_cycles_per_sample * (double)(pcr->num_samples);

    /* Compute dot and cross product,
     * used for FLL later as well
     * as to determine if we had a bitflip */
    cross = (IPM1 * QP) - (IP * QPM1);
    dot = (IPM1 * IP) + (QPM1 * QP);

    // Correct for bitflips in phase, this might
    // not work with low SNR where we get ~180deg change
    // in phase due to noise...
    if (acqhdr->codes_per_databit > 0) {
      // This method works if PLL
      // bitsign=IP>0?+1:-1;
      // if the dot product is positive, no bit flip
      // if it is negative, we had a bitflip
      bitsign *= dsign(dot);
    }

    if (IP != 0) {
      carrier_residual_rad = atan(QP / IP);
    } else {
      carrier_residual_rad = 0;
    }

    carrier_residual_s =
        carrier_residual_rad / (2.0 * M_PI) / acqhdr->carrier_frequency;

    // This is the accumulated phase at the END of the interval
    // FIXME: If we start counter rotation with a phase offset,
    // this is not reflected here in 'time', should add start
    // phase at some point
    carrier_phase_cycles += (nominal_phase_cycles - pcr_out->integrated_phase);

    // Accumulated phase at end of interval (in seconds)
    carrier_phase_s = carrier_phase_cycles / acqhdr->carrier_frequency;

    // Compute at center of interval the phase
    carrier_phase_s_center =
        carrier_phase_s -
        ((((nominal_phase_cycles - pcr_out->integrated_phase) /
           acqhdr->carrier_frequency)) /
         2.0);

    // Time tags lie at center of correlation interval
    receivertime_s_center =
        (samples_skipped + sample_counter + ((pcr->num_samples - 1) / 2.0)) /
        datahdr->fs;

    // Keep track of elapsed time
    elapsed_time_s = sample_counter / datahdr->fs;

    // If we're reading backwards this is actually negative time
    if (read_backwards) {
      receivertime_s_center *= -1.0;
    }

    // This is the offset of our integration center
    // to the ACTUAL center of the epoch. For instance,
    // this value is subtracted from GRAIL NowTime(), NowTime()
    // is reported at the CENTER of the bit. tau_delta_tx computes
    // the offset then from the center of the bit to where we made our
    // observation
    // (spacecraft_nowtime_center) - (model_time_center)
    // Another way to view this: we compute the offset of the
    // integration center (modeltime_s_center) and the epoch center
    // (which will be the model_time evaluated at the bit center)
    tau_delta_tx =
        // We do it funny here, first compute after, say 20 epocs
        // what the txtime is, then subtract out 1/2 a epoch
        // to center it at our epoch
        (((double)(epoch_counter + 1) * nominal_chips_per_epoch -
          (nominal_chips_per_epoch / 2.0)) /
         (dir_sign * cod_freq_orig)) -
        modeltime_s_center;

    sample_counter += pcr->num_samples;
    epoch_counter += 1;

    /* After about 20 FLL updates
     * we should be ready to do PLL */
    if (epoch_counter == fll_loops && fll_loops != 0) {
      do_pll = 1;
    }

    carrier_observable = (carrier_phase_s_center - carrier_residual_s);
    code_observable = -(modeltime_s_center - receivertime_s_center);

    // Save carrier-code in
    cmc[cmc_counter] = (code_observable - carrier_observable);
    cmc_time[cmc_counter] = receivertime_s_center;
    // Only start doing CMC after about 5 seconds of tracking
    if (elapsed_time_s > 5.0) {
      cmc_counter++;
    }
    // Approx every 10 sec do a fit
    if (cmc_counter >= CMC_FIT_SEC * epocs_per_sec && elapsed_time_s < 30.0) {
#ifdef FITTER
      (void)fitter(cmc_time, cmc, CMC_FIT_SEC * epocs_per_sec, 1, &chisq, coef,
                   cmc_time[0], 0);
#else
      coef[0] = 0;
      coef[1] = 0;
      coef[2] = 0;
      coef[3] = 0;
      chisq = 0;
#endif
      // Compute the rate of code-carrier in Hz
      cmc_rate_hz = acqhdr->carrier_frequency * coef[1];
      current_cmc = cmc_rate_hz;
      // fprintf(stderr,"CMC rate is %g Hz\n",cmc_rate_hz);
      //  If we have data bits, make sure we're not
      //  on a false lock (eg +-25 Hz off for GPS)
      if (datarate_hz > 0) {
        if (fabs(cmc_rate_hz) > (0.2 * (datarate_hz / 2.0)) &&
            fabs(cmc_rate_hz) < (1.8 * (datarate_hz / 2.0))) {
          // If we already fixed cmc but the rate is still off,
          // print out error
          if ((int)cmc_fix_hz != 0) {
            fprintf(stderr, "ERROR: CMC not fixed, %g Hz!\n", cmc_fix_hz);
            cmc_error = 1;
          } else {
            fprintf(stderr, "CMC rate is %g Hz, fixing tracking\n",
                    cmc_rate_hz);
            cmc_fix_hz = -cmc_rate_hz;
          }
        }
      } else {
        if (fabs(cmc_rate_hz) > 5.0) {
          // If we already fixed cmc but the rate is still off,
          // print out error
          if ((int)cmc_fix_hz != 0) {
            fprintf(stderr, "ERROR: CMC not fixed (no bits), %g Hz!\n",
                    cmc_fix_hz);
            cmc_error = 1;
          } else {
            fprintf(stderr,
                    "CMC rate is %g Hz, fixing tracking (no data bits)\n",
                    cmc_rate_hz);
            cmc_fix_hz = -cmc_rate_hz;
          }
        }
      }
      cmc_counter = 0;
    }

    if (cmc_counter >= CMC_FIT_SEC * epocs_per_sec) {
      cmc_counter = 0;
    }

    // write data to stdout -----------------------------
    fprintf(stdout,
            "% 12.12f % 12.2lf % 12.2lf % 12.2lf % 12.2lf % 12.2lf % 12.2lf",
            receivertime_s_center, IE, QE, IP, QP, IL, QL);
    fprintf(stdout,
            "% 15lf % 15lf % 20.15lf % 20.15lf % 9.3lf % 8d % 5.1lf % 6.1lf\n",
            dir_sign * (pcr->fc - acqhdr->intermediate_frequency),
            (acqhdr->coderate - cod_freq_orig),

            // Integrated phase
            carrier_observable,

            // Pseudorange = modeltime(tR) - receivertime(tR)
            // tR = time at center of integration interval
            code_observable,

            // This is the model time  (time as transmitted by satellite) for
            // the actual center of integration window,
            // not center of epoch. This way it aligns with our other
            // observables
            (tau_delta_tx * 1.0e9),
            // pcr->num_samples, vsnr, noise);
            pcr->num_samples, vsnr, cmc_rate_hz);
    // pcr->num_samples, vsnr, noise);

    car_freq_adj_old = car_freq_adj;

    // FLL/PLL -------------------------------------------------------
    // Compute phase error (Prompt correlator)
    if (do_pll) {
      carrier->phase_error = atan(QP / IP) / (2.0 * M_PI);
      carrier->freq_error = 0;
    } else {
      // If this is the first time, just say we have zero error
      if (IPM1 == 0 && QPM1 == 0) {
        carrier->phase_error = 0;
        carrier->freq_error = 0;
        freq_update = 0;
      } else {
        freq_update = (freq_update + 1) % 2;
        // Only update every other time
        // since the freq needs to be the
        // same for 2 accumulations to
        // estimate freq.
        if (freq_update == 1) {
          // FLL assumes that vector does NOT stradle
          // a bit transition. For weak signal
          // detection where we have to integrate
          // over 1 bit period, this will happen
          // all the time. Thus we correct for bit flips
          // here at the expense of halfing our pull-in range
          //  IE. for 20ms (50Hz) integrations our pull in should
          //  +- 25Hz, but with this bit correction it is
          //  +-12.5 Hz.
          /* correct for bit flips */
          cross *= dsign(dot);
          dot *= dsign(dot);
          // Compute error in cycles/second
          // If we have a bit flip
          freq_error = (atan2(cross, dot)) /
                       (modeltime_s_center - old_modeltime_s_center) /
                       (2.0 * M_PI);
          carrier->freq_error = freq_error;
          carrier->phase_error = 0;
          car_freq_adj = loop_run3(carrier);
        }
      }
    }
    // Compute code error

    code->phase_error =
        (sqrt((IE * IE) + (QE * QE)) - sqrt((IL * IL) + (QL * QL))) /
        (sqrt((IE * IE) + (QE * QE)) + sqrt((IL * IL) + (QL * QL)));

    // Compute adjustment frequencies for PLL and DLL
    if (!disable_pll_feedback && do_pll)
      car_freq_adj = loop_run3(carrier);

    if (!disable_dll_feedback)
      cod_freq_adj = loop_run3(code);

    // Update phase and frequency --------------------------------
    pcr->start_phase = pcr_out->phase_next;
    pcr->fc = car_freq_orig + car_freq_adj + cmc_fix_hz;

    // Change integration length by one sample -------------------
    // to keep us aligned with code epoch
    // (if necessary!)
    // FIXME: For glonass P-code, codelength should
    // probably be codelength*numperiods to properly
    // adjust every 1msec or 4msec instead of every 1sec
    chip_offset = -fmod(corr.chips_next - corr.start_chips_orig,
                        (double)(acqhdr->codelength));
    if (chip_offset > (double)acqhdr->codelength / 2.0) {
      chip_offset -= acqhdr->codelength;
    } else if (chip_offset < (double)(-acqhdr->codelength / 2.0)) {
      chip_offset += acqhdr->codelength;
    }

    if (chip_offset < 0)
      pcr->num_samples = samples_to_use + 1;
    else
      pcr->num_samples = samples_to_use - 1;

    IPM1 = IP;
    QPM1 = QP;
    old_modeltime_s_center = modeltime_s_center;
    corr.start_chips = corr.chips_next;
    // FIXME: This is bad practice, changing acqusition frequencies!
    acqhdr->coderate = cod_freq_orig + (initial_code_dop + cod_freq_adj);
    acqhdr->sample_number = 0; // Don't skip samples for next iter.
    // Whenever coderate changes, we have to update parameters
    calc_prn_params(acqhdr);

    free(pcr_out->data);
    free(pcr_out);
    free(pdata);

    // Lock detect.
    // - VSNR is good
    // - Skip check if feedback is disabled
    // - code-carrier rate is good
    if ((((vsnr < TRACKING_THRESHOLD) && (accum_count >= PROMPT_WINDOW_SIZE)) &&
         (disable_pll_feedback == 0 && disable_dll_feedback == 0)) ||
        (cmc_error)) {
      fprintf(stderr, "WARNING: Loss of lock, vsnr = %0.1f cmc_error is %d\n",
              vsnr, cmc_error);
      signal_locked = 0;
    }
  }

alldone:
  if (pdata != NULL)
    free(pdata);
  if (datahdr != NULL)
    free(datahdr);
  if (pcr != NULL)
    free(pcr);
  if (carrier != NULL)
    free(carrier);
  if (code != NULL)
    free(code);
  if (pprn != NULL)
    free(pprn);
  if (acqhdr != NULL)
    free(acqhdr);
  if (cmc != NULL)
    free(cmc);
  if (cmc_time != NULL)
    free(cmc_time);

  //  if (corr.data!=NULL) free(corr.data);

  return 0;
}
