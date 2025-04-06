/* tone tracker
 *
 *
 * $Id: tone_track.c 120 2013-07-24 00:14:14Z esterhui $
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

#define HDR_DIR ("/home/esterhui/.signalTracker/")

#define PROMPT_WINDOW_SIZE (100)
#define TRACKING_THRESHOLD (0.0) // 2 sigma [V/V]

void displayUsage(char *progname) {
  fprintf(stderr, "%s [options] < inputfile > outputfile\n\n", progname);
  fprintf(stderr, "Reads IF data formatted with prsr_parse from stdin and "
                  "writes I/Q accumulations to stdout\n\n");
  fprintf(stderr, "OPTIONS\n");
  fprintf(stderr, "\t-a seconds,vsnr:   Assist PLL by first doing a FFT to get "
                  "frequency\n");
  fprintf(stderr, "\t                  FFT will use 'seconds' of data, eg 0.1 "
                  "will use 100ms\n");
  fprintf(stderr,
          "\t                  FFT will use 'vsnr' as threshold, eg 6 V/V\n");
  fprintf(stderr, "\t-f tx,dop     :   Nominal TX frequency, and\n");
  fprintf(stderr,
          "\t              :   doppler estimate in Hz eg. 1575420000,2000\n");
  fprintf(
      stderr,
      "\t-l phase,freq :   phase/freq loop bandwidth in (Hz), 18,5 default\n");
  fprintf(stderr, "\t-n numsamp    :   Samples to integrate over\n");
  fprintf(stderr, "\t-d            :   Disable feedback\n");
  fprintf(stderr,
          "\t-i numaccum   :   Number of integrations to perform, 0=inf\n");
  fprintf(stderr, "\t-m channel    :   If multiple channels in data, use "
                  "channel m, default m=0\n");
  fprintf(stderr, "\t-s seektime   :   Seek this many seconds into data\n");
  fprintf(stderr,
          "\t                  before starting tracking. Default=0.0 \n");
  fprintf(stderr, "\t-c dop_freq   :   If -a (FFT) option used, will reject "
                  "tones with doppler\n");
  fprintf(stderr,
          "\t              :   greather than dop_freq offset from IF.\n");
  fprintf(stderr, "\t-q            :   Swap I/Q\n");
  fprintf(
      stderr,
      "\t-w headerfile :   Use this header file (describing data format)\n");

  return;
}

int main(int argc, char **argv) {
  int numsamp = 1000;
  double integration_time_seconds = 0;
  double seektime = 0;
  double car_freq_orig = 0;
  int chan_to_use = 0;
  double nominal_freq = 0, if_freq = 0, dop_freq = 0;
  unsigned int numaccum = 0;
  double IP = 0, QP = 0;
  double car_freq_adj = 0;
  double freq_error = 0;
  int signal_locked = 1;
  double nominal_phase_cycles;
  int do_pll = 1;
  double carrier_phase_s = 0, carrier_phase_s_center, carrier_phase_cycles = 0;
  double carrier_observable = 0;
  double receivertime_s_center = 0, receivertime_s_center_old = 0;
  double carrier_residual_s = 0, carrier_residual_rad;
  double nominal_cycles_per_sample;
  unsigned long long samples_skipped;
  double prompt_window[PROMPT_WINDOW_SIZE];
  double vsnr = 0, noise = 0;
  double rss_prompt = 0;
  double IPM1 = 0, QPM1 = 0, cross = 0, dot = 0;
  unsigned long long int sample_counter = 0;
  unsigned int accum_count = 0;
  double delta_time = 0;
  double lbw_phase = 40.0, lbw_freq = 5.0;

  complex_double_t *pdata = NULL;
  data_header_t *phdr = NULL;
  acq_parameters_t acqhdr;
  acq_parameters_t *ahdr = &acqhdr;

  cr_t *pcr = NULL;
  cr_t *pcr_out = NULL;

  complex_double_t accum;
  int disable_feedback = 0;

  loop_t *carrier = (loop_t *)malloc(sizeof(loop_t));
  pcr = (cr_t *)malloc(sizeof(cr_t));

  /* --- Scan input arguments --- */
  int c, rc;

  opterr = 0;

  int flipiq = 0;
  int do_fftacquire = 0;
  float fftacquire_seconds = 0.001, fftacquire_threshold = 6;
  double doppler_cutoff_hz = 0;
  int headerspecified = 0;

  char hdr_file[512];

  while ((c = getopt(argc, argv, "qhda:n:k:m:f:i:l:c:s:w:")) != -1) {
    switch (c) {
    case 'd':
      disable_feedback = 1;
      break;
    case 'q':
      flipiq = 1;
      break;
    case 'f':
      rc = sscanf(optarg, "%lf,%lf", &nominal_freq, &dop_freq);
      if (rc != 2) {
        fprintf(stderr, "Couldn't parse -f flag\n");
        fprintf(stderr, "should look like 1575420000,1000\n");
        return -1;
      }
      break;
    case 's':
      sscanf(optarg, "%lf", &seektime);
      break;
    case 'l':
      rc = sscanf(optarg, "%lf,%lf", &lbw_phase, &lbw_freq);
      if (rc != 2) {
        fprintf(stderr, "Couldn't parse -l flag\n");
        fprintf(stderr, "should look like loop_gain,loop_bandwidth\n");
        fprintf(stderr, "For instance -l 0.707,10 (10 Hz bandwidth)\n");
        return -1;
      }
      break;
    case 'i':
      sscanf(optarg, "%d", &numaccum);
      break;
    case 'c':
      sscanf(optarg, "%lf", &doppler_cutoff_hz);
      doppler_cutoff_hz = fabs(doppler_cutoff_hz);
      break;
    case 'n':
      sscanf(optarg, "%d", &numsamp);
      break;
    case 'm':
      sscanf(optarg, "%d", &chan_to_use);
      break;
    case 'w':
      headerspecified = 1;
      strncpy(hdr_file, optarg, 510);
      break;
    case 'a':
      rc = sscanf(optarg, "%f,%f", &fftacquire_seconds, &fftacquire_threshold);
      if (rc != 2) {
        fprintf(stderr, "Couldn't parse -a flag\n");
        fprintf(stderr, "should look like 0.020,30\n");
        return -1;
      }
      do_fftacquire = 1;
      break;
    case 'h':
      displayUsage(argv[0]);
      exit(0);
    default:
      displayUsage(argv[0]);
      exit(0);
    }
  }

  // Parse data header
  FILE *fid_datahdr;

  if (!headerspecified) {
    sprintf(hdr_file, "%s/data_hdr_chan%02d.txt", HDR_DIR, chan_to_use);
    fprintf(stderr, "WARNING: Using external header file %s\n", hdr_file);
  }
  if ((fid_datahdr = fopen(hdr_file, "r")) == NULL) {
    perror("Can't open header file");
  }
  phdr = parseDataHeader(fid_datahdr);
  if (phdr == NULL)
    return -1;

  // Init PLL stuff
  integration_time_seconds = 1.0 / phdr->fs * numsamp; // Compute integ time

  loop_init(carrier, integration_time_seconds, lbw_phase, lbw_freq);

  if_freq = nominal_freq - phdr->rf_freq;

  // Fill acquisition header
  ahdr->fs = phdr->fs;
  ahdr->carrier_frequency = nominal_freq;
  ahdr->intermediate_frequency = if_freq;
  ahdr->doppler = if_freq + dop_freq;
  ahdr->codes_per_databit = 0;
  fprintf(stderr, "seektime is %lf, fs is %lf\n", seektime, ahdr->fs);
  ahdr->sample_number = (uint64_t)roundl(seektime * ahdr->fs);

  samples_skipped = ahdr->sample_number;

  ahdr->carrier_phase = 0;
  ahdr->vsnr_threshold = fftacquire_threshold;
  ahdr->doppler_search_space_ss = doppler_cutoff_hz;

  /* ---- FFT to refine solution ---- */
  if (do_fftacquire) {
    int NFFT = (int)round(phdr->fs * fftacquire_seconds);
    pdata = parseDataPayload(stdin, phdr, ahdr->sample_number, NFFT,
                             chan_to_use, flipiq);
    samples_skipped += NFFT;
    fft_acquire_tone(pdata, NFFT, ahdr);
    if (ahdr->found) {
      fprintf(stderr,
              "FFT Acquire found tone at %0.1f Hz with 1sVSNR of %0.1f V/V "
              "(%0.1f dB-Hz)\n",
              ahdr->doppler, ahdr->vsnr_1s, 20.0 * log10(ahdr->vsnr_1s));
    } else {
      fprintf(stderr,
              "No tones found above threshold of %0.1f V/V with integration of "
              "%f seconds\n",
              ahdr->vsnr_threshold, fftacquire_seconds);
      free(phdr);
      free(pcr);
      return -1;
    }
  }

  pcr->num_samples = numsamp;

  pcr->fs = phdr->fs;
  pcr->fc = ahdr->doppler;
  pcr->start_phase = ahdr->carrier_phase;

  // print header -------------------------------------------------
  printAccumHeader(stdout, phdr, ahdr, numsamp, 0);

  fprintf(stdout, "# timetag:s, IP, QP,");
  fprintf(stdout, "Doppler:Hz, Phase:s, ");
  fprintf(stdout, "samplesused, vsnr:V/V, noise\n");

  // Save off the original frequency, need this for
  // PLL to work properly
  car_freq_orig = pcr->fc;

  // This is the time in seconds that a perfect integration (over n-epocs
  // will take)

  nominal_cycles_per_sample = ahdr->intermediate_frequency / ahdr->fs;

  // Account for the start offset of the prompt correlator,
  // if we are off by, say 0.5 chips, the start time is not zero
  // seconds!

  double epoch_counter = 0; // Counts how many integrations we've done

  // ------- Loop through all data and run acq. on blocks
  while (((pdata = parseDataPayload(stdin, phdr, ahdr->sample_number,
                                    pcr->num_samples, chan_to_use, flipiq)) !=
          NULL) &&
         signal_locked && (accum_count < numaccum || numaccum == 0)) {

    // Advanced the data pointer
    pcr->data = pdata;

    // counter rotation ------------------------------------------
    pcr_out = counterRotate(pcr);

    // accumulation ----------------------------------------------
    accum = accumulate(pcr_out->data, pcr_out->num_samples);

    accum_count++;

    IP = accum.real;
    QP = accum.imag;

    // Keep a 'moving window' of 100 samples of rss(IP,QP)
    rss_prompt = sqrt(IP * IP + QP * QP);

    prompt_window[(accum_count - 1) % PROMPT_WINDOW_SIZE] = rss_prompt;

    if (accum_count >= PROMPT_WINDOW_SIZE) {
      vsnr = compute_snr(rss_prompt, prompt_window, PROMPT_WINDOW_SIZE, &noise);
    }

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

    if (IP != 0) {
      carrier_residual_rad = atan(QP / IP);
    } else {
      carrier_residual_rad = 0;
    }

    carrier_residual_s =
        carrier_residual_rad / (2.0 * M_PI) / ahdr->carrier_frequency;

    // This is the accumulated phase at the END of the interval
    // FIXME: If we start counter rotation with a phase offset,
    // this is not reflected here in 'time', should add start
    // phase at some point
    carrier_phase_cycles += (nominal_phase_cycles - pcr_out->integrated_phase);

    // Accumulated phase at end of interval (in seconds)
    carrier_phase_s = carrier_phase_cycles / ahdr->carrier_frequency;

    // Compute at center of interval the phase
    carrier_phase_s_center =
        carrier_phase_s -
        ((((nominal_phase_cycles - pcr_out->integrated_phase) /
           ahdr->carrier_frequency)) /
         2.0);

    // Time tags lie at center of correlation interval
    receivertime_s_center =
        (samples_skipped + sample_counter + ((pcr->num_samples - 1) / 2.0)) /
        phdr->fs;

    sample_counter += pcr->num_samples;
    epoch_counter += 1;

    carrier_observable = (carrier_phase_s_center - carrier_residual_s);

    // write data to stdout -----------------------------
    fprintf(stdout, "% 12.12f % 12.2lf % 12.2lf", receivertime_s_center, IP,
            QP);
    fprintf(stdout, "% 15lf % 20.15lf % 8d % 5.1lf % 6.1lf\n",
            (pcr->fc - ahdr->intermediate_frequency),

            // Integrated phase
            carrier_observable,
            // pcr->num_samples, vsnr, noise);
            pcr->num_samples, vsnr, freq_error);

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
      } else {
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
        /* FIXME, taking out
        cross*=dsign(dot);
        dot*=dsign(dot);
        */
        // Compute error in cycles/second
        // If we have a bit flip
        delta_time = receivertime_s_center - receivertime_s_center_old;
        freq_error = (atan2(cross, dot)) / (delta_time) / (2.0 * M_PI);
        carrier->freq_error = freq_error;
      }
    }

    // Compute adjustment frequencies for PLL and DLL
    if (!disable_feedback)
      car_freq_adj = loop_run3(carrier);

    // Update phase and frequency --------------------------------
    pcr->start_phase = pcr_out->phase_next;
    pcr->fc = car_freq_orig + car_freq_adj;

    IPM1 = IP;
    QPM1 = QP;
    // FIXME: This is bad practice, changing acqusition frequencies!
    ahdr->sample_number = 0; // Don't skip samples for next iter.

    // Lock detect.
    // - VSNR is good
    // - Skip check if feedback is disabled
    // - code-carrier rate is good
    if ((((vsnr < TRACKING_THRESHOLD) && (accum_count >= PROMPT_WINDOW_SIZE)) &&
         (disable_feedback == 0))) {
      fprintf(stderr, "WARNING: Loss of lock, vsnr = %0.1f\n", vsnr);
      signal_locked = 0;
    }

    receivertime_s_center_old = receivertime_s_center;

    free(pcr_out->data);
    free(pcr_out);
    free(pdata);
  }

  fprintf(stderr, "Done processing\n");
  free(phdr);
  free(pcr);

  return 0;
}
