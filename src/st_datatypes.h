#ifndef __ST_DATATYPES_H
#define __ST_DATATYPES_H

#include <inttypes.h>

#include "armadillo"

typedef std::vector<double> stdvec;

#define DATA_BLOCKSIZE 524288
#define PRN_LENGTH 1023

using namespace arma;

typedef struct {
  float real;
  float imag;
} complex_double_t;

/* -------------------------------------------------------------------------- */
/** @brief Parameters used for acquisition stored in the struct
 *
 * Some of these parameters are filled in by read_prn()
 */
/* -------------------------------------------------------------------------- */
typedef struct {
  double fs;                /**< Sample frequency in Hz */
  double coderate;          /**< Code rate in Hz */
  double carrier_frequency; /**< Signal (true) TX Carrier frequency in Hz */
  double intermediate_frequency;  /**< IF at sampled of carrier in Hz*/
  uint32_t codelength;            /**< Code length, eg 1023 */
  unsigned int codes_per_databit; /**< Number of PN codes in
                                  1 databit */
  double samples_per_code_period; /**< How many samples per code
                                    perdiod we have */
  double samples_per_chip;        /**< How many samples per chip */
  double chips_per_sample;        /**< How many chips per sample */
  uint32_t bytes_per_code_period; /**< How many integer bytes per
                                    code period we have */
  uint32_t prn;                   /**< What PRN code this is */
  char description[160];          /**< ASCII description */
  double vsnr_threshold;          /**< Raw detection Voltage SNR */
  double doppler_search_space_ss; /**< Single sided search space, Hz */

  // ---- Post acquisition params
  int num_lags;
  double lag_spacing_s;   /**< lag spacing in seconds */
  double *mag_vs_lag;     /**< rss(Iaccum,Qaccum) vs vs. lag of doppler bin
                               where the peak magnitude was detected */
  double *phase_vs_lag;   /**< atan2(Qaccum,Iaccum) vs. lag
                             where the peak magnitude was detected */
  double noise_floor;     /**< An estimate of the noise floor, measured
                            in 1/integration_period bandwidth */
  double *rss_vs_lag_I;   /**< atan2(Qaccum,Iaccum) vs. lag  */
  double *rss_vs_lag_Q;   /**< atan2(Qaccum,Iaccum) vs. lag  */
  uint8_t found;          /**< Signal found? 1=found 0=not found */
  uint64_t sample_number; /**< At what sample number code delay = 0 chips
                           */
  double magnitude;       /**< Peak output of rss(Iaccum,Qaccum); */
  double vsnr_1s;         /**< 1s Voltage SNR (estimated by computing
                                the ratio of the  rss(correlation peak)
                                to an average over 0.25*codelength noise
                                lags and using the integration time
                                to scale to 1s */
  double vsnr;            /**< Raw Voltage SNR (estimated by computing
                                   the ratio of the  rss(correlation peak)
                                   to an average over 0.25*codelength noise
                                   lags */
  double code_phase;      /**< Phase of code at sample 0 */
  double carrier_phase;   /**< Phase of carrier at sample 0 in rad */
  double
      doppler; /**< Doppler of detected signal (this includes the IF) in Hz */
  double doppler_bin_size_hz; /**< Size of doppler bins searched in Hz */
} acq_parameters_t;

/* -------------------------------------------------------------------------- */
/** @brief Header for IF data
 *
 */
/* -------------------------------------------------------------------------- */
typedef struct {
  double fs;                      /**< Sample frequency in Hz */
  double rf_freq;                 /**< RF Sky frequency (eg 1575.112e6) in Hz */
  unsigned int year;              /**< time tag - year */
  unsigned int doy;               /**< time tag - day of year */
  unsigned int sec;               /**< time tag - second of day */
  int data_time_offset;           /**< in nano-seconds */
  int bits_per_sample;            /**< bits per sample*/
  int num_channels;               /* number of IF channels in data file */
  unsigned int data_size_samples; /**< Payload data length in samples */
  size_t data_size_bytes;         /**< Payload data length in bytes */
  char description[160];          /**< One line description */
} data_header_t;

typedef struct {
  data_header_t *header;
  complex_double_t *data;
} chan_t;

typedef struct {
  complex_double_t *data;
  unsigned int num_samples;
  double fs;
  double fc;
  double start_phase;      /**< Phase to start rotation with */
  double integrated_phase; /**< Integrated phase for interval in cycles*/
  double phase_next;       /**< Phase for next sample */
} cr_t;

typedef struct {
  complex_double_t *data;
  unsigned int number;     /**< Number of correlators */
  double start_chips;      /**< First correlator start position in chips */
  double spacing_chips;    /**< Correlator spacing in cips */
  double chips_next;       /**< Filled in by prnMultiplyAccumulate */
  double start_chips_orig; /**< Original start chips */
  double integrated_phase; /**< Integrated phase over interval */
} corr_t;

typedef struct {
  double phase_error;
  double freq_error;

  double T; // integration time in seconds

  /* 3rd order PLL */
  double Wn3p;
  double Wn2p;
  double Wnp;
  double a3;
  double b3;

  /* 2nd order FLL */
  double Wnf;
  double Wn2f;
  double a2;

  /* integrators */
  double acc_delay;
  double vel_delay;
} loop_t;

typedef struct {
  vec timetag_s;
  vec ip;
  vec qp;
  vec doppler_hz;
  vec phase_s;
  vec nsamp;
  double fs; // sample rate in Hz
} st_accum_tone;

typedef struct {
  int antenna_number;
  double rf_frequency_hz;
  /* ---- Amplitude related ----- */
  uint32_t num_accum; // number of accumulations
  double rss_mean;
  double rss_std;
  double integ_time_s;
  double vsnr1s;
  double cn0;
  /* ---- Phase related ----- */
  double doppler_mean_hz;
  double doppler_std_hz;
  double phase_ns;        // phase of quad fit (mean phase)
  double phase_rate;      // phase rate of quad fit
  double phase_rate_rate; // phase rate rate of quad fit
  double phase_std_ns;    // phase scatter (after removing poly)
  double phase_vsnr1s;    // phase snr
} stats_signal;

typedef struct {
  int antenna_number;
  double rf_frequency_hz;
  uint32_t num_samples; // number of IF samples used
  double ratio[32];     // ratio of number of ones to total samples counted,
                        // should be 0.5
} stats_ifsamples;

#endif
