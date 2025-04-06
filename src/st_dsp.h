#ifndef __ST_DSP_H__
#define __ST_DSP_H__

/**
 * @file dsp.h
 * @brief Prototypes for signaltracking DSP
 * @author Stephan Esterhuizen
 * @date 2011-08-23
 */

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "armadillo"
#include "st_datatypes.h"
#include <fftw3.h>

/* -------------------------------------------------------------------------- */
/** @brief resample_prncode
 *
 * This function resamples PRN code from one sample per chip to N-samples
 * per chip. Where N is fs/coderate
 *
 * @param prn_in
 * @param codelength Code length in integer chips (length of prn_in vector)
 * @param numsamples number of desired samples at rate fs
 * @param fs Sample rate in Hz
 * @param coderate Code rate in Hz
 *
 * @return Pointer to +1 -1 int8_t resampled vector containing numsamples items
 */
/* -------------------------------------------------------------------------- */
int8_t *resample_prncode(uint8_t *prn_in, uint32_t codelength,
                         uint32_t numsamples, double fs, double coderate);

/* -------------------------------------------------------------------------- */
/** @brief fft_acquire
 *
 * This software uses a hybrid FFT acquisition method, here is the scoop:
 *  In order to acquire a PN code, we have to resolve both doppler and
 *  code delay uncertainties. If we happen to multyply with a correctly
 *  delayed PN code, the GRAIL signal will still contain
 *      - Carrier + doppler + data bits
 *  The FFT of the above described signal will be effectively a tone modulated
 *  with, say, 50Hz of data. This is how this acquisition algorithm works:
 *      1) Multiply,say, 1ms of sampled data with all possible
 *      combinations of PN code. If we have a short PN sequence (eg 1023 chips)
 *      this will be 2046 multiplies * 4000 samples (assuming 4000 samples
 *      in 1ms)
 *      2) Now take each of these 2046 vectors, FFT them and find a peak. As
 *      long as a data bit transition didn't occur and the signal is present,
 *      a peak will be evident at the correct doppler and code delay.
 *
 * @param pbuf_data
 * @param pbuf_prncode
 * @param numsamples
 * @param acq
 *
 * @return - 1 on success
 */
/* -------------------------------------------------------------------------- */
int fft_acquire(complex_double_t *pbuf_data, int8_t *pbuf_prncode,
                int numsamples, acq_parameters_t *acq);

/* -------------------------------------------------------------------------- */
/** @brief counterRotate
 *
 * Mixes data with reference complex sinusoid, returns structure with
 * counter-rotated data.
 *
 * @param pcr
 *
 * @return
 */
/* -------------------------------------------------------------------------- */
cr_t *counterRotate(cr_t *pcr);

/* -------------------------------------------------------------------------- */
/** @brief accumulate
 *
 *  Computes the sums of the complex data array, returns sums
 *
 * @param data
 * @param num_samples
 *
 * @return complex_double_t containing accumulation is returned
 */
/* -------------------------------------------------------------------------- */

complex_double_t accumulate(complex_double_t *data, unsigned int num_samples);

/* -------------------------------------------------------------------------- */
/** @brief loop_init
 *
 *  Initializes a few constants, etc for the carrier PLL
 * @param carrier
 * @param integration_time
 *
 * @return
 */
/* -------------------------------------------------------------------------- */
int loop_init(loop_t *params, double integration_time, double BW_phase,
              double BW_freq);

/* -------------------------------------------------------------------------- */
/** @brief pll_run
 *
 *  Given accumulations and the costas_t structure with parameters, this
 *  function returns the frequency adjustment necessary for the next
 *  integration
 *
 * @param params
 *
 * @return Frequency adjustment (double)
 */
/* -------------------------------------------------------------------------- */
double loop_run3(loop_t *params);

/* -------------------------------------------------------------------------- */

// Obtained from:
// http://www.strchr.com/standard_deviation_in_one_pass
/* -------------------------------------------------------------------------- */
/** @brief compute_snr Given a amplitude value and a window of numsamp
 * amplitude values, it first computes the scatter of amplitude in the window
 * and then divides the ampl value by the scatter.
 *
 * @param ampl The S part in Signal/Noise
 * @param ampl_window The scatter of this becomes the N part in S/N
 * @param numsamp Number of samples to use when computing scatter
 * @param *noise The value of the noise floor will be written to the ptr
 *
 * @return Signal to Noise ratio (ampl/std(ampl_window))
 */
/* -------------------------------------------------------------------------- */
double compute_snr(double ampl, double ampl_window[], int numsamp,
                   double *noise);

int prnMultiplyAccumulate(corr_t *corr, cr_t *pcr, uint8_t *pprn,
                          acq_parameters_t *acq);

int serial_acquire(complex_double_t *pbuf_data, uint8_t *pprn, int numsamples,
                   double initial_doppler_hz, float ss_doppler_span_hz,
                   acq_parameters_t *acq);

int fft2_acquire(complex_double_t *pbuf_data, uint8_t *pprn, int numsamples,
                 double initial_doppler_hz, float ss_doppler_span_hz,
                 acq_parameters_t *acq);

// Returns the sign (+1 or -1)
// of a double number
double dsign(double number);

// Given array of timetags and values, does a lsq fit
double fitter(double *timetag_s, double *value, int nelements, int order,
              double *chisq, double *coef, double eval_at_time, int verbose);

// Finds first tone in data
int fft_acquire_tone(complex_double_t *pbuf_data, int numsamples,
                     acq_parameters_t *acq);

// Given I/Q vectors, finds statistics
// rf_frequency is nominal frequency of carrier in cyc/second (Hz)
stats_signal *computeStatistics(st_accum_tone *d, double rf_frequency,
                                int antenna);

// Removes polynomial from data
stdvec *removePoly(double *time_s, double *data, uint32_t num_elements,
                   double *coef, int order);

#endif
