#ifndef __ST_IO_H__
#define __ST_IO_H__

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <unistd.h>

#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include <sys/time.h>
#include <sys/types.h>

#include <arpa/inet.h>

#include "st_datatypes.h"

#define DATA_BLOCKSIZE 524288

#define TEMP_BLOCKSIZE 524288

/* States for reading header & data */
#define READ_HEADER 0
#define READ_DATA 1

#define READ_TIMEOUT                                                           \
  (0.1) /* Timeout after this many seconds when no samples                     \
         show up */
#define NO_DATA_SLEEP                                                          \
  (0.02) /* Sleep this long if no data available before                        \
            trying again */

/* -------------------------------------------------------------------------- */
/** @brief read_prncode
 *
 * @param filename
 * @param acq
 *
 * @return
 */
/* -------------------------------------------------------------------------- */
uint8_t *read_prncode(char *filename, acq_parameters_t *acq, int reverse);

int calc_prn_params(acq_parameters_t *acq);

/* -------------------------------------------------------------------------- */
/** @brief prsr_parse
 *
 * Parses 1 frame of PRSR data (usually 1 second of data)
 *
 * @param fd File descriptor (file or socket)
 * @param chan_header
 * @param chan_data
 *
 * @return
 */
/* -------------------------------------------------------------------------- */
chan_t *prsr_parse(int fd);

data_header_t *parseDataHeader(FILE *fid);

complex_double_t *parseDataPayload(FILE *fid, data_header_t *phdr,
                                   unsigned long long seeksamples,
                                   unsigned int numsamples, int chan,
                                   int flipiq);

/* -------------------------------------------------------------------------- */
/** @brief printDataHeader
 *
 * @param pheader
 * @param fid
 *
 * @return
 */
/* -------------------------------------------------------------------------- */
int printDataHeader(FILE *fid, data_header_t *pheader);

/* -------------------------------------------------------------------------- */
/** @brief printAcquisitionHeader
 *
 * @param fid
 * @param hpre
 * @param hpost
 *
 * @return
 */
/* -------------------------------------------------------------------------- */
int printAcquisitionHeader(FILE *fid, acq_parameters_t *hdr);

/* -------------------------------------------------------------------------- */
/** @brief parseAcquisitionHeader
 *
 * @param fid
 *
 * @return
 */
/* -------------------------------------------------------------------------- */
acq_parameters_t *parseAcquisitionHeader(FILE *fid);

/* -------------------------------------------------------------------------- */
/** @brief prnMultiplyAccumulate
 *
 *  Multiplies data in pcr with data in pprn according to corr parameters.
 *  and accumulates
 *
 * @param corr
 * @param pcr
 * @param pprn
 * @param acq
 *
 * @return
 */
/* -------------------------------------------------------------------------- */
int prnMultiplyAccumulate(corr_t *corr, cr_t *pcr, uint8_t *pprn,
                          acq_parameters_t *acq);

/* -------------------------------------------------------------------------- */
/** @brief skipBytes
 *
 * Skip this many bytes in a stream. We use this instead of fseek, since fseek
 * will fail on a PIPE
 *
 * @param fid
 * @param bytes
 *
 * @return 0 on error, 1 on success
 */
/* -------------------------------------------------------------------------- */
int skipBytes(FILE *fid, long long int bytes);

/* -------------------------------------------------------------------------- */
/** @brief input_timeout
 *
 * Blocks on file descriptor filedes for seconds or until activity, this should
 * be used when doing reads for data, but it looks like it immediately
 * returns when doing reads from pipes and here files. I'll leave it in here
 * for now.
 *
 * Obtained from :
 *  http://www.gnu.org/s/hello/manual/libc/Waiting-for-I_002fO.html
 *
 * @param filedes File descriptor fileno(FILE*)
 * @param seconds Will time out after this many seconds
 *
 * @return 0 on timeout, -1 on error, 1 on file descriptor ready
 */
/* -------------------------------------------------------------------------- */
int input_timeout(int filedes, unsigned int seconds);

/* -------------------------------------------------------------------------- */
/** @brief printCorrelWaveform
 *
 * @param fid
 * @param num_lags
 * @param lag_s
 * @param lag_phase
 * @param lag_magnitude
 * @param datahdr
 * @param acqhdr
 *
 * @return
 */
/* -------------------------------------------------------------------------- */
int printCorrelWaveform(FILE *fid, int num_lags, double *lag_s,
                        double *lag_magnitude, double *lag_phase,
                        data_header_t *datahdr, acq_parameters_t *acqhdr);

void printAccumHeader(FILE *fid, data_header_t *datahdr,
                      acq_parameters_t *acqhdr, double integ_periods,
                      double lag_spacing_chips);

int makeargv(char *string, char *argv[], int argvsize);

void printStatistics(FILE *fid, std::vector<stats_signal *> stats);

/**
 * Reads signal tracker data, (raw, tone, or PR)
 */
mat *parse_st_data(const char *filename, double *fs);

/**
 * Reads signal tracker tone data
 */
st_accum_tone *parse_st_accum_tone(const char *filename, int skip);

#define I_TIMETAG (0)
#define I_IP (1)
#define I_QP (2)
#define I_DOP_HZ (3)
#define I_PHASE_SEC (4)
#define I_NSAMP (5)

void printStatistics(FILE *fid, std::vector<stats_signal *> stats);

#endif
