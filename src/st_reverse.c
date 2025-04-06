/* st_reverse
 *
 *  Reverses a ST data file. IE, if ST_8BIT_IQ, byte1,2 will become
 *  byte 99,byte100. This is used for playing files backwards.
 *  The timetag is also properly modifed
 *
 * $Id: st_reverse.c 89 2012-06-06 20:56:26Z esterhui $
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

void displayUsage(char *progname) {
  fprintf(stderr, "%s [-s -n] < inputfile > outputfile\n\n", progname);
  fprintf(stderr, "Reads ST data from stdin and writes I/Q accumulations in "
                  "reverse order to stdout\n\n");
  fprintf(stderr, "$Revision: 75 $\n\n");
  fprintf(stderr, "OPTIONS\n");
  fprintf(stderr, "\t-s seektime    :  Seek this many seconds into file\n");
  fprintf(stderr, "\t-n seconds        Save this many seconds of data\n");
  fprintf(stderr, "\t                  Default is all data, or -1\n");

  return;
}

int main(int argc, char **argv) {

  int i, j;
  data_header_t *phdr = NULL;
  void *pdata = NULL;
  FILE *fid = stdin;
  int8_t *pbuf8 = NULL;

  double seek_seconds = 0, num_seconds = -1;
  long seek_sample_number, first_data_sample, seek_bytes;
  unsigned long samples_read = 0;
  size_t rc_read, rc;
  size_t one_sec_samples = 0;
  int bytes_per_sample_pair = 1;

  /* --- Scan input arguments --- */
  int c;

  opterr = 0;

  while ((c = getopt(argc, argv, "hn:s:")) != -1) {
    switch (c) {
    case 'n':
      sscanf(optarg, "%lf", &num_seconds);
      break;
    case 's':
      sscanf(optarg, "%lf", &seek_seconds);
      seek_seconds = round(seek_seconds);
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
  phdr = parseDataHeader(fid);
  if (phdr == NULL)
    return -1;

  one_sec_samples = phdr->fs;
  first_data_sample = ftell(fid);

  // Allocate space for data
  if (phdr->bits_per_sample == 1) {
    bytes_per_sample_pair = 1;
    pdata = malloc(one_sec_samples * bytes_per_sample_pair);
  } else if (phdr->bits_per_sample == 8) {
    bytes_per_sample_pair = 2;
    pdata = malloc(one_sec_samples * bytes_per_sample_pair);
  }

  pbuf8 = (int8_t *)pdata;

  fprintf(stderr, "First sample starts at [%lf]\n", seek_seconds * phdr->fs);

  phdr->data_size_samples = num_seconds * phdr->fs;
  phdr->sec += seek_seconds;
  phdr->rf_freq *= -1.0;

  printDataHeader(stdout, phdr);

  while (1) {

    // Here we add +1 so the first sample is truly on the
    // integer second, otherwise the first sample will
    // be integer_second - 1/fs;
    seek_sample_number = (seek_seconds * phdr->fs - samples_read) + 1;

    seek_bytes =
        ((seek_sample_number - one_sec_samples) * bytes_per_sample_pair);
    rc = fseek(fid, seek_bytes + first_data_sample, SEEK_SET);
    if (rc == -1) {
      perror("Seek failed");
      break;
    }

    rc_read = fread(pdata, bytes_per_sample_pair, one_sec_samples, fid);

    if (rc_read != one_sec_samples) {
      fprintf(stderr, "rc_read=%ld, requested=%ld\n", rc_read, one_sec_samples);
      perror("Read failed");
      break;
    }

    samples_read += rc_read;

    // Now loop backwards through samples and write to stdout
    for (i = 0; i < one_sec_samples; i++) {
      j = ((one_sec_samples - 1) * bytes_per_sample_pair) -
          (i * bytes_per_sample_pair);
      rc = fwrite(pbuf8 + j, 1, bytes_per_sample_pair, stdout);
      if (rc != bytes_per_sample_pair) {
        perror("fwrite failed");
        break;
      }
    }

    fprintf(stderr, "samples read is %lf/%lf\n", samples_read / phdr->fs,
            num_seconds);
    if (samples_read / phdr->fs >= num_seconds) {
      break;
    }
  }

  fclose(fid);
  free(pdata);
  free(phdr);

  return 0;
}
