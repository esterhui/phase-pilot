/* tone_process.c
 *
 *
 * $Id: tone_process.c 11 2011-08-28 22:52:11Z esterhui $
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
  fprintf(stderr, "%s [-a -n] < inputfile > outputfile\n\n", progname);
  fprintf(stderr, "Reads IF data formatted with prsr_parse from stdin and "
                  "writes counter-rotated data to stdout\n\n");
  fprintf(stderr, "OPTIONS\n");
  fprintf(stderr, "\t-a acqfile    :   File containing initial parameters\n");
  fprintf(stderr, "\t-n numsamp    :   Samples to integrate over\n");

  return;
}

int main(int argc, char **argv) {
  size_t samples_to_use = 0, samples_processed = 0;
  char *acqfile = NULL;
  int numsamp = 1000;

  complex_double_t *pdata = NULL;
  data_header_t *phdr = NULL;

  cr_t *pcr = NULL;
  cr_t *pcr_out = NULL;

  pcr = (cr_t *)malloc(sizeof(cr_t));

  /* --- Scan input arguments --- */
  int c;

  opterr = 0;

  while ((c = getopt(argc, argv, "ha:n:")) != -1) {
    switch (c) {
    case 'a':
      acqfile = optarg;
      break;
    case 'n':
      sscanf(optarg, "%d", &numsamp);
      break;
    case 'h':
      displayUsage(argv[0]);
      exit(0);
    default:
      displayUsage(argv[0]);
      exit(0);
    }
  }

  samples_to_use = numsamp;
  pcr->num_samples = numsamp;
  pcr->fs = 3.2e6;
  pcr->fc = 308e3;
  pcr->start_phase = 0;

  // ------- Loop through all data and run acq. on blocks
  while (1) {
    // Read one full frame of data
    phdr = parseDataHeader(stdin);
    if (phdr == NULL)
      break;
    pdata = parseDataPayload(stdin, phdr);
    if (pdata == NULL)
      break;

    assert(phdr->data_size_samples > samples_to_use);

    samples_processed = 0;
    while (samples_processed < phdr->data_size_samples) {
      // Advanced the data pointer
      pcr->data = pdata + samples_processed;

      // Do the counter rotation
      pcr_out = counterRotate(pcr);
      assert(pcr_out != NULL);

      // Write data to stdout as doubles
      if (fwrite(pcr_out->data, sizeof(*pcr_out->data), pcr_out->num_samples,
                 stdout) != pcr_out->num_samples) {
        perror("Error writing samples");
      }

      // Move the start phase over for next batch
      pcr->start_phase = pcr_out->phase_next;
      samples_processed += samples_to_use;

      free(pcr_out->data);
      free(pcr_out);
    }

    if (!phdr)
      free(phdr);
    if (!pdata)
      free(pdata);
  }

  free(pcr);

  return 0;
}
