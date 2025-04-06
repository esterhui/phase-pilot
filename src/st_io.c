#include "st_io.h"

#ifdef COMPUTE_STATISTICS
using namespace arma;
#endif
using namespace std;

uint8_t *read_prncode(char *filename, acq_parameters_t *acq, int reverse) {
  size_t rc = 1;
  char line[160];
  char temp[160];
  uint8_t *buf = NULL, *buf2;
  unsigned int i;

  FILE *fid;
  fid = fopen(filename, "r");

  if (fid == NULL) {
    perror(NULL);
    return NULL;
  }

  // First read the header --------------

  // Get the #-----
  fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "#-%s-\n", temp);
  // Get PRN
  fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "PRN%*[ ]:%*[ ]%d", &acq->prn);
  // Get chipping rate
  fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "Chipping Rate%*[ ]:%*[ ]%lf", &acq->coderate);
  // Get carrier frequency
  fgets(line, sizeof(line), fid);
  rc &=
      sscanf(line, "Carrier Frequency%*[ ]:%*[ ]%lf", &acq->carrier_frequency);
  // Get intermediate frequency
  // fgets(line, sizeof(line), fid);
  // rc&=sscanf(line,"Intermediate Frequency%*[ ]:%*[ ]%lf",
  //&acq->intermediate_frequency);
  // Get codes per databit
  fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "Codes per databit%*[ ]:%*[ ]%d", &acq->codes_per_databit);
  // Get Length
  fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "Length%*[ ]:%*[ ]%d", &acq->codelength);
  // Get Description
  fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "Description%*[ ]:%*[ ]%s", acq->description);
  // Get the #-----
  fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "#-%s-\n", temp);

  if (rc == 0) {
    fprintf(stderr, "Couldn't parse PRN file header\n");
    return NULL;
  }

  // ---------------------------

  // Init some acq stuff
  calc_prn_params(acq);

  // Allocate space for PRN
  buf = (uint8_t *)malloc(acq->codelength);
  assert(buf != NULL);

  rc = fread(buf, 1, acq->codelength, fid);

  if (rc != acq->codelength) {
    fprintf(
        stderr,
        "read_prncode: Out of data, expected %d got %d, cheers mate <=====\n",
        acq->codelength, (int)rc);
    free(buf);
    return NULL;
  }

  // If we're reading backwards in time, we should reverse the PRN code
  if (reverse) {
    buf2 = (uint8_t *)malloc(acq->codelength);
    for (i = 0; i < acq->codelength; i++) {
      buf2[i] = buf[(acq->codelength - 1) - i];
    }
    free(buf);
    buf = buf2;
  }

  fclose(fid);

  return buf;
}

int calc_prn_params(acq_parameters_t *acq) {
  acq->samples_per_chip = acq->fs / acq->coderate;
  acq->chips_per_sample = acq->coderate / acq->fs;
  acq->samples_per_code_period = (acq->samples_per_chip * acq->codelength);
  acq->bytes_per_code_period = (uint32_t)acq->samples_per_code_period;
  return 1;
}

data_header_t *parseDataHeader(FILE *fid) {
  size_t rc = 1;
  char *frc;
  char line[160];
  char temp[160];

  data_header_t *hdr = (data_header_t *)malloc(sizeof(data_header_t));

  if (fid == NULL) {
    perror(NULL);
    free(hdr);
    return NULL;
  }

  // First read the header --------------

  // Get the #-----
  frc = fgets(line, sizeof(line), fid);

  if (frc == NULL) {
    fprintf(stderr, "No header, done\n");
    fseek(fid, 0L, SEEK_SET);
    free(hdr);
    return NULL;
  }

  rc &= sscanf(line, "#-%s-\n", temp);
  // Get Date
  fgets(line, sizeof(line), fid);
  rc &=
      sscanf(line, "Date%*[ ]:%*[ ]%u %u %u", &hdr->year, &hdr->doy, &hdr->sec);
  // Get offset from utc (in ns)
  fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "UTC Offset%*[ ]:%*[ ]%d [ns]", &hdr->data_time_offset);
  // Get sample rate
  fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "Sample rate%*[ ]:%*[ ]%lf [MHz]", &hdr->fs);
  hdr->fs *= 1e6; /* Convert to Hz */
  // Get RF sky frequency
  fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "RF Frequency%*[ ]:%*[ ]%lf [Hz]", &hdr->rf_freq);
  // Get Length
  fgets(line, sizeof(line), fid);
  rc &=
      sscanf(line, "Number of samples%*[ ]:%*[ ]%ud", &hdr->data_size_samples);
  // Get Bits per sample
  fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "Bits per sample%*[ ]:%*[ ]%d", &hdr->bits_per_sample);
  // Get Number of channels
  fgets(line, sizeof(line), fid);
  int rc2 = sscanf(line, "Channels%*[ ]:%*[ ]%d", &hdr->num_channels);
  if (rc2 == 0) {
    fprintf(stderr, "Channel not found, assuming 1\n");
    hdr->num_channels = 1;
  } else {
    fgets(line, sizeof(line), fid);
  }
  // Get Description
  rc &= sscanf(line, "Description%*[ ]:%*[ ]%s", hdr->description);
  // Get the #-----
  fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "#-%s-\n", temp);
  // Get the #-----

  if (rc == 0) {
    fseek(fid, 0L, SEEK_SET);
    free(hdr);
    return NULL;
  }

  // ---------------------------
  if (hdr->bits_per_sample > 1) {
    // Compute size in bytes, we divide by 8 to convert to bytes and
    // multiply by 2 assuming complex real/imag data.
    hdr->data_size_bytes =
        hdr->bits_per_sample * hdr->data_size_samples / 8 * 2;
  } else if (hdr->bits_per_sample == 1) {
    hdr->data_size_bytes = hdr->data_size_samples;
  } else {
    fprintf(stderr, "bits_per_sample of [%d] not supported\n",
            hdr->bits_per_sample);
  }

  return hdr;
}

complex_double_t *parseDataPayload(FILE *fid, data_header_t *phdr,
                                   unsigned long long seeksamples,
                                   unsigned int numsamples, int chan,
                                   int flipiq) {
  size_t samples_to_read = 0, samples_left = 0, samples_read = 0;
  size_t rc_read = 0, prev_rc_read = 0;
  size_t i;
  double timeout_counter = 0;
  size_t bytes_per_sample = 1;
  int isTRIG_RAW = 0;
  int isPX1500_REAL = 0;
  int isGPS1A_REAL = 0;
  int isGNURADIO_COMPLEX = 0;
  int isGNURADIO_COMPLEX_SHORT = 0;
  int isReal = 0;

  complex_double_t *pdata = NULL;

  /* ---- Our buffers ------ */
  void *pbuf_data = NULL;
  int8_t *pbuf8 = NULL;
  int16_t *pbuf16 = NULL;
  int32_t *pbuf32 = NULL;
  float *pbuf32f = NULL;
  double *pbuf64f = NULL;

  // Check if this is TriG data
  if (strstr(phdr->description, "TRIG_RAW") != NULL) {
    isTRIG_RAW = 1;
  }
  // Check if this is PX1500 real-only 8-bit data
  else if (strstr(phdr->description, "PX1500_8BIT") != NULL) {
    isPX1500_REAL = 1;
    isReal = 1;
  } else if (strstr(phdr->description, "GPS1A_8BIT") != NULL) {
    isGPS1A_REAL = 1;
    isReal = 1;
  } else if (strstr(phdr->description, "GNURADIO_COMPLEX_SHORT") != NULL) {
    isGNURADIO_COMPLEX_SHORT = 1;
  } else if (strstr(phdr->description, "GNURADIO_COMPLEX") != NULL) {
    isGNURADIO_COMPLEX = 1;
  }

  samples_left = numsamples;
  // Sanity check, we only support 1 and 8 bits per sample currently
  // assert((phdr->bits_per_sample==1) || (phdr->bits_per_sample==8));
  assert((phdr->bits_per_sample == 1) || (phdr->bits_per_sample == 8) ||
         (phdr->bits_per_sample == 16) || (phdr->bits_per_sample == 32) ||
         (phdr->bits_per_sample == 64));

  if (isReal) {
    // 1 = IQ, 8= bits_to_bytes
    bytes_per_sample = (1 * phdr->bits_per_sample * phdr->num_channels) / 8;
  } else {
    // 2 = IQ, 8= bits_to_bytes
    bytes_per_sample = (2 * phdr->bits_per_sample * phdr->num_channels) / 8;
  }
  // fprintf(stderr,"bytes per sample is %zd\n",bytes_per_sample);

  // Allocate space for data
  pbuf_data = malloc(numsamples * bytes_per_sample); // 2 byte per I/Q pair

  pbuf8 = (int8_t *)pbuf_data;
  pbuf16 = (int16_t *)pbuf_data;
  pbuf32 = (int32_t *)pbuf_data;
  pbuf32f = (float *)pbuf_data;
  pbuf64f = (double *)pbuf_data;

  // fprintf(stderr,"skipping %zd bytes\n",bytes_per_sample*seeksamples);
  if (skipBytes(fid, bytes_per_sample * seeksamples) < 0) {
    fprintf(stderr, "Seek into file failed");
    goto alldone;
  }

  // Allocate real buffer here, since 1byte of input data
  // is 2 bytes of complex int8 data, do times 2 here on length.
  pdata = (complex_double_t *)malloc(numsamples * sizeof(complex_double_t));
  assert(pdata != NULL);

  // ------ Read the data
  while (samples_left > 0) {
    if (samples_left > DATA_BLOCKSIZE) {
      samples_to_read = DATA_BLOCKSIZE;
    } else {
      samples_to_read = samples_left;
    }

    // fseek(fid,0L,SEEK_SET);
    rc_read = fread(pbuf_data, bytes_per_sample, samples_to_read, fid);

    // fprintf(stderr,"chan[%d] %08x\n",chan,pbuf32[1]);

    if (isGNURADIO_COMPLEX && phdr->bits_per_sample == 32) {
      // fprintf(stderr,"PX1500 data\n");
      for (i = 0; i < rc_read; i++) {
        pdata[samples_read + i].real = 10000 * pbuf32f[2 * i + 0];
        pdata[samples_read + i].imag = 10000 * pbuf32f[2 * i + 1];
      }

    } else if (isGNURADIO_COMPLEX && phdr->bits_per_sample == 64) {
      // SMAP data, double, but convert to complex
      for (i = 0; i < rc_read; i++) {
        pdata[samples_read + i].real = (float)(10000.0 * pbuf64f[2 * i + 0]);
        pdata[samples_read + i].imag = (float)(10000.0 * pbuf64f[2 * i + 1]);
      }
    } else if (isGNURADIO_COMPLEX_SHORT && phdr->bits_per_sample == 16) {
      for (i = 0; i < rc_read; i++) {
        pdata[samples_read + i].real = (float)pbuf16[2 * i + 0];
        pdata[samples_read + i].imag = (float)pbuf16[2 * i + 1];
      }
    } else if (isPX1500_REAL && phdr->bits_per_sample == 8) {
      // fprintf(stderr,"PX1500 data\n");
      for (i = 0; i < rc_read; i++) {
        pdata[samples_read + i].real =
            (double)(int8_t)((uint8_t)pbuf8[i] - 128);
        pdata[samples_read + i].imag = 0;
      }
    } else if (isGPS1A_REAL && phdr->bits_per_sample == 8) {
      // fprintf(stderr,"GPS1A data\n");
      for (i = 0; i < rc_read; i++) {
        pdata[samples_read + i].real = (double)pbuf8[i];
        pdata[samples_read + i].imag = 0;
      }
    } else if (phdr->bits_per_sample == 8) {
      for (i = 0; i < rc_read; i++) {
        pdata[samples_read + i].real = (double)(int8_t)(pbuf16[i] & 0xff);
        pdata[samples_read + i].imag =
            (double)(int8_t)((pbuf16[i] >> 8) & 0xff);
      }
    } else if (phdr->bits_per_sample == 1) {
      for (i = 0; i < rc_read; i++) {
        if (flipiq == 1) { // Invert IQ
          if (phdr->num_channels == 16) {
            // fprintf(stderr,"%d\n",((pbuf32[i]>>(2*chan))&0x1)==0x1?+1:-1);
            pdata[samples_read + i].imag =
                ((pbuf32[i] >> (2 * chan)) & 0x1) == 0x1 ? +1 : -1;
            pdata[samples_read + i].real =
                ((pbuf32[i] >> (2 * chan)) & 0x2) == 0x2 ? +1 : -1;
          } else if (phdr->num_channels == 8) {
            // fprintf(stderr,"%d\n",((pbuf32[i]>>(2*chan))&0x1)==0x1?+1:-1);
            pdata[samples_read + i].imag =
                ((pbuf16[i] >> (2 * chan)) & 0x1) == 0x1 ? +1 : -1;
            pdata[samples_read + i].real =
                ((pbuf16[i] >> (2 * chan)) & 0x2) == 0x2 ? +1 : -1;
          } else if (isTRIG_RAW) {
            pdata[samples_read + i].imag =
                ((pbuf8[i] >> (chan)) & 0x1) == 0x1 ? +1 : -1;
            pdata[samples_read + i].real =
                ((pbuf8[i] >> (4 + chan)) & 0x1) == 0x1 ? +1 : -1;
          } else {
            pdata[samples_read + i].imag =
                ((pbuf8[i] >> (2 * chan)) & 0x1) == 0x1 ? +1 : -1;
            pdata[samples_read + i].real =
                ((pbuf8[i] >> (2 * chan)) & 0x2) == 0x2 ? +1 : -1;
          }
        } else if (flipiq == 0) {
          if (phdr->num_channels == 16) {
            pdata[samples_read + i].real =
                ((pbuf32[i] >> (2 * chan)) & 0x1) == 0x1 ? +1 : -1;
            pdata[samples_read + i].imag =
                ((pbuf32[i] >> (2 * chan)) & 0x2) == 0x2 ? +1 : -1;
          } else if (phdr->num_channels == 8) {
            // fprintf(stderr,"%d\n",((pbuf32[i]>>(2*chan))&0x1)==0x1?+1:-1);
            pdata[samples_read + i].real =
                ((pbuf16[i] >> (2 * chan)) & 0x1) == 0x1 ? +1 : -1;
            pdata[samples_read + i].imag =
                ((pbuf16[i] >> (2 * chan)) & 0x2) == 0x2 ? +1 : -1;
          } else if (isTRIG_RAW) {
            pdata[samples_read + i].real =
                ((pbuf8[i] >> (chan)) & 0x1) == 0x1 ? +1 : -1;
            pdata[samples_read + i].imag =
                ((pbuf8[i] >> (4 + chan)) & 0x1) == 0x1 ? +1 : -1;
          } else {
            pdata[samples_read + i].real =
                ((pbuf8[i] >> (2 * chan)) & 0x1) == 0x1 ? +1 : -1;
            pdata[samples_read + i].imag =
                ((pbuf8[i] >> (2 * chan)) & 0x2) == 0x2 ? +1 : -1;
          }
        }
        if (flipiq == 2) { // Trig SCIP format (I,I,I,I,Q,Q,Q,Q)
          if (phdr->num_channels == 4) {
            pdata[samples_read + i].real =
                ((pbuf8[i] >> (chan)) & 0x1) == 0x1 ? +1 : -1;
            pdata[samples_read + i].imag =
                ((pbuf8[i] >> (4 + chan)) & 0x1) == 0x1 ? +1 : -1;
          } else {
            fprintf(stderr, "Only 4 channels supported for SciP format\n");
            exit(1);
          }
        }
      }
    }

    samples_read += rc_read;
    samples_left -= rc_read;

    if (rc_read < 1) {
      // Should do some type of select() here, but when I do,
      // it just immediately returns for pipes and here files
      usleep(NO_DATA_SLEEP * 1e6);

      if (prev_rc_read < 1)
        timeout_counter += (NO_DATA_SLEEP);
      if (timeout_counter > READ_TIMEOUT) {
        fprintf(stderr, "INFO: No more data, exiting\n");
        goto alldone;
      }
    }

    // Save off old rc_read, we use this to increment
    // a timeout counter.
    prev_rc_read = rc_read;
  }

  // Done with temp buf
  free(pbuf_data);

  return pdata;

alldone:
  free(pdata);
  free(pbuf_data);
  return NULL;
}

int printDataHeader(FILE *fid, data_header_t *pheader) {
  fprintf(fid, "#---------------------------------------------\n");
  fprintf(fid, "Date              :   %d %d %d\n", pheader->year, pheader->doy,
          pheader->sec);
  fprintf(fid, "UTC Offset        :   %d [ns]\n", pheader->data_time_offset);
  fprintf(fid, "Sample rate       :   %f [MHz]\n", pheader->fs / 1.0e6);
  fprintf(fid, "RF Frequency      :   %lf\n", pheader->rf_freq);
  fprintf(fid, "Number of samples :   %d\n", pheader->data_size_samples);
  fprintf(fid, "Bits per sample   :   %d\n", pheader->bits_per_sample);
  fprintf(fid, "Description       :   %s\n", pheader->description);
  fprintf(fid, "#---------------------------------------------\n");
  fflush(fid);

  return 1;
}

int printAcquisitionHeader(FILE *fid, acq_parameters_t *hdr) {
  fprintf(fid, "#---------------------------------------------\n");
  fprintf(fid, "PRN                :   %d\n", hdr->prn);
  fprintf(fid, "Sample Number      :   %lld\n", hdr->sample_number);
  fprintf(fid, "Code Phase         :   %lf [chips]\n", hdr->code_phase);
  fprintf(fid, "Doppler Frequency  :   %lf [Hz]\n", hdr->doppler);
  fprintf(fid, "Carrier Phase      :   %lf [rad]\n", hdr->carrier_phase);
  fprintf(fid, "CNo                :   %0.1lf [dB-Hz]\n",
          20.0 * log10(hdr->vsnr_1s));
  fprintf(fid, "1sVSNR             :   %0.0f [V/V]\n", hdr->vsnr_1s);
  fprintf(fid, "VSNR               :   %0.1f [V/V]\n", hdr->vsnr);
  fprintf(fid, "Magnitude          :   %0.1f\n", hdr->magnitude);
  fprintf(fid, "Noise              :   %0.1f\n", hdr->noise_floor);
  fprintf(fid, "#---------------------------------------------\n");
  fflush(fid);

  return 1;
}

acq_parameters_t *parseAcquisitionHeader(FILE *fid) {
  size_t rc = 1;
  char *frc;
  char line[160];
  char temp[160];

  acq_parameters_t *hdr = (acq_parameters_t *)malloc(sizeof(acq_parameters_t));

  if (fid == NULL) {
    fprintf(stderr, "Acquisition file not provided, use -a\n");
    free(hdr);
    return NULL;
  }

  // First read the header --------------

  // Get the #-----
  frc = fgets(line, sizeof(line), fid);

  if (frc == NULL) {
    fprintf(stderr, "No header, done\n");
    free(hdr);
    return NULL;
  }

  rc &= sscanf(line, "#-%s-\n", temp);
  // Get  PRN
  frc = fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "PRN%*[ ]:%*[ ]%d", &hdr->prn);
  // Get Sample Number
  frc = fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "Sample Number%*[ ]:%*[ ]%lld", &hdr->sample_number);
  // Get Code Phase
  frc = fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "Code Phase%*[ ]:%*[ ]%lf [chips]", &hdr->code_phase);
  // Get Carrier Freq
  frc = fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "Doppler Frequency%*[ ]:%*[ ]%lf [Hz]", &hdr->doppler);
  // Get Carrier Phase
  frc = fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "Carrier Phase%*[ ]:%*[ ]%lf [rad]", &hdr->carrier_phase);
  // Get CNo
  frc = fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "CNo%*[ ]:%*[ ]%lf [Hz]", &hdr->magnitude);
  // Get 1sVSNR
  frc = fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "1sVSNR%*[ ]:%*[ ]%lf [V/V]", &hdr->vsnr_1s);
  // Get VSNR
  frc = fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "VSNR%*[ ]:%*[ ]%lf [V/V]", &hdr->vsnr);
  // Get Magnitude
  frc = fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "Magnitude%*[ ]:%*[ ]%lf", &hdr->magnitude);
  // Get Noise
  frc = fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "Noise%*[ ]:%*[ ]%lf", &hdr->noise_floor);
  // Get the #-----
  frc = fgets(line, sizeof(line), fid);
  rc &= sscanf(line, "#-%s-\n", temp);

  if (rc == 0) {
    fprintf(stderr, "Couldn't parse Acquisition file header\n");
    free(hdr);
    return NULL;
  }

  return hdr;
}

int skipBytes(FILE *fid, long long int bytes) {
  int rc = 0, prev_rc = 0;
  long long int bytesToRead, bytesLeft;
  float timeout_counter = 0;
  int64_t rc_long;

  // Allocate a temp buffer
  void *tempBuf = malloc(TEMP_BLOCKSIZE);

  bytesLeft = bytes;

  // If no seek necessary, will just return
  if (bytes == 0) {
    free(tempBuf);
    return 1;
  }

  // First try to use fseek, if it fails, read byte-for-byte
  if ((rc_long = lseek64(fileno(fid), bytes, SEEK_CUR)) == bytes) {
    free(tempBuf);
    return 1;
  } else {
    perror("Seek failed : ");
    fprintf(
        stderr,
        "Fseek DID NOT work (rc=%lld), seeking %lld bytes manually ***** \n",
        rc_long, bytes);
  }

  while (bytesLeft > 0) {
    if (bytesLeft >= TEMP_BLOCKSIZE) {
      bytesToRead = TEMP_BLOCKSIZE;
    } else {
      bytesToRead = bytesLeft;
    }
    rc = fread(tempBuf, 1, bytesToRead, fid);
    bytesLeft -= rc;
    if (rc < 1) {
      usleep(NO_DATA_SLEEP * 1e6);
      if (prev_rc < 1)
        timeout_counter += (NO_DATA_SLEEP);
      if (timeout_counter > READ_TIMEOUT) {
        fprintf(stderr, "INFO: No more data, exiting\n");
        break;
      }
    }
    prev_rc = rc;
  }

  free(tempBuf);

  if (bytes == 0)
    rc = 1;
  return rc == 0 ? -1 : 1;
}

int printCorrelWaveform(FILE *fid, int num_lags, double *lag_s,
                        double *lag_magnitude, double *lag_phase,
                        data_header_t *datahdr, acq_parameters_t *acqhdr) {
  int i;

  fprintf(fid, "# $Revision: 124 $, Data Date:%d %d %d,", datahdr->year,
          datahdr->doy, datahdr->sec);
  fprintf(fid, "UTC Offset:%d ns,", datahdr->data_time_offset);
  fprintf(fid, "Fs: %f Hz,", datahdr->fs);
  fprintf(fid, "Samples skipped:%lld\n", acqhdr->sample_number);

  fprintf(fid, "# Lag:s, Magnitude, Phase:ns\n");

  for (i = 0; i < num_lags; i++) {
    fprintf(fid, "% 9.9f % 9.1f % 14.12f\n", lag_s[i], lag_magnitude[i],
            lag_phase[i]);
  }

  return 1;
}

void printAccumHeader(FILE *fid, data_header_t *datahdr,
                      acq_parameters_t *acqhdr, double integ_periods,
                      double lag_spacing_chips) {
  fprintf(fid, "# $Revision: 124 $, Data Date:%d %d %d,", datahdr->year,
          datahdr->doy, datahdr->sec);
  fprintf(fid, "UTC Offset:%d ns,", datahdr->data_time_offset);
  fprintf(fid, "Fs: %f Hz,", datahdr->fs);
  fprintf(fid, "Samples skipped:%lld,", acqhdr->sample_number);
  fprintf(fid, "Accumultions per databit:%lf,",
          acqhdr->codes_per_databit / integ_periods);
  // Here we multiply by 2 because lag_spacing_chips is the
  // spacing between E,P,L, thus E-to-L is 2x that
  fprintf(fid, "E-L Spacing: %0.3f Chips\n", lag_spacing_chips * 2.0);
}

#include <ctype.h>

int makeargv(char *string, char *argv[], int argvsize) {
  char *p = string;
  int i;
  int argc = 0;

  for (i = 0; i < argvsize; i++) {
    /* skip leading whitespace */
    while (isspace(*p))
      p++;

    if (*p != '\0')
      argv[argc++] = p;
    else {
      argv[argc] = 0;
      break;
    }

    /* scan over arg */
    while (*p != '\0' && !isspace(*p))
      p++;
    /* terminate arg: */
    if (*p != '\0' && i < argvsize - 1)
      *p++ = '\0';
  }

  return argc;
}

#ifdef COMPUTE_STATISTICS
mat *parse_st_data(const char *filename, double *fs, int skip) {
  mat *D = new mat;
  string line1, line2;
  int rc;

  ifstream stfile;
  stfile.open(filename, ios::in);
  if (!stfile.is_open()) {
    cout << "Unable to open file : " << filename << endl;
    delete D;
    return NULL;
  }

  // Read the two header lines
  getline(stfile, line1);
  getline(stfile, line2);

  // Find the FS
  rc = sscanf(line1.c_str(), "%*[^,], %*[^,], %*[^,], Fs: %lf", fs);
  if (rc != 1) {
    /*
    cout << "Error, couldn't load sample rate from line 1 of file"
        << filename  << endl;
        */
    stfile.close();
    delete D;
    return NULL;
  }

  bool status = D->load(stfile);

  if (status == false) {
    // cout << "Error, couldn't load ascii data from " << filename << endl;
    stfile.close();
    delete D;
    return NULL;
  }

  size_t rows = D->n_rows;
  size_t cols = D->n_cols;

  // Get new matrix where we cut off the first 'skip
  // elements
  // D=D(skip:end,0:end)
  *D = D->submat(span(skip, rows - 1), span(0, cols - 1));

  stfile.close();

  return D;
}

st_accum_tone *parse_st_accum_tone(const char *filename, int skip) {
  double fs;
  st_accum_tone *d = new st_accum_tone;

  // Read the thing
  mat *D = parse_st_data(filename, &fs, skip);

  if (D == NULL) {
    delete d;
    return NULL;
  }

  /* break matrix up into stuff we want */
  d->timetag_s = D->col(I_TIMETAG);
  d->ip = D->col(I_IP);
  d->qp = D->col(I_QP);
  d->doppler_hz = D->col(I_DOP_HZ);
  d->phase_s = D->col(I_PHASE_SEC);
  d->nsamp = D->col(I_NSAMP);
  d->fs = fs;

  return d;
}

void printStatistics(FILE *fid, std::vector<stats_signal *> stats) {
  int numstats = stats.size();
  fprintf(fid, "---------------------------------------------------------------"
               "----------\n");
  fprintf(fid, "|Ant, |  CN0  | Ampl*|  Ampl  | Phase  |ampl/phase| Velocity | "
               "Phase*   |\n");
  fprintf(fid, "|chan |(dB-Hz)| scat | 1sVSNR | 1sVSNR | SNR (dB) |  (ns/s)  | "
               "scat (ps)|\n");
  fprintf(fid, "|--------------------------------------------------------------"
               "---------|\n");
  for (int i = 0; i < numstats; i++) {
    if (stats[i]->cn0 >= 0)
      fprintf(fid,
              "|%2d,%2d|%7.1lf|%6.0lf|%8.0lf|%8.0lf|%10.1lf|%10.3lf|%10.1lf|\n",
              stats[i]->antenna_number, (i % 4), stats[i]->cn0,
              stats[i]->rss_std, stats[i]->vsnr1s, stats[i]->phase_vsnr1s,
              20 * log10(stats[i]->vsnr1s / stats[i]->phase_vsnr1s),
              stats[i]->phase_rate, stats[i]->phase_std_ns * 1e3);
    else
      fprintf(fid,
              "|%2d,%2d|  .    |   .  |    .   |   .    |      .   |     .    "
              "|      .   |\n",
              stats[i]->antenna_number, (i % 4));
    if ((i + 1) % 4 == 0) {
      fprintf(fid, "|----------------------------------------------------------"
                   "-------------|\n");
    }
  }
}
#endif
