/* fft_acquire.c
 *
 * This program reads I/Q data and performs acquisition
 * on PRN-based sequences. This software uses a hybrid FFT acquisition method,
 * here is the scoop:
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
 * $Id: fft_acquire.c 125 2013-10-16 01:10:55Z esterhui $
 * 
 */ 
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include <assert.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <fftw3.h>

#include "st_datatypes.h"
#include "st_dsp.h"
#include "st_io.h"

#define FINE_BINS_SS (20)  // Single sided fine bins to search
#define HDR_DIR ("/home/esterhui/.signalTracker/")

// Only good to SNRs of about 200 V/V with 10ms integration times

void displayUsage(char *progname) 
{
    fprintf(stderr,"%s [-p -t] < inputfile > outputfile\n\n",progname);
    fprintf(stderr,"Reads PRSR data formatted with prsr_parse from stdin ");
    fprintf(stderr,"and writes acquisition parameters to stdout\n\n");
    fprintf(stderr,"OPTIONS\n");
    fprintf(stderr,"\t-p prnfile    :   File containing PRN to correlate against\n");
    fprintf(stderr,"\t-e            :   Enhanced doppler search (fine doppler steps)\n");
    fprintf(stderr,"\t-t intnum     :   Number of code periods to use in acquisition,\n");
    fprintf(stderr,"\t                  eg 1,2,3, etc. 1 is default\n");
    fprintf(stderr,"\t-s seektime   :   Seek this many seconds into data\n");
    fprintf(stderr,"\t                  before starting acq. Default=0.0 \n");
    fprintf(stderr,"\t-c correl.txt :   Writes correlation waveform to this file\n");
    fprintf(stderr,"\t-m channel    :   If multiple channels in data, use channel m, default m=0\n");
    fprintf(stderr,"\t-i doppler    :   Initial estimate of doppler in Hz\n");
    fprintf(stderr,"\t-r dop_range  :   Will search +- this range in Hz (+-35000 Hz default)\n");
    fprintf(stderr,"\t-y acq_type   :   Type of acquisition to use:\n");
    fprintf(stderr,"\t                  0 - Hybrid FFT Acquire\n");
    fprintf(stderr,"\t                  1 - Serial Acquire\n");
    fprintf(stderr,"\t                  2 - FFT/IFFT Acquire\n");
    fprintf(stderr,"\t-v vsnr_thres :   VSNR detection threshold, default=6.0:\n");
    fprintf(stderr,"\t-a            :   Only print acquired signals\n");
    fprintf(stderr,"\t-q            :   Flip IQ\n");
    fprintf(stderr,"\t-w            :   Repeated acquisition until found (off by default)\n");

    return;
}

int main(int argc, char **argv) 
{
    size_t samples_to_use=0;
    double integration_periods=1;
    acq_parameters_t acq, acq_old, acq_max;
    char *prnfile=NULL;
    char *correlfile=NULL;
    FILE *fid_correl=NULL;
    FILE *fid_datahdr=NULL;
    int chan_to_use=0;
    double initial_doppler_hz=0;
    int read_backwards=0;
    int dop=0;
    acq_old.vsnr=0;
    acq_old.doppler=0;
    int repeated_acquisition=0;

    double vsnr_threshold=6.0;
    double doppler_range_ss_hz=8000;

    unsigned long long seek_samples=0, num_acquisitions=0, samples_processed=0;
    unsigned long long samples_processed_old=0;
    double seektime=0;
    int fft_type=2;
    int use_old_samples_processed=1;

    complex_double_t* pdata=NULL;
    complex_double_t* pdata_old=NULL;
    data_header_t* phdr=NULL;

    /* ---- Our buffers ------ */
    uint8_t *pprn=NULL;
    int8_t *pprn_sampled=NULL;

    /* --- Scan input arguments --- */
    int c;
    int i;
     
    opterr = 0;
    int finedoppler=0;
    int flipiq=0;
    int printAcquired=0;

    memset(&acq, 0, sizeof(acq));
    memset(&acq_old, 0, sizeof(acq_old));
    memset(&acq_max, 0, sizeof(acq_max));

    while ((c = getopt (argc, argv, "waheqp:t:s:f:c:m:y:i:v:r:")) != -1) {
        switch (c) {
            case 'a':
                printAcquired=1;
                break;
            case 'p':
                prnfile=optarg;
                break;
            case 'e':
                finedoppler=1;
                break;
            case 'w':
                repeated_acquisition=1;
                break;
            case 'q':
		flipiq=1;
                break;
            case 'c':
                correlfile=optarg;
                break;
            case 't':
                sscanf(optarg,"%lf",&integration_periods);
                break;
            case 's':
                sscanf(optarg,"%lf",&seektime);
                break;
            case 'i':
                sscanf(optarg,"%lf",&initial_doppler_hz);
                break;
            case 'm':
                sscanf(optarg,"%d",&chan_to_use);
                break;
            case 'v':
                sscanf(optarg,"%lf",&vsnr_threshold);
                break;
            case 'y':
                sscanf(optarg,"%d",&fft_type);
                break;
            case 'r':
                sscanf(optarg,"%lf",&doppler_range_ss_hz);
                break;
            case 'h':
                displayUsage(argv[0]);
                exit(0);
            default:
                displayUsage(argv[0]);
                exit(0);
        }
    }



    // Read the first header to get some info
    fprintf(stderr,"Trying alternate header file\n");
    char hdr_file[512];
    sprintf(hdr_file,"%s/data_hdr_chan%02d.txt",HDR_DIR,chan_to_use);
    //sprintf(hdr_file,"%s/data_hdr_chan%02d.txt",HDR_DIR,2);
    if ((fid_datahdr=fopen(hdr_file,"r"))==NULL) {
        perror("Can't open alternate header file");
    }
    phdr=parseDataHeader(fid_datahdr);

    
    assert(phdr!=NULL);
    // Fill in sample rate for the acquisition stuff
    acq.fs=phdr->fs;
    acq.vsnr_threshold = vsnr_threshold;


    if (phdr->rf_freq<0) {
        read_backwards=1;
    }
    else read_backwards=0;

    // Read PRN code once
    if ((pprn=read_prncode(prnfile,&acq,read_backwards))==NULL) {
        fprintf(stderr,"Couldn't read PRN file: %s\n",prnfile);
        goto alldone;
    }

    assert(pprn!=NULL);

    if (correlfile!=NULL) {
        fid_correl=fopen(correlfile,"wb");
        if (fid_correl==NULL) {
            perror(NULL);
            goto alldone;
        }
    }

        // FIXME: Hack for old PRSR files with incorrect RF_FREQ setting
    if ((phdr->rf_freq<1e9 || phdr->rf_freq>3e9) && phdr->rf_freq > 0) {
        phdr->rf_freq=2032e6;
        fprintf(stderr,"WARNING, changing rf_freq to %lf Hz\n",
            phdr->rf_freq);
    }



    if (phdr->rf_freq<0) {
        acq.intermediate_frequency=-(phdr->rf_freq+acq.carrier_frequency);
    }
    else {
        acq.intermediate_frequency=(acq.carrier_frequency-phdr->rf_freq);
    }
    fprintf(stderr,"Nominal IF is %lf Hz\n", acq.intermediate_frequency);
    fprintf(stderr,"Search range is +- %lf kHz\n", doppler_range_ss_hz/1000.0);

    // Allocate enough for one code period at rate FS
    samples_to_use=(size_t)round(acq.samples_per_code_period*integration_periods);

    // Sample PRN code at Fs, produce twice the length so 
    // when we start indexing into the sampled PN code, we
    // can always grab a whole code_length
    pprn_sampled=resample_prncode(pprn, acq.codelength, 2*samples_to_use,
        acq.fs, acq.coderate);
    assert(pprn_sampled!=NULL);


    // If user asked to seek into data file, correctly initialize things here
    seek_samples = seektime*acq.fs;
    samples_processed = seek_samples;
    //fprintf(stderr,"samples processed is %lld\n",samples_processed);

    // ------- Loop through all data and run acq. on blocks, we always
    // For now if we don't find any data, seek 1 second into the file
    // and try again.
    while((pdata=parseDataPayload(stdin,phdr,seek_samples,
                samples_to_use,chan_to_use,flipiq))!=NULL) {

        switch (fft_type) {
            case 0 :
                // Hybrid FFT method
                fft_acquire(pdata, pprn_sampled,samples_to_use,&acq);
                break;
            case 1 :
                // Time-domain method
                serial_acquire(pdata, pprn,samples_to_use,initial_doppler_hz,
                    doppler_range_ss_hz,&acq);
                break;
            case 2 :
                // FFT/IFFT method
                fft2_acquire(pdata, pprn,samples_to_use,initial_doppler_hz,
                    doppler_range_ss_hz,&acq);
                break;
        }

        num_acquisitions++;

        // If we're integrating over exactly a bit, 
        // skip by 0.5 integration period to ensure
        // the next integration might not have a bit flip
        if (integration_periods==acq.codes_per_databit) {
            seek_samples=(int)round(samples_to_use*0.5);
        }
        else {
            // Otherwise don't skip any samples, just
            // continue where we left off
            seek_samples=0;
        }

        if (num_acquisitions>1) {
            samples_processed_old=samples_processed;
            samples_processed+=samples_to_use+seek_samples;
        }
        // Do one more acquisition (and don't print the header)
        if (num_acquisitions%2==1) {
            // Update number of samples used
            acq_old=acq;
            free(pdata_old);
            pdata_old=pdata;
            continue;
        }
        else  {
            /*
                fprintf(stderr,"vnsr [%0.3f/%0.3f] dop[%0.1f/%0.1f] <====\n",
                    acq_old.vsnr,acq.vsnr,
                    acq_old.doppler-acq.intermediate_frequency,
                    acq.doppler-acq.intermediate_frequency
                    );
                    */
            if (acq_old.vsnr>acq.vsnr){
                acq_max=acq_old;
                free(pdata);
                pdata=pdata_old;
                // When we fix the acq.samples below,
                // use this number since that's
                // where our seek was
                use_old_samples_processed=1;
            }
            else {
                use_old_samples_processed=0;
                acq_max=acq;
            }
        }

        // Quit here if not found
        if (acq_max.vsnr<vsnr_threshold) {
            if (printAcquired!=1) {
                // fft_acquire doens't know about the seek_samples
                // offset in the data we gave it, add the offset here
                if (use_old_samples_processed){
                    acq_max.sample_number+=samples_processed_old;
                }
                else{
                    acq_max.sample_number+=samples_processed;
                }
                printAcquisitionHeader(stdout,&acq_max);
            }
            free(pdata);
            if (pdata!=pdata_old) {
                free(pdata_old);
            }
            pdata=NULL;
            pdata_old=NULL;
            // Stop if we're not doing repeated acquisitions
            if (repeated_acquisition==0)
                break;
            // If we are, go ahead Jim
            else {
                // Clear up some stuff first though
                free(acq.mag_vs_lag);
                free(acq.phase_vs_lag);
                free(acq.rss_vs_lag_I);
                free(acq.rss_vs_lag_Q);
                free(acq_old.mag_vs_lag);
                free(acq_old.phase_vs_lag);
                free(acq_old.rss_vs_lag_I);
                free(acq_old.rss_vs_lag_Q);
                continue;
            }
        }

        // Now move the initial doppler to get finer resolution
        if (finedoppler) {
            // Save off initial doppler;
            double init_dop=acq_max.doppler-acq.intermediate_frequency;
            double binhalf=acq_max.doppler_bin_size_hz/2.0;
            for (dop=-binhalf; dop<binhalf; dop+=(binhalf/FINE_BINS_SS)) {
                fft2_acquire(pdata, pprn,samples_to_use,init_dop+dop,
                    0,&acq);
                // Look for best one!
                /*
                fprintf(stderr,"vnsr [%0.3f] dop[%0.1f] <====\n",
                    acq.vsnr, acq.doppler-acq.intermediate_frequency);
                    */
                if (acq.vsnr > acq_max.vsnr) {
                    acq_max=acq;
                }
            }
        }

        free(pdata);


        fprintf(stderr,"vnsr [%0.3f] dop[%0.1f] <====\n",
            acq_max.vsnr, acq_max.doppler-acq_max.intermediate_frequency);

        // Again, add samples processed, fft doesn't know
        // how many samples we are in
        acq_max.sample_number+=samples_processed;
        printAcquisitionHeader(stdout,&acq_max);

        if (fid_correl!=NULL) {
            double current_lag=0;
            double *lag_s;
            
            // Number of lags is equal to number of samples in one code
            // period, we basically correlate with each sample lag.
            // FIXME: this lag calculation is done once here and once
            // in st_dsp.c, should consolidate!
            int num_lags=acq.num_lags;
            lag_s=(double*)malloc(num_lags*sizeof(double));
            for (i=0; i < num_lags; i++) {
                lag_s[i]=current_lag;
                current_lag+=acq.lag_spacing_s;
            }
            printCorrelWaveform(fid_correl, num_lags, lag_s, acq.mag_vs_lag, 
               acq.phase_vs_lag, phdr, &acq);

            free(lag_s);
            
        }


        // After 2nd acq. if not found, quit.
        break;
    }

    free(acq.mag_vs_lag);
    free(acq.phase_vs_lag);
    free(acq.rss_vs_lag_I);
    free(acq.rss_vs_lag_Q);

    free(acq_old.mag_vs_lag);
    free(acq_old.phase_vs_lag);
    free(acq_old.rss_vs_lag_I);
    free(acq_old.rss_vs_lag_Q);

alldone:
    //if (acq.mag_vs_lag!=NULL) free(acq.mag_vs_lag);
    //if (acq.phase_vs_lag!=NULL) free(acq.phase_vs_lag);
    free(phdr);
    free(pprn_sampled);
    free(pprn);

    return 0;
}
