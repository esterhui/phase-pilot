#include "st_dsp.h"
#include <math.h>
#include <time.h>
#include <vector>

using namespace std;

#ifdef FITTER
    #include <gsl/gsl_multifit.h>
#endif


// When using mean of rss of 'noise correlators', this
// yields a noise floor slightly higher than truth,
// here we adjust the noise by the correct ammount to
// present a realistic VSNR
// Experimentally determined with this matlab program:
// --------
// numtrails = 1000;
//
// mymean=zeros(1,numtrails);
//
// for ind=1:numtrails,
//     I=randn(1,1000*1000);
//     Q=randn(1,1000*1000);
//     mymean(ind)=mean(sqrt(Q.*Q+I.*I));
// end 
// mean(mymean)
//
#define FFT_NOISE_SCALE (1.0/1.25329)
// Function below 
// From: http://stackoverflow.com/questions/824118/why-is-floor-so-slow
#define PSEUDO_FLOOR( V ) ((V) >= 0 ? (int)(V) : (int)((V) - 1))

// Resamples complex vector x by p/q
// This uses "stephan's horrible resampling" algorithm,
// which is basically a sample_hold (if p/q > 1) or
// a CIC(sum&dump) combined with a sample hold if p/q < 1
fftw_complex* resample(fftw_complex *x, int xlen, int xlen_new) 
{
    fftw_complex *xnew, *xnew_sum_dump;
    long int j=0,k=0;
    double i=0;

    double step_size=xlen/(double)xlen_new;
    int step_size_floor=(int)floor(step_size);

    double frac_l=0;
    double frac_r=0;


    int sum_dump_len = (int)ceil(xlen/(double)step_size_floor);

    xnew = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * xlen_new );

    xnew_sum_dump = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) *
        sum_dump_len );

    // Clear out new array first
    for (j=0; j < xlen_new; j++) {
        xnew[j][0]=0;
        xnew[j][1]=0;
    }

    for (j=0; j < sum_dump_len; j++) {
        xnew_sum_dump[j][0]=0;
        xnew_sum_dump[j][1]=0;
    }

    // If we are giving more samples than we have, it's easy,
    // do simple sample-hold interpolation
    j=0;
    if (xlen_new>xlen) {
        while(j < xlen_new) {
            xnew[j][0]=x[((uint32_t)floor(i))][0];
            xnew[j][1]=x[((uint32_t)floor(i))][1];
            i+=step_size;
            j++;
        }
    }
    else if (xlen_new<xlen) {
        j=0;
        double samples_given=0;
        for (k=0; k < xlen; k++) {
            // Done with this accum?
            if (samples_given+1>step_size) {
                // frac_l is how much we give the 
                // current accumulation
                frac_l=step_size-samples_given;
                // frac_r is how much of this sample
                // the next accum gets
                frac_r=1.0-frac_l;
                
                // Current accum gets some of the next sample
                xnew[j][0]+=frac_l*x[k][0];
                xnew[j][1]+=frac_l*x[k][1];
                // Next accum gets the rest of the sample
                xnew[j+1][0]+=frac_r*x[k][0];
                xnew[j+1][1]+=frac_r*x[k][1];
                j=j+1;

                // Count the fractional of the current sample
                // we gave the new accum
                samples_given=frac_r;
                //fprintf(stdout,"samp=%d, accum=%d, samples_given=%f (F)\n",k,j,samples_given);

            }
            else {
                xnew[j][0]+=x[k][0];
                xnew[j][1]+=x[k][1];
                samples_given+=1;
                //fprintf(stdout,"samp=%d, accum=%d, samples_given=%f\n",k,j,samples_given);
            }
        }
    }
    /*
    else if (xlen_new<xlen) {
        for (j=0; j < sum_dump_len; j++) {
            for (k=0; k < step_size_floor; k++) {
                // break out early if we try to sum/dump past 
                // input vector length
                if(j*step_size_floor+k >= xlen) break;
                xnew_sum_dump[j][0]+=x[j*step_size_floor+k][0];
                xnew_sum_dump[j][1]+=x[j*step_size_floor+k][1];
            }
        }
        step_size = sum_dump_len/(double)xlen_new;
        i=0; j=0;
        while(j < xlen_new) {
            xnew[j][0]=xnew_sum_dump[((uint32_t)floor(i))][0];
            xnew[j][1]=xnew_sum_dump[((uint32_t)floor(i))][1];
            i+=step_size;
            j++;
        }
    }
    */
    else {
        xnew=x;
    }

    free(xnew_sum_dump);

    return xnew;
}

int8_t* resample_prncode(uint8_t *prn_in, uint32_t codelength, 
    uint32_t numsamples, double fs, double coderate)
{
    double step_size=coderate/fs; // Chips per sample
    double i=0;
    uint32_t j=0;
    int8_t *prn_out=NULL;

    /*
    fprintf(stderr,"prn resample step size is %lf\n",step_size);
    */

    // Allocate space for resampled PRN, we do 2* here so 
    // we can index, say 3 chips in and still grab 1023 chips
    // w/o a problem (easy way to do offsetting)
    prn_out=(int8_t*)malloc(numsamples);

    // Sanity checks
    assert(prn_in!=NULL);
    assert(prn_out!=NULL);

    while(j < numsamples) {
        prn_out[j]=(prn_in[((uint32_t)floor(i))%codelength])==0?-1:1;
        i+=step_size;
        j++;
    }
    // Make sure we asigned enough samples
    assert(j==numsamples);
    return prn_out;
}

// FFT acquisition
int fft_acquire(complex_double_t *pbuf_data,
                int8_t *pbuf_prncode,
                int numsamples,
                acq_parameters_t *acq)
{
    int i,j;
    int n;
    double mag, oldmag=0, phase=0;
    int iold=0,jold=0;
    double noise_mean=0;
    double frequency_step=0;

    double integrations_per_second = acq->fs/numsamples;

    // Number of lags is equal to number of samples in one code
    // period, we basically correlate with each sample lag.
    // FIXME: this lag calculation is done once here and once
    // in fft_acquire.c, should consolidate!
    int num_lags=(int) ceil(acq->samples_per_code_period);

    int lag_skip=(int)round(num_lags/acq->codelength/2.0);

    if (lag_skip < 1)
        lag_skip=1;

    num_lags=(int)ceil(num_lags/lag_skip);
    

    acq->num_lags=num_lags;
    acq->lag_spacing_s=lag_skip/acq->coderate;
    // Allocate space for correlation waveform
    acq->mag_vs_lag=(double*)malloc( (num_lags)*sizeof(double));
    acq->phase_vs_lag=(double*)malloc( (num_lags)*sizeof(double));
    acq->rss_vs_lag_I=(double*)malloc( (num_lags)*sizeof(double));
    acq->rss_vs_lag_Q=(double*)malloc( (num_lags)*sizeof(double));


    /*
    fprintf(stderr,"lag skip is %d\n",lag_skip);
    fprintf(stderr,"num_lags is %d\n",num_lags);
    */

    // This could be huge, but let's store a 2D delay dop map anyways
    // delay_dop_map[delay][dop]
    double **delay_dop_map_mag;
    double **delay_dop_map_phase;

    delay_dop_map_mag=(double**)malloc(num_lags*sizeof(double*));
    delay_dop_map_phase=(double**)malloc(num_lags*sizeof(double*));
    for (i=0; i<num_lags; i++) {
        delay_dop_map_mag[i]=(double*)malloc(numsamples*sizeof(double));
        delay_dop_map_phase[i]=(double*)malloc(numsamples*sizeof(double));
    }

    fftw_complex *in;
    fftw_complex *out;
    fftw_plan plan_forward;
    

    acq->doppler=0;
    acq->sample_number=0;
    acq->magnitude=0;
    acq->code_phase=0;
   
    // Setup FFT plan
    n= numsamples;

    // Compute how big our doppler bins are
    frequency_step=(acq->fs/(double)n);

    acq->doppler_bin_size_hz=frequency_step;

    //nc = ( n / 2 ) + 1;
    out = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * n );
    in = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * n );
    plan_forward = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // j = different model code delays, for each code delay
    // we multiply all the input data and then FFT
    for (j=0; j < num_lags; j++) {
        // Before we FFT the data, multiply input data with
        // model code delay
        for (i=0; i < numsamples; i++) {
            in[i][0]=(pbuf_data[i].real)*pbuf_prncode[i+j*lag_skip];      // I
            in[i][1]=(pbuf_data[i].imag)*pbuf_prncode[i+j*lag_skip];    // Q
        }
        // Now FFT
        fftw_execute ( plan_forward );
        
        // Loop through all dopplers, but skip i=0 since that's
        // the DC bin, skip last bin too for good measure?
        for ( i = 1; i < n; i++ ) {
            // Compute vector magnitude sqrt(I^2+Q^2)
            mag=sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1]);
            phase=atan2(out[i][1],out[i][0]);
            // Save into our 3d structure
            delay_dop_map_mag[j][i]=mag;
            delay_dop_map_phase[j][i]=phase;

            //if (i==3260) printf("%d %0.1f\n",j,mag);
            //if (i==2052) printf("%d %0.1f\n",j,mag);
            if (mag > oldmag){ //j==1142, i==359
                iold=i;
                jold=j;
                oldmag=mag;
                // Since this is a double sided FFT, 
                // doppler is from -fs/2 to fs/2
                acq->doppler=(i*acq->fs/(double)n);
                // Complex FFT scales from 0 to fs,
                // but fs/2 to fs is really -fs/2 to 0
                // so let's fix it here.
                if (acq->doppler > acq->fs/2.0)
                    acq->doppler-=acq->fs/2.0;
                acq->code_phase=acq->codelength
                    -(j*lag_skip/acq->samples_per_code_period
                      *(double)(acq->codelength));
                acq->sample_number=round(
                    acq->samples_per_code_period-j*lag_skip);
                acq->magnitude=mag;
            }
        }
    }

    // Now save off correlation waveform corresponding to the doppler
    // where the max peak was detected.
    for (j=0; j< num_lags; j++) {
        acq->mag_vs_lag[j]=delay_dop_map_mag[j][iold];
        acq->phase_vs_lag[j]=delay_dop_map_phase[j][iold];
        //fprintf (stdout, "%d %f \n",j,acq->mag_vs_lag[j]);
    }

    // Now compute the mean noise floor (sum over 1/2 of the lags, but
    // start the lag summing 0.25*total_lags from the detected peak to
    // ensure we don't sum some of the peak)
    for (j=0; j< num_lags/2; j++) {
        noise_mean+=acq->mag_vs_lag[((jold+j+num_lags/4)%num_lags)];
    }
    acq->noise_floor=noise_mean/(num_lags/2);

    // Compute 1sVSNR
    // 1sVSNR = magnitude/noise_floor * sqrt(integrations per second).
    // The last term scales the SNR to a 1sVSNR
    acq->vsnr = acq->magnitude/acq->noise_floor;
    acq->vsnr_1s = acq->vsnr*sqrt(integrations_per_second);

    // 6-sigma
    if (acq->vsnr >= acq->vsnr_threshold) {
        acq->found=1;
    }
    else {
        acq->found=0;
    }

    acq->carrier_phase=0; // FIXME, should make estimate at some point

    // --------- output for other programs -----------
    /*
    fprintf (stdout, "i/j/n            : %d %d %d\n",iold,jold,n);
    fprintf (stdout, "PRN              : %02d\n",acq->prn);
    fprintf (stdout, "Doppler          : %0.2f Hz\n",acq->doppler);
    fprintf (stdout, "Sample Offset    : %d\n",acq->sample_number);
    fprintf (stdout, "Code Phase       : %0.2f chips\n",acq->code_phase);
    fprintf (stdout, "Peak             : %0.1f\n",acq->magnitude);
    */
    /*
    fprintf (stdout,"%02d\t%0.2f\t%d\t%0.2f\t%0.1f\n",
        acq->prn,p.doppler,p.sample_number,p.code_phase,
        p.magnitude);
    fflush(stdout);
    */

    for (i=0; i<num_lags; i++) {
        free(delay_dop_map_mag[i]);
        free(delay_dop_map_phase[i]);
    }
    free(delay_dop_map_mag);
    free(delay_dop_map_phase);

    fftw_destroy_plan(plan_forward);
    fftw_free(out);
    fftw_free(in);

    return acq->found;
}

// Searial acquisition
int serial_acquire(complex_double_t *pbuf_data,
                uint8_t *pprn,
                int numsamples,
                double initial_doppler_hz,
                float ss_doppler_span_hz, 
                acq_parameters_t *acq)
{
    int i,j;
    double mag, oldmag=0, phase=0;
    int jold=0;
    double noise_mean=0;
    double bin_size_hz=0;
    int num_bins = 0;
    int cur_bin=0, cur_sign=1;
    double I,Q, cur_doppler_hz;

    cr_t CR_in;
    cr_t *CR_out;
    corr_t corr;

    double lag_spacing_chips = 1/10.0;
    double integrations_per_second = acq->fs/numsamples;

    int num_lags=(int)ceil(acq->codelength/lag_spacing_chips);
    fprintf(stderr,"num_lags is %d\n",num_lags);

    bin_size_hz = integrations_per_second/2.0;

    acq->doppler_bin_size_hz=bin_size_hz;

    //fprintf(stderr,"bin size is %0.1lf Hz\n",bin_size_hz);
    num_bins = (int)round(2.0*ss_doppler_span_hz/bin_size_hz);
    //fprintf(stderr,"num bins is %d\n",num_bins);

    acq->doppler=0;
    acq->sample_number=0;
    acq->magnitude=0;
    acq->code_phase=0;
   
    CR_in.fs=acq->fs;
    CR_in.start_phase=0;
    CR_in.num_samples=numsamples;
    CR_in.data=pbuf_data;

    corr.number=1;
    corr.start_chips=0; 
    corr.start_chips_orig=0;
    corr.spacing_chips=0;

    acq->num_lags=num_lags;
    acq->lag_spacing_s = lag_spacing_chips / acq->coderate;
    // Allocate space for correlation waveform
    acq->mag_vs_lag=(double*)malloc( (num_lags)*sizeof(double));
    acq->phase_vs_lag=(double*)malloc( (num_lags)*sizeof(double));
    acq->rss_vs_lag_I=(double*)malloc( (num_lags)*sizeof(double));
    acq->rss_vs_lag_Q=(double*)malloc( (num_lags)*sizeof(double));

    // j = different model code delays, for each code delay
    // we multiply all the input data and then FFT
    for (i=0; i < num_bins; i++) {
        cur_doppler_hz = initial_doppler_hz+cur_sign*cur_bin*bin_size_hz;
        fprintf(stderr,"Searching doppler %0.1f Hz\n",cur_doppler_hz);
        CR_in.fc=acq->intermediate_frequency+cur_doppler_hz;
        // FIXME: Should free data somewhere
        CR_out=counterRotate(&CR_in);
        // Flip signs
        cur_sign*=-1;
        // Increment doppler bin if we're done
        if (i%2==0) {
            cur_bin++;
        }

        // Now process each lag
        for (j=0; j < num_lags; j++) {

            // Setup the new correlation point
            corr.start_chips = j*lag_spacing_chips;
            
            // Correlate
            prnMultiplyAccumulate(&corr,CR_out,pprn,acq);
            
            // Get I/Q
            I=corr.data[0].real;
            Q=corr.data[0].imag;
            
            // Compute vector magnitude sqrt(I^2+Q^2)
            mag=sqrt(I*I+Q*Q);
            phase=atan2(Q,I);
            // Save into our 3d structure
            acq->mag_vs_lag[j]=mag;
            acq->phase_vs_lag[j]=phase;
            acq->rss_vs_lag_I[j]=sqrt(I*I);
            acq->rss_vs_lag_Q[j]=sqrt(Q*Q);

            if (mag > oldmag){ //j==1142, i==359
                jold=j;
                oldmag=mag;
                // Since this is a double sided FFT, 
                // doppler is from -fs/2 to fs/2
                acq->doppler=CR_in.fc;
                /*
                acq->code_phase=acq->codelength -(corr.start_chips);
                acq->sample_number=round(acq->samples_per_code_period*
                    (acq->code_phase/acq->codelength));
                    */
                acq->code_phase=(corr.start_chips);
                acq->sample_number=round(acq->samples_per_code_period*(
                    (acq->code_phase/acq->codelength)));
                acq->magnitude=mag;
            }
        }

        // Now compute the mean noise floor (sum over 1/2 of the lags, but
        // start the lag summing 0.25*total_lags from the detected peak to
        // ensure we don't sum some of the peak)
        noise_mean=0;
        for (j=0; j< num_lags/2; j++) {
            noise_mean+=acq->rss_vs_lag_I[((jold+j+num_lags/4)%num_lags)];
            //fprintf(stderr,"%f\n",acq->rss_vs_lag[((jold+j+num_lags/4)%num_lags)]);
        }
        acq->noise_floor=noise_mean/(num_lags/2);

        // Compute 1sVSNR
        // 1sVSNR = magnitude/noise_floor * sqrt(integrations per second).
        // The last term scales the SNR to a 1sVSNR
        acq->vsnr = acq->magnitude/acq->noise_floor;
        acq->vsnr_1s = acq->vsnr*sqrt(integrations_per_second);

        // 6-sigma
        if (acq->vsnr >= acq->vsnr_threshold) {
            acq->found=1;
        }
        else {
            acq->found=0;
        }

        acq->carrier_phase=0; // FIXME, should make estimate at some point

        if (acq->found) break;
    }

    return acq->found;
}

fftw_complex* complexMultiply(fftw_complex a, fftw_complex b)
{
    fftw_complex *temp;
    temp = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * 1 );

    //  (a+bi) (c+di) = (ac-bd) + i(bc+ad)
    (*temp)[0]=a[0]*b[0]-(a[1]*b[1]); // Real part = a.real*b.real-a.imag*b.imag
    (*temp)[1]=a[1]*b[0]+(a[0]*b[1]); // Imag part = a.imag*b.real+a.real*b.imag

    return temp;
}

int fft2_acquire(complex_double_t *pbuf_data,
                uint8_t *pprn,
                int numsamples,
                double initial_doppler_hz,
                float ss_doppler_span_hz, 
                acq_parameters_t *acq)
{
    int i,j;
    fftw_complex *x,*x_resampled, *x_fft;
    fftw_complex *prn_resampled, *prn_fft;
    fftw_complex *prn_resampled2;
    int8_t *prn_resampled_int;
    double I,Q,mag,phase,oldmag=0;
    int jold=0;
    complex_double_t *pbuf_data_conj=NULL;

    double integrations_per_second = acq->fs/numsamples;

    // First compute FFT length
    float numperiods=numsamples/acq->samples_per_code_period;
    float fftN=2.0*acq->codelength*numperiods;
    int fft_length=pow(2,(int)ceil(log2(fftN)));

    /*
    fprintf(stderr,"FFT length is %d\n",fft_length);
    fprintf(stderr,"code length is %d\n",acq->codelength);
    fprintf(stderr,"periods is %f\n",numperiods);
    fprintf(stderr,"initial dop is  %f\n",initial_doppler_hz);
    */


    float frequency_step = floor(acq->coderate/acq->codelength/ 
        numperiods);

    acq->doppler_bin_size_hz=frequency_step;

    int frequency_bins = (int)(2*round(ss_doppler_span_hz/frequency_step)+1);

    /*
    fprintf(stderr,"frequency step is %0.1f Hz\n",frequency_step);
    fprintf(stderr,"frequency bins is %d\n",frequency_bins);
    */

    double samples_per_chip = fft_length/(numperiods*acq->codelength);

    // one_fft_length is how many non-repeating lags we have
    // if integration time is 2x code period, we stack 2 PRN codes
    // together, thus when we do fft acquire there will be two peaks
    // one_fft_length tells the software to only search the first
    // unique segment.

    int one_fft_length;
    if (numperiods>=1) {
        one_fft_length = (int)round(samples_per_chip*acq->codelength);
    }
    else {
        one_fft_length = (int)round(samples_per_chip*acq->codelength*numperiods);
    }

    /*
    fprintf(stderr,"samples per chip is %lf\n",samples_per_chip);
    */

    int num_lags = one_fft_length;

    /*
    fprintf(stderr,"number of lags is %d\n",num_lags);
    */

    acq->num_lags=one_fft_length;
    if (numperiods>=1)
        acq->lag_spacing_s=acq->codelength/(double)acq->num_lags/acq->coderate;
    else
        acq->lag_spacing_s=acq->codelength*numperiods/(double)acq->num_lags/acq->coderate;
    // Allocate space for correlation waveform
    acq->mag_vs_lag=(double*)malloc( (num_lags)*sizeof(double));
    acq->phase_vs_lag=(double*)malloc( (num_lags)*sizeof(double));
    acq->rss_vs_lag_I=(double*)malloc( (num_lags)*sizeof(double));
    acq->rss_vs_lag_Q=(double*)malloc( (num_lags)*sizeof(double));

    cr_t CR_in;
    cr_t *CR_out;

    // ---------------------------------------------------------------------- 
    // ----------------- Remove nominal IF from samples ---------------------
    // ---------------------------------------------------------------------- 

    pbuf_data_conj=(complex_double_t*)malloc(
        numsamples*(sizeof(complex_double_t)));
    for (i=0; i < numsamples; i++) {
        pbuf_data_conj[i].imag=-pbuf_data[i].imag;
        pbuf_data_conj[i].real=pbuf_data[i].real;
        //fprintf(stdout,"%d %lf %lf\n",i,pbuf_data[i].real, pbuf_data[i].imag);
    }

    CR_in.fs=acq->fs;
    CR_in.start_phase=0;
    CR_in.num_samples=numsamples;
    CR_in.data=pbuf_data_conj;

    CR_in.fc=-(acq->intermediate_frequency+initial_doppler_hz);

    CR_out=counterRotate(&CR_in);

    // ---------------------------------------------------------------------- 
    // ---------------- Resample to fft_length & fft ------------------------
    // ---------------------------------------------------------------------- 

    x = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * numsamples );
    x_fft = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * fft_length );

    // Convert to fftw_complex
    for (i=0; i < numsamples; i++) {
        x[i][0]=CR_out->data[i].real;
        x[i][1]=CR_out->data[i].imag;
        //fprintf(stdout,"%d %lf %lf\n",i,x[i][0], x[i][1]);
    }

    free(CR_out->data);
    free(CR_out);


    x_resampled=resample(x, numsamples, fft_length);
    assert(x_resampled!=NULL);

    /*
    for (i=0; i < numsamples; i++) {
        fprintf(stdout,"%d %lf %lf\n",i,x[i][0], x[i][1]);
    }
    */

    // And FFT the resampled data
    fftw_plan p;
    p = fftw_plan_dft_1d(fft_length, x_resampled, x_fft, FFTW_FORWARD,
        FFTW_ESTIMATE);
    fftw_execute(p);

    /*
    for (i=0; i < fft_length; i++) {
        fprintf(stdout,"%d %lf %lf\n",i,x_fft[i][0], x_fft[i][1]);
    }
    */

    // ---------------------------------------------------------------------- 
    // ---------------- Now resample & FFT prn code -------------------------
    // ---------------------------------------------------------------------- 
    /*
    prn_resampled_int=resample_prncode(pprn, acq->codelength, 
        fft_length, samples_per_chip*acq->coderate, acq->coderate);
        */
    prn_resampled_int=resample_prncode(pprn, acq->codelength, 
        numsamples, acq->fs, acq->coderate);


    //fprintf(stderr,"numsamp is %d, fft_length is %d\n",numsamples,fft_length);
    prn_resampled= (fftw_complex*)fftw_malloc( 
            sizeof ( fftw_complex ) * numsamples );
    prn_fft = (fftw_complex*)fftw_malloc ( 
        sizeof ( fftw_complex ) * fft_length );

    // Convert to complex (imagnary part = 0)
    //for (i=0; i < fft_length; i++) {
    for (i=0; i < numsamples; i++) {
        prn_resampled[i][0]=(double)prn_resampled_int[i];
        prn_resampled[i][1]=0;
        //fprintf(stdout,"%d %lf %lf\n",i,prn_resampled[i][0], prn_resampled[i][1]);
    }
    free(prn_resampled_int);

    prn_resampled2=resample(prn_resampled, numsamples, fft_length);

    /*
    for (i=0; i < fft_length; i++) {
        fprintf(stdout,"%d %lf %lf %lf\n",i,x_resampled[i][0], x_resampled[i][1],prn_resampled2[i][0]);
    }
    */


    fftw_destroy_plan(p);
    //p = fftw_plan_dft_1d(fft_length, prn_resampled, prn_fft, FFTW_FORWARD,
    p = fftw_plan_dft_1d(fft_length, prn_resampled2, prn_fft, FFTW_FORWARD,
        FFTW_ESTIMATE);
    fftw_execute(p);

    // Now take conjugate of prn_fft 
    // conj(x) = real(x)-j*imag(x)
    for (i=0; i < fft_length; i++) {
        prn_fft[i][1]*=-1.0;
    }

    
    /*
    for (i=0; i < fft_length; i++) {
        fprintf(stdout,"%d %lf %lf\n",i,prn_fft[i][0], prn_fft[i][1]);
    }
    */
    // ---------------------------------------------------------------------- 
    // ---------------- Multiply in freq domain & ifft ----------------------
    // ---------------------------------------------------------------------- 

    double cur_doppler_hz;
    int cur_bin=0, cur_sign=1;
    fftw_complex *temp;

    fftw_complex *fft_correlation= (fftw_complex*)fftw_malloc( 
            sizeof ( fftw_complex ) * fft_length );
    fftw_complex *time_correlation= (fftw_complex*)fftw_malloc( 
            sizeof ( fftw_complex ) * fft_length );

    
    int cur_fft_bin=0;

    for (i=0; i <  frequency_bins; i++) {
        cur_doppler_hz = initial_doppler_hz+cur_sign*cur_bin*frequency_step;
        //fprintf(stderr,"Searching doppler %0.1f Hz\n",cur_doppler_hz);

        // Multiply prn_fft with x_fft point by point
        for (j=0; j < fft_length; j++) {
            cur_fft_bin=(j-(cur_sign*cur_bin))%fft_length;
            if (cur_fft_bin < 0) {
                cur_fft_bin = fft_length+cur_fft_bin;
            }
            temp=complexMultiply(prn_fft[j],x_fft[cur_fft_bin]);
            fft_correlation[j][0]=(*temp)[0];
            fft_correlation[j][1]=(*temp)[1];
            /*
            fprintf(stdout,"%d %lf %lf\n",j,fft_correlation[j][0],
                fft_correlation[j][1]);
                */
            free(temp);
        }

        // Now ifft ---------
        fftw_destroy_plan(p);
        p = fftw_plan_dft_1d(fft_length, fft_correlation, time_correlation, 
            FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_execute(p);

        /*
        for (j=0; j < fft_length; j++) {
            fprintf(stdout,"%d %lf %lf\n",j,time_correlation[j][0],
                time_correlation[j][1]);
        }
        */
        for (j=0; j < one_fft_length; j++) {
            I=time_correlation[j][0];
            Q=time_correlation[j][1];

            mag=sqrt(I*I+Q*Q);
            phase=atan2(Q,I);
            // Save into our 3d structure
            acq->mag_vs_lag[j]=mag;
            acq->phase_vs_lag[j]=phase;

            if (mag > oldmag){
                oldmag=mag;
                jold=j;
                // Since this is a double sided FFT, 
                // doppler is from -fs/2 to fs/2
                acq->doppler=acq->intermediate_frequency+cur_doppler_hz;
                if (numperiods>=1)
                    acq->code_phase=j/(double)one_fft_length*acq->codelength;
                else
                    acq->code_phase=j/(double)one_fft_length*acq->codelength*numperiods;
                acq->sample_number=round(acq->samples_per_code_period*(
                    (acq->code_phase/acq->codelength)));
                acq->magnitude=mag/fft_length;
            }
        }

        // Flip signs
        cur_sign*=-1;
        // Increment doppler bin if we're done
        if (i%2==0) {
            cur_bin++;
        }

        // Now compute the mean noise floor (sum over 1/2 of the lags, but
        // start the lag summing 0.25*total_lags from the detected peak to
        // ensure we don't sum some of the peak)
        double noise_mean=0;
        for (j=0; j< num_lags*0.75; j++) {
            noise_mean+=acq->mag_vs_lag[((jold+j+num_lags/8)%num_lags)];
        }
        acq->noise_floor=noise_mean/(num_lags*0.75)/fft_length*FFT_NOISE_SCALE;

        // Compute 1sVSNR
        // 1sVSNR = magnitude/noise_floor * sqrt(integrations per second).
        // The last term scales the SNR to a 1sVSNR
        if (acq->noise_floor <=0) {
            acq->vsnr = 0;
        }
        else {
            acq->vsnr = acq->magnitude/acq->noise_floor;
        }
        acq->vsnr_1s = acq->vsnr*sqrt(integrations_per_second);

        if (acq->vsnr >= acq->vsnr_threshold) {
            acq->found=1;
        }
        else {
            acq->found=0;
        }

        if (acq->found) break;
    }
    
    // FIXME: many other mallocs not freed here
    fftw_destroy_plan(p);
    free(pbuf_data_conj);
    free(x);
    free(x_fft);
    free(prn_fft);
    free(prn_resampled);
    free(prn_resampled2);
    free(fft_correlation);
    free(time_correlation);
    free(x_resampled);

    return acq->found;
}

cr_t* counterRotate(cr_t* pcr)
{
    unsigned int i;
    double phase=0;
    float Ireal, Iimag, Qreal, Qimag;

    // Compute Sample period T in seconds
    double Ts=1.0/pcr->fs;
    double t=0;
    double start_phase=0;

    cr_t* pcr_out=(cr_t*)malloc(sizeof(cr_t));
    pcr_out->data=(complex_double_t*)malloc(
        pcr->num_samples*sizeof(complex_double_t));

    pcr_out->fs=pcr->fs;
    pcr_out->fc=pcr->fc;
    pcr_out->start_phase=pcr->start_phase;
    pcr_out->num_samples=pcr->num_samples;
    
    // Build up phase and mix
    double fc2pi=pcr->fc*2.0*M_PI; // 2*pi*fc
    for (i=0; i < pcr->num_samples; i++) {
        phase=(fc2pi*t)+pcr->start_phase;
        if (i==0) start_phase=phase;

        t+=Ts; // Increment time

        Ireal=(cos(phase))*(pcr->data[i].real);
        Iimag=(-sin(phase))*(pcr->data[i].real);
        Qreal=(sin(phase))*(pcr->data[i].imag);
        Qimag=(cos(phase))*(pcr->data[i].imag);

        pcr_out->data[i].real=Ireal+Qreal;
        pcr_out->data[i].imag=Iimag+Qimag;
    }

    // Compute integrated phase in cycles
    pcr_out->integrated_phase=pcr->fc*Ts*pcr->num_samples;
    
    // Next phase will be current phase plus 1 clock tick
    // FIXME, for many hour data sets, we might run out of phase resolution
    // since phase here keeps incrementing
    pcr_out->phase_next=phase+(2.0*M_PI*pcr->fc*Ts);

    return pcr_out;
}


complex_double_t accumulate(complex_double_t* data, unsigned int num_samples)
{
    unsigned int i;
    complex_double_t accum;
    accum.real=0;
    accum.imag=0;

    for (i=0; i < num_samples; i++){
        //printf("[%d %lf %lf]\n",i,data[i].real,data[i].imag);
        accum.real+=data[i].real;
        accum.imag+=data[i].imag;
    }

    return accum;
}

int loop_init(loop_t* params,double integration_time, double BW_phase, double BW_freq)
{
    params->T=integration_time;

    params->acc_delay=0.0;
    params->vel_delay=0.0;
    params->freq_error=0.0;
    params->phase_error=0.0;

    /* First compute PLL multipliers */
    /* From GPS: Principles and applications "Loop Filters" */
    params->Wnp     = BW_phase/0.7845; // [rad/sec]
    params->Wn2p    = params->Wnp*params->Wnp;
    params->Wn3p    = params->Wnp*params->Wnp*params->Wnp;
    params->a3      = 1.1;
    params->b3      = 2.4;

    /* Now compute FLL multipliers */
    params->Wnf     = BW_freq/0.53;
    params->Wn2f    = params->Wnf*params->Wnf;
    params->a2      = sqrt(2)*params->Wnf;

    return 1;
}


// 3rd order loop, tracks freq & acceleration
double loop_run3(loop_t* params)
{
    double freq_adj;

    // ---- Setup phase inputs ----
    double p1=params->phase_error*params->Wn3p*params->T;
    double p2=params->phase_error*params->a3*params->Wn2p;
    double p3=params->phase_error*params->b3*params->Wnp;
    
    // ---- Setup freq inputs ----
    double f1=params->freq_error*params->a2*params->Wnf;
    double f2=params->freq_error*params->Wn2f*params->T;
    
    // ---------- Acceleration accumulation ----------

    double acc1 = f2+p1+params->acc_delay; // First summation
    double acc2 = 0.5*(acc1+params->acc_delay); // 2nd accumulation
    params->acc_delay = acc1; // Delay block
    double acc3 = f1 + acc2 + p2; // 3rd accumulation

    // ---------- Velocity accumulation ----------
    double vel1 = (params->T*acc3)+params->vel_delay;
    double vel2 = 0.5*(vel1+params->vel_delay);
    params->vel_delay = vel1; // Delay block

    // Compute frequency adjustment necessary (in Hz)
    freq_adj = vel2+p3;

    return freq_adj;
}

int prnMultiplyAccumulate(corr_t* corr, cr_t* pcr, uint8_t* pprn, 
    acq_parameters_t* acq)
{
    unsigned int samplenum;
    unsigned int corrnum;
    unsigned int prn_ind;
    int chip;

    // Reserve space for all correlator phases
    double* corr_phase=(double*)malloc(corr->number*sizeof(double));


    // Always mod
    corr->start_chips=fmod(corr->start_chips,acq->codelength);

    // Account for negative phase
    if (corr->start_chips<0) corr->start_chips+=acq->codelength;

    // Allocate space for correlator output
    corr->data=(complex_double_t*)malloc(sizeof(*corr->data)*corr->number);
    for (corrnum=0; corrnum < corr->number; corrnum++){
        corr->data[corrnum].real=0;
        corr->data[corrnum].imag=0;
        
        // We use corr_phase as our index into the pprn vector
        // if the phase in chips is 23.3
        // we will be at chip[23] fraction 0.3
        corr_phase[corrnum]=fmod(acq->codelength-
            fmod(corr->start_chips+corrnum*corr->spacing_chips,
            acq->codelength),acq->codelength);
    }

    for (corrnum=0; corrnum < corr->number; corrnum++){
        for (samplenum=0; samplenum < pcr->num_samples; samplenum++){
            // prn_ind is the integer chips to index into
            // our (eg. 1023 chip) prn code
            prn_ind=PSEUDO_FLOOR(corr_phase[corrnum]);

            // Convert 1 0 to +1 -1
            chip=pprn[prn_ind]==1?+1:-1;
            corr->data[corrnum].real+=(pcr->data[samplenum].real*chip);
            corr->data[corrnum].imag+=(pcr->data[samplenum].imag*chip);

            // Advance phase by 1 step
            corr_phase[corrnum]+=acq->chips_per_sample;
            if (corr_phase[corrnum] > acq->codelength) {
                corr_phase[corrnum]-=acq->codelength;
            }
        }
    }

    // The next phase should be this
    corr->chips_next = -corr_phase[0];

    // Keep track of integrated phase/chips here for the interval
    corr->integrated_phase=acq->chips_per_sample*pcr->num_samples;

    free(corr_phase);

    return 1;
}

// Obtained from:
// http://www.strchr.com/standard_deviation_in_one_pass
double compute_snr(double amplitude, double ampl_window[], 
                int numsamp, double *noise) {
    int i;
    if(numsamp == 0)
        return 0.0;
    double sum = 0;
    for(i = 0; i < numsamp; i++)
       sum += ampl_window[i];
    double mean = sum / numsamp;
    double sq_diff_sum = 0;
    for(i = 0; i < numsamp; i++) {
       double diff = ampl_window[i] - mean;
       sq_diff_sum += diff * diff;
    }
    double variance = sq_diff_sum / numsamp;

    /*
    fprintf(stderr,"var=%lf, [%lf,%lf,%lf...%lf]\n",variance,ampl_window[0],
        ampl_window[1],ampl_window[2],ampl_window[numsamp]);
        */
    *noise=sqrt(variance);

    return mean/(*noise);
}

double dsign(double number)
{
    return(number>=0?+1:-1);
}


#ifdef FITTER
double fitter(double *timetag_s, double *value, int nelements, 
    int order, double *chisq, double *coef, double eval_at_time, int verbose)
{
    int i=0, n, j=0;
    double fitted_value;
    double xi, yi;
    gsl_matrix *X, *cov;
    gsl_vector *y, *w, *c;

    // order of 2 actually has 3 coefficients,
    // eg x^2, x^1, x^0.
    int gsl_order = order+1;

    // Sanity check, make sure the time
    // we asked it to evaluate lies in the
    // data we gave it
    if ((eval_at_time) < (timetag_s[0])) {
        return NAN;
    }
    if((eval_at_time) > (timetag_s[nelements-1])){
        return NAN;
    }

    n = nelements;

    X = gsl_matrix_alloc (n, gsl_order);
    y = gsl_vector_alloc (n);
    w = gsl_vector_alloc (n);

    c = gsl_vector_alloc (gsl_order);
    cov = gsl_matrix_alloc (gsl_order, gsl_order);

    for (i = 0; i < n; i++)
    {
        xi=timetag_s[i];
        yi=value[i];

        for(j=0; j < gsl_order; j++) {
            gsl_matrix_set(X, i, j, pow(xi, j));
        }

        gsl_vector_set (y, i, yi);
    }

    gsl_multifit_linear_workspace * work 
        = gsl_multifit_linear_alloc (n, gsl_order);
    gsl_multifit_linear (X, y, c, cov,
        chisq, work);
    gsl_multifit_linear_free (work);

#define C(i) (gsl_vector_get(c,(i)))

    if (verbose) {
        fprintf (stderr,"----------------\n");
        fprintf (stderr,"    # x[eval] = %g x[start] = %g x[end]= %g\n", 
            eval_at_time,timetag_s[0],timetag_s[n-1]);
        fprintf (stderr,"    # chisq = %g [npts=%d]\n", *chisq,n);
        if (order > 2) {
            fprintf (stderr,
                "    # best fit: Y = %g + %g X + %g X^2 + %g X^3\n", 
            C(0), C(1), C(2), C(3));
        }
        if (order > 1) {
            fprintf (stderr,
                "    # best fit: Y = %g + %g X + %g X^2\n", 
            C(0), C(1), C(2));
        }
        else {
            fprintf (stderr,"    # best fit: Y = %g + %g X \n", 
            C(0), C(1));
        }
    }

    /* Get coefficients */
    for (i=0; i <= order; i++) {
        coef[i]=C(i);
    }

    /* Now compute the value at requested timetag */
    if (order==1) {
        fitted_value = C(0)+C(1)*eval_at_time;
    }
    else {
        fitted_value = C(0)+C(1)*eval_at_time+C(2)*pow(eval_at_time,2);
    }
    // FIXME: > 2nd order still seems to screw up

    gsl_matrix_free (X);
    gsl_vector_free (y);
    gsl_vector_free (w);
    gsl_vector_free (c);
    gsl_matrix_free (cov);

    return fitted_value;
}
#endif

int fft_acquire_tone(complex_double_t *pbuf_data,
                int numsamples,
                acq_parameters_t *acq)
{
    int i;
    int n;
    double mag=0, oldmag=0, oldmag_min=INFINITY;
    double frequency_step=0;
    double temp_dop=0;

    double integrations_per_second = acq->fs/numsamples;

    fftw_complex *in;
    fftw_complex *out;
    fftw_plan plan_forward;
    
    acq->doppler=0;
    acq->sample_number=0;
    acq->magnitude=0;
    acq->code_phase=0;
   
    // Setup FFT plan
    n= numsamples;

    // Compute how big our doppler bins are
    frequency_step=(acq->fs/(double)n);

    acq->doppler_bin_size_hz=frequency_step;

    //nc = ( n / 2 ) + 1;
    out = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * n );
    in = (fftw_complex*)fftw_malloc ( sizeof ( fftw_complex ) * n );
    plan_forward = fftw_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for (i=0; i < numsamples; i++) {
        in[i][0]=(pbuf_data[i].real);      // I
        in[i][1]=(pbuf_data[i].imag);    // Q
    }
    // Now FFT
    fftw_execute ( plan_forward );
        
    acq->noise_floor=0;
    // Loop through all dopplers, but skip i=0 since that's
    // the DC bin, skip last bin too for good measure?
    for ( i = 1; i < n-1; i++ ) {
        // Compute vector magnitude sqrt(I^2+Q^2)
        mag=sqrt(out[i][0]*out[i][0]+out[i][1]*out[i][1]);
        acq->noise_floor+=mag;

        // doppler is from -fs/2 to fs/2
        temp_dop=(i*acq->fs/(double)n);
        // Complex FFT scales from 0 to fs,
        // but fs/2 to fs is really -fs/2 to 0
        // so let's fix it here.
        if (temp_dop > acq->fs/2.0)
            temp_dop-=acq->fs;

        //fprintf(stderr,"%lf %f\n",temp_dop,mag);

        if (mag > oldmag){ //j==1142, i==359
            // Since this is a double sided FFT, 
            // doppler is from -fs/2 to fs/2
            temp_dop=(i*acq->fs/(double)n);
            // Complex FFT scales from 0 to fs,
            // but fs/2 to fs is really -fs/2 to 0
            // so let's fix it here.
            if (temp_dop > acq->fs/2.0)
                temp_dop-=acq->fs;

            // Only save if not DC and not +-fs/2 and
            // within doppler search space (or ignore if
            // search space is 0)
            if ((temp_dop != 0) &&
                (temp_dop < acq->fs/2.0 && temp_dop > -acq->fs/2.0) &&
                (acq->doppler_search_space_ss==0||(fabs(
                    temp_dop-acq->intermediate_frequency)<
                    acq->doppler_search_space_ss))) {
                //fprintf(stderr,"saving mag %f at doppler %f\n",mag,temp_dop);
                acq->magnitude=mag;
                oldmag=mag;
                acq->doppler=temp_dop;
            }
        }
    }
    acq->noise_floor/=n;

    // Compute 1sVSNR
    // 1sVSNR = magnitude/noise_floor * sqrt(integrations per second).
    // The last term scales the SNR to a 1sVSNR
    acq->vsnr = acq->magnitude/acq->noise_floor;
    acq->vsnr_1s = acq->vsnr*sqrt(integrations_per_second);

    //fprintf(stderr,"VSNR is %0.1f, thresh is %f\n",acq->vsnr,acq->vsnr_threshold);
    // 6-sigma and do not track anything that looks like DC
    if (acq->vsnr >= acq->vsnr_threshold) {
        acq->found=1;
    }
    else {
        acq->found=0;
    }

    acq->carrier_phase=0; // FIXME, should make estimate at some point

    // --------- output for other programs -----------
    /*
    fprintf (stdout, "i/n              : %d %d\n",iold,n);
    fprintf (stdout, "Doppler          : %0.2f Hz\n",acq->doppler);
    fprintf (stdout, "Peak             : %0.1f\n",acq->magnitude);
    fprintf (stdout, "Noise            : %0.1f\n",acq->noise_floor);
    */

    fftw_destroy_plan(plan_forward);
    fftw_free(out);
    fftw_free(in);

    return acq->found;
}

#ifdef COMPUTE_STATISTICS
stats_signal* computeStatistics(st_accum_tone *d, double rf_frequency, int antenna)
{
    vec rss;

    stats_signal *stats = new stats_signal;

    stats->antenna_number=antenna;
    stats->rf_frequency_hz=rf_frequency;

    // If no data or not enough samples
    if (d==NULL || d->timetag_s.size()<1) {
        // Indicate things are bad
        stats->cn0=NAN;
        return stats;
    }

    rss=sqrt(square(d->ip)+square(d->qp));

    double fs=d->fs;
    double numsamp_mean=mean(d->nsamp);
    double integ_time_s=numsamp_mean/fs;

    double rss_std_expected=sqrt(numsamp_mean);
    double rss_mean=mean(rss);
    double vsnr1s=rss_mean/rss_std_expected*sqrt(1/integ_time_s);
    double cn0=20.0*log10(vsnr1s)-3.0;

    if (d->timetag_s.size()>1) {
        stats->rss_std=stddev(rss);
        stats->doppler_std_hz=stddev(d->doppler_hz);
    }
    else {
        stats->rss_std=0;
        stats->doppler_std_hz=0;
    }

    stats->rss_mean=rss_mean;
    stats->integ_time_s=integ_time_s;
    stats->vsnr1s=vsnr1s;
    stats->cn0=cn0;

    stats->doppler_mean_hz=mean(d->doppler_hz);

    /* ------ Determine phase ---------- */
    double chisq;
    double coef[3];
    stdvec timetag_s_vector,phase_s_vector;

    // convert to normal vectors
    timetag_s_vector=conv_to<stdvec>::from(d->timetag_s);
    phase_s_vector=conv_to<stdvec>::from(d->phase_s);

    int verbose=0;
    int order=2;
    uint32_t num_elements=timetag_s_vector.size();

    stats->num_accum=num_elements;

    if (num_elements>3) {
        // run fitter to get coefficients
        fitter(timetag_s_vector.data(),phase_s_vector.data(),
            num_elements,order,&chisq,coef,NAN,verbose);

        // Now remove a polynomial from the phase data
        stdvec *phase_vector=removePoly(timetag_s_vector.data(),phase_s_vector.data(),
            num_elements,coef,order);
        // Convert back to good ol' arma magic
        d->phase_s=conv_to<vec>::from(*phase_vector);

        stats->phase_ns=coef[0]*1e9; // ns
        stats->phase_rate=coef[1]*1e9; // ns/s
        stats->phase_rate_rate=coef[2]*1e9; // ns/s/s
        stats->phase_std_ns=stddev(d->phase_s)*1e9; 
        // Compute phase noise, JPL Pub 95-6, eq 7.14
        double cycle_length_ns=1.0/rf_frequency*1e9;
        double phase_scatter_cyc=stats->phase_std_ns/cycle_length_ns;
        stats->phase_vsnr1s=1.0/(2.0*M_PI*phase_scatter_cyc)*sqrt(1/integ_time_s);;
    }

    return stats;
}

// Loop through data and remove polynomial described by coeff, return
// a new vector
stdvec* removePoly(double *time_s, double *data, uint32_t num_elements, double *coef,  int order)
{

    stdvec *newdata = new stdvec;
    double fitted_value;
    double t;
    for (uint32_t i=0; i < num_elements; i++) {
        // Time in seconds;
        t=time_s[i];
        if (order==1) {
            fitted_value = coef[0]+coef[1]*t;
        }
        else if (order==2){
            fitted_value = coef[0]+coef[1]*t+coef[2]*pow(t,2);
        }
        else {
            fitted_value = coef[0]+coef[1]*t+coef[2]*pow(t,2)
                +coef[3]*pow(t,3);
        }

        newdata->push_back(data[i]-fitted_value);
    }

    return newdata;
}
#endif
