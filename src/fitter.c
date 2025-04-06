/* Fitter - fits a quadratic to given data
 *
 *
 * $Id: fitter.c 110 2013-04-07 22:27:03Z esterhui $
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

//#include "st_datatypes.h"
#include "st_dsp.h"
#include "st_io.h"

#define MAXCOLUMNS (20)

using namespace std;

void displayUsage(char *progname) 
{
    fprintf(stderr,"%s [-o -r -h] < inputfile > outputfile\n\n",progname);
    fprintf(stderr,"Reads timestamp,value files and returns a fit\n");
    fprintf(stderr, "$Revision: 78 $\n\n");
    fprintf(stderr,"OPTIONS\n");
    fprintf(stderr,"\t-o order      :   Order of the fit, eg -o2 for quadratic fit\n");
    fprintf(stderr,"\t-s starttime      :   Request timetags to start at this time. eg -s 1.75 will produce timetags centered at 1.75 seconds\n");
    fprintf(stderr,"\t-r output_rate:   Output data rate in Hz\n");
    fprintf(stderr,"\t-v            :   Verbose output (print statistics)\n");
    fprintf(stderr,"\t-h            :   This help menu\n");

    return;
}

#define BLOCK_SIZE (128)

int main(int argc, char **argv) 
{

    int fit_order=2;
    float output_rate_hz=1;
    double eval_time=0;
    double chisq=0,chisq_sum[MAXCOLUMNS];
    double coef[10];
    unsigned long int num_fits[MAXCOLUMNS];
    double start_time=0, end_time=0, start_eval_time;
    double observation_interval=0;
    int c;
    int i[MAXCOLUMNS],j=0;
    int pts_min[MAXCOLUMNS], pts_max[MAXCOLUMNS];
    int verbose=0;
    double *xi[MAXCOLUMNS], *yi[MAXCOLUMNS];
    double fit_value;


    char *av[MAXCOLUMNS];
    int ac=0,ac_old=0;
    double array[MAXCOLUMNS];


    /* --- Scan input arguments --- */

    while ((c = getopt (argc, argv, "hvn:o:r:s:")) != -1) {
        switch (c) {
            case 'o':
                sscanf(optarg,"%d",&fit_order);
                break;
            case 'v':
                verbose=1;
                break;
            case 'r':
                sscanf(optarg,"%f",&output_rate_hz);
                break;
            case 's':
                sscanf(optarg,"%lf",&start_eval_time);
                break;
            case 'h':
                displayUsage(argv[0]);
                exit(0);
            default:
                fprintf(stderr,"argc is %d, c is %c\n",argc,c);
                displayUsage(argv[0]);
                exit(0);
        }
    }

    observation_interval=1/output_rate_hz;

    if (verbose) {
        fprintf(stderr,"Fitting Interval : %8.3f [s]\n",observation_interval);
    }

    // For now just allocate all columns
    for(j = 0; j < MAXCOLUMNS; j++) {
        num_fits[j]=0;
        i[j]=0;
        chisq_sum[j]=0;
        pts_max[j]=0;
        pts_min[j]=1e9;
        xi[j]=(double*)malloc(sizeof(double)*BLOCK_SIZE);
        yi[j]=(double*)malloc(sizeof(double)*BLOCK_SIZE);
    }

    char line[1024];

    while(fgets(line,sizeof(line),stdin)!=NULL)
    {
        
        //printf("processing line %d, %f %f\n",i[1],array[0],array[1]);
        // Skip if it's a comment (and print out)
        if (line[0]=='#') {
            fprintf(stdout,"%s",line);
            continue;
        }

        // Figure out how many columns we have (save in ac)
        ac = makeargv(line, av, MAXCOLUMNS);

        // Make sure number of columns didn't change
        if (ac_old>0) {
            assert(ac==ac_old);
        }
        ac_old=ac;


        // Convert columns to float
        for(j = 0; j < ac; j++) {
            array[j] = atof(av[j]);
        }

        // If this is the first time we're computing evaluation time
        // do some special stuff
        if (eval_time==0) {
            // See where we should start taking fitting data.
            // if observation_interval is 1sec, we
            // take data from 0.5 - 1.5 seconds. but if
            // our data starts at, say, 0.6 seconds we skip
            // and start at 1.5 - 2.5 seconds
            // Only compute start time if we're not explicitly given one
            if (start_time==0 && start_eval_time==0) {
                start_time=ceil(array[0]/observation_interval)*
                    observation_interval
                    +observation_interval/2.0;
                start_eval_time=start_time;
            }
            // Else we already gave start time as command line
            // now adjust to be actually at beginning and not
            // user-requested center
            else {
                start_time=start_eval_time-observation_interval/2.0;
            }

            // Here we figure out what the closest n*obs_time
            // is to still get timetags to align properly
            if (array[0]+observation_interval/2.0 > start_time) {
                start_time=start_eval_time+(observation_interval*
                    ceil((array[0]-start_eval_time)/observation_interval)
                    )-observation_interval/2.0;
            }
            eval_time=start_time+observation_interval/2.0;
            end_time=start_time+observation_interval;
            if (verbose) {
                fprintf(stderr,"Start time : %8.6f [s]\n",start_time);
                fprintf(stderr,"End time   : %8.6f [s]\n",end_time);
                fprintf(stderr,"Current time   : %8.6f [s]\n",array[0]);
            }
        }

        /* We have enough points, do the fit */
        if ((array[0]) >= (end_time)) {
            if (verbose) {
                fprintf(stderr,"Enough points for fit\n");
            }
            int print_timetag=0;
            for (j=1; j < ac; j++) {

                // Only fit for phase/range, the rest just compute average
                if (j==9 || j==10) {
                    fit_value=fitter(xi[j],yi[j],i[j],fit_order,
                        &chisq,coef, eval_time,verbose);
                }
                else {
                    // ---- Just compute the average -----
                    fit_value=0;
                    for (int count=0; count < i[j]; count++) {
                        if (j==1 || j==3 || j==5) {
                            fit_value+=fabs(yi[j][count]);
                        }
                        else { 
                            fit_value+=yi[j][count];
                        }

                    }
                    fit_value/=i[j];
                }

                // Throw out fit if bad points
                if (fpclassify(fit_value)==FP_NAN) {
                    // Request to restart things
                    eval_time=0;
                }
                else {
                    num_fits[j]++;
                    if (print_timetag==0) {
                        fprintf(stdout,"%10.6lf",eval_time);
                        print_timetag=1;
                    }
                    if (j==9 || j==10)  {
                        fprintf(stdout," %0.15lf", fit_value);
                    }
                    else  {
                        fprintf(stdout," %0.3lf", fit_value);
                    }
                
                    /* Keep track of min/max points */
                    if (i[j] > pts_max[j]) {
                        pts_max[j]=i[j];
                    }
                    if (i[j] < pts_min[j]) {
                        pts_min[j]=i[j];
                    }
                    chisq_sum[j]+=chisq;
                }
                i[j]=0;
            }

            // Move on if this was a bad fit.
            if (eval_time==0) {
                continue;
            }

            start_time=end_time;
            end_time=start_time+observation_interval;
            fprintf(stdout,"\n");
            eval_time=start_time+observation_interval/2.0;
        }

        /* Store in array for fitting if times match, if
         * not anough memory available, grow buffer in multiples 
         * of BLOCKSIZE*/
        if ((array[0]) >= (start_time) 
            && (array[0]) < (end_time)) {
            for (j=1; j < ac; j++) {
                xi[j][i[j]]=array[0];
                yi[j][i[j]]=array[j];
                i[j]++;
                // we rolled over, reallocate more memory
                if (((i[j]%BLOCK_SIZE)!=((i[j]-1)%BLOCK_SIZE)) && i[j]>0) {
                    xi[j]=(double*)realloc(xi[j],sizeof(double)*
                        ((i[j]/BLOCK_SIZE)+1)*BLOCK_SIZE);
                    yi[j]=(double*)realloc(yi[j],sizeof(double)*
                        ((i[j]/BLOCK_SIZE)+1)*BLOCK_SIZE);
                    assert(xi[j]!=NULL);
                    assert(yi[j]!=NULL);
                }
            }
        }
    }
    
    // For now just allocate all columns
    if (verbose) {
        fprintf(stdout,"# column number, num_fits, pts_min, pts_max, mean chisq\n");
    }
    for(j = 1; j < ac; j++) {
        if (verbose) {
            fprintf(stdout,"# %02d %4ld %3d %3d %g\n",j,
                num_fits[j],pts_min[j],pts_max[j],chisq_sum[j]/num_fits[j]);
        }

        free (xi[j]);
        free (yi[j]);
    }
    
    return 0;
}

