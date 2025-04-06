

# Set up the environment
# for the tadder 750 Cross-Compiler
ifndef  environmentComplete
ifdef   TRIG_INSTRUMENT
ifdef   ELDK_5_1
	GCC                   = /opt/eldk-5.1/powerpc/sysroots/i686-eldk-linux/usr/bin/powerpc-linux/powerpc-linux-gcc
	GPP                   = /opt/eldk-5.1/powerpc/sysroots/i686-eldk-linux/usr/bin/powerpc-linux/powerpc-linux-g++
	LD					  = /opt/eldk-5.1/powerpc/sysroots/i686-eldk-linux/usr/bin/powerpc-linux/powerpc-linux-ld
	environmentString     = "TRIG_INSTRUMENT ELDK_5_1 750 Cross-Compiler Environment"
	INCLUDE				  =-I. -I../fftw-3.3.3/api/ -L../fftw-3.3.3/.libs -I../armadillo-3.800.2/include -L../armadillo-3.800.2/lib -I../gsl-1.15 -L../gsl-1.15/.libs -L../gsl-1.15/cblas/.libs
else
	GCC                   = /opt/ppc_6xx/usr/bin/ppc_7xx-gcc
	GPP                   = /opt/ppc_6xx/usr/bin/ppc_7xx-g++
	LD					  =/opt/ppc_6xx/usr/bin/ppc_7xx-ld
	environmentString     = "TRIG_INSTRUMENT ppc_7xx 750 Cross-Compiler Environment"
endif
	GCC_CFLAGS            = -fsigned-char
	GCC_LFLAGS            =
	GPP_CFLAGS            = -fsigned-char
	GPP_LFLAGS            =
	environmentComplete   = yes
endif
endif

# Set up for no environment variable being available
ifndef environmentComplete
	INCLUDE=-I.
	GCC                   = /usr/bin/gcc
	GCC_CFLAGS            =
	GCC_LFLAGS            =
	GPP                   = /usr/bin/g++
	GPP_CFLAGS            =
	GPP_LFLAGS            =
	LD					  = /usr/bin/ld
	environmentComplete = No
	environmentString = "No Environment setup"
endif



CC=$(GPP)
CFLAGS=-g -Wall -O2  -Wno-unused-result
ifdef FITTER
	LIBS=-lm libst.so -lfftw3 -lgsl -lgslcblas
else
	#LIBS=-lm libst.so -lfftw3 
	#LIBS=-lm libst.so -lfftw3 -larmadillo -lgfortran -llapack -lblas -lgsl -lgslcblas
	LIBS=-lm libst.so -lfftw3 -larmadillo -llapack -lblas -lgsl -lgslcblas
endif
CC_OBJS=st_dsp.o st_io.o

#all: fft_acquire tone_track prn_track st_reverse 
#all: fft_acquire tone_track prn_track st_checkaccum fitter
all: fft_acquire tone_track prn_track fitter

libst.so: ${CC_OBJS}
	$(LD) -g -Ur ${CC_OBJS} -o libst.so


$(CC_OBJS) : ${@:.o=.c} ${@:.o=.h}
	    $(CC) $(CFLAGS) $(@:.o=.c) -c $(INCLUDE) $(DEFINES)

fft_acquire: fft_acquire.c st_dsp.o st_io.o
	#$(CC) $(CFLAGS) fft_acquire.c -o fft_acquire st_io.o st_dsp.o -lm -lfftw3   $(INCLUDE)
	$(CC) $(CFLAGS) fft_acquire.c -o fft_acquire $(LIBS)   $(INCLUDE)

fitter: fitter.c libst.so
	$(CC) $(CFLAGS) fitter.c -o fitter $(LIBS) $(INCLUDE)

tone_process: tone_process.c libst.so
	$(CC) $(CFLAGS) tone_process.c -o tone_process $(LIBS) $(INCLUDE)

tone_track: tone_track.c libst.so
	$(CC) $(CFLAGS) tone_track.c -o tone_track $(LIBS) $(INCLUDE)

prn_track: prn_track.c libst.so
	$(CC) $(CFLAGS) prn_track.c -o prn_track $(LIBS) $(INCLUDE)

st_reverse: st_reverse.c libst.so
	$(CC) $(CFLAGS) st_reverse.c -o st_reverse $(LIBS) $(INCLUDE)

st_checkaccum: st_checkaccum.cpp libst.so
	$(GPP) $(CFLAGS) st_checkaccum.cpp -o st_checkaccum $(LIBS) $(INCLUDE)

md: data/header1e.txt /mnt/bigdrive/esterhui/work/grail/prsr/data/20110803/test_prn1e.bin
	cat data/header1e.txt > data/data1e.bin
	cat /mnt/bigdrive/esterhui/work/grail/prsr/data/20110803/test_prn1e.bin >> data/data1e.bin

clean:
	\rm -f *.o fft_acquire tone_track prn_track st_reverse st_checkaccum fitter *.so *.pyc
