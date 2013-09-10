#this is ispc 1.3.1dev
ISPC=ispc
ICPC=icpc
GXX=g++
KNCH=knc.h

ifeq ($(PROFILE),1)
	ICPC+=-DPROFILE
endif

all: radix_avx radix_mic


radix_avx: radix.cpp
	$(ICPC) -O3 -o radix_avx radix.cpp -xavx -openmp  -g -ansi-alias -funroll-all-loops

radix_mic: radix.cpp 
	$(ICPC) -O3 -o radix_mic radix.cpp -mmic -openmp  -g -ansi-alias -funroll-all-loops -fargument-noalias 

clean:
	/bin/rm -rf radix_avx radix_mic

clean_all: clean
	/bin/rm -rf *~

run:
