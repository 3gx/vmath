#this is ispc 1.3.1dev
ISPC=ispc
ICPC=icpc
GXX=g++
KNCH=knc.h

ifeq ($(ALIGN),1)
	ALIGNED = --opt=force-aligned-memory
endif 

all: exclusive_scan_avx exclusive_scan_mic

exclusive_scan_avx: exclusive_scan.cpp exclusive_scan.ispc
	/bin/mkdir -p objs
	$(ISPC) -O3 --arch=x86-64 --target=avx exclusive_scan.ispc -o objs/exclusive_scan.o -h objs/exclusive_scan_ispc.h $(ALIGNED)
	$(ICPC) -O3 -o exclusive_scan_avx exclusive_scan.cpp objs/exclusive_scan.o -I./objs -xavx -openmp  -no-prec-div -no-prec-sqrt -g


exclusive_scan_mic: exclusive_scan.cpp exclusive_scan.ispc
	/bin/mkdir -p objs
	$(ISPC) -O3 exclusive_scan.ispc --emit-c++ --target=generic-16 -o exclusive_scan_zmm.cpp -h objs/exclusive_scan.h --c++-include-file=$(KNCH) $(ALIGNED)
	$(ICPC) -O3 -o exclusive_scan_mic exclusive_scan.cpp  exclusive_scan_zmm.cpp -I./objs -mmic -openmp  -no-prec-div -no-prec-sqrt -no-opt-prefetch -restrict -g


clean:
	/bin/rm -rf exclusive_scan_avx exclusive_scan_mic objs exclusive_scan_zmm.cpp 

clean_all: clean
	/bin/rm -rf *~

run:
