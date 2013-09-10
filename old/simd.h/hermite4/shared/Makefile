#this is ispc 1.3.1dev
ISPC=ispc
ICPC=icpc
GXX=g++
KNCH=knc.h
HERMITE4=./hermite4struct.ispc_h

ifeq ($(FLAT),1)
	HERMITE4=./hermite4flat.ispc_h
endif

ifeq ($(FAST),1)
  FASTMATH=--opt=fast-math
endif

all: hermite4_avx hermite4_mic

hermite4_avx: hermite4.cpp $(HERMITE4)
	/bin/mkdir -p objs
	ln -sf $(HERMITE4) ./hermite4.ispc_h
	$(ISPC) -O3 --arch=x86-64 --target=avx hermite4.ispc -o objs/hermite4_avx_ispc.o -h objs/hermite4_ispc.h $(FASTMATH)
	$(ICPC) -O3 -o hermite4_avx hermite4.cpp objs/hermite4_avx_ispc.o -I./objs -xavx -openmp  -no-prec-div -no-prec-sqrt -restrict -g


hermite4_mic: hermite4.cpp $(HERMITE4)
	/bin/mkdir -p objs
	ln -sf $(HERMITE4) ./hermite4.ispc_h
	$(ISPC) -O3 hermite4.ispc --emit-c++ --target=generic-16 -o hermite4_zmm.cpp -h objs/hermite4_ispc.h --c++-include-file=$(KNCH) $(FASTMATH)
	$(ICPC) -O3 -o hermite4_mic hermite4.cpp  hermite4_zmm.cpp -I./objs -mmic -openmp  -no-prec-div -no-prec-sqrt -no-opt-prefetch -restrict -g


clean:
	/bin/rm -rf hermite4_avx hermite4_mic objs hermite4_zmm.cpp hermite4.ispc_h

clean_all: clean
	/bin/rm -rf *~

run:
	GOMP_CPU_AFFINITY=[8-15] OMP_NUM_THREADS=4 ./hermite4_avx
