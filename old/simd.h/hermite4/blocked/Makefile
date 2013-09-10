#this is ispc 1.3.1dev
ISPC=ispc
ICPC=icpc
KNCH=knc.h

ifeq ($(FAST),1)
  FASTMATH=--opt=fast-math
endif
 
all: hermite4_avx hermite4_mic hermite6_avx hermite6_mic hermite8_avx hermite8_mic

hermite4_avx: compile.cpp  parameter.cpp hermite4.h hermite4.ispc nbodysystem.h
	/bin/mkdir -p objs
	$(ISPC) -O3 --arch=x86-64 --target=avx hermite4.ispc -o objs/hermite4_avx_ispc.o -h objs/hermite4_ispc.h $(FASTMATH)
	$(ICPC) -O3 -o hermite4_avx compile.cpp parameter.cpp objs/hermite4_avx_ispc.o -I./objs -xavx -openmp  -no-prec-div -no-prec-sqrt -g -DFOURTH 


hermite4_mic: compile.cpp  parameter.cpp hermite4.h hermite4.ispc nbodysystem.h
	/bin/mkdir -p objs
	$(ISPC) -O3 hermite4.ispc --emit-c++ --target=generic-16 -o hermite4_zmm.cpp -h objs/hermite4_ispc.h --c++-include-file=$(KNCH) $(FASTMATH) 
	$(ICPC) -O3 -o hermite4_mic compile.cpp parameter.cpp  hermite4_zmm.cpp -I./objs -mmic -openmp  -no-prec-div -no-prec-sqrt -no-opt-prefetch -g -DFOURTH 

hermite6_avx: compile.cpp  parameter.cpp hermite6.h hermite6.ispc nbodysystem.h
	/bin/mkdir -p objs
	$(ISPC) -O3 --arch=x86-64 --target=avx hermite6.ispc -o objs/hermite6_avx_ispc.o -h objs/hermite6_ispc.h $(FASTMATH)
	$(ICPC) -O3 -o hermite6_avx compile.cpp parameter.cpp objs/hermite6_avx_ispc.o -I./objs -xavx -openmp  -no-prec-div -no-prec-sqrt -g -DSIXTH

hermite6_mic: compile.cpp  parameter.cpp hermite6.h hermite6.ispc nbodysystem.h
	/bin/mkdir -p objs
	$(ISPC) -O3 hermite6.ispc --emit-c++ --target=generic-16 -o hermite6_zmm.cpp -h objs/hermite6_ispc.h --c++-include-file=$(KNCH) $(FASTMATH) 
	$(ICPC) -O3 -o hermite6_mic compile.cpp parameter.cpp  hermite6_zmm.cpp -I./objs -mmic -openmp  -no-prec-div -no-prec-sqrt -no-opt-prefetch -g -DSIXTH

hermite8_avx: compile.cpp  parameter.cpp hermite8.h hermite8.ispc nbodysystem.h
	/bin/mkdir -p objs
	$(ISPC) -O3 --arch=x86-64 --target=avx hermite8.ispc -o objs/hermite8_avx_ispc.o -h objs/hermite8_ispc.h $(FASTMATH)
	$(ICPC) -O3 -o hermite8_avx compile.cpp parameter.cpp objs/hermite8_avx_ispc.o -I./objs -xavx -openmp  -no-prec-div -no-prec-sqrt -g -DEIGHTH

hermite8_mic: compile.cpp  parameter.cpp hermite8.h hermite8.ispc nbodysystem.h
	/bin/mkdir -p objs
	$(ISPC) -O3 hermite8.ispc --emit-c++ --target=generic-16 -o hermite8_zmm.cpp -h objs/hermite8_ispc.h --c++-include-file=$(KNCH) $(FASTMATH) 
	$(ICPC) -O3 -o hermite8_mic compile.cpp parameter.cpp  hermite8_zmm.cpp -I./objs -mmic -openmp  -no-prec-div -no-prec-sqrt -no-opt-prefetch -g -DEIGHTH


clean:
	/bin/rm -rf hermite4_avx hermite4_mic objs hermite4_zmm.cpp hermite6_avx hermite6_mic hermite6_zmm.cpp hermite8_avx hermite8_mic hermite8_zmm.cpp

clean_all: clean
	/bin/rm -rf *~

