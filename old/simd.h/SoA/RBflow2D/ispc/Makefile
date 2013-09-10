CXX = g++
CXX = $(GCC47)/bin/g++ -O4  -Wall -g  -I$(HOME)/usr/include

LD = $(CXX)

ifeq ($(OMP), 1)
	OMPFLAGS  = -fopenmp -D__OMP__
endif
ifeq ($(OMPGLIBC), 1)
	OMPFLAGS += -D_GLIBCXX_PARALLEL
endif

OFLAGS = -O3 -g -Wall
ifeq ($(GDB), 1)
	OFLAGS = -O0 -g -Wall
endif
	
#OFLAGS  += -funroll-all-loops
CXXFLAGS = $(OFLAGS) -Wstrict-aliasing=2 $(OMPFLAGS)

ifeq ($(AVX), 1)
	CXXFLAGS += -mavx -mtune=corei7-avx -march=corei7-avx
else
	CXXFLAGS += -msse4
endif

ifeq ($(RELEASE), 1)
	CXXFLAGS += -DNDEBUG
endif

ifeq ($(SINGLE),1)
	CXX += -DSINGLE
endif

LDFLAGS = 

SRCPATH = ./
SRC = Integrate.cpp rtc.cpp
OBJ = $(SRC:%.cpp=%.o) 

LIB = rbflow

RM = /bin/rm

all: $(LIB) 

$(LIB): $(OBJ) 
	/bin/rm -rf lib$@.a
	/bin/rm -rf tmp_iscp
	mkdir tmp_iscp
	cd tmp_iscp; ar x ../librbflow_ispc.a
	ar qv lib$@.a $^  tmp_iscp/*.o
	/bin/rm -rf tmp_iscp
	ranlib lib$@.a

%.o: $(SRCPATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	/bin/rm -rf *.o lib$(LIB).a

clean_all:
	/bin/rm -rf *.o lib$(LID).a *~

$(OBJ): FArray.h \
 	Integrate.h \
 	divg.h \
 	hdnl.h \
 	invtr2.h \
 	invtr3.h \
 	invtrro.h \
 	invtr_all.h \
 	update_qp.h \
 	phcalcTP.h \
 	prcalc.h \
 	solq2.h \
 	solq3.h \
 	solro.h\
	tridiag.h \
 	updvp.h \
 	solT.h \
	hdnl_all.h \
	sse_fp32.h \
	sse_fp64.h \
	avx_fp32.h \
	avx_fp64.h

