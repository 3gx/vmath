CXX = g++  
# CXX = $(GCC47)/bin/g++ -std=gnu++11 -O4  -Wall -g  -I$(HOME)/usr/include  

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
	ar qv lib$@.a $^
	ranlib lib$@.a

%.o: $(SRCPATH)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	/bin/rm -rf *.o lib$(LIB).a

clean_all:
	/bin/rm -rf *.o lib$(LID).a *~

$(OBJ): FArray.h \
 	Integrate.h \
	tridiag.h \
 	solT.h \
 	hdnl.h \
 	invtr.h \
 	solq2.h \
 	solq3.h \
 	solro.h\
 	divg.h \
 	phcalcTP.h \
 	update_qp.h \
	AVX/avx.h AVX/avx_fp32.h AVX/avx_fp64.h \
	SSE/sse.h SSE/sse_fp32.h SSE/sse_fp64.h \
	../dims.h

