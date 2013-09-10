
ifeq ($(OSX), 1) 
CXX=g++-mp-4.7
CXXFLAGS=-Iobjs/ -O3 -m64 -mavx -I../../../ -Wall
LDFLAGS= -lstdc++
CPP_S=$(addprefix objs/, $(CPP_SRC:.cpp=.s))
LD = clang
else
CXX=g++
CXXFLAGS=-Iobjs/ -O3 -m64 -mavx -I../../../ -Wall
CPP_S=$(addprefix objs/, $(CPP_SRC:.cpp=.s))
LDFLAGS=
LD=$(CXX)
endif

ifeq ($(OMP),1)
	CXXFLAGS += -fopenmp
	LDFLAGS += -fopenmp
ifeq ($(OSX),1)
  LDFLAGS += -L/opt/local/lib/gcc47 -lgomp	
endif
endif

default: $(EXAMPLE)

all: $(EXAMPLE)

.PHONY: dirs clean

dirs:
	/bin/mkdir -p objs/

objs/%.cpp objs/%.s: dirs

clean:
	/bin/rm -rf objs *~ $(EXAMPLE) 

$(EXAMPLE): $(CPP_S) 
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

objs/%.s: %.cpp dirs 
	$(CXX) $< $(CXXFLAGS) -S -o $@

objs/%.s: ../%.cpp dirs
	$(CXX) $< $(CXXFLAGS) -S -o $@



