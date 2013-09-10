ISPC_LIB=rbflow_ispc

ISPC_SRC = init.ispc hdnl.ispc invtr.ispc
ISPC_TARGETS=sse4


TASK_CXX=tasksys.cpp
TASK_LIB=-lpthread
TASK_OBJ=objs/tasksys.o

CXX=g++
CXXFLAGS=-Iobjs/ -O3 -m64
LIBS=-lm $(TASK_LIB) -lstdc++
ISPC=ispc -O3 --arch=x86-64 $(ISPC_FLAGS)
#ISPC_OBJS=$(addprefix objs/, $(ISPC_SRC:.ispc=)_ispc.o)
#ISPC_HEADER=objs/$(ISPC_SRC:.ispc=_ispc.h)
ISPC_OBJS=$(ISPC_SRC:%.ispc=objs/%_ispc.o)
ISPC_HEADER=objs/$(ISPC_SRC:.ispc=_ispc.h)

all: dirs $(ISPC_LIB)

# .PHONY: dirs clean

dirs:
	/bin/mkdir -p objs/

objs/%.cpp objs/%.o objs/%.h: dirs

clean:
	/bin/rm -rf ./objs *~ lib$(ISPC_LIB).a 

$(ISPC_LIB): $(ISPC_OBJS) $(TASK_OBJ)
	/bin/rm -rf lib$@.a
	ar qv lib$@.a $^ 
	ranlib lib$@.a

$(TASK_OBJ): tasksys.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

objs/%.o: %.cpp dirs $(ISPC_HEADER)
	$(CXX) $< $(CXXFLAGS) -c -o $@

objs/%.o: ../%.cpp dirs
	$(CXX) $< $(CXXFLAGS) -c -o $@

objs/$(ISPC_LIB).o: objs/%_ispc.h

objs/%_ispc.h objs/%_ispc.o: %.ispc
	$(ISPC) --target=$(ISPC_TARGETS) $< -o objs/$*_ispc.o -h objs/$*_ispc.h
