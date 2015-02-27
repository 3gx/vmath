
CXX=icpc -m64
LD=icpc -xavx -m64
CXXFLAGS = -std=c++11  -I../../vmath 
CXXFLAGS+= -O3
CXXFLAGS+= -xavx

default: simple


clean:
	/bin/rm -rf *~ *.o simple

simple: simple.o
	$(LD) -o $@ $<

simple.o: simple.cpp 
	$(CXX) $(CXXFLAGS) -c -o $@ $<

