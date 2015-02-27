
CXX=icpc -m64
LD=icpc -m64
CXXFLAGS = -std=c++11  -I../../vmath 
CXXFLAGS+= -O3 -ansi-alias -fp-model fast -inline-forceinline -inline-level=1 
CXXFLAGS+= -xavx

default: mandelbrot 


clean:
	/bin/rm -rf *~ *.o mandelbrot

mandelbrot: mandelbrot.o mandelbrot_serial.o
	$(LD) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

