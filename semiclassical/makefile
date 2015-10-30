include ~/makefile

% :: %.cpp
	$(CXX) -I$(MYCPP) $< -o $@ -litpp



fc=gfortran
FFLAGS=-O3 -msse4  -funroll-loops --param max-unroll-times=4

# wigner: wigner.f90 torus_mac.o plotdistColor.o
wigner: wigner.f90  torus_mac.o
	$(fc) -o wigner wigner.f90 torus_mac.o  $(FFLAGS) -llapack -lm -lfftw3 
torus_mac.o: torus_mac.f90
	$(fc) -c torus_mac.f90 $(FFLAGS) -llapack -lm -lfftw3 

%.o :: %.cpp
	$(CXX) -c -I$(MYCPP) $<  -litpp

test.o :: test.cpp 
	$(CXX) -c -I$(MYCPP) -I$(MYCPP)/semiclassical  -I$(MYCPP)/cpp  -I$(PHASESPACE) $<  -litpp

test :: test.o torus_mac.o
	$(CXX) -I$(MYCPP) -I$(MYCPP)/semiclassical $^ -o $@ -litpp  -llapack -lm -lfftw3


CFLAGS=-O3 -lm -lfftw3
onestateeco.out : onestateeco.cc semiclassical.cpp
	$(CXX) -I . onestateeco.cc -o onestateeco.out $(CFLAGS)  -litpp


