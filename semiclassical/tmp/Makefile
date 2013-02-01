#LIBCFLAGS=-O3 -c -funroll-all-loops
#CFLAGS=-O3 -funroll-all-loops
#LIBCFLAGS=-O3 -c -mcpu=i686 -funroll-all-loops
#CFLAGS=-O3 -lm -mcpu=i686 -funroll-all-loops
#LIBCFLAGS=-O3 -c -mcpu=i386 
#CFLAGS=-O3 -lm -mcpu=i386 
#GCC=/usr/bin/gcc

LIBCFLAGS=-O3 -c 
#-m32  
#CFLAGS=-O3 -lm -I/usr/local/fftw/include -L/usr/local/fftw/lib -lfftw3
CFLAGS=-O3 -lm -lfftw3
#-m32	
GCC=g++ -w

ejemplo ::
	make onestateeco.out
	./onestateeco.out 1. 0.4 0 150 0.3 0.4 1 3 2.

##los programas que usan qstate en lugar de quantum, no estan limitados por
##la dimension, i.e. no son "estados de qubits"

all: complex.o qstate.o random.o onestateeco.out strip 
# all: complex.o quantum.o qstate.o random.o onestateeco.out eco2.out eco_pos.out eco_har.out loopitmap.out maps_eco.out maps2eco.out strip 

strip : 
	strip  onestateeco.out 

random.o : random.cc random.h
	$(GCC) $(LIBCFLAGS) random.cc

# quantum.o : quantum.cc quantum.h gcd.h
# 	$(GCC) $(LIBCFLAGS) quantum.cc

qstate.o : qstate.cc qstate.h 
# qstate.o : qstate.cc qstate.h gcd.h
	$(GCC) $(LIBCFLAGS) qstate.cc
	
onestateeco.out : onestateeco.cc qstate.o 
	$(GCC) onestateeco.cc qstate.o -o onestateeco.out $(CFLAGS)  -litpp

# eco2.out : eco2.cc quantum.o complex.o random.o 
# 	$(GCC) eco2.cc quantum.o complex.o random.o -o eco2.out $(CFLAGS)

# eco_pos.out : eco_pos.cc quantum.o complex.o random.o 
# 	$(GCC) eco_pos.cc quantum.o complex.o random.o -o eco_pos.out $(CFLAGS)
# 
# eco_har.out : eco_har.cc quantum.o complex.o random.o 
# 	$(GCC) eco_har.cc quantum.o complex.o random.o -o eco_har.out $(CFLAGS)
# 	
# loopitmap.out : loopitmap.cc quantum.o complex.o random.o 
# 	$(GCC) loopitmap.cc quantum.o complex.o random.o -o loopitmap.out $(CFLAGS)
# 
# maps_eco.out : maps_eco.cc qstate.o complex.o random.o 
# 	$(GCC) maps_eco.cc qstate.o complex.o random.o -o maps_eco.out $(CFLAGS)
# 	
# maps2eco.out : maps2eco.cc qstate.o complex.o random.o 
# 	$(GCC) maps2eco.cc qstate.o complex.o random.o -o maps2eco.out $(CFLAGS)
# 
clean: 
	rm -f *.o 
	rm -f *.out
	rm -f *~


