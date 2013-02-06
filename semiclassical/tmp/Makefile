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

onestateeco.out : onestateeco.cc 
	$(GCC) onestateeco.cc -o onestateeco.out $(CFLAGS)  -litpp

clean: 
	rm -f *.o 
	rm -f *.out
	rm -f *~


