
#INSTRUCCION = ifort -c -fast
INSTRUCCION = ifort -c -O3
INSTRUCCION = ifort -c 
COMPILER=ifort
LIBRARIES = -llapack95 -llapack_i386 -lblas_i386
MODULESDIR= -I/home/users/carlosp/lib/lapack95_modules


rmt_generator.o :: rmt_generator.f90 math.o linear.o heap.o
	$(COMPILER) -ipo -c -O3 -I/opt/intel/mkl/8.1.1/interfaces/blas95/ $(MODULESDIR) rmt_generator.f90


linear.o :: math.o linear.f90
	$(COMPILER) -ipo -c -O3 $(MODULESDIR) linear.f90 

heap.o :: heap.f90
	$(COMPILER) -ipo -c -O3 heap.f90

clean   ::
	rm -f *dvi  *aux *log *~ *.o *bbl *blg core.* *.out *.toc *.bm *.end *.fgx *.tbx *.mod

the_evolutors.o :: the_evolutors.f90 makefile bit_manipulation.mod
	$(COMPILER) -ipo -c -O3 the_evolutors.f90

bit_manipulation.mod :: math.o

math.o :: math.f90
	$(COMPILER) -ipo -c -ipo -O3 math.f90

chain_symmetries.o ::. chain_symmetries.f90 math.mod
	$(COMPILER) -ipo -c -O3 chain_symmetries.f90

math.mod : math.o

quantumF90.o :: quantum.f90 
	$(COMPILER) -ipo -c -O3 quantumF90.f90

quantum.o :: quantum.f90
	$(INSTRUCCION) -ipo definitions1.f90 math.f90  diagonalize.f90 newevolutor.f90 quantum.f90


diagonalize.o :: diagonalize.f90
	$(INSTRUCCION) -ipo math.f90 diagonalize.f90

newevolutor.o :: newevolutor.f90 definitions1.f90 math.f90
	$(INSTRUCCION) -ipo definitions1.f90 math.f90 newevolutor.f90


definitions1.o :: definitions1.f90 
	$(INSTRUCCION) -ipo definitions1.f90 

libcuantica.a :: math.o definitions1.o diagonalize.o newevolutor.o quantum.o
	ar rcs libcuantica.a math.o definitions1.o diagonalize.o newevolutor.o quantum.o

all :: libcuantica.a


modules.tgz :: 
	tar -cvzf modules.tgz *90 *.f makefile




new_testhqrii.o: new_testhqrii.f
	   ifort -ipo -c -O1 -tpp7 -pad new_testhqrii.f

aleat.o: aleat.f
	   ifort -ipo -c -O3 -tpp7 aleat.f

new_householder.o: new_householder.f
	   ifort -ipo -c -O2 -tpp7 -axP new_householder.f

qr.o     : qr.f
	   ifort -ipo -c -O3 -tpp7 -axP -fp -pad qr.f

inverse_iter.o: inverse_iter.f
	   ifort -ipo -c -O3 -tpp7 -axP -fp -pad inverse_iter.f

back_transformation.o: back_transformation.f
	   ifort -ipo -c -O2 -tpp7 -axP  -pad -align back_transformation.f

vhqrii_modules.o: vhqrii_modules.f
	   ifort -ipo -c -O3 -tpp7 -axP -ip -pad vhqrii_modules.f
