include ~/makefile

#$(OPTIONSONE)



varios.o :: varios.f90 dpt/dpt.o math.o
	$(COMPILER) $(NOLINK) $(OPTFLAGS) -Idpt varios.f90

variosNI.o :: varios.f90 dpf/dpfNI.o mathNI.o
	$(COMPILER) -o variosNI.o $(NOLINK) -Idpf varios.f90


modules.tgz :: 
	tar -cvzf modules.tgz *90 *.f makefile

new_householder.o : new_householder.f
	$(COMPILERNEWHOUSEHOLDER) new_householder.f
new_householderNI.o : new_householder.f
	$(COMPILER) $(NOLINK) -o new_householderNI.o new_householder.f

qr.o     : qr.f
	$(COMPILERQR) qr.f
qrNI.o     : qr.f
	$(COMPILER) $(NOLINK) -o qrNI.o qr.f

inverse_iter.o: inverse_iter.f
	$(COMPILERINVERSEITER) inverse_iter.f
inverse_iterNI.o: inverse_iter.f
	$(COMPILER) $(NOLINK) -o inverse_iterNI.o inverse_iter.f

back_transformation.o: back_transformation.f
	$(COMPILERBACKTRANSFORMATION) back_transformation.f
back_transformationNI.o: back_transformation.f
	$(COMPILER) $(NOLINK) -o back_transformationNI.o back_transformation.f

vhqrii_modules.o: vhqrii_modules.f
	   $(COMPILERVHQRIIMODULES) vhqrii_modules.f
vhqrii_modulesNI.o: vhqrii_modules.f
	$(COMPILER) $(NOLINK) -o vhqrii_modulesNI.o vhqrii_modules.f

all_nonmkl :: new_householder.o qr.o inverse_iter.o back_transformation.o vhqrii_modules.o quantumF90NI.o  quantumF90.o quantum.o diagonalize.o newevolutor.o definitions1.o heap.o the_evolutors.o the_evolutorsNI.o math.o mathNI.o chain_symmetries.o chain_symmetriesNI.o 


all_mkl :: rmt_generator.o rmt_generatorNI.o linear.o linearNI.o 
all :: all_nonmkl all_mkl


advanced_math.o :: advanced_math.f90 math.o
advanced_mathNI.o :: advanced_math.f90 mathNI.o
chain_symmetries.o ::. chain_symmetries.f90 math.o
chain_symmetriesNI.o ::. chain_symmetries.f90 mathNI.o
concurrence.o :: math.o concurrence.f90
concurrenceNI.o :: math.o concurrence.f90
diagonalize.o :: diagonalize.f90 math.o
diagonalizeNI.o :: diagonalize.f90 mathNI.o
exact_ps.o :: exact_ps.f90 ps_rmt.o
exact_psNI.o :: exact_ps.f90 ps_rmtNI.o
linear.o :: linear.f90 math.o
linearNI.o :: linear.f90 mathNI.o
newevolutor.o :: newevolutor.f90 definitions1.o math.o
newevolutorNI.o :: newevolutor.f90 definitions1.o math.o
quantum.o :: quantum.f90 math.o definitions1.o diagonalize.o newevolutor.o
quantumNI.o :: quantum.f90 mathNI.o definitions1NI.o diagonalizeNI.o newevolutorNI.o
quantumF90.o :: quantumF90.f90 math.o
quantumF90NI.o :: quantumF90.f90 mathNI.o 
quantumNI.o :: quantum.f90 mathNI.o definitions1NI.o diagonalizeNI.o newevolutorNI.o
rmt_generator.o :: rmt_generator.f90 math.o linear.o heap.o
rmt_generatorNI.o :: rmt_generator.f90 mathNI.o linearNI.o heapNI.o
the_evolutors.o :: the_evolutors.f90  math.o
the_evolutorsNI.o :: the_evolutors.f90 mathNI.o


%.o :: %.f
	$(FC) $(NOLINK) $(FCOPTIMIZE) $(MKLMODULES) $<
%NI.o :: %.f 
	$(FC) $(NOLINK) -o $@  $(MKLMODULES) $<

%.o :: %.f90 
	$(FC) $(NOLINK) $(FCOPTIMIZE) $(MKLMODULES) $<



%NI.o :: %.f90 
	$(FC) $(NOLINK) $(FCWARNINGS) -o $@  $(MKLMODULES) $<

dpt/dpt.o :: dpt/dpt.f90
	cd dpt; $(FC) $(NOLINK) $(FCOPTIMIZE) dpt.f90
dpf/dpf.o :: dpf/dpf.f90
	cd dpf; $(FC) $(NOLINK) $(FCOPTIMIZE) dpf.f90
dpt/dptNI.o :: dpt/dpt.f90
	cd dpt; $(FC) -o dptNI.o $(NOLINK) dpt.f90
dpf/dpfNI.o :: dpf/dpf.f90
	cd dpf; $(FC) -o dpfNI.o $(NOLINK) dpf.f90

# build_rmt_hams.o :: build_rmt_hams.f90
# 	$(COMPILER) $(NOLINK) $(OPTFLAGS) -I$(BLASNFMKL) -I$(LAPACKNFMKL) build_rmt_hams.f90
# build_rmt_hamsNI.o :: build_rmt_hams.f90
# 	$(COMPILER) -o build_rmt_hamsNI.o $(NOLINK) -I$(BLASNFMKL) -I$(LAPACKNFMKL) build_rmt_hams.f90
# build_rmt_hams.o build_rmt_hamsNI.o 
# 	build_rmt_hams.o build_rmt_hamsNI.o 


all_modules :: math.o mathNI.o linear.o linearNI.o heap.o heapNI.o  \
	rmt_generator.o rmt_generatorNI.o ps_rmt.o ps_rmtNI.o \
	filehandling.o filehandlingNI.o exact_ps.o exact_psNI.o \
	concurrence.o concurrenceNI.o  the_evolutors.o the_evolutorsNI.o \
	advanced_math.o advanced_mathNI.o chain_symmetries.o chain_symmetriesNI.o \
	quantumF90NI.o quantumF90.o quantumNI.o quantum.o \
	newevolutorNI.o newevolutor.o \
 	diagonalizeNI.o diagonalize.o \
	definitions1NI.o definitions1.o 
	
# new_householder.o : new_householder.f
# new_householderNI.o : new_householder.f
# qr.o     : qr.f
# qrNI.o     : qr.f
# inverse_iter.o: inverse_iter.f
# inverse_iterNI.o: inverse_iter.f
# back_transformation.o: back_transformation.f
# back_transformationNI.o: back_transformation.f
# vhqrii_modules.o: vhqrii_modules.f
# vhqrii_modulesNI.o: vhqrii_modules.f
# dpt/dpt.o :: dpt/dpt.f90
# dpf/dpf.o :: dpf/dpf.f90
# dpt/dptNI.o :: dpt/dpt.f90
# dpf/dpfNI.o :: dpf/dpf.f90
# varios.o :: varios.f90 dpt/dpt.o math.o
# variosNI.o :: varios.f90 dpf/dpfNI.o mathNI.o

