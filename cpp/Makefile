#CXXFLAGS=-O
CXXFLAGS=`itpp-config --cflags`
LIBS=`itpp-config --static --libs`
-include ~/makefile
# Configuration details. Adjust to your needs with
# LOGNAME and, for example, ifeq ($(HOSTNAME),jesus)
# LDLIBS = -litpp
NVCCFLAGS= -arch=sm_13
INCLUDES := -I. -I ../

default:: clean test_spins

all:: clean documentation test_spins test_spins_fast test_memory_time test_spins2d test_rmt basic_libs.tgz testing test_dev_random

documentation ::
	@doxygen ./Doxyfile
# register:: register.cpp
# 	$(CXX) -I../ $(CXXFLAGS) -o $@ $@.cpp $(LIBS)

test_memory_time:: test_memory_time.cpp
	$(CXX)  -I../ -Wall $(PATHTLCLAP)  $(CXXFLAGS) -o $@ $@.cpp $(LIBS)

# test_rmt:: test_rmt.cpp
# 	$(CXX)  -Wall $(INCLUDETCLAP) $(CXXFLAGS) -o $@ $@.cpp -litpp

test_spins2d:: test_spins2d.cpp
# Test if nvcc is available if so, then
# 	nvcc -D INCLUDECUDA=TRUE -I . -arch=sm_13 $(INCLUDECUDA) $(INCLUDETCLAP) $(CXXFLAGS) -o $@ $@.cu -litpp
# else use the usual compilation line, with some preprocessing
# 	$(CXX) -I . $(INCLUDETCLAP) $(CXXFLAGS)
	$(CXX) -I../  $(CXXFLAGS) $< -Wall -o $@ $(LIBS)

test_spins_fast :: test_spins.cpp
	$(CXX) -O3  -march=native -fomit-frame-pointer -I../ -I . $(INCLUDETCLAP) $(CXXFLAGS) -o $@ test_spins.cpp $(LIBS)

test_spins:: test_spins.cpp
	$(CXX)  -Wall -I../ -I . $(INCLUDETCLAP) $(CXXFLAGS) -o $@ $@.cpp $(LIBS)

test_dev_random:: test_dev_random.cpp
	$(CXX)  -Wall -I../ -I . $(INCLUDETCLAP) $(CXXFLAGS) -o $@ $@.cpp $(LIBS)

#test_dev_random_kanbalam :: test_dev_random
#	bsub  -oo salida -eo error -q pruebas -n 5 srun ./test_dev_random

testing :: testing.cpp cfp_math.cpp itpp_ext_math.cpp
	$(CXX) -I../ -I . -o $@ $@.cpp $(LIBS)
	#$(CXX)  -o $@ -I $(MYCPP) -litpp $@.cpp

% :: %.cpp
	$(CXX) $< -Wall -o $@ -I $(MYCPP) -litpp

test_rmt:: test_rmt.cpp
	$(CXX)  -I../ -Wall $(PATHTLCLAP)  $(CXXFLAGS) -o $@ $@.cpp $(LIBS)

basic_libs.tgz :: cfp_math.cpp cfp_math.h dev_random.cpp itpp_ext_math.cpp itpp_ext_math.h purity_RMT.cpp purity_RMT.h RMT.cpp RMT.h spinchain.cpp
	tar -cvzf $@ $^

.PHONY: clean

clean::
	rm -f test_spins
	rm -f test_spins2d
	rm -f test_memory_time
	rm -f test_rmt
	rm -f test_spins_fast
	rm -f testing
	rm -f test_dev_random
	rm -f basic_libs.tgz
	rm -rf doc
	echo Clean done
