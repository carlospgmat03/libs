#CXXFLAGS=-O
CXXFLAGS=`itpp-config --cflags`
LIBS=`itpp-config --static --libs`
-include ~/makefile
# Configuration details. Adjust to your needs with
# LOGNAME and, for example, ifeq ($(HOSTNAME),jesus)
# LDLIBS = -litpp
NVCCFLAGS= -arch=sm_13
INCLUDES := -I. -I ../

documentation ::
	doxygen
register :: register.cpp quantum_state.hpp my_random.hpp
	$(CXX) $(CXXFLAGS) -o $@ $@.cpp

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

test_spins:: test_spins.cpp
	$(CXX)  -Wall -I../ -I . $(INCLUDETCLAP) $(CXXFLAGS) -o $@ $@.cpp $(LIBS)

test_dev_random_kanbalam :: test_dev_random
	bsub  -oo salida -eo error -q pruebas -n 5 srun ./test_dev_random

testing :: testing.cpp ../cpp/cfp_math.cpp $(MYCPP)/cpp/itpp_ext_math.cpp
	$(CXX)  -o $@ -I $(MYCPP) -litpp $@.cpp

% :: %.cpp
	$(CXX) $< -Wall -o $@ -I $(MYCPP) -litpp

test_rmt:: test_rmt.cpp RMT.cpp

basic_libs.tgz :: cfp_math.cpp cfp_math.h dev_random.cpp itpp_ext_math.cpp itpp_ext_math.h purity_RMT.cpp purity_RMT.h RMT.cpp RMT.h spinchain.cpp
	tar -cvzf $@ $^
