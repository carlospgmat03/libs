#CXXFLAGS=-O
include ~/makefile

documentation ::
	doxygen 
register :: register.cpp quantum_state.hpp my_random.hpp
	$(CXX) $(CXXFLAGS) -o $@ $@.cpp

test_memory_time:: test_memory_time.cpp
	$(CXX)  -Wall $(PATHTLCLAP)  $(CXXFLAGS) -o $@ $@.cpp -litpp

# test_rmt:: test_rmt.cpp
# 	$(CXX)  -Wall $(INCLUDETCLAP) $(CXXFLAGS) -o $@ $@.cpp -litpp

test_spins:: $(MYCPP)/test_spins.cpp $(MYCPP)/spinchain.cpp $(MYCPP)/itpp_ext_math.cpp
	$(CXX)  -Wall -I . $(INCLUDETCLAP) $(CXXFLAGS) -o $@ $@.cpp -litpp

test_dev_random_kanbalam :: test_dev_random
	bsub  -oo salida -eo error -q pruebas -n 5 srun ./test_dev_random

testing :: testing.cpp ../cpp/cfp_math.cpp $(MYCPP)/cpp/itpp_ext_math.cpp
	$(CXX)  -o $@ -I $(MYCPP) -litpp $@.cpp

% :: %.cpp  $(MYCPP)/itpp_ext_math.cpp
	g++ $< -o $@ -I $(MYCPP) -litpp

test_rmt:: test_rmt.cpp RMT.cpp

basic_libs.tgz :: cfp_math.cpp cfp_math.h dev_random.cpp itpp_ext_math.cpp itpp_ext_math.h purity_RMT.cpp purity_RMT.h RMT.cpp RMT.h spinchain.cpp 
	tar -cvzf $@ $^
