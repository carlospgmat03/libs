#include <iostream>
#include <tclap/CmdLine.h>
#include "RMT.cpp"
// #include <purity_RMT.cpp>
using namespace std;

//typedef double wp;
//typedef long double wp;
// #include <dev_random.cpp>
//#include <my_random.hpp>

//using namespace random_variables;

int main(int argc, char* argv[])
{
// 	Random semilla_uran;
	int q=13;
	itpp::cmat test_1, test_2;
        test_1 = RMT::RandomGUEDeltaOne(1<<q);
        system("sleep 20");
        
//         test_2 = RMT::RandomGUEDeltaOne(1<<q);
//         test_3 = test_1 * test_2;
//         test_2.set_size(1<<q, 1<<q);
// 	itpp::imat test(50);
//   	cout << PurityRMT::QubitEnvironmentHamiltonian(3,0.) <<endl;
// 	cout << RMT::FlatSpectrumGUE(5,.1) <<endl;
	return 0;
}


