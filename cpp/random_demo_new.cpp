#include <iostream>
using namespace std;

//typedef double wp;
typedef long double wp;
#include <my_random.hpp>
//#include <my_random.hpp>

using namespace random_variables;

int main()
{
  // global_random_number_generator.seed(42);
//	QuantumState psi(4);
//	psi.print( );

#include <my_random.hpp>
  cout.precision(10);
    cout.width(40);
  for(int i = 0; i < 2; i++)
    cout << SampleNormal (1, 0)<<endl;
  cout << "ahora lo bueno\n";
  for(int i = 0; i < 10; i++)
    cout << ComplexSampleNormal (1, 10) << ComplexSampleNormal () <<endl;
  
  return 0;
}
