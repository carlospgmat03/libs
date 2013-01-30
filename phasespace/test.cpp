//  Para meter el estado en cuestion usamos
//  cut -f3- prueba_estado_entrada.txt| ./test --i1 50 -o test_wigner_transformation
//
// {{{ Include and TCLAP
#include <iostream>
#include <tclap/CmdLine.h>
#include <itpp/itbase.h>
// #include <RMT.cpp>
// #include <purity_RMT.cpp>
#include <phasespace.cpp>
using namespace std;
using namespace itpp;
// using namespace RMT;

//typedef double wp;
//typedef long double wp;
#include <dev_random.cpp>
//#include <my_random.hpp>

//using namespace random_variables;
TCLAP::CmdLine cmd("Command description message", ' ', "0.1");
TCLAP::ValueArg<string> optionArg("o","option", "Option"
    ,false,"nichts", "string",cmd);
TCLAP::ValueArg<unsigned int> seed("s","seed",
    "Random seed [0 for urandom]",false, 243243,"unsigned int",cmd);
// "Random seed [0 for urandom]",false, 938475975,"unsigned int",cmd);
TCLAP::ValueArg<int> dimension("n","dimension",
    "Dimension of the system",false, 4,"int",cmd);
TCLAP::SwitchArg no_general_report("","no_general_report",
    "Print the general report", cmd);
  TCLAP::ValueArg<int> i1("","i1", "First integer argument",false, 1,"int",cmd);
// }}}
std::complex<double>     Im=std::complex<double>(0,1);

int main(int argc, char* argv[]) { // {{{
  cout.precision(17); cerr.precision(17);
  // {{{ initial definitions
  cmd.parse( argc, argv );
  int error=0;
  string option=optionArg.getValue();
  cout.precision(17); cerr.precision(17);
  // }}}
  // {{{ Set seed for random
  unsigned int semilla=seed.getValue();
  if (semilla == 0){
    Random semilla_uran; semilla=semilla_uran.strong();
  } 
  RNG_reset(semilla);
  // }}}
  // {{{ Report on the screen
  if(!no_general_report.getValue()){
    cout << "#linea de comando: "; 
    for(int i=0;i<argc;i++){ 
      cout <<argv[i]<<" " ;
    } cout << endl ;
    cout << "#semilla = " << semilla << endl; 
    error += system("echo \\#hostname: $(hostname)");
    error += system("echo \\#comenzando en: $(date)");
    error += system("echo \\#uname -a: $(uname -a)");
    error += system("echo \\#working dir: $(pwd)");
  }
  // }}}
  if //{{{ option == loquesea
    (option == "test_wigner_transformation"){ // {{{
      int d=i1.getValue();
      cvec state(d);
      double x,y;
      for (int i=0; i<d; i++){
        cin >> x  >> y ;
        state(i)=x+std::complex<double>(0,1)*y;
        //         cout << rho(i,j) << endl;
        //       cout << "moco, i=" << i << endl;
      }
      mat w=phasespace::wigner(state);
//       cout << w(0,0) << " " << w(0,1) << endl;
      for (int i=0; i<2*d; i++){
        for (int j=0; j<2*d; j++){
          cout << w(i,j) << endl;
        }
      }
      return 0;
    } // }}}
  //}}}
} // }}}


