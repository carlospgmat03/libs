// {{{ Include and TCLAP
#include <iostream>
#include <tclap/CmdLine.h>
#include "cfp_math.cpp"
#include "itpp_ext_math.cpp"
#include "RMT.cpp"
#include "purity_RMT.cpp"
#include "dev_random.cpp"

using namespace std;
using namespace itpp;
using namespace itppextmath;
using namespace RMT;
using namespace cfpmath;
//typedef double wp;
//typedef long double wp;

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
// }}}
int main(int argc, char* argv[]) { // {{{
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
    (option == "get_gse_spectrum"){ // {{{
      vec eig=eig_sym(RandomGSEDeltaOne(dimension.getValue()));
      for (int i=0; i<eig.size(); i++){  cout << eig(i) << endl;}
//       cout << eig << endl;
      return 0;
    } // }}}
  else if (option == "get_gse_flat_spectrum"){ // {{{
      vec eig=FlatSpectrumGSE(dimension.getValue());
      for (int i=0; i<eig.size(); i++){  cout << eig(i) << endl;}
//     cout << eig <<endl;
    return 0;
  } // }}}
  else if (option == "get_gse_flat_spectrum"){ // {{{
      vec eig=FlatSpectrumGSE(dimension.getValue());
      for (int i=0; i<eig.size(); i++){  cout << eig(i) << endl;}
//     cout << eig <<endl;
    return 0;
  } // }}}
  //}}}
} // }}}


