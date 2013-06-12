 // {{{ includes,  namespaces, TCLAP (command line options) and OpenMP (parallel processing).
// Compilation: g++ -o Lambda Lambda.cpp -litpp -fopenmp
#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
#include "cfp_math.cpp"
#include "itpp_ext_math.cpp"
#include "RMT.cpp"
#include <tclap/CmdLine.h>
#include "dev_random.cpp"
#include <cstdio>
#include <ctime>

#define N 50000

using namespace std;
using namespace itpp;
using namespace itppextmath;
using namespace cfpmath;

  TCLAP::CmdLine cmd("Command description message", ' ', "0.1");
  TCLAP::ValueArg<string> optionArg("o","option", "Option"
      ,false,"nichts", "string",cmd);
  TCLAP::ValueArg<unsigned int> seed("s","seed",
      "Random seed [0 for urandom]",false, 243243,"unsigned int",cmd);
  // "Random seed [0 for urandom]",false, 938475975,"unsigned int",cmd);
  TCLAP::ValueArg<int> dimension("n","dimension",
      "Dimension of the system",false, 4,"int",cmd);
  TCLAP::ValueArg<int> time_steps("t","time_steps",
      "Number of time steps",false, 20,"int",cmd);
  TCLAP::ValueArg<double> perturbation("","perturbation",
      "Scalar epsilon that gives the perturbation in H= epsilon*H_0 + V",false, 0.01,"double",cmd);
  TCLAP::ValueArg<double> total_time("","total_time",
      "Total time of evolution",false, 10.,"double",cmd);
  //           "Number of particles per dimension",false, 9,"int",cmd);
  TCLAP::SwitchArg no_general_report("","no_general_report",
      "Print the general report", cmd);
// }}}
void test_fast_LambdaMatrix(){ // {{{
    int n = dimension.getValue();
    double t=10.;
//     test_UDiagUdagger();
//     test_trAB();
//     test_multiply_by_sigma_leftmost_qubit();
//     abort();
    cmat H = RMT::RandomGUE(n)/sqrt(double(n)),W, U;
    vec eigenvalues;
//     LambdaMatrix(H, t);
    eig_sym(H,eigenvalues,W); // Calculate the eigenvalues and eigenvectors
    cout << "Valor " << norm(LambdaMatrixQubitLeft(H, t) -LambdaMatrixQubitLeft(W, eigenvalues, t))<< endl;
  return;
} //}}}
int main(int argc, char* argv[]){ // {{{
  // {{{ initial definitions
  freopen ("output.out", "w", stdout);//write to file with standard I/O commands
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
    (option == "get_all_matrix_elements"){ // {{{
      //     test_fast_LambdaMatrix();
      //     abort(); 
      //
      // #linea de comando: ./Lambda -o get_all_matrix_elements 
      // #semilla = 243243
      // #hostname: papa
      // #comenzando en: Tue May 24 11:25:44 CDT 2011
      // #uname -a: Linux papa 2.6.34-gentoo-r6 #28 SMP Tue Mar 8 12:13:06 CST 2011 i686 Intel(R) Core(TM) i5 CPU 750 @ 2.67GHz GenuineIntel GNU/Linux
      // #working dir: /home/carlosp/investigacion/nacho
      // [[-0.025112053964205398 0.021834274403129628 0.057318461906268249]
      //  [-0.03851518928118873 0.018383377993690227 0.051855234390517896]
      //  [0.1378363982382266 -0.083457252349165484 -0.12048641972925504]]
      // #terminando: Tue May 24 11:25:44 CDT 2011
      //
      //
      //
      int n = dimension.getValue();
      double t=10.;
      cmat H = RMT::RandomGUE(n)/sqrt(double(n)),W;
      mat l;
      vec eigenvalues;
      eig_sym(H,eigenvalues,W); // Calculate the eigenvalues and eigenvectors
      for (int i=0; i<=time_steps.getValue(); i++){
        t=i*total_time.getValue()/time_steps.getValue();
        l = LambdaMatrixQubitLeft(W, eigenvalues, t) ;
        cout << t << " " ;
        for (int col=0; col<l.cols(); col++){ for (int row=0; row<l.rows(); row++){
          cout << l(row,col) << " " ;
        }}
        cout << endl;
      }
    // }}}
   } else if  (option == "get_all_matrix_elements_perturbation"){ // {{{
      //
      // #linea de comando: ./Lambda -o get_all_matrix_elements _perturbation
      //
	int epsilon;
	double max;
	//double mark_med[101];
	#pragma omp parallel for ordered schedule(dynamic) //Option for parallel processing. Optimized for 'for' loops. Each core is doing a different instance of the loop.
	for(epsilon = 0; epsilon < 60; epsilon +=2){
		int n = dimension.getValue();
		double t=10., ep=double(epsilon)/(n*10);;
		int cont;		
		cmat H;
		cmat W;
		mat l;
		vec eigenvalues, alpha=zeros(time_steps.getValue() + 1);
		time_t start, end;

		for(cont = 0; cont < N; cont++){
			H = ep*kron(eye(2),diag(RMT::RandomPUEspectrum(n/2))) + (RMT::RandomPUE(n)/sqrt(double(n))); //Calculates N random hamiltonians
			eig_sym(H,eigenvalues,W); // Calculates the eigenvalues and eigenvectors (complexity O(n^3))
			for (int i=0; i<=time_steps.getValue(); i++){
				t=i*total_time.getValue()/time_steps.getValue();
				l = LambdaMatrixQubitLeft(W, eigenvalues, t); //Calculates the components of the quantum chanel at time t (complexity O(n^3))
				l /= N; //Calculate average components
				alpha(i) = alpha(i) + (l(0,0) + l(1,1) + l(2,2))/3;
			}
		}

		#pragma omp ordered //This line is necessary to get the output in the correct order
		{

		for(int i = 0; i <= time_steps.getValue(); i++)
			printf("%lf ", alpha(i));
		cout << endl << endl;
		}
	}
  } else if (option == "nichts") { // {{{
  // }}}
  } else { // {{{
    cout << "Error en la opcion. Mira, esto es lo que paso: "
      << optionArg.getValue() << endl;
  } // }}}
// }}}
  // {{{ Final report
  if(!no_general_report.getValue()){
    error += system("echo \\#terminando:    $(date)");
  }
  // }}}
  return 0;
} // }}}
