// {{{ Include and TCLAP
#include <iostream>
#include <tclap/CmdLine.h>
#include "cfp_math.cpp"
#include "itpp_ext_math.cpp"
#include "spinchain.cpp"
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
TCLAP::ValueArg<string> optionArg("o","option", "Option" ,false,"nichts", "string",cmd);
TCLAP::ValueArg<string> channelArg("","ChannelType", "The Channel type" ,false,"RandomFuzzyMeasurement", "string",cmd);
TCLAP::ValueArg<unsigned int> seed("s","seed",
    "Random seed [0 for urandom]",false, 243243,"unsigned int",cmd);
// "Random seed [0 for urandom]",false, 938475975,"unsigned int",cmd);
TCLAP::ValueArg<int> qubits("q","qubits", "Qubits in the system",false, 2,"int",cmd);
TCLAP::ValueArg<int> sample_size("","sample_size",
    "Sample size",false, 20,"int",cmd);
TCLAP::ValueArg<int> integer1("","i1", "One integer to pass, for example first index of block",false, 0,"int",cmd);
TCLAP::ValueArg<int> integer2("","i2", "Second integer to pass, for example second index of block",false, 1,"int",cmd);
TCLAP::ValueArg<int> integer3("","i3", "Third integer to pass, for example intersection bits",false, 1,"int",cmd);
TCLAP::ValueArg<int> dimension("d","dimension", "Dimension of the system",false, 4,"int",cmd);
TCLAP::ValueArg<double> swap_probability("","swap_probability", "Some probability",false, 0.2,"double",cmd);
TCLAP::ValueArg<double> betaTCLAP("","beta", "Inverse temperature for montecarlo",false, 1.,"double",cmd);
TCLAP::ValueArg<double> delta_stepTCLAP("","delta", "Step size",false, 0.01,"double",cmd);
TCLAP::ValueArg<double> polarization_z("","r_z", "Polarization in z direction",false, 0.2,"double",cmd);
TCLAP::ValueArg<double> errorTCLAP("","max_tolerance", "Maximum error tolerated",false, 0.001,"double",cmd);
TCLAP::SwitchArg no_general_report("","no_general_report","Print the general report",cmd);
// }}}

struct FuzzyElement {  // {{{
    double p; 
    itpp::Array<itpp::ivec> permutation;
}; 
Array<FuzzyElement> ChainFuzzyMeasurement(int q, double p){ //{{{
  Array<FuzzyElement> ChannelChain(q+1);
  ChannelChain(0).p = p;
  ChannelChain(0).permutation.set_size(0);

  for (int i=0; i<q; i++){
    ChannelChain(i+1).p = (1.-p)/q;
    ChannelChain(i+1).permutation.set_size(1);
    ChannelChain(i+1).permutation(0)=vec_2(i,(i+1)%q);
  }
  return ChannelChain;
} // }}}
Array<FuzzyElement> AllToAllFuzzyMeasurement(int q, double p){ //{{{
  Array<FuzzyElement> Channel((q*(q-1)/2));
  Channel(0).p = p;
  Channel(0).permutation.set_size(0);

  double p2=(1.-p)/(Channel.size()-1);
  int k=1;
  for (int i=0; i<q-1; i++){
    for (int j=i+1; j<q; j++){
      Channel(k).p = p2;
      Channel(k).permutation.set_size(1);
      Channel(k).permutation(0)=vec_2(i,j);
    }
  }
  return Channel;
} // }}}
Array<FuzzyElement> TrivialFuzzyMeasurement(int q){ //{{{
  Array<FuzzyElement> ChannelChain(1);
  ChannelChain(0).p = 1;
  ChannelChain(0).permutation.set_size(0);
  return ChannelChain;
} // }}}
Array<FuzzyElement> SingleFullRotationFuzzyMeasurement(int q, double p){ //{{{
  Array<FuzzyElement> Channel(2);
  cout << "Back on black" << endl;
  Channel(0).p = 1-p;
  Channel(0).permutation.set_size(0);
  Channel(1).p = p;
  cout << "Back on black 5" << endl;
  Channel(1).permutation.set_size(1);
  cout << "Back on black 6" << endl;
  Channel(1).permutation(0).set_size(q);
  cout << "Back on black 7" << endl;
  Channel(1).permutation(0) = itpp::linspace_fixed_step(1,q-2);
  cout << "Back on black 9" << endl;
  Channel(1).permutation(0)(q-1) = 0;
  cout << "Back on black 10" << endl;
  return Channel;
} // }}}
Array<FuzzyElement> RandomFuzzyMeasurement(int q, double p){ // {{{
//   cout << factorial(q) << endl;
  Array<FuzzyElement> ChannelChain(factorial(q));
  double total_p=0;
  ivec line_permutation(q);
  ChannelChain(0).p = p;
  for (int i=0; i<q; i++){ line_permutation(i)=i; }
  ChannelChain(0).permutation= ChangeNotationPermutation(line_permutation);
  for (int i=1; i<factorial(q); i++){
    ChannelChain(i).p = itpp::randu();
    total_p += ChannelChain(i).p;
//       cout <<  ChangeNotationPermutation(perline) << endl;
    line_permutation=next_permutation(line_permutation);
    ChannelChain(i).permutation=  ChangeNotationPermutation(line_permutation);
  }
  for (int i=1; i<factorial(q); i++){
    ChannelChain(i).p = ChannelChain(i).p *(1-p)/total_p;
//     total_p += ChannelChain(i).p;
//       cout <<  ChangeNotationPermutation(perline) << endl;
//     line_permutation=next_permutation(line_permutation);
//     ChannelChain(i).permutation=  ChangeNotationPermutation(line_permutation);
  }
  return ChannelChain;
} // }}}
Array<FuzzyElement> RandomFuzzyMeasurementTwo(int q, double p){ // {{{
//   cout << factorial(q) << endl;
  Array<FuzzyElement> ChannelChain((q*(q-1))/2 + 1);
  double total_p=0;
  int k=1;
  ivec line_permutation(q);
  ChannelChain(0).p = p;
  for (int i=0; i<q; i++){ line_permutation(i)=i; }
  ChannelChain(0).permutation= ChangeNotationPermutation(line_permutation);
  for (int i=0; i<q-1; i++){
    for (int j=i+1; j<q; j++){
      ChannelChain(k).p = itpp::randu();
      total_p += ChannelChain(k).p;
      //       cout <<  ChangeNotationPermutation(perline) << endl;
      //     line_permutation=next_permutation(line_permutation);
      ChannelChain(k).permutation= SingleTwoPermutation(i, j);
      k++;
    }
  }
  for (int i=1; i<ChannelChain.size(); i++){
    ChannelChain(i).p = ChannelChain(i).p *(1-p)/total_p;
//     total_p += ChannelChain(i).p;
//       cout <<  ChangeNotationPermutation(perline) << endl;
//     line_permutation=next_permutation(line_permutation);
//     ChannelChain(i).permutation=  ChangeNotationPermutation(line_permutation);
  }
  return ChannelChain;
} // }}}
Array<FuzzyElement> RandomFuzzyMeasurementChain(int q, double p){ // {{{
//   cout << factorial(q) << endl;
  Array<FuzzyElement> ChannelChain(q + 1);
  double total_p=0;
  int k=1;
  ivec line_permutation(q);
  ChannelChain(0).p = p;
  for (int i=0; i<q; i++){ line_permutation(i)=i; }
  ChannelChain(0).permutation= ChangeNotationPermutation(line_permutation);
  for (int i=1; i<q+1; i++){
      ChannelChain(i).p = itpp::randu();
      total_p += ChannelChain(i).p;
      ChannelChain(i).permutation= SingleTwoPermutation(i-1, i%q);
  }
  for (int i=1; i<ChannelChain.size(); i++){
    ChannelChain(i).p = ChannelChain(i).p *(1-p)/total_p;
//     total_p += ChannelChain(i).p;
//       cout <<  ChangeNotationPermutation(perline) << endl;
//     line_permutation=next_permutation(line_permutation);
//     ChannelChain(i).permutation=  ChangeNotationPermutation(line_permutation);
  }
  return ChannelChain;
} // }}}
double total_probability(Array<FuzzyElement>& Channel){ // {{{
  int n=Channel.size();
  double p=0;
  for (int i=0; i<n; i++){
    p+= Channel(i).p;
  }
  return p;
} // }}}
void NormalizeMeasurement(Array<FuzzyElement>& Channel){ // {{{
//   cout << factorial(q) << endl;
  double p0=total_probability(Channel);
  for (int i=0; i<Channel.size(); i++){
    Channel(i).p = Channel(i).p/p0;
  }
  return ;
} // }}}
std::ostream& operator << (std::ostream& o, const FuzzyElement& a){ // {{{
    o << "permutation: " << a.permutation << "\tp: " << a.p ;
    return o;
} // }}}
std::ostream& operator << (std::ostream& o, const itpp::Array<FuzzyElement>& a){ // {{{
  for (int i=0; i<a.size(); i++){
    o << i << "\t: " << a(i) << std::endl;
  }
  return o;
} // }}}
Array<FuzzyElement> ChannelSelect(int q, double p, string ChannelType){ // {{{
//   Array<FuzzyElement> Channel;
  if (ChannelType == "ChainFuzzyMeasurement") return ChainFuzzyMeasurement(q, p);
  if (ChannelType == "AllToAllFuzzyMeasurement") return AllToAllFuzzyMeasurement(q, p);
  if (ChannelType == "TrivialFuzzyMeasurement") return TrivialFuzzyMeasurement(q);
  if (ChannelType == "RandomFuzzyMeasurement") return RandomFuzzyMeasurement(q, p);
  if (ChannelType == "RandomFuzzyMeasurementTwo") return RandomFuzzyMeasurementTwo(q, p);
  if (ChannelType == "RandomFuzzyMeasurementChain") return RandomFuzzyMeasurementChain(q, p);
  if (ChannelType == "Rotation") return SingleFullRotationFuzzyMeasurement(q, p);
  cerr << "No se encontro el tipo en ChannelSelect" << endl;
  abort();
} // }}}
 // }}}

itpp::cmat coarse_graining_two_positions(itpp::cvec state, double swap_probability){ // {{{
        itpp::cvec swaped_state;
        swaped_state = state;
        apply_swap(swaped_state, 0, 1);
        return (1-swap_probability)*partial_trace(state, 2) 
           + (swap_probability)*partial_trace(swaped_state, 2);
} // }}}

void apply_sigma_x_like(itpp::cvec& state, int target_level_0, int target_level_1, double theta){// {{{

  std::complex<double> I(0,1);
  double cos_theta  = cos(-theta);
  double sin_theta  = sin(-theta);

  itpp::cmat m(2,2);
  itpp::cvec mini_state(2);
  m(0,0) = cos_theta;
  m(1,1) = cos_theta;
  m(0,1) = I*sin_theta;
  m(1,0) = I*sin_theta;
  mini_state(0)=state(target_level_0);
  mini_state(1)=state(target_level_1);
  mini_state = m*mini_state;
  state(target_level_0)=mini_state(0);
  state(target_level_1)=mini_state(1);
//    std::cout << "state" << state << "\n";
return;
}// }}}
void apply_sigma_y_like(itpp::cvec& state, int target_level_0, int target_level_1, double theta){// {{{
  std::complex<double> I(0,1);
  double cos_theta  = cos(-theta);
  double sin_theta  = sin(-theta);

  itpp::cmat m(2,2);
  itpp::cvec mini_state(2);
  m(0,0) = cos_theta;
  m(1,1) = cos_theta;
  m(0,1) = sin_theta;
  m(1,0) = -sin_theta;
  mini_state(0)=state(target_level_0);
  mini_state(1)=state(target_level_1);
  mini_state = m*mini_state;
  state(target_level_0)=mini_state(0);
  state(target_level_1)=mini_state(1);
//    std::cout << "state" << state << "\n";
return;
}// }}}
void apply_sigma_z_like(itpp::cvec& state, int target_level, double theta){// {{{
  std::complex<double> I(0,1);
//    std::cout << "state antes " << state << "\n";
  state(target_level)= state(target_level)*exp(-I*theta);
//    std::cout << "state despues " << state << "\n";
return;
}// }}}
void apply_variationONSITE(itpp::cvec& state, int number_gell_man_like, double epsilon){// {{{
  int d=state.size();
  ivec pars=IntegerDigits(number_gell_man_like,d, 2);
//   cout << "number_gell_man_like=" << number_gell_man_like << ", pars=" << pars << endl; 
//   cout << "number_gell_man_like=" << number_gell_man_like << endl; 
  if (pars(0) > pars(1)){
    apply_sigma_x_like(state,pars(0),pars(1),epsilon);
  }
  else if(pars(0) == pars(1)){
    apply_sigma_z_like(state,pars(0),epsilon);
  }
  else if(pars(0) < pars(1)){
    apply_sigma_y_like(state,pars(0),pars(1),epsilon);
  }
}//}}}
itpp::cvec apply_variation(itpp::cvec& stateorig, int number_gell_man_like, double epsilon){// {{{
  cvec state=stateorig;
//    std::cout << "state antes " << state << "\n";
  apply_variationONSITE(state, number_gell_man_like, epsilon);
//    std::cout << "state despues" << state << "\n";
  return state;
}//}}}
cmat GellManLikeHamiltonian(vec parameters){ // {{{

  std::complex<double> I(0,1);
  int d=isqrt(parameters.size());
  cmat H(d,d); H=0.;
  ivec pars;
  for (int number_gell_man_like=0; number_gell_man_like< d*d; number_gell_man_like++){
    pars=IntegerDigits(number_gell_man_like,d, 2);
    if (pars(0) > pars(1)){
      H(pars(0), pars(1))+=parameters(number_gell_man_like);
      H(pars(1), pars(0))+=parameters(number_gell_man_like);
    }
    else if(pars(0) == pars(1)){
      H(pars(0), pars(0))+=parameters(number_gell_man_like);
    }
    else if(pars(0) < pars(1)){
      H(pars(0), pars(1))+=-I*parameters(number_gell_man_like);
      H(pars(1), pars(0))+=I*parameters(number_gell_man_like);
    }

  }
  return H;
} // }}}

double distance_operators(const itpp::cmat& A,const itpp::cmat& B){ // {{{
return norm(A-B,"fro");
} // }}}
double f(vec alpha, cvec& state, cmat& target_state, double swap_probability){ // {{{

  std::complex<double> I(0,1);
  cmat H = GellManLikeHamiltonian(-alpha);
  cvec psi_varied = exponentiate(I*H)*state;
  cmat rho_test=coarse_graining_two_positions(psi_varied, swap_probability);
  return distance_operators(target_state,rho_test);
} // }}}
double f_single_component(double alpha_i, int i, cvec& state, cmat& target_state, double swap_probability){ // {{{
  std::complex<double> I(0,1);
//   cmat H = GellManLikeHamiltonian(-alpha);
  cvec psi_varied = apply_variation(state, i, alpha_i);
  cmat rho_test=coarse_graining_two_positions(psi_varied, swap_probability);
  return distance_operators(target_state,rho_test);
} // }}}
vec gradient_CG_old(cvec& state, double epsilon, cmat& target_state, double swap_probability){ // {{{
  int d=state.size();
  vec grad(d*d); 
  cvec psi_plus, psi_minus;
  cmat rho_test;
  cmat rho_plus, rho_minus;


//     cout << epsilon  << endl;
//   rho_test=coarse_graining_two_positions(target_state, swap_probability);

  for (int i=0; i<d*d; i++){
    psi_plus = apply_variation(state, i, epsilon);
    rho_plus = coarse_graining_two_positions(psi_plus, swap_probability);
//     cout << state  << endl;
//     cout << psi_plus  << endl;
    psi_minus=apply_variation(state, i, -epsilon);
    rho_minus = coarse_graining_two_positions(psi_minus, swap_probability);
//   dis=distance_operators(target_state,rho_test);
    grad(i)=distance_operators(rho_plus, target_state)-distance_operators(rho_minus,target_state);
    grad(i)=grad(i)/(2*epsilon) ; 
//     if (i==2){
// 
//       cout << "hola, i=" << i << endl;
//       cout << "state    = " << state << endl;
//       cout << "psi_plus = " << psi_plus << endl;
//       cout << "rho_plus = " << rho_plus << endl;
//       cout << "d(rhp_plus,target_state) " << distance_operators(rho_plus, target_state)<< endl;
//       cout << "grad("<<i << ")=" << grad(i) << endl << endl; 
//     }
//     cout << grad(i) << endl;
//     abort();
  }
  return grad;
} // }}}
vec gradient_CG(cvec& state, double epsilon, cmat& target_state, double swap_probability){ // {{{
  int d=state.size();
  vec grad(d*d); 
  for (int i=0; i<d*d; i++){
    grad(i)=(f_single_component(epsilon, i, state, target_state, swap_probability)- 
     f_single_component(-epsilon, i, state, target_state, swap_probability))/(2*epsilon);
  }
  return grad;
} // }}}

double backtracking_line_search_old(cvec& psi,cmat& target_state, double swap_probability){ // {{{

  cout << "esta obsoleta, porque luego se usa otro gradiente y eso hace pedos." << endl; 
  abort();
  cout << "entrando a backtracking_line_search" << endl; 
  int dim=psi.size(); 
  int counter=0;

  double f_x = f(zeros(dim*dim), psi, target_state, swap_probability);
  double alpha_j = f_x; // En los algortmos de wiki se recomienda poner uno, pero es un overestimate. Ya tenenos una idea del error con la f_x

  double epsilon_gradient=min(0.1, f_x/10);

  vec grad = gradient_CG(psi, epsilon_gradient, target_state, swap_probability);

  vec p = -grad;

  double c=0.5, tau=0.5,m = -dot(p,p);

  double t=-c*m, f_xpap = f(alpha_j*p, psi, target_state, swap_probability);

//   cout << "f_x: " << f_x << endl;
//   cout << "f_xpap: " << f_xpap << endl;
//   cout << "t: " << t << endl;
  int counter_max=10;

  while (f_x -f_xpap < alpha_j*t && counter<counter_max){
    alpha_j = tau* alpha_j;
    f_xpap = f(alpha_j*p, psi, target_state, swap_probability);
    counter++;
//     cout << std::fixed <<"f_xpap=" << std::setprecision(5) << f_xpap 
//       << ", Delta f="<< std::setprecision(5) << f_x -f_xpap
//       << ", alpha_j*t=" << std::setprecision(5) <<alpha_j*t   << endl;
    //       cout << f_xpap << endl;
    //       cout << f_x -f_xpap << ", " << alpha_j*t   << endl;
  }
  if (counter >= counter_max){
    cerr << "Pedo en backtracking_line_search 90834752 " << endl; 
    cerr << "counter: " << counter << endl; 
    abort();
  }

  cout << "saliendo de backtracking_line_search. Error: " << f_x -f_xpap << endl; 
  return alpha_j;
} // }}}
cvec evolve_state_gradient(cvec& psi, vec& grad, double& alpha_j){ // {{{
  std::complex<double> I(0,1);
  cmat H = GellManLikeHamiltonian(-alpha_j*grad);
  return exponentiate(-I*H)*psi;
} // }}}
double backtracking_line_search_and_update(cvec& psi,cmat& target_state, double swap_probability){ // {{{

//   cout << "entrando a backtracking_line_search" << endl; 
  int dim=psi.size(); 
  int counter=0, counter_max=30;

  double f_x = f(zeros(dim*dim), psi, target_state, swap_probability);
  double alpha_j = 1.; // En los algortmos de wiki se recomienda poner uno, pero es un overestimate. Ya tenenos una idea del error con la f_x
//   double alpha_j = f_x; // En los algortmos de wiki se recomienda poner uno, pero es un overestimate. Ya tenenos una idea del error con la f_x

  double epsilon_gradient=min(0.00001, f_x/100);
  vec grad = gradient_CG(psi, epsilon_gradient, target_state, swap_probability);
  vec p = -grad;
  double c=0.5, tau=0.5,m = -dot(p,p);
  double t=-c*m, f_xpap = f(alpha_j*p, psi, target_state, swap_probability);

//   cout << "f_x: " << f_x << endl << "f_xpap: " << f_xpap << endl 
//     << "Delta f: " << f_x - f_xpap << endl << "t: " << t << endl 
//     << "alpha_j*t: " << alpha_j*t<< endl;


  while (f_x -f_xpap < alpha_j*t && counter<counter_max){
    alpha_j = tau* alpha_j;
    f_xpap = f(alpha_j*p, psi, target_state, swap_probability);
    counter++;
//     cout << std::fixed <<"f_xpap=" << std::setprecision(5) << f_xpap 
//       << ", Delta f="<< std::setprecision(5) << f_x -f_xpap
//       << ", alpha_j*t=" << std::setprecision(5) <<alpha_j*t   << endl;
    //       cout << f_xpap << endl;
//     cout << "En backtracking_line_search: " << f_x -f_xpap << ", " << alpha_j*t   << endl;
  }
  if (counter >= counter_max){
    cerr << "Pedo en backtracking_line_search 90834752 " << endl; 
    cerr << "counter: " << counter << endl; 
    abort();
  }
  psi = evolve_state_gradient(psi, grad, alpha_j);

//   cout << "saliendo de backtracking_line_search. Error: " << f_x -f_xpap << endl; 
//   return alpha_j;
  return f(zeros(dim*dim), psi, target_state, swap_probability);
} // }}}

double update_state_backtracking(cvec& psi, cmat& target_state, double swap_probability){ // {{{
  cout << "entrando a update_state_backtracking" << endl;
  std::complex<double> I(0,1);
  double alpha_j= backtracking_line_search_old(psi, target_state, swap_probability);
  double epsilon_gradient = min(0.1, abs(alpha_j)/100);
//   vec grad = gradient_CG(psi, 0.001, target_state, swap_probability);
  vec grad = gradient_CG(psi, epsilon_gradient, target_state, swap_probability);
  cout << "update_state_backtracking alpha_j=" << alpha_j<< endl;
  cout << "epsilon_gradient=" << epsilon_gradient<< endl;
//   cout << "grad=" << grad << endl;
  psi = evolve_state_gradient(psi, grad, alpha_j);
  double error;
  error=distance_operators(target_state,
        coarse_graining_two_positions(psi, swap_probability));
  cout << "Error en update_state_backtracking: " << error << endl;
  cout << "saliendo de update_state_backtracking" << endl;
  abort();
  return error;
  //   return;
} // }}}
void cool_down_gradient(cvec& psi, cmat& target_state, double swap_probability, int counter_max_cooldown, double max_tolerance){ // {{{

  int dim=psi.size(); 
  double error = f(zeros(dim*dim), psi, target_state, swap_probability);
  int counter=0;
//   cout << "Inicio cool_down_gradient, psi=" << psi << endl;
//   cout << "En cool_down_gradient, max_tolerance=" << max_tolerance<< endl;
  while (error > max_tolerance && counter < counter_max_cooldown){
    error = backtracking_line_search_and_update(psi,target_state, swap_probability);
//     cout << "En cool_down_gradient, error=" << error << endl;
    counter++ ; 
  }
  if (counter >= counter_max_cooldown){
    cerr << "# Pedo rutina  cool_down_gradient 243987545" << endl; 
    //       cerr << "counter: " << counter << endl; 
    cerr << "# " <<std::setprecision(5)<< psi << endl;
    //         cout << std::setprecision(5)<< std::fixed <<  psi << endl;
    abort();
    cerr << "# Error final =" <<error<< endl;
  } 
//   cout << "Fin cool_down_gradient, psi=" << psi << endl;
  return;
  } // }}}
void montecarlo_step(cvec& psi, cmat& target_state, double& swap_probability, // {{{
    double delta, double beta, double& error){ 

  int dim=psi.size(); 
  vec grad = randn(dim*dim);
  cvec psip = evolve_state_gradient(psi, grad, delta);
  double errorp=distance_operators(target_state,
      coarse_graining_two_positions(psip,swap_probability));
  if (errorp < error){
    psi = psip;
    error=errorp;
  } else {
    if (randu() < exp(-beta*(errorp - error))){
      // cout << " SI da el salto." << endl;
      psi = psip;
      error=errorp;
    } else {
      //           cout << " No da el salto." << endl;
    }
  }
  return; 
} // }}}


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
  if //{{ option == loquesea
    (option == "empty"){ // {{{
      return 0;
    } // }}}
  else if (option == "exploregradientrace"){ // {{{
    std::complex<double> I(0,1);
    int dim = 4;
    cvec psi= RandomState(dim);
    cmat target_state=(eye(2) + sigma(vec_3(0.,0.0,0.5)))/2;
    double swap_probability=0.3, epsilon=0.1;
    for (int i=0; i<100000; i++){
      gradient_CG(psi, epsilon, target_state, swap_probability);
    }
    return 0;
  } // }}}
  else if (option == "test_f_vs_f_single"){ // {{{
    std::complex<double> I(0,1);
    int dim = dimension.getValue();
    cvec psi= RandomState(dim);
    cmat target_state=partial_trace(RandomState(dim), 2);
    vec alpha(dim*dim);
    double alpha_i=0.1, f1, f2, error=0., swap_probability=0.3;
    for (int i=0; i < dim*dim; i++){
      alpha=0.;
      alpha(i) = alpha_i;
      f1 =  f(alpha, psi, target_state, swap_probability);
      f2 = f_single_component(alpha_i, i, psi, target_state, swap_probability);
      error += abs(f1 - f2);
    }
    cout << "error total: " << error << endl;
    return 0;
  } // }}}
  else if (option == "compare_new_old_gradient"){ // {{{
    std::complex<double> I(0,1);
//     int dim = 4;
    int dim = dimension.getValue();
    cvec psi= RandomState(dim);
    cmat target_state=partial_trace(RandomState(dim), 2);
    double swap_probability=0.3, epsilon=0.1;
     cout <<  norm(gradient_CG(psi, epsilon, target_state, swap_probability)-
      gradient_CG_old(psi, epsilon, target_state, swap_probability)) << endl ;
    return 0;
  } // }}}
  else if (option == "exploregradientdirection"){ // {{{
    std::complex<double> I(0,1);
    //     cmat H;
    int dim = dimension.getValue(), direcciones_raras = 0;
    vec grad, alpha;
    double swap_probability=0.3, epsilon=0.001;
    double dis, learning_rate=0.1, distance_initial;

    cmat target_state=partial_trace(RandomState(dim), 2);
    cvec psi_0= RandomState(dim);

    distance_initial= f(zeros(dim*dim), psi_0, target_state, swap_probability);

    grad = gradient_CG(psi_0, epsilon, target_state, swap_probability);
    alpha = learning_rate*grad/norm(grad);
    double distancia_direccion_optima= f(-alpha, psi_0, target_state, swap_probability);

    cout << "distancia inicial=" << distance_initial << endl;
    cout << "distancia optima=" << distancia_direccion_optima << endl;
    //     double dis;

    for (int i=0; i<sample_size.getValue(); i++){
      grad=itpp::randn(dim*dim);
      alpha=learning_rate*grad/norm(grad);
      dis= f(alpha, psi_0, target_state, swap_probability);
      //       cout <<  dis << endl;
      if (dis < distancia_direccion_optima){
        cout <<  dis << endl;
        direcciones_raras +=1;
      }
    }
    cout << "Direcciones raras encontradas=" << direcciones_raras << endl;

    return 0;
  } // }}}
  else if (option == "testgellmanfastoperations"){ // {{{
    std::complex<double> I(0,1);
      cvec psi, psi_varied_slow, psi_varied_fast ;
//       int dim = 7;
    int dim = dimension.getValue();
      psi= RandomState(dim);
      psi= 0.;
      psi(0)= 1.;

      cmat H;

      vec parameters(dim*dim);
      double epsilon=0.1, error=0.;


      for (int number_gell_man_like=0; number_gell_man_like<dim*dim; number_gell_man_like++){
//       int number_gell_man_like=1;{
//         cout << "number_gell_man_like: " << number_gell_man_like << endl;
        parameters=0.;
        psi_varied_fast = apply_variation(psi, number_gell_man_like, epsilon);
        parameters( number_gell_man_like) = epsilon;
        H = GellManLikeHamiltonian(parameters);
//         cout << H << endl;
        psi_varied_slow = exponentiate(-I*H)*psi;

//         cout << psi << endl;
//         cout << psi_varied_fast << endl;
//         cout << psi_varied_slow << endl;
//         cout << norm(psi-psi_varied_fast) << ", "
//            << norm(psi_varied_slow-psi_varied_fast) << endl;
        error+=norm(psi_varied_slow-psi_varied_fast) ;
//         abort();
      }

        cout << "Error: " <<error << endl;

      return 0;
    } // }}}
  else if (option == "test_adan"){ // {{{
    std::complex<double> I(0,1);
    cvec psi, psi_varied;
    cmat H;
    int dim = 4;
    double dis, learning_rate;
    cmat target_state, rho_test;
    target_state=(eye(2) + sigma(vec_3(0.,0.0,0.5)))/2;

    psi= RandomState(dim);
    psi(0)=(1.+I)/sqrt(2); psi(1)=1-I; psi(2)=1-I/sqrt(2); psi(3)=1+I; psi=psi/norm(psi);
    //       apply_sigma_x_like_2_levels(psi, 0, 2, 0.1);

    //       cvec psi_plus= apply_variation(psi, 0, 0.1);
    cout << psi << endl;
    //       cout << psi_plus << endl;
    //       abort();
    double swap_probability=0.3, epsilon=0.1;
    rho_test=coarse_graining_two_positions(psi, swap_probability);
    dis=distance_operators(target_state,rho_test);
    vec grad;

    grad = gradient_CG(psi, epsilon, target_state, swap_probability);
    cout << "Gradiente: " << grad << endl;
      return 0;
    } // }}}
  else if (option == "exploregradientvaryrate"){ // {{{
    std::complex<double> I(0,1);
    cvec psi, psi_varied;
    cmat H;
    int dim = 4;
    double dis, learning_rate;
    cmat target_state, rho_test;
    target_state=(eye(2) + sigma(vec_3(0.,0.0,0.5)))/2;

    psi= RandomState(dim);

    // cout << psi << endl;
    // cout << psi_plus << endl;
    // abort();
    double swap_probability=0.3, epsilon=0.1;
    rho_test=coarse_graining_two_positions(psi, swap_probability);
//     cout << "rho_test" << rho_test << endl;
//     cout << "target_state" << target_state << endl;
    dis=distance_operators(target_state,rho_test);
//     cout << target_state - rho_test << endl; 
//     cout << norm( target_state - rho_test , "fro" )<< endl; 
//     cout << norm( rho_test-target_state   , "fro" )<< endl; 
//     cout <<  "diferencia" << rho_test-target_state<< endl; 
//     rho_test << endl;
    vec grad;

    grad = gradient_CG(psi, epsilon, target_state, swap_probability);
//     cout << "Gradiente: " << grad << endl;
//     cout << "norma target_state:" << norm(target_state, "fro" ) << endl;
//     cout << "norma rho_test:" << norm(rho_test, "fro" ) << endl;
//     cout << "distancia entre distance_operators(target_state,rho_test)=" << dis << endl;
//     abort(); 
    grad=grad/norm(grad);

    double max_learning_rate=1.;
    int number_points=10;
    H = GellManLikeHamiltonian(grad);
    for (int i=0; i<number_points; i++){
      learning_rate=i*max_learning_rate/number_points;
      psi_varied = exponentiate(I*learning_rate*H)*psi;
//       cout << "psi_varied: " << psi_varied << endl;
      rho_test=coarse_graining_two_positions(psi_varied, swap_probability);
      dis=distance_operators(target_state,rho_test);
      cout <<  "learning_rate: " <<  setprecision(2) << learning_rate << ", " << setprecision(5) <<  dis << endl;
    }


//       cout <<  grad << endl;
//       cout <<  coarse_graining_two_positions(psi, swap_probability) << endl;


//          rho = coarse_graining_all_qubits(psi);
//          cout << Purity(rho) << endl;
      return 0;
    } // }}}
  else if (option == "explorelinesearch"){ // {{{
    std::complex<double> I(0,1);
    int dim = dimension.getValue();
    double learning_rate;
    cmat rho_test;
    cmat target_state=partial_trace(RandomState(dim), 2);
    target_state=(eye(2) + sigma(vec_3(0.,0.0,0.5)))/2;

    cvec psi= RandomState(dim);
//     psi=0.; psi(0)=1; psi(1)=1; psi=psi/norm(psi);
//    psi(0)=0.79495; psi(1)=-0.08747; psi(2)=-0.47415; psi(3)=-0.36821;psi=psi/norm(psi);
//    psi(0)=0.76077; psi(1)=-0.07081; psi(2)=-0.48252; psi(3)=-0.42824;psi=psi/norm(psi);
//     psi(0) =-0.34155+0.06873*I; psi(1)=0.04518+0.15559*I; psi(2)=0.72485-0.56136*I; psi(3)=0.08904-0.06246*I; psi/norm(psi);

//     cout << std::setprecision(5) << "psi = " << psi << endl;
//     cout << std::setprecision(5) << "target_state = " << target_state << endl;

    double swap_probability=0.3, epsilon_gradient=0.1;


    double error_inicial = distance_operators(target_state,coarse_graining_two_positions(psi, swap_probability));
    double error=1;
    int counter =0, counter_max=40; 
    cout << "error_inicial=" << error_inicial << endl << endl;
    while (error > 0.00001 && counter < counter_max){
//       error=  update_state_backtracking(psi, target_state, swap_probability);

      error = backtracking_line_search_and_update(psi,target_state, swap_probability);


      cout << "Error iterado=" << error << endl;
//       cout << std::setprecision(5)<< std::fixed <<  psi << endl;
      counter++ ; 
//       return 1;
    }
    if (counter >= counter_max){
      cerr << "Pedo en opcion explorelinesearch o983457209 " << endl; 
      cerr << "counter: " << counter << endl; 
      abort();
    }

    cout << std::setprecision(5)<< std::fixed <<  psi << endl;


    return 0;
  } // }}}
  else if (option == "get_convergent_states"){ // {{{
    std::complex<double> I(0,1);
    int dim = dimension.getValue();
//     cmat target_state=partial_trace(RandomState(dim), 2);
    cmat target_state=(eye(2) + sigma(vec_3(0.,0.0,polarization_z.getValue())))/2;

    cvec psi;
//     double max_tolerance=0.01;
    double max_tolerance=errorTCLAP.getValue();
// TCLAP::ValueArg<double> error("","error", "Maximum error tolerated",false, 0.001,"double",cmd);

    double p=swap_probability.getValue(), error;
    int counter, counter_max=40; 
    for (int i=0; i< sample_size.getValue(); i++){
      error=1; counter=0; psi=RandomState(dim);
//    error_inicial=distance_operators(target_state,coarse_graining_two_positions(psi,p));
      while (error > max_tolerance && counter < counter_max){
        error = backtracking_line_search_and_update(psi,target_state, p);
        counter++ ; 
      }
      if (counter >= counter_max){
        cerr << "# Pedo en opcion get_slow_convergent_states 249385723985 " << endl; 
        //       cerr << "counter: " << counter << endl; 
        cerr << "# " <<std::setprecision(5)<< psi << endl;
//         cout << std::setprecision(5)<< std::fixed <<  psi << endl;
        //       abort();
        cerr << "# Error final =" <<error<< endl;
      } else {
        cout << std::setprecision(8)<< psi << endl;
      }
    }
    return 0;
  } // }}}
  else if (option == "get_states_bruto"){ // {{{
    std::complex<double> I(0,1);
    int dim = dimension.getValue();
    //     cmat target_state=partial_trace(RandomState(dim), 2);
    cmat target_state=(eye(2) + sigma(vec_3(0.,0.0,polarization_z.getValue())))/2;

    cvec psi;
    //     double max_tolerance=0.01;
    double max_tolerance=errorTCLAP.getValue();
    // TCLAP::ValueArg<double> error("","error", "Maximum error tolerated",false, 0.001,"double",cmd);

    double p=swap_probability.getValue(), error;
    int counter=0, counter_max=sample_size.getValue();
    int counter_total=0;
    while (counter < counter_max){
      psi=RandomState(dim);
      counter_total++;
      error= distance_operators(target_state,
          coarse_graining_two_positions(psi,p));
      if (error < max_tolerance){
        cout << std::setprecision(8)<< psi << endl;
        counter++;
      }
    }
    cout << "# Numero de estados probados:  " << counter_total << endl;
    return 0;
  } // }}}
  else if (option == "get_relative_states"){ // {{{
    std::complex<double> I(0,1);
    int dim = dimension.getValue();
    //     cmat target_state=partial_trace(RandomState(dim), 2);
    cmat target_state=(eye(2) + sigma(vec_3(0.,0.0,polarization_z.getValue())))/2;

    cvec psi;
    //     double max_tolerance=0.01;
    double max_tolerance=errorTCLAP.getValue();
    // TCLAP::ValueArg<double> error("","error", "Maximum error tolerated",false, 0.001,"double",cmd);

    double p=swap_probability.getValue(), error;
    int counter=0, counter_max=sample_size.getValue();
    int counter_total=0;
    while (counter < counter_max){
      psi=RandomState(dim);
      counter_total++;
      error= distance_operators(target_state,
          coarse_graining_two_positions(psi,p));
      if (error < max_tolerance){
        cout << std::setprecision(8)<< psi << endl;
        counter++;
      }
    }
    cout << "# Numero de estados probados:  " << counter_total << endl;
    return 0;
  } // }}}
  else if (option == "test_get_states_bruto"){ // {{{
    std::complex<double> I(0,1);
    int dim = dimension.getValue();
    //     cmat target_state=partial_trace(RandomState(dim), 2);
    cmat target_state=(eye(2) + sigma(vec_3(0.,0.0,polarization_z.getValue())))/2;

    cmat rho_tmp,r2,  rho(dim,dim);
    rho=0.;
    cvec psi, p2;
    //     double max_tolerance=0.01;
    double max_tolerance=errorTCLAP.getValue();
    // TCLAP::ValueArg<double> error("","error", "Maximum error tolerated",false, 0.001,"double",cmd);

    double p=swap_probability.getValue(), error, e2;
    int counter=0, counter_max=sample_size.getValue();
    int counter_total=0;
    while (counter < counter_max){
      psi=RandomState(dim);
      counter_total++;
      error= distance_operators(target_state,
          coarse_graining_two_positions(psi,p));
      if (error < max_tolerance){
        cout << "Error inicial:" << error << endl;

        p2 = psi; 
        spinchain::apply_magnetic_kick(p2, vec_3(0.,0.,2*pi*randu()));
        e2 = distance_operators(target_state, coarse_graining_two_positions(p2,p));
        cout << "Error al aplicar una rotacion en z por un angulo raro:" << error << endl;


        rho_tmp = Proyector(psi);



        spinchain::apply_magnetic_kick(p2, vec_3(0.,0.,.43));

//         p2 = psi; p2(1)=-psi(1);p2(2)=-psi(2); 
//         p2 = psi; p2(1)=-I*psi(1);p2(2)=I*psi(2); 
        r2 = Proyector(p2);
        e2 = distance_operators(target_state, coarse_graining_two_positions(p2,p));
//         r2 = rho_tmp;
//         r2(1,2) = -rho_tmp(1,2);
//         r2(2,1) = -rho_tmp(2,1);
        rho+=rho_tmp;
        counter++;
        cout << std::setprecision(4)<< rho_tmp(1,2) << ", "
             <<  rho(1,2)/counter << ", "
             <<  error << ", " << e2 << endl;
//         cout << std::setprecision(8)<< psi << endl;
//         cout << rho_tmp << endl;
//         cout << r2 << endl;
        abort();
      }
    }
    cout << "# Numero de estados probados:  " << counter_total << endl;
    return 0;
  } // }}}
  else if (option == "metropolis_errors"){ // {{{
    std::complex<double> I(0,1);
//     int dim = dimension.getValue();
    int dim = 4;
    double p=swap_probability.getValue() ;
    cmat rho=(eye(2) + sigma(vec_3(0.,0.0,polarization_z.getValue())))/2;
    cvec psi = RandomState(dim), psip ; 
    double errorp, error=distance_operators(rho,coarse_graining_two_positions(psi,p));
    //     double max_tolerance=0.01;
    double delta=delta_stepTCLAP.getValue(), beta=betaTCLAP.getValue();
    // TCLAP::ValueArg<double> error("","error", "Maximum error tolerated",false, 0.001,"double",cmd);
// TCLAP::ValueArg<double> betaTCLAP("","beta", "Inverse temperature for montecarlo",false, 1.,"double",cmd);
// TCLAP::ValueArg<double> delta_stepTCLAP("","delta", "Step size",false, 0.01,"double",cmd);

    int counter=0, counter_max=sample_size.getValue();
    vec grad;
    while (counter < counter_max){
      grad = randn(dim*dim);
      psip = evolve_state_gradient(psi, grad, delta);
      errorp=distance_operators(rho,coarse_graining_two_positions(psip,p));
      cout << error  << endl;
//       cout << "error:" << error << ", errorp:"<< errorp << endl;
      if (errorp < error){
        psi = psip;
        error=errorp;
      } else {
        if (randu() < exp(-beta*(errorp - error))){
//           cout << " SI da el salto." << endl;
          psi = psip;
          error=errorp;

        } else {
//           cout << " No da el salto." << endl;
        }
      }
      counter++;
    }
    return 0;
  } // }}}
  else if (option == "metropolis"){ // {{{
    int dim = 4;
    double p=swap_probability.getValue() ;
    cmat rho=(eye(2) + sigma(vec_3(0.,0.0,polarization_z.getValue())))/2;
    cmat rho_montecarlo(4,4);
    cvec psi = RandomState(dim) ; 

    int counter_max_cooldown=100;
    double max_tolerance_cool=0.01;
    double delta=delta_stepTCLAP.getValue(), beta=betaTCLAP.getValue();
    double error;

    error=distance_operators(rho,coarse_graining_two_positions(psi,p));
//     cout << "Error inicial:" << error << endl;
//     cout << error << endl;
//     cout << "Psi antes de cooldown:" << psi << endl;
    // Cool down
    cool_down_gradient(psi, rho, p, counter_max_cooldown, max_tolerance_cool);
    error=distance_operators(rho,coarse_graining_two_positions(psi,p));
//     cout << "Error despues de cooldown:" << error << endl;
//     cout << "Psi despues de cooldown:" << psi << endl;
//     cout << error << endl;

    // Relajacion
    int relaxation_time=40;
//       cout << "Antes de la relajacion" << endl;
    for (int t=0; t< relaxation_time; t++){
//       cout << "en la relajacion t=" << t << endl;
      montecarlo_step(psi, rho, p, delta, beta, error);
//       cout << error << endl;
    }
    // Start counting
    int evolution_time=400000;
    rho_montecarlo=0.;
    for (int t=0; t< evolution_time; t++){
      montecarlo_step(psi, rho, p, delta, beta, error);
      rho_montecarlo+= Proyector(psi);
//       cout << error << endl;
    }
    rho_montecarlo= rho_montecarlo/evolution_time;
    cout << rho_montecarlo << endl;

    return 0;
  } // }}}
  else if (option == "coarse_graining_two_positions_vary_probability"){ // {{{
      int dim=4;
      cvec psi;
      cmat rho_i, rho_f;

      for (int i=0; i<sample_size.getValue(); i++){  
        psi= RandomState(dim);
        rho_i = Proyector(psi);
        rho_f = coarse_graining_two_positions(psi, swap_probability.getValue());
        cout << Concurrence(rho_i) << " "
             << Concurrence(rho_f) << " "
             << Purity(rho_f) << " "
             << endl;
      }
//       cout << partial_trace_qubits(rhob, 2) << endl;
//     cout << eig <<endl;
    return 0;
  } // }}}
  //}}
} // }}}


