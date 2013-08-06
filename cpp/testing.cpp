 // {{{ includes,  namespaces, and TCLAP (command line options)
#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
#include "cpp/cfp_math.cpp"
#include "cpp/itpp_ext_math.cpp"
#include "RMT.cpp"
#include <tclap/CmdLine.h>

#include "dev_random.cpp"
// using namespace itpp;
// using namespace itppextmath;
// using namespace cfpmath;
  TCLAP::CmdLine cmd("Command description message", ' ', "0.1");
  TCLAP::ValueArg<string> optionArg("o","option", "Option" ,false,"nichts", "string",cmd);
  TCLAP::ValueArg<unsigned int> seed("s","seed", "Random seed [0 for urandom]",false, 243243,"unsigned int",cmd);
  TCLAP::ValueArg<int> i1("","i1", "First integer argument",false, 1,"int",cmd);
  TCLAP::ValueArg<int> i2("","i2", "Second integer argument",false, 2,"int",cmd);
  TCLAP::ValueArg<int> i3("","i3", "Third integer argument",false, 1,"int",cmd);
  TCLAP::ValueArg<double> d1("","d1", "First double argument",false, 0.,"double",cmd);
  TCLAP::SwitchArg no_general_report("","no_general_report", "Print the general report", cmd);
// }}}
using namespace std;  
using namespace cfpmath;  
using namespace itppextmath;  

using namespace itpp;  
int main(int argc, char* argv[]){
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
  itpp::RNG_reset(semilla);
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
    (option == "test_insert_bit"){ // {{{
      cout << cfpmath::insert_bit(i1.getValue(),i2.getValue(),i3.getValue())<<std::endl;
    // }}}
  } else if (option == "test_extend_qubit_operator_coherent_with_sigma") { // {{{
    int qubits=i1.getValue(), j, coded_positions;
    double total_error=0.;
    ivec sigma_array(qubits);
    cmat sigma_long, sigma_short;
    for (int i=0; i<integer_pow(4,qubits) ;i++){
      j=i;
      for (int digit=0; digit< qubits;digit++){ sigma_array(digit)=j%4; j=j>>2; }
      sigma_long = sigma(sigma_array);
      coded_positions = bvec_to_int(sigma_array!=0);
      sigma(sigma_array(sigma_array!=0));
      sigma_short = extend_qubit_operator(sigma(sigma_array(sigma_array!=0)), coded_positions, qubits);
      total_error+= norm(sigma_long-sigma_short) ;
    }
    cout << total_error << endl;
  // }}}
  } else if (option == "test_rotate_bits") { // {{{
    cout << rotate_bits(i1.getValue(), i2.getValue()) << endl ;
  // }}}
  } else if (option == "test_partial_trace_mixed") { // {{{
    int qubits=i1.getValue(); int d=pow_2(qubits); 
    cmat rho(d,d);
    double x,y;
    for (int i=0; i<d; i++){
      for (int j=0; j<d; j++){
        cin >> x  >> y ;
        rho(i,j)=x+std::complex<double>(0,1)*y;
//         cout << rho(i,j) << endl;
      }
    }
//     cout << rho << "\n";
//     cvec state=RandomState(d); rho=Proyector(state);
    cmat rho_A=partial_trace_qubits(rho, i2.getValue());
//     cout << rho_A << "\n";
    int da=rho_A.cols();
    for (int i=0; i<da; i++){
      for (int j=0; j<da; j++){
        cout << real(rho_A(i, j)) << " " << imag(rho_A(i,j)) <<endl;
      }
    }
  // }}}
  } else if (option == "test_reverse_bits") { // {{{
    cout << reverse_bits(i1.getValue(), i2.getValue()) << endl ;
  // }}}
  } else if (option == "test_extract_digits") { // {{{
    int length=i1.getValue();
    int maxn=pow_2(length);
    int n1a, n1b, n2a, n2b;
    int error=0;
    for (int n=0; n<maxn; n++){
      for (int nwhich=0; nwhich<maxn; nwhich++){
        extract_digits(n, n1a, n2a, nwhich);
        extract_digits(n, length, n1b, n2b, nwhich);
        error += abs(n1a-n1b) + abs(n2a-n2b );
//         cout << n1a-n1b << ", " << n2a-n2b << n1a<< n1b << n2a << n2b << endl;
      }
    }
    cout << "The error is " << error << endl; 
  // }}}
  } else if (option == "test_partial_trace") { // {{{
    int q=i1.getValue();
    int d=pow_2(q);
    int mask=i2.getValue();;
    cmat rho;
    double x,y;
    cvec state(d);
    for (int i=0; i<d; i++){
      cin >> x  >> y ;
      state(i)=x+std::complex<double>(0,1)*y;
    }
//     cout << "mask " << mask << endl;
//     cout << "estado " << state << endl;
    rho= partial_trace_qubits(state,mask);
//     cout << rho << endl;
//     abort(); 
    for (int i=0; i<rho.cols(); i++){
      for (int j=0; j<rho.cols(); j++){
        cout << real(rho(i,j)) << " " << imag(rho(i,j)) <<endl;
      }
    }


  // }}}
  } else if (option == "test_all_bit_rotations") { // {{{
    ivec res=all_bit_rotations(i1.getValue(), i2.getValue());
    for (int i=0; i<res.size(); i++){ 
      cout <<res(i) << " ";
    }
    cout << endl ;
//     all_bit_rotations(i1.getValue(), i2.getValue()) ;
  // }}}
  } else if (option == "test_GetState_mixed") { // {{{
    int da=i1.getValue(), db=i2.getValue();
    int d=da*db;
    double x,y, t=d1.getValue();
    itpp::vec lambda;
    cvec state(d), state_out;
    cmat H(d,d), eigenvectors, rho_out;
    for (int i=0; i<d; i++){
      cin >> x  >> y ;
      state(i)=x+std::complex<double>(0,1)*y;
    }
    for (int i=0; i<d; i++){
      for (int j=0; j<d; j++){
        cin >> x  >> y ;
        H(i,j)=x+std::complex<double>(0,1)*y;
      }
    }
//     state_out= H*state;
    
    eig_sym(H,lambda, eigenvectors);
    rho_out = GetState(eigenvectors, lambda, state, t, da);
    for (int i=0; i<da; i++){
      for (int j=0; j<da; j++){
        cout << real(rho_out(i, j)) << " " << imag(rho_out(i,j)) <<endl;
      }
    }
//     for (int i=0; i<rho.cols(); i++){
//       for (int j=0; j<rho.cols(); j++){
//         cout << real(rho(i,j)) << " " << imag(rho(i,j)) <<endl;
//       }
//     }

//     all_bit_rotations(i1.getValue(), i2.getValue()) ;
  // }}}
  } else if (option == "test_GetState") { // {{{
    int da=i1.getValue(), db=i2.getValue();
    int d=da*db;
    double x,y, t=d1.getValue();
    itpp::vec lambda;
    cvec state(d), state_out;
    cmat H(d,d), eigenvectors;
    for (int i=0; i<d; i++){
      cin >> x  >> y ;
      state(i)=x+std::complex<double>(0,1)*y;
    }
    for (int i=0; i<d; i++){
      for (int j=0; j<d; j++){
        cin >> x  >> y ;
        H(i,j)=x+std::complex<double>(0,1)*y;
      }
    }
//     state_out= H*state;
    
    eig_sym(H,lambda, eigenvectors);
    state_out = GetState(eigenvectors, lambda, state, t);
    for (int i=0; i<state_out.size(); i++){
      cout << real(state_out(i)) << " " << imag(state_out(i)) <<endl;
    }
//     for (int i=0; i<rho.cols(); i++){
//       for (int j=0; j<rho.cols(); j++){
//         cout << real(rho(i,j)) << " " << imag(rho(i,j)) <<endl;
//       }
//     }

//     all_bit_rotations(i1.getValue(), i2.getValue()) ;
  // }}}
  } else if (option == "test_arbitrary_1_qubit_gate") { // {{{
    int q=i1.getValue();
    int d=pow_2(q);
    cvec state(d), state_1,state_2;
    double total_error = 0.;
    
    state=RandomState(d); 
//     state = 0.; state(1)=1.;

//     Array<cmat> SingleQubitGates(0);

//     SingleQubitGates=concat(sigma(0), SingleQubitGates);
//     SingleQubitGates=concat(sigma(1), SingleQubitGates);
//     SingleQubitGates=concat(sigma(2), SingleQubitGates);
//     SingleQubitGates=concat(sigma(3), SingleQubitGates);
//       
    // combinaciones de single qubit gates, leugo controles. 

//     for (int i=0; i<SingleQubitGates.size(); i++){
    for (int t=0; t<q; t++){
      // Hadamard gate {{{
      state_1 = state;
      state_2 = state;
      apply_gate(state_1, 1<<t, to_cmat(hadamard_matrix()));
      apply_hadamard(state_2, t);
      total_error += norm(state_1-state_2);
      if (total_error > 0.000000000001){
        cout << total_error << endl;
        cout << "Estado original=" << state << endl;
        cout << "state_1=" << state_1 << endl;
        cout << "state_2=" << state_2 << endl;
        cout << "h=" << to_cmat(hadamard_matrix()) << endl;
        abort();
      }
      // }}}
      // The four sigmas {{{
      for (int i=0; i<4; i++){
        state_1 = state;
        state_2 = state;
        apply_gate(state_1, 1<<t, sigma(i));
        apply_sigma(state_2, i, t);
        total_error += norm(state_1-state_2);
        if (total_error > 0.000000000001){
          cout << total_error << endl;
          cout << "(i, t)=(" << i << "," << t << ")" << endl;
          cout << "Estado original=" << state << endl;
          cout << "state_1=" << state_1 << endl;
          cout << "state_2=" << state_2 << endl;
          abort();
        }
      } // }}}
    }
    cout << "Error total = " << total_error << endl; 
      
  // }}}
  } else if (option == "test_arbitrary_2_qubit_gate") { // {{{
//     cout << "in option test_arbitrary_2_qubit_gate 1" << endl;
    int q=i1.getValue();
    int d=pow_2(q);
    cvec state(d), state_1,state_2;
    double total_error = 0.;
    Array<cmat> SingleQubitGates(5);
    SingleQubitGates(0)=sigma(0);
    SingleQubitGates(1)=sigma(1);
    SingleQubitGates(2)=sigma(2);
    SingleQubitGates(3)=sigma(3);
    SingleQubitGates(4)=to_cmat(hadamard_matrix());
    int ns=SingleQubitGates.size();
    cmat cu;
//     cout << "in option test_arbitrary_2_qubit_gate 2" << endl;
    state=RandomState(d); 
//     state = 0.; state(1)=1.;
    // las compuertas controladas, c-Pauli y c-R_k

    // All single qubit gates combined  {{{
    for (int t1=0; t1<q; t1++){ for (int t2=0; t2<t1; t2++){ // 
//       cout << "adsfa " << t1 << endl;
      for (int i1=0; i1<ns; i1++){ for (int i2=0; i2<ns; i2++){
//         cout << "adsfa (t1, t2, i1, i2)=(" << t1 << "," << t2 << "," << i1 << ","  << i2 << ")"<< endl;
        state_1 = state;
        state_2 = state;
        apply_gate(state_1, 1<<t1, SingleQubitGates(i1));
        apply_gate(state_1, 1<<t2, SingleQubitGates(i2));
        apply_gate(state_2, (1<<t1)+(1<<t2), TensorProduct(SingleQubitGates(i1),SingleQubitGates(i2)));
        total_error += norm(state_1-state_2);
        if (total_error > 0.000000000001){ // {{{
          cout << total_error << endl;
          cout << "sigma("<<i1<<") en "<< t1 << ", sigma("<<i2<<") en "<< t2 << endl;
          cout << "Estado original=" << state << endl;
          cout << "state_1=" << state_1 << endl;
          cout << "state_2=" << state_2 << endl;
          abort();
        } //}}}
      } } 
    } }// }}}
    // Swap gates {{{
    for (int t1=0; t1<q; t1++){ for (int t2=0; t2<t1; t2++){ 
      if (t1==t2){
        continue;
      }
//       cout << t1 << t2 << endl;
      state_1 = state;
      state_2 = state;
      apply_swap(state_1, t1, t2);
      apply_gate(state_2, (1<<t1)+(1<<t2),to_cmat(swap_matrix()));
      total_error += norm(state_1-state_2);
        if (total_error > 0.000000000001){ // {{{
          cout << total_error << endl;
//           cout << "sigma("<<i1<<") en "<< t1 << ", sigma("<<i2<<") en "<< t2 << endl;
          cout << "Estado original=" << state << endl;
          cout << "state_1=" << state_1 << endl;
          cout << "state_2=" << state_2 << endl;
          abort();
        } //}}}
    } }// }}}
    // control not gates {{{
//     cout << "Working on the cnot gate" << endl;
    for (int c=0; c<q; c++){ for (int t=0; t<q; t++){ 
      if (c==t){
        continue;
      }
//       cout << "c=" << c << ", t=" << t << endl;
      for (int i=0; i<ns; i++){ 
        state_1 = state; state_2 = state;
//         apply_control_u(state, c, t, SingleQubitGates(i));
//         cu = control_u_matrix( SingleQubitGates(i));
        apply_cnot(state_1, c, t);
        cu = to_cmat(cnot_matrix());
        if (c<t){
          cu=swap_matrix()*cu*swap_matrix();
        }
//         cout << "cu=" << cu << endl;
        apply_gate(state_2, (1<<c)+(1<<t),cu);
        total_error += norm(state_1-state_2);
        if (total_error > 0.000000000001){ // {{{
          cout << total_error << endl;
          //           cout << "sigma("<<i1<<") en "<< t1 << ", sigma("<<i2<<") en "<< t2 << endl;
          cout << "Control, target = " << c << ", " << t  << endl;
          cout << "Estado original=" << state << endl;
          cout << "state_1=" << state_1 << endl;
          cout << "state_2=" << state_2 << endl;
          abort();
        } //}}}
      }
    } }// }}}
    // control pauli gates {{{
//     cout << "Working on the cnot gate" << endl;
    for (int c=0; c<q; c++){ for (int t=0; t<q; t++){ 
      if (c==t){
        continue;
      }
//       cout << "c=" << c << ", t=" << t << endl;
      for (int i=0; i<ns; i++){ 
        state_1 = state; state_2 = state;
        apply_control_u(state_1, c, t, SingleQubitGates(i));
        cu = control_u_matrix( SingleQubitGates(i));
//         cout << cu << endl;
        if (c<t){
          cu=swap_matrix()*cu*swap_matrix();
        }
//         cout << "cu=" << cu << endl;
        apply_gate(state_2, (1<<c)+(1<<t),cu);
        total_error += norm(state_1-state_2);
        if (total_error > 0.000000000001){ // {{{
          cout << total_error << endl;
          //           cout << "sigma("<<i1<<") en "<< t1 << ", sigma("<<i2<<") en "<< t2 << endl;
          cout << "Control, target = " << c << ", " << t  << endl;
          cout << "Estado original=" << state << endl;
          cout << "state_1=" << state_1 << endl;
          cout << "state_2=" << state_2 << endl;
          abort();
        } //}}}
      }
    } }// }}}
    cout << "Error total = " << total_error << endl; 
      
  // }}}
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
}
