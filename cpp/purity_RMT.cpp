// include files {{{
#ifndef  purity_RMT_VARIOUS
#define purity_RMT_VARIOUS
#include <RMT.cpp>
// #include <purity_RMT.h>
#include <QuantumStateITPP.cpp>
#include <cmath>
#include <cfp_math.cpp>
#include <itpp_ext_math.cpp>
// }}}
namespace PurityRMT{
  itpp::Mat<std::complex<double> > SpectatorHamiltonian(int dim_env, int dim_cou, int dim_wit, double lambda, double percet_out=0.1){ // {{{
    int dimc=dim_cou*dim_wit;
    itpp::Mat<std::complex<double> > tmp(dimc*dim_env,dimc*dim_env);
    tmp=kron(itpp::eye_c(dimc),
        itpp::diag(itpp::to_cvec(RMT::FlatSpectrumGUE(dim_env,percet_out))))
      + lambda*kron(itpp::eye_c(dim_wit),
          RMT::RandomGUE(dim_cou*dim_env,"sigma_offdiag=1"));
    return tmp;
  } // }}}
  itpp::Mat<std::complex<double> > TwoQubitSpectatorHamiltonianQQE(int dim_env, double lambda, double percet_out=0.1, double split_qubit=0.){// {{{
    // The structure of the Hilbert space is 
    // mcH=C_2 \otimes C_2 mcH_\env 
    itpp::Mat<std::complex<double> > tmp(4*dim_env,4*dim_env);
    tmp=kron(itpp::eye_c(4),
        itpp::diag(itpp::to_cvec(RMT::FlatSpectrumGUE(dim_env,percet_out))))
      + lambda*kron(itpp::eye_c(2),RMT::RandomGUE(2*dim_env,"sigma_offdiag=1"));
    if (abs(split_qubit)>0.){
      itpp::cvec internal(2);
      internal(0)=split_qubit/2;
      internal(1)=-split_qubit/2;
      tmp=tmp+kron(kron(itpp::eye_c(2),itpp::diag(internal)),itpp::eye_c(dim_env));
    }
    return tmp;
  } // }}}
  itpp::Mat<std::complex<double> > TwoQubitSpectatorHamiltonianEQQ(int dim_env, double lambda, double percet_out=0.1, double split_qubit=0.){// {{{
    // The structure of the Hilbert space is 
    // mcH=mcH_\env \otimes C_2 \otimes C_2
    // H= H_{env} \otimes 1_4 + \lambda V_{\env,1} \otimes 1_2
    itpp::Mat<std::complex<double> > tmp(4*dim_env,4*dim_env);
    tmp=kron(itpp::diag(itpp::to_cvec(RMT::FlatSpectrumGUE(dim_env,percet_out))),
        itpp::eye_c(4))
      + lambda*kron(RMT::RandomGUE(2*dim_env,"sigma_offdiag=1"),itpp::eye_c(2));
    if (abs(split_qubit)>0.){
      itpp::cvec internal(2);
      internal(0)=split_qubit/2;
      internal(1)=-split_qubit/2;
      tmp=tmp+kron(kron(itpp::eye_c(dim_env),itpp::diag(internal)),itpp::eye_c(2));
    }
    return tmp;
  } // }}}
  itpp::Mat<std::complex<double> > TwoQubitCommonEnvironment(int dim_env, itpp::vec lambda){// {{{
    // The structure of the Hilbert space is 
    // H=H_\env \otimes H_1 \otimes H_0
    //           ACA VOY
    // 		itpp::Mat<std::complex<double> > tmp(4*dim_env,4*dim_env);
    double percet_out=0.1;
    int q=cfpmath::log_base_2(dim_env)+2;
    int n_e=4*(dim_env-1), n_e0=n_e+1, n_e1=n_e+2;

    return itppextmath::extend_qubit_operator(itpp::diag(RMT::FlatSpectrumGUE(dim_env,percet_out)), n_e, q) // H_\env
      + lambda(0)*itppextmath::extend_qubit_operator(RMT::RandomGUE(2*dim_env,"sigma_offdiag=1"),n_e0, q)  // V_\env,0
      + lambda(1)*itppextmath::extend_qubit_operator(RMT::RandomGUE(2*dim_env,"sigma_offdiag=1"),n_e1, q); // V_\env,1

    //                 return tmp;
  } // }}}
  itpp::Mat<std::complex<double> > QubitEnvironmentHamiltonian(int dim_env, double lambda, double percet_out=0.1){ // {{{
    return kron(itpp::eye_c(2),itpp::diag(itpp::to_cvec (RMT::FlatSpectrumGUE(dim_env,percet_out))))
      + lambda*RMT::RandomGUE(2*dim_env,"sigma_offdiag=1");
  } // }}}
  void EspectatorEigenvaluesBell(){ // {{{
    int dim_env=4;
    double lambda=0.01, max_time=1.5;
    int NumberHamiltonians=2, NumberInitialConditions=1, time_steps=5;
#ifdef ASK
    clog << "Inserte dim_env                 "; 
    cin >> dim_env;                 cout<<endl;
    clog << "Inserte lambda                  "; 
    cin >> lambda;                  cout<<endl;
    clog << "Inserte max_time                "; 
    cin >> max_time;                cout<<endl;
    clog << "Inserte time_steps              "; 
    cin >> time_steps;              cout<<endl;
    clog << "Inserte NumberHamiltonians      "; 
    cin >> NumberHamiltonians;      cout<<endl;
    clog << "Inserte NumberInitialConditions "; 
    cin >> NumberInitialConditions; cout<<endl;
#endif
    QuantumStateITPP psi_0=QuantumStateITPP("Bell"), psi;
    itpp::Mat<std::complex<double> > H,U,U_dagger,rho;
    itpp::Vec<double> eigenvalues_rho, eigenvalues;
    double t;
    for (int i_hamiltonians=0;i_hamiltonians<NumberHamiltonians;
        i_hamiltonians++){
      H=TwoQubitSpectatorHamiltonianQQE(dim_env, lambda);
      itpp::eig_sym(H, eigenvalues, U);
      U_dagger=itpp::hermitian_transpose(U);
      for (int i_initial_conditions=0;
          i_initial_conditions<NumberInitialConditions;
          i_initial_conditions++){
        psi_0=StateTensor(QuantumStateITPP("Bell"),
            QuantumStateITPP(dim_env,"Random"));
        // Aca paso a la base que diagonaliza 
        // el hamiltoniano
        psi_0=U_dagger*psi_0;
        for (int i_time=0; i_time<time_steps; i_time++){
          t=cfpmath::linear_interval(i_time, 
              time_steps, 0., max_time);
          psi=evolve_with_phases(psi_0, 
              eigenvalues, t);
          psi=U*psi;
          rho=psi.PartialTrace(4);
          itpp::eig_sym(rho,eigenvalues_rho);
          std::cout << t<<" ";
          for (int i=0;i<4;i++){
            std::cout << eigenvalues_rho(i)<<" ";
          }
          std::cout << std::endl;
        }
      }
    }
  } // }}}
  void Espectator2QubitsAvRhoArguments(const int dim_env, const double max_time, const double lambda, const int NumberHamiltonians, const int NumberInitialConditions, const int time_steps, const std::string state_type, double const split_qubit=0.){ // {{{
    QuantumStateITPP psi_0, psi;
    itpp::Mat<std::complex<double> > H,U,U_dagger,rho;
    itpp::Array<itpp::cmat> rho_average(time_steps);
    itpp::Vec<double> eigenvalues_rho, eigenvalues;
    double t;
    for (int i=0;i<time_steps;i++){rho_average(i)=itpp::zeros_c(4,4);}
    for (int i_hamiltonians=0;i_hamiltonians<NumberHamiltonians; i_hamiltonians++){
      H=TwoQubitSpectatorHamiltonianQQE(dim_env, lambda,0.2, split_qubit);
      itpp::eig_sym(H, eigenvalues, U);
      U_dagger=itpp::hermitian_transpose(U);
      for (int i_initial_conditions=0; i_initial_conditions<NumberInitialConditions; i_initial_conditions++){
        psi_0=StateTensor(QuantumStateITPP(4,state_type),
            QuantumStateITPP(dim_env,"Random"));
        psi_0=U_dagger*psi_0;
        for (int i_time=0; i_time<time_steps; i_time++){
          t=cfpmath::linear_interval(i_time, 
              time_steps, 0., max_time);
          psi=U*evolve_with_phases(psi_0, 
              eigenvalues, t);
          rho=psi.PartialTrace(4);
          rho_average(i_time)+=rho;
        }
      }
    }
    for (int i_time=0; i_time<time_steps; i_time++){
      rho_average(i_time)/=(NumberHamiltonians*NumberInitialConditions);
      std::cout<<cfpmath::linear_interval(i_time, time_steps, 0., max_time)<<" ";
      itppextmath::PrintCompactHermitian(rho_average(i_time));
      std::cout<<std::endl;
    }
  } // }}}
  void TestTwoQubitSpectatorHamiltonian(const int dim_env, const double lambda, const int NumberHamiltonians, double const split_qubit=0.){ // {{{
    itpp::ivec state;
    std::ifstream seed_file2("/tmp/ss2.txt");seed_file2>>state;seed_file2.close();
    itpp::vec eigen_env;
    RNG_set_state (state);
    std::cout <<"El malo!!"<<std::endl;
    eigen_env=RMT::FlatSpectrumGUE(dim_env,0.2);
  } // }}}
  void EspectatorEigenvalues2QubitsArguments(const int dim_env, const double max_time, const double lambda, const int NumberHamiltonians, const int NumberInitialConditions, const int time_steps, const std::string state_type, double const split_qubit=0.){ // {{{
    QuantumStateITPP psi_0, psi;
    itpp::Mat<std::complex<double> > H,U,U_dagger,rho;
    itpp::Vec<double> eigenvalues_rho, eigenvalues;
    double t;
    for (int i_hamiltonians=0;i_hamiltonians<NumberHamiltonians; 
        i_hamiltonians++){
      H=TwoQubitSpectatorHamiltonianQQE(dim_env, lambda,0.2, split_qubit);
      itpp::eig_sym(H, eigenvalues, U);
      U_dagger=itpp::hermitian_transpose(U);
      for (int i_initial_conditions=0;
          i_initial_conditions<NumberInitialConditions;
          i_initial_conditions++){
        psi_0=StateTensor(QuantumStateITPP(4,state_type),
            QuantumStateITPP(dim_env,"Random"));
        psi_0=U_dagger*psi_0;
        for (int i_time=0; i_time<time_steps; i_time++){
          t=cfpmath::linear_interval(i_time, 
              time_steps, 0., max_time);
          psi=U*evolve_with_phases(psi_0, 
              eigenvalues, t);
          rho=psi.PartialTrace(4);
          itpp::eig_sym(rho,eigenvalues_rho);
          std::cout << t<<" ";
          for (int i=0;i<4;i++){
            std::cout << eigenvalues_rho(i)<<" ";
          }
          std::cout<<itppextmath::Purity(rho)<<" "<<itppextmath::Concurrence(rho);
          itppextmath::PrintCompactHermitian(rho);
          std::cout<<itppextmath::vonNeumann(eigenvalues_rho)<<std::endl;
        }
      }
    }
  } // }}}
  void SingleQubitArguments(const int dim_env, const double max_time, const double lambda, const int NumberHamiltonians, const int NumberInitialConditions, const int time_steps, const std::string state_type){ // {{{
    QuantumStateITPP psi_0, psi;
    itpp::Mat<std::complex<double> > H,U,U_dagger,rho;
    itpp::Vec<double> eigenvalues_rho, eigenvalues;
    double t;
    for (int i_hamiltonians=0;i_hamiltonians<NumberHamiltonians; 
        i_hamiltonians++){
      H=QubitEnvironmentHamiltonian(dim_env, lambda);
      itpp::eig_sym(H, eigenvalues, U);
      U_dagger=itpp::hermitian_transpose(U);
      for (int i_initial_conditions=0;
          i_initial_conditions<NumberInitialConditions;
          i_initial_conditions++){
        psi_0=StateTensor(QuantumStateITPP(2,state_type),
            QuantumStateITPP(dim_env,"Random"));
        psi_0=U_dagger*psi_0;
        for (int i_time=0; i_time<time_steps; i_time++){
          t=cfpmath::linear_interval(i_time, 
              time_steps, 0., max_time);
          psi=U*evolve_with_phases(psi_0, 
              eigenvalues, t);
          rho=psi.PartialTrace(2);
          itpp::eig_sym(rho,eigenvalues_rho);
          std::cout << t<<" ";
          for (int i=0;i<2;i++){
            std::cout << eigenvalues_rho(i)<<" ";
          }
          std::cout<<itppextmath::Purity(rho);
          itppextmath::PrintCompactHermitian(rho);
          std::cout<<std::endl;
        }
      }
    }
  } // }}}
  void EspectatorEigenvalues2Qubits(){ // {{{
    int dim_env=4;
    double lambda=0.01, max_time=1.5;
    int NumberHamiltonians=2, NumberInitialConditions=1, time_steps=5;
    std::string state_type="Bell";
    state_type="Separable";
    std::clog << "Entrando en la opcion para una condicion inicial arbitraria hh"
      << std::endl;
#ifdef ASK
    std::clog << "Inserte dim_env                 ";
    std::cin >> dim_env;                 std::cout<<std::endl;
    std::clog << "Inserte lambda                  ";
    std::cin >> lambda;                  std::cout<<std::endl;
    std::clog << "Inserte max_time                "; 
    std::cin >> max_time;                std::cout<<std::endl;
    std::clog << "Inserte time_steps              "; 
    std::cin >> time_steps;              std::cout<<std::endl;
    std::clog << "Inserte NumberHamiltonians      "; 
    std::cin >> NumberHamiltonians;      std::cout<<std::endl;
    std::clog << "Inserte NumberInitialConditions "; 
    std::cin >> NumberInitialConditions; std::cout<<std::endl;
    std::clog << "Inserte el tipo de estado       ";
    std::cin >> state_type;               std::cout<<std::endl;
#endif
    EspectatorEigenvalues2QubitsArguments(dim_env,
        max_time, lambda, NumberHamiltonians, NumberInitialConditions,
        time_steps, state_type);
  } // }}}
  double Time_Estimate_Purity(const double purity, const double lambda, const double Heisenberg_time=2*itpp::pi){ // {{{
    return (1-purity)/(3*lambda*lambda*Heisenberg_time);
  } // }}}
  void get_time_info(const int dim_env,  const double lambda){ // {{{
    itpp::Mat<std::complex<double> > H, U, U_dagger;
    QuantumStateITPP psi_0, psi;

    double required_purity=0.9;
    itpp::Vec<double> eigenvalues;

    H=TwoQubitSpectatorHamiltonianEQQ(dim_env, lambda);
    psi_0=QuantumStateITPP(4*dim_env,"BellRandom");

    double t_m1=0,t=Time_Estimate_Purity(required_purity, lambda),t_new;
    double P_m1=1,P,P_new,epsilon_purity=0.0000001;
    unsigned counter_loop=0,max_counter=10;
    itpp::eig_sym(H, eigenvalues, U); U_dagger=itpp::hermitian_transpose(U);
    psi_0=U_dagger*psi_0;
    P=itppextmath::Purity( ( U*evolve_with_phases(psi_0,eigenvalues, t)).PartialTrace(4));
    while (abs(P-required_purity)>epsilon_purity ){
      t_new=t_m1+(t-t_m1)*(P_m1-required_purity)/(P_m1-P);
      P_new=itppextmath::Purity( ( U*evolve_with_phases(psi_0,eigenvalues, t_new)).PartialTrace(4));
      std::cout << t_new<<" "<<P_new<< std::endl;
      t_m1=t;t=t_new;P_m1=P;P=P_new;
      counter_loop++;
      if (counter_loop>max_counter) {
        std::cerr<<"Paila, muchas iteraciones get_time_info \n"; 
        exit(1);
      }
    }
  } // }}}
}
#endif // purity_RMT_VARIOUS
