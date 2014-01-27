#ifndef  SPIN_CHAIN
#define SPIN_CHAIN
#include "itpp_ext_math.cpp"

//usar namespace conflictúa con los headers de cuda
//using namespace itpp;
namespace spinchain{ // {{{ .h
void apply_magnetic_kick(itpp::cvec&, itpp::vec, int);
void apply_magnetic_kick(itpp::cvec&, itpp::vec);
void apply_ising_z(itpp::cvec&, itpp::vec&);
void apply_ising_z(itpp::cvec&, double);
void apply_ising_z(itpp::cvec&, double, int, int);
void apply_ising_z_spectator(itpp::cvec&, double, double);
void apply_ising_z_spectator_double_chain(itpp::cvec&, double, double, double, int);
void apply_ising_z_vinayak(itpp::cvec&, double, double);
void apply_ising_z_common_environment_chain(itpp::cvec&, double, double);
void apply_common_environment_chain(itpp::cvec&, double, double, itpp::vec);
void apply_chain(itpp::cvec& state, double J, itpp::vec magnetic_field);
void apply_ising_star(itpp::cvec& , double , double );
void apply_ising_star_double(itpp::cvec& , double , double, double, int);
void apply_ising_star_double_closed(itpp::cvec& , double , double, double, int);
void apply_ising_star_most(itpp::cvec& , double , double );
void apply_kick_star(itpp::cvec& , itpp::vec , itpp::vec );
void apply_kick_spectator_double_chain(itpp::cvec& , itpp::vec , itpp::vec, itpp::vec, int);
void apply_kick_star_most(itpp::cvec& , itpp::vec , itpp::vec );

itpp::cvec project_base_state(int, int, int);
itpp::cvec apply_external_reflection(itpp::cvec&);
itpp::cvec project_state_vertical_momentum(int , itpp::cvec& , int );
itpp::cvec apply_vertical_external_reflection(itpp::cvec& , int );
} // }}}
namespace spinchain{ // {{{
  // Advanced building blocks {{{
  itpp::cvec apply_chain_spit_state(itpp::cvec state,itpp::vec magnetic_field, double J){// {{{
    itpp::cvec tmp=state;
    apply_chain(tmp, J, magnetic_field);
    return tmp;
  } //}}}
  void apply_star_most(itpp::cvec& state, double J, itpp::vec magnetic_field, double J_interaction, itpp::vec local_magnetic_field){// {{{
    // This topology is the usual closed chain, but with a qubit (the q-1 qubit) coupled
    // uniformely to all members
    //     *   *
    //  *         *
    //       *     
    //  *         *
    //     *   *
    apply_ising_star_most(state, J, J_interaction);
    apply_kick_star_most(state, magnetic_field, local_magnetic_field);
    return;
  } //}}}
  void apply_star(itpp::cvec& state, double J, itpp::vec magnetic_field, double J_interaction, itpp::vec local_magnetic_field){// {{{
    // This topology is the usual closed chain, but with a qubit (the 0th qubit) coupled
    // uniformely to all members
    //     *   *
    //  *         *
    //       *     
    //  *         *
    //     *   *
    apply_ising_star(state, J, J_interaction);
    apply_kick_star(state, magnetic_field, local_magnetic_field);
    return;
  } //}}}
  void apply_star_double(itpp::cvec& state, double J, itpp::vec magnetic_field, double Jcoupling1,double Jcoupling2, itpp::vec local_magnetic_field,int qubits_env1){// {{{
   
    apply_ising_star_double(state, J, Jcoupling1,Jcoupling2,qubits_env1);
    apply_kick_star(state, magnetic_field, local_magnetic_field);
    return;
  } //}}}
  void apply_star_double_closed(itpp::cvec& state, double J, itpp::vec magnetic_field, double Jcoupling1,double Jcoupling2, itpp::vec local_magnetic_field,int qubits_env1){// {{{
   
    apply_ising_star_double_closed(state, J, Jcoupling1,Jcoupling2,qubits_env1);
    apply_kick_star(state, magnetic_field, local_magnetic_field);
    return;
  } //}}}
  void apply_dephasing_chain(itpp::cvec& state, double Jenv, itpp::vec magnetic_field_env,  double  J_interaction_qubit_env, double Delta){// {{{
// ApplyDephasingChain[psi0_, Delta_, Jenv_, benv_, Jinteraction_] 
//  statenew = ApplyMagnetickKick[statenew, {0, 0, Delta/2}, 0];
//  statenew = ApplyMagnetickKickStarEnvironment[statenew, benv];
//  (*U2 interno del env, las ising y la interaccion con el medio*)
//  statenew = ApplyIsingStar[statenew, Jenv, Jinteraction];
    itpp::vec bz(3);
    bz=0.;
    bz(2)=Delta/2;
    apply_star(state, Jenv, magnetic_field_env, J_interaction_qubit_env, bz);
    return;
  } //}}}
  void apply_chain(itpp::cvec& state, double J, itpp::vec magnetic_field){// {{{
    apply_ising_z(state, J);
    apply_magnetic_kick(state, magnetic_field);
    return;
  } //}}}
  void apply_spectator(itpp::cvec& state, double Jenv, double Jcoupling, itpp::vec magnetic_field){// {{{
    apply_ising_z_spectator(state, Jenv, Jcoupling);
    apply_magnetic_kick(state, magnetic_field);
    return;
  } //}}}
//esto tambien:
  void apply_spectator_double_chain(itpp::cvec& state, double Jenv, double Jcoupling1, double Jcoupling2,int qubits_env1, itpp::vec magnetic_field_central, itpp::vec magnetic_field_env1, itpp::vec magnetic_field_env2){// {{{
    apply_ising_z_spectator_double_chain(state, Jenv, Jcoupling1,Jcoupling2,qubits_env1);
    apply_kick_spectator_double_chain(state, magnetic_field_central, magnetic_field_env1, magnetic_field_env2, qubits_env1);
    return;
  } //}}}

  void apply_vinayak(itpp::cvec& state, double Jenv, double Jcoupling, itpp::vec magnetic_field){// {{{
    apply_ising_z_vinayak(state, Jenv, Jcoupling);
    apply_magnetic_kick(state, magnetic_field);
    return;
  } //}}}
  void apply_common_environment_chain(itpp::cvec& state, double Jenv, double Jcoupling, itpp::vec magnetic_field){// {{{
    apply_magnetic_kick(state, magnetic_field);
    apply_ising_z_common_environment_chain(state, Jenv, Jcoupling);
    return;
  } //}}}
  void apply_ising_z_vinayak(itpp::cvec& state, double Jenv, double Jcoupling){// {{{
    //   std::cout << "la j es " << J <<" y la constante es " << expmij << std::endl;
    int qubits=cfpmath::log_base_2(state.size());
    if (qubits%2 != 0){
      std::cerr << "Not an even chain, qubits="<<qubits<<" must be an even number" << std::endl;
      abort();
    }
    itpp::vec J(qubits); 
    J=Jenv; 
    //     J(0)=1.;
    //     J(1)=0.;
    //     J(2)=1.;
    //     J(3)=0.;
    J(qubits/2-1)=Jcoupling; J(qubits-1)=0.;
    apply_ising_z(state, J);
    return;
  } // }}}
  void apply_ising_z_spectator(itpp::cvec& state, double Jenv, double Jcoupling){// {{{
    //   std::cout << "la j es " << J <<" y la constante es " << expmij << std::endl;
    int qubits=cfpmath::log_base_2(state.size());
    itpp::vec J(qubits); 
    J=Jenv; 
    J(0)=0.; J(1)=Jcoupling; J(qubits-1)=0.;
    apply_ising_z(state, J);
    return;
  } // }}}
//intento mio:
 void apply_ising_z_spectator_double_chain(itpp::cvec& state, double Jenv, double Jcoupling1,double Jcoupling2,int qubits_env1){// {{{
    int qubits=cfpmath::log_base_2(state.size());
    itpp::vec J(qubits);
    J=Jenv;
    J(0)=0.; J(1)=Jcoupling1; J(qubits_env1+1)=Jcoupling2; J(qubits-1)=0.;
    apply_ising_z(state, J);
    return;
  } // }}}

  void apply_ising_z_common_environment_chain(itpp::cvec& state, double Jenv, double Jcoupling){// {{{
    //   std::cout << "la j es " << J <<" y la constante es " << expmij << std::endl;
    int qubits=cfpmath::log_base_2(state.size());
    itpp::vec J(qubits); 
    J=Jenv; 
    J(0)=0.; J(1)=Jcoupling; J(qubits-1)=Jcoupling;
    //     std::cerr << J << std::endl;
    apply_ising_z(state, J);
    return;
  } // }}}
  // }}}
  // Intermediate building blocks {{{
  void apply_kick_star_most(itpp::cvec& state, itpp::vec magnetic_field, itpp::vec local_magnetic_field){// {{{
    int qubits=cfpmath::log_base_2(state.size());
    // Kick in the central qubit
//     std::cout << "En apply_kick_star_most, local_magnetic_field="<<local_magnetic_field << std::endl;
    apply_magnetic_kick(state,local_magnetic_field,qubits-1);
    // kick in the rest
    for (int Position=0; Position<qubits-1; Position++){
      apply_magnetic_kick(state,magnetic_field,Position);
    }
    return;
  } // }}}
  void apply_kick_star(itpp::cvec& state, itpp::vec magnetic_field, itpp::vec local_magnetic_field){// {{{
    int qubits=cfpmath::log_base_2(state.size());
    // Kick in the central qubit
    apply_magnetic_kick(state,local_magnetic_field,0);
    // kick in the rest
    for (int Position=1; Position<qubits; Position++){
      apply_magnetic_kick(state,magnetic_field,Position);
    }
    return;
  } // }}}
  void apply_kick_spectator_double_chain(itpp::cvec& state, itpp::vec magnetic_field_central,itpp::vec magnetic_field_env1,
itpp::vec magnetic_field_env2, int qubits_env1){// {{{
    int qubits=cfpmath::log_base_2(state.size());
    
    // Kick in the central system: a Bell pair
    apply_magnetic_kick(state,magnetic_field_central,0);
    apply_magnetic_kick(state,magnetic_field_central,1);
    // kick in the near environment
    for (int Position=2; Position<qubits_env1+2; Position++){
      apply_magnetic_kick(state,magnetic_field_env1,Position);
    }
    // kick in the far environment
    for (int Position=qubits_env1+2; Position<qubits; Position++){
      apply_magnetic_kick(state,magnetic_field_env2,Position);
    }

    return;
  } // }}}

  void apply_magnetic_kick(itpp::cvec& state, itpp::vec magnetic_field){// {{{
    int qubits=cfpmath::log_base_2(state.size());
    for (int Position=0; Position<qubits; Position++){
      apply_magnetic_kick(state,magnetic_field,Position);
    }
    return;
  } // }}}
  void apply_ising_star_most(itpp::cvec& state, double J, double J_interaction){// {{{
//     std::cout << "En apply_ising_star_most J_interaction="<<J_interaction << std::endl;
    int qubits=cfpmath::log_base_2(state.size());
    // Interaction of the central qubit with the rest
    for (int Position=0; Position<qubits-1; Position++){
      apply_ising_z(state, J_interaction, qubits-1, Position);
    }
    for (int Position=0; Position<qubits-1; Position++){
      apply_ising_z(state, J, Position, (Position+1)%(qubits-1));

    }
    return;
  } // }}}
  void apply_ising_star(itpp::cvec& state, double J, double J_interaction){// {{{
//     std::cout << "En apply_ising_star J_interaction="<<J_interaction << std::endl;
    int qubits=cfpmath::log_base_2(state.size());
    // Interaction of the central qubit with the rest
//     std::cout << "Step 1 " << state << std::endl;
    for (int Position=1; Position<qubits; Position++){
//       std::cout << "Apply J_interaction between " <<  itpp::vec_2(0, Position) << std::endl ;
      apply_ising_z(state, J_interaction, 0, Position);
    }
//     std::cout << "Step 2 " << state << std::endl;
    for (int Position=1; Position<qubits-1; Position++){
//       std::cout << "Apply J  between " <<  itpp::vec_2(Position, (Position+1)%qubits) << std::endl ;
      apply_ising_z(state, J, Position, (Position+1)%qubits);

    }
//     std::cout << "Step 3 " << state << std::endl;
    apply_ising_z(state, J, qubits-1, 1);
//     std::cout << "Step 4 " << state << std::endl;
    return;
  } // }}}



void apply_ising_star_double(itpp::cvec& state, double J, double Jcoupling1,double Jcoupling2, int qubits_env1){// {{{
   
    int qubits=cfpmath::log_base_2(state.size());
    
      apply_ising_z(state, Jcoupling1, 0, 1);
  
    for (int Position=1; Position<qubits-1; Position++){


      apply_ising_z(state, J, Position, (Position+1)%qubits);

    }
    apply_ising_z(state, Jcoupling2, qubits_env1, qubits_env1+1);
    apply_ising_z(state, Jcoupling2, qubits-1, 1);
    return;
  } // }}}

 void apply_ising_star_double_closed(itpp::cvec& state, double J, double Jcoupling1,double Jcoupling2,int qubits_env1){// {{{
//     std::cout << "En apply_ising_star J_interaction="<<J_interaction << std::endl;
    int qubits=cfpmath::log_base_2(state.size());
    // Interaction of the central qubit with fthe rest
//     std::cout << "Step 1 " << state << std::endl;
    for (int Position=1; Position<qubits_env1; Position++){
//       std::cout << "Apply J_interaction between " <<  itpp::vec_2(0, Position) << std::endl ;
      apply_ising_z(state, Jcoupling1, 0, Position);
    }
//     std::cout << "Step 2 " << state << std::endl;
    for (int Position=1; Position<qubits_env1-1; Position++){
//       std::cout << "Apply J  between " <<  itpp::vec_2(Position, (Position+1)%qubits) << std::endl ;
      apply_ising_z(state, J, Position, (Position+1)%qubits_env1);

    }
//     std::cout << "Step 3 " << state << std::endl;
    apply_ising_z(state, J, qubits_env1-1, 1);
//     std::cout << "Step 4 " << state << std::endl;
   
    for (int Position=1; Position<qubits-qubits_env1; Position++){
//       std::cout << "Apply J_interaction between " <<  itpp::vec_2(0, Position) << std::endl ;
      apply_ising_z(state, Jcoupling2, Position, Position+qubits_env1);
    }
    
    for (int Position=qubits_env1+1; Position<qubits-qubits_env1-1; Position++){
//       std::cout << "Apply J  between " <<  itpp::vec_2(Position, (Position+1)%qubits) << std::endl ;
      apply_ising_z(state, J, Position, (Position+1)%qubits);

    }
    
    apply_ising_z(state, J, qubits-1, qubits_env1+1);
//     std::cout << "Step 4 " << state << std::endl;
   


 return;
  } // }}}


  void apply_ising_z(itpp::cvec& state, double J){// {{{
    int qubits=cfpmath::log_base_2(state.size());
    itpp::vec Jv(qubits);
    Jv=J;
    apply_ising_z(state, Jv);
    return;
  } // }}}
  void apply_ising_z(itpp::cvec& state, itpp::vec& J){// {{{
    //   std::cout << "la j es " << J <<" y la constante es " << expmij << std::endl;
    int qubits=cfpmath::log_base_2(state.size());
    for (int i=0; i<qubits; i++){
      apply_ising_z(state, J(i), i, (i+1)%qubits);
    } 
    return;
  } // }}}
  // }}}
  // Symmetries in the homogeneous case {{{
  class CompactSymmetricBaseMember{ // {{{
    public:
      int generator;
      int k; //Sector of the symmetry
      int qubits; //Sector of the symmetry
      bool degenerate;
      bool sign; // 
      // We have a generator, |n>, to which we apply P_k.
      // This should not be a ket with zero length. 
      // Consider product
      // x=<n|R P_k |n> (R is the external reflection operator, R |ijkl> = |lkji>)
      // 
      // if x=0, the generator is degenerate, (and the boolean has 
      // a true value) and we use both 
      // P_k |n> \pm K R P_k |n> (K is the antiunitary symmetry)
      // 
      // if x \ne 0, then on of the P_k |n> \pm K R P_k |n>
      // is not null and can be the basis. 
      // The sign is + if sign=true and - if sign=false

  }; // }}}
  std::ostream &operator<<(std::ostream &os, const CompactSymmetricBaseMember &psi){ // {{{
    os << "(Generator " << psi.generator << "; symmetry sector " << psi.k << "; qubits " 
      << psi.qubits << "; degenerate " << psi.degenerate << ", sign " << psi.sign << ")" ;
    return os;
  } // }}}
  itpp::Array<CompactSymmetricBaseMember> build_rotationally_symmetric_base_states_compact(int qubits, int sector){ //{{{
    //     std::cout << "Si buenas 1 " << std::endl;
    itpp::Array<CompactSymmetricBaseMember> basis_states;
    CompactSymmetricBaseMember basis_state;
    int dim=cfpmath::pow_2(qubits);
    itpp::Vec<bool> Consider(dim);
    Consider=true;
    int degeneration_period;
    itpp::ivec all_n_rotated, reversed_all_n_rotated;
    //     std::cout << "Si buenas 2 dim=" << dim << std::endl;
    int J;
    for (int n=0; n<dim; n++){
//       std::cout <<"studying n=" << n<< std::endl;
      J=cfpmath::primitive_period_bit_rotation(n, qubits); 
      if (sector*J%qubits!=0){
        Consider(n)=false;
//         std::cout << "faaaaaaaaaaaalse" << std::endl;
      }
      //       std::cout << "Si buenas A n=" << n << std::endl;
      if(Consider(n)){
        all_n_rotated=itppextmath::all_bit_rotations(n, qubits);
        //         std::cout << "all_n_rotated=" <<  all_n_rotated << std::endl;
        reversed_all_n_rotated.set_size(all_n_rotated.size());
        //         std::cout << "Si buenas before for n=" << n << std::endl;
        for (int i=0; i<all_n_rotated.size(); i++){
          //           std::cout << "Estudiando el caso n=" << n << ", i=" <<i << std::endl;
          reversed_all_n_rotated(i)=cfpmath::reverse_bits(all_n_rotated(i), qubits);
          //           std::cout << "all_n_rotated("<<i<<")="<< all_n_rotated(i) <<
          //             ", reversed_all_n_rotated("<<i<<")="<< reversed_all_n_rotated(i) <<std::endl;
          Consider(all_n_rotated(i))=false;
          Consider(reversed_all_n_rotated(i))=false;
        }
//         std::cout << "Si buenas B n=" << n << std::endl;
        basis_state.k = sector;
        basis_state.generator = n;
        basis_state.qubits = qubits;
        if(itppextmath::Contains(all_n_rotated, cfpmath::reverse_bits(n,qubits))){// entonces <n|R P_k |n> != 0 {{{
          basis_state.degenerate=true;
          degeneration_period=itpp::min_index(abs(reversed_all_n_rotated - n)) ;
          if (((2*degeneration_period*sector%qubits)==0) && (((degeneration_period*sector)%qubits) != 0)){

            basis_state.sign=false;
          } else {
            basis_state.sign=true;
          }
//           abort();
          basis_states = itpp::concat(basis_states,basis_state);
//           std::cout <<"Just added " << basis_state.generator << std::endl; 
        } // }}}
        else {// entonces <n|R P_k |n> == 0 // {{{
          basis_state.degenerate=false;
          basis_state.sign=true;
          basis_states = itpp::concat(basis_states,basis_state);
          basis_state.sign=false;
          basis_states = itpp::concat(basis_states,basis_state);
          // so we add both states with signs
          // entonces <n|R P_k |n> = 0

        } //}}}
      }
    }
    return basis_states;
  } // }}}
  itpp::Array<CompactSymmetricBaseMember> build_rotationally_symmetric_base_states_compact(int qubits){ //{{{
    //     std::cout << "Si buenas 1 " << std::endl;
    itpp::Array<CompactSymmetricBaseMember> basis_states(0);
    for (int k=0; k<qubits; k++){
      basis_states= itpp::concat(basis_states,build_rotationally_symmetric_base_states_compact(qubits, k));
    }
    return basis_states;
  } // }}}
  itpp::cvec DecodeCompactRotationallySymetricBasisState(CompactSymmetricBaseMember encoded_state){ // {{{
    int q=encoded_state.qubits, n=encoded_state.generator, k=encoded_state.k;
//     std::cout << "Decoding the state given by q="<<q<<", n="<<n<<", k="<<k
//       << ", degenerate=" <<encoded_state.degenerate 
//       << ", sign=" <<encoded_state.sign 
//       << std::endl;
    //     if (CompactSymmetricBaseMember.degenerate){
    //     if (true ){
    itpp::cvec state_tmp=project_base_state(k, n, q);
//     std::cout << "norma de estado temporal " << norm(state_tmp) << std::endl
//       << "Estado: " << state_tmp << std::endl;
    if ( encoded_state.sign ){
      itpp::cvec state=state_tmp+conj(apply_external_reflection(state_tmp));
//       std::cout << "norma de estado (x) " << norm(state) << std::endl;
      return state/norm(state);
    } else {
      itpp::cvec state=state_tmp-conj(apply_external_reflection(state_tmp));
//       std::cout << "norma de estado (y) " << norm(state) << std::endl;
//       abort(); 
      return state/norm(state);
    }
  } // }}}
  itpp::cvec apply_external_reflection(itpp::cvec& state_in){ // {{{
    // The operator is defined as
    // R|i_0 i_1 i_2 ... i_{n-1}> = |i_{n-1} ... i_1 i_0>
    itpp::Vec<bool> reflected(state_in.size());
    itpp::cvec state=state_in;
    int n_reflected; 
    std::complex<double> tmp;
    reflected=false;
    int qubits=cfpmath::log_base_2(state.size());
    for (int n=0; n<state.size(); n++){
      n_reflected=cfpmath::reverse_bits(n, qubits);
      if (!reflected(n_reflected)){
        tmp=state(n);
        state(n)=state(n_reflected);
        state(n_reflected)=tmp;
        reflected(n)=true;
        reflected(n_reflected)=true;
      }
    }
    return state;
  } // }}}
  itpp::cvec apply_rotation(itpp::cvec& state_in, int power){ // {{{
    // The operator is defined as
    // R^k|i_{n-1}  ... i_1 i_0> = |i_{k-1} ... i_0 i_{n-1} ... i_k>
    // i. e. is a rotation to the right of the bits. 
    int d=state_in.size();
    itpp::cvec state(d);
    int qubits=cfpmath::log_base_2(d);
    for (int n=0; n<d; n++){
      state(cfpmath::rotate_bits(n, qubits, power))=state_in(n);
    }
    return state;
  } // }}}
  itpp::cvec apply_rotation(itpp::cvec& state_in){ // {{{
    // The operator is defined as
    // R|i_{n-1}  ... i_1 i_0> = |i_0 i_{n-1} ... i_1>
    // i. e. is a rotation to the right of the bits. 
    int d=state_in.size();
    itpp::cvec state(d);
    int qubits=cfpmath::log_base_2(d);
    for (int n=0; n<d; n++){
      state(cfpmath::rotate_bits(n, qubits))=state_in(n);
    }
    return state;
  } // }}}
  itpp::cvec project_base_state(int k, int base_state, int qubits){ // {{{
    // la idea es que agarro un n particular Veo si lo debo considerar. 
    // luego entonces marco los que no debo considerar porque son ciclos del man
    // luego reviso si proyecta a 0. 
    int J=cfpmath::primitive_period_bit_rotation(base_state, qubits); 
    itpp::cvec state(cfpmath::pow_2(qubits)); 
    state=0.;
    if (k*J%qubits!=0){ 
      return state; 
    }

//     std::complex<double> Imag=std::complex<double>(0,1);
    int n_rotated=base_state;
    for (int j=0; j<J; j++){
      state(n_rotated)=exp(-std::complex<double>(0,1)*2.*itpp::pi*double(j*k)/double(qubits));
      n_rotated=cfpmath::rotate_bits(n_rotated, qubits); 
    }
    return state/norm(state); 
    //Evalute if the state is cero
  } // }}}
  itpp::cvec project_state(int k, itpp::cvec& state_in){ // {{{
    // Creo que la formula general es 
    // P_k = \sum_{j=0}^L \varphi_{j,k} T^j
    int d=state_in.size();
    int q=cfpmath::log_base_2(d);
    itpp::cvec state(d);
    state=state_in;
    for (int i=1; i<q; i++){
      state+= exp(-2.*itpp::pi*std::complex<double>(0,1)*double(i*k/double(q)))* apply_rotation(state_in, i);
//       state+= exp(-2*(double(2)*Imag)*k*i/double(q)) * apply_rotation(state_in, i);
//       state+= exp(-2*itpp::pi*Imag*k*i/double(q)) * apply_rotation(state_in, i);
//       state+= exp(*k*i/double(q)) * apply_rotation(state_in, i);
    }
    return state/norm(state); 
    //Evalute if the state is cero
  } // }}}
  // }}}
  // Symmetries in the 2D case {{{
  class CompactSymmetric2DBaseMember{ // {{{
    public:
      // This are the equivalents to the generators. The provide
      // the basic elements of the total wave function
      itpp::Array<CompactSymmetricBaseMember*> generators; 
      // Symmetry Sector in the vertical direction
      int k_v;
      // Sign, to see if we consider P_k |n> \pm K R P_k |n>
      bool sign; 
      int vertical_dimension( ){ // {{{
        return generators.size();
      } // }}}
      int horizontal_dimension( ){ // {{{
        return (*(generators(0))).qubits;
      } // }}}
      int vertical_momentum( ){ // {{{
        return k_v;
      } // }}}
      int horizontal_momentum( ){ // {{{
        int total_k_h=0;
        for (int i=0; i<generators.size(); i++){
          total_k_h += (*(generators(i))).k;
        }
        return total_k_h%horizontal_dimension();
      } // }}}
  }; // }}}
  CompactSymmetric2DBaseMember CreateCompactSymmetric2DBaseMemberFromGenerator(CompactSymmetricBaseMember g, itpp::Array<CompactSymmetricBaseMember>& basis_states_horizontal){ // {{{
//     std::cout <<  "Getting in CreateCompactSymmetric2DBaseMemberFromGenerator\n" ;
    int nh = (basis_states_horizontal(0)).qubits;
//     int q=g.generators.size();
//     int nv = cfpmath::log_base_2(basis_states_horizontal.size());
    int q=g.qubits;
    int nv=q/nh;
//     std::cout << "CreateCompactSymmetric2DBaseMemberFromGenerator nh, nv= " << nh << ", " << nv << std::endl;
//     abort();
      
    CompactSymmetric2DBaseMember basis_2d_state; 
    basis_2d_state.generators.set_size(nv);
    int mask, basis_number;
    itpp::ivec  horizontal_basis_state_numbers(nv);
//     std::cout <<  "In CreateCompactSymmetric2DBaseMemberFromGenerator jklaljkdfa\n" ;
    for (int i_row=0; i_row < nv; i_row++){
      mask = ((cfpmath::pow_2(nh)-1)<<(nh*i_row)); 
      basis_number = (g.generator & mask) >> (nh*i_row);
      basis_2d_state.generators(i_row)=&basis_states_horizontal(basis_number);
    } // cout << "Smaller generators = " << horizontal_basis_state_numbers << endl;
    basis_2d_state.sign = g.sign;
    basis_2d_state.k_v  = g.k%(nv);

//     std::cout <<  "Getting out CreateCompactSymmetric2DBaseMemberFromGenerator\n" ;
    return basis_2d_state; 

  } // }}}
  void build_rotationally_symmetric_2D_base(int nh, int nv, itpp::Mat<itpp::Array<CompactSymmetric2DBaseMember> >& AllMatrixElements, itpp::Array<CompactSymmetricBaseMember>& basis_states_horizontal){ //{{{

    int q= nv*nh; int d=cfpmath::pow_2(q);
    itpp::Array<CompactSymmetricBaseMember> basis_states_2d;
    basis_states_horizontal=build_rotationally_symmetric_base_states_compact(nh);
    basis_states_2d=build_rotationally_symmetric_base_states_compact(q);
    CompactSymmetric2DBaseMember basis_2d_state; 
    AllMatrixElements.set_size(nh, nv);
    int kv, kh;
    for (int i=0; i< d; i++){
      basis_2d_state = CreateCompactSymmetric2DBaseMemberFromGenerator(basis_states_2d(i), basis_states_horizontal);
      kh=basis_2d_state.horizontal_momentum();
      kv=basis_2d_state.vertical_momentum();
      AllMatrixElements(kh, kv)=itpp::concat(AllMatrixElements(kh, kv), basis_2d_state);
//       itpp::cvec state2d = DecodeCompactSymmetric2DBaseMember(basis_2d_state);
    }
    return;
    
  } // }}}
  itpp::cvec DecodeCompactSymmetric2DBaseMember(CompactSymmetric2DBaseMember g){ // {{{
//     std::cout << "Gettin in DecodeCompactSymmetric2DBaseMember" << std::endl;

    int nh = g.horizontal_dimension(), nv = g.vertical_dimension();
//     std::cout << "DecodeCompactSymmetric2DBaseMember " << nh << ", " << nv << std::endl;
    itpp::Array<CompactSymmetricBaseMember> basis_states;
    basis_states=build_rotationally_symmetric_base_states_compact(nh);
 
    itpp::Array<itpp::cvec> statesbasic;
    statesbasic.set_size(nv);
    for (int  i_row=0; i_row < nv; i_row++){
      statesbasic(i_row)=DecodeCompactRotationallySymetricBasisState(*g.generators(i_row));
//       total_k_horizontal += basis_states(horizontal_basis_state_numbers(i_row)).k; 
    }
//     std::cout << "DecodeCompactSymmetric2DBaseMember 094832" << std::endl;
    itpp::cvec prestate, KRprestate;
    prestate=itppextmath::TensorProduct(statesbasic);
//     std::cout << "DecodeCompactSymmetric2DBaseMember 6554654 prestate.size()=" << prestate.size() << std::endl;
    prestate=project_state_vertical_momentum(g.k_v, prestate, nh);
//     std::cout << "DecodeCompactSymmetric2DBaseMember 87979987" << std::endl;
    KRprestate = apply_vertical_external_reflection(prestate, nh);
//     std::cout << "DecodeCompactSymmetric2DBaseMember 7856r41kj" << std::endl;
    KRprestate = itpp::conj(KRprestate);
//     KRprestate = itpp::conj(apply_vertical_external_reflection(prestate, nh));
//     std::cout << "Gettin out DecodeCompactSymmetric2DBaseMember" << std::endl;
    if (g.sign){
      return prestate + KRprestate;
    } else {
      return prestate - KRprestate;
    }
  } // }}}
  itpp::cvec apply_vertical_rotation(itpp::cvec& state_in, int horizontal_dimension){ // {{{
    // the bits are ordered as follows
    //
    // 8   9  10  11
    // 4   5   6   7
    // 0   1   2   3
    //
    // In this case, horizontal dimension is 4.
    // The above state gets transformed into 
    //
    // 4   5   6   7
    // 0   1   2   3
    // 8   9  10  11
    //
//     std::cout << "Entrando a apply_vertical_rotation" << std::endl;
    int d=state_in.size();
    itpp::cvec state(d);
    int qubits=cfpmath::log_base_2(d);
    for (int n=0; n<d; n++){
      state(cfpmath::rotate_bits(n, qubits, horizontal_dimension))=state_in(n);
    }
//     std::cout << "Saliendo de apply_vertical_rotation" << std::endl;
    return state;
  } // }}}
  itpp::cvec apply_horizontal_rotation(itpp::cvec& state_in, int horizontal_dimension){ // {{{
    // the bits are ordered as follows
    //
    // 8   9  10  11
    // 4   5   6   7
    // 0   1   2   3
    //
    // In this case, horizontal dimension is 4.
    // The above state gets transformed into 
    //
    // 11  8   9  10 
    //  7  4   5   6 
    //  3  0   1   2 
    //
    // 9  10  11   8   
    // 5   6   7   4   
    // 1   2   3   0   
    //
    // For example for a  3x2 |3> goes into |5> :
    // 0 0 0   => 0 0 0 
    // 1 1 0   => 1 0 1
    int d=state_in.size();
    itpp::cvec state(d);
    int qubits=cfpmath::log_base_2(d);
    for (int n=0; n<d; n++){
      state(cfpmath::apply_horizontal_rotation(n, qubits, horizontal_dimension))=state_in(n);
    }
    return state;
  } // }}}
  itpp::cvec apply_horizontal_rotation(itpp::cvec& state_in, int horizontal_dimension, int power){ // {{{
    // the documentation is found in the function without the power, that is
    //
    // itpp::cvec apply_horizontal_rotation(itpp::cvec&, int )
    //
    itpp::cvec state, tmp_state;
    state = state_in;
    for (int n=0; n<power; n++){
      tmp_state = apply_horizontal_rotation(state, horizontal_dimension);
      state=tmp_state;
    }
    return state;
  } // }}}
  itpp::cvec project_state_horizontal_momentum(int k, itpp::cvec& state_in, int horizontal_dimension){ // {{{
    // Creo que la formula general es 
    // P_k = \sum_{j=0}^L \varphi_{j,k} T^j
    itpp::cvec state;
    state=state_in;
    std::complex<double> phase;
    for (int i=1; i<horizontal_dimension; i++){
      phase = exp(-2.*itpp::pi*std::complex<double>(0,1)*double(i*k)/double(horizontal_dimension));
//       std::cout << i << ", " << phase << ", " << double(i*k)/double(q) << std::endl;
      state+= phase*apply_horizontal_rotation(state_in, horizontal_dimension, i);
    }
    return state/norm(state); 
  } // }}}
  itpp::cvec apply_vertical_rotation(itpp::cvec& state_in, int horizontal_dimension, int power){ // {{{
    // the bits are ordered as follows
    //
    // 8   9  10  11
    // 4   5   6   7
    // 0   1   2   3
    //
    itpp::cvec state, tmp_state;
    state = state_in;
    for (int n=0; n<power; n++){
      tmp_state = apply_vertical_rotation(state, horizontal_dimension);
      state=tmp_state;
    }
    return state;
  } // }}}
  itpp::cvec apply_vertical_external_reflection(itpp::cvec& state_in, int horizontal_dimension){ // {{{
//     std::cout << "In apply_vertical_external_reflection 0" << std::endl;
    itpp::Vec<bool> reflected(state_in.size());
    itpp::cvec state=state_in;
    int n_reflected; 
    std::complex<double> tmp;
    reflected=false;
    int qubits=cfpmath::log_base_2(state.size());
    for (int n=0; n<state.size(); n++){
//       n_reflected=cfpmath::reverse_bits(n, qubits);
      n_reflected = cfpmath::apply_vertical_external_reflection(n, qubits, horizontal_dimension);
//       std::cout << "In apply_vertical_external_reflection n=" << n 
//         <<", n_reflected=" << n_reflected << std::endl;
      if (!reflected(n_reflected)){
        tmp=state(n);
        state(n)=state(n_reflected);
        state(n_reflected)=tmp;
        reflected(n)=true;
        reflected(n_reflected)=true;
      }
    }
//     std::cout << "In apply_vertical_external_reflection 999999999999999999999" << std::endl;
    return state;
  } // }}}
  itpp::cvec project_state_vertical_momentum(int k, itpp::cvec& state_in, int horizontal_dimension){ // {{{
    // Creo que la formula general es 
    // P_k = \sum_{j=0}^L \varphi_{j,k} T^j
//     std::cout << "Entrando a project_state_vertical_momentum" << std::endl;
    itpp::cvec state;
    state=state_in;
//     std::cout << "En project_state_vertical_momentum 424234234 state_in.size()=" << state_in.size()<< std::endl;
    int qubits=cfpmath::log_base_2(state_in.size());
//     std::cout << "En project_state_vertical_momentum 64646234" << std::endl;
    int vertical_dimension = qubits/horizontal_dimension;

    std::complex<double> phase;
    for (int i=1; i<vertical_dimension; i++){
      phase = exp(-2.*itpp::pi*std::complex<double>(0,1)*double(i*k)/double(vertical_dimension));
//       std::cout << i << ", " << phase << ", " << double(i*k)/double(q) << std::endl;
      state+= phase*apply_vertical_rotation(state_in, horizontal_dimension, i);
    }
//     std::cout << "Saliendo de project_state_vertical_momentum" << std::endl;

    return state/norm(state); 
  } // }}}
  // }}}
  // Matrices for different evolution operators{{{
  itpp::cmat MatrixForIsing_star(int qubits, double J, itpp::vec magnetic_field, double J_interaction, itpp::vec local_magnetic_field, int sector){// {{{

    itpp::Array<CompactSymmetricBaseMember> basis_states;
    basis_states=build_rotationally_symmetric_base_states_compact(qubits-1, sector);
    itpp::cvec state_r, state_l;
    itpp::Array<itpp::cvec> state_q(2); state_q(0)=itpp::to_cvec(itpp::vec_2(1.,0.)); state_q(1)=itpp::to_cvec(itpp::vec_2(0.,1.));
    int d_sector=basis_states.size();
    itpp::cmat U(2*d_sector,2*d_sector);
    for (int i=0; i<d_sector; i++){ for (int iq=0; iq<2; iq++){
      state_r=itppextmath::TensorProduct(DecodeCompactRotationallySymetricBasisState(basis_states(i)),state_q(iq));
      apply_star(state_r, J, magnetic_field,J_interaction, local_magnetic_field);
      for (int j=0; j<d_sector; j++){for (int jq=0; jq<2; jq++){
        state_l=itppextmath::TensorProduct(DecodeCompactRotationallySymetricBasisState(basis_states(j)),state_q(jq));
        U(2*j+jq, 2*i+iq) =  itpp::dot(conj(state_l),state_r) ;
      }}
    }}
    return U;
  } // }}}
  itpp::cmat MatrixForIsing_star_most(int qubits, double J, itpp::vec magnetic_field, double J_interaction, itpp::vec local_magnetic_field, int sector){// {{{

    itpp::Array<CompactSymmetricBaseMember> basis_states;
    basis_states=build_rotationally_symmetric_base_states_compact(qubits-1, sector);
    itpp::cvec state_r, state_l;
    itpp::Array<itpp::cvec> state_q(2); state_q(0)=itpp::to_cvec(itpp::vec_2(1.,0.)); state_q(1)=itpp::to_cvec(itpp::vec_2(0.,1.));
    int d_sector=basis_states.size();
    itpp::cmat U(2*d_sector,2*d_sector);
//     std::cout << "Hola cabronsito, en la rutina problematica" << std::endl;
    for (int i=0; i<d_sector; i++){ for (int iq=0; iq<2; iq++){
      state_r=itppextmath::TensorProduct(state_q(iq), DecodeCompactRotationallySymetricBasisState(basis_states(i)));
//       std::cout << state_r << std::endl;
      apply_star_most(state_r, J, magnetic_field,J_interaction, local_magnetic_field);
//       std::cout << state_r << std::endl << std::endl;
      for (int j=0; j<d_sector; j++){for (int jq=0; jq<2; jq++){
        state_l=itppextmath::TensorProduct(state_q(jq), DecodeCompactRotationallySymetricBasisState(basis_states(j)));
        U(j+d_sector*jq, i+d_sector*iq) =  itpp::dot(conj(state_l),state_r) ;
//         U(2*j+jq, 2*i+iq) =  dot(conj(state_l),state_r) ;
      }}
    }}
//     std::cout << "Hola cabronsito, en la rutina problematica" << U<< std::endl;
//     abort();
    return U;
  } // }}}
  itpp::cmat MatrixForIsing_star_most(int qubits, double J, itpp::vec magnetic_field, double J_interaction, itpp::vec local_magnetic_field){// {{{
    int d=cfpmath::pow_2(qubits);
    itpp::cvec state(d);
    itpp::cmat U(d,d);
    U=0.;
    for (int i=0; i<d; i++){
      state=0.;
      state(i)=1.;
      apply_star_most(state, J, magnetic_field, J_interaction, local_magnetic_field);
      U.set_col(i,state);
    }
    return U;
  } // }}}
  itpp::cmat MatrixForIsing_star(int qubits, double J, itpp::vec magnetic_field, double J_interaction, itpp::vec local_magnetic_field){// {{{
    int d=cfpmath::pow_2(qubits);
    itpp::cvec state(d);
    itpp::cmat U(d,d);
    U=0.;
    for (int i=0; i<d; i++){
      state=0.;
      state(i)=1.;
      apply_star(state, J, magnetic_field, J_interaction, local_magnetic_field);
      U.set_col(i,state);
    }
    return U;
  } // }}}
  // }}}
  // Basic building blocks of evoution operators {{{
  void apply_ising_z(itpp::cvec& state, double J, int position1, int position2){// {{{
    std::complex<double> expmij=exp(-std::complex<double>(0,1)*J);
    std::complex<double> exppij=exp( std::complex<double>(0,1)*J);
    //   std::cout << "la j es " << J <<" y la constante es " << expmij << std::endl;
    for (int i=0; i<state.size(); i++){
      if(cfpmath::test_bit(i, position1) == cfpmath::test_bit(i, position2)){
        state(i)= expmij*state(i);
      } else {
        state(i)= exppij*state(i);
      }
    }
    return;
  } // }}}
  void apply_magnetic_kick(itpp::cvec& state, itpp::vec magnetic_field, int position){// {{{
    itpp::ivec pos(2);
    itpp::cvec moco(2);
    if (norm(magnetic_field)<0.000000000000001) return;
    for (int i=0; i<state.size()/2; i++){
      pos(0)=cfpmath::merge_two_numbers(0,i,cfpmath::pow_2(position));
      pos(1)=cfpmath::merge_two_numbers(1,i,cfpmath::pow_2(position));
      moco= itppextmath::exp_minus_i_b_sigma(magnetic_field)*state(pos);
      state(pos(0))= moco(0);
      state(pos(1))= moco(1);
    }
    return;
  } // }}}
  // }}}
  // Testing routines {{{
  itpp::cmat MatrixForIsingZ(int i, int j, int total){// {{{
    int mini=std::min(i,j), maxi=std::max(i,j);
    //     std::cout << "hola papi isingZ\n";

    itpp::cmat tmp(1,1);    tmp=1;
    tmp=itppextmath::TensorProduct(tmp,itpp::eye_c(cfpmath::pow_2(mini)));
    //     std::cout << "hola papi 1 "<< tmp.rows()<<" isingZ\n";
    tmp=itppextmath::TensorProduct(tmp,itppextmath::sigma(3));
    //     std::cout << "hola papi 2 "<< tmp.rows()<<" isingZ\n";
    tmp=itppextmath::TensorProduct(tmp,itpp::eye_c(cfpmath::pow_2(abs(i-j)-1)));
    //     std::cout << "hola papi 3 "<< tmp.rows()<<" isingZ\n";
    tmp=itppextmath::TensorProduct(tmp,itppextmath::sigma(3));
    //     std::cout << "hola papi 4 "<< tmp.rows()<<" isingZ\n";
    tmp=itppextmath::TensorProduct(tmp,itpp::eye_c(cfpmath::pow_2(total-maxi-1)));
    //     std::cout << "hola papi 5 "<< tmp.rows()<<" isingZ\n";

    // 
    //     std::cout << "hola papi isingZ xx \n";
    //     std::cout << "min="<<mini<<", max="<<maxi<<"\n";
    // //     return TensorProduct(itpp::eye_c(pow_2( min(i,j) )))
    return tmp; 
  } // }}}
  itpp::cmat MatrixForIsing_open_chain(double J, int total){// {{{
    itpp::cmat tmp(cfpmath::pow_2(total),cfpmath::pow_2(total));
    tmp=0.;
    for (int i=1; i<total; i++){
      tmp=tmp+MatrixForIsingZ(i-1, i, total);
    }
    return J*tmp;
  } // }}}
  itpp::cmat MatrixForIsing_chain(double J, int total){// {{{
    itpp::cmat tmp=MatrixForIsingZ(0, total-1, total);
    for (int i=1; i<total; i++){
      tmp=tmp+MatrixForIsingZ(i-1, i, total);
    }
    return J*tmp;
  } // }}}
  itpp::cmat Inefficient_magnetic(itpp::vec b, int total){// {{{
    itpp::cmat tmp(cfpmath::pow_2(total),cfpmath::pow_2(total));
    tmp=0.;
    for (int i=0; i<total; i++){
      tmp=tmp+itppextmath::sigma(b, i, total);
    }
    return tmp;
  } // }}}
  // }}}
  // Unsorted
  // itpp::vec eigenvalues(itpp::vec MagenticField, double Ising, int Dimension, std::string type_h){ // {{{
  //   itpp::cmat h= hamiltonian(MagenticField,Ising, Dimension, type_h);
  //   return eig_sym(h);
  // } // }}}
  // Unsorted 
} // }}}
#endif // SPIN_CHAIN
