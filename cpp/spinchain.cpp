#ifndef  SPIN_CHAIN
#define SPIN_CHAIN
#include <itpp_ext_math.cpp>

using namespace itpp;
namespace spinchain{ // {{{ .h
void apply_magnetic_kick(cvec&, vec, int);
void apply_magnetic_kick(cvec&, vec);
void apply_ising_z(cvec&, vec&);
void apply_ising_z(cvec&, double);
void apply_ising_z(cvec&, double, int, int);
void apply_ising_z_spectator(cvec&, double, double);
void apply_ising_z_vinayak(cvec&, double, double);
void apply_ising_z_common_environment_chain(cvec&, double, double);
void apply_common_environment_chain(cvec&, double, double, vec);
void apply_chain(cvec& state, double J, vec magnetic_field);
void apply_ising_star(cvec& , double , double );
void apply_ising_star_most(cvec& , double , double );
void apply_kick_star(cvec& , vec , vec );
void apply_kick_star_most(cvec& , vec , vec );
cvec project_base_state(int, int, int);
cvec apply_external_reflection(itpp::cvec&);
} // }}}
namespace spinchain{ // {{{
  // Advanced building blocks
  cvec apply_chain_spit_state(cvec state,vec magnetic_field, double J){// {{{
    cvec tmp=state;
    apply_chain(tmp, J, magnetic_field);
    return tmp;
  } //}}}
  void apply_star_most(cvec& state, double J, vec magnetic_field, double J_interaction, vec local_magnetic_field){// {{{
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
  void apply_star(cvec& state, double J, vec magnetic_field, double J_interaction, vec local_magnetic_field){// {{{
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
  void apply_chain(cvec& state, double J, vec magnetic_field){// {{{
    apply_ising_z(state, J);
    apply_magnetic_kick(state, magnetic_field);
    return;
  } //}}}
  void apply_spectator(cvec& state, double Jenv, double Jcoupling, vec magnetic_field){// {{{
    apply_ising_z_spectator(state, Jenv, Jcoupling);
    apply_magnetic_kick(state, magnetic_field);
    return;
  } //}}}
  void apply_vinayak(cvec& state, double Jenv, double Jcoupling, vec magnetic_field){// {{{
    apply_ising_z_vinayak(state, Jenv, Jcoupling);
    apply_magnetic_kick(state, magnetic_field);
    return;
  } //}}}
  void apply_common_environment_chain(cvec& state, double Jenv, double Jcoupling, vec magnetic_field){// {{{
    apply_magnetic_kick(state, magnetic_field);
    apply_ising_z_common_environment_chain(state, Jenv, Jcoupling);
    return;
  } //}}}
  void apply_ising_z_vinayak(cvec& state, double Jenv, double Jcoupling){// {{{
    //   std::cout << "la j es " << J <<" y la constante es " << expmij << std::endl;
    int qubits=cfpmath::log_base_2(state.size());
    if (qubits%2 != 0){
      std::cerr << "Not an even chain, qubits="<<qubits<<" must be an even number" << std::endl;
      abort();
    }
    vec J(qubits); 
    J=Jenv; 
    //     J(0)=1.;
    //     J(1)=0.;
    //     J(2)=1.;
    //     J(3)=0.;
    J(qubits/2-1)=Jcoupling; J(qubits-1)=0.;
    apply_ising_z(state, J);
    return;
  } // }}}
  void apply_ising_z_spectator(cvec& state, double Jenv, double Jcoupling){// {{{
    //   std::cout << "la j es " << J <<" y la constante es " << expmij << std::endl;
    int qubits=cfpmath::log_base_2(state.size());
    vec J(qubits); 
    J=Jenv; 
    J(0)=0.; J(1)=Jcoupling; J(qubits-1)=0.;
    apply_ising_z(state, J);
    return;
  } // }}}
  void apply_ising_z_common_environment_chain(cvec& state, double Jenv, double Jcoupling){// {{{
    //   std::cout << "la j es " << J <<" y la constante es " << expmij << std::endl;
    int qubits=cfpmath::log_base_2(state.size());
    vec J(qubits); 
    J=Jenv; 
    J(0)=0.; J(1)=Jcoupling; J(qubits-1)=Jcoupling;
    //     std::cerr << J << std::endl;
    apply_ising_z(state, J);
    return;
  } // }}}
  // Intermediate building blocks
  void apply_kick_star_most(cvec& state, vec magnetic_field, vec local_magnetic_field){// {{{
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
  void apply_kick_star(cvec& state, vec magnetic_field, vec local_magnetic_field){// {{{
    int qubits=cfpmath::log_base_2(state.size());
    // Kick in the central qubit
    apply_magnetic_kick(state,local_magnetic_field,0);
    // kick in the rest
    for (int Position=1; Position<qubits; Position++){
      apply_magnetic_kick(state,magnetic_field,Position);
    }
    return;
  } // }}}
  void apply_magnetic_kick(cvec& state, vec magnetic_field){// {{{
    int qubits=cfpmath::log_base_2(state.size());
    for (int Position=0; Position<qubits; Position++){
      apply_magnetic_kick(state,magnetic_field,Position);
    }
    return;
  } // }}}
  void apply_ising_star_most(cvec& state, double J, double J_interaction){// {{{
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
  void apply_ising_star(cvec& state, double J, double J_interaction){// {{{
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
  void apply_ising_z(cvec& state, double J){// {{{
    int qubits=cfpmath::log_base_2(state.size());
    vec Jv(qubits);
    Jv=J;
    apply_ising_z(state, Jv);
    return;
  } // }}}
  void apply_ising_z(cvec& state, vec& J){// {{{
    //   std::cout << "la j es " << J <<" y la constante es " << expmij << std::endl;
    int qubits=cfpmath::log_base_2(state.size());
    for (int i=0; i<qubits; i++){
      apply_ising_z(state, J(i), i, (i+1)%qubits);
    }
    return;
  } // }}}
  // Symmetries in the homogeneous case
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
  Array<CompactSymmetricBaseMember> build_rotationally_symmetric_base_states_compact(int qubits, int sector){ //{{{
    //     std::cout << "Si buenas 1 " << std::endl;
    Array<CompactSymmetricBaseMember> basis_states;
    CompactSymmetricBaseMember basis_state;
    int dim=cfpmath::pow_2(qubits);
    Vec<bool> Consider(dim);
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
          basis_states = concat(basis_states,basis_state);
//           std::cout <<"Just added " << basis_state.generator << std::endl; 
        } // }}}
        else {// entonces <n|R P_k |n> == 0 // {{{
          basis_state.degenerate=false;
          basis_state.sign=true;
          basis_states = concat(basis_states,basis_state);
          basis_state.sign=false;
          basis_states = concat(basis_states,basis_state);
          // so we add both states with signs
          // entonces <n|R P_k |n> = 0

        } //}}}
      }
    }
    return basis_states;
  } // }}}
  Array<CompactSymmetricBaseMember> build_rotationally_symmetric_base_states_compact(int qubits){ //{{{
    //     std::cout << "Si buenas 1 " << std::endl;
    Array<CompactSymmetricBaseMember> basis_states(0);
    for (int k=0; k<qubits; k++){
      basis_states= concat(basis_states,build_rotationally_symmetric_base_states_compact(qubits, k));
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
    cvec state_tmp=project_base_state(k, n, q);
      if ( encoded_state.sign ){
        cvec state=state_tmp+conj(apply_external_reflection(state_tmp));
        return state/norm(state);
      } else {
        cvec state=state_tmp-conj(apply_external_reflection(state_tmp));
        return state/norm(state);
      }
  } // }}}
  cvec apply_external_reflection(itpp::cvec& state_in){ // {{{
    // The operator is defined as
    // R|i_0 i_1 i_2 ... i_{n-1}> = |i_{n-1} ... i_1 i_0>
    Vec<bool> reflected(state_in.size());
    cvec state=state_in;
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
  itpp::cvec apply_rotation(itpp::cvec& state_in){ // {{{
    // The operator is defined as
    // R|i_{n-1}  ... i_1 i_0> = |i_0 i_{n-1} ... i_1>
    // i. e. is a rotation to the right of the bits. 
    int d=state_in.size();
    cvec state(d);
    int qubits=cfpmath::log_base_2(d);
    for (int n=0; n<d; n++){
      state(cfpmath::rotate_bits(n, qubits))=state_in(n);
    }
    return state;
  } // }}}
  cvec project_base_state(int k, int base_state, int qubits){ // {{{
    // la idea es que agarro un n particular Veo si lo debo considerar. 
    // luego entonces marco los que no debo considerar porque son ciclos del man
    // luego reviso si proyecta a 0. 
    int J=cfpmath::primitive_period_bit_rotation(base_state, qubits); 
    cvec state(cfpmath::pow_2(qubits)); 
    state=0.;
    if (k*J%qubits!=0){ 
      return state; 
    }

    std::complex<double> Imag=std::complex<double>(0,1);
    int n_rotated=base_state;
    for (int j=0; j<J; j++){
      state(n_rotated)=exp(-2*pi*Imag*j*k/qubits);
      n_rotated=cfpmath::rotate_bits(n_rotated, qubits); 
    }
    return state/norm(state); 
    //Evalute if the state is cero
  } // }}}
  cmat MatrixForIsing_star(int qubits, double J, vec magnetic_field, double J_interaction, vec local_magnetic_field, int sector){// {{{

    Array<CompactSymmetricBaseMember> basis_states;
    basis_states=build_rotationally_symmetric_base_states_compact(qubits-1, sector);
    cvec state_r, state_l;
    Array<cvec> state_q(2); state_q(0)=to_cvec(vec_2(1.,0.)); state_q(1)=to_cvec(vec_2(0.,1.));
    int d_sector=basis_states.size();
    cmat U(2*d_sector,2*d_sector);
    for (int i=0; i<d_sector; i++){ for (int iq=0; iq<2; iq++){
      state_r=itppextmath::TensorProduct(DecodeCompactRotationallySymetricBasisState(basis_states(i)),state_q(iq));
      apply_star(state_r, J, magnetic_field,J_interaction, local_magnetic_field);
      for (int j=0; j<d_sector; j++){for (int jq=0; jq<2; jq++){
        state_l=itppextmath::TensorProduct(DecodeCompactRotationallySymetricBasisState(basis_states(j)),state_q(jq));
        U(2*j+jq, 2*i+iq) =  dot(conj(state_l),state_r) ;
      }}
    }}
    return U;
  } // }}}
  cmat MatrixForIsing_star_most(int qubits, double J, vec magnetic_field, double J_interaction, vec local_magnetic_field, int sector){// {{{

    Array<CompactSymmetricBaseMember> basis_states;
    basis_states=build_rotationally_symmetric_base_states_compact(qubits-1, sector);
    cvec state_r, state_l;
    Array<cvec> state_q(2); state_q(0)=to_cvec(vec_2(1.,0.)); state_q(1)=to_cvec(vec_2(0.,1.));
    int d_sector=basis_states.size();
    cmat U(2*d_sector,2*d_sector);
//     std::cout << "Hola cabronsito, en la rutina problematica" << std::endl;
    for (int i=0; i<d_sector; i++){ for (int iq=0; iq<2; iq++){
      state_r=itppextmath::TensorProduct(state_q(iq), DecodeCompactRotationallySymetricBasisState(basis_states(i)));
//       std::cout << state_r << std::endl;
      apply_star_most(state_r, J, magnetic_field,J_interaction, local_magnetic_field);
//       std::cout << state_r << std::endl << std::endl;
      for (int j=0; j<d_sector; j++){for (int jq=0; jq<2; jq++){
        state_l=itppextmath::TensorProduct(state_q(jq), DecodeCompactRotationallySymetricBasisState(basis_states(j)));
        U(j+d_sector*jq, i+d_sector*iq) =  dot(conj(state_l),state_r) ;
//         U(2*j+jq, 2*i+iq) =  dot(conj(state_l),state_r) ;
      }}
    }}
//     std::cout << "Hola cabronsito, en la rutina problematica" << U<< std::endl;
//     abort();
    return U;
  } // }}}
  cmat MatrixForIsing_star_most(int qubits, double J, vec magnetic_field, double J_interaction, vec local_magnetic_field){// {{{
    int d=cfpmath::pow_2(qubits);
    cvec state(d);
    cmat U(d,d);
    U=0.;
    for (int i=0; i<d; i++){
      state=0.;
      state(i)=1.;
      apply_star_most(state, J, magnetic_field, J_interaction, local_magnetic_field);
      U.set_col(i,state);
    }
    return U;
  } // }}}
  cmat MatrixForIsing_star(int qubits, double J, vec magnetic_field, double J_interaction, vec local_magnetic_field){// {{{
    int d=cfpmath::pow_2(qubits);
    cvec state(d);
    cmat U(d,d);
    U=0.;
    for (int i=0; i<d; i++){
      state=0.;
      state(i)=1.;
      apply_star(state, J, magnetic_field, J_interaction, local_magnetic_field);
      U.set_col(i,state);
    }
    return U;
  } // }}}
  // Basic building blocks
  void apply_ising_z(cvec& state, double J, int position1, int position2){// {{{
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
  void apply_magnetic_kick(cvec& state, vec magnetic_field, int position){// {{{
    ivec pos(2);
    cvec moco(2);
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
  // Testing
  cmat MatrixForIsingZ(int i, int j, int total){// {{{
    int mini=std::min(i,j), maxi=std::max(i,j);
    //     std::cout << "hola papi isingZ\n";

    cmat tmp(1,1);    tmp=1;
    tmp=itppextmath::TensorProduct(tmp,eye_c(cfpmath::pow_2(mini)));
    //     std::cout << "hola papi 1 "<< tmp.rows()<<" isingZ\n";
    tmp=itppextmath::TensorProduct(tmp,itppextmath::sigma(3));
    //     std::cout << "hola papi 2 "<< tmp.rows()<<" isingZ\n";
    tmp=itppextmath::TensorProduct(tmp,eye_c(cfpmath::pow_2(abs(i-j)-1)));
    //     std::cout << "hola papi 3 "<< tmp.rows()<<" isingZ\n";
    tmp=itppextmath::TensorProduct(tmp,itppextmath::sigma(3));
    //     std::cout << "hola papi 4 "<< tmp.rows()<<" isingZ\n";
    tmp=itppextmath::TensorProduct(tmp,eye_c(cfpmath::pow_2(total-maxi-1)));
    //     std::cout << "hola papi 5 "<< tmp.rows()<<" isingZ\n";

    // 
    //     std::cout << "hola papi isingZ xx \n";
    //     std::cout << "min="<<mini<<", max="<<maxi<<"\n";
    // //     return TensorProduct(eye_c(pow_2( min(i,j) )))
    return tmp; 
  } // }}}
  cmat MatrixForIsing_open_chain(double J, int total){// {{{
    cmat tmp(cfpmath::pow_2(total),cfpmath::pow_2(total));
    tmp=0.;
    for (int i=1; i<total; i++){
      tmp=tmp+MatrixForIsingZ(i-1, i, total);
    }
    return J*tmp;
  } // }}}
  cmat MatrixForIsing_chain(double J, int total){// {{{
    cmat tmp=MatrixForIsingZ(0, total-1, total);
    for (int i=1; i<total; i++){
      tmp=tmp+MatrixForIsingZ(i-1, i, total);
    }
    return J*tmp;
  } // }}}
  cmat Inefficient_magnetic(vec b, int total){// {{{
    cmat tmp(cfpmath::pow_2(total),cfpmath::pow_2(total));
    tmp=0.;
    for (int i=0; i<total; i++){
      tmp=tmp+itppextmath::sigma(b, i, total);
    }
    return tmp;
  } // }}}
  // Various
  // vec eigenvalues(vec MagenticField, double Ising, int Dimension, std::string type_h){ // {{{
  //   cmat h= hamiltonian(MagenticField,Ising, Dimension, type_h);
  //   return eig_sym(h);
  // } // }}}
  // Unsorted 
} // }}}
#endif // SPIN_CHAIN

