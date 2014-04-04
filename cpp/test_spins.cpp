// {{{ includes and using
#include <iostream>
#include <tclap/CmdLine.h>
// #include <RMT.cpp>
// #include <purity_RMT.cpp>
#include "itpp_ext_math.cpp"
#include "spinchain.cpp"
using namespace std;
using namespace itpp;
using namespace itppextmath;
using namespace cfpmath;
using namespace spinchain;
// }}}
// {{{ tlacp
  TCLAP::CmdLine cmd("Command description message", ' ', "0.1");

  TCLAP::ValueArg<string> optionArg("o","option", "Option" ,false,"nichts", "string",cmd);
  TCLAP::ValueArg<int> qubits("q","qubits", "Number of qubits",false, 4,"int",cmd);
  TCLAP::ValueArg<int> i1("","i1", "Integer parameter",false, 4,"int",cmd);
  TCLAP::ValueArg<int> i2("","i2", "Integer parameter",false, 4,"int",cmd);
  TCLAP::ValueArg<int> i3("","i3", "Integer parameter",false, 4,"int",cmd);
  TCLAP::ValueArg<int> position("","position", "The position of something",false, 1,"int",cmd);
  TCLAP::ValueArg<int> position2("","position2", "The position of something",false, 3,"int",cmd);
  TCLAP::ValueArg<double> Ising("","ising_z", "Ising interaction in the z-direction",false, 1.4,"double",cmd);
  TCLAP::ValueArg<double> coupling("","coupling", "Ising interaction for the coupling",false, 0.01,"double",cmd);
  TCLAP::ValueArg<double> Delta("","Delta", "Energy separation",false, 1,"double",cmd);
  TCLAP::ValueArg<double> bx("","bx", "Magnetic field in x direction",false, 1.4,"double",cmd);
  TCLAP::ValueArg<double> by("","by", "Magnetic field in y direction",false, 1.4,"double",cmd);
  TCLAP::ValueArg<double> bz("","bz", "Magnetic field in z direction",false, 1.4,"double",cmd);
// }}}
std::complex<double> Im(0,1);

int main(int argc, char* argv[]) { //{{{
// 	Random semilla_uran;
// 	itpp::RNG_reset(semilla_uran.strong());
//   	cout << PurityRMT::QubitEnvironmentHamiltonian(3,0.) <<endl;
// 	cout << RMT::FlatSpectrumGUE(5,.1) <<endl;
//
  cmd.parse( argc, argv );
  cout.precision(16);

  string option=optionArg.getValue();
  if (option=="test_kick_single_spin"){ // {{{// {{{

    int dim=pow_2(qubits.getValue());
    cvec state(dim);
    double x,y;
    ifstream myReadFile;
    myReadFile.open("/tmp/estado.dat");
    for (int i=0; i<dim; i++){
      myReadFile >> x >> y ;
      state(i) = complex<double>(x,y) ;
    }
    myReadFile.close();
    
    vec b(3); b(0)=bx.getValue(); b(1)=by.getValue(); b(2)=bz.getValue();
    apply_magnetic_kick(state,b,position.getValue());
    for (int i=0; i<dim; i++){
      cout << real(state(i)) << " " << real(-Im*state(i)) << endl;
    }//}}}
  } else if(option=="test_kick") { // {{{
    int dim=pow_2(qubits.getValue());
    cvec state(dim);
    double x,y;
    ifstream myReadFile;
    myReadFile.open("/tmp/estado.dat");
    for (int i=0; i<dim; i++){
      myReadFile >> x >> y ;
      state(i) = complex<double>(x,y) ;
    }
    myReadFile.close();
    
    vec b(3); b(0)=bx.getValue(); b(1)=by.getValue(); b(2)=bz.getValue();
    apply_magnetic_kick(state,b);
    for (int i=0; i<dim; i++){
      cout << real(state(i)) << " " << real(-Im*state(i)) << endl;
    }
    // }}}
  } else if(option=="test_hadamard") { // {{{
    int dim=pow_2(qubits.getValue());
    cvec state(dim);
    double x,y;
    for (int i=0; i<dim; i++){
      cin >> x >> y ;
      state(i) = complex<double>(x,y) ;
    }
    apply_hadamard(state,position.getValue());
    for (int i=0; i<dim; i++){
      cout << real(state(i)) << " " << real(-Im*state(i)) << endl;
    }
    // }}}
  } else if(option=="test_timing") {// {{{
    int dim=pow_2(qubits.getValue());
    cvec state(dim);
    state=RandomState(dim);
    double x,y;
        
    vec b(3); b(0)=bx.getValue(); b(1)=by.getValue(); b(2)=bz.getValue();

    for (int t=0; t<i1.getValue(); t++){
    apply_chain(state,Ising.getValue(),b);
    }
//     for (int i=0; i<dim; i++){
//       cout << real(state(i)) << " " << real(-Im*state(i)) << endl;
//     }//}}}
  } else if(option=="test_apply_chain") {// {{{
    int dim=pow_2(qubits.getValue());
    cvec state(dim);
    double x,y;
    ifstream myReadFile;
    myReadFile.open("/tmp/estado.dat");
    for (int i=0; i<dim; i++){
      myReadFile >> x >> y ;
      state(i) = complex<double>(x,y) ;
//       cout << x <<", "<<y<<endl;
    }
    myReadFile.close();
        
    vec b(3); b(0)=bx.getValue(); b(1)=by.getValue(); b(2)=bz.getValue();

    apply_chain(state,Ising.getValue(),b);
    for (int i=0; i<dim; i++){
      cout << real(state(i)) << " " << real(-Im*state(i)) << endl;
    }//}}}
  } else if(option=="test_kicked_ising_chain_eduardo") {// {{{
    int dim=pow_2(10);
    vec b(3);
    b(0)=1.;
    b(1)=0.;
    b(2)=1.;
    double J=1.;
    cvec state(dim);
    cvec state_out_eduardo(dim);
    double x,y;
    ifstream myReadFile;
    myReadFile.open("/tmp/c.dat");
//     myReadFile.open("/home/carlosp/Desktop/salida.dat");
    for (int i=0; i<dim; i++){
      myReadFile >> x >> y ;
      state(i) = complex<double>(x,y) ;
//       cout << i << ", " << state(i) << endl;
//       cout << x <<", "<<y<<endl;
    }
    for (int i=0; i<dim; i++){
      myReadFile >> x >> y ;
      state_out_eduardo(i) = complex<double>(x,y) ;
//       cout << x <<", "<<y<<endl;
    }
    for (int i=0; i<dim; i++){
//       state_out_eduardo(i) = complex<double>(x,y) ;
//       cout << i<< "; " << state(i) <<", "<<state_out_eduardo(i)<<endl;
    }
    myReadFile.close();
apply_chain(state,J,b);
// norm(state,state_out_eduardo)
// cout <<state<<endl<<endl<<state_out_eduardo;
cout <<norm(state) << endl;
cout <<norm(state_out_eduardo) << endl;
cout <<norm(state-state_out_eduardo) << endl;
    
//     apply_ising_z(state,J,b);
    for (int i=0; i<dim; i++){
//       cout << real(state(i)) << " " << real(-Im*state(i)) << endl;
    }//}}}
  } else if(option=="test_ising_single") {// {{{
    int dim=pow_2(qubits.getValue());
    cvec state(dim);
    double x,y;
    ifstream myReadFile;
    myReadFile.open("/tmp/estado.dat");
    for (int i=0; i<dim; i++){
      myReadFile >> x >> y ;
      state(i) = complex<double>(x,y) ;
//       cout << x <<", "<<y<<endl;
    }
    myReadFile.close();
    
    apply_ising_z(state,Ising.getValue(),position.getValue(), position2.getValue());
    for (int i=0; i<dim; i++){
      cout << real(state(i)) << " " << real(-Im*state(i)) << endl;
    }//}}}
  } else if(option=="test_dephasing_chain") {// {{{
    int dim=pow_2(qubits.getValue());
    cvec state(dim);
    double x,y;
    ifstream myReadFile;
    myReadFile.open("/tmp/estado.dat");
    for (int i=0; i<dim; i++){ myReadFile >> x >> y ; state(i) = complex<double>(x,y) ; }
    myReadFile.close();
    vec b_env(3); b_env(0)=bx.getValue(); b_env(1)=by.getValue(); b_env(2)=bz.getValue();
    double J_interaction_qubit_env=coupling.getValue();
    apply_dephasing_chain(state,Ising.getValue(),b_env, J_interaction_qubit_env, Delta.getValue());
    for (int i=0; i<dim; i++){
      cout << real(state(i)) << " " << real(-Im*state(i)) << endl;
    }//}}}
  } else if(option=="test_apply_common_environment_chain") {// {{{
    int dim=pow_2(qubits.getValue());
    cvec state(dim);
    double x,y;
    ifstream myReadFile;
    myReadFile.open("/tmp/estado.dat");
    for (int i=0; i<dim; i++){ myReadFile >> x >> y ; state(i) = complex<double>(x,y) ; }
    myReadFile.close();
    vec b(3); b(0)=bx.getValue(); b(1)=by.getValue(); b(2)=bz.getValue();
    apply_common_environment_chain(state,Ising.getValue(),coupling.getValue(), b);
    for (int i=0; i<dim; i++){
      cout << real(state(i)) << " " << real(-Im*state(i)) << endl;
    }//}}}
  } else if(option=="test_reduced_MatrixForIsing_star") {// {{{
    int q=qubits.getValue();
    int dim=pow_2(q);
    vec b(3); b(0)=bx.getValue(); b(1)=by.getValue(); b(2)=bz.getValue();
    vec local_b(3); local_b=0.;
    double  J=itpp::randu(), J_interaction=0.;
//     cmat U=MatrixForIsing_star(q, J, b, J_interaction, local_b);
    cvec state_env, state, state_q;
    cmat rho;
    state_env=RandomState(dim);
    state_q=RandomState(2);
    state_q=0.;
    state_q(0)=1.;
    state=TensorProduct(state_env,state_q);
    apply_star(state,J, b, J_interaction, local_b);
    rho= partial_trace_qubits(state,1);
    cout << rho << " "  << endl;
  //}}}
  } else if(option=="test_MatrixForIsing_star") {// {{{
    int q=qubits.getValue();
    int dim=pow_2(q);
    vec b=itpp::randu(3), local_b=itpp::randu(3);double  J=itpp::randu(), J_interaction=itpp::randu();
    cmat U=MatrixForIsing_star(q, J, b, J_interaction, local_b);
    cvec state, state_matrix;
    state=RandomState(dim); state_matrix=state;
    apply_star(state,J, b, J_interaction, local_b);
    cout << norm(U*state_matrix - state) << " "  << endl;
  //}}}
  } else if(option=="test_extended_base_orthonormality") {// {{{
    cout << "quiero probar la unitariedad en la base computacional"
      << " del operador de estrella y verificar la completez de la "
      << "base rotacional extendida a un qubit" << endl;
    int q=qubits.getValue(); // con i1=2 hay problemas
    int d=pow_2(q);
    vec b=itpp::randu(3), local_b=itpp::randu(3);
//     double  J=itpp::randu(), J_interaction=itpp::randu(); J=0.; J_interaction=0.;
    b=1.; local_b=0.; 
    b=vec_3(pi/2.,0.,0.);
    Array<cvec> state_q(2);
    cvec state_r, state_l;
    state_q(0)=to_cvec(vec_2(1.,0.)); state_q(1)=to_cvec(vec_2(0.,1.));
//     cmat U=MatrixForIsing_star(q, J, b, J_interaction, local_b);
    cmat U(d,d);
    Array<CompactSymmetricBaseMember> basis_states=build_rotationally_symmetric_base_states_compact(q);
    for (int i=0; i<d/2; i++){ for (int iq=0; iq<2; iq++){
      state_r=TensorProduct(DecodeCompactRotationallySymetricBasisState(basis_states(i)),state_q(iq));
      cout << "norm(state_r)=" << norm(state_r) << endl;
      for (int j=0; j<d/2; j++){for (int jq=0; jq<2; jq++){
        state_l=TensorProduct(DecodeCompactRotationallySymetricBasisState(basis_states(j)),state_q(jq));
        U(2*j+jq, 2*i+iq) =  dot(conj(state_l),state_r) ;
      }}
    }}
    
    cout << "Primero debemos escribir eficientemente la U en forma de bloque. Luego "
      <<" calcular el canal cuantico y verificar que da lo mismo" << endl;
    cout << "Is U=1 " << test_unit_matrix(U)<< endl;
//     cout << "Is U unitary? " << test_unitarity(U)<< endl;

    // }}}
  } else if(option=="test_start_unitarity") {// {{{
    cout << "quiero probar la unitariedad en la base computacional"
      << " del operador de estrella y verificar la completez de la "
      << "base rotacional extendida a un qubit" << endl;
    int q=qubits.getValue(); // con i1=2 hay problemas
    int d=pow_2(q);
    vec b=itpp::randu(3), local_b=itpp::randu(3);double  J=itpp::randu(), J_interaction=itpp::randu();
    Array<cvec> state_q(2); state_q(0)=to_cvec(vec_2(1.,0.)); state_q(1)=to_cvec(vec_2(0.,1.));
    cvec state_r, state_l;
    cmat U(d,d);
    Array<CompactSymmetricBaseMember> basis_states=build_rotationally_symmetric_base_states_compact(q-1);
    for (int i=0; i<d/2; i++){ for (int iq=0; iq<2; iq++){
      state_r=TensorProduct(DecodeCompactRotationallySymetricBasisState(basis_states(i)),state_q(iq));
      apply_star(state_r, J, b,J_interaction, local_b);
      // cout << "norm(state_r)=" << norm(state_r) << endl;
      // cout << "state_r=" <<state_r << endl;
      for (int j=0; j<d/2; j++){for (int jq=0; jq<2; jq++){
        state_l=TensorProduct(DecodeCompactRotationallySymetricBasisState(basis_states(j)),state_q(jq));
        // cout << "state_l=" <<state_l << endl;
        U(2*j+jq, 2*i+iq) =  dot(conj(state_l),state_r) ;
      }}
    }}
    
    cout << "Is U unitary? " << test_unitarity(U)<< endl;

    // }}}
  } else if(option=="test_start_equivalence") {// {{{

    vec b=itpp::randu(3), local_b=itpp::randu(3);double  J=itpp::randu(), J_interaction=itpp::randu();
    cvec state_least, state_most, state_most_least;
    state_least=RandomState(pow_2(qubits.getValue()));
    state_most=apply_rotation(state_least);
    apply_star_most(state_most, J, b, J_interaction, local_b);
    apply_star(state_least, J, b, J_interaction, local_b);
    state_most_least = apply_rotation(state_least);
    cout << "Are the two operators equivalen? " << norm(state_most_least-state_most) << endl;

    // }}}
  } else if(option=="test_project_base_state") {// {{{
    Array<CompactSymmetricBaseMember> basis_states, tmp_basis;
    cvec state_l, state_r, state;
    int q=qubits.getValue();
    int d=pow_2(q);cmat U(d,d);
    double error=0.;
    // orthonormality in each sector {{{
    for (int k=0; k<q; k++){
      basis_states=build_rotationally_symmetric_base_states_compact(q, k);
      int size_space=basis_states.size();
      cmat U(size_space, size_space);
      U=0.;
      for (int i=0; i<size_space; i++){
        state_l=DecodeCompactRotationallySymetricBasisState(basis_states(i));
        //       cout << "Norma state_l=" << norm(state_l) << endl;
        //       abort(); 
        for (int j=0; j<size_space; j++){
          state_r=DecodeCompactRotationallySymetricBasisState(basis_states(j));
          //         cout << "norma(state_r)=" << norm(state_r) << endl;
          //         if (j==1) abort(); 
          U(i,j) =  dot(conj(state_l),state_r) ;
          //         cout << dot(conj(state_l),state_r) << endl;
        }
      }
      //     cout << U << endl;
      //     cout << "ortonormalidad del sector k=" << k << endl;
      error += norm(eye_c(size_space)-U); 
    }
    cout <<"error de la ortonormalidad en casa sector es " << error << endl;
    // }}}
    basis_states.set_size(0); // {{{ Total dimension must be ok
    ivec sizes_k(q);
    for (int k=0; k<q; k++){
//       tmp_basis =  build_rotationally_symmetric_base_states_compact(q, k);
      basis_states= build_rotationally_symmetric_base_states_compact(q, k);
      sizes_k(k)=basis_states.size();
    }
//     cout <<"dimension de cada sistema  " << sizes_k  << endl;
    error+=abs(sum(sizes_k)-pow_2(q));
    cout <<"error de la dimension total " << error << endl;
    // }}}
    basis_states.set_size(0); // {{{ Orthonormality in the whole hilbert space
    U=0.;

    for (int k=0; k<q; k++){
      basis_states= concat(basis_states,build_rotationally_symmetric_base_states_compact(q, k));
    }
    for (int i=0; i<d; i++){
      state_l=DecodeCompactRotationallySymetricBasisState(basis_states(i));
      //       cout << "Norma state_l=" << norm(state_l) << endl;
      //       abort(); 
      for (int j=0; j<d; j++){
        state_r=DecodeCompactRotationallySymetricBasisState(basis_states(j));
        //         cout << "norma(state_r)=" << norm(state_r) << endl;
        //         if (j==1) abort(); 
        U(i,j) =  dot(conj(state_l),state_r) ;
        //         cout << dot(conj(state_l),state_r) << endl;
      }
    }
    //     cout << U << endl;
    //     cout << "ortonormalidad del sector k=" << k << endl;
    error += norm(eye_c(d)-U); 
    cout <<"error de la ortonormalidad en todo el espacio es " << error << endl;
    // }}}
    // in each sector the eigenvectors are eigenvectors {{{
    for (int k=0; k<q; k++){
      basis_states=build_rotationally_symmetric_base_states_compact(q, k);
      int size_space=basis_states.size();
      for (int i=0; i<size_space; i++){
        state=DecodeCompactRotationallySymetricBasisState(basis_states(i));
        error+=norm(apply_rotation(state)-exp((2.*pi*Im*k)/q)*state);
//         cout << "k= "<< k << ", i=" << i << "... error=" << error << endl;
//         if( k==0 && i==0){
//           cout << norm(state) << endl; 
//           cout << norm(apply_rotation(state)) << endl; 
//           cout << norm(exp((2.*pi*Im*k)/q)*state) << endl; 
//           cout << "1." << state << endl; 
//           cout << "2." << apply_rotation(state) << endl;
//           cout << "3." << exp((2.*pi*Im*k)/q)*state << endl; 
//           cout << (apply_rotation(state))-(exp((2.*pi*Im*k)/q)*state) << endl;
// //           abort();
//         }
//         abort(); 
      }
      //     cout << U << endl;
      //     cout << "ortonormalidad del sector k=" << k << endl;
//       error += norm(eye_c(size_space)-U); 
    }
    cout <<"error de los eigenvalores  " << error << endl;
    // }}}
    basis_states.set_size(0); // {{{ Completeness in the whole hilbert space
    U=0.;
    for (int k=0; k<q; k++){
      basis_states= concat(basis_states,build_rotationally_symmetric_base_states_compact(q, k));
    }
    for (int i=0; i<d; i++){
      state_l=DecodeCompactRotationallySymetricBasisState(basis_states(i));
      U+=Proyector(state_l);
    }
    error += norm(eye_c(d)-U); 
    cout <<"error de la completes en todo el espacio es " << error << endl;
    // }}}
    std::cout << "error total=" << error << std::endl;
    //}}}
  } else if(option=="test_horizontal_proyector") {// {{{
    cvec Ppsi, state;
    double error=0.;
    // PArece que en todas las dimensiones funciona 
    //  ./test_spins -o test_horizontal_proyector --i1 3 --i2 3 --i3 5
    int nv=i2.getValue(), nh=i3.getValue(); int q= nv*nh; int d=pow_2(q);
    cout << "Test to see if projector works." << "Grid " << nh << "x" << nv  << endl;
    std::complex<double> eigen_phase;
    for (int k=0; k<nh; k++){
      state = RandomState(d);
      eigen_phase = exp(2.*itpp::pi*std::complex<double>(0,1)*(double(k)/nh));
      Ppsi = project_state_horizontal_momentum(k, state, nh);
      error += norm (eigen_phase*Ppsi - apply_horizontal_rotation(Ppsi, nh));
      error += abs(norm(Ppsi)-1);
    }
    cout << "Error = " << error << endl;
    //}}}
  } else if(option=="test_vertical_proyector") {// {{{
    cvec Ppsi, state;
    double error=0.;
    // PArece que en todas las dimensiones funciona 
    //  ./test_spins -o test_horizontal_proyector --i1 3 --i2 3 --i3 5
    int nv=i2.getValue(), nh=i3.getValue(); int q= nv*nh; int d=pow_2(q);
    cout << "Test to see if vertical projector works." << "Grid " << nh << "x" << nv  << endl;
    std::complex<double> eigen_phase;
    for (int k=0; k<nh; k++){
      state = RandomState(d);
      eigen_phase = exp(2.*itpp::pi*std::complex<double>(0,1)*(double(k)/nv));
      Ppsi = project_state_vertical_momentum(k, state, nh);
      error += norm (eigen_phase*Ppsi - apply_vertical_rotation(Ppsi, nh));
      error += abs(norm(Ppsi)-1);
    }
    cout << "Error = " << error << endl;
    //}}}
  } else if(option=="test_create_base_2d") {// {{{
    Array<CompactSymmetricBaseMember> basis_states, tmp_basis, basis_states_many_body;
    Array<cvec> statesbasic;
    CompactSymmetricBaseMember g;
    cvec state_l, state_h, state_r, state, prestate;
    int q, d;
    double error=0.;

    int nv=3, nh=4; q= nv*nh;d=pow_2(q);
    int tv_sector=0;

    // Generacion de los estados base con los que se construira el resto de cosas. {{{
    // Aca genero los estados que llamare phi_i con i=0,...,15 (2^nh-1) Esos son los que serviran de base
    // para construir los otros. E
    basis_states=build_rotationally_symmetric_base_states_compact(nh);
    // Ahora quiero construir uno en particular, del sector con momento vertical igual a 1.
    basis_states_many_body=build_rotationally_symmetric_base_states_compact(q);
    for (int i=0; i<basis_states_many_body.size(); i++){
//       cout << i <<"; " << basis_states_many_body(i) << endl;
    }
    // }}}
    g=basis_states_many_body(31);
    cout << g << endl;


    int mask;
    ivec  horizontal_basis_state_numbers(nv);
    for (int i_row=0; i_row < nv; i_row++){
      mask = ((pow_2(nh)-1)<<(nh*i_row)); 
//       cout << "Generator=" << g.generator <<", Mask=" << mask << ", number=" << ( (g.generator & mask) >> (nh*i_row)) << endl; 
      horizontal_basis_state_numbers(i_row) = (g.generator & mask) >> (nh*i_row);

    } // cout << "Smaller generators = " << horizontal_basis_state_numbers << endl;
//     for (int ib=0; ib< pow_2(nh) ; ib++){
//       std::cout << "basis small =" << basis_states(ib) << std::endl;
//     }
    
    statesbasic.set_size(nv);
    int total_k_horizontal=0;
    std::complex<double> phase;
    for (int  i_row=0; i_row < nv; i_row++){
      statesbasic(i_row)=DecodeCompactRotationallySymetricBasisState(basis_states(horizontal_basis_state_numbers(i_row)));
      total_k_horizontal += basis_states(horizontal_basis_state_numbers(i_row)).k; 
    }
    total_k_horizontal = total_k_horizontal % nh;
    phase = exp(2.*itpp::pi*std::complex<double>(0,1)*(double(total_k_horizontal)/nh));
    prestate=TensorProduct(statesbasic);
    cout << "Test if eigenstate of horizontal rotation " 
      << "Proportionality constant " << proportionality_constant(prestate, apply_horizontal_rotation(prestate, nh)) 
      << " (error=" << proportionality_test(prestate, apply_horizontal_rotation(prestate, nh)) << ")" << endl; 
    cout << "A ver " << norm(phase*prestate - apply_horizontal_rotation(prestate, nh)) << endl; 

//     Aca voy, ahora toca hacer las traslaciones en forma vertical. Ver que el momento total de lo anterior es igual a la suma de momentos. 


    cout << "See that the general projection operator works the same on basis states." << endl; 
    std::cout << "error total=" << error << std::endl;

    cvec bi_symmetric_state; 

    bi_symmetric_state = project_state_vertical_momentum(tv_sector, prestate, nh);
    cout << "* Test if eigenstate of horizontal rotation " 
      << "Proportionality constant " << proportionality_constant(bi_symmetric_state, apply_horizontal_rotation(bi_symmetric_state, nh)) 
      << " (error=" << proportionality_test(bi_symmetric_state, apply_horizontal_rotation(bi_symmetric_state, nh)) << ")" << endl; 
    cout << "* Test if eigenstate of vertical rotation " 
      << "Proportionality constant " << proportionality_constant(bi_symmetric_state, apply_vertical_rotation(bi_symmetric_state, nh)) 
      << " (error=" << proportionality_test(bi_symmetric_state, apply_vertical_rotation(bi_symmetric_state, nh)) << ")" << endl; 

    cout << "Test to see if projector works" << endl;
    int k=0;
    state = RandomState(d);
    phase = exp(-2.*itpp::pi*std::complex<double>(0,1)*(double(k)/nh));
    state_h = project_state_horizontal_momentum(k, state, nh);
    cout << norm (state_h - apply_horizontal_rotation(state_h, nh))  << ", " << norm(state_h) << endl;
    //}}}
  } else {// {{{
    vec b(3);
    b(0)=-0.404;b(1)=0.0844;b(2)=0.361;
    cout <<   exp_minus_i_b_sigma(b) << endl;
    // b = {-0.404, 0.084, 0.361};
    // InputForm[MatrixExp[-I b.Table[Pauli[i], {i, 3}]]]
    // {{0.8534308186484069 - 0.34318420526987775*I, 
    // -0.07985449651709063 + 0.3840621022964837*I}, 
    // {0.0798544965170907 + 0.3840621022964837*I, 
    // 0.8534308186484068 + 0.3431842052698777*I}}
  }//}}}
//}}}
  return 0;
}//}}}
