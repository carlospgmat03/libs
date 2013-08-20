// {{{ includes and using
#include <iostream>
#include <tclap/CmdLine.h>
// #include <RMT.cpp>
// #include <purity_RMT.cpp>
#include "itpp_ext_math.cpp"
#include "spinchain.cpp"
#include "cuda_routines.cu"

using namespace std;
// using namespace itpp;
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
    itpp::cvec state(dim);
    double x,y;
    ifstream myReadFile;
    myReadFile.open("/tmp/estado.dat");
    for (int i=0; i<dim; i++){
      myReadFile >> x >> y ;
      state(i) = complex<double>(x,y) ;
    }
    myReadFile.close();
    
    itpp::vec b(3); b(0)=bx.getValue(); b(1)=by.getValue(); b(2)=bz.getValue();
    apply_magnetic_kick(state,b,position.getValue());
    for (int i=0; i<dim; i++){
      cout << real(state(i)) << " " << real(-Im*state(i)) << endl;
    }//}}}
  } else if(option=="test_project_base_state") {// {{{
    itpp::Array<CompactSymmetricBaseMember> basis_states, tmp_basis;
    itpp::cvec state_l, state_r, state;
    int q=qubits.getValue();
    int d=pow_2(q);itpp::cmat U(d,d);
    double error=0.;
    // orthonormality in each sector {{{
    for (int k=0; k<q; k++){
      basis_states=build_rotationally_symmetric_base_states_compact(q, k);
      int size_space=basis_states.size();
      itpp::cmat U(size_space, size_space);
      U=0.;
      for (int i=0; i<size_space; i++){
        state_l=DecodeCompactRotationallySymetricBasisState(basis_states(i));
        //       cout << "Norma state_l=" << norm(state_l) << endl;
        //       abort(); 
        for (int j=0; j<size_space; j++){
          state_r=DecodeCompactRotationallySymetricBasisState(basis_states(j));
          //         cout << "norma(state_r)=" << norm(state_r) << endl;
          //         if (j==1) abort(); 
          U(i,j) =  itpp::dot(conj(state_l),state_r) ;
          //         cout << itpp::dot(conj(state_l),state_r) << endl;
        }
      }
      //     cout << U << endl;
      //     cout << "ortonormalidad del sector k=" << k << endl;
      error += norm(itpp::eye_c(size_space)-U); 
    }
    cout <<"error de la ortonormalidad en casa sector es " << error << endl;
    // }}}
    basis_states.set_size(0); // {{{ Total dimension must be ok
    itpp::ivec sizes_k(q);
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
      basis_states= itpp::concat(basis_states,build_rotationally_symmetric_base_states_compact(q, k));
    }
    for (int i=0; i<d; i++){
      state_l=DecodeCompactRotationallySymetricBasisState(basis_states(i));
      //       cout << "Norma state_l=" << norm(state_l) << endl;
      //       abort(); 
      for (int j=0; j<d; j++){
        state_r=DecodeCompactRotationallySymetricBasisState(basis_states(j));
        //         cout << "norma(state_r)=" << norm(state_r) << endl;
        //         if (j==1) abort(); 
        U(i,j) =  itpp::dot(conj(state_l),state_r) ;
        //         cout << itpp::dot(conj(state_l),state_r) << endl;
      }
    }
    //     cout << U << endl;
    //     cout << "ortonormalidad del sector k=" << k << endl;
    error += norm(itpp::eye_c(d)-U); 
    cout <<"error de la ortonormalidad en todo el espacio es " << error << endl;
    // }}}
    // in each sector the eigenvectors are eigenvectors {{{
    for (int k=0; k<q; k++){
      basis_states=build_rotationally_symmetric_base_states_compact(q, k);
      int size_space=basis_states.size();
      for (int i=0; i<size_space; i++){
        state=DecodeCompactRotationallySymetricBasisState(basis_states(i));
        error+=norm(apply_rotation(state)-exp((2.*itpp::pi*Im*double(k))/double(q))*state);
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
//       error += norm(itpp::eye_c(size_space)-U); 
    }
    cout <<"error de los eigenvalores  " << error << endl;
    // }}}
    basis_states.set_size(0); // {{{ Completeness in the whole hilbert space
    U=0.;
    for (int k=0; k<q; k++){
      basis_states= itpp::concat(basis_states,build_rotationally_symmetric_base_states_compact(q, k));
    }
    for (int i=0; i<d; i++){
      state_l=DecodeCompactRotationallySymetricBasisState(basis_states(i));
      U+=Proyector(state_l);
    }
    error += norm(itpp::eye_c(d)-U); 
    cout <<"error de la completes en todo el espacio es " << error << endl;
    // }}}
    std::cout << "error total=" << error << std::endl;
    //}}}
  } else if(option=="test_horizontal_proyector") {// {{{
    itpp::cvec Ppsi, state;
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
    itpp::cvec Ppsi, state;
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
    itpp::Array<CompactSymmetricBaseMember> basis_states, tmp_basis, basis_states_many_body;
    itpp::Array<itpp::cvec> statesbasic;
    CompactSymmetricBaseMember g;
    itpp::cvec state_l, state_h, state_r, state, prestate;
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
    itpp::ivec  horizontal_basis_state_numbers(nv);
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

    itpp::cvec bi_symmetric_state; 

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
  } else if(option=="test_commutator") {// {{{
    itpp::cvec Ppsi, state;
    double error=0.;
    // PArece que en todas las dimensiones funciona 
    //  ./test_spins -o test_horizontal_proyector --i1 3 --i2 3 --i3 5
    int nv=i2.getValue(), nh=i3.getValue(); int q,  d;
    cout << "Test to see if vertical projector works." << "Grid " << nh << "x" << nv  << endl;
    itpp::vec b(3);
    b(0)=1.42; b(1)=2.534; b(2)=3.78;
    double J=1.2;
//     b=0.00001;
//     b(0)=2*itpp::pi; b(1)=0.; b(2)=0.;
//     J=itpp::pi/4;

    itpp::cvec state_0 , s1, s2, s3, s4, s5;
//     state_0=0.;
//     state_0(i1.getValue())=1.;

    for (int nh=2; nh<5; nh++){ for (int nv=2; nv<5; nv++){
    	q=nv*nh;
	d=pow_2(q);
	state_0 = RandomState(d);
	// Conmutador [U, Tv] {{{
	s1=apply_vertical_rotation(state_0, nh);
	itppcuda::apply_floquet( s1, J, b);
	s2=state_0;
	itppcuda::apply_floquet(s2, J, b);
	s3=apply_vertical_rotation(s2, nh);
	cout << "Error en [U, Tv]=" << norm (s3 - s1) << endl;
	error += norm (s3 - s1);
	// }}}
	// Conmutador [U, Th] {{{
	s1=apply_horizontal_rotation(state_0, nh);
	itppcuda::apply_floquet2d( s1, J, b, nh);
	s2=state_0;
	itppcuda::apply_floquet2d(s2, J, b, nh);
	s3=apply_horizontal_rotation(s2, nh);
	error += norm (s3 - s1);
	cout << "Error en [U, Th]=" << norm (s3 - s1) << endl;
	// }}}
    } }
    cout << "Error total = " << error << endl;
    //}}}
  } else {// {{{
    itpp::vec b(3);
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
