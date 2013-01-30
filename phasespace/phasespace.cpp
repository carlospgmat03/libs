// include files {{{
// #include <iostream>
#ifndef PHASESPACE
#define PHASESPACE
// #include <itpp_ext_math.h>
// #include <cfp_math.cpp>
// #include <purity_RMT.cpp>
#include <itpp/itbase.h>
// #include <itpp/stat/misc_stat.h>
// }}}

extern "C" { 
// It seems that a good deal of information can be obtained through, say 
// nm torus_mac.o
// 
// Algunos ejemplos de como llamar ritunas de fortran en c++ estan en las librerias de cpp, de ahi se llama
// desde c++ con las rutinas nativas de itpp. Los archivos utiles son 
// /home/carlosp/investigacion/nacho/itppstrip/itpp-4.2/itpp/base/algebra/eigen.cpp
// /home/carlosp/investigacion/nacho/itppstrip/itpp-4.2/itpp/base/algebra/lapack.h
void __phase_space_routines_MOD_wigner_from_state_tmp(std::complex<double> *phi); 
void __phase_space_routines_MOD_wigner_from_state(int& n, std::complex<double> *phi, double *wigner); 
}

namespace phasespace{ // {{{ .h
} // }}}
namespace phasespace{ // {{{
  itpp::mat wigner(itpp::cvec& state){ // {{{
    int n=state.size();
      itpp::mat w(2*n, 2*n);
//       __phase_space_routines_MOD_wigner_from_state_tmp(state._data());
//     std::cout << "in wigner, n=" << n << std::endl;
//       __phase_space_routines_MOD_wigner_from_state(n, state._data());
      __phase_space_routines_MOD_wigner_from_state(n, state._data(), w._data());
//       std::cout <<"Maricon" << w << std::endl;

      return w;
  } //}}}
} // }}}
#endif // PHASESPACE

