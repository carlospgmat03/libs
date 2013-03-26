#ifndef ITPP_EXT_MATH_H
#define ITPP_EXT_MATH_H
#include <itpp/itbase.h>
namespace itppextmath{ // {{{ .h
template <class Num_T>  double test_symmetry(const itpp::Mat<Num_T >);
double test_real_symmetric(const itpp::cmat);
itpp::cmat sigma(int );
itpp::cvec OperatorToVectorPauliBasis(itpp::cmat );
itpp::cmat multiply_by_sigma_leftmost_qubit(const itpp::cmat& A, int sigma_label);
} // }}}


#endif
