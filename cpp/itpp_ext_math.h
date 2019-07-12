#ifndef ITPP_EXT_MATH_H
#define ITPP_EXT_MATH_H
#include <itpp/itbase.h>
namespace itppextmath{ // {{{ .h
template <class Num_T>  double test_symmetry(const itpp::Mat<Num_T >);
double test_real_symmetric(const itpp::cmat);
itpp::cmat sigma(int );
itpp::cvec OperatorToVectorPauliBasis(itpp::cmat );
itpp::cmat multiply_by_sigma_leftmost_qubit(const itpp::cmat&, int);
itpp::ivec IntegerDigits(int, int, int);
itpp::Vec<bool> int_to_bool(const int, const int );
int bvec_to_int(itpp::Vec<bool>);
} // }}}


#endif
