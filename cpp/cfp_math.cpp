#ifndef CFP_MATH_VARIOUS
#define CFP_MATH_VARIOUS

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "cfp_math.h"

#ifndef M_LOG2E
#define M_LOG2E 1.44269504088896340736 //log2(e)
#endif  

namespace cfpmath{// {{{

  bool even_number_of_ones_base_2(const int );
}// }}}
namespace cfpmath{// {{{
// http://stackoverflow.com/questions/3064926/how-to-write-log-base2-in-c-c
inline long double log2(const long double x){ // {{{
    return  log(x) * M_LOG2E;
} //}}}
double h_function(double x){// {{{
  if (x<=0.) return 0. ;
  if (x>=1.) return 0. ;
  return -x*log2(x);
}// }}}
double linear_interval(int iterator, int max_iterator, double min, double max){// {{{
  return min + iterator*(max-min)/(max_iterator-1);
}// }}}
double next_point_secant(double tnm1, double tnm2, double fnm1, double fnm2){// {{{
 return tnm1-fnm1*(tnm1-tnm2)/(fnm1-fnm2); 
  
}// }}}
// Bitwise functions and power 2 related// {{{
bool parity_sum_digits_base_2(const int number_in){// {{{
  //! This function returns the parity of the sum of the digits base 2
  /*! It is useful for fermions to find to which block (even or odd) a particular
   * base member 
   *  belongs to. 
   *  \param number_in Is the base number which we find if it belongs t
   *  o the even sector. 
   *
   *  Example:
   *  \verbatim
   parity_sum_digits_base_2(3)  = parity_sum_digits_base_2(101)  = 0
   parity_sum_digits_base_2(13) = parity_sum_digits_base_2(1101) = 1 \endverbatim
   */
  return !even_number_of_ones_base_2(number_in);
}// }}}
int BitCount(const unsigned int u) {// {{{
  // taken from http://tekpool.wordpress.com/category/bit-count/
  unsigned int uCount;

  uCount = u
    - ((u >> 1) & 033333333333)
    - ((u >> 2) & 011111111111);
  return
    ((uCount + (uCount >> 3))
     & 030707070707) % 63;
}// }}}
inline bool are_two_bits_on_one(const int n, const int position_1, const int position_2){// {{{
  if (! (n & pow_2(position_1)) ) return false;
  if (! (n & pow_2(position_2)) ) return false;
  return true;
}// }}}
inline int set_bit(const int a, const int bit){// {{{
  return (a | pow_2(bit));// debo usar el or
}// }}}
inline bool test_bit(const int a, const int bit){// {{{
  return (a & pow_2(bit));// debo usar el or
}// }}}
inline int remove_bit(const int n, const int bit_position){// {{{
  int n_first_bits = n & (pow_2(bit_position)-1);
  int n_second_bits= ((n - n_first_bits)  & ~pow_2(bit_position))>> 1;
  return n_first_bits+ n_second_bits;
}// }}}
inline int insert_bit(const int n, const int bit_position, const int bit){// {{{
  int n_first_nm1_bits = n & (pow_2(bit_position)-1);
  int n_with_0 = ((n - n_first_nm1_bits) << 1)+ n_first_nm1_bits;
  if (bit==0) {
    return n_with_0;
  } else if (bit==1) {
    return n_with_0 + (1<<bit_position);
  } else {
    std::cerr << "Error en insert_bit" << std::endl;
    abort();
  }
}// }}}
int merge_two_numbers(const int a, const int b, const int mask){// {{{
  //! Merge two integers bit by bit. 
  /*! In this routine we merge two numbers. It is useful for doing the tensor
   *  product. 'a' and 'b' are the input number wheares as usual mask
   *  indicates the position. 
   *  \param a Is the first number
   *  \param b The second number to be merged 
   *  \param mask the number enconding the places of a
   *
   *  Example:
   *  \verbatim
   *  mask               = 0 0 1 0 1 0 0 1 = 41
   *  a                  =     1   0     1 = 5
   *  b                  = 1 0   0   1 0   = 18
   *  merge_two_numbers  = 1 0 1 0 0 1 0 1 = 165
   *   \endverbatim
   */
  int n1=a, n2=b, j=0, tmp=0;
  for (;!( n1==0 && n2==0 );){
    if(test_bit(mask,j)){
      if (test_bit(n1,0)) tmp=set_bit(tmp,j); 
      n1 >>= 1;
    }else{
      if (test_bit(n2,0)) tmp=set_bit(tmp,j); 
      n2 >>= 1;
    }
    j++;
  }
  return tmp;
}// }}}
void extract_digits(int n, int& n1, int& n2, int which){// {{{
  //   This routine takes numin and puts two numbers
  //   n1out and n2out which result from the digits
  //   of numin that are marked with the 1 bits
  //   of the number nwhich
  //   exambple
  //   nwhich=   0 1 0 1 0 0 1 = 41
  //   numin=    0 1 1 0 1 1 1 = 55
  //   n1out=      1   0     1 = 5
  //   n2out=    0   1   1 1   = 7
  n1=0; n2=0;
  int counter_1=0, counter_2=0;        
  int nrot=n;
  for (int j=0; nrot != 0 ; j++){
    if(test_bit(which,j)){
      if (test_bit(n,j)) n1=set_bit(n1, counter_1);
      counter_1++;
    } else {
      if (test_bit(n,j)) n2=set_bit(n2, counter_2);
      counter_2++;
    }
    nrot= nrot>>1; 
  }
}// }}}
void extract_digits(int n, int length, int& n1, int& n2, int which){// {{{
  //   This routine takes numin and puts two numbers
  //   n1out and n2out which result from the digits
  //   of numin that are marked with the 1 bits
  //   of the number nwhich
  //   exambple
  //   nwhich=   0 1 0 1 0 0 1 = 41
  //   numin=    0 1 1 0 1 1 1 = 55
  //   n1out=      1   0     1 = 5
  //   n2out=    0   1   1 1   = 7
  n1=0; n2=0;
  int counter_1=0, counter_2=0;        
  for (int j=0; j<length; j++){
    if(test_bit(which,j)){
      if (test_bit(n,j)) n1=set_bit(n1, counter_1);
      counter_1++;
    } else {
      if (test_bit(n,j)) n2=set_bit(n2, counter_2);
      counter_2++;
    }
  }
}// }}}
bool even_number_of_ones_base_2(const int number_in){// {{{
  //! This function returns true if the number has an even number of ones in its representation in base 2.
  /*! It is useful for fermions to find to which block (even or odd) a particular base member 
   *  belongs to. 
   *  \param number_in Is the base number which we find if it belongs to the even sector. 
   *
   *  Example:
   *  \verbatim
   even_number_of_ones_base_2(3)  = even_number_of_ones_base_2(101)  = true
   even_number_of_ones_base_2(13) = even_number_of_ones_base_2(1101) = false \endverbatim
   *  see also http://www.dreamincode.net/code/snippet761.htm
   */
  bool tmp=true;
  int  n=number_in;
  while(n != 0){
    if(n & 1) tmp =!tmp ; 
    n >>= 1 ;
  }
  return tmp;
}// }}}
bool same_parity_sum_digits_base_2(const int n1, const int n2){// {{{
  return (even_number_of_ones_base_2(n1) == 
      even_number_of_ones_base_2(n2));
}// }}}
int pow_2(const int n){// {{{
  //! This function returns 2^n
  /*! \param n the number to which 2 is going to be powered
   *
   *  Example:
   *  \verbatim
   pow_2(3)  = 8 \endverbatim
   */
  if (n==0) return 1;
  return 1 << n;
}// }}}
bool is_integer_power_2(unsigned i){ // {{{
  // http://www.cprogramming.com/snippets/show.php?tip=10&count=30&page=0
  return !((i-1) & i); 
}// }}}
int log_base_2(const int n){// {{{
  //     if (! is_integer_power_2(n)){std::cout<<"Error en log_base_2";abort();}
  if (! is_integer_power_2(n)){std::cout<<"Error en log_base_2\n";abort();}
  int tmp=0, n_tmp=n;
  while(tmp<50){
    if (n_tmp==1) return tmp;
    tmp++;
    n_tmp>>=1;
  }
  abort();
}// }}}
void test_log_base_2(){// {{{
  for (int i=1; i<9999999; i*=2){
    std::cout << "log_2 ("<<i <<")="<< log_base_2(i) << "\n";
  }
}// }}}
bool parity_order_inversion_base_2(const int number_in){// {{{
  //! This function returns the number of 1,1 swaps required to invert a number
  /*! Consider a number 29 =0 1 1 1 0 1
   *                        a b c d e f
   *
   *  If we want to invert the cables, we would obtain a number of intersections:
   *
   *                        f e d c b a
   *                        1 0 1 1 1 0
   *  So there is a number of crossings. The question is how many ones appear in the 
   *  crossing? For example, consider 4 bits:
   *
   *  0 0 0 0 => Require no crossing (Even parity)
   *  0 0 0 1 => Require no crossing (Even parity)
   *  0 0 1 0 => Again no crossing   (Even parity)
   *  0 0 1 1 => Require one crossing(Odd parity)
   */
  return ((cfpmath::BitCount(number_in)%4)>1); 
}// }}}
// 
// 
// 
// 
// 
// 
// 
// 
// 
// 
// int rotate_bits(int number_in, int starting_digit, int ending_digit, int power){ // {{{
//   // We want to rotate an "inner" part of the number. Consider
//   // number_in = 0  0  1  1  0  1  0  1  1  0  0 = 428
//   // starting digit =3
//   // ending digit = 8
//   // Digits to rotate 
//   // number_in = 0  0  1  1  0  1  0  1  1  0  0 = 428
//   // digits    = 10 9  8  7  6  5  4  3  2  1  0 
//   // number    =       x  x  x  x  x  x          = mask
//   //           =       1  1  0  1  0  1          = 53
//   //  power = 3, so we have to rotatet it to the left 3 times. 
//   //  rotated bits =   1  0  1  1  1  0          = 46
//   // new number= 0  0  1  0  1  1  1  0  1  0  0 = 372
//   int n_fixed, n_to_rotate, n_rotated;
//   int mask=pow_2(ending_digit+1)-pow_2(ending_digit);
//   extract_digits(number_in, 
// } // }}}
int rotate_bits(int number_in, int size_register){ // {{{
// We rotate the bits of a number, but only the first "size_register" bits
int mask=(1<<(size_register+1))-1;
return 
(number_in&(~mask)) + ((number_in&mask)/2) + ((number_in%2)<<(size_register-1));
} // }}}
int rotate_bits(int number_in, int size_register, int power){ // {{{
// We rotate the bits of a number, but only the first "size_register" bits
int n_out=number_in;
return n_out;
} // }}}
int primitive_period_bit_rotation(int n, int size_register){ // {{{
// Calculates the minumum J such that
// T^j n = n
// where T is the "rotate_bits" function, n is an 
// integer and size_register is the size of the representation
//
int primitive_period=1, n_rotated=rotate_bits(n, size_register);
while (n_rotated != n){
  n_rotated=rotate_bits(n_rotated, size_register);
  primitive_period++;
}
return primitive_period;

} // }}}
int reverse_bits(int number_in, int size_register){ // {{{
  // We rotate the bits of a number, but only the first "size_register" bits
  int number_refelected=0;
  for (int i=0; i<size_register; i++){

    if(test_bit(number_in, i)){
      number_refelected=set_bit(number_refelected,size_register-i-1);
    }
  }
  return number_refelected;
} // }}}
//// }}}
// Number theoretical// {{{
  int maximum_prime_power_divisor(const int number, const int base){// {{{
    int power=0;
    while (number%integer_pow(base,power)==0){
      power++;
//       std::cout << "base="<< base << ", number="<< number << ", integer_pow(base,power)=" << integer_pow(base,power) << std::endl;
    }
    return power-1;
  }// }}}
  bool is_integer_power(unsigned i, unsigned base){ // {{{
    // http://www.cprogramming.com/snippets/show.php?tip=10&count=30&page=0
    return !((i-1) & i); 
  }// }}}
  int coarse_grain(const int number, const int coarse_grain_parameter){// {{{
    return number - number%coarse_grain_parameter;
  }// }}}
  int integer_pow(const int base, const int power){// {{{
    if (power==0 && base == 0){
      std::cout << "Algun lio en integer_pow weuwrtafds\n"; abort();
    }
    if (power < 0){
      std::cout << "Algun lio en integer_pow 385oosdiur\n"; abort();
    }
    int tmp=1;
    for (int i=0; i<power; i++){
      tmp*=base;
    }
    return tmp;
  }// }}}
  bool test_int_log(const int base, const int n){// {{{
    //! Test if an integer is a power of another. 
    /*! Check if \f$ n=base^{m} \f$ with \f$ m \f$ some   element of lhe even and odd lists are identical. 
     * 000,010. This elements are even (e) and odd(o). So the structure is completly determined by this 
     * par even odd.  * so the element we must consider is the upper corner, namely 0000,0101. The reduced elements are
     * 000,010. This elements are even (e) and odd(o). So the structure is completly determined by this 
     * par even odd. The following table documents the structure
     *
     *  \param position the position of the new qubit. 
     */
    int test=1, result=0;
    for (; test < n; test = base*test){ result++; }
    return (test==n);
//     if (test ==n){
//       return true;
//     } else {
//       return false;
//     }
  }// }}}
  int integer_part_log(const int base, const int n){// {{{
    int test=1, result=0;
//     std::cout << "integer_part_log " << base << ", " << n << std::endl;
    for (; test <= n; test = base*test){  result++ ; }
    return result;
  }// }}}
  int int_log(const int base, const int n){// {{{
    int test=1, result=0;
    for (; test < n; test = base*test){ result++; }
    if (test ==n){
      return result;
    } else {
      std::cerr << "Error en int_log. base=" << base << ", n=" << n << std::endl;
      abort();
    }
  }// }}}
  int translate_coupled_base_to_decoupled(const int &number_in){// {{{

//     return number_in /2; 
    return number_in >> 1; 
  }// }}}
  long isqrt (long x){// {{{
    /* Integer square root by Halleck's method, with Legalize's speedup */
    /* from: http://home.utah.edu/~nahaj/factoring/isqrt.c.html*/
    long   squaredbit, remainder, root;

    if (x<1) return 0;

    /* Load the binary constant 01 00 00 ... 00, where the number
     * of zero bits to the right of the single one bit
     * is even, and the one bit is as far left as is consistant
     * with that condition.)
     */
    squaredbit  = (long) ((((unsigned long) ~0L) >> 1) & 
        ~(((unsigned long) ~0L) >> 2));
    /* This portable load replaces the loop that used to be 
     * here, and was donated by  legalize@xmission.com 
     */

    /* Form bits of the answer. */
    remainder = x;  root = 0;
    while (squaredbit > 0) {
      if (remainder >= (squaredbit | root)) {
        remainder -= (squaredbit | root);
        root >>= 1; root |= squaredbit;
      } else {
        root >>= 1;
      }
      squaredbit >>= 2; 
    }

    return root;
  }// }}}
}// }}}
#endif // CFP_MATH_VARIOUS
