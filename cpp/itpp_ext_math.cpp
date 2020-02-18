// include files {{{
#ifndef  ITPP_EXT_MATH_VARIOUS
#define ITPP_EXT_MATH_VARIOUS
#include "itpp_ext_math.h"
#include "cfp_math.cpp"
#include <algorithm>    // std::next_permutation, std::sort

// #include <purity_RMT.cpp>
#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
// }}}
namespace itppextmath{ // {{{
  itpp::mat rotation_matrix(double theta){ // {{{
    itpp::mat R(2,2);
    R(0,0)=cos(theta);
    R(1,1)=R(0,0);
    R(1,0)=sin(theta);
    R(0,1)=-R(1,0);
    return R;
  } //}}}
  double sum_positive_derivatives(const itpp::vec& f){ // {{{
    double tmp=0;
    for (int i=0; i<f.size()-1; i++){
      if(f(i+1)>f(i)){
        tmp+=f(i+1)-f(i);
      }
    }
    return tmp;
  } // }}}
  // number theoretical
template <class Num_T> int counter_elements_close_fixed_value(itpp::Vec<Num_T> v, Num_T x, double epsilon){ //{{{
int counter=0;
for (int i=0; i< v.size(); i++){
  if( abs(v(i)-x) < epsilon) counter++;
}
return counter;

} // }}}
template <class Num_T> itpp::bmat zero_non_zero(itpp::Mat<Num_T> m, double epsilon){ //{{{
  int mr =  m.rows();
  int mc =  m.cols();
  itpp::bmat z(mr,mc);
  z=1;
//   std::cout << "hola (mr,mc)=(" << mr << "," << mc << ")" << std::endl;
  for (int i=0; i< mr; i++){
    for (int j=0; j< mc; j++){
//       std::cout << "hola " << i << ", " << j << std::endl;
      if(abs(m(i,j))< epsilon) z(i,j)=0;
    }
  }
return z;
} // }}}
// Inquiry {{{
template <class Num_T> int locateLargestAbs(const itpp::Vec<Num_T>& e){ // {{{
  int position=0;
  for (int i=0; i<e.size(); i++){
    if (abs(e(i)) > abs(e(position))){
      position=i;
    }
  }
  return position;
} // }}}
double compare(const itpp::cmat& A, const itpp::cmat& B){ // {{{
  return itpp::norm(A - B);
} //}}}
template <class Num_T> bool Contains(const itpp::Mat<Num_T>& v, const itpp::Vec<Num_T>& e){// {{{
  for (int i=0; i<v.rows();i++){ if (e==v.get_row(i)) return true; }
  return false;
} //}}}
template <class Num_T> bool Contains(const itpp::Array<Num_T >& v, const Num_T& e){// {{{
  for (int i=0; i<v.size();i++){ if (e==v(i)) return true; }
  return false;
} //}}}
template <class Num_T> bool Contains(const itpp::Array<Num_T >& BigArray, const itpp::Array<Num_T >& SmallArray){// {{{
  for (int i=0; i<SmallArray.size();i++){ if ( !Contains(BigArray, SmallArray(i))) return false; }
  return true;
} //}}}
template <class Num_T> bool Contains(const itpp::Vec<Num_T >& v, const Num_T& e){// {{{
  for (int i=0; i<v.size();i++){ if (e==v(i)) return true; }
  return false;
} //}}}
template <class Num_T> bool AreEqual(const itpp::Array<Num_T >& a1, const itpp::Array<Num_T >& a2){// {{{
//! Routine to check if two arrays are equal.
/*! It must be superseeded by the operator ==, but I still dont know how
 */
  if (a1.size() != a2.size()) return false;
  for (int i=0; i<a1.size();i++){ if (a1(i)!=a2(i)) return false; }
  return true;
} //}}}
template <class Num_T> bool AreAllEqual(const itpp::Vec<Num_T >& v){// {{{
  if (v.size() == 0) return true;
  Num_T e=v(0);
  for (int i=1; i<v.size();i++){ if (e!=v(i)) return false; }
  return true;
} //}}}
bool DoTheyIntersect(const itpp::Array< itpp::ivec >& Array1, const itpp::Array< itpp::ivec >& Array2){// {{{
//! Do intersection of two arrays is empty.
/*! It treat them as sets.
 */
  for (int i=0; i<Array1.size(); i++){
    if ( itppextmath::Contains(Array2, Array1(i) )) {return true;}
  }
  return false;
} //}}}
double test_similar_complex_list(itpp::cvec& z1, itpp::cvec& z2){ // {{{
  if (z1.size() != z2.size()){
    std::cerr << "En test_similar_complex_list, tamaños diferentes" << std::endl;
    return z1.size() + z2.size();
  }
  double error=0.;//, e;
  int dim = z1.size(), min_pos; 
  itpp::cvec zt = z2;
  itpp::Vec<bool> checked(dim);
  itpp::vec errors;
  for (int i=0; i<dim; i++){
    errors = abs(zt-z1(i));
//     min_index(errors);
    error += min(errors, min_pos);
//     error+=e;
    zt.del(min_pos);
  }
// std::complex<double> z;
  return error;

} // }}}
template <class Num_T> bool AreEqual_modulo_order(const itpp::Array<Num_T >& a1, const itpp::Array<Num_T >& a2){// {{{
//! Routine to check if two arrays are equal.
/*! It must be superseeded by the operator ==, but I still dont know how
 */
  if (a1.size() != a2.size()) return false;
  for (int i=0; i<a1.size();i++){ if ( !Contains(a1, a2(i))) {return false; }}
  return true;
} //}}}
template <class Num_T> int total_elements_depth_2(const itpp::Array<Num_T >& v){// {{{
  int result = 0;
  for (int i=0; i<v.size(); i++){
    result += v(i).size();
  }
  return result;
} //}}}
template <class Num_T> Num_T Last(const itpp::Array<Num_T >& v){// {{{
  return v(v.length()-1);
} //}}}
template <class Num_T> int Position_FirstIntersection(const itpp::Vec<Num_T >& v, const Num_T& e){// {{{
  for (int i=0; i<v.size();i++){ if (e==v(i)) return i; }
  std::cerr << "Not found in Position_FirstIntersection\n";
  abort();
  return -1;
} //}}}
template <class Num_T> int Position_FirstIntersection(const itpp::Array<Num_T >& v, const Num_T& e){// {{{
  for (int i=0; i<v.size();i++){ if (e==v(i)) return i; }
  std::cerr << "Not found in Position_FirstIntersection\n";
  abort();
  return -1;
} //}}}
template <class Num_T> double compare(const itpp::Array<Num_T >& a1, const itpp::Array<Num_T >& a2){// {{{
//! Routine to check if two arrays are equal.
/*! It must be superseeded by the operator ==, but I still dont know how
 */
  if (a1.size() != a2.size()) return -1;
  double x=0;
  for (int i=0; i<a1.size();i++){
    x+= compare(a1(i),a2(i));
  }
  return x;
} //}}}
template <class Num_T> Num_T proportionality_constant(const itpp::Vec<Num_T >& a1, const itpp::Vec<Num_T >& a2){// {{{
//! Attempt to get a proportionality constant between two vectors
/*!
 */
  if (a1.size() != a2.size()){
    std::cerr << "Not even wrong" << std::endl;
    abort();
  }
  int position_largest = locateLargestAbs(a1);
  return a1(position_largest)/a2(position_largest);
} //}}}
template <class Num_T> double proportionality_test(const itpp::Vec<Num_T >& a1, const itpp::Vec<Num_T >& a2, Num_T proportionality_constant){// {{{
//! Routine to check if two vectors are proportional.
/*!
 */
  if (a1.size() != a2.size()) return -1;
  return norm(a1-proportionality_constant*a2)/norm(a1);
} //}}}
template <class Num_T> double proportionality_test(const itpp::Vec<Num_T >& a1, const itpp::Vec<Num_T >& a2){// {{{
//! Routine to check if two vectors are proportional.
/*!
 */
  if (a1.size() != a2.size()) return -1;
  int position_largest = locateLargestAbs(a1);
//   std::cout << "position_largest " << position_largest <<", a1(p)=" << a1(position_largest) << std::endl;
  Num_T proportionality_constant=a1(position_largest)/a2(position_largest);
  return proportionality_test(a1, a2, proportionality_constant);
} //}}}
// }}}
// Reordering, replacing, inequalities {{{
itpp::cvec ReorderLastTwoSubsystems(itpp::cvec& psi, int n_middle, int n_last){ // {{{
  // the idea is to have a state psi_abc and to reorder it to get
  // the state psi_acb. n_middle is the dimension of the middle Hilbert space
  // and n_last the dimension of the last Hilbert space

  itpp::cvec psi_reorderd(psi.size());
  std::complex<double> tmp;
//   std::cout << "en ReorderLastTwoSubsystems. n_middle="<<n_middle<<", n_last="<<n_last << std::endl;
  int nb=n_middle, nc=n_last;
  int na=psi.size()/(nb*nc);
  int index_1, index_2;
  for (int ia=0; ia<na; ia++){ for (int ib=0; ib<nb; ib++){ for (int ic=0; ic<nc; ic++){

    index_1 = ia*nb*nc+ ib*nc+ ic;
    index_2 = ia*nb*nc+ ic*nb+ ib;
//     std::cout << "["<<ia<<ib<<ic<<"] Poniendo el original "
//     << index_1<<" en la posicion "<<index_2 << std::endl;
    psi_reorderd(index_2) = psi(index_1);
//     psi(index_1)= psi(index_2);
//     psi(index_2)= tmp;
  } } }
  return psi_reorderd;
} // }}}
void ReorderTwoSubsystems(itpp::cvec& psi, int n_last){ // {{{
  int n_middle=psi.size()/n_last;
  psi = ReorderLastTwoSubsystems(psi, n_middle, n_last);
  return;
} // }}}
bool vector_greater_than(const itpp::ivec& v, const itpp::ivec& w){ // {{{
  if (v.size() > w.size()){
    return true;
  } else if( v.size() < w.size()){
    return false;
  } else {
    for (int i=0; i< v.size(); i++){
      if ( v(i)>w(i)) {
        return true;
      } else if ( v(i)<w(i)) {
        return false;
      }
    }
    return false;
  }

} //}}}
template <class Num_T> itpp::Array<Num_T >  Replace(const itpp::Array<Num_T > v, Num_T old, Num_T nuevo){ // {{{
  itpp::Array<Num_T > tmp=v;
  for (int i=0; i<v.size(); i++){
    if (v(i)== old){
      tmp(i)=nuevo;
    }
  }
  //   return;
  return tmp;
} //}}}
itpp::Array<itpp::ivec> sort(const itpp::Array<itpp::ivec>& input){ // {{{
  //     http://en.wikipedia.org/wiki/Bubble_sort
  itpp::Array<itpp::ivec> tmp(0), result(0);
  itpp::ivec moco;
  bool swaped;
  for (int size=0; result.size() != input.size() ;size++){
    //       std::cerr << "size " << size <<"\n";
    for (int i = 0; i< input.size(); i++){
      if (input(i).size() == size){tmp = itpp::concat(tmp,input(i));}
    }
    swaped=true;
    //       std::cerr << " analizing size " << size <<"\n";
    while ( (swaped) && tmp.size() >1 ){
      swaped=false;
      //         std::cerr << " New pass " <<"\n";
      for (int i=0; i<tmp.size() -1; i++){
        //           std::cerr << "Analizing positions " << i <<" and "<<i+1 <<"\n";
        if ( vector_greater_than(tmp(i),tmp(i+1)) ){
          //             std::cerr << tmp(i) << " ?>?" << tmp(i+1)  <<"\n";
          moco=tmp(i); tmp(i)=tmp(i+1); tmp(i + 1)=moco;
          swaped=true;
        }
      }

    }
    result = itpp::concat(tmp,result);
    tmp.set_size(0);
  }
  return result;
} //}}}
template <class Num_T> itpp::Array<Num_T > Flatten(const itpp::Array<itpp::Array<Num_T > >& vct){ // {{{
  itpp::Array<Num_T > tmp(0);
  for (int i=0; i<vct.size(); i++){
    tmp=itpp::concat(tmp,vct(i));
  }
  return tmp;
} //}}}
template <class Num_T> itpp::Array<Num_T > Shuffle(const itpp::Array<Num_T >& vct){ // {{{
  // from http://www.algoblog.com/2007/06/04/permutation/
  itpp::Array<Num_T > v = vct;
  for (int i=1; i<v.size(); i++){ v.swap(i,itpp::randi(0,i)); }
  return v;
} //}}}
template <class Num_T> void swap(itpp::Vec<Num_T >& v, const int p1, const int p2){ // {{{
//   std::cout << "Entrando en el swap. v=" << v << ", p1=" << p1 << ", p2=" << p2 << std::endl;
  Num_T x;
  x=v(p1);
  v(p1)=v(p2);
  v(p2)=x;
//   std::cout << "saliendo  del swap. v=" << v << ", p1=" << p1 << ", p2=" << p2 << std::endl;
  return ;
} //}}}
template <class Num_T> itpp::Array<Num_T > CircularRotateRight(const itpp::Array<Num_T >& v){ // {{{
  itpp::Array<Num_T > tmp=v;
  tmp.shift_right(Last(tmp));
  return tmp;
} //}}}
template <class Num_T> itpp::Array<Num_T > CircularRotateRight(const itpp::Array<Num_T >& v, const int n){ // {{{
  int n_true=n%v.size();
  if (n_true == 0){return v;}
  if(n_true<0){n_true += v.size();}
  itpp::Array<Num_T > tmp=v;
  for (int i=0; i< n_true; i++){ tmp = CircularRotateRight(tmp); }
  return tmp;
} //}}}
template <class Num_T> itpp::Array<Num_T > CircularRotateRight(const itpp::Array<Num_T >& v, const int i_left, const int i_right){ // {{{
  itpp::Array<Num_T > tmp=v(i_left, i_right);
  tmp.shift_right(Last(tmp));
  return tmp;
} //}}}
template <class Num_T> itpp::Array<Num_T > CircularRotateLeft(const itpp::Array<Num_T >& v){ // {{{
  itpp::Array<Num_T > tmp=v;
  tmp.shift_left(v(0));
  return tmp;
} //}}}
itpp::cmat Reorder_state_tensor_form(itpp::cvec vector,int which){ // {{{
  //   std::cout << "@Reorder_state_tensor_form -2" << std::endl;
  int dim1 = cfpmath::pow_2(cfpmath::BitCount(which));
  //   std::cout << "@Reorder_state_tensor_form -1" << std::endl;
  int dim2 = vector.size()/dim1;
  //   std::cout << "@Reorder_state_tensor_form 0" << std::endl;
//   int qubits=cfpmath::log_base_2(vector.size());
  itpp::cmat out(dim2,dim1); out=0.;
  int col, row;
  //   std::cout << "@Reorder_state_tensor_form 1" << std::endl;
  for (int j=0; j<vector.size(); j++){
    cfpmath::extract_digits(j,col,row,which);
    out(row,col)=vector(j);
  }
  return out;
} // }}}
// }}}
// set manipulation and permutations {{{
// itpp::Array<itpp::ivec> all_permutations(vec p){
itpp::ivec next_permutation(const itpp::ivec& v){ // {{{

  int n = v.size();
  itpp::ivec w(n);
  int myints[n];
  for (int i=0; i<n; i++){
      myints[i] = v(i);
  } 
  std::next_permutation(myints,myints+n);
  for (int i=0; i<n; i++){
      w(i)=myints[i];
  } 

return w;
} // }}}
itpp::Array<itpp::ivec> SingleTwoPermutation(int i, int j){ // {{{
  itpp::Array<itpp::ivec> permutation(1);
  permutation(0)=itpp::vec_2(i,j);
  return permutation;
} // }}}
itpp::Array<itpp::ivec> ChangeNotationPermutation(const itpp::ivec& cycle){ // {{{
  // https://www.youtube.com/watch?v=MpKG6FmcIHk
  // el algoritmo comienza por 
//   cout << "En la rutina jodida  ChangeNotationPermutation cycle=" << cycle <<endl; 
  int q=cycle.size();
  itpp::Array<itpp::ivec> permutation_cycle_notation(0);
  itpp::Vec<bool> checked(q);
  checked = false;
  int current_position=0;

  itpp::ivec current_cycle;

//   cout << "2 En la rutina jodida  ChangeNotationPermutation " <<endl; 
  while ( current_position <q){
    current_cycle.set_size(0);
//     cout << "En el while. current_position=" << current_position << endl;
    while (!checked(current_position)){
      current_cycle = concat(current_cycle,current_position);
      checked(current_position) = true;
      current_position = cycle(current_position);
//       cout << "En el while while. current_position=" << current_position << endl;
    }
//     cout << "Sali del while interno. current_cycle=" << current_cycle<< endl;
    if(current_cycle.size()>1){
      permutation_cycle_notation = concat(permutation_cycle_notation, current_cycle);
    }
    current_position++;
  }
  return permutation_cycle_notation;
} // }}}
void all_permutations(int q){ // {{{
// for (int i=0, p.size(); i++){
// p(i) concatenado con las permutaciones de p(sin la i)
// }
// https://stackoverflow.com/questions/32678349/all-permutations-c-with-vectorint-and-backtracking
 int myints[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};

//         vector<int> v = { 1, 2, 3, 4 };
//   std::sort (myints,myints+3);

//   std::cout << "The 3! possible permutations with 3 elements:\n";
  do {
//     std::cout << myints << '\n';
    for (int i=0; i<q; i++){
      std::cout << myints[i] <<" " ;
    } 
    std::cout << std::endl;
//     std::cout << myints[0] << ' ' << myints[1] << ' ' << myints[2] << '\n';
  } while ( std::next_permutation(myints,myints+q) );

  std::cout << "After loop: " << myints[0] << ' ' << myints[1] << ' ' << myints[2] << '\n';

//   itpp::ivec v="1 2 3";
// //     std::vector<int> v = {1,2,3};
//     do {
// //         print(v);
// 
// 
//     std::cout << v << std::endl;
// //     for (int e : v) {
// //         std::cout << " " << e;
// //     }
// //     std::cout << std::endl;
// 
//     } while (std::next_permutation(v.begin(), v.end()));
return;
} // }}}
template <class Num_T> void interchange_permutations_2line(itpp::Vec<Num_T>& b, const itpp::ivec &p){ // {{{
  // This routine makes a permutation on vector b, coded with vector p. 
  // If 
  // b = [45 89 22 31 23 76]
  // and
  // p = [5 3 2 1 0 4]
  //
  // Then the element b(0)=45 will go to position p(0)=5
  // element b(1)=89 will go to position p(1)=3
  // etc
  // so the permuted b will be 
  // b -> [23 31 22 89 76 45]
  if(b.size() != p.size()){
    std::cerr << "interchange_permutations() in cfpmath: dimension mismatch" << std::endl;
    abort();
  }
  //Num_T temp;
  itpp::ivec q=p;
  int qtmp;
  for (int k = 0; k < q.size(); k++) {
    while (q(k) != k){
      qtmp=q(k);
      swap(q,k,qtmp);
      swap(b,k,qtmp);
    }
  }
} 
// }}}
template <class Num_T> void interchange_permutations_cycle(itpp::Vec<Num_T>& b, const itpp::ivec &c){ // {{{
  // This routine makes a permutation on vector b, coded with cycle c. 
  // If 
  // b = [45 89 22 31 23 76]
  // and
  // c = [4 1 2]
  //
  // b = [45 23 89 31 22 76]
  // so that what is in position c(i) goes to position c(i+1)
  // e.g. what is in position c(0)=4, that is p(c(0))=p(4) goes to position c(1)=1
  // e.g. what is in position c(1)=1,  goes to position c(2)=2
  // e.g. what is in position c(2)=2,  goes to position c(3)=c(0)=4
//   Num_T temp;
//   itpp::ivec q=p;
//   int qtmp;
  itpp::Vec<Num_T> b_sub=b.get(c);
  b_sub.shift_right(b_sub(b_sub.size()-1)) ;
  for (int k=0; k<c.size(); k++){
    b(c(k))=b_sub(k);
  }
} 
// }}}
int interchange_bits_2line(const int& n_in, const itpp::ivec &p){ // {{{
  // See documentaiton for 
  // void interchange_permutations_2line 
  // and for 
  // interchange_bits_cycle
  int n = n_in;
  itpp::Vec<bool> nbits=int_to_bool(n, p.size());
//   std::cout << "Control 1 " << nbits << ", " << c << std::endl;
  interchange_permutations_2line(nbits, p);
//   std::cout << "Control 2 " << std::endl;
  n = bvec_to_int(nbits);
//   std::cout << "Control 3 " << std::endl;
  return n;
} 
// }}}
int interchange_bits_cycle(const int& n_in, const itpp::ivec &c){ // {{{
  // This routine makes a permutation on integer n, coded with cycle c. 
  // If 
  // n = 16 = 010000 = [0 0 0 0 1 0 0 ...]
  // and
  // c = [3 4 1 2]
  //
  // b = [0 1 0 0 0 ... ] = 00010 = 2
  int n = n_in;
  int size = cfpmath::integer_part_log(2,n);
  if (max(c)+1 > size){ size=max(c)+1; 
//     std::cout << size << std::endl; 
  }
  itpp::Vec<bool> nbits=int_to_bool(n, size);
//   std::cout << "Control 1 " << nbits << ", " << c << std::endl;
  interchange_permutations_cycle(nbits, c);
//   std::cout << "Control 2 " << std::endl;
  n = bvec_to_int(nbits);
//   std::cout << "Control 3 " << std::endl;
  return n;
} 
// }}}
int interchange_bits_permutation(const int& n_in, const itpp::Array<itpp::ivec> &c){ // {{{
  // This routine makes a permutation on integer n, coded with cycle c. 
  // If 
  // n = 16 = 010000 = [0 0 0 0 1 0 0 ...]
  // and
  // c = [3 4 1 2]
  //
  // b = [0 1 0 0 0 ... ] = 00010 = 2
  int n = n_in;
  for (int i=0; i< c.size(); i++){
    n = interchange_bits_cycle(n, c(i));
  }
  return n;
} 
// }}}
template <class Num_T> itpp::Mat<Num_T > ArrayToMatrix(const itpp::Array<itpp::Vec<Num_T> >& A){// {{{
//   itpp::Mat<Num_T > tmp;
  itpp::Mat<Num_T > tmp(A.size(),A(0).size());
  for (int i=0; i<A.size(); i++){
    tmp.set_row(i,A(i));
  }
  return tmp; 
  
}// }}}
template <class Num_T> itpp::Mat<Num_T > OuterSum(const itpp::Vec<Num_T >& A, const itpp::Vec<Num_T >& B){// {{{
  itpp::Mat<Num_T > tmp(A.size(), B.size());
  for (int iA=0; iA<A.size(); iA++){
    for (int iB=0; iB<A.size(); iB++){
      tmp(iA, iB)=A(iA)+B(iB);
    }
  }
  return tmp; 
}// }}}
template <class Num_T> itpp::Array<Num_T > Outer(const itpp::Array<Num_T >& A, const itpp::Array<Num_T >& B){// {{{
  itpp::Array<Num_T > tmp(0);
  for (int i=0; i<A.size(); i++){
    tmp=itpp::concat(tmp, Elementwiseconcat(A(i),B));
  }
  return tmp; 
  
}// }}}
template <class Num_T> bool unique_elements(const Num_T & A){// {{{
  return A.size() == Union(A).size();
}// }}}
template <class Num_T> Num_T Union(const Num_T & A){// {{{
  if (A.size() == 0){
    return A;
  }
  Num_T tmp(1);
  tmp(0) = A(0);
  for (int i = 1; i< A.size(); i++){
    if (!Contains(tmp, A(i))){ tmp = itpp::concat(tmp, A(i)); }
  }
  return tmp;
}// }}}
template <class Num_T> itpp::Array<Num_T > Union(const itpp::Array<Num_T >& A, const itpp::Array<Num_T >& B){// {{{
  if (A.size() == 0){
    return Union(B);
  }
  itpp::Array<Num_T > tmp;
  tmp.set_size(1); 
  tmp(0) = A(0);
  for (int i = 1; i< A.size(); i++){
    if (!Contains(tmp, A(i))){ tmp = itpp::concat(tmp, A(i)); }
  }
  for (int i = 0; i< B.size(); i++){
    if (!Contains(tmp, B(i))){ tmp = itpp::concat(tmp, B(i)); }
  }
  return tmp;
}// }}}
template <class Num_T> itpp::Array<Num_T > del(const itpp::Array<Num_T >& A, const int i){// {{{
  int sA = A.size();
  if (i==0){
    if (sA==1){
      itpp::Array<Num_T > tmp(0);
      return tmp;
    } else {
      return A(1,sA-1);
    }
  } else if (i==sA -1) {
    return A(0,i-1);
  } else {
    return itpp::concat(A(0,i-1),A(i+1,sA-1));
  }
}// }}}
template <class Num_T> itpp::Array<Num_T > minus(const itpp::Array<Num_T >& A, const Num_T & B){// {{{
  itpp::Array<Num_T > A_Union=Union(A);
  if (!Contains(A_Union,B)){ 
    return A_Union;
  } else {
    int pi=Position_FirstIntersection(A_Union,B);
    return del(A_Union,pi);
  }
}// }}}
template <class Num_T> itpp::Array<Num_T > minus(const itpp::Array<Num_T >& A, const itpp::Array<Num_T >& B){// {{{
  itpp::Array<Num_T > A_union=Union(A), tmp;
  for (int i = 0; i< A_union.size(); i++){
    if (!Contains(B, A_union(i))){ tmp = itpp::concat(tmp, A_union(i)); }
  }
  return tmp ;
}// }}}
bool EmptyIntersection(const itpp::Array< itpp::ivec >& Array1, const itpp::Array< itpp::ivec >& Array2){// {{{
//! Stablish if 2 Arrays have an intersection
/*! It treat them as sets.  
 */
  for (int i=0; i<Array1.size(); i++){
    if ( itppextmath::Contains(Array2, Array1(i) )){
      return false;
    }
  } 
  return true;
}// }}}
itpp::Array< itpp::ivec > Intersection(const itpp::Array< itpp::ivec >& Array1, const itpp::Array< itpp::ivec >& Array2){// {{{
//! Get intersection of two arrays. 
/*! It treat them as sets.  
 */
  itpp::Array<itpp::ivec> tmp; 
  for (int i=0; i<Array1.size(); i++){
    if ( itppextmath::Contains(Array2, Array1(i) )){
      tmp=itpp::concat(tmp,Array1(i));
    }
  } 
  return Union(tmp);
}// }}}
template <class Num_T> itpp::Array<Num_T > RemoveItem(const itpp::Array<Num_T >& v, const int index){// {{{
  if (index < 0 || index > v.size()-1){
    std::cerr << "Something wrong in RemoveItem\n";
    abort();
  }
  if (v.size()==1){return itpp::Array<Num_T >(0);}
  if (index == 0){ return v( 1, v.size()-1);}
  if (index == v.size()-1){ return v(0, v.size()-2);}
  return itpp::concat(v(0, index -1), v(index + 1, v.size()-1));
}// }}}
template <class Num_T> itpp::Array<Num_T > RemoveItem(const itpp::Array<Num_T >& v, const itpp::ivec& indices){// {{{
  itpp::Array<Num_T > tmp = v;
  itpp::ivec indices_sorted = indices;
  itpp::sort(indices_sorted);
  for (int i =indices.size()-1; i>=0; i--){ tmp=RemoveItem(tmp,indices(i)); }
  return tmp;
}// }}}
itpp::Mat<int>  BulkPositionUnion2PositionSets(const itpp::Mat<int>& A, const itpp::Mat<int>& B){// {{{
  return itpp::concat_vertical(A,B);
}// }}}
// }}}
// Linear Algebra {{{
// Vectorization of density matrices. 
itpp::vec vectorization_density_matrix(const itpp::cmat& rho){ // {{{
  // We use the gell man matrices, ordered as in https://blog.suchideas.com/2016/04/sun-gell-mann-matrices-in-mathematica/
  int dim=rho.cols();
  itpp::vec rho_vec(dim*dim-1);
//   itpp::vec rho_vec(dim*dim);
  int vector_index=0;
  for (int i=1; i<dim; i++){
    for (int j=0; j<i; j++){
      rho_vec(vector_index)=2*real(rho(j,i));
      rho_vec(vector_index+(dim*(dim-1)/2))=-2*imag(rho(j,i));
//       cout << vector_index << ", " <<vector_index+(dim*(dim-1)/2) <<  endl;
      vector_index++;
    }
  }
  itpp::vec rho_diag = real(itpp::diag(rho));
  double x,y;

  y=sqrt((2.*(dim))/(dim+1));
  for (int k=1; k<dim; k++){
    y=-sqrt((2.*k)/(k+1));
    x=-y/k;
//     reconstruyendo 
//     std::cout << "k="<<k<<";  x="<<x<<";    y="<<y<<std::endl;
    rho_vec(k+dim*(dim-1)-1) = x*sum(rho_diag.get(0,k-1)) + rho_diag(k)*y;
  }
//   rho_vec(dim*dim-1)=
  return rho_vec;

} // }}}
itpp::cmat devectorization_density_matrix(const itpp::vec& rho_vec){ // {{{
  // We use the gell man matrices, ordered as in https://blog.suchideas.com/2016/04/sun-gell-mann-matrices-in-mathematica/
  int dim=cfpmath::isqrt(rho_vec.size()+1);
  std::complex<double> Im(0.,1.);
  itpp::cmat rho(dim,dim);
//   itpp::vec rho_vec(dim*dim);
  int vector_index=0;
  for (int i=1; i<dim; i++){
    for (int j=0; j<i; j++){
      rho(j,i)=(rho_vec(vector_index) - Im*rho_vec(vector_index+(dim*(dim-1)/2)))/2.;
      rho(i,j)=conj(rho(j,i));
//       rho_vec(vector_index)=2*real(rho(j,i));
//       rho_vec(vector_index+(dim*(dim-1)/2))=2*imag(rho(j,i));
//       cout << vector_index << ", " <<vector_index+(dim*(dim-1)/2) <<  endl;
      vector_index++;
    }
  }
//   itpp::vec rho_diag = real(itpp::diag(rho));
  double y, x, rho_parcial=0.;
  int r_index;

  for (int k=dim-1; k > 0;  k--){
    r_index = dim*dim - (dim-k)-1;
    x = sqrt(2./(k*k+k));
    y = -k*x;
//     std::cout << "k="<<k<<";  x="<<x<<";    y="<<y<<std::endl;
    rho(k,k) = (rho_vec(r_index) - x*( 1- rho_parcial))/(y-x);
    rho_parcial += real(rho(k,k));
//     std::cout << "r_index =" << r_index << std::endl;
//     rho_vec(k+dim*(dim-1)) = (sum(rho_diag.get(0,k)) - rho_diag(k+1))*y;
  }
  rho(0,0)=1-rho_parcial;
//   rho_vec(dim*dim-1)=
  return rho;

} // }}}
//
// Traces and Partial traces 
itpp::cmat partial_trace(const itpp::cvec& state, int dim_leave){ // {{{
  // This partial trace allows to leave a space of arbitrary dimension. Let
  // H = H_a H_b be th ehilbert space. "i" labels H_a and "j" the H_b
  // then
  // psi = \sum \alpha_ij |ij>
  // so
  // rho_AB = \sum_{iji'j'} \alpha_{ij} \alpha_{i'j'}^* |ij><i'j'|
  // rho_B  = \sum_{ji'j'} \alpha_{ij} \alpha_{ij'}^* |j><j'|
  itpp::cmat rho_B(dim_leave,dim_leave );
//   std::cout << "@partial_trace_qubits 2" << itpp::norm(reordered_state) <<std::endl;
  for (int j =0; j< dim_leave; j++){
    for (int jp =0; jp< dim_leave; jp++){
      rho_B(j,jp)=0;
      for (int i =0; i< state.size()/dim_leave; i++){
        rho_B(j,jp)+=conj(state(i*dim_leave +jp) )*state(i*dim_leave +j);
      }
    }
  }
//   std::cout << itpp::trace(rho_B) << std::endl;
//   abort();
  return rho_B;

} // }}}

/**
 * \brief Partial trace of a complex state.
 *
 * Calculates the partial trace of the specified composite state, obtaining
 * the subspace indicated with the second parameter.
 *
 * The position of the qubits that define the desired subspace is specified in
 * binary form in the variable `which`, for example, `0000011`, with `00000`
 * the elements of the environment and `11` are the elements of the desired
 * subspace, the variable is written in its decimal value.
 *
 * rho^A = tr_B (rho^AB)
 *
 * \param state [itpp::cvec] composite state in the form of complex vector.
 * \param which [int] position of desired subspace qubits, specified in binary, witten in decimal
 *
 * \return resulting matrix [itpp::cmat], with complex elements.
*/
itpp::cmat partial_trace_qubits(itpp::cvec state, int which){ // {{{
//   std::cout << "@partial_trace_qubits 1" << std::endl;
  int size_rho=cfpmath::pow_2(cfpmath::BitCount(which));
//   std::cout << "@partial_trace_qubits size_rho=" << size_rho
//     << ", which = " << which << std::endl;
  itpp::cmat reordered_state = Reorder_state_tensor_form(state, which), rho(size_rho,size_rho );
//   std::cout << "@partial_trace_qubits 2" << itpp::norm(reordered_state) <<std::endl;
  for (int i1 =0; i1< size_rho; i1++){
    for (int i2 =0; i2< size_rho; i2++){
      rho(i1,i2)=itpp::conj(reordered_state.get_col(i2))*reordered_state.get_col(i1);
    }
  }
//   std::cout << itpp::trace(rho) << std::endl;
//   abort();
  return rho;

} // }}}

/**
 * \brief Partial trace of a density complex matrix.
 *
 * Calculates the partial trace of the specified density matrix, obtaining
 * the subspace indicated with the second parameter.
 *
 * The position of the qubits that define the desired subspace is specified in
 * binary form in the variable `which`, for example, `0000011`, with `00000`
 * the elements of the environment and `11` are the elements of the desired
 * subspace, the variable is written in its decimal value.
 *
 * rho^A = tr_B (rho^AB)
 *
 * \param rho [itpp::cmat] density complex matrix.
 * \param which [int] position of desired subspace qubits, specified in binary, witten in decimal
 *
 * \return resulting matrix [itpp::cmat], with complex elements.
*/
itpp::cmat partial_trace_qubits(const itpp::cmat& rho, int which){ // {{{
//   std::cout << "@partial_trace_qubits 1" << std::endl;
  int DimensionB=cfpmath::pow_2(cfpmath::BitCount(which));

//   std::cout << "@partial_trace_qubits DimensionB " << DimensionB << std::endl;
  itpp::cmat rho_out(DimensionB, DimensionB);
  int TotalDimension=rho.cols();
  int DimensionA=TotalDimension/DimensionB;
  std::complex<double> tmp;
  int ab1, ab2;
  for (int b1 =0; b1< DimensionB; b1++){
    for (int b2 =0; b2< DimensionB; b2++){
      tmp=0.;
      for (int a=0; a<DimensionA; a++){
        ab1=cfpmath::merge_two_numbers(b1, a, which);
        ab2=cfpmath::merge_two_numbers(b2, a, which);
        tmp+=rho(ab1, ab2);
      }
      rho_out(b1,b2)=tmp;
    }
  }
//   std::cout << itpp::trace(rho) << std::endl;
//   abort();
  return rho_out;

} // }}}
std::complex<double> trAB(const itpp::cmat& A, const itpp::cmat& B){ // {{{
  // Calculates tr A.B
  //
  std::complex<double> tr=0;
  int ni=A.cols(); // 2
  int nj=A.rows(); // 1
  if (ni!=B.rows() || nj !=B.cols()){
    std::cerr << "Error en trAB" << std::endl;
    abort();
  }
  for (int i=0; i<ni; i++){
    for (int j=0; j<nj; j++){
      tr+=B(i,j)*A(j,i);
//       std::cout << A(i,j)*B(j,i) <<
    }
  }
  return tr;
} // }}}
// Tensor Products
template <class Num_T> itpp::Mat<Num_T> postpend_tensor_identity(const itpp::Mat<Num_T>& Matrix, const int dim_identity){ // {{{
  itpp::Mat<Num_T> tmp(Matrix.rows()*dim_identity, Matrix.cols()*dim_identity);
  tmp.zeros();
  Num_T x;
  int pos_cols, pos_rows;
  for (int i_rows=0; i_rows <  Matrix.rows(); i_rows++){
    for (int i_cols=0; i_cols <  Matrix.cols(); i_cols++){
      x=Matrix(i_rows, i_cols);
      pos_cols=dim_identity*i_cols;
      pos_rows=dim_identity*i_rows;
      for (int i_diag=0; i_diag<dim_identity; i_diag++){tmp(pos_rows+i_diag, pos_cols+i_diag)= x;}
    }
  }
  return tmp;
} // }}}
template <class Num_T> itpp::Mat<Num_T> extend_qubit_operator(const itpp::Mat<Num_T>& Matrix, const int encoded_nontrivial_positions, const int total_number_of_qubits){ // {{{

  itpp::Mat<Num_T> tmp(cfpmath::pow_2(total_number_of_qubits),cfpmath::pow_2(total_number_of_qubits));
  tmp=0.;
  int number_of_nontrivial_positions = cfpmath::BitCount(encoded_nontrivial_positions);
  int i_col_total, j_col_total;
  for (int i=0; i<cfpmath::pow_2(total_number_of_qubits-number_of_nontrivial_positions); i++) {
    for (int i_col=0; i_col<cfpmath::pow_2(number_of_nontrivial_positions); i_col++){
      i_col_total = cfpmath::merge_two_numbers(i_col, i, encoded_nontrivial_positions);
      for (int j_col=0; j_col<cfpmath::pow_2(number_of_nontrivial_positions); j_col++){
        j_col_total = cfpmath::merge_two_numbers(j_col, i, encoded_nontrivial_positions);
        tmp(i_col_total, j_col_total)=Matrix(i_col, j_col);
      }
    }
  }
  return tmp;

}
// }}}
template <class Num_T> itpp::Mat<Num_T> single_qubit_operator_tensor_identity(const itpp::Mat<Num_T>& Matrix, const int position_operator, const int total_number_of_qubits){ // {{{
  itpp::Mat<Num_T> tmp(cfpmath::pow_2(total_number_of_qubits),cfpmath::pow_2(total_number_of_qubits));
  int i0,i1;
  tmp=0.;
  for (int i=0; i<cfpmath::pow_2(total_number_of_qubits-1); i++) {
    i0=cfpmath::insert_bit(i, position_operator, 0);
    i1=cfpmath::insert_bit(i, position_operator, 1);
    tmp(i0,i0)=Matrix(0,0);
    tmp(i0,i1)=Matrix(0,1);
    tmp(i1,i0)=Matrix(1,0);
    tmp(i1,i1)=Matrix(1,1);
  }
  return tmp;

} // }}}
template <class Num_T> itpp::Mat<Num_T> prepend_tensor_identity(const itpp::Mat<Num_T>& Matrix, const int dim_identity){ // {{{
  int rows_mat= Matrix.rows(), cols_mat= Matrix.cols() ;
  itpp::Mat<Num_T> tmp(rows_mat*dim_identity, cols_mat*dim_identity);
  tmp.zeros();
  for (int i_diag=0; i_diag<dim_identity; i_diag++){
    tmp.set_submatrix(rows_mat*i_diag, cols_mat*i_diag, Matrix);
  }
  return tmp;
} // }}}

/**
 * \brief Tensor product of three matrices.
 *
 * Calculates the tensor product of three specified matrices, generating another
 * matrix with the result.
 *
 * Matrix1 ⊗ Matrix2 ⊗ Matrix3
 *
 * \param Matrix1 [itpp::Mat<Num_T>] matrix 1 with data type elements Num_T.
 * \param Matrix2 [itpp::Mat<Num_T>] matrix 2 with data type elements Num_T.
 * \param Matrix3 [itpp::Mat<Num_T>] matrix 3 with data type elements Num_T.
 *
 * \return resulting matrix [itpp::Mat<Num_T>], with same data type as in the parameters.
*/
template <class Num_T> itpp::Mat<Num_T> TensorProduct(const itpp::Mat<Num_T>& Matrix1, const itpp::Mat<Num_T>& Matrix2, const itpp::Mat<Num_T>& Matrix3){ // {{{
  return TensorProduct(TensorProduct(Matrix1,Matrix2),Matrix3);
} // }}}

/**
 * \brief Tensor product of two vectors.
 *
 * Calculates the tensor product of two specified vectors, generating another
 * vector with the result.
 *
 * Vec1 ⊗ Vec2
 *
 * \param Vec1 [itpp::Vec<Num_T>] vector 1 with data type elements Num_T.
 * \param Vec2 [itpp::Vec<Num_T>] vector 2 with data type elements Num_T.
 *
 * \return resulting vector [itpp::Vec<Num_T>], with same data type as in the parameters.
*/
template <class Num_T> itpp::Vec<Num_T> TensorProduct(const itpp::Vec<Num_T>& Vec1, const itpp::Vec<Num_T>& Vec2){ // {{{
  itpp::Vec<Num_T> tmp(Vec1.size()*Vec2.size());
  tmp.zeros();
  for (int i1=0; i1<Vec1.size(); i1++){
    for (int i2=0; i2<Vec2.size(); i2++){
//       std::cout << i1 << ", " << i2 << ", " <<i1*Vec1.size() + i2 << std::endl;
      tmp(i1*Vec2.size() + i2 )=Vec1(i1)*Vec2(i2);
//       std::cout << Vec1(i1)*Vec2(i2) << std::endl;
    }
  }
//       std::cout << "Hola cabron " << std::endl;
  return tmp;
} // }}}

/**
 * \brief Tensor product of given states.
 *
 * Calculates the tensor product of all the states specified in the array
 * received, making the product in the order in which these states are found
 * within the array.
 *
 * States(0) ⊗ States(1) ⊗ States(2) ⊗ ... ⊗ States(States.size())
 *
 * \param States [itpp::Array<itpp::Vec<Num_T>>] set of states in an array.
 *
 * \return resulting vector [itpp::Vec<Num_T>], with same data type as in the parameters.
*/
template <class Num_T> itpp::Vec<Num_T> TensorProduct(const itpp::Array<itpp::Vec<Num_T> >& States){ // {{{
  itpp::Vec<Num_T> tmp=States(0);
  for (int i=1; i<States.size(); i++){
    tmp=TensorProduct(tmp,States(i));
  }
  return tmp;
} // }}}

/**
 * \brief Tensor product of two vectors.
 *
 * Calculates the tensor product of two specified vectors, generating another
 * vector with the result. First vector is of type itpp::cvec (complex),
 * while the second is of type itpp::vec (double).
 *
 * Vec1 ⊗ Vec2
 *
 * \param Vec1 [itpp::cvec] vector 1 with data type elements complex<double>.
 * \param Vec2 [itpp::vec] vector 2 with data type elements double.
 *
 * \return resulting vector [itpp::cvec], with data type complex<double>.
*/
itpp::cvec TensorProduct(const itpp::cvec& Vec1, const itpp::vec& Vec2){ // {{{
  return TensorProduct(Vec1, to_cvec(Vec2));
  // return Vec1;
} // }}}

/**
 * \brief Tensor product of two matrices.
 *
 * Calculates the tensor product of two specified matrices, generating another
 * matrix with the result.
 * Uses the itpp::kron method.
 *
 * Matrix1 ⊗ Matrix2
 *
 * \param Matrix1 [itpp::Mat<Num_T>] matrix 1 with data type elements Num_T.
 * \param Matrix2 [itpp::Mat<Num_T>] matrix 2 with data type elements Num_T.
 *
 * \return resulting matrix [itpp::Mat<Num_T>], with same data type as in the parameters.
*/
template <class Num_T> itpp::Mat<Num_T> TensorProduct(const itpp::Mat<Num_T>& Matrix1, const itpp::Mat<Num_T>& Matrix2){ // {{{
//   cerr << "There is the itpp::kron function that does exactly what I wanted here" << endl ;
//   abort();
//   itpp::Mat<Num_T> tmp(Matrix1.rows()*Matrix2.rows(), Matrix1.cols()*Matrix2.cols());
//   tmp.zeros();
//   for (int i_row=0; i_row< Matrix1.rows(); i_row++){
//     for (int i_col=0; i_col< Matrix1.cols(); i_col++){
//       tmp.set_submatrix(i_row*Matrix2.rows(),i_col*Matrix2.cols(), Matrix1(i_row,i_col)*Matrix2);
//     }
//   }
  return itpp::kron(Matrix1, Matrix2);
} // }}}

/**
 * \brief Tensor product of given matrices.
 *
 * Calculates the tensor product of all the matrices specified in the array
 * received, making the product in the order in which these matrices are found
 * within the array.
 *
 * Matrices(0) ⊗ Matrices(1) ⊗ Matrices(2) ⊗ ... ⊗ Matrices(Matrices.size())
 *
 * \param Matrices [itpp::Array<itpp::Mat<Num_T>>] set of matrices in an array.
 *
 * \return resulting matrix [itpp::Mat<Num_T>], with same data type as in the parameters.
*/
template <class Num_T> itpp::Mat<Num_T> TensorProduct(const itpp::Array<itpp::Mat<Num_T> >& Matrices){ // {{{
  itpp::Mat<Num_T> tmp=Matrices(0);
  for (int i=1; i<Matrices.size(); i++){
    tmp=itpp::kron(tmp,Matrices(i));
  }
  return tmp;
} // }}}

/**
 * \brief Tensor power of given state (complex vector).
 *
 * Calculates the tensor product as many times as specified with the second argument.
 *
 * state ⊗ state ⊗ state ⊗ ... ⊗ state (`qub` times)
 *
 * \param state [itpp::cvec] state in the form of complex vector.
 * \param qub [int] power, how many times the tensor product will be applied to the state.
 *
 * \return resulting vector [itpp::cvec], with complex elements.
*/
itpp::cvec TensorPow(const itpp::cvec state, int qub){ //{{{
  itpp::cvec newstate;
  newstate=state;
  for(int i=0;i<qub-1;i++){
    newstate=TensorProduct(newstate,state);
  }
  return newstate;
} //}}}

/**
 * \brief Tensor power of given complex matrix.
 *
 * Calculates the tensor product as many times as specified with the second argument.
 *
 * Matrix ⊗ Matrix ⊗ Matrix ⊗ ... ⊗ Matrix (`qub` times)
 *
 * \param state [itpp::cvec] matrix to calculate the tensor product.
 * \param qub [int] power, how many times the tensor product will be applied to the matrix.
 *
 * \return resulting matrix [itpp::cmat], with complex elements.
*/
itpp::cmat TensorPow(const itpp::cmat Matrix, int qub){ //{{{
  itpp::cmat newMatrix;
  newMatrix=Matrix;
  for(int i=0;i<qub-1;i++){
    newMatrix=TensorProduct(newMatrix,Matrix);
  }
  return newMatrix;
} //}}}

// Exponentiation
itpp::mat exponentiate_real_symmetric(const itpp::mat& Matrix){ // {{{
// template <class Num_T> itpp::Mat<Num_T> exponentiate(itpp::Mat<Num_T> Matrix){
  //Calculates exp(H)
  itpp::mat tmp;
  itpp::vec eigenvalues;
  itpp::mat V;
  bool exito = itpp::eig_sym(Matrix, eigenvalues, V);
  double how_symmetric = test_symmetry(Matrix);
  if (exito==false || how_symmetric > 10e-12){
    std::cerr << "Algun pedo en exponentiate_real_symmetric\n";
    abort();
  }
  tmp = V*itpp::diag(itpp::exp(eigenvalues))* V.hermitian_transpose();
  return tmp;
} // }}}
itpp::cmat exponentiate(const itpp::cmat& Matrix){ // {{{
//   std::cout << "@exponentiate, 1" << std::endl;
// template <class Num_T> itpp::Mat<Num_T> exponentiate(itpp::Mat<Num_T> Matrix){
  //Calculates exp(H)
  double how_real_symmetric = test_real_symmetric(Matrix);
  if (how_real_symmetric < 10e-14){
    return itpp::to_cmat(exponentiate_real_symmetric(itpp::real(Matrix)));
  }

  itpp::cmat tmp;
  itpp::cvec eigenvalues;
  itpp::cmat V;
  bool exito = itpp::eig(Matrix, eigenvalues, V);
  double revision;
  revision = itpp::norm(Matrix -
      V*itpp::diag(eigenvalues)* V.hermitian_transpose())/itpp::norm(Matrix);
//   std::cout << "@exponentiate, 50. Error es " << revision << std::endl;
  if (revision > 10e-12 || !exito){
    std::cerr << "Matrix=" << Matrix << std::endl;
    std::cerr << "Eigenvalores=" << eigenvalues << std::endl;
    std::cerr << "Ortonormalidad de eigenvectores=" << itpp::norm(V*V.hermitian_transpose()-itpp::eye_c(Matrix.size())) << std::endl;
    std::cerr << "Ortonormalidad de eigenvectores, matriz completa=\n" << abs(V*V.hermitian_transpose()) << std::endl;
    for (int i=0; i<Matrix.rows(); i++){
      std::cerr << "Ecuacion de eigengectores sirve? i=" << i << ", " << norm(Matrix*V.get_col(i)-eigenvalues(i)*V.get_col(i)) << std::endl;
    }
    std::cerr << "Algun pedo en exponentiate! \nRevision=" << revision
      <<"Normas, "<<itpp::norm(Matrix) << ", " << exito << std::endl;
    std::cerr << "La rutina de lapack tiene pedos cuando esta muy sparse la matriz,\n"
      <<"por ejemplo un elemento del modelo de free fermions con 2D en primeros vecinos" << std::endl;
    abort();
  }
  tmp = V*itpp::diag(itpp::exp(eigenvalues))* V.hermitian_transpose();
//   std::cout << "@exponentiate, 100" << std::endl;
  return tmp;
} // }}}
itpp::cmat exponentiate_nonsym(const itpp::cmat& Matrix){ // {{{
// template <class Num_T> itpp::Mat<Num_T> exponentiate(itpp::Mat<Num_T> Matrix){
  //Calculates exp(H) for NON-HERMITIAN matrix
  double how_real_symmetric = test_real_symmetric(Matrix);
  if (how_real_symmetric < 10e-14){
    return itpp::to_cmat(exponentiate_real_symmetric(itpp::real(Matrix)));
  }

  itpp::cmat tmp;
  itpp::cvec eigenvalues;
  itpp::cmat V;
  itpp::cmat invV;
  bool exito = itpp::eig(Matrix, eigenvalues, V);
  bool exito2 = itpp::inv(V,invV);
  double revition;
  revition = itpp::norm(Matrix -
      (V*itpp::diag(eigenvalues)* invV));
  if (revition > 10e-12 || !exito || !exito2){
    std::cerr << "Algun xxx pedo en exponentiate " << revition << ", " << exito << std::endl;
    std::cerr << "La rutina de lapack tiene pedos cuando esta muy sparse la matriz,\n"
      <<"por ejemplo un elemento del modelo de free fermions con 2D en primeros vecinos" << std::endl;
    abort();
  }
  tmp = V*itpp::diag(itpp::exp(eigenvalues))* invV;
  return tmp;
} // }}}
itpp::cmat exponentiate(const itpp::mat& Matrix){ // {{{
// template <class Num_T> itpp::Mat<Num_T> exponentiate(itpp::Mat<Num_T> Matrix){
  //Calculates exp(H)
  double how_symmetric = test_symmetry(Matrix);
  if (how_symmetric < 10e-14){
    return itpp::to_cmat(exponentiate_real_symmetric(Matrix));
  }

  itpp::cmat tmp;
  itpp::cvec eigenvalues;
  itpp::cmat V;
  bool exito = itpp::eig(Matrix, eigenvalues, V);
  double revition;
  revition = itpp::norm(Matrix -
      V*itpp::diag(eigenvalues)* V.hermitian_transpose());
  if (revition > 10e-12 || !exito){
    std::cerr << "Algun pedo en exponentiate" << revition << ", " << exito << std::endl;
    std::cerr << "La rutina de lapack tiene pedos cuando esta muy sparse la matriz,\n"
      <<"por ejemplo un elemento del modelo de free fermions con 2D en primeros vecinos" << std::endl;
    abort();
  }
  tmp = V*itpp::diag(itpp::exp(eigenvalues))* V.hermitian_transpose();
  return tmp;
} // }}}
// Others
void applyUDiagUdagger(const itpp::cmat& U, const itpp::cvec& E,  itpp::cvec& psi){ // {{{
  // Calculates U.diag(E).U^\dagger.psi with diag(E) a diagonal matrix
  psi = U.hermitian_transpose()*psi;
  psi = elem_mult(E, psi);
  psi = U*psi;
  return ;
} // }}}
itpp::cmat UDiagUdagger(const itpp::cmat& U, const itpp::cvec& E){ // {{{
  // Calculates U.diag(E).U^\dagger with diag(E) a diagonal matrix
// template <class Num_T> itpp::Mat<Num_T> exponentiate(itpp::Mat<Num_T> Matrix){
  //Calculates exp(H)
  int n=E.length();
  itpp::cmat Result(n,n);
  std::complex<double> A;

  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      A=0.;
      for (int k=0; k<n; k++){
        A+=U(i,k)*E(k)*std::conj(U(j,k));
      }
      Result(i,j)=A;
    }
  }
  return Result;
} // }}}
itpp::cmat UDiagUdagger(const itpp::cmat& U, const itpp::vec& E){ // {{{
return UDiagUdagger(U, to_cvec(E));
} // }}}
itpp::cmat InteractionPicture(const itpp::cmat& H0, const itpp::cmat& As, double t){ // {{{
  std::cout << "Warning, using a slow routine asdfasd" << std::endl;
  // Returns exp[i t H0] As exp[-i t H0]
  // acording to wikipedia, this is the standard sign
  // http://en.wikipedia.org/wiki/Interaction_picture
  itpp::cmat U;
  std::complex<double> I(0.,1.);
  U = exponentiate(-I*t*H0);
//   std::cout << "She was all dressed up in black" << std::endl;
  itpp::cmat tmp;
  return itpp::hermitian_transpose(U)*As*U;
//   return itpp::hermitian_transpose(U)*As*U;
} // }}}
itpp::cmat InteractionPicture(const itpp::cmat& V0, const itpp::cmat& As){ // {{{
  std::cout << "Warning, using a slow routinen 42384938274" << std::endl;
  // Returns V0^\dagger As V0
  // acording to wikipedia, this is the standard sign
  // http://en.wikipedia.org/wiki/Interaction_picture
  return itpp::hermitian_transpose(V0)*As*V0;
//   return itpp::hermitian_transpose(U)*As*U;
} // }}}
template <class Num_T> itpp::Mat<Num_T> Proyector(const itpp::Vec<Num_T>& psi){ // {{{
  int n=psi.size();
  itpp::Mat<Num_T> p(n,n);
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      p(i,j) = psi(i)*conj(psi(j));
    }
  }
  return p;
} // }}}
template <class Num_T> itpp::Mat<Num_T> DirectSum(const itpp::Array<itpp::Mat<Num_T> >& Matrices){ // {{{

  itpp::Mat<Num_T> tmp=Matrices(0);
  for (int i=1; i<Matrices.size(); i++){
    tmp=DirectSum(tmp, Matrices(i));
  }
  return tmp;
} // }}}
template <class Num_T> itpp::Mat<Num_T> DirectSum(const itpp::Mat<Num_T>& Matrix1, const itpp::Mat<Num_T>& Matrix2){ // {{{
  itpp::Mat<Num_T> tmp(Matrix1.rows()+ Matrix2.rows(), Matrix1.cols()+ Matrix2.cols());
  tmp.zeros();
  tmp.set_submatrix(0,0, Matrix1);
  tmp.set_submatrix(Matrix1.rows(), Matrix1.cols(), Matrix2);
  return tmp;
} // }}}
template <class Num_T> itpp::Vec<Num_T> get_column(const itpp::Array<itpp::Vec<Num_T> >& A, const int element){ // {{{
  itpp::Vec<Num_T > tmp(A.size());
  for (int i = 0; i< A.size(); i++){ tmp(i) = A(i)(element); }
  return tmp;
} //}}}
// Multiplication
std::complex<double> expected_value(const itpp::cmat& A, const itpp::cvec& psi){ // {{{
  return itpp::dot(conj(psi), A*psi);
} // }}}
template <class Num_T> itpp::Mat<Num_T> AtimesDiagB(const itpp::Mat<Num_T>& A, const itpp::Vec<Num_T>& B){ // {{{
  itpp::Mat<Num_T> tmp(A.cols(),A.rows());
  for (int i=0; i<tmp.cols(); i++){
    tmp.set_col(i,B(i)*A.get_col(i));
  }
  return tmp;
} // }}}
// template <class Num_T> itpp::Vec<Num_T> DiagonalTimesVector(const itpp::Vec<Num_T>& DiagMatrix, const itpp::Vec<Num_T>& Psi){ // {{{
// //   itpp::Mat<Num_T> tmp(A.cols(),A.rows());
//   return;Psi = DiagMatrix* Psi;
//   return;
// } // }}}
// }}}
// Reporting {{{
void PrintCompactHermitian(itpp::Mat<std::complex<double> >& rho){ // {{{
	int dim=rho.cols();
	std::cout <<" ";
	for (int i =0; i< dim;i++){
		std::cout <<real(rho(i,i))<<" ";
// 		std::cout << "x"<<rho(i,i)<<" ";
		for (int j =i+1; j< dim;j++){
			std::cout <<real(rho(i,j))<<" "<<imag(rho(i,j)) << " ";
		}
	}
} // }}}
// }}}
// Bit manipulation {{{
itpp::Array< itpp::vec> numbers_sorted_by_bits_on(int number_of_bits){ // {{{
  int q=number_of_bits;
  itpp::Array< itpp::vec> numeros_por_bits_on(q+1);
  itpp::ivec entradas_metidas(q+1);
  entradas_metidas=0;
  int size_i, bits_on;
  int dim=cfpmath::pow_2(q);
  // Calculo de numeros con k qubits
  for (int i=0; i<=q; i++){
    size_i = cfpmath::factorial(q)/(cfpmath::factorial(i)*cfpmath::factorial(q-i));
    numeros_por_bits_on(i).set_size(size_i);
//         cout << i<<", "<< size_i << ", " << dim << endl;
  }
  for (int i=0; i<dim; i++){
    bits_on = cfpmath::BitCount(i);
//         cout << i<<", "<< bits_on << endl;
    numeros_por_bits_on(bits_on)(entradas_metidas(bits_on)) = i;
    entradas_metidas(bits_on)++;
  }
  return numeros_por_bits_on;
}// }}}
itpp::ivec IntegerDigits(int n,int base=2, int length=0){// {{{
  int number=n;
  int size=length;
  if (length==0){
    size = cfpmath::integer_part_log(base,n);
  }
  itpp::ivec result(size);
  result=0;

  int counter=0;
  while (n!=0 && counter<size){
    result(counter) = number%base;
    number=(number - result(counter))/base;
    counter ++;
  }
  return result;
}// }}} Matrix=
itpp::ivec IntegerDigits(const int n, itpp::ivec base){// {{{
  int number=n;
  itpp::ivec result(base.size());
  result=0;

  for (int counter=0; counter<base.size(); counter++){
    result(counter) = number%(base(counter));
    number=(number - result(counter))/base(counter);
  }
  return result;
}// }}} Matrix=
int bvec_to_int(itpp::Vec<bool> bool_array){// {{{
  int tmp=0;
  for (int i=0; i<bool_array.size(); i++){
    if (bool_array(i)){
      tmp= tmp | (1<<i);
    }
  }
  return tmp;
}// }}}
itpp::Vec<bool> int_to_bool(const int n, const int total_number_bits){// {{{
  itpp::Vec<bool> tmp(total_number_bits);
  for (int i=0; i<total_number_bits; i++){ tmp(i)=cfpmath::test_bit(n,i); }
  return tmp;
}// }}}
itpp::ivec all_n_bit_numbers(int bits_on, int total_number_bits){// {{{
  itpp::ivec tmp(0);
  for (int i=0; i<cfpmath::pow_2(total_number_bits); i++){
//     std::cout << i <<" "<<  cfpmath::BitCount(i)  << "\n";
    if (cfpmath::BitCount(i)==bits_on){
      tmp = itpp::concat(tmp,i);
//       std::cout << i <<" "<<  cfpmath::BitCount(i)  << " siiii\n";
    }
  }
  return tmp;
}// }}}
itpp::ivec all_bit_rotations(int n, int size_register){ //{{{
  itpp::ivec tmp(size_register);
  tmp(0)=n;
  int candidate;
  for( int i =0 ; i<size_register-1 ; i++){
    candidate=cfpmath::rotate_bits(tmp(i), size_register);
    if(candidate == n){
      tmp.set_size(i+1, true);
      break;
    }
//     std::cout << "i="<< i << std::endl;
    tmp(i+1)=candidate;
  }
//   return ;
  return tmp;
} // }}}
// }}}
// Functional Element manipulation, on all elements {{{
template <class Num_T> Num_T  Chop(const Num_T& A, double epsilon=1e-12){ // {{{
  if (abs(A)<epsilon){
    return 0;
  } else {
    return A;
  }
} // }}}
template <class Num_T> itpp::Vec<Num_T > Chop(const itpp::Vec<Num_T>& A, double epsilon=1e-12){ // {{{
  itpp::Vec<Num_T > tmp(A.size());
  for (int i=0; i<A.size(); i++){
    tmp(i)=Chop(A(i));
  }
  return tmp;
} // }}}
template <class Num_T> itpp::Mat<Num_T > Chop(const itpp::Mat<Num_T>& A, double epsilon=1e-12){ // {{{
  itpp::Mat<Num_T > tmp(A.rows(),A.cols());
  for (int r=0; r<A.rows(); r++){
    for (int c=0; c<A.cols(); c++){
      tmp(r,c)=Chop(A(r,c),epsilon);
    }
  }
  return tmp;
} // }}}
template <class Num_T, class Cosa> itpp::Array<Num_T > ElementwiseAddInFirstElement(const Cosa joda, const itpp::Array<Num_T>& A){ // {{{
  itpp::Array<Num_T> tmp(A.size());
  for (int i=0; i<A.size(); i++){
    tmp(i)=itpp::concat(joda, A(i));
  }
  return tmp;
} // }}}
template <class Num_T> itpp::Array<Num_T > ElementwiseRemoveFirstElement(const itpp::Array<Num_T>& A){ // {{{
  itpp::Array<Num_T > tmp(A.size());
  for (int i=0; i<A.size(); i++){
    tmp(i)=A(i)(1,A(i).size()-1);
  }
  return tmp;
} // }}}
template <class Num_T> itpp::Array<Num_T > Elementwiseconcat(const Num_T& A, const itpp::Array<Num_T >& B){ // {{{
  itpp::Array<Num_T > tmp(B.size());
  for (int i=0; i<B.size(); i++){
    tmp(i)=itpp::concat(A,B(i));
  }
  return tmp;
} // }}}
itpp::ivec Modulize(const itpp::ivec& intvec,const itpp::ivec n){ // {{{
  itpp::Vec<int> tmp(intvec.size());
  for (int i =0;i<intvec.size();i++){
    tmp(i)=intvec(i)%n(i);
    if (tmp(i)<0){
      tmp(i)=tmp(i)+n(i);
    }
  }
  return tmp;
} // }}}
itpp::ivec Modulize(const itpp::ivec& intvec,const int n){ // {{{
  itpp::Vec<int> tmp(intvec.size());
  for (int i =0;i<intvec.size();i++){
    tmp(i)=intvec(i)%n;
    if (tmp(i)<0){
      tmp(i)=tmp(i)+n;
    }
  }
  return tmp;
} // }}}
itpp::Array<itpp::ivec> Modulize(const itpp::Array<itpp::ivec>& intvec,const itpp::ivec n){ // {{{
  itpp::Array<itpp::ivec> tmp(intvec.size());
  for (int i =0;i<intvec.size();i++){
    tmp(i)=Modulize(intvec(i), n);
  }
  return tmp;
} // }}}
itpp::Array<itpp::ivec> Modulize(const itpp::Array<itpp::ivec>& intvec,const int n){ // {{{
  itpp::Array<itpp::ivec> tmp(intvec.size());
  for (int i =0;i<intvec.size();i++){
    tmp(i)=Modulize(intvec(i), n);
  }
  return tmp;
} // }}}
itpp::Mat<int> Modulize(const itpp::Mat<int>& set_of_vectors,const int n){ // {{{
  itpp::Mat<int> tmp(set_of_vectors.rows() ,set_of_vectors.cols() );
  for (int i=0; i < set_of_vectors.rows(); i++){
    tmp.set_row(i,Modulize(set_of_vectors.get_row(i),n) );
  }
  return tmp;
} // }}}
itpp::Array<itpp::ivec> DeScale(const itpp::Array<itpp::ivec>& set_of_vectors,const int scale_factor){ // {{{
  itpp::Array<itpp::ivec> tmp(set_of_vectors.size());
  for (int i=0; i < set_of_vectors.size(); i++){
    tmp(i)=set_of_vectors(i)/scale_factor;
  }
  return tmp;
} // }}}
itpp::Array<itpp::ivec> Scale(const itpp::Array<itpp::ivec>& set_of_vectors,const int scale_factor){ // {{{
  itpp::Array<itpp::ivec> tmp(set_of_vectors.size());
  for (int i=0; i < set_of_vectors.size(); i++){
    tmp(i)=scale_factor*set_of_vectors(i);
  }
  return tmp;
} // }}}
itpp::Array<itpp::ivec> Displace(const itpp::Array<itpp::ivec>& set_of_vectors,const itpp::ivec displacement_vector){ // {{{
  itpp::Array<itpp::ivec> tmp(set_of_vectors.size());
  for (int i=0; i < set_of_vectors.size(); i++){
    tmp(i)=set_of_vectors(i)+displacement_vector;
  }
  return tmp;
} // }}}
itpp::Mat<int> Displace(const itpp::Mat<int>& set_of_vectors,const itpp::Vec<int> displacement_vector){ // {{{
  itpp::Mat<int> tmp(set_of_vectors.rows() ,set_of_vectors.cols() );
  for (int i=0; i < set_of_vectors.rows(); i++){
    tmp.set_row(i,set_of_vectors.get_row(i)+displacement_vector);
  }
  return tmp;
} // }}}
itpp::Array<itpp::ivec> DisplaceScaleModulize( const itpp::Array<itpp::ivec>& set_of_vectors,const itpp::ivec displacement_vector, const int scale_factor, const int n){ // {{{
  itpp::Array<itpp::ivec> tmp(set_of_vectors.size());
  for (int i=0; i < set_of_vectors.size(); i++){
    tmp(i)=set_of_vectors(i);
    tmp(i)=scale_factor*tmp(i);
    tmp(i)=tmp(i)+displacement_vector;
    tmp(i)=Modulize(tmp(i),n);
  }
  return tmp;

} // }}}
template <class Num_T> itpp::Array<Num_T*> to_pointers_1(itpp::Array<Num_T>& v){ // {{{
  itpp::Array<Num_T*> result(v.size());
  for (int i=0; i<v.size(); i++){
    result(i)=&(v(i));
  }
  return result;
} // }}}
template <class Num_T> itpp::Array<itpp::Array<Num_T*> > to_pointers_2(itpp::Array<itpp::Array<Num_T> >& v){ // {{{
  itpp::Array<itpp::Array<Num_T*> > result(v.size());
  for (int i=0; i<v.size(); i++){
    result(i).set_size(v(i).size());
    for (int j=0; j<v(i).size(); j++){
      result(i)(j)=&(v(i)(j));
    }
  }
  return result;
} // }}}
// }}}
// 2D grid specialized functions {{{
itpp::ivec GetBlock(const itpp::Array<itpp::ivec>& positions, const int length_cell, const int length_grid){
  if (positions.size() == 1 && positions(0)==itpp::vec_2(0,0)){
    // std::cerr << "Poaca papa, positions.size()="<< positions.size() << "\n";
    return itpp::vec_2(0,0);
  }
  if (positions.size() != 4){
    std::cerr << "Something wrong in GetBlock sdlkjhflakshf, positions.size()="<< positions.size()
      << ", length_cell=" << length_cell<<"\n";
    abort();
  }
  itpp::ivec position_block(positions(0).size() ), las;
  int cell_size;
  for (int i=0; i< positions(0).size(); i++){
    las =Union(get_column(positions,i));
    itpp::sort(las);
//     las =itppextmath::sort(Union(get_column(positions,i)));
    if (las.size() != 2){
      std::cerr << "Something wrong in GetBlock o8e7dkjhf\n";
      abort();
    }
    if ( las(0) == 0){
      if( cfpmath::test_int_log(length_cell, las(1))){
        position_block(i)=0;
        cell_size= las(1)-las(0);
      } else if ( cfpmath::test_int_log(length_cell, length_grid - las(1)) ) {
        cell_size = length_grid - las(1);
        position_block(i) = length_grid / cell_size -1;
      }
    } else {
      cell_size = las(1) - las(0);
      position_block(i)=las(0)/cell_size;
    }
  }
  return position_block;
}
void get_base_position_and_delta(const  itpp::Array<itpp::ivec > Positions, int& delta, itpp::ivec& base_position){
  int candidate_delta;
  base_position=Positions(0);
  for (int i=1; i< Positions.size(); i++){
    if (Positions(i)(0) < base_position(0)) base_position(0) = Positions(i)(0);
    if (Positions(i)(1) < base_position(1)) base_position(1) = Positions(i)(1);
  }
  delta = -1;
  for (int i=0; i< Positions.size() ; i++){
    candidate_delta= itpp::max(Positions(i)-base_position );
    if (candidate_delta > 0 && delta < 0) delta = candidate_delta;
    if (candidate_delta < delta)          delta = candidate_delta;
  }
  return;
}
itpp::Array<itpp::ivec> rotationaly_order_vector_renormalized_lattices_3(const itpp::Array<itpp::ivec >& A, const int length_mera){
  // This routine identifies the canonical order of a set of points
  // that me a valid set of points after a renormalization step in a
  // 3^n x 3^n or 2*3^n x 2*3^n lattice
  if (A.size() != 4){
    std::cerr << "@rotationaly_order_vector_renormalized_lattices_3. Error, the size must be 4. Actual size: "
      << A.size() << std::endl;
    abort();
  }

  int Delta;
  itpp::ivec x = Union(get_column(A,0));sort(x);
  itpp::ivec y = Union(get_column(A,1));sort(y);
  itpp::ivec base(2);base=itpp::vec_2(x(0), y(0));
  Delta=x(1)-x(0);

  if (x(0)==0 && !cfpmath::test_int_log(3, Delta)){
      Delta = cfpmath::maximum_prime_power_divisor(x(1),3);
      base(0)=x(1);
  }
  if (y(0)==0 && !cfpmath::test_int_log(3, y(1)-y(0))){
      base(1)=y(1);
  }
  itpp::Array<itpp::ivec> DesiredPositions(4);
  DesiredPositions(0) = base;
  DesiredPositions(1) = base+itpp::vec_2(0,Delta);
  DesiredPositions(2) = base+itpp::vec_2(Delta,Delta);
  DesiredPositions(3) = base+itpp::vec_2(Delta,0) ;
  DesiredPositions = Modulize(DesiredPositions, length_mera);
  if (!AreEqual_modulo_order(DesiredPositions, A)){
    std::cerr << "@rotationaly_order_vector_renormalized_lattices_3. Error, las posiciones no coinciden:"
      << "DesiredPositions: " << DesiredPositions << std::endl
      << "A: " << A << std::endl;
  }
  return DesiredPositions;
}
itpp::Array<itpp::ivec> rotational_order(const itpp::Array<itpp::ivec >& A){
  itpp::ivec pos_x = get_column(A,0), pos_y = get_column(A,1);
  itpp::ivec base = itpp::vec_2(itpp::min(pos_x), itpp::min(pos_y));
  int Delta =  itpp::max(pos_x) - itpp::min(pos_x);
  itpp::Array<itpp::ivec> DesiredPositions(4);
  DesiredPositions(0) = base;
  DesiredPositions(1) = base+itpp::vec_2(0,Delta);
  DesiredPositions(2) = base+itpp::vec_2(Delta,Delta);
  DesiredPositions(3) = base+itpp::vec_2(Delta,0) ;
  return DesiredPositions;
}
// }}}
// Fermionic specialized functions {{{
void UpDownMostIndexFermionicMatrix(itpp::cmat& MatrixEven, itpp::cmat& MatrixOdd, const std::string direction, const std::string position){
  if ( (position != "most" && position != "least") ||
      (direction != "outgoing_to_incomming" && direction != "incomming_to_outgoing")){
    std::cerr << "Error in DownMostSignificantIndexFermionicMatrix, xuweor98\n"
      << "direction:"<< direction << "\n"
      << "position:"<< position << "\n"
      << "( ):"<< (direction != "outgoing_to_incomming" && direction != "incomming_to_outgoing") << "\n"
      << "(position != most && position != least):"<< (position != "most" && position != "least") << "\n"
      ;
    abort();
  }
  int d_in  = MatrixEven.cols(), d_out = MatrixEven.rows();
  if (d_in != MatrixOdd.cols() || d_out != MatrixOdd.rows()){
    std::cerr << "Error in DownMostSignificantIndexFermionicMatrix, o2u234\n";
    abort();
  }
  itpp::cmat new_even, new_odd;
  if (direction == "outgoing_to_incomming"){
    new_even.set_size(d_out/2,2*d_in);
    new_odd.set_size(d_out/2,2*d_in);
    if (position=="most"){
      new_even.set_submatrix(0,0, MatrixEven(0, d_out/2 -1, 0, d_in-1));
      new_even.set_submatrix(0,d_in, MatrixOdd(d_out/2, d_out -1, 0, d_in-1));
      new_odd.set_submatrix(0,0, MatrixOdd(0, d_out/2 -1, 0, d_in-1));
      new_odd.set_submatrix(0,d_in, MatrixEven(d_out/2, d_out -1, 0, d_in-1));
    }
    if (position=="least"){
      bool parity_i, parity_j;
      std::complex<double> matrix_element;
      for (int i=0; i < d_out; i++){
        for (int j=0; j < 2*d_in; j++){
          parity_i = cfpmath::parity_sum_digits_base_2(i);
          parity_j = cfpmath::parity_sum_digits_base_2(j);
          if (parity_j == 0 ){
            matrix_element = MatrixEven(i,j/2);
          } else {
            matrix_element = MatrixOdd(i,j/2);
          }
          if (parity_i == 0 ){
            new_even(i/2,j) = matrix_element;
          } else {
            new_odd(i/2,j) = matrix_element;
          }
        }
      }
    }
  }
  if (direction == "incomming_to_outgoing"){
    new_even.set_size(d_out*2,d_in/2);
    new_odd.set_size(d_out*2,d_in/2);
    if (position=="most"){
      new_even.set_submatrix(0,0, MatrixEven(0, d_out -1, 0, d_in/2-1));
      new_even.set_submatrix(d_out, 0, MatrixOdd(0, d_out -1, d_in/2, d_in-1));
      new_odd.set_submatrix(0,0, MatrixOdd(0, d_out -1, 0, d_in/2 -1));
      new_odd.set_submatrix(d_out, 0, MatrixEven(0, d_out -1, d_in/2, d_in-1));
    }
    if (position=="least"){
      bool parity_i, parity_j;
      std::complex<double> matrix_element;
      for (int i=0; i < 2*d_out; i++){
        for (int j=0; j < d_in; j++){
          parity_i = cfpmath::parity_sum_digits_base_2(i);
          parity_j = cfpmath::parity_sum_digits_base_2(j);
          if (parity_i == 0 ){
            matrix_element = MatrixEven(i/2,j);
          } else {
            matrix_element = MatrixOdd(i/2,j);
          }
          if (parity_j == 0 ){
            new_even(i,j/2) = matrix_element;
          } else {
            new_odd(i,j/2) = matrix_element;
          }
        }
      }
    }
  }
  MatrixEven = new_even;
  MatrixOdd = new_odd;
  return;
}
void UpDownMostIndexFermionicMatrix(itpp::cmat& MatrixEven, itpp::cmat& MatrixOdd, const int places_most, const int places_least){
//! Routine to move the matrices so as to simulate moving an index up or down.
/*! This routine works at the level of matrices as for qubits I have not a good
 * idea of what would mean when two fermions with the same name are comming
 * in (or out). The positive sign indicates clockwise movement.
 */
  if (places_most>0){
    for (int i=0; i< places_most; i++){
      UpDownMostIndexFermionicMatrix(MatrixEven, MatrixOdd, "incomming_to_outgoing", "most");
    }
  }
  if (places_most<0){
    for (int i=0; i> places_most; i--){
      UpDownMostIndexFermionicMatrix(MatrixEven, MatrixOdd, "outgoing_to_incomming", "most");
    }
  }
  if (places_least>0){
    for (int i=0; i< places_least; i++){
      UpDownMostIndexFermionicMatrix(MatrixEven, MatrixOdd, "outgoing_to_incomming", "least");
    }
  }
  if (places_least<0){
    for (int i=0; i> places_least; i--){
      UpDownMostIndexFermionicMatrix(MatrixEven, MatrixOdd, "incomming_to_outgoing", "least");
//       std::cerr << "hello " << i << "\n";
    }
  }
  return;
}
// }}}
// Quantum Info, and advanced functions {{{
// Some important matrices of operations and states
itpp::cmat sigma(itpp::ivec is){ // {{{
  int qubits=is.size();
  itpp::Array<itpp::cmat> sigmas(qubits);
  if (qubits==0){

    return itpp::mat_1x1( std::complex<double> (1,0));
  }
  for (int i=0; i<qubits; i++){
    sigmas(qubits-i-1)=sigma(is(i));
  }
  return TensorProduct(sigmas);
} // }}}
itpp::cmat sigma(int i){ // {{{
  std::complex<double> imaginary(0,1), c0(0.,0.);
  switch (i) {
    case 0:
      return itpp::to_cmat(itpp::mat_2x2(1, 0, 0, 1));
    case 1:
      return itpp::to_cmat(itpp::mat_2x2(0, 1, 1, 0));
    case 2:
      return itpp::mat_2x2(c0, -imaginary, imaginary, c0);
    case 3:
      return itpp::to_cmat(itpp::mat_2x2(1, 0, 0, -1));
    default:
      std::cout << "Sigma not found. i=" << i << "\n";
      abort();
  }
} // }}}
itpp::cmat sigma(itpp::vec b){ // {{{
  if (b.size()==3){
    return b(0)*sigma(1)+b(1)*sigma(2)+ b(2)*sigma(3);
  } else {
    std::cout << "Sigma not found 89p790as7dfu.\n";
    abort();
  }
} // }}}
itpp::cmat sigma(itpp::vec b, int position, int total){ // {{{
  itpp::cmat tmp = sigma(b);
  if (position==total-1){
    tmp=TensorProduct(tmp,itpp::eye_c(cfpmath::pow_2(total-1)));
  } else if (position==0){
    tmp=TensorProduct(itpp::eye_c(cfpmath::pow_2(total-1)),tmp);
  } else if (0<position && position<total-1){
    tmp=TensorProduct(tmp,itpp::eye_c(cfpmath::pow_2(position)));
    tmp=TensorProduct(itpp::eye_c(cfpmath::pow_2(total-position-1)),tmp);
  } else {
    std::cout << "Sigma not found asdfa \n";
    abort();
  }
  return tmp;
} // }}}
itpp::mat hadamard_matrix(){// {{{
  itpp::mat tmp(2,2);
  tmp=1;
  tmp(1,1)=-1;
  tmp=tmp/sqrt(2.);
  return tmp;
} //}}}
itpp::mat swap_matrix(){// {{{
  itpp::mat tmp(4,4);
  tmp=0.;
  tmp(0,0)=1;
  tmp(1,2)=1;
  tmp(2,1)=1;
  tmp(3,3)=1;
  return tmp;
} //}}}
itpp::mat cnot_matrix(){// {{{
  itpp::mat tmp(4,4);
  tmp=0.;
  tmp(0,0)=1;
  tmp(1,1)=1;
  tmp(2,3)=1;
  tmp(3,2)=1;
  return tmp;
} //}}}
itpp::cmat control_u_matrix(itpp::cmat& u){// {{{
  itpp::cmat tmp(4,4);
  tmp=0.;
  tmp(0,0)=1;
  tmp(1,1)=1;
  tmp.set_submatrix(2,2,u);
  return tmp;
} //}}}
itpp::ivec diagonal_sigma_z(int encoded_positions, int qubits){ // {{{
  itpp::ivec tmp(cfpmath::pow_2(qubits));
  tmp=1;
//   std::cout << "Encoded positions = " << encoded_positions << std::endl;
  for (int i=0; i<cfpmath::pow_2(qubits); i++){
//     std::cout << "i="<<i<<", i&encoded_positions = " << (i&encoded_positions) << std::endl;
    if (cfpmath::parity_sum_digits_base_2(i&encoded_positions)){
      tmp(i)=-1;
    }
  }
  return tmp;
} // }}}

/**
 * \brief Generates a Werner state.
 *
 * Generates the density matrix of a werner state with a given alpha value.
 *
 * \param alpha [const double] value between 0 and 1.
 *
 * \return generated Werner state within a matrix [itpp::Mat].
*/
itpp::Mat<std::complex<double> > Werner(const double alpha){// {{{
	itpp::Mat<std::complex<double>  > tmp(4,4);
 	tmp=0.;
	tmp(0,0)=tmp(3,3)=(2-alpha)/4;
	tmp(1,1)=tmp(2,2)=alpha/4;
	tmp(0,3)=tmp(3,0)=(1-alpha)/2;
	return tmp;
}// }}}

itpp::cvec StatisticallyNormalRandomState(unsigned int dim){ // {{{
  itpp::cvec tmp=itpp::randn_c(dim);
  return tmp/sqrt(dim);
} // }}}
/**
 * \brief Generates a random state.
 *
 * Given a specified dimension, it generates a normalized random state.
 *
 * \param dim [unsigned int] random state dimension.
 *
 * \return randomly generated state within a vector [itpp::cvec].
*/
itpp::cvec RandomState(unsigned int dim){ // {{{
  itpp::cvec tmp=itpp::randn_c(dim);
  return tmp/norm(tmp);
} // }}}
itpp::cvec RandomSeparableState(unsigned int q){ // {{{
//   std::cout << "entrando a RandomSeparableState, q= " << q <<std::endl;
  int dim=cfpmath::pow_2(q);
  itpp::cmat z=itpp::randn_c(q,2);
//   z(0,0)=0;
//   z(0,1)=1;
//   z(1,1)=0;
//   z(2,1)=0;
//   z(3,1)=0;
//   z(1,0)=0;
//   z(2,0)=0;
//   z(3,0)=1;
  itpp::cvec psi(dim);
  std::complex<double> c;
//   std::cout << "z=" << z << std::endl;

  for (int i=0; i< dim; i++){
    c=1;
//     std::cout << "i=" << i << std::endl;
    for (unsigned int b=0; b<q; b++){
//        std::cout << "test result:" << cfpmath::test_bit(i,b) << std::endl;
      if (cfpmath::test_bit(i,b)){
        c*=z(b,0);
//         std::cout << "In 0 " << b << std::endl;
      } else {
        c*=z(b,1);
//         std::cout << "In 1 " << b << std::endl;
      }
    }
    psi(i)=c;
  }
//   psi= psi/norm(psi);
//   std::cout << "norma = " << norm(psi) << std::endl;
  return psi/norm(psi);
} // }}}
/**
 * \brief Generates a basis state.
 *
 * Given a dimension and a position for the basis, it generates a basis state.
 *
 * \param dim [unsigned int] basis state dimension.
 * \param basis_number [unsigned int] basis state dimension.
 *
 * \return generated basis state within a vector [itpp::cvec].
*/

itpp::cvec BasisState(unsigned int dim, unsigned int basis_number){ // {{{
  itpp::cvec tmp(dim);
  tmp.zeros();
  tmp(basis_number)=1;
  return tmp;
} // }}}

/**
 * \brief Generates a bell state.
 *
 * Generates a bell state with specified dimension, default is 4.
 *
 * \param dim [int] bell state dimension.
 *
 * \return generated bell state within a vector [itpp::vec].
*/
itpp::vec BellState(int dim=4){// {{{
  itpp::vec tmp(dim);
  tmp.zeros();
  tmp(0)=1;
  tmp(dim-1)=1;
  return tmp/sqrt(2.);
}// }}}
itpp::cvec BlochToQubit(double theta, double phi){// {{{
  itpp::cvec tmp(2);
  std::complex<double> I(0,1);
  tmp.zeros();
  tmp(0)=cos(theta/2);
  tmp(1)=sin(theta/2)*exp(I*phi);
  return tmp;
}// }}}
itpp::vec BellState(double theta){// {{{
  int dim=4;
  itpp::vec tmp(dim);
  tmp.zeros();
  tmp(0)=cos(theta);
  tmp(dim-1)=sin(theta);
  return tmp;
}// }}}

/**
 * \brief Applies the exponential function.
 *
 * From sakurai eq (3.2.44) we have that
 * exp(- i \sigma \cdot n \phi/2) =
 *  cos(\phi/2) \openone -  i sin(\phi/2) \sigma \cdot n
 *
 * exp(- i \sigma \cdot n \theta) =
 *  cos(\theta) \openone -  i sin(\theta) \sigma \cdot n
 *
 * \param b [itpp::vec]
 *
 * \return [itpp::cmat].
*/
itpp::cmat exp_minus_i_b_sigma(itpp::vec b){// {{{
  double theta=itpp::norm(b);
//   std::cout << theta << std::endl;
  std::complex<double> I(0,1);
  itpp::vec n=b/theta;
  return cos(theta)*sigma(0)-I*sin(theta)*sigma(n);
}// }}}

itpp::cmat magnetic(itpp::vec b, int total){// {{{
  itpp::cmat tmp(cfpmath::pow_2(total),cfpmath::pow_2(total));
  tmp=0.;
  for (int i=0; i<total; i++){
    tmp=tmp+sigma(b, i, total);
  }
  return tmp;
}// }}}
// Apply operations
template <class Num_T> void apply_sigma_x(itpp::Vec<Num_T>& state, int target_bit){// {{{

  int pos_1, pos_2;
  Num_T x;
  for (int i=0; i<state.size()/2; i++){
    pos_1 = cfpmath::merge_two_numbers(0, i, cfpmath::pow_2(target_bit));
    pos_2 = cfpmath::set_bit(pos_1,target_bit);
//     std::cout << "i="<< i << ", pos_1="<< pos_1 << ", pos_2="<< pos_2 << ", target_bit=" << target_bit << "\n";
    x=state(pos_1);
    state(pos_1)=state(pos_2);
    state(pos_2)=x;
//     std::cout << "state(pos_1)="<< state(pos_1) << ", state(pos_2)="<< state(pos_2)  << "\n";
  }
//    std::cout << "state" << state << "\n";
return;
}// }}}
template <class Num_T> void apply_sigma_y(itpp::Vec<Num_T>& state, int target_bit){// {{{
  int pos_1, pos_2;
  std::complex<double> Im(0.,1.);
  Num_T x;
  for (int i=0; i<state.size()/2; i++){
    pos_1 = cfpmath::merge_two_numbers(0, i, cfpmath::pow_2(target_bit));
    pos_2 = cfpmath::set_bit(pos_1,target_bit);
//     std::cout << "i="<< i << ", pos_1="<< pos_1 << ", pos_2="<< pos_2 << "\n";
    x=state(pos_1);
    state(pos_1)=-Im*state(pos_2);
    state(pos_2)=Im*x;
//     std::cout << "state(pos_1)="<< state(pos_1) << ", state(pos_2)="<< state(pos_2)  << "\n";
  }
//   std::cout << "state" << state << "\n";
}// }}}
template <class Num_T> void apply_sigma_z(itpp::Vec<Num_T>& state, int target_bit){// {{{
  int pos;
  Num_T x;
  for (int i=0; i<state.size()/2; i++){
    pos = cfpmath::merge_two_numbers(1, i, cfpmath::pow_2(target_bit));
//     pos = cfpmath::set_bit(target_bit,target_bit);
    state(pos)=-state(pos);
  }
//   std::cout << "state" << state << "\n";
}// }}}
// void apply_sigma(itpp::cvec state, itpp::vec b, int PositionQubit){ // {{{
//
//  abort();
// } // }}}
void apply_sigma(itpp::cvec& state, int Pauli_i, int PositionQubit){ // {{{
//   std::cout << "Por entrar en el lio PositionQubit=" << PositionQubit << ", Pauli_i=" << Pauli_i <<  std::endl;
  if (Pauli_i==1){
    apply_sigma_x(state, PositionQubit);
  } else if (Pauli_i==2) {
    apply_sigma_y(state, PositionQubit);
  } else if (Pauli_i==3) {

    apply_sigma_z(state, PositionQubit);
  }
  return;
} // }}}
void apply_sum_sigma_x(itpp::cvec& state){ // {{{
  // Applies the operator \sum_i \sigma_x^i to state
  itpp:: cvec state_orig=state, new_state=itpp::zeros_c(state.size());
  int qubits=cfpmath::log_base_2(state.size());
  for(int i=0;i<qubits;i++){
    apply_sigma_x(state_orig,i);
    new_state=new_state+state_orig;
    state_orig=state;
  }
  state=new_state;
  return;
} // }}}
void apply_gate(itpp::cvec& state, int nwhich, itpp::cmat gate){// {{{
  // This gate has to be applied with care. Notice for example that
  // with nwhich=3=11, the gate XY and YX are different. We are assuming
  // that the user has taken care of ordering the incoming gate
  // carefully. For example, if you want to apply H in qubit 0 and
  // Y gate in qubit 2, you sould set
  // nwhich = 1 0 1 = 5
  // and then set gate = YH
  // since HY will give you an incorrect result. Moreover, if you want
  // to apply a controled gate, it is not good to put
  // 1 (oplus) U
  // if the control qubit is ?????? I dont know, but there is a correct
  // and an incorrect way of putting it. !! see the testing routine
  // to set it up-.
  //
  int qs=cfpmath::BitCount(nwhich);
  int qs2=cfpmath::pow_2(qs);
  if ( (gate.rows() != qs2) || (gate.rows() != qs2)){
    std::cerr << "Error apply_gate" << std::endl;
    std::cerr << "gate.rows()=" << gate.rows() << std::endl;
    std::cerr << "qs2=" << qs2 << std::endl;
    std::cerr << "nwhich=" << nwhich << std::endl;
    abort();

  }
  itpp::cvec moco;
  itpp::ivec pos(qs2);
  for (int i=0; i<state.size()/qs2; i++){
    for (int j=0; j<qs2; j++){
      pos(j)=cfpmath::merge_two_numbers(j,i,nwhich);
    }
//     std::cout << "pos=" << pos << std::endl;
    moco=gate*state(pos);
//     std::cout << "moco=" << moco << std::endl;
    for (int j=0; j<qs2; j++){
      state(pos(j))=moco(j);
    }
//     std::cout << "state(pos)=" << state(pos) << std::endl;

  }
//   std::cout << "gate" << gate << std::endl;
//   std::cout << "este deberia estar rotado" << state << std::endl;
  return;
}// }}}
void apply_hadamard(itpp::cvec& state, int position){// {{{
  itpp::ivec pos(2);
  itpp::cvec moco;
  itpp::mat h=hadamard_matrix();
  for (int i=0; i<state.size()/2; i++){
    pos(0)=cfpmath::merge_two_numbers(0,i,cfpmath::pow_2(position));
    pos(1)=cfpmath::merge_two_numbers(1,i,cfpmath::pow_2(position));
    moco= h*state(pos);
    state(pos(0))= moco(0);
    state(pos(1))= moco(1);
  }
  return;
} // }}}
itpp::cmat multiply_by_sigma_leftmost_qubit(const itpp::cmat& A, int sigma_label){ // {{{
if (sigma_label==0){
  return A;
}
int n=A.cols();
itpp::cmat B(n,n);
if (sigma_label==1) {
  B.set_submatrix(0  ,0,A.get_rows(n/2,n-1  ));
  B.set_submatrix(n/2,0,A.get_rows(0,  n/2-1));
  return B;
} else if (sigma_label==2) {
  std::complex<double> I(0.,1.);
  B.set_submatrix(0  ,0,-I*A.get_rows(n/2,n-1  ));
  B.set_submatrix(n/2,0, I*A.get_rows(0,  n/2-1));

//   std::cout << "En la rutina A=" << A.get_rows(n/2,n-1) << std::endl;
//   std::cout << "En la rutina B=" << B(0,n/2-1,0,n-1) << std::endl;
  return B;

} else if (sigma_label==3) {
  B.set_submatrix(0  ,0, A.get_rows(0,n/2-1));
  B.set_submatrix(n/2,0,-A.get_rows(n/2,n-1));
//   B(n/2,n-1,0,n-1)=-A.get_rows(n/2,n-1);
//   B(0,n/2-1,0,n-1)=A.get_rows(0,n/2-1);
  return B;
} else {
  std::cerr << "sigma_label=" << sigma_label<<", not valid. aborting"<<std::endl;
  abort();
}


} // }}}
void apply_inverse_Rk(itpp::cvec& state, int k, int position1, int position2){// {{{
  std::complex<double> phi=exp(-std::complex<double>(0,1)*2.*itpp::pi/double(cfpmath::pow_2(k)));
  int mask = cfpmath::pow_2(position1) + cfpmath::pow_2(position2);
  int n;
  for (int j=0; j<state.size()/4; j++){
    n=cfpmath::merge_two_numbers(3, j, mask);
    state(n)=phi*state(n);
  }
  return;
} // }}}
void apply_swap(itpp::cvec& state, int position1, int position2){// {{{
  int mask = cfpmath::pow_2(position1) + cfpmath::pow_2(position2);
  int n1, n2;
  for (int j=0; j<state.size()/4; j++){
    n1=cfpmath::merge_two_numbers(2, j, mask);
    n2=cfpmath::merge_two_numbers(1, j, mask);
    swap(state, n1, n2);
  }
  return;
} // }}}
void apply_swap(itpp::cmat& rho, int position1, int position2){// {{{
  int mask = cfpmath::pow_2(position1) + cfpmath::pow_2(position2);
  int n1, n2;
  for (int j=0; j<rho.cols()/4; j++){
    n1=cfpmath::merge_two_numbers(2, j, mask); 
    n2=cfpmath::merge_two_numbers(1, j, mask); 
    rho.swap_rows(n1, n2);
    rho.swap_cols(n1, n2);
//     swap(state, n1, n2);
  }
  return;
} // }}}
void apply_cnot(itpp::cvec& state, int control, int target){// {{{
  int mask = cfpmath::pow_2(control) + cfpmath::pow_2(target);
  int n1, n2;
  for (int j=0; j<state.size()/4; j++){
    if(control < target){
      n1=cfpmath::merge_two_numbers(1, j, mask);
      n2=cfpmath::merge_two_numbers(3, j, mask);
    } else if (control > target) {
      n1=cfpmath::merge_two_numbers(2, j, mask);
      n2=cfpmath::merge_two_numbers(3, j, mask);
    }
    swap(state, n1, n2);
  }
  return;
} // }}}
void apply_control_u(itpp::cvec& state, int control, int target, itpp::cmat& u){// {{{
  int mask = cfpmath::pow_2(control) + cfpmath::pow_2(target);
  itpp::cvec moco;
//   int n1, n2;
  itpp::ivec pos(2);
  for (int j=0; j<state.size()/4; j++){
    if(control < target){
      pos(0)=cfpmath::merge_two_numbers(1, j, mask);
      pos(1)=cfpmath::merge_two_numbers(3, j, mask);
    } else if (control > target) {
      pos(0)=cfpmath::merge_two_numbers(2, j, mask);
      pos(1)=cfpmath::merge_two_numbers(3, j, mask);
    }
    moco=u*state(pos);
    for (int j=0; j<2; j++){
      state(pos(j))=moco(j);
    }
  }
  return;
} // }}}
// Quantum information functions on states
double InverseParticipationRatio(const itpp::cvec& psi){ // {{{
  double tmp=0.;
  for (int i=0; i<psi.size();i++){
    tmp+=std::pow(abs(psi(i)),4);
  }
  return tmp;

} // }}}

/**
 * \brief Calculates the concurrence.
 *
 * Calculates the concurrence of the specified density matrix.
 *
 * C = max{0, L1, L2, L3, L4}
 *
 * where Li are the eigenvalues of
 *
 * sqrt(rho (sigma_y ⊗ sigma_y) rho* (sigma_y ⊗ sigma_y))
 *
 * \param rho [itpp::Mat<std::complex<double>] density matrix.
 *
 * \return resulting value [double], concurrence of the density matrix.
*/
double Concurrence(const itpp::Mat<std::complex<double>> rho){ // {{{
	itpp::Mat<std::complex<double> > rho_tilde(4,4);
	for (unsigned int i1=0;i1<4;i1++){
		for (unsigned int i2=0;i2<4;i2++){
			rho_tilde(i1,i2)=conj(rho(3-i1,3-i2));
			if( ((cfpmath::BitCount(i1)+cfpmath::BitCount(i2) )%2) == 1){
				rho_tilde(i1,i2)=-rho_tilde(i1,i2);
			}
		}
	}
	rho_tilde=rho*rho_tilde;
	itpp::Vec<double> eigenvalues=real(eig(rho_tilde));
	if (max(eigenvalues)<0){
		return 0.;
	}
	else{
		return 2*sqrt(max(eigenvalues))-sum(sqrt(abs(eigenvalues)));
	}
} // }}}
double ConcurrenceFromPure(const itpp::cvec psi){ // {{{
	return 2*abs(psi(1)*psi(2)-psi(0)*psi(3));
} // }}}
double vonNeumann(const itpp::vec& lambda){// {{{
	double tmp=0;
	for (int i=0; i<lambda.size(); i++){tmp+=cfpmath::h_function(lambda(i));}
// 	return tmp/log(lambda.size());
	return tmp;
}// }}}
double vonNeumann(const itpp::cmat rho){// {{{
  itpp::Vec<double> eigenvalues=real(eig(rho));
  return vonNeumann(eigenvalues);
} // }}}

/**
 * \brief Calculates the purity.
 *
 * Calculates the purity of the specified density matrix.
 *
 * P(rho) = Tr rho^2
 *
 * \param rho [itpp::Mat<std::complex<double>] density matrix.
 *
 * \return resulting value [double], purity of the density matrix.
*/
double Purity(const itpp::Mat<std::complex<double> >& rho){// {{{
	double P=0;
	for (int i =0; i< rho.rows(); i++){
		P += abs(itpp::dot (rho.get_col(i),rho.get_row(i)));
	}
	return P;
}// }}}
/**
 * \brief Calculates the purity.
 *
 * Calculates the purity of the specified eigenvalues. Applying the dot product
 * to the eigenvalues received
 *
 * eigenvalues (dot) eigenvalues
 *
 * \param eigenvalues [itpp::vec] in the form of a vector.
 *
 * \return resulting value [double], purity of the eigenvalues.
*/
double Purity(const itpp::vec& eigenvalues){// {{{
        return itpp::dot(eigenvalues,eigenvalues);
}// }}}

// Others
itpp::mat LambdaMatrixQubitRight(itpp::cmat U){ // {{{ Quantum channel for a single qubit, slow
  itpp::mat Lambda(3,3);
  int n=U.cols();
//   U = exponentiate(-I*t*H);
  for (int j=1; j<=3; j++){
    for (int k=1; k<=3; k++){
      Lambda(j-1,k-1)=real(itpp::trace(
            itpp::kron(itpp::eye_c(n/2), sigma(j))*U*itpp::kron(itpp::eye_c(n/2), sigma(k))*itpp::hermitian_transpose(U)
            ))/n;
    }
  }
//   cout << Lambda(0,0) << endl;
  return Lambda;
} //}}}
itpp::mat LambdaMatrixQubitLeft(itpp::cmat U){ // {{{ Quantum channel for a single qubit, slow
  itpp::mat Lambda(3,3);
  int n=U.cols();
//   U = exponentiate(-I*t*H);
  for (int j=1; j<=3; j++){
    for (int k=1; k<=3; k++){
      Lambda(j-1,k-1)=real(itpp::trace(
            itpp::kron(sigma(j),itpp::eye_c(n/2))*U*itpp::kron(sigma(k),itpp::eye_c(n/2))*itpp::hermitian_transpose(U)
            ))/n;
    }
  }
//   cout << Lambda(0,0) << endl;
  return Lambda;
} //}}}
itpp::mat LambdaMatrixQubitLeft(itpp::cmat H, double t){ // {{{ Quantum channel for a single qubit, slow
  itpp::cmat W, U;
  std::complex<double> I(0.,1.);
  U = exponentiate(-I*t*H);
  return LambdaMatrixQubitLeft(U);
} //}}}
itpp::mat LambdaMatrixQubitLeft(itpp::cmat W, itpp::vec E, double t){ // {{{Quantum channel for a single qubit, eigenphases & Eigenenergies
  std::complex<double> I(0.,1.);
  itpp::cmat U_dagger=UDiagUdagger(W, exp(-I*t*E));
  itpp::cmat U=hermitian_transpose(U_dagger);
//   std::cout << U << std::endl;
  itpp::Array<itpp::cmat> sigmaU(3), sigmaUdagger(3);
  for (int j=1; j<=3; j++){
    sigmaU(j-1)      =multiply_by_sigma_leftmost_qubit(U,j);
    sigmaUdagger(j-1)=multiply_by_sigma_leftmost_qubit(U_dagger,j);
  }
  itpp::mat Lambda(3,3);
  int n=E.length();
  for (int j=1; j<=3; j++){
    for (int k=1; k<=3; k++){
      Lambda(j-1,k-1)=real(trAB(sigmaU(j-1),sigmaUdagger(k-1)))/n;
    }
  }
//   cout << Lambda(0,0) << endl;
  return Lambda;
} //}}}
itpp::Array<itpp::cvec> PrepareInitialStatesForSingleQubitTomography(){ //{{{

    itpp::Array<itpp::cvec> InitialStates(4);
    InitialStates(0)=itpp::to_cvec(itpp::vec_2(1.,0.)); //|0>
    InitialStates(1)=itpp::to_cvec(itpp::vec_2(0,1)); //|1>
    InitialStates(2)=itpp::to_cvec(itpp::vec_2(1,1))/sqrt(2); // |x+>=(|0>+|1>)/sqrt(2)
    InitialStates(3)=itpp::vec_2(std::complex<double>(1.,0.),std::complex<double>(0.,1.))/sqrt(2); // |y+>=(|0>+i|1>)/sqrt(2)
    return InitialStates;
} // }}}
itpp::mat TomographySingleQubit(itpp::Array<itpp::cmat>& rhos_f){ // {{{
  // The idea is to get lambda. We must first get the vector form
  // of each guy.
  itpp::Array<itpp::cvec> rhos_vectors(4);
  for (int i=0; i<4; i++){
    rhos_vectors(i)=OperatorToVectorPauliBasis(rhos_f(i));
  }
  itpp::mat Lambda(4,4);
  Lambda=0.;
  Lambda(0,0)=1;
  Lambda.set_col(0, real(rhos_vectors(0)+rhos_vectors(1))/2); // Lambda_i0 = (r^f_1 + r^f_2)/2
  Lambda.set_col(3, real(rhos_vectors(0)-rhos_vectors(1))/2); // Lambda_i3 = (r^f_1 - r^f_2)/2
  Lambda.set_col(1, real(rhos_vectors(2))-Lambda.get_col(0));  // Lambda_i1 = (r^f_1 - r^f_2)/2
  Lambda.set_col(2, real(rhos_vectors(3))-Lambda.get_col(0));  // Lambda_i2 = (r^f_1 - r^f_2)/2

  return Lambda;

} // }}}
itpp::cvec OperatorToVectorPauliBasis(itpp::cmat rho){ // {{{
  itpp::cvec r(4);
  for (int i=0; i<4; i++){
    r(i)=itpp::trace(sigma(i)*rho);
  }
  return r;
} // }}}
// }}}
// Quantum Evolution {{{
itpp::cvec GetState(itpp::cmat& U, itpp::vec& eigenvalues, itpp::cvec& psi_0, double t){ // {{{
  itpp::cmat U_dagger = itpp::hermitian_transpose(U);
  itpp::cvec psi_0_prime = U_dagger*psi_0;
  itpp::cvec psi_t_prime = itpp::elem_mult(exp(-eigenvalues * t * std::complex<double>(0.,1.)),psi_0_prime);
  itpp::cvec psi_t = U* psi_t_prime;

  //   cout << "In GetState, test_unitarity(U)=" <<  test_unitarity(U) << endl;
  //   cout << "In GetState, test_unitarity(U_dagger)=" <<  test_unitarity(U_dagger) << endl;
  //   cout << "In GetState, norm(psi_0)=" <<  norm(psi_0) << endl;
  //   cout << "In GetState, norm(psi_0_prime)=" <<  norm(psi_0_prime) << endl;
  //   cout << "In GetState, norm(psi_t)=" <<  norm(psi_t) << endl;
  //   abort();
  return psi_t;

} // }}}
itpp::cmat GetState(itpp::cmat& U, itpp::vec& eigenvalues, itpp::cvec& psi_0, double t, int dimension_central_system){ // {{{
  int d=eigenvalues.size();
  if ( U.cols()!= d || U.rows()!= d || psi_0.size() != d){
    std::cerr << "Error en GetState (for mixed)\n"
      << "U.cols()=" << U.cols() << ", U.rows()=" << U.rows()
      << ", eigenvalues.size()="<< eigenvalues.size() << ", psi_0.size()=" << psi_0.size()
      << std::endl;
    abort();
  }

  itpp::cvec tmp=GetState(U, eigenvalues, psi_0,t);
  itpp::cmat tmp_r = partial_trace(tmp, dimension_central_system);
  return partial_trace(GetState(U, eigenvalues, psi_0,t), dimension_central_system);
} // }}}
//  }}}
// Testing {{{
template <class Num_T>  double test_symmetry(const itpp::Mat<Num_T > H){ // {{{
 return itpp::norm(H-H.transpose() );
} // }}}
template <class Num_T>  bool sorted(const itpp::Vec<Num_T > v){ // {{{ Test if it is sorted, from smallest to greatest
  int n=v.size();
  if (n==0){return true;}
  for (int i=1; i<n; i++){
    if(v(i)<v(i-1)) return false;
  }
  return true;
} // }}}
void test_GetState_pure(){ // {{{
//   So the idea is to have a Hamiltonian, its eigenvectors and eigenvalues,
//   time and check if the evolution is good.
  std::cout << "In test_GetState" << std::endl;
  int dim=50;
  double t=itpp::randu();
//   t=1.;
//   t=itpp::pi/2;
  std::complex<double> I(0,1);
  itpp::cmat H, eigenvectors, eih;
  itpp::cvec psi_0, psi_slow, psi_fast;
  itpp::vec lambda;

  H=itpp::randn_c(dim,dim); H=H+H.hermitian_transpose();
//   H(0,0)=0.; H(1,1)=0.; H(0,1)=0.+1.*I; H(1,0)=std::conj(H(0,1));
//   H(0,1)=2.;
//   H(1,0)=2.;
  psi_0=itpp::randn_c(dim); psi_0=psi_0/norm(psi_0);
//   psi_0=0.; psi_0(0)=1.;
  eih=exponentiate(-I*t*H);
//   std::cout << "H=" << H << std::endl << "eih=" << eih << std::endl;
  psi_slow = exponentiate(-I*t*H)*psi_0;
  eig_sym(H,lambda, eigenvectors);
  psi_fast = GetState(eigenvectors, lambda, psi_0, t);
  std::cout << "1-Producto interior="
    << 1.-itpp::dot(conj(psi_fast), psi_fast) << std::endl;
//   std::cout << psi_slow << std::endl<< psi_fast << std::endl;

  abort();
  return;
//  return itpp::norm(H-H.transpose() );
} // }}}
void test_GetState_mixed(){ // {{{
//   So the idea is to have a Hamiltonian, its eigenvectors and eigenvalues,
//   time and check if the evolution is good.
//   std::cout << "In test_GetState_mixed" << std::endl;
  int dim_a=5, dim_b=2;
  int dim=dim_a*dim_b;
  double t=itpp::randu();
//   t=itpp::pi/2;
  std::complex<double> I(0,1);
  itpp::cmat H, eigenvectors, eih, rho_fast;
  itpp::cvec psi_0, psi_slow, psi_fast;
  itpp::vec lambda;

  H=itpp::randn_c(dim,dim); H=H+H.hermitian_transpose();
//   H(0,0)=0.; H(1,1)=0.; H(0,1)=0.+1.*I; H(1,0)=std::conj(H(0,1));
//   H(0,1)=2.;
//   H(1,0)=2.;
  psi_0=itpp::randn_c(dim); psi_0=psi_0/norm(psi_0);
//   psi_0=0.; psi_0(0)=1.;
//   eih=exponentiate(-I*t*H);
//   std::cout << "H=" << H << std::endl << "eih=" << eih << std::endl;
//   psi_slow = exponentiate(-I*t*H)*psi_0;
  eig_sym(H,lambda, eigenvectors);
//   psi_fast = GetState(eigenvectors, lambda, psi_0, t);
  rho_fast = GetState(eigenvectors, lambda, psi_0, t, dim_a);
//   std::cout << "tr rho=" << itpp::trace(rho_fast) << std::endl;
//   std::cout << psi_slow << std::endl<< psi_fast << std::endl;

//   abort();
  return;
//  return itpp::norm(H-H.transpose() );
} // }}}
double test_reality(const itpp::cmat& H){ // {{{
 return itpp::norm(itpp::imag(H));
} // }}}
double test_real_symmetric(const itpp::cmat H){ // {{{
 return test_symmetry(H) + test_reality(H);
} // }}}
double test_unit_matrix(const itpp::cmat U){ // {{{
  if (U.cols() != U.rows()){
    std::cout << "La matriz no es cuadrada\n";
    abort();
    return -1;
  }
 return itpp::norm(U - itpp::eye_c(U.rows() ));
} // }}}
double test_unitarity(const itpp::cmat U){ // {{{
  if (U.cols() != U.rows()){
    std::cout << "La matriz no es cuadrada\n";
    abort();
    return -1;
  }
 return itpp::norm(U*U.hermitian_transpose() - itpp::eye_c(U.rows() ))
          + itpp::norm(U.hermitian_transpose()*U - itpp::eye_c(U.rows() ));
} // }}}
void  test_multiply_by_sigma_leftmost_qubit(){ // {{{
  int q=4;
  int n=cfpmath::pow_2(q);
//   n=30;
  itpp::cmat test_matrix;
  itpp::Array<itpp::cmat> test_matrix_results(4);
  itpp::Array<itpp::cmat> test_matrix_results_long(4);
  itpp::randn_c (n, n, test_matrix);
  double error=0.;
  for (int i =0; i<4; i++){
    test_matrix_results(i)=multiply_by_sigma_leftmost_qubit(test_matrix,i);
  }
  test_matrix_results_long(0)=test_matrix;
  test_matrix_results_long(1)=sigma(itpp::vec("1 0 0"),q-1,q)*test_matrix;
  test_matrix_results_long(2)=sigma(itpp::vec("0 1 0"),q-1,q)*test_matrix;
  test_matrix_results_long(3)=sigma(itpp::vec("0 0 1"),q-1,q)*test_matrix;
//   std::cout << "test_matrix_results_long:" << std::endl << test_matrix_results_long(2) << std::endl;
//   std::cout << "test_matrix_results:" << std::endl << test_matrix_results(2) << std::endl;
//   std::cout << "Carga y control 10" << std::endl;
  for (int i =0; i<4; i++){
//     std::cout << "Carga y control i=" << i << ", error acumulado=" << error <<  std::endl;
//     std::cout << test_matrix_results(i) << std::endl;
    error+=norm(test_matrix_results(i)-test_matrix_results_long(i));
  }
  std::cout << "Error in multiply_by_sigma_leftmost_qubit: " << error<< std::endl;
  return;
} // }}}
void test_UDiagUdagger(){ // {{{
  int n=5;
  itpp::cmat U; itpp::randn_c(n, n, U);
  itpp::cvec e; itpp::randn_c(n, e);
  itpp::cmat W1 = UDiagUdagger(U, e);
  itpp::cmat W2 = U*diag(e)*hermitian_transpose(U);
  std::cout << "error en uno y otro metodo=" << norm(W1-W2) << std::endl;
  return;
} // }}}
void test_trAB(){ // {{{
  std::complex<double> r1, r2;
  int n1=15, n2=8;
  itpp::cmat A, B; itpp::randn_c(n1, n2, A);itpp::randn_c(n2, n1, B);

//   B=1;
//   A=3;
  r1 = trAB(A, B);
  r2 = itpp::trace(A*B);
//   std::cout << "Test result A*B " << A*B << std::endl;
//   std::cout << "Test result " << r1 << ", " <<  r2 << std::endl;
  std::cout << "Test result " << abs(r1-r2) << std::endl;
} // }}}
//  }}}
// Not yet ordered {{{

template <class Num_T> itpp::Mat<Num_T > Array_to_Matrix(const itpp::Array<itpp::Vec<Num_T > >& v){ // {{{
//   std::cout << "Lo que entra:" << v << std::endl;
  int n = v.length();
  itpp::Mat<Num_T > tmp(n, v(0).size());
  for (int i=0; i<n; i++){
    tmp.set_row(i,v(i));
  }
  return tmp;
} //}}}
template <class Num_T> itpp::Array<itpp::Vec<Num_T> > Matrix_to_Array(const itpp::Mat<Num_T >& v){ // {{{
  int n = v.rows();
  itpp::Array<itpp::Vec<Num_T> > tmp(n);
  for (int i=0; i<n; i++){
    tmp(i) = v.get_row(i);
  }
  return tmp;
} //}}}
template <class Num_T> itpp::Vec<Num_T > my_mod(const itpp::Vec<Num_T >& v){ // {{{
itpp::Vec<Num_T > tmp(v.size());
  for (int i=0; i<v.size(); i++){
    tmp(i) = cfpmath::my_mod(v(i));
  }
  return tmp;
} //}}}
template <class Num_T> itpp::Array<Num_T > my_mod(const itpp::Array<Num_T >& v){ // {{{
itpp::Array<Num_T > tmp(v.size());
  for (int i=0; i<v.size(); i++){
    tmp(i) = cfpmath::my_mod(v(i));
  }
  return tmp;
} //}}}
itpp::cvec read_state(int d){ // {{{
  itpp::cvec psi(d);
  double x,y;
  for (int j=0; j<d; j++){
    std::cin >> x  >> y ;
    psi(j)=x+std::complex<double>(0,1)*y;
  }
  return psi;
} // }}}
itpp::Vec<std::string> split_string(std::string sentence){// {{{
  // http://stackoverflow.com/questions/236129/c-how-to-split-a-string
  std::string word;
  std::istringstream iss(sentence);
  int number_words=0;
  while( iss >> word ){number_words++;}
  itpp::Vec<std::string> splitted(number_words);
  number_words=0;
  std::istringstream iss2(sentence);
//   iss(sentence);
  while( iss2 >> word ){splitted(number_words)=word; number_words++;}
  return splitted;
}// }}}
itpp::Mat<std::complex<double> > Hermitian_From_Compact(const itpp::vec& compact_rho){// {{{
	int dim=cfpmath::isqrt(compact_rho.size()), counter=0;
	itpp::Mat<std::complex<double>  > tmp(dim,dim);
	for (int i=0;i<dim;i++){
		tmp(i,i)=compact_rho(counter); counter++;
		for (int j=i+1;j<dim;j++){
// 			tmp(i,j)=compact_rho(counter)+compact_rho(counter+1);
			tmp(i,j)=std::complex<double>(compact_rho(counter),compact_rho(counter+1));
			tmp(j,i)=conj(tmp(i,j));counter+=2;
		}
	}
	return tmp;
}// }}}
itpp::Array< itpp::ivec > LinSpace(const int i_i, const int i_f, const int step=1){// {{{
  int index=0;
  itpp::Array< itpp::ivec > tmp((i_f-i_i)/step + 1);
  for (int i=i_i; i<= i_f; i+=step){
//     std::cerr << "Creating the index " << index << "  with value " << i << "\n";
    tmp(index)=itpp::vec_1(i);
    index++;
  }
//   std::cerr << tmp  << " es este\n";
  return tmp;
}// }}}
template <class Num_T> itpp::ivec waist(const itpp::Array<Num_T >& A){// {{{
  itpp::ivec tmp(A.size());
  for (int i = 0; i< A.size(); i++){
    tmp(i) = A(i).size();
  }
  return tmp;
}// }}}
// }}}
// Pretty printing {{{
void standard_printing(itpp::cvec& c){ // {{{

  for (int i=0; i<c.size(); i++){
    std::cout <<  real(c(i)) << std::endl;
    std::cout <<  imag(c(i)) << std::endl;
  }
  return;
} // }}}
void standard_printing(itpp::cmat& c){ // {{{

  for (int i=0; i<c.rows(); i++){
    for (int j=0; j<c.cols(); j++){
      std::cout <<  real(c(i,j)) << std::endl;
      std::cout <<  imag(c(i,j)) << std::endl;
    }
  }
  return;
} // }}}
void standard_printing(itpp::vec& c){ // {{{

  for (int i=0; i<c.size(); i++){
    std::cout <<  c(i) << std::endl;
//     std::cout << i << "  " <<  c(i) << std::endl;
  }
  return;
} // }}}
void pretty_printing_1(itpp::cvec& c){ // {{{

  for (int i=0; i<c.size(); i++){
    std::cout << i << " " << c(i) << std::endl;
  }
  return;
} // }}}
// }}}
} // }}}
#endif // ITPP_EXT_MATH_VARIOUS
