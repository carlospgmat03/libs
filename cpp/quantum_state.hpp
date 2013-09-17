#include <valarray>
#include <iostream>
#include <complex>
#include <cmath> 
#include "my_random.hpp"
#include <string>

// http://msdn2.microsoft.com/en-us/library/7a40a0e3(VS.80).aspx
/** \brief A class to hold discreate quantum states
 *
 * This class allows to define a quantum state with arbitrary finite
 * dimension. Its can normalize, randomize, and print
 */
// Sugiere David que use un parametro de plantilla:
// template<class wp> class QuantumState { }
class QuantumState
{
	public:
	/// The coeficients, complex numbers
	std::valarray<std::complex<wp> > coefficients ;
	/// The constructor. Initialized to a random state by default
	QuantumState(unsigned int dimension, std::string state_type="Random"){
		coefficients.resize(dimension);
		if (state_type=="Random"){
			InitializeRandom( ) ;
		} else {
			std::cout << "Illegal type initializing a QuantumState";
			exit(1);
		}
	}
	/// Method to have the size of the quantum state
	unsigned int size( ) const{		
		return coefficients.size() ;
	}
	/// Method to get the norm of the quantum state
	wp norm( ) const{		
		wp r=0;
		for (int i=0; i< size();i++){
			// remplazar aca pow por multiplicar la expresion por si
			// misma o usar la de blitz
			r+= std::pow( abs(coefficients[i]), 2 );
		}
		return std::sqrt(r) ;
	}
	/// A method to reinitialize to a random state
	void InitializeRandom( ){
		for (int i=0; i < coefficients.size( ); i++ ){
			coefficients[i]=random_variables::ComplexSampleNormal();
		}
		coefficients/= norm() ;
	}
	/// prints the state in a nice way
	void print( ){
		for (int i=0; i < coefficients.size( ); i++ ){
			std::cout << "[" << i << "]=" << coefficients[i] << std::endl ;
		}
	}
};
///Dot product for two different states
std::complex<wp> dot_product(  const QuantumState& psi_left, const QuantumState& psi_right){
	std::complex<wp> r;
	r=0;
	for (int i=0; i < psi_left.size( ) ; i++){
		r+= conj(psi_left.coefficients[i]) * psi_right.coefficients[i];
	}
	return r ;
}


