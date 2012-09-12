#ifndef QuantumStateITPP_CLASS
#define QuantumStateITPP_CLASS


// #include <QuantumStateITPP.h>
template <class StateType> StateType StateTensor(const StateType& psi_1, const StateType& psi_2){
	StateType tmp(psi_1.size()*psi_2.size());
	for (int i =0; i<psi_1.size(); i++){
		tmp.coefficients.set_subvector(i*psi_2.size(),psi_1.coefficients(i)*psi_2.coefficients);
	}
	return tmp;
}
class QuantumStateITPP{
	public:
	itpp::Vec<std::complex<double> > coefficients;
	QuantumStateITPP(unsigned int dimension, std::string state_type="Random"){
		coefficients.set_size(dimension);
		if (state_type=="Random"){
			InitializeRandom( ) ;
		}
		else if (state_type=="Separable"){
			InitializeQubitSeparable( );	
		}
		else if (state_type=="BellRandom"){
			BellRandom( );	
		}
		else if (state_type=="Bell"){
			coefficients=0.;
			coefficients(0)=1;
			coefficients(dimension-1)=1;
			normalize();
		} else {
			std::cerr  << "Illegal type initializing a QuantumState";
			exit(1);
		}
	}
	QuantumStateITPP(std::string state_type="Bell"){
		if (state_type=="Bell"){
			coefficients.set_size(4);
			coefficients="1 0 0 1";
			normalize();
		} else {
			std::cerr  << "Illegal type initializing a QuantumState";
			exit(1);
		}
	}
	/// Method to have the size of the quantum state
	int size( ) const{		
		return coefficients.size() ;
	}
	void normalize( ){
		coefficients/= norm() ;		
	}
	double norm( ) const{		
		double r=0;
		for (int i=0; i< size();i++){
			r+= std::pow( abs(coefficients(i)), 2 );
		}
		return std::sqrt(r) ;
	}
	/// A method to reinitialize to a random state
	void InitializeRandom( ){
		coefficients=itpp::randn_c(size());
		normalize( ) ;
	}
	void InitializeQubitSeparable( ){
		QuantumStateITPP onequbit(2);
		QuantumStateITPP tmp(1);
		tmp.coefficients=1;
		for (int q=1; q<size() ; q*=2){
			onequbit.InitializeRandom();
			tmp=StateTensor(tmp,onequbit);
		}
		coefficients=tmp.coefficients;
		return;
	}
	void BellRandom( ){
		coefficients=(StateTensor(QuantumStateITPP(4,"Bell"),
				QuantumStateITPP(size()/4,"Random"))).coefficients;
		return;
	}
	void set(const QuantumStateITPP state_in){
		coefficients.set_size(state_in.size() );
		coefficients=state_in.coefficients;
	}
	itpp::Mat<std::complex<double> > PartialTrace(const int size_H_1){
  		int size_H_2=size() /size_H_1;
		itpp::Mat<std::complex<double> > rho(size_H_1,size_H_1);
		itpp::Vec<std::complex<double> > vec_row(size_H_2),vec_col(size_H_2);
		if (size_H_2*size_H_1 != size()){
			std::cerr  << "error PartialTrace" ; 
			exit(1);
		}
		for (int j_row=0 ; j_row< size_H_1 ; j_row++){
		vec_row=coefficients.mid(size_H_2*j_row,size_H_2);
			for (int j_col=j_row; j_col< size_H_1 ; j_col++){
				vec_col=coefficients.mid(size_H_2*j_col,size_H_2);
				rho(j_row,j_col)=dot(itpp::conj(vec_row),vec_col);
    				rho(j_col,j_row)= conj(rho(j_row,j_col));  
			}
		}
		return rho ;
	}
};
std::ostream &operator<<(std::ostream &os, const QuantumStateITPP  &psi){
os << psi.coefficients;return os;}
QuantumStateITPP operator*(const itpp::Mat<std::complex<double> > &H, const QuantumStateITPP& psi){
	QuantumStateITPP tmp;
	tmp.coefficients=H*psi.coefficients;
	return tmp;
}
QuantumStateITPP evolve_with_phases(QuantumStateITPP &psi, itpp::Vec<double> eigenvalues, double t){
// 	psi=exp(-I*time*eigenvalues) * psi_0
	QuantumStateITPP tmp;
	tmp.coefficients = elem_mult(itpp::exp(eigenvalues * t * std::complex<double>(0.,1.)),psi.coefficients) ;
	return tmp;
}
#endif // QuantumStateITPP
