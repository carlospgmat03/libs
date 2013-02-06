#include "semiclassical.cpp"
int main(){ // {{{
	double q0=0.3,p0=0.4;
	int dim=50, tmax=10;
	double alpha_up=4;
	double Dt=2.4048;
	double sigma=40;
	itpp::cvec ai, bi;
	std::cout.precision(16);
	alpha_up=alpha_up/(4.*M_PI*M_PI);
	double normalizedDt=Dt/(2.*M_PI*dim); //

	semiclassical::initialize_fftplans(dim);
	semiclassical::inicio_estado(ai,q0,p0, sigma, dim);
	bi=ai;
	for(int t=0;t<=tmax;t++){
		semiclassical::kick_std(alpha_up,ai);
		semiclassical::kick_std(alpha_up+normalizedDt,bi);		
	} 
	for (int i=0; i<dim; i++){
		std::cout << real(ai(i)) << " " << imag(ai(i)) << std::endl;
	}
	return 0;
} // }}}
