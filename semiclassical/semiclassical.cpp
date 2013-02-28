#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
#include <fftw3.h>
namespace semiclassical{ // {{{
// FFT initialize  {{{
//********************************************//
  fftw_complex *in;
  fftw_complex *out;
  fftw_plan planb;
  fftw_plan planf;
  //
//********************************************//
void initialize_fftplans(int dim){ 
	in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*dim);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*dim);
	planf=fftw_plan_dft_1d(dim,in,in,FFTW_FORWARD,FFTW_ESTIMATE);
	planb=fftw_plan_dft_1d(dim,out,out,FFTW_BACKWARD,FFTW_ESTIMATE);
        return;
} // }}}
void copyF(itpp::cvec &a,fftw_complex *in){ // {{{
  int nx=size(a);
int i,j;
double xi,xj;
for(i=0;i<nx;i++)
   {
	   xi=std::real(a(i));
	   xj=std::imag(a(i));
	   in[i][0]=xi;
	   in[i][1]=xj;
   }	
return;
} // }}}
void copyB(itpp::cvec &a,fftw_complex *out){ // {{{
  int nx=a.size();
int i,j;
double xi,xj;
for(i=0;i<nx;i++)
     {
	   xi=out[i][0];
	   xj=out[i][1];
	   a(i)=std::complex<double>(xi,xj);
     }
return;
} // }}}
void FFT(itpp::cvec &a,int dir){ // {{{
int nq_tot,nq2,N,nx,ny;
nq_tot=log2(a.size());
nq2=nq_tot/2;
N=1<<(nq_tot);
nx=ny=1<<(nq2);
double norma=norm(a);
if(dir==0)
  {
    copyF(a,in);
    fftw_execute(planf);
    copyB(a,in);
  }
if(dir==1)
  {
    copyF(a,out);
    fftw_execute(planb);
    copyB(a,out);
  }

a=a/norm(a);
return;
} //}}}
itpp::cvec coherent_state(double q0,double p0, double sigma, int N){ // {{{
//  int 
  itpp::cvec a;
  a.set_size(N);
 a(0)=1.;
//  int N=a.size();
	int jmax=8;
	double arg=2*M_PI*N;
	std::complex<double> sumc;
// 	double sigma=coefsigma;
	double tip,tiq;
	tip=tiq=0.;
	for(int i=0;i<N;i++){
				sumc=std::complex<double>(0.,0.);
				for(int jj=-jmax;jj<=jmax;jj++){
					double dx=jj+q0-((double)i+tiq)/((double)N);
				//	double ex=arg*(-0.5*pow(dx/sigma,2.)-complex(0.,p0*dx-jj*tip/((double)N)));
				double ex=arg*-0.5*pow(dx/sigma,2.);
				double ex1=arg*p0*dx-jj*tip/((double)N);
					sumc+=exp(ex)*std::complex<double>(cos(ex1),sin(ex1));		
				}
 				a(i)=sumc;
// 				a(i)=1.;
	}
// 	a=a/itpp::norm(a);
	a=a/itpp::norm(a);
    return a;
} // }}}
itpp::cvec coherent_state(double q0,double p0, int N){ // {{{
  return coherent_state(q0, p0, 1., N);
} //}}}
// Harper map {{{
/***********************************/
void U_x_har(double k,itpp::cvec &a){ // {{{
  int l,N=a.size();
  double x;
  double theta;
  itpp::cvec phases(N);
  for(l=0;l<N;l++){
    theta=2*M_PI*l/((double)N);
    x=2*M_PI*N*k*cos(theta);
    phases[l]=std::complex<double>(cos(x),sin(x));
  }
  a=elem_mult(a,phases);
  return ;
} // }}}
void U_p_har(double k, itpp::cvec &a){ // {{{
  int l,N=a.size();
  double x,theta;
  itpp::cvec phases(N);
  for(l=0;l<N;l++){
    theta=2*M_PI*l/((double)N);
    x=2*M_PI*N*k*cos(theta);
    phases[l]=std::complex<double>(cos(x),sin(x));
  }
  a=elem_mult(a,phases);
  return ;
} // }}}
void kick_har(double kx, double kp, itpp::cvec &a){ // {{{
  U_x_har(kx,a);
  FFT(a,1);
  U_p_har(kp,a);
  FFT(a,0); 
} // }}}
// }}}
// Sawtooth map {{{
//  q'=q+p'
//  p'=p+K(q mod 1 -0.5)
void U_p_saw(double T,itpp::cvec &a){ // {{{
  double x;
  int l, N=a.size();
  itpp::cvec phases(N);
  for(l=0;l<N;l++){
    x=(-2*M_PI*(0.5*l*l/((double)N)));
    phases[l]=std::complex<double>(cos(x),sin(x));
  }
  a=elem_mult(a,phases);
} // }}}
void U_x_saw(double T, itpp::cvec &a){ // {{{
  int l,N=a.size();
  double x;
  itpp::cvec phases(N);
  if (N%2 != 0){
    std::cerr << "The dimension must be even. In this case, dimension=" 
      << N << std::endl<<"Aborting @U_x_std " << std::endl;
    abort();
  }
  for(l=0;l<N;l++){
    x=2*M_PI*(T)*(0.5*(l-N/2)*(l-N/2)/((double)N)); 
    phases[l]=std::complex<double>(cos(x),sin(x));
  }
  a=elem_mult(a,phases);
  return ;
} // }}}
void kick_saw(double T,itpp::cvec &a){ // {{{
  double tn=2*itpp::pi*T;
//   std::cout << "En el saw, T=" << T << std::endl; 
//   std::cout << "En el saw, 1 psi(0)=" << a(0) << std::endl; 
//   U_x_saw(T,a);
//   FFT(a,1);
//   U_p_saw(0.,a);
//   FFT(a,0); 
//   std::cout << "En el saw, T=" << T << std::endl; 
//   std::cout << "En el saw, 1 psi(0)=" << a(0) << std::endl; 
  U_x_saw(tn,a);
//   std::cout << "En el saw, 2 psi(0)=" << a(0) << std::endl; 
  FFT(a,1);
//   std::cout << "En el saw, 3 psi(0)=" << a(0) << std::endl; 
  U_p_saw(0.,a);
//   std::cout << "En el saw, 4 psi(0)=" << a(0) << std::endl; 
  FFT(a,0); 
//   std::cout << "En el saw, 5 psi(0)=" << a(0) << std::endl; 
} // }}}
// }}}
// standard map {{{
// p'=p+k sin(q) mod(2 pi) transicion en k~0.96
// q'=q+p'
void U_p_std(double T,itpp::cvec &a){ // {{{
  int l,N;
  double x;double T1=2.0; // de acuerdo a valores tipicos usados con dima
  N=a.size();
  if (N%2 != 0){
    std::cerr << "The dimension must be even. In this case, dimension=" << N 
      << std::endl<<"Aborting @U_p_std " << std::endl;
    abort();
  }
  itpp::cvec phases(N);
  for(l=0;l<N;l++){
    // Correccion de nacho, el menos no va!!!
    // x=-M_PI*(double)l*(double)l/((double)N);
    x=M_PI*(double)l*(double)l/((double)N);
    if(T>=0.){
      phases[l]=std::complex<double>(cos(x),sin(x));
    }else{
      phases[l]=std::complex<double>(cos(x),-sin(x));
    }
  }
  a=elem_mult(a,phases);
} // }}}
void U_x_std(double k, itpp::cvec &a){ //{{{
  int l,N;
  double x;
  double theta;
  // std::complex<double> *phases;
  double alpha2=0.;
  N=size(a);
  if (N%2 != 0){
    std::cerr << "The dimension must be even. In this case, dimension=" << N 
      << std::endl<<"Aborting @U_x_std " << std::endl;
    abort();
  }
  itpp::cvec phases(N);
  for(l=0;l<N;l++){
    theta=2*M_PI*l/((double)N);
    // Correccion de nacho, el menos no va!!!
    // 	x=(-2*M_PI*N*k*cos(theta));
    x=(2*M_PI*N*k*cos(theta));
    phases(l)=std::complex<double>(cos(x),sin(x));
  }
  a=elem_mult(a,phases);
  return ;
} // }}}
void kick_std(double T,itpp::cvec &a){ // {{{
//         std::cout << "En kick_std, 1 a(0)=" << a(0) << std::endl;
	U_x_std(T,a);
//         std::cout << "En kick_std, 2 a(0)=" << a(0) << std::endl;
	FFT(a,1);
	U_p_std(2.0,a);
	FFT(a,0);
} // }}}
// }}}
void print_state_p(itpp::cvec &a,char *filename){ // {{{
FILE *fp;
double pp;


fp=fopen(filename,"w");
for(int i=0;i<a.size();i++){
    pp=(double)i;
	fprintf(fp,"%lf\t %.14e\t %.14e\t %.14e\n",pp,abs(a(i))*abs(a(i)),real(a(i)),imag(a(i)));
   }	
 //fprintf(fp,"%lf\t %20.18lf \n",(double)NN,absquad(a(0))); //real_teil(a(0)),imag_teil(a(0)));
fprintf(fp,"\n");
fclose(fp);
return;
} // }}}
} // }}}
