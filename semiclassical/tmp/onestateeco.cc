// #include <itpp_ext_math.h>
// #include <cfp_math.cpp>
#include <itpp/itbase.h>
#include <itpp/stat/misc_stat.h>
#define EXP2(x) (1<<(x))
// #include "qstate.h"
#include <fftw3.h>
// #include <stdlib.h>
using namespace std;
// double px0,py0,x0,Y0;
// double fx,Fx,Fy,gammax,gammay;
double alpha_up;

//********************************************//
  fftw_complex *in;
  fftw_complex *out;
  fftw_plan planb;
  fftw_plan planf;
//********************************************//
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
N=EXP2(nq_tot);
nx=ny=EXP2(nq2);
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
void put_coh(itpp::cvec& a,double p0,double q0, double sigma){ // {{{
//  int 
 a(0)=1.;
 int N=a.size();
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
    return;
} // }}}
// standard map
// p'=p+k sin(q) mod(2 pi) transicion en k~0.96
// q'=q+p'
void U_p_std(double T,itpp::cvec &a){ // {{{
    int l,N;
    double x;double T1=2.0; // de acuerdo a valores tipicos usados con dima
    N=a.size();
    itpp::cvec phases(N);
    for(l=0;l<N;l++){
	x=-M_PI*(double)l*(double)l/((double)N);
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
double x,k1;
double theta,fase;
// std::complex<double> *phases;
k1=1.1;
double del=0.9;
double alpha1=1.;
double alpha2=0.;
fase=0.0; // *rand_double();
    N=size(a);
//     phases=new std::complex<double>[N+5];
    itpp::cvec phases(N);
    for(l=0;l<N;l++){
	theta=2*M_PI*l/((double)N);
	if(k>alpha_up){
	x=(-2*M_PI*N*alpha_up*(alpha1*cos(theta)+alpha2*cos(2.*theta)/2.))+
	(-2*M_PI*N*(k-alpha_up)*cos(theta));
	}else{	
	x=(-2*M_PI*N*k*(alpha1*cos(theta)+alpha2*cos(2.*theta)/2.));}
	//(N*k1*cos((double)theta)+N*k*(sin(theta)-0.5*sin(2*theta)));
	//x=-2*T*theta;
	phases(l)=std::complex<double>(cos(x),sin(x));
// 	if (l==10) cout << l << ", " << phases[l] << "; a(l)="<< a(l) << ", a(l)*phases[l]=" << a(l)*phases[l] <<  endl;
    }
//     cout << "moco moco " << phases << "moco moco " << endl;
    a=elem_mult(a,phases);
//     cout << "a ver " << a(10) << endl;
//     delete[] phases;
return ;
} // }}}
void kick_std(double T,itpp::cvec &a){ // {{{
	U_x_std(T,a);
// 	cout << "En kick_std, a(0)=" << a(0) << endl;
	FFT(a,1);
	U_p_std(2.0,a);
	FFT(a,0);
} // }}}
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
void inicio_estado(itpp::cvec &a, double q0, double p0, double sigma){ // {{{
	put_coh(a,p0,q0, sigma);
	a=a/itpp::norm(a);
	return;
} // }}}
int main(){ // {{{
	double q0=0.3,p0=0.4;
	int nqb=50, tmax=10;
	alpha_up=4;
	double Dt=2.4048;
	double sigma=40;
	itpp::cvec ai(nqb), bi(nqb);

	int dimq=nqb;
	alpha_up=alpha_up/(4.*M_PI*M_PI);
	Dt=Dt/(2.*M_PI*dimq); //

	//Dt=Dt/(2.*M_PI*dimq);
	// fftw plan initialize ---------------------------------------- {{{
	in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*dimq);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*dimq);
	planf=fftw_plan_dft_1d(dimq,in,in,FFTW_FORWARD,FFTW_ESTIMATE);
	planb=fftw_plan_dft_1d(dimq,out,out,FFTW_BACKWARD,FFTW_ESTIMATE);
	//-------------------------------------------------------------- }}}
	inicio_estado(ai,q0,p0, sigma);
// 	std::cout << ai << std::endl;
	bi=ai;
	for(int t=0;t<=tmax;t++){
// 		cout << a(0) << ", " << ai(0) << endl; 
		kick_std(alpha_up,ai);
// 		cout << a(0) << ", " << ai(0) << endl; 
// 		abort(); 
		kick_std(alpha_up+Dt,bi);		
	} //end t loop
	std::cout << ai << std::endl;
// 	std::cout << bi << std::endl;
	return 0;
} // }}}
