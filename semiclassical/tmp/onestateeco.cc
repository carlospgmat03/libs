/*************************************************/
//
#include "complex.h"
#include "random.h"
#include "qstate.h"
#include "math.h"
#include <time.h>
#include <fftw3.h>

int mapa;
double px0,py0,x0,Y0;
double aa,V0,AA,Dt,tmax,dim;
double fx,Fx,Fy,gammax,gammay;
double fy,coefsigma,sat_t;
time_t seconds00,seconds0,seconds1,duracion;
double w=1.5;
int posi;
double alpha_up;
double rleft,rright;
double ma,mb;

//********************************************//
  fftw_complex *in;
  fftw_complex *out;
  fftw_plan planb;
  fftw_plan planf;
//********************************************//

void copyF(q_state &a,fftw_complex *in){ // {{{
  int nx=a.N();
int i,j;
double xi,xj;
for(i=0;i<nx;i++)
   {
	   xi=real_teil(a(i));
	   xj=imag_teil(a(i));
	   in[i][0]=xi;
	   in[i][1]=xj;
   }	
return;
} // }}}
void copyB(q_state &a,fftw_complex *out){ // {{{
  int nx=a.N();
int i,j;
double xi,xj;
for(i=0;i<nx;i++)
     {
	   xi=out[i][0];
	   xj=out[i][1];
	   a(i)=complex(xi,xj);
     }
return;
} // }}}
void normalize(q_state &a,double norm){ // {{{
int N=a.N();
for(int i=0;i<N;i++){
		a(i)*=norm;
	}
return;	
} // }}}
void FFT(q_state &a,int dir){ // {{{
int nq_tot,nq2,N,nx,ny;
nq_tot=a.L();
nq2=nq_tot/2;
N=EXP2(nq_tot);
nx=ny=EXP2(nq2);
double norma=a.norm();
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

a.normalize();
normalize(a,norma);
return;
} //}}}
void put_point_pos(q_state &a,double q0){ // {{{
    double p0,sigma2;
    int i,num,N1;
double x=a.N()/8.;
 N1=a.N();
    a.put_to_zero();
    sigma2=(double)a.N()/x;
    p0=0; //a.N()/2.;
   // q0=2.0; //M_PI; //a.N()/2.; //M_PI; //2.0;
    // a.gauss_add(p0,q0,sigma2);
    //a(N1/2)=1.0; 
	posi=(int)(q0*N1);

    a(posi)=1.0;
    //a(3*N1/4)=1.0;
    a.normalize();
} // }}}
void put_coh(q_state &a,double p0,double q0){ // {{{
    int p,N;
	N=a.N();
	int jmax=8;
	double arg=2*M_PI*N;
	complex sumc;
	double sigma=1.*coefsigma;
	double tip,tiq;
	tip=tiq=0.;
	for(int i=0;i<N;i++){
				sumc=complex(0.,0.);
				for(int jj=-jmax;jj<=jmax;jj++){
					double dx=jj+q0-((double)i+tiq)/((double)N);
				//	double ex=arg*(-0.5*pow(dx/sigma,2.)-complex(0.,p0*dx-jj*tip/((double)N)));
				double ex=arg*-0.5*pow(dx/sigma,2.);
				double ex1=arg*p0*dx-jj*tip/((double)N);
					sumc+=exp(ex)*complex(cos(ex1),sin(ex1));		
				}
				a(i)=sumc;
	}
	a.normalize();
    return;
} // }}}

// standard map
// p'=p+k sin(q) mod(2 pi) transicion en k~0.96
// q'=q+p'
void U_p_std(double T,q_state &a){ // {{{
    int l,N;
    double x;double T1=2.0; // de acuerdo a valores tipicos usados con dima
    complex *phases;
    N=a.N();
int N2=N/2;
    phases=new complex[N+5];
    for(l=0;l<N;l++){
      //x=(-(E[l]+beta*absquad(a(l)))*T);
	x=-M_PI*(double)l*(double)l/((double)N);
	if(T>=0.){
		phases[l]=complex(cos(x),sin(x));
	}else{
		phases[l]=complex(cos(x),-sin(x));
	}
    }
    a.phase_mult(phases);
    delete[] phases;
} // }}}
void U_x_std(double k, q_state &a){ //{{{
int l,N;
double x,k1;
double theta,fase;
complex *phases;
k1=1.1;
load_seed();
double del=0.9;
double alpha1=1.;
double alpha2=0.;
fase=0.0; // *rand_double();
save_seed();
    N=a.N();
    phases=new complex[N+5];
    for(l=0;l<N;l++){
	theta=2*M_PI*l/((double)N);
	if(k>alpha_up){
	x=(-2*M_PI*N*alpha_up*(alpha1*cos(theta)+alpha2*cos(2.*theta)/2.))+
	(-2*M_PI*N*(k-alpha_up)*cos(ma*theta));
	}else{	
	x=(-2*M_PI*N*k*(alpha1*cos(theta)+alpha2*cos(2.*theta)/2.));}
	//(N*k1*cos((double)theta)+N*k*(sin(theta)-0.5*sin(2*theta)));
	//x=-2*T*theta;
	phases[l]=complex(cos(x),sin(x));
    }
    a.phase_mult(phases);
    delete[] phases;
return ;
} // }}}
void kick_std(double T,q_state &a){ // {{{
	if(T>=0.){	
		//~ printf("normal std\n");
		  U_x_std(T,a);
		  FFT(a,1);
		  U_p_std(2.0,a);
		  FFT(a,0);
	  }else{
		  //~ printf("inversa std\n");
	  	 FFT(a,1);
   		 U_p_std(T,a);
		 FFT(a,0);
		 U_x_std(T,a);
	  }

 //		printf("std\t %lf\n",T);
    } // }}}
void print_state_p(q_state &a,char *filename){ // {{{
FILE *fp;
double pp;
int NN=a.N();
double N2=(double)NN/2.;


fp=fopen(filename,"w");
for(int i=0;i<NN;i++){
  /*  
      if(i<=N2)pp=(double)i;
      if(i>N2)pp=-(NN-(double)i);
  */
    pp=(double)i;
	fprintf(fp,"%lf\t %.14e\t %.14e\t %.14e\n",pp,absquad(a(i)),real_teil(a(i)),imag_teil(a(i)));
   }	
 //fprintf(fp,"%lf\t %20.18lf \n",(double)NN,absquad(a(0))); //real_teil(a(0)),imag_teil(a(0)));
fprintf(fp,"\n");
fclose(fp);
return;
} // }}}
void inicio_estado(q_state &a, double q0, double p0){ // {{{
	if(coefsigma<0){
		if(coefsigma>-2){
		put_point_pos(a,q0);
		}else{
		put_point_pos(a,p0);
		FFT(a,1);
		}
	}else{
		put_coh(a,p0,q0);
	}
	a.normalize();
	    return;
} // }}}
// 
int main(int argc,char **argv){
  int N,i,cont;
  double t4=10000.0;
  double t6=1000000.0;
  double t8=100000000.0;
  double rel_err,eps,q0,p0,lam,
		dfiddet,pureza,sx,sy,ds;
  double dPR,dPtheta,r_in,r_out,E0;
  FILE *fp;
  FILE *fg;
  FILE *ff;
//  double gammax,gammay,Fx,Fy;
double ri,pri,PR,theta,ptheta,ftheta;
double FAm1,dFA,FIDm1;
double *En;
double *AVprob;
double *r_av;
double *AVLOG;
double *AVD_FID;
complex AVFID[20000];
double *AVdFID;
double *AVIPR;
double *AVIPR_B;
//

///parametros a entrar: px0,py0,x0,y0,aa(ancho del pot) , V0
seconds00=time(NULL);
//load_seed();	
int nqb=6;
ma=1.;mb=1.;
double tmax;
tmax=100000;
 double W=4.;
double beta=1.;
double gamma=0.01;
double Mt;
complex FA;
//OJO que nqb es la dimension aca
  // use above default values if argc=1
  if(argc<10 && argc>1) 
    fehler("USAGE (8 parametros):\n K \n deltaK \n tmax \n nqb (dimension) \n q0 \n p0 \n ma (if map is cat then it's a of pcat, if map is std then it's armonic) \n mapa (1:cat,2:har,3:std,4:saw) \n T_EST(<0 pos,>0 squeeze)");
  if(argc>=10){
    if(sscanf(argv[1],"%lf",&alpha_up)==0) 
      fehler("USAGE (8 parametros):\n K \n deltaK \n tmax \n nqb (dimension) \n q0 \n p0 \n ma (if map is cat then it's a of pcat, if map is std then it's armonic) \n mapa (1:cat,2:har,3:std,4:saw) \n T_EST(<0 pos,>0 squeeze)");
    if(sscanf(argv[2],"%lf",&Dt)==0) 
      fehler("USAGE (8 parametros):\n K \n deltaK \n tmax \n nqb (dimension) \n q0 \n p0 \n ma (if map is cat then it's a of pcat, if map is std then it's armonic) \n mapa (1:cat,2:har,3:std,4:saw) \n T_EST(<0 pos,>0 squeeze)");
    if(sscanf(argv[3],"%lf",&tmax)==0) 
      fehler("USAGE (8 parametros):\n K \n deltaK \n tmax \n nqb (dimension) \n q0 \n p0 \n ma (if map is cat then it's a of pcat, if map is std then it's armonic) \n mapa (1:cat,2:har,3:std,4:saw) \n T_EST(<0 pos,>0 squeeze)");
    if(sscanf(argv[4],"%d",&nqb)==0) 
      fehler("USAGE (8 parametros):\n K \n deltaK \n tmax \n nqb (dimension) \n q0 \n p0 \n ma (if map is cat then it's a of pcat, if map is std then it's armonic) \n mapa (1:cat,2:har,3:std,4:saw) \n T_EST(<0 pos,>0 squeeze)");
    if(sscanf(argv[5],"%lf",&q0)==0) 
      fehler("USAGE (8 parametros):\n K \n deltaK \n tmax \n nqb (dimension) \n q0 \n p0 \n ma (if map is cat then it's a of pcat, if map is std then it's armonic) \n mapa (1:cat,2:har,3:std,4:saw) \n T_EST(<0 pos,>0 squeeze)");
    if(sscanf(argv[6],"%lf",&p0)==0) 
      fehler("USAGE (8 parametros):\n K \n deltaK \n tmax \n nqb (dimension) \n q0 \n p0 \n ma (if map is cat then it's a of pcat, if map is std then it's armonic) \n mapa (1:cat,2:har,3:std,4:saw) \n T_EST(<0 pos,>0 squeeze)");
    if(sscanf(argv[7],"%lf",&ma)==0) 
      fehler("USAGE (8 parametros):\n K \n deltaK \n tmax \n nqb (dimension) \n q0 \n p0 \n ma (if map is cat then it's a of pcat, if map is std then it's armonic) \n mapa (1:cat,2:har,3:std,4:saw) \n T_EST(<0 pos,>0 squeeze)");
	if(sscanf(argv[8],"%d",&mapa)==0) 
      fehler("USAGE (8 parametros):\n K \n deltaK \n tmax \n nqb (dimension) \n q0 \n p0 \n ma (if map is cat then it's a of pcat, if map is std then it's armonic) \n mapa (1:cat,2:har,3:std,4:saw) \n T_EST(<0 pos,>0 squeeze)");
  	if(sscanf(argv[9],"%lf",&coefsigma)==0) 
      fehler("USAGE (8 parametros):\n K \n deltaK \n tmax \n nqb (dimension) \n q0 \n p0 \n ma (if map is cat then it's a of pcat, if map is std then it's armonic) \n mapa (1:cat,2:har,3:std,4:saw) \n T_EST(<0 pos,>0 squeeze)");

  }
  mb=ma;
//  alpha_up=0.5;
  q_state a(nqb);

printf("ma=%lf\n",ma);
lam=log((2+ma*ma+sqrt(ma*ma*(4+ma*ma)))/2.);

q_state b(nqb);
q_state c(nqb);

duracion=clock();
int dimq=a.N();
//Dt=Dt*512/(double)a.N();
double kout=alpha_up;
// if(mapa==3){
alpha_up=alpha_up/(4.*M_PI*M_PI);
double Dtout=Dt;
Dt=Dt/(2.*M_PI*dimq); //

//Dt=Dt/(2.*M_PI*dimq);
// fftw plan initialize ---------------------------------------- {{{
in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*dimq);
out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*dimq);
 planf=fftw_plan_dft_1d(dimq,in,in,FFTW_FORWARD,FFTW_ESTIMATE);
 planb=fftw_plan_dft_1d(dimq,out,out,FFTW_BACKWARD,FFTW_ESTIMATE);
//-------------------------------------------------------------- }}}
inicio_estado(a,q0,p0);
print_state_p(a,"state_ini.txt");
b=a;
for(int t=0;t<=(int)tmax;t++){
	kick_std(alpha_up,a);
	kick_std(alpha_up+Dt,b);		
} //end t loop
print_state_p(a,"state_a.txt");
print_state_p(b,"state_b.txt");
return 0;
}


