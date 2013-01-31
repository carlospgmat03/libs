/*
 * Copyright (C) 2003 Klaus Frahm <frahm@irsamc.ups-tlse.fr>
 * Quantware MIPS Center, Laboratoire de Physique Theorique
 * University Paul Sabatier, Toulouse III
 * 118, route de Narbonne, 31062 Toulouse Cedex 4 - FRANCE
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public Licens
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-
 *
 */

#include "qstate.h"
#include "random.h"
#include <math.h>

// #include "complex.cc"



int quant_init_flag=0;
complex *e_itheta;


// initializes the global variables
void q_state::global_initialize(){
    int i;
    double phase;

    //    printf("MAIN INITIALIZATION !!!\n");
    quant_init_flag=1;
    e_itheta=new complex[MAX_QUBITS+10];
    
    e_itheta[0]=-1.0;
    e_itheta[1]=I;
    e_itheta[2]=complex(M_SQRT1_2,M_SQRT1_2);
    phase=M_PI/8.0;
    for(i=3;i<=MAX_QUBITS;i++){
	e_itheta[i]=complex(cos(phase),sin(phase));
	phase/=2.0;
    }
}

// initializes the pointer
void q_state::initialize(int ndim){
    if(!quant_init_flag) global_initialize();
    if(ndim<0 || ndim>MAX_N){
	printf("ndim  =  %d   is a ",ndim);
	fehler("N too big !!");
    }
    nbqubits=(int)(log(1.*ndim)/log(2.));
    nbstates=ndim;
    nbpbits=(nbqubits+1)/2; // default value for number p-qubits
    coeffs=new complex[nbstates+5];
}

// puts all coeffs to 0
void q_state::put_to_zero(){
    int i;

    for(i=0;i<nbstates;i++) coeffs[i]=0;
}


// default constructor, all qubits put to 0
q_state::q_state(int ndim){
    initialize(ndim);
    put_to_zero();
    coeffs[0]=1.0;
}


// constructor, with "coeffs" given by the vector "vec"
q_state::q_state(int ndim,const complex *vec){
    int i;

    initialize(ndim);
    for(i=0;i<nbstates;i++) coeffs[i]=vec[i];
}


// constructor, with coeffs[pos]=1 and coeffs[i]=0 if i!=pos
q_state::q_state(int ndim,int pos){
    initialize(ndim);
    put_to_zero();
    coeffs[pos]=1.0;
}


// copy construcor
q_state::q_state(const q_state &state){
    int i;

    initialize(state.nbqubits);
    nbpbits=state.nbpbits;
    for(i=0;i<nbstates;i++) coeffs[i]=state.coeffs[i];
}


// destructor
q_state::~q_state(){
    delete[] coeffs;
}

// assignment operator
q_state& q_state::operator=(q_state &a){
    if(this!=&a){
	// adapt the vector size
	if(nbqubits!=a.nbqubits){
	    delete[] coeffs;
	    initialize(a.nbqubits);
	}
	nbpbits=a.nbpbits;
	for(int i=0;i<nbstates;i++)
	    coeffs[i]=a.coeffs[i];
    }
    return *this;
}
//assignment operator
q_state& q_state::operator+(q_state &a){
    if(this!=&a){
	// adapt the vector size
	if(nbqubits!=a.nbqubits){
	    delete[] coeffs;
	    initialize(a.nbqubits);
	}
	nbpbits=a.nbpbits;
	for(int i=0;i<nbstates;i++)
	    coeffs[i]=coeffs[i]+a.coeffs[i];
    }
    return *this;
}
//assignment operator
void q_state::const_mult(complex phi){
	for(int i=0;i<nbstates;i++)
	    coeffs[i]=phi*coeffs[i];
}
// reverse qubit order
void q_state::reverse_state(){
    int i,x,j,k;
    complex temp;

    for(i=0;i<nbstates;i++){
	j=0; x=i;
	for(k=0;k<nbqubits;k++){
	    j=(j<<1) | (x&1);
	    x=x>>1;
	}
	if(i<j){
	    temp=coeffs[i];
	    coeffs[i]=coeffs[j];
	    coeffs[j]=temp;
	}
    }
}

// reverse qubit order
// restricted for qubits i with:  nmin <= i < nmax 
void q_state::reverse_state(int nmin,int nmax){
    int i,x,j,k,mask;
    complex temp;

    mask=0;
    for(k=nmax;k<nbqubits;k++) mask=(mask<<1) | 1;
    for(k=nmin;k<nmax;k++) mask=(mask<<1);
    //    mask=mask<<(nmax-nmin);
    for(k=0;k<nmin;k++)	mask=(mask<<1) | 1;

    for(i=0;i<nbstates;i++){
	j=0; x=i>>nmin;
	for(k=nmin;k<nmax;k++){
	    j=(j<<1) | (x&1);
	    x=x>>1;
	}
	j=j<<nmin;
	j=j | (mask&i);
	if(i<j){
	    temp=coeffs[i];
	    coeffs[i]=coeffs[j];
	    coeffs[j]=temp;
	}
    }
}



// put coeffs to a complex vector
void q_state::put(complex *a){
    int i;

    for(i=0;i<nbstates;i++) a[i]=coeffs[i];
}

// get coeffs from a complex vector
void q_state::get(complex *a){
    int i;

    for(i=0;i<nbstates;i++) coeffs[i]=a[i];
}

// set value of "nbpbits"
void q_state::put_pbits(int p){
    if(p<0) p=0;
    if(p>nbqubits) p=nbqubits;
    nbpbits=p;
}






// multiplication with an arbitrary phase-vector
void q_state::phase_mult(complex *phases){
    int l;
    double a,b;

    for(l=0;l<nbstates;l++){
      //      coeffs[l]*=phases[l];
      a=coeffs[l].real()*phases[l].real()-coeffs[l].imag()*phases[l].imag();
      b=coeffs[l].real()*phases[l].imag()+coeffs[l].imag()*phases[l].real();
      coeffs[l]=complex(a,b);
    }
}
    
// effective difference between two qubit-states
double q_state::max_diff(const q_state &a){
    double dmax,diff;
    int l;

    if(nbqubits!=a.nbqubits) return 9.999E+99;
    dmax=0.0;
    for(l=0;l<nbstates;l++){
	diff=abs(coeffs[l]-a.coeffs[l]);
	if(diff>dmax) dmax=diff;
    }
    return dmax;
}

// norm of difference of two qubit-states
double q_state::norm_diff(const q_state &a){
    double sum;
    int i;

    sum=0.0;
    if(nbstates!=a.nbstates) return 9.999E+99;
    for(i=0;i<nbstates;i++){
	sum+=absquad(coeffs[i]-a.coeffs[i]);
    }
    return sqrt(sum);
}

// simple output
void q_state::print(char *message){
    int i;

    printf("%s\n----------------\n",message);
    for(i=0;i<nbstates;i++)
	printf("%5d  %16.8lf  +  i * %16.8lf\n",
	       i,coeffs[i].real(),coeffs[i].imag());
    printf("\n");
}

// file output, textfile-version
void q_state::print_state(char *file_name){
    FILE *fp;
    int i;

    fp=fopen(file_name,"w");
    fprintf(fp,"##%d %d\n",nbqubits,nbpbits);
    for(i=0;i<nbstates;i++){
	fprintf(fp,"%25.17le %25.17le\n",coeffs[i].real(),coeffs[i].imag());
    }
    fclose(fp);
}

// file output, textfile-version
void q_state::print_state_real(char *file_name){
    FILE *fp;
    int i;

    fp=fopen(file_name,"w");
  //  fprintf(fp,"%d %d\n",nbqubits,nbpbits);
    for(i=0;i<nbstates;i++){
	fprintf(fp,"%d %25.17le \n",i,coeffs[i].real());
    }
    fclose(fp);
}

// prob. file output, textfile-version
void q_state::print_prob(char *file_name){
    FILE *fp;
    int i;

    fp=fopen(file_name,"w");
  //  fprintf(fp,"%d %d\n",nbqubits,nbpbits);
    for(i=0;i<nbstates;i++){
	fprintf(fp,"%d %25.17le \n",i,absquad(coeffs[i]));
    }
    fclose(fp);
}

// overload of previous fucntion...
// prob. suming over last nbr qubits.
// file output, textfile-version
void q_state::print_prob(char *file_name,int nbr){
    FILE *fp;
    int i,j,num,dim,dim1;
    double suma;
num=nbqubits-nbr;
dim=EXP2(num);
dim1=EXP2(nbr);
    fp=fopen(file_name,"w");
  //  fprintf(fp,"%d %d\n",nbqubits,nbpbits);
    for(i=0;i<dim;i++){
      suma=0.;
       for(j=0;j<dim1;j++){
	suma+=real_teil(coeffs[i<<nbr|j]*conj(coeffs[i<<nbr|j]));
		}
	fprintf(fp,"%d %25.17le \n",i,suma);
    }
    fclose(fp);
}

// overload of previous fucntion...
// prob. suming over last nbr qubits.
// file output, textfile-version
//also prints epsilon column
void q_state::print_prob(char *file_name,int nbr,double eps){
    FILE *fp;
    int i,j,num,dim,dim1;
    double suma;
num=nbqubits-nbr;
dim=EXP2(num);
dim1=EXP2(nbr);
    fp=fopen(file_name,"w");
  //  fprintf(fp,"%d %d\n",nbqubits,nbpbits);
    for(i=0;i<dim;i++){
      suma=0.;
       for(j=0;j<dim1;j++){
	suma+=real_teil(coeffs[i<<nbr|j]*conj(coeffs[i<<nbr|j]));
		}
	fprintf(fp,"%9.7le %d %25.17le \n",eps,i,suma);
    }
    fclose(fp);
}
//
// couting peaks after applying shor gate seqcuence
// (Hadamard+Modular exp+FT)
int q_state::peak_count(int nbr){
    int i,j,count;
    double suma;
    double *buff;
    if(nbr>nbqubits)
	fehler("2nd register is too big!!!");
    int nz=EXP2(nbqubits-nbr);
    int nr=EXP2(nbr);
    buff= new double[nz + 5];
    for(i=0;i<nz;i++){
      suma=0.;
       for(j=0;j<nr;j++){
	suma+=real_teil(coeffs[i<<nbr|j]*conj(coeffs[i<<nbr|j]));
		}
	buff[i]=suma;
    }
    
    count=0;
    if(buff[0]>buff[1])count++;
    for(i=1;i<nz-1;i++){
	if(buff[i-1]<buff[i]&&buff[i+1]<buff[i])count++;
      }
    if(buff[nz-2]<buff[nz-1])count++;
   delete[] buff;
  return count;
}

// file input, textfile-version
void q_state::read_state(char *file_name){
    FILE *fp;
    int i,L,p;
    double x,y;

    fp=fopen(file_name,"r");
    if(fp==NULL) fehler("Inputfile not found !!");
    fscanf(fp,"%d%d",&L,&p);
    // adapt the vector size
    if(nbqubits!=L){
	delete[] coeffs;
	initialize(L);
    }
    put_pbits(p);
    for(i=0;i<nbstates;i++){
	fscanf(fp,"%lf%lf",&x,&y);
	coeffs[i]=complex(x,y);
    }
    fclose(fp);
}

// file input, textfile-version, 2nd version allowing for 
// double numbers for L and p
void q_state::read_state2(char *file_name){
    FILE *fp;
    int i,L,p;
    double x,y,LL,pp;

    fp=fopen(file_name,"r");
    if(fp==NULL) fehler("Inputfile not found !!");
    //    fscanf(fp,"%d%d",&L,&p);
    fscanf(fp,"%lf%lf",&LL,&pp);
    L=(int)(LL+0.1);
    p=(int)(pp+0.1);
    // adapt the vector size
    if(nbqubits!=L){
	delete[] coeffs;
	initialize(L);
    }
    put_pbits(p);
    for(i=0;i<nbstates;i++){
	fscanf(fp,"%lf%lf",&x,&y);
	coeffs[i]=complex(x,y);
    }
    fclose(fp);
}

// file output, raw-version
void q_state::save_state(char *file_name){
    FILE *fp;
    int i;

    fp=fopen(file_name,"w");
    fwrite(coeffs,sizeof(complex),nbstates,fp);
    fclose(fp);
}

// file input, raw-version
void q_state::load_state(char *file_name){
    FILE *fp;
    int i;

    fp=fopen(file_name,"r");
    if(fp==NULL) fehler("Inputfile for raw-loading not found !!");
    fread(coeffs,sizeof(complex),nbstates,fp);
    // test if file is in Big-Endian format
    if(isnan(norm())){
      printf("Converting to Little-Endian format.\n");
      little_endian_convert();
    }

    fclose(fp);
}

// converting from and to little Endian format
void q_state::little_endian_convert(void){
  char *t,tmp;
  int i,j;

  for(i=0;i<nbstates;i++){
    t=(char*)(&(coeffs[i].real()));
    for(j=0;j<4;j++){
      tmp=t[j];
      t[j]=t[7-j];
      t[7-j]=tmp;
    }
    t=(char*)(&(coeffs[i].imag()));
    for(j=0;j<4;j++){
      tmp=t[j];
      t[j]=t[7-j];
      t[7-j]=tmp;
    }
  }
}


// norm-value
double q_state::norm(void){
    double sum;
    int i;
    
    sum=0.0;
    for(i=0;i<nbstates;i++) sum+=absquad(coeffs[i]);
    return sqrt(sum);
}

// IPR-value
double q_state::IPR(void){
    double sum,x;
    int i;
    
    sum=0.0;
    for(i=0;i<nbstates;i++){
	x=absquad(coeffs[i]);
	sum+=x*x;
    }
    return 1.0/sum;

}

// normalization
void q_state::normalize(void){
    double fak;
    int i;

    fak=1.0/norm();
    for(i=0;i<nbstates;i++) coeffs[i]*=fak;
}

// mean-value
double q_state::mean(void){
    double mean;
    int i,N2,x;

    N2=nbstates/2;
    mean=0.0;
    for(i=0;i<nbstates;i++){
      x=i;
     //       if(x>N2) x-=nbstates;
      mean+=((double)x*absquad(coeffs[i]));
    }
    return mean;
}

// variance
double q_state::variance(void){
    double vv,mm,xx,L2,L;
    int i;

    mm=mean();
    vv=0.0;
    L=(double)nbstates;
    L2=L/2.0;
    for(i=0;i<nbstates;i++){
	xx=(double)i-mm;
	if(xx<-L2){ xx+=L;}
	if(xx>L2){xx-=L;}
	vv+=(xx*xx*absquad(coeffs[i]));
    }
    return vv;
}

// inverse participation ration
double q_state::ipr(void){
  double sum,x;
  int i;

  sum=0.0;
  for(i=0;i<nbstates;i++){
    x=absquad(coeffs[i]);
    sum+=x*x;
  }
  return 1.0/sum;
}
//
// measure nr-th qubit in a quatum state
int q_state::measure(int nr){
  int i,j;
  double prob[2];
  double x;
  // direct return with bad value if bit index out of range
  if(nr<0 || nr>=nbqubits) return -1; 
  // calculate the two probabilites
  prob[0]=0.0;
  prob[1]=0.0;
  for(i=0;i<nbstates;i++){
    j=(i>>nr)&1;
    prob[j]+=absquad(coeffs[i]);
  }
  // draw the random qubit value: 
  // "0" with probability prob[0]
  // "1" with probability prob[1]
  if((x=rand_double())<prob[0]) j=0; else j=1;
    // projection, put coeffs[i] to zero if nr-th bit != j
    // for(i=0;i<nbstates;i++){
    //   if(((i>>nr)&1)!=j){
    //     coeffs[i]=0;
    //   }
    // }
    // normalization after projection
    //normalize();
    // return measured qubit value
  return j;
}

// effective volume of a state, i.e. number of sites 
// with |psi(p)| > cut_value/sqrt(N)
int q_state::volume(double cut_value){
  int i,number;

  number=0;
  cut_value=cut_value*cut_value/(double)nbstates;
  for(i=0;i<nbstates;i++){
    if(absquad(coeffs[i])>cut_value) number++;
  }
  return number;
}
/**/
// fidelity
double q_state::fidelity(const q_state &a){
    complex sum;
    int i;

    sum=0.0;
    if(nbstates!=a.nbstates) return 0.0;
    for(i=0;i<nbstates;i++){
	sum+=conj(coeffs[i])*a.coeffs[i];
    }
    return absquad(sum);
}
///////////////////////
complex q_state::fid_amp(const q_state &a){
    complex sum;
    int i;

    sum=0.0;
    if(nbstates!=a.nbstates) return 0.0;
    for(i=0;i<nbstates;i++){
	sum+=conj(coeffs[i])*a.coeffs[i];
    }
    return sum;
}

// provides a random phase e^{ix} with -eps/2 <= x < eps/2

complex rand_phase(double eps){
    double x=random_rectangular_old(eps);
    return complex(cos(M_PI*x),sin(M_PI*x));
}
//**********************************************************************/

//*********************************************************
// Pauli gates
//********************************************************
void q_state::global_phase(complex phase){
   complex tt;
    for(int l=0;l<nbstates;l++){
	tt=coeffs[l]*phase;
	coeffs[l]=tt;
    }
}

//------------------------------
