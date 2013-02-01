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

int quant_init_flag=0;
std::complex<double> *e_itheta;
std::complex<double> Im(0,1);

// initializes the global variables
void q_state::global_initialize(){
    int i;
    double phase;

    //    printf("MAIN INITIALIZATION !!!\n");
    quant_init_flag=1;
    e_itheta=new std::complex<double>[MAX_QUBITS+10];
    
    e_itheta[0]=-1.0;
    e_itheta[1]=std::complex<double>(0.1);
    e_itheta[2]=std::complex<double>(M_SQRT1_2,M_SQRT1_2);
    phase=M_PI/8.0;
    for(i=3;i<=MAX_QUBITS;i++){
	e_itheta[i]=std::complex<double>(cos(phase),sin(phase));
	phase/=2.0;
    }
}

// initializes the pointer
void q_state::initialize(int ndim){
    if(!quant_init_flag) global_initialize();
    if(ndim<0 || ndim>MAX_N){
// 	std::cerr << "Algo mal" << std::endl;
// 	abort();
// 	printf("ndim  =  %d   is a ",ndim);
// 	fehler("N too big !!");
    }
    nbqubits=(int)(log(1.*ndim)/log(2.));
    nbstates=ndim;
    nbpbits=(nbqubits+1)/2; // default value for number p-qubits
    coeffs=new std::complex<double>[nbstates+5];
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
// multiplication with an arbitrary phase-vector
void q_state::phase_mult(std::complex<double> *phases){
    int l;
    double a,b;

    for(l=0;l<nbstates;l++){
      //      coeffs[l]*=phases[l];
      a=coeffs[l].real()*phases[l].real()-coeffs[l].imag()*phases[l].imag();
      b=coeffs[l].real()*phases[l].imag()+coeffs[l].imag()*phases[l].real();
      coeffs[l]=std::complex<double>(a,b);
    }
}
    

// norm-value
double q_state::norm(void){
    double sum;
    int i;
    
    sum=0.0;
    for(i=0;i<nbstates;i++) sum+=abs(coeffs[i])*abs(coeffs[i]);
    return sqrt(sum);
}

// normalization
void q_state::normalize(void){
    double fak;
    int i;

    fak=1.0/norm();
    for(i=0;i<nbstates;i++) coeffs[i]*=fak;
}


