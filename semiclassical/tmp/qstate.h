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

#ifndef QUANTUM

#define QUANTUM

#include "complex.h"

#define EXP2(x) (1<<(x))
#define DEFAULT_NBQUBITS 10 
#define MAX_QUBITS 28
#define MAX_N EXP2(MAX_QUBITS)
#define MAX_PHASES 1000
#define STAT_ITER_NUM 2

inline int min(int x,int y) {return (x>y)? y:x;};

inline int max(int x,int y) {return (x<y)? y:x;};

inline double min(double x,double y) {return (x>y)? y:x;};

inline double max(double x,double y) {return (x<y)? y:x;}


class q_state{
    int nbqubits;  // number of qubits
    int nbstates;  // number of states = 2^nbqubits
    int nbpbits;   // number of qubits used for p values in classical case
    complex *coeffs; // complex coefficients for each state

    // initializes the global variables
    void global_initialize();
    // initializes the pointer
    void initialize(int nbq);

 public:
complex& operator()(int i){
#ifdef RANGE_CHECK
    if ((i<0) || (i>=nbstates)){
      printf("Indizies: %5d  %5d\n",i);
      fehler("Bereichsfehler!");
    }
#endif
    return coeffs[i];
  }
    // puts all coeffs to 0
    void put_to_zero();
    // --> Constructors and Destructor
    // default constructor, all qubits put to 0
    q_state(int nbq=DEFAULT_NBQUBITS);
    // constructor, with "coeffs" given by the vector "vec"
    q_state(int nbq,const complex *vec);
    // constructor, with coeffs[pos]=1 and coeffs[i]=0 if i!=pos
    q_state(int nbq,int pos);
    // copy construcor
    q_state(const q_state &state);
    // destructor
    ~q_state();
    // assignment operator
    q_state& operator=(q_state &a);
    // assignment operator
    q_state& operator+(q_state &a);
    // assignment operator
    void const_mult(complex phi);

    // --> diverse quantum-gate operations
    // reverse qubit order
    void reverse_state();
    // reverse qubit order
    // restricted for qubits i with:  nmin <= i < nmax 
    void reverse_state(int nmin,int nmax);
    // binary reverse of first and second half of the qubits
    void bin_reverse_state();
       // put coeffs to a complex vector
    void put(complex *a);
    // get coeffs from a complex vector
    void get(complex *a);
    // set value of "nbpbits"
    void put_pbits(int p);
    // diverse administrative or other tools
    // multiplication with arbitrary phase-vector
    void phase_mult(complex *phases);
    // effective difference between two qubit-states
    double max_diff(const q_state &a);
    // norm of difference of two qubit-states
    double norm_diff(const q_state &a);
    // acces to nbstates
    int N(void){ return nbstates; }
    // acces to nbqubits
    int L(void){ return nbqubits; }
    // acces to nbpbits
    int Lp(void){ return nbpbits; }
    // simple output to console
    void print(char *message);
    // file output, textfile-version
    void print_state(char *file_name);
    // file output (real part only), textfile-version
    void print_state_real(char *file_name);
    // file output (prob), textfile-version
    void print_prob(char *file_name);
    // same as above, setting the leftomost nbr qubits to 0.
    void print_prob(char *file_name,int nbr);
    // overload of previous fucntion...
// prob. suming over last nbr qubits.
// file output, textfile-version
//also prints epsilon column
    void print_prob(char *file_name,int nbr,double eps);
    // couting peaks after applying shor gate seqcuence
    // (Hadamard+Modular exp+FT)
    int peak_count(int nbr);
    // file input, textfile-version
    void read_state(char *file_name);
    // file input, textfile-version, 2nd version allowing for 
    // double numbers for L and p
    void read_state2(char *file_name);
    // file output, raw-version
    void save_state(char *file_name);
    // file input, raw-version
    void load_state(char *file_name);
    // converting from and to little Endian format
    void little_endian_convert(void);
    // norm-value
    double norm(void);
    // IPR-value
    double IPR(void);
    // normalization
    void normalize(void);
    // mean-value
    double mean(void);
    // variance
    double variance(void);
    // inverse participation ration
    double ipr(void);
    // measure nr-th qubit in a quatum state
    int measure(int nr);
    // effective volume of a state, i.e. number of sites 
    // with |psi(p)| > cut_value/sqrt(N)
    int volume(double cut_value);
    // fidelity
    double fidelity(const q_state &a);
    // fidelity ampliytude
    complex fid_amp(const q_state &a);
    
//Pauli ops
    //global phase
    void global_phase(complex phase);
    //Pauli Z operator on qubit m  OLD 
   
};

extern void print(complex);
extern double random_rectangular_old(double w);
// provides a random phase e^{ix} with -eps/2 <= x < eps/2
complex rand_phase(double eps); 
extern unsigned int seed_wert_1;
extern unsigned int seed_wert_2;

#endif /* !QUANTUM */
