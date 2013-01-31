#ifndef COMPLEX 

#define COMPLEX

#define I complex(0,1)

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define zahl complex

inline void fehler(char *s){
  puts(s);
  exit(1);
}

class complex{
	double re,im;
public:
	complex(double r=0,double i=0){ re=r; im=i; }
	double& real(){ return re; }
	double& imag(){ return im; }
	complex conjug(){ return complex(re,-im); }
	complex times_ii(){ return complex(-im,re); }
	complex operator+(complex a){
	  return complex(re+a.re,im+a.im);
	}
	complex operator-(complex a){
	  return complex(re-a.re,im-a.im);
	}
	complex operator-(){
	  return complex(-re,-im);
	}
	complex operator*(complex a){
	  return complex(re*a.re-im*a.im,re*a.im+im*a.re);
	}
	complex operator/(complex a){
	  double r=a.re*a.re+a.im*a.im;
	  return complex((re*a.re+im*a.im)/r,(im*a.re-re*a.im)/r);
	}
	complex operator+(double a){
	  return complex(re+a,im);
	}
	complex operator-(double a){
	  return complex(re-a,im);
	}
	complex operator*(double b){
	  return complex(re*b,im*b);
	}
	complex operator/(double b){
	  return complex(re/b,im/b);
	}
	friend complex operator+(double b,complex a){
	  return complex(a.re+b,a.im);
	}
	friend complex operator-(double b,complex a){
	  return complex(b-a.re,-a.im);
	}
	friend complex operator*(double b,complex a){
	  return complex(a.re*b,a.im*b);
	}
	friend complex operator/(double b,complex a){
	  double fak=b/(a.re*a.re+a.im*a.im);
	  return complex(a.re*fak,-a.im*fak);
	}
	int operator==(complex a){
	  return ((re==a.re) && (im==a.im));
	}
	int operator!=(complex a){
	  return ((re!=a.re) || (im!=a.im));
	}
	complex& operator+=(complex a){
	  re+=a.re; im+=a.im;
	  return *this;
	}
	complex& operator-=(complex a){
	  re-=a.re; im-=a.im;
	  return *this;
	}
	complex& operator*=(complex a){
	  double tmp_re;
	  tmp_re=re*a.re-im*a.im; 
	  im=re*a.im+im*a.re;
	  re=tmp_re;
	  return *this;
	}
	complex& operator/=(complex a){
	  double tmp_re,r;

	  r=a.re*a.re+a.im*a.im;
	  tmp_re=(re*a.re+im*a.im)/r; 
	  im=(im*a.re-re*a.im)/r;
	  re=tmp_re;
	  return *this;
	}
	complex& operator+=(double a){
	  re+=a;
	  return *this;
	}
	complex& operator-=(double a){
	  re-=a;
	  return *this;
	}
	complex& operator*=(double a){
	  re*=a; im*=a;
	  return *this;
	}
	complex& operator/=(double a){
	  re/=a; im/=a;
	  return *this;
	}
        friend complex lese();
	friend complex csqrt(complex);
	friend void print(complex);
	friend double absquad(complex a){ return a.re*a.re+a.im*a.im; }
	friend double abs(complex  a){ return sqrt(absquad(a)); }
	friend double skalar_produkt(complex a,complex b){
		return a.re*b.re+a.im*b.im;
	}
};

inline complex conj(complex x){
  return x.conjug();
}

inline double real_teil(complex x){
  return x.real();
}

inline double imag_teil(complex x){
  return x.imag();
}

/*
  extern complex operator+(complex,complex);
  extern complex operator-(complex,complex);
  extern complex operator+(complex,complex);
  extern complex operator+(complex,complex);
  extern int operator==(complex,complex);
*/
extern complex csqrt(complex);
extern complex lese(void);
extern void print(complex);

#endif
