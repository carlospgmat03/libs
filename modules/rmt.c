#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"
#include "rmt.h"

/***************************************************************/
/* Generators of GUE, GOE, CUE, HOE                            */
/* March, 2001                                                 */
/***************************************************************/


/* Gaussian random number from NR */
/* uses ran1() from NR */
double gasdev(long *idum)
{
	extern double ran1(long *idum);
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*ran1(idum)-1.0;
			v2=2.0*ran1(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


void GOE(int n, double **o, double sigma, long *iseed) {
  /* real symmetric Gaussian random matrix */
  /* P(H) \propto \exp{-Tr(H^2)/2sigma^2}. */
  /* output in o[0..n-1][0..n-1]. */
  int i,j;
  double tmp;

  /* upper triangle;sqrt(2) is because of lower el.*/
  tmp=sigma/sqrt(2);
  for(j=1;j<n;j++) for(i=0;i<j;i++) o[i][j]=gasdev(iseed)*tmp;
  /* diagonal */
  for(i=0;i<n;i++) o[i][i]=gasdev(iseed)*sigma;

  /* fill dependent elements - lower triangle */
  for(i=1;i<n;i++) for(j=0;j<i;j++) o[i][j]=o[j][i];

  return ;
}


void GUE(int n, double **h, double sigma, long *iseed, enum trace traceless) {
  /* hermitean Gaussian random matrix. */
  /* P(H) \propto \exp{-Tr(H^2)/2/sigma^2}. */
  /* <|H_ij|^2>=sigma^2 */
  /* output in h[0...n-1][0..2(n-1)+1]. */ 
  int i,j;
  double tmp;

  /* diagonal - real */
  tmp=sigma/sqrt(2);
  for(i=0;i<n;i++) {
    h[i][2*i]=gasdev(iseed)*sigma;
    h[i][2*i+1]=0; 
  }
  /* upper triangle - complex */
  for(j=1;j<n;j++) for(i=0;i<j;i++) {
    h[i][2*j]=gasdev(iseed)*tmp;
    h[i][2*j+1]=gasdev(iseed)*tmp;
  }

  /* fill dependent elements - lower triangle */
  for(i=1;i<n;i++) for(j=0;j<i;j++) {
    h[i][2*j]=h[j][2*i];
    h[i][2*j+1]=-h[j][2*i+1];
  }
  if (traceless==TRACELESS) { /* substract TrH /n */
    tmp=0;
    for(i=0;i<n;i++) tmp+=h[i][2*i];
    tmp=tmp/n;
    for(i=0;i<n;i++) h[i][2*i]-=tmp;
  }

  return ;
}


void HOE(int n, double **o, long *iseed)
     /* RANDOM ORTOGONAL matrix - from O(N), */
     /* by eigenvectors of a symmetric real matrix-GOE */
     /* that is in turn obtained by P(H)\propto\exp{-Tr(H^2)/2}. */
{
  int i,j,t1,t2;
  double *d,*e,tmp;
  
  d=(double *)malloc(sizeof(double)*n);
  e=(double *)malloc(sizeof(double)*n);
  
  /* get GOE into o[][] */
  GOE(n,o,1,iseed);

  /* shift indices of o[][] */
  for(i=0;i<n;i++) o[i]--;
  
  /* eigenvalues of real symmetric matrix in a */
  tred2(o-1,n,d-1,e-1,1);
  
  tqli(d-1,e-1,n,o-1,1);
  /* eigenvalues are in d[], eigenvec. in a[][] */

  /* free memory - we don't need eigenvalues */
  free(e);
  free(d);

  /* shift o[] back */
  for(i=0;i<n;i++) o[i]++;

  return;

  /* finaly, just for sure, permute eigenvectors - eig. is in column */
  /* n/2 of transpositions */
  for(i=0;i<n/2;i++) {
    /* two random indices */
    t1=rand() % n;
    t2=rand() % n;
    /* switch */
    for(j=0;j<n;j++) {
      tmp=o[j][t1];
      o[j][t1]=o[j][t2];
      o[j][t2]=tmp;
    }
  }
}


void CUE(int n, double **o, long *iseed)
     /* RANDOM UNITARY matrix - from U(N), */
     /* by eigenvectors of a compelex hermitean matrix-GUE */
     /* that is in turn obtained by P(H)\propto\exp{-Tr(H^2)/2}. */
{
  int i,j,t1,t2;
  double *re,*im,*wr,*fv1,*fv2,*fm1,*xr,*xi;
  double co,si,angle,tmp;
  int matz,ierr;

  /* Horrible memory management (uses 3x) !! */
  /* GUE matrix */
  re=(double *)malloc(sizeof(double)*n*n);
  im=(double *)malloc(sizeof(double)*n*n);
  /* eigenvectors */
  xr=(double *)malloc(sizeof(double)*n*n);
  xi=(double *)malloc(sizeof(double)*n*n);
  /* eigenvalues */
  wr=(double *)malloc(sizeof(double)*n);
  
  /* get the GUE matrix */
  GUE(n,o,1,iseed,TRACE);
  
  /* rewrite matrix in FORTRAN - column first order */
  for(i=0;i<n;i++) for(j=0;j<n;j++) {
    re[j*n+i]=o[i][2*j];
    im[j*n+i]=o[i][2*j+1];
  }

  /* workspace */
  fv1=(double *)malloc(sizeof(double)*n);
  fv2=(double *)malloc(sizeof(double)*n);
  fm1=(double *)malloc(sizeof(double)*n*2);
  matz=1; 

  /* solve using fortran routine */
  ch_(&n,&n,re,im,wr,&matz,xr,xi,fv1,fv2,fm1,&ierr);
  if (ierr) {
    fprintf(stderr,"ch returned error code %d.\n",ierr);
    exit(1);
  }

  /* rewrite xr xi into o[][] and multiply each eigvec. with */
  /* random phase - because of some solver conventions. */

  for(j=0;j<n;j++) {
    angle=2*M_PI*((double)rand()/RAND_MAX-0.5); /* random phase */
    //    angle=0;
    co=cos(angle);
    si=sin(angle);
    for(i=0;i<n;i++) { /* j runs over diff. eigenvectors */
      o[i][2*j]=xr[i+n*j]*co-xi[i+n*j]*si;
      o[i][2*j+1]=xi[i+n*j]*co+xr[i+n*j]*si;
    }
  }

  /* finaly, just for sure, permute eigenvectors - eig. is in column */
  /* n/2 of transpositions */
  for(i=0;i<n/2;i++) {
    /* two random indices */
    t1=rand() % n;
    t2=rand() % n;
    /* switch */
    for(j=0;j<n;j++) {
      tmp=o[j][2*t1];
      o[j][2*t1]=o[j][2*t2];
      o[j][2*t2]=tmp;
      tmp=o[j][2*t1+1];
      o[j][2*t1+1]=o[j][2*t2+1];
      o[j][2*t2+1]=tmp;
    }
  }

  /* free tmp. memory */
  free(re);
  free(im);
  free(xr);
  free(xi);
  free(wr);
  free(fv1);
  free(fv2);
  free(fm1);
}


void COE(int n, double **o, long *iseed)
     /* COE matrix from CUE matrix U as U^T U */
{
  int i,j,k;
  double **tmp;
  
  tmp=(double**)malloc(sizeof(double*)*n);
  for(i=0;i<n;i++) tmp[i]=(double*)malloc(sizeof(double)*2*n);

  /* get CUE */
  CUE(n,tmp,iseed);

  /* U^T U */
  for(i=0;i<n;i++) for(j=0;j<n;j++) {
    o[i][2*j]=o[i][2*j+1]=0;
    for(k=0;k<n;k++) {
      o[i][2*j]+=tmp[k][2*i]*tmp[k][2*j]-tmp[k][2*i+1]*tmp[k][2*j+1];
      o[i][2*j+1]+=tmp[k][2*i]*tmp[k][2*j+1]+tmp[k][2*i+1]*tmp[k][2*j];
    }
  }
  
  for(i=0;i<n;i++) free(tmp[i]);
  free(tmp);
}



void tracelessU(int n, double **U, long *iseed)
     /* In U[][] returns traceless unitary matrix of size n. */
     /* U=V D V+, where D is diagonal with l,-l pairs. */
{
  double **tmp,**tmp1,*tmpp;
  int i,j;
  
  tmp=(double**)malloc(sizeof(double*)*n);
  for(i=0;i<n;i++) tmp[i]=(double*)malloc(sizeof(double)*2*n);

  tmp1=(double**)malloc(sizeof(double*)*n);
  for(i=0;i<n;i++) tmp1[i]=(double*)malloc(sizeof(double)*2*n);

  /* get n/2 random real eigenvalues */
  tmpp=(double *)malloc(sizeof(double)*n);
  for(i=0;i<n/2;i++) {
    tmpp[i]=(double)rand()/RAND_MAX;
    tmpp[n-1-i]=M_PI+tmpp[i];
  }
  /* shuffle a bit */
  for(i=0;i<n;i++) {
    j=rand() % n;
    /* switch i-th and j-th */
    if (i!=j) {
      tmpp[j]+=tmpp[i];
      tmpp[i]=tmpp[j]-tmpp[i];
      tmpp[j]-=tmpp[i];
    }
  }
  /* get random unitary */
  CUE(n,tmp,iseed);

  /* multiply by diagonal - exponential  */
  for(i=0;i<n;i++) for(j=0;j<n;j++) {
    tmp1[i][2*j]=cos(tmpp[j])*tmp[i][2*j]-sin(tmpp[j])*tmp[i][2*j+1];
    tmp1[i][2*j+1]=cos(tmpp[j])*tmp[i][2*j+1]+sin(tmpp[j])*tmp[i][2*j];
  }
  /* multiply by U[][] */
  adjungate_matrix(n,tmp);
  cmultiply(n,tmp1,tmp,U);

  for(i=0;i<n;i++) free(tmp[i]);
  free(tmp);
  for(i=0;i<n;i++) free(tmp1[i]);
  free(tmp1);
  free(tmpp);
}
