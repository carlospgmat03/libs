Program semiagam
! Constantes y definicion de variables  {{{
  parameter(nmax=510)     ! maximum dimension of Hilbert space 500
  parameter(ndist=2*nmax)    ! maximum dimension of distributions
  real*8 pr1,pr2,tip,tiq,dpr1,dpr2,rnum1,rnum2,beta,hpr1,hpr2
  real*4 spl(1:2*nmax,1:2*nmax)
  complex*16 sigma,phi(0:nmax-1),phi2(0:nmax-1),phi3(0:nmax-1)
  complex*16 rho1(0:nmax-1,0:nmax-1),rho2(0:nmax-1,0:nmax-1)
  complex*16 rho3(0:nmax-1,0:nmax-1)
  complex*16 wigner(0:ndist-1,0:ndist-1),cmunu(0:ndist-1,0:ndist-1),temp
  complex*16 cmn(0:ndist-1,0:ndist-1)
  complex*16 hus(0:ndist-1,0:ndist-1),hus2(0:ndist-1,0:ndist-1)
  complex*16 chord(0:ndist-1,0:ndist-1)
  equivalence (wigner(0,0),hus(0,0))!,chord(0,0))
  complex*16 work1(2*ndist+16),work2(2*ndist+16),work3(2*ndist+16)
  complex*16 psptrace,norm1,norm2,d1,d2
  complex*16 overlap,echo,entro2,echo1,echo2
  real*8 slope,entro1,q,p,q2,p2,pi,dospi
  real*8 entropy,epsp,epsq,entr1,entr2
  character*3 type1,type2,map,rep,noise,cene,ckmax
  character*5 col,cpr1,cepsp,cpi,cqi
  character*30 title,sub,ty
  character*25 fnombre,aname
  character*1 test2,dummy3,evtype,ans1
  character*2 stateop,var1
  character*6 reptype
  Real tr(6)
  integer direction,fsteps
  Integer,Parameter::ndim=16000
  real(8) xx(ndim),yy(ndim),xe(ndim),ye(ndim)

  pi=4.d0*datan(1.d0)
  tip=0.d0;tiq=0.d0
  ! }}}
! Read the state {{{
  n=0
  do
   read(*,*,end=10)ra, rr, xr,xi
!    read(*,*,end=10)xr,xi
   phi(n)=dcmplx(xr,xi)
   n=n+1
  end do
10 continue
! }}}
! PRocesamiento de la informacion {{{
  !call init_state(n,nr,q,p,tip,tiq,type1,phi)  
  call state2kirk(n,phi,tip,tiq,rho1,nmax,work1)
  call kirk2wig(1,n,rho1,nmax,nr,wigner,ndist,work1)
  nx=n*nr ;ny=nx
  wigner=wigner/real(nx,8)
  call zfft2d(1,nx,nx,wigner,ndist,work1)
  call zfft2d(-1,nx,nx,wigner,ndist,work1)
  wigner=wigner/(nx**2)
! }}}
! EScritura {{{
    do i=0,nx-1
      do j=0,nx-1
        write(*,*) i,j,real(wigner(i,j),8)
      end do
    end do
  ! }}}
 contains
!! subroutines {{{
! }}}
end  Program semiagam
