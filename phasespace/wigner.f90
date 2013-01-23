Program semiagam
use phase_space_routines
! Constantes y definicion de variables  {{{
  parameter(nmax=510)     ! maximum dimension of Hilbert space 500
  parameter(ndist=2*nmax)    ! maximum dimension of distributions
  complex(kind(1d0)) ::  phi(0:nmax-1)
  complex(kind(1d0)), allocatable ::  phi_a(:)
  complex*16 rho1(0:nmax-1,0:nmax-1)
  complex*16 wigner(0:ndist-1,0:ndist-1)
  real(kind(1d0)), allocatable :: wigner_tst(:,:)
  complex*16 work1(2*ndist+16)

  pi=4.d0*datan(1.d0)
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
 allocate(phi_a(n), wigner_tst(2*n,2*n))
 phi_a=phi(0:n-1)
!  print*, "hola"
!  call wigner_from_state_2(phi_a,wigner_tst)
 print*, size(phi_a)
 call wigner_from_state(phi_a,wigner_tst)
!  print*, "adios"
! Procesamiento de la informacion {{{
  !call init_state(n,nr,q,p,tip,tiq,type1,phi)  
!   print*, "A la entrada de state2kirk, main", n,phi(0),rho1(0,0),nmax,work1(1)
  call state2kirk(n,phi,0d0,0d0,rho1,nmax,work1)
!   print*, "En el main  , salida de state2kik", rho1(0,0)
!   print*,n,rho1(0,0), nmax, nr, wigner(0,0), ndist, work1(1)
  call kirk2wig(1,n,rho1,nmax,nr,wigner,ndist,work1)
!   print*, "Main", wigner(0,0)
  nx=n*nr ;ny=nx
  wigner=wigner/real(nx,8)
!   print*, "En la main  , salida de kirk2wig", wigner(0,0)
  call zfft2d(1,nx,nx,wigner,ndist,work1)
  call zfft2d(-1,nx,nx,wigner,ndist,work1)
!   print*, "En la main  , salida de zfft2d", wigner(0,0), nx**2
  wigner=wigner/(nx**2)
! }}}
! EScritura {{{
        write(*,*) real(wigner(0,0),8),wigner_tst(1,1)
!         write(*,*) sum(abs(wigner-wigner_tst))
    tmp=0d0
    do i=0,nx-1
      do j=0,nx-1
         tmp=tmp+abs(real(wigner(i,j),8)-wigner_tst(i+1,j+1))
!         write(*,*) i,j,real(wigner(i,j),8),wigner_tst(i,j)
      end do
    end do
    print*, tmp
  ! }}}
 deallocate (phi_a, wigner_tst)
 contains
!! subroutines {{{
subroutine wigner_from_state_2(phi,wigner) !{{{
  implicit none
  complex(kind(1d0)), intent(in)   :: phi(:)
  complex(kind(1d0)), allocatable  :: wigner_tmp(:,:)
  real(kind(1d0)), intent(out)  :: wigner(:,:)
  complex(kind(1d0)), allocatable  :: rho1(:,:), work1(:)
  integer n, nr
  n=size(phi)
  allocate (rho1(n,n), work1(4*n+16), wigner_tmp(2*n,2*n))
  if ((size(wigner,1).ne.2*n).or.(size(wigner,2).ne.2*n)) then
    print*,"Tamanos no compatibles en la rutina wigner_from_state"
    print*,size(phi), size(wigner,1), size(wigner,2)
  endif
  wigner_tmp=0d0
  call state2kirk(n,phi,0d0,0d0,rho1,n,work1)
  call kirk2wig(1,n,rho1,n,nr,wigner_tmp,2*n,work1)
  wigner_tmp=wigner_tmp/(2*n)
  call zfft2d(1,2*n,2*n,wigner_tmp,2*n,work1)
  call zfft2d(-1,2*n,2*n,wigner_tmp,2*n,work1)
  wigner_tmp=wigner_tmp/(4*n*n)
  deallocate(rho1, work1)
  wigner=real(wigner_tmp)
  deallocate(wigner_tmp)
end subroutine wigner_from_state_2 ! }}}

! }}}
end  Program semiagam
