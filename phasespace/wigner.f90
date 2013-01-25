Program semiagam
use phase_space_routines
! Constantes y definicion de variables  {{{
  parameter(nmax=510)     ! maximum dimension of Hilbert space 500
  complex(kind(1d0)) ::  phi(0:nmax-1)
  complex(kind(1d0)), allocatable ::  phi_a(:)
  real(kind(1d0)), allocatable :: wigner(:,:)
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
 allocate(phi_a(n), wigner(2*n,2*n))
 phi_a=phi(0:n-1)
 call wigner_from_state(phi_a,wigner)
! test {{{
    do i=1,2*n
      do j=1,2*n
         print*, wigner(i,j)
!          tmp=tmp+abs(real(wigner(i,j),8)-wigner_tst(i+1,j+1))
      end do
    end do
  ! }}}
 deallocate (phi_a, wigner)
end  Program semiagam
