module concurrencemod
  use bit_manipulation, only: bits_on_one
  use mkl95_precision, only: wp => dp
  use mkl95_lapack, only : geev
  implicit none
contains
  function concurrence(rho)
    implicit none
    real(kind(1d0))                  :: concurrence
    complex(kind(1d0)), intent(in)   :: rho(0:,0:)
    complex(kind(1d0))               :: rho_tilde(0:3,0:3)
    real(kind(1d0))                  :: eigen_values(0:3)
    complex(kind(1d0))               :: ceigen_values(0:3)
    integer                          :: j1,j2
    if ((size(rho,1).ne.4).or.(size(rho,2).ne.4)) then
      print*,"algun pedo en concurrence",size(rho,1),size(rho,2)
      stop
    endif
    do j1=0,3
       do j2=0,3
          rho_tilde(j1,j2)=conjg(rho(3-j1,3-j2))
          if (mod(bits_on_one(j1)+bits_on_one(j2),2)==1) rho_tilde(j1,j2)=-rho_tilde(j1,j2) 
       end do
    end do
    rho_tilde=matmul(rho,rho_tilde)
    call geev(rho_tilde,ceigen_values); eigen_values=real(ceigen_values)
    if ( maxval(eigen_values)<0) then; Concurrence=0
    else; Concurrence=2*sqrt((maxval(eigen_values) ))- sum(sqrt(abs(eigen_values)))
    endif
  end function Concurrence
end module concurrencemod

