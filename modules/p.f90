program test
  use circular_ensembles
  implicit none
  integer          :: n
  real(kind(1d0)), allocatable :: sp(:)
  complex(kind(1d0)), allocatable :: mp(:,:)
  real                           :: time_1,time_2
  n=4
  allocate(sp(n),mp(n,n))
  call cpu_time(time_1)
  call generate_COE_matrix(mp); 
  call  cpu_time(time_2)
  deallocate(sp,mp)
end program test
