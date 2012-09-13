program test
  use GSE_generator
  implicit none
  integer          :: n
  real(kind(1d0)), allocatable :: sp(:)
  complex(kind(1d0)), allocatable :: mp(:,:)
  real                           :: time_1,time_2
  n=4
  allocate(sp(n),mp(n,n))
  call generate_COE_matrix(mp); 
  call LA_SYEV(mp,sp)
  print*,sp
  deallocate(sp,mp)
end program test

