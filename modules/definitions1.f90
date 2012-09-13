module definitions
!  use math
  IMPLICIT NONE
  integer :: q
  integer, parameter :: type_Random=1,type_Bell_Random=3
  integer, parameter :: type_Random_Random=4,type_Base=2, type_file=5, type_Desentangled=6
  integer, parameter :: type_Desentangled_random=7,type_superpos_random=8
  type system_car_scalar
     integer :: qbits
     real(kind(1.d0)) :: J
     real(kind(1.d0)) :: h(3)
  end type system_car_scalar
  type system_car_Rudi
     integer :: qbits
     real(kind(1.d0)) :: J
     real(kind(1.d0))  :: h1
     real(kind(1.d0))  :: h2
  end type system_car_Rudi
  type system_car
     integer :: qbits
     real(kind(1.d0)), pointer :: J(:)
     real(kind(1.d0)) :: h(3)
  end type system_car
  type system_carvJ
     integer :: qbits
     real(kind(1.d0)), pointer :: J(:)
     real(kind(1.d0)) :: h(3)
  end type system_carvJ
  type system_carvJvB
     integer :: qbits
     real(kind(1.d0)), pointer :: J(:)
     real(kind(1.d0)), pointer :: h(:,:)
  end type system_carvJvB
  type system_carvJvBvM
     integer :: qbits
     real(kind(1.d0)), pointer :: J(:,:)
     real(kind(1.d0)), pointer :: h(:,:,:)
  end type system_carvJvBvM
  interface showvaluesparameters
     module procedure showvaluesparameters1
     module procedure showvaluesparameters2
  end interface
contains

  subroutine showvaluesparameters1(sc)
    type(system_carvJ), intent(in) :: sc
    integer                      :: j1
    open(99,file="/dev/stderr")
    write(99,*)" vary J, homogeneous h"
    write(99,*)" qubits ",sc%qbits
    write(99,*)"jotas"
    do j1=0,size(sc%j)-1
       write(99,*)"  j(",j1,")=",sc%j(j1)
    end do
    write(99,*)"Magnetic Field",sc%h
    close(99)
  end subroutine showvaluesparameters1

  subroutine showvaluesparameters2(sc)
    type(system_carvJvB), intent(in) :: sc
    integer                      :: j1
    open(99,file="/dev/stderr")
    write(99,*)" vary J, vary h"
    write(99,*)" qubits ",sc%qbits
    write(99,*)"jotas"
    do j1=0,size(sc%j)-1
       write(99,*)"  j(",j1,")=",sc%j(j1)
    end do
    write(99,*)"Magnetic Field"
    do j1=0,size(sc%j)-1
       write(99,*)"  h(",j1,")=",sc%h(:,j1)
    end do
    close(99)
  end subroutine showvaluesparameters2

end module definitions
