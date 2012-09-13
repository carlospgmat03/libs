module filehandling
  implicit none
contains
  subroutine Getreal_test_sort(archivo,spectrum,sorted)
    implicit none
    real(kind(1d0)),intent(out) :: spectrum(:)
    logical,intent(out)         :: sorted
    real (kind(1d0))            :: Last_Element_Tested
    character(len=*)            :: archivo
    integer                     :: j1,size_jodo
    size_jodo=size(spectrum);  open(unit=80,file=trim(archivo),status="old");
    read(80,*)Spectrum(1); Last_Element_Tested=Spectrum(1);sorted=.true.
    do j1=2,size_jodo; 
      read(80,*)Spectrum(j1); if (Spectrum(j1)<Last_Element_Tested)sorted=.false.; Last_Element_Tested=Spectrum(j1)
    end do; close(unit=80)
  end subroutine Getreal_test_sort

  function Get_LineNumber(archivo)
    implicit none
    integer            :: Get_LineNumber
    character          :: ch
    character(len=*)   :: archivo
    integer            :: ios
    open(unit=80,file=trim(archivo),status="old")
    Get_LineNumber=0
    do
       read(80,*,iostat=ios)ch
       if(ios.ne.0) exit
       Get_LineNumber=Get_LineNumber+1
    end do
    close(unit=80)
  end function Get_LineNumber

  subroutine GetReals(archivo,array)
    implicit none
    integer                       :: size_a
    real(kind(1d0)),intent(out)   :: array(:,:)
    character(len=*)              :: archivo
    integer                       :: j1
    size_a=size(array,1)
    open(unit=80,file=trim(archivo),status="old")
    do j1=1,size_a
       read(80,*)array(j1,:)
    end do
    close(unit=80)
  end subroutine GetReals


  subroutine GetReal_1(archivo,array)
   implicit none
    integer                       :: size_a
    real(kind(1d0)),intent(out)   :: array(:)
    character(len=*)              :: archivo
    integer                       :: j1
    size_a=Get_LineNumber(archivo)
    open(unit=80,file=trim(archivo),status="old")
    do j1=1,size_a
       read(80,*)array(j1)
    end do
    close(unit=80)
  end subroutine GetReal_1
end module filehandling


