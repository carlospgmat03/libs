module math
  implicit none
  real(kind(1d1)), parameter :: pi=3.14159265358979d0
  complex(kind(1d1)), parameter :: I=(0,1d0)
  interface randominteger
    module procedure randominteger_scalar
    module procedure randominteger_1d
    module procedure randominteger_2d
  end interface
  interface randominteger_pm
    module procedure randominteger_pm_scalar
    module procedure randominteger_pm_1d
    module procedure randominteger_pm_2d
  end interface
  interface norma
     module procedure norma_real
     module procedure norma_complex
  end interface
  interface mean
     module procedure mean_real
     module procedure mean_complex
  end interface
  interface MatrixTrace
     module procedure MatrixTrace_real
     module procedure MatrixTrace_complex
  end interface
  interface RandomGaussian
     module procedure RandomGaussian_scalar_vary
     module procedure RandomGaussian_scalar
  end interface
  interface arg
     module procedure arg_scalar
     module procedure arg_array
  end interface
  interface standard_deviation
     module procedure standard_deviation_real
     module procedure standard_deviation_complex
  end interface
  interface dagger
     module procedure dagger_real
     module procedure dagger_complex
  end interface
contains
  subroutine swap_values_complex(a,b)
    complex(kind(1d0)), intent(inout) :: a,b
    complex(kind(1d0))     :: z
    z=a;a=b;b=z
  end subroutine swap_values_complex
  function IdMat_real(size_matrix)
    integer, intent(in) :: size_matrix
    integer             :: j
    real                :: IdMat_real(size_matrix,size_matrix)
    IdMat_real=0d0;
    do j=1,size_matrix; IdMat_real(j,j)=1d0; enddo
  end function IdMat_real
  function IdMat_complex(size_matrix)
    integer, intent(in) :: size_matrix
    integer             :: j
    complex             :: IdMat_complex(size_matrix,size_matrix)
    IdMat_complex=0d0;
    do j=1,size_matrix; IdMat_complex(j,j)=1d0; enddo
  end function IdMat_complex

  function RandomInteger_pm_scalar(n)
    ! it gives a random integer between -n and n with uniform distribution
    implicit none
    integer :: RandomInteger_pm_scalar
    integer, intent(in)                :: n
    real(kind(1d0))                    :: v
    call random_number(v); v=(2d0*n+1)*v-n
    RandomInteger_pm_scalar=floor(v)
  end function RandomInteger_pm_scalar
  subroutine RandomInteger_pm_1d(n,integer_array)
    implicit none
    integer,intent(out)          :: Integer_array(:)
    integer, intent(in)          :: n
    integer                      :: j_s
    do j_s=1,size(Integer_array); Integer_Array(j_s)=RandomInteger_pm_scalar(n); enddo
  end subroutine RandomInteger_pm_1d
  subroutine RandomInteger_pm_2d(n,integer_array)
    implicit none
    integer,intent(out)          :: integer_array(:,:)
    integer, intent(in)          :: n
    integer                      :: j_s
    do j_s=1,size(Integer_array,1); call  RandomInteger_pm_1d(n,integer_array(j_s,:)); enddo
  end subroutine RandomInteger_pm_2d
  function RandomInteger_scalar(n)
    ! it gives a random integer between 0 and n-1 with uniform distribution
    implicit none
    integer :: RandomInteger_scalar
    integer, intent(in) :: n
    real(kind(1d0))                    :: v
    call random_number(v)
    RandomInteger_scalar=floor(v*n)
  end function RandomInteger_scalar
  subroutine RandomInteger_1d(n,integer_array)
    implicit none
    integer,intent(out)          :: Integer_array(:)
    integer, intent(in)          :: n
    integer                      :: j_s
    do j_s=1,size(Integer_array); Integer_Array(j_s)=RandomInteger_scalar(n); enddo
  end subroutine RandomInteger_1d
  subroutine RandomInteger_2d(n,integer_array)
    implicit none
    integer,intent(out)          :: integer_array(:,:)
    integer, intent(in)          :: n
    integer                      :: j_s
    do j_s=1,size(Integer_array,1); call  RandomInteger_1d(n,integer_array(j_s,:)); enddo
  end subroutine RandomInteger_2d
  function random_permutation(n)
    ! produces a random permutation , i.e. if n=3, produces, at some sample iteration
    !           1           3           2
    ! from http://www.cise.ufl.edu/~cop5536/hw2-05fa.html
    integer, intent(in)      :: n
    integer                  :: random_permutation(0:n-1),n_before,n_tmp,n_i
    integer                  :: rpg(0:n-1)
    do n_i=0,n-1;random_permutation(n_i)=n_i+1; enddo
    rpg=random_permutation_generator(n)
    do n_i=0,n-1
       n_before=rpg(n_i);  n_tmp=random_permutation(n_i);
       random_permutation(n_i)=random_permutation(n_before); random_permutation(n_before)=n_tmp
    end do
  end function random_permutation
  function random_permutation_generator(n)
    ! produces the sequence that performes a random permutation, see
    ! from http://www.cise.ufl.edu/~cop5536/hw2-05fa.html
    integer, intent(in)      :: n
    integer                  :: random_permutation_generator(0:n-1),n_i
    do n_i=0,n-1; random_permutation_generator(n_i)=RandomInteger(n-n_i)+n_i;  end do
  end function random_permutation_generator
  function test_multiplicity(state1,state2)
    complex(kind(1d0))            :: test_multiplicity
    complex(kind(1d0))            :: state1(:),state2(:)
    integer                       :: position(1)
    position=maxloc(abs(state2)); test_multiplicity=state1(position(1))/state2(position(1))
  end function test_multiplicity
  function test_parallel(state1,state2)
    real(kind(1d0))               :: test_parallel
    complex(kind(1d0))            :: state1(:),state2(:)
    test_parallel=abs(dot_product(state1,state2))/(norma(state2)*norma(state1))
  end function test_parallel
  function dagger_complex(matrix)
    complex(kind(1d0)),intent(in) :: matrix(:,:)
    complex(kind(1d0))            :: dagger_complex(size(matrix(1,:)),size(matrix(1,:)))
    if (size(matrix(1,:)).ne.size(matrix(:,1))) then
       print*,"algun pedo por ahi en dagger_complex",size(matrix(1,:)),size(matrix(:,1))
       stop
    end if
    dagger_complex = conjg(transpose(matrix))
  end function dagger_complex
  function dagger_real(matrix)
    real(kind(1d0)),intent(in) :: matrix(:,:)
    real(kind(1d0))            :: dagger_real(size(matrix(1,:)),size(matrix(1,:)))
    if (size(matrix(1,:)).ne.size(matrix(:,1))) then
       print*,"algun pedo por ahi en dagger_complex",size(matrix(1,:)),size(matrix(:,1))
       stop
    end if
    dagger_real = transpose(matrix)
  end function dagger_real
  function Rotation2D(theta)
    implicit none
    real(kind(1d0))             :: Rotation2D(2,2)
    real(kind(1d0)),intent(in)  :: theta
    Rotation2D(1,1)=cos(theta);    Rotation2D(2,2)=Rotation2D(1,1)
    Rotation2D(1,2)=-sin(theta);   Rotation2D(2,1)=-Rotation2D(1,2)
  end function Rotation2D
  function MatrixTrace_real(matrix)
    implicit none
    real(kind(1d0))  :: MatrixTrace_real
    real(kind(1d0)), intent(in)::  matrix(:,:)
    integer                       :: j
    if (size(matrix(1,:)).ne.size(matrix(:,1))) then
       print*,"Error en MatrixTrace R"
       stop
    end if
    MatrixTrace_real=0d0
    do j=1,size(matrix(1,:))
       MatrixTrace_real=MatrixTrace_real+matrix(j,j)
    end do
  end function MatrixTrace_real
  function MatrixTrace_complex(matrix)
    implicit none
    complex (kind(1d0))  :: MatrixTrace_complex
    complex(kind(1d0)), intent(in)::  matrix(:,:)
    integer                       :: j
    if (size(matrix(1,:)).ne.size(matrix(:,1))) then
       print*,"Error en MatrixTrace C"
       stop
    end if
    MatrixTrace_complex=0d0
    do j=1,size(matrix(1,:))
       MatrixTrace_complex=MatrixTrace_complex+matrix(j,j)
    end do
  end function MatrixTrace_complex
  function Distribution_Moment(x,k)
    implicit none
    real(kind(1d0)), intent(in)        :: x(:),k
    real(kind(1d0))                    :: x_bar
    real(kind(1d0))                    :: Distribution_Moment
    !print*,"entrnado en distribution moment"
    x_bar=mean_real(x)
    Distribution_Moment=mean_real( (x-x_bar)**k )
  end function Distribution_Moment
  subroutine angulo_de_sincos(seno, coseno, angulo, error)
    implicit none
    real (kind(1d0)), intent(in)    :: seno, coseno
    real (kind(1d0)), intent(out)   :: angulo, error
    real (kind(1d0))                :: xseno, xcoseno
    error=sqrt(seno**2+coseno**2)
    xseno=seno/error
    xcoseno=coseno/error
    error=1d0-error
    angulo=acos(coseno)
    if (seno<0d0) angulo=2*pi-angulo
  end subroutine angulo_de_sincos
  subroutine angulo_de_sincos2(seno, coseno, angulo, error)
    implicit none
    real (kind(1d0)), intent(in)    :: seno, coseno
    real (kind(1d0)), intent(out)   :: angulo, error
    real (kind(1d0))                :: xseno, xcoseno
    error=sqrt(seno**2+coseno**2)
    xseno=seno/error
    xcoseno=coseno/error
    error=1d0-error
    if (abs(xcoseno)<=abs(xseno)) then
       angulo=acos(xcoseno)
       if (seno<0d0) angulo=2*pi-angulo
    else
       angulo=asin(xseno)
       if (coseno<0)  angulo=pi-angulo
       if (angulo<0)  angulo=2*pi+angulo
    end if
    !write(*,*)error,seno,coseno
  end subroutine angulo_de_sincos2
  subroutine get_angulo_de_sincos_array(senos,cosenos,angulos,errores)
    real(kind(1d0)), intent(in)       :: senos(:), cosenos(:)
    real(kind(1d0)), intent(out)      :: angulos(:), errores(:)
    real(kind(1d0)), allocatable      :: abs_senos(:)
    real(kind(1d0))                   :: cosp,x,sx
    integer                           :: j,pos(1),s
    s=size(senos)
    if(((s.ne.size(cosenos)).or.(s.ne.size(errores))).or.(s.ne.size(angulos))) then
       print*,"algo mal en get_angulo_de_sincos_array"; stop
    end if
    allocate(abs_senos(s)); abs_senos=abs(senos)
    do j=1,s
       cosp=cosenos(j); if (cosp>1d0)  cosp=1d0; if (cosp<-1d0) cosp=-1d0
       x=acos(cosp); sx=sin(x); pos=minloc(abs(abs_senos-sx))
       call angulo_de_sincos2(senos(pos(1)), cosenos(j) , angulos(j) , errores(j)); 
    end do
    deallocate(abs_senos)
  end subroutine get_angulo_de_sincos_array
  subroutine get_angulo_de_sincos_array2(senos,cosenos,angulos,errores)
    real(kind(1d0)), intent(in)       :: senos(:), cosenos(:)
    real(kind(1d0)), intent(out)      :: angulos(:), errores(:)
    real(kind(1d0)), allocatable      :: abs_senos(:)
    logical(kind(1d0)), allocatable   :: non_used(:)
    real(kind(1d0))                   :: cosp,x,sx
    integer                           :: j,pos(1),s
    s=size(senos)
    if(((s.ne.size(cosenos)).or.(s.ne.size(errores))).or.(s.ne.size(angulos))) then
       print*,"algo mal en get_angulo_de_sincos_array"; stop
    end if
    allocate(abs_senos(s),non_used(s)); abs_senos=abs(senos); non_used=.true.
    do j=1,s
       cosp=cosenos(j); if (cosp>1d0)  cosp=1d0; if (cosp<-1d0) cosp=-1d0
       x=acos(cosp); sx=sin(x); pos=minloc(abs(abs_senos-sx),MASK=non_used)
       call angulo_de_sincos2(senos(pos(1)), cosenos(j) , angulos(j) , errores(j)); 
       non_used(pos)=.false.
    end do
    deallocate(abs_senos,non_used)
  end subroutine get_angulo_de_sincos_array2
  subroutine get_angulo_de_sincos_array3(senos,cosenos,angulos,errores_f)
    real(kind(1d0)), intent(in)  :: senos(:),cosenos(:)
    real(kind(1d0)), intent(out) :: angulos(:),errores_f(:)
    real(kind(1d0))              :: errores(size(senos))
    real(kind(1d0))              :: senos_sq(size(senos)),senos_tmp(size(senos)),cosenos_tmp(size(senos))
    integer                      :: dimtot,j,j1,positions(size(senos)),pos_cos,pos_sin
    dimtot=size(senos)
    senos_sq=senos**2
    senos_tmp=senos; cosenos_tmp=cosenos
    do j=1,dimtot
      do j1=1,dimtot+1-j;
        call locmins(cosenos_tmp(j1)**2,senos_sq(:dimtot+1-j),positions(j1),errores(j1))
      enddo
      pos_cos=minval(minloc(errores(:dimtot+1-j)))
      pos_sin=positions(pos_cos)
      call angulo_de_sincos2(senos_tmp(pos_sin), cosenos_tmp(pos_cos) , angulos(j) , errores_f(j))
      cosenos_tmp(pos_cos:dimtot+1-j)=cshift(cosenos_tmp(pos_cos:dimtot+1-j),1)
      senos_tmp(pos_sin:dimtot+1-j)=cshift(senos_tmp(pos_sin:dimtot+1-j),1)
      senos_sq(pos_sin:dimtot+1-j)=cshift(senos_sq(pos_sin:dimtot+1-j),1)
    enddo
    contains
    subroutine locmins(coseno_sq,sinarray_sq,posicion,error)
      real(kind(1d0)), intent(in)  :: coseno_sq,sinarray_sq(:)
      real(kind(1d0)), intent(out) :: error
      integer, intent(out)         :: posicion
      posicion=minval(minloc(abs(1-sqrt(coseno_sq+sinarray_sq))))
      error=abs(1-sqrt(coseno_sq+sinarray_sq(posicion)))
    end subroutine locmins
  end subroutine get_angulo_de_sincos_array3


  subroutine get_angulo_de_sincos_array2_m(senos,cosenos,angulos,errores)
    real(kind(1d0)), intent(in)       :: senos(:), cosenos(:)
    real(kind(1d0)), intent(out)      :: angulos(:), errores(:)
    call get_angulo_de_sincos_array2(senos,cosenos,angulos,errores)
    angulos=2*pi-angulos
  end subroutine get_angulo_de_sincos_array2_m
  function sort_real(x)
    real (kind(1d0)),intent(in)      :: x(:)
    real (kind(1d0))                 :: sort_real(size(x))
    integer                          :: s,j1,j2
    s=size(x)
    sort_real=x
    do j1=1,s-1
       do j2=1,s-1
          if (sort_real(j2)<sort_real(j2+1)) then 
             sort_real(j2:j2+1)=sort_real(j2+1:j2:-1)
          end if
       enddo
    enddo
  end function sort_real
  function standard_deviation_real(z)
    real(kind(1d0)), intent(in)     :: z(:)
    real(kind(1d0))                 :: standard_deviation_real
    real(kind(1d0))                 :: x3
    integer                         :: size_z
    size_z=size(z)
    x3=mean(z**2) - mean(z)**2
    if (x3<=0) then 
       standard_deviation_real=0d0
    else
       standard_deviation_real=sqrt(size_z*x3/(size_z-1d0))
    end if
  end function standard_deviation_real
  function standard_deviation_complex(z)
    complex(kind(1d0)), intent(in)  :: z(:)
    real(kind(1d0))                 :: standard_deviation_complex
    real(kind(1d0))                 :: x3
    integer                         :: size_z
    size_z=size(z)
    x3=mean(abs(z)**2) - abs(mean(z))**2
    if (x3<=0) then 
       standard_deviation_complex=0d0
    else
       standard_deviation_complex=sqrt(size_z*x3/(size_z-1d0))
    end if
  end function standard_deviation_complex
  function mean_complex(z)
    complex(kind(1d0)), intent(in)  :: z(:)
    complex(kind(1d0))              :: mean_complex
    
    mean_complex=sum(z)/size(z)
  end function mean_complex
  function mean_real(z)
    real(kind(1d0)), intent(in)  :: z(:)
    real(kind(1d0))              :: mean_real
    mean_real=sum(z)/size(z)
  end function mean_real
  subroutine control(j1)
    implicit none
    integer :: j1(:)
    open(unit=99,file="control.dat",status="unknown")
    write(99,*)j1
    close(99)
  end subroutine control
  function arg_array(z)
    implicit none
    complex(kind(1d0)), intent(in)  :: z(:)
    real(kind(1d0))                 :: arg_array(size(z))
    integer                         :: j1
    do j1=1,size(z)
       arg_array(j1)=arg_scalar(z(j1))
    end do
  end function arg_array
  function arg_scalar(z)
    implicit none
    real(kind(1d0))                 :: arg_scalar
    complex(kind(1d0)), intent(in)  :: z
    complex(kind(1d0))              :: zn
    if (abs(z)==0d0) then
       arg_scalar=0d0
    else
       zn=z/abs(z)
       arg_scalar=acos(real(zn,kind(1d0)))*sign(1d0,aimag(zn))
      !if (arg_scalar<0d0) arg_scalar=arg_scalar+2*pi
    end if
  end function arg_scalar
  function RandomGaussian_scalar_vary(A,x0)
    ! Statistical theories of spectra: fluctuations
    ! C. E. Porter and Rosenzwweig
    ! pg. 279
    implicit none
    complex(kind(1d0))                 :: RandomGaussian_scalar_vary
    real(kind(1d0)), intent(in)        :: a,x0
    real(kind(1d0))                    :: u,v
    call random_number(u); call random_number(v)
    u=u+0.00000000000000001d0
    RandomGaussian_scalar_vary= &
         (A*sqrt(-2*log(u))*cos(2*pi*v)+x0)+I*(A*sqrt(-2*log(u))*sin(2*pi*v)+x0)
  end function RandomGaussian_scalar_vary
  function RandomGaussian_scalar()
    implicit none
    complex(kind(1d0))                 :: RandomGaussian_scalar
    RandomGaussian_scalar= RandomGaussian_scalar_vary(1d0,0d0)
  end function RandomGaussian_scalar
  function norma_complex(state)
    implicit none
    real(kind(1d0)) :: norma_complex
    complex(kind(1d0)), intent(in) :: state(:)
    norma_complex=sqrt(real(dot_product(state,state),kind(1d0)))
  end function norma_complex
  function norma_real(state)
    implicit none
    real(kind(1d0)) :: norma_real
    real(kind(1d0)), intent(in) :: state(:)
    norma_real=sqrt(real(dot_product(state,state),kind(1d0)))
  end function norma_real
  function MatrixTraceSquare(matrix)
    implicit none
    complex (kind(1d0))  :: MatrixTraceSquare
    complex(kind(1d0)), intent(in)::  matrix(:,:)
    integer                       :: j
    if (size(matrix(1,:)).ne.size(matrix(:,1))) then
       print*,"Error en MatrixTrace"
       stop
    end if
    MatrixTraceSquare=0d0
    do j=1,size(matrix(1,:))
      !MatrixTraceSquare=MatrixTraceSquare+dot_product(matrix(j,:),matrix(:,j))
       MatrixTraceSquare=MatrixTraceSquare+dot_product(conjg(matrix(j,:)),matrix(:,j))
    end do
  end function MatrixTraceSquare
  function IntegerLog(number,base)
    implicit none
    integer               :: IntegerLog
    integer, intent(in)   :: number,base    
    IntegerLog=0
    IntegerLog=floor(log(real(number))/log(real(base)))
    if (base**IntegerLog.ne.number) IntegerLog=IntegerLog+1
    if (base**IntegerLog.ne.number) then
       print*,"error en IntegerLog"
       stop
    end if
    return
  end function IntegerLog
  function IntegerLogBase2(number)
    implicit none
    integer               :: IntegerLogBase2
    integer, intent(in)   :: number
    select case (number)
       case (1);         IntegerLogBase2=0
       case (2);         IntegerLogBase2=1
       case (4);         IntegerLogBase2=2
       case (8);         IntegerLogBase2=3
       case (16);        IntegerLogBase2=4
       case (32);        IntegerLogBase2=5
       case (64);        IntegerLogBase2=6
       case (128);       IntegerLogBase2=7
       case (256);       IntegerLogBase2=8
       case (512);       IntegerLogBase2=9
       case (1024);      IntegerLogBase2=10
       case (2048);      IntegerLogBase2=11
       case (4096);      IntegerLogBase2=12
       case (8192);      IntegerLogBase2=13
       case (16384);     IntegerLogBase2=14
       case (32768);     IntegerLogBase2=15
       case (65536);     IntegerLogBase2=16
       case (131072);    IntegerLogBase2=17
       case (262144);    IntegerLogBase2=18
       case (524288);    IntegerLogBase2=19
       case (1048576);   IntegerLogBase2=20
       case (2097152);   IntegerLogBase2=21
       case (4194304);   IntegerLogBase2=22
       case (8388608);   IntegerLogBase2=23
       case (16777216);  IntegerLogBase2=24
       case (33554432);  IntegerLogBase2=25
       case (67108864);  IntegerLogBase2=26
       case (134217728); IntegerLogBase2=27
       case default; print*,"error en averiguar el numero en la base apropiada"; stop
    end select
    return
  end function IntegerLogBase2
  subroutine single_integer_set_seed(seed)
    integer, intent(in)  :: seed
    integer              :: number_integers_fortran_seed,maxnumber,i
    integer, allocatable :: full_seed(:)
    call random_seed(SIZE = number_integers_fortran_seed)
    maxnumber=Huge(seed)
    allocate(full_seed(number_integers_fortran_seed))
    call random_seed(put=spread(seed,1,number_integers_fortran_seed))
    do i=1,number_integers_fortran_seed
      full_seed(i)=randominteger(maxnumber)
    enddo
    call random_seed(put=full_seed)
    deallocate(full_seed)
  end subroutine single_integer_set_seed
end module math
module bit_manipulation
  implicit none
contains
  function AddZeros(IntegerIn,Position1,Position2)
    !    IntegerIn      0  0  a4  a3  a2  a1  a0
    !    Position1=0                           X
    !    Position1=3               X
    !    AddZeros      a4  a3  a2  0  a1  a0   0
    implicit none
    integer,intent(in) :: IntegerIn,Position1,Position2
    integer            :: AddZeros
    integer            :: P1,P2
    if(Position1<Position2) then;   P2=Position1; P1=Position2
    else;                           P1=Position1; P2=Position2
    end if
    AddZeros=ishft(IntegerIn,2)
    AddZeros=ishftc(AddZeros,-1,P1+1)
    AddZeros=ishftc(AddZeros,-1,P2+1)
  end function AddZeros
  function AddOneZero(IntegerIn,Position)
    !    IntegerIn      0  0  a4  a3  a2  a1  a0
    !    Position=2                   X
    !    AddZeros       0 a4  a3  a2  0   a1  a0
    implicit none
    integer,intent(in) :: IntegerIn,Position
    integer            :: AddOneZero
    AddOneZero=ishft(IntegerIn,1)
    AddOneZero=ishftc(AddOneZero,-1,Position+1)
  end function AddOneZero
  function bit_inversion(number_in,qubits)
    integer, intent(in)   :: number_in,qubits
    integer               :: bit_inversion
    logical               :: bits(qubits)
    bits=get_bits(number_in,qubits)
    bits=bits(qubits:1:-1)
    bit_inversion=from_bits(bits)
  end function bit_inversion
  function from_bits(bits)
    implicit none
    logical, intent(in) :: bits(0:)
    integer             :: from_bits
    integer             :: j1
    from_bits=0
    do j1=0,size(bits)-1
       if (bits(j1)) from_bits=from_bits+2**j1
    end do
  end function from_bits
  function get_bits(num,length)
    implicit none
    integer, intent(in)      :: num,length
    logical                  :: get_bits(0:length-1)
    integer                  :: j1
    do j1=0,length-1
       get_bits(j1)=  btest(num,j1)
    end do
  end function get_bits
  function bits_on_one(num)
    !     This function returns the number of bits set
    !     to one on 
    implicit none
    integer   :: bits_on_one
    integer, intent(in)   :: num
    integer               :: j
    if (num<0) then
       print*,"error en funcion bits_on_one con num=",num
    end if
    bits_on_one=0
    if (num.ne.0) then 
       do j=0,floor(log(real(num)) /log(2d0))
          if (btest(num,j)) bits_on_one=bits_on_one+1
       enddo
    end if
    return
  end function bits_on_one
  subroutine exdig(number_in,L,n1out,n2out,nwhich)
    !     This routine takes numin and puts two numbers
    !     n1out and n2out which result from the digits
    !     of numin that are marked with the 1 bits
    !     of the number nwhich
    !     exambple
    !     nwhich=   0 1 0 1 0 0 1 = 41
    !     numin=    0 1 1 0 1 1 1 = 55
    !     n1out=      1   0     1 = 5
    !     n2out=    0   1   1 1   = 7
    implicit none
    integer, intent(in)  :: L,nwhich,number_in
    integer, intent(out) :: n1out,n2out
    integer              :: j,numin, L1,L2
    n1out=0
    n2out=0
    numin=number_in
    L1=bits_on_one(nwhich)
    L2=L-L1
    do j=0,L-1
       if (btest(nwhich,j)) then
          if (btest(numin,0)) n1out=ibset(n1out,0)
          n1out=ishftc(n1out,-1,L1)
       else
          if (btest(numin,0)) n2out=ibset(n2out,0)
          n2out=ishftc(n2out,-1,L2)
       endif
       numin=ishftc(numin,-1,L)
    end do
    return
  end subroutine exdig
  function merge_two_integers(a,b,nwhich)
    ! in this routine we merge two numbers. It is useful for doing the tensor
    ! product. 'a' and 'b' are the input number wheares as usual nwhich
    ! indicates the position. Example
    ! nwhich             = 0 0 1 0 1 0 0 1 = 41
    ! a                  =     1   0     1 = 5
    ! b                  = 1 0   0   1 0   = 18
    ! merge_two_integers = 1 0 1 0 0 1 0 1 = 165
    integer, intent(in)      :: a,b,nwhich
    integer                  :: merge_two_integers
    integer                  :: num_1,num_2,j
    num_1=a; num_2=b; j=0;     merge_two_integers=0
    if(bits_on_one(a).gt.bits_on_one(nwhich)) stop "ojo creo que aca hay un loop potencial"
    do while(num_1.ne.0.or.num_2.ne.0)
      if (btest(nwhich,j)) then; 
        if (btest(num_1,0)) merge_two_integers=ibset(merge_two_integers,j); 
        num_1=ishft(num_1,-1)
      else; 
        if (btest(num_2,0)) merge_two_integers=ibset(merge_two_integers,j); 
        num_2=ishft(num_2,-1); 
      endif
      j=j+1
    enddo
  end function merge_two_integers

end module bit_manipulation
