module mylinearalgebra
  use mkl95_precision, only: wp => dp
  use mkl95_lapack, only : la_syev=>syev,la_heev=>heev
  !use la_precision, only: wp => dp
  !use f95_lapack, only : la_syev, la_heev
  !!!! OJO LA RUTINA LA_SYEVR FALLA PARA NUMERO GRANDE DE QUBITS=8 !  use f95_lapack, only : la_syevr
  use math
  use bit_manipulation
  implicit none
  interface tensor_routine
     module procedure TensorProduct_real_routine
     module procedure TensorProduct_complex_routine
  end interface
  interface tensor
     module procedure MatrixTensorProduct_real
     module procedure MatrixTensorProduct_complex
  end interface
  interface Add_operator_some_qubits
    ! We Perform the following operation:
    ! H=H+ I_{n} \otimes V
    ! n indicates the space in which the identity is acting
    !   for example if the idendity acts on the first qubit and third qubit then
    !   n= 0 0 0 0 1 0 1=5
      ! What I want to do is the tensor product of 2 matrices corresponding to qubits.
    ! The qubits that correspond to space A (and only those) must be set on 1 in nwhich
    ! Example, Let H_A be the identity matrix for one qubit, and space A beeing
    ! the one of qubit 3. Then nwhich=(0 0 0 1 0 0 0)base2=8
    module procedure Add_operator_some_qubits_real
    module procedure Add_operator_some_qubits_complex
  end interface
  interface diagonal_from_vector
     module procedure diagonal_from_vector_real
     module procedure diagonal_from_vector_complex
  end interface
  interface set_identity
    module procedure set_identity_real
    module procedure set_identity_complex
    module procedure identity_double_complex_nk
  end interface
  interface identity
    module procedure identity_double_real
    module procedure identity_double_complex
  end interface
contains
  function Matrix_I_Exp_complex(H)
    implicit none
    complex(kind(1d0)), dimension(:,:)                       :: H
    complex(kind(1d0)), dimension(size(H(1,:)),size(H(1,:))) :: Matrix_I_Exp_complex
    complex(kind(1d0)), dimension(size(H(1,:)),size(H(1,:))) :: U
    real(kind(1d0)), dimension(size(H(1,:)))                 :: eigenvalues
    U=H;   call la_heev(U,eigenvalues,JOBZ='V')
    Matrix_I_Exp_complex=matmul(U,matmul(diagonal_from_vector(exp(I*eigenvalues)),conjg(Transpose(U))))
  end function Matrix_I_Exp_complex
  function Matrix_I_Exp_real(H)
    implicit none
    real(kind(1d0)), dimension(:,:)                          :: H
    complex(kind(1d0)), dimension(size(H(1,:)),size(H(1,:))) :: Matrix_I_Exp_real
    real(kind(1d0)), dimension(size(H(1,:)),size(H(1,:)))    :: U
    real(kind(1d0)), dimension(size(H(1,:)))                 :: eigenvalues
    U=H;   call la_syev(U,eigenvalues,JOBZ='V')
    Matrix_I_Exp_real=matmul(U,matmul(diagonal_from_vector(exp(I*eigenvalues)),Transpose(U)))
  end function Matrix_I_Exp_real
  subroutine Add_operator_some_qubits_real(H,V,nwhich)
    real(kind(1d0)),intent(in)    :: V(0:,0:)
    real(kind(1d0)),intent(inout) :: H(0:,0:)
    integer,intent(in)            :: nwhich

    integer                          :: j1,j2,jI,dim_I,position_x,position_y
    dim_I=2**bits_on_one(nwhich)
    if ((size(H,1).ne.size(H,2)).or.(size(V,1).ne.size(V,2)).or.(size(V,1)* dim_I.ne.size(H,1))) then
      print*,"pedo en Add_Hamiltonian_Last_Qubit_real",size(H,1),size(H,2),size(V,1),size(V,2),dim_I,nwhich; stop
    endif

    do j1=0,size(V,1)-1
      do j2=0,size(V,1)-1
        do jI=0,dim_I-1
          position_x=merge_two_integers(jI,j1,nwhich)
          position_y=merge_two_integers(jI,j2,nwhich)
          !print*,j1,j2,jI,position_x
          H(position_x,position_y)=H(position_x,position_y)+V(j1,j2)
        enddo
      enddo
    enddo
  end subroutine Add_operator_some_qubits_real
  subroutine Add_operator_some_qubits_complex(H,V,nwhich)
    complex(kind(1d0)),intent(in)    :: V(0:,0:)
    complex(kind(1d0)),intent(inout) :: H(0:,0:)
    integer,intent(in)               :: nwhich

    integer                          :: j1,j2,jI,dim_I,position_x,position_y
    dim_I=2**bits_on_one(nwhich)
    if ((size(H,1).ne.size(H,2)).or.(size(V,1).ne.size(V,2)).or.(size(V,1)* dim_I.ne.size(H,1))) then
      print*,"pedo en Add_Hamiltonian_Last_Qubit_complex",size(H,1),size(H,2),size(V,1),size(V,2),dim_I,nwhich; stop
    endif

    do j1=0,size(V,1)-1
      do j2=0,size(V,1)-1
        do jI=0,dim_I-1
          position_x=merge_two_integers(jI,j1,nwhich)
          position_y=merge_two_integers(jI,j2,nwhich)
          !print*,j1,j2,jI,position_x
          H(position_x,position_y)=H(position_x,position_y)+V(j1,j2)
        enddo
      enddo
    enddo
  end subroutine Add_operator_some_qubits_complex
  function MatrixTensorProduct_real(H_a,H_b,nwhich)
    ! What I want to do is the tensor product of 2 matrices corresponding to qubits.
    ! The qubits that correspond to space A (and only those) must be set on 1 in nwhich
    ! Example, Let H_A be the identity matrix for one qubit, and space A beeing
    ! the one of qubit 3. Then nwhich=(0 0 0 1 0 0 0)base2=8
    integer, intent(in)         :: nwhich
    real(kind(1d0)), intent(in) :: H_a(0:,0:),H_b(0:,0:)
    real(kind(1d0))             :: MatrixTensorProduct_real(0:size(H_a(1,:))*size(H_b(1,:))-1,0:size(H_a(1,:))*size(H_b(1,:))-1)
    integer                     :: total_size, j1, j2, j1a, j1b, j2a, j2b, qubits
    total_size=size(H_a(0,:))*size(H_b(0,:))
    qubits=IntegerLog(total_size,2)
    do j1=0,total_size-1
       do j2=0,total_size-1
          call exdig(j1,qubits,j1a,j1b,nwhich)
          call exdig(j2,qubits,j2a,j2b,nwhich)
          MatrixTensorProduct_real(j1,j2)=H_a(j1a,j2a)*H_b(j1b,j2b)
       enddo
    enddo
  end function MatrixTensorProduct_real
  subroutine TensorProduct_real_routine(H_a,H_b,nwhich,H_out)
    integer, intent(in)            :: nwhich
    real(kind(1d0)), intent(in) :: H_a(0:,0:),H_b(0:,0:)
    real(kind(1d0)), intent(out):: H_out(0:,0:)
    !MatrixTensorProduct_complex(0:size(H_a(1,:))*size(H_b(1,:))-1,0:size(H_a(1,:))*size(H_b(1,:))-1)
    integer                        :: total_size, j1, j2, j1a, j1b, j2a, j2b, qubits
    total_size=size(H_a,1)*size(H_b,1)
    if((size(H_a,1).ne.size(H_a,2)).or.(Size(H_b,1).ne.size(H_b,2))&
      .or.(total_size.ne.size(H_out,1)).or.(total_size.ne.size(H_out,2))) then
      print*,"pedo en TensorProduct_real_routine",size(H_a,1),size(H_b,1),size(H_a,2),size(H_b,2),size(H_out,1),size(H_out,2)
      stop
    end if
    qubits=IntegerLogBase2(total_size)
    do j1=0,total_size-1
       do j2=0,total_size-1
          call exdig(j1,qubits,j1a,j1b,nwhich)
          call exdig(j2,qubits,j2a,j2b,nwhich)
          H_out(j1,j2)=H_a(j1a,j2a)*H_b(j1b,j2b)
       enddo
    enddo
  end subroutine TensorProduct_real_routine
  subroutine TensorProduct_complex_routine(H_a,H_b,nwhich,H_out)
    integer, intent(in)            :: nwhich
    complex(kind(1d0)), intent(in) :: H_a(0:,0:),H_b(0:,0:)
    complex(kind(1d0)), intent(out):: H_out(0:,0:)
    !MatrixTensorProduct_complex(0:size(H_a(1,:))*size(H_b(1,:))-1,0:size(H_a(1,:))*size(H_b(1,:))-1)
    integer                        :: total_size, j1, j2, j1a, j1b, j2a, j2b, qubits
    print*,"hola tpc"
    total_size=size(H_a,1)*size(H_b,1)
    if((size(H_a,1).ne.size(H_a,2)).or.(size(H_b,1).ne.size(H_b,2))&
      .or.(total_size.ne.size(H_out,1)).or.(total_size.ne.size(H_out,2))) then
      print*,"pedo en TensorProduct_complex_routine",size(H_a,1),size(H_b,1),size(H_a,2),size(H_b,2),size(H_out,1),size(H_out,2)
      stop
    end if
    qubits=IntegerLogBase2(total_size)
    do j1=0,total_size-1
       do j2=0,total_size-1
          call exdig(j1,qubits,j1a,j1b,nwhich)
          call exdig(j2,qubits,j2a,j2b,nwhich)
          H_out(j1,j2)=H_a(j1a,j2a)*H_b(j1b,j2b)
       enddo
    enddo
  end subroutine TensorProduct_complex_routine
  function MatrixTensorProduct_complex(H_a,H_b,nwhich)
    integer, intent(in)         :: nwhich
    complex(kind(1d0)), intent(in) :: H_a(0:,0:),H_b(0:,0:)
    complex(kind(1d0))             :: MatrixTensorProduct_complex(0:size(H_a(1,:))*size(H_b(1,:))-1,&
      0:size(H_a(1,:))*size(H_b(1,:))-1)
    integer                     :: total_size, j1, j2, j1a, j1b, j2a, j2b, qubits
    !print*,"Hola puto, que estres"
    total_size=size(H_a(0,:))*size(H_b(0,:))
    qubits=IntegerLog(total_size,2)
    do j1=0,total_size-1
       do j2=0,total_size-1
          call exdig(j1,qubits,j1a,j1b,nwhich)
          call exdig(j2,qubits,j2a,j2b,nwhich)
          !print*,j1,j2,j1a,j1b,j2a,j2b
          MatrixTensorProduct_complex(j1,j2)=H_a(j1a,j2a)*H_b(j1b,j2b)
       enddo
    enddo
  end function MatrixTensorProduct_complex
  function identity_double_complex_nk(size_matrix)
    integer, intent(in)           :: size_matrix
    integer                       :: j
    complex(kind(1d0))            :: identity_double_complex_nk(size_matrix,size_matrix)
    identity_double_complex_nk=0d0
    do j=1,size_matrix
       identity_double_complex_nk(j,j)=1d0
    end do
  end function identity_double_complex_nk
  function identity_double_complex(size_matrix,number_kind)
    integer, intent(in)           :: size_matrix
    complex(kind(1d0))            :: number_kind
    integer                       :: j
    complex(kind(1d0))            :: identity_double_complex(size_matrix,size_matrix)
    identity_double_complex=0d0
    if (number_kind==1d0) then
    endif
    do j=1,size_matrix
       identity_double_complex(j,j)=1d0
    end do
  end function identity_double_complex
  function identity_double_real(size_matrix)
    integer, intent(in)           :: size_matrix
    integer                       :: j
    real(kind(1d0))               :: identity_double_real(size_matrix,size_matrix)
    identity_double_real=0d0
    do j=1,size_matrix
       identity_double_real(j,j)=1d0
    end do
  end function identity_double_real
  subroutine set_identity_real(matrix)
    real(kind(1d0)), intent(out)    :: matrix(:,:)
    integer                         :: size_matrix,j
    matrix=0d0; do j=1,size_matrix;matrix(j,j)=1d0; end do
  end subroutine set_identity_real
  subroutine set_identity_complex(matrix)
    complex(kind(1d0)), intent(out) :: matrix(:,:)
    integer                         :: size_matrix,j
    matrix=0d0; do j=1,size_matrix;matrix(j,j)=1d0; end do
  end subroutine set_identity_complex
  function diagonal_from_vector_real(vector)
    real(kind(1d0))           :: vector(:)
    real(kind(1d0))           :: diagonal_from_vector_real(size(vector),size(vector))
    integer                   :: j1
    diagonal_from_vector_real=0d0
    do j1=1,size(vector)
       diagonal_from_vector_real(j1,j1)=vector(j1)
    end do
  end function diagonal_from_vector_real
  function diagonal_from_vector_complex(vector)
    complex(kind(1d0))        :: vector(:)
    complex(kind(1d0))        :: diagonal_from_vector_complex(size(vector),size(vector))
    integer                   :: j1
    diagonal_from_vector_complex=0d0
    do j1=1,size(vector)
       diagonal_from_vector_complex(j1,j1)=vector(j1)
    end do
  end function diagonal_from_vector_complex
end module mylinearalgebra
