module GUE_generator
  !  use la_precision, only: wp => dp
  !  use f95_lapack, only : la_heev
  use math, only: RandomGaussian,pi,I
  use mylinearalgebra, only: diagonal_from_vector,la_heev
  implicit none
contains
  subroutine generate_GUE_matrix(matrix)
    implicit none
    integer                         :: nsize,j1,j2
    complex(kind(1d0)), intent(out) :: matrix(:,:)
    nsize=size(matrix(1,:))
    do j1=1,nsize
       matrix(j1,j1)=real(RandomGaussian(1d0,0d0))
       do j2=j1+1,nsize
          matrix(j1,j2)=RandomGaussian(1/sqrt(2d0),0d0)
          matrix(j2,j1)=conjg(matrix(j1,j2))
       end do
    end do
  end subroutine generate_GUE_matrix  
  subroutine generate_GUE_normal(matrix)
    ! Aca los eigenvalores van de -N/Sqrt[Pi] a N/Sqrt[Pi]. En el centro de semicirculo la densidad
    ! es de Sqrt[Pi]/2
    integer                         :: nsize
    complex(kind(1d0)), intent(out) :: matrix(:,:)
    nsize=size(matrix(1,:))
    call generate_GUE_matrix(matrix)
    matrix=matrix*sqrt(nsize/(4*pi))
  end subroutine generate_GUE_normal
  subroutine generate_GUE_unfolded_normal_matrix(Hamiltonian)
    complex(kind(1d0)), intent(out)        :: Hamiltonian(:,:)
    real(kind(1d0)), allocatable           :: Spectrum(:)
    complex(kind(1d0)), allocatable        :: H_tmp(:,:),H_tmp2(:,:)
    integer                                :: sizesystem, j1
    sizesystem=size(Hamiltonian(:,1))
    allocate(Spectrum(sizesystem),H_tmp(sizesystem,sizesystem),H_tmp2(sizesystem,sizesystem))
    call generate_GUE_normal(Hamiltonian)
    call la_heev(Hamiltonian,spectrum,JOBZ='V')
    do j1=1,sizesystem
       Spectrum(j1)=approxcumulative(Spectrum(j1),sizesystem)
    end do
    do j1=1,sizesystem;H_tmp(j1,:)=Spectrum(j1)*conjg(Hamiltonian(:,j1));enddo;H_tmp2=Hamiltonian
    Hamiltonian=Matmul(H_tmp2,H_tmp)
    deallocate(Spectrum,H_tmp,H_tmp2)
  end subroutine generate_GUE_unfolded_normal_matrix
  subroutine generate_GUE_unfolded_normal_spectrum(spectrum)
    real(kind(1d0)), intent(out)           :: spectrum(:)
    complex(kind(1d0)), allocatable        :: Hamiltonian(:,:)
    integer                                :: sizesystem, j1
    sizesystem=size(spectrum)
    allocate(Hamiltonian(sizesystem,sizesystem))
    call generate_GUE_normal(Hamiltonian)
    call la_heev(Hamiltonian,spectrum)
    do j1=1,sizesystem
       Spectrum(j1)=approxcumulative(Spectrum(j1),sizesystem)
    end do
    deallocate(Hamiltonian)
  end subroutine generate_GUE_unfolded_normal_spectrum
  function approxcumulative(Energy,size_spectrum)
    ! se asume que entra un espectro con densidad de semicirculo, ordenado con minimo en -N/sqrt(pi) y maximo en N/sqrt(pi)
    ! a la salida se encuentra un espectro con mean level spacing igual a uno
    real(kind(1d0))                :: approxcumulative
    real(kind(1d0)), intent(in)    :: energy
    integer, intent(in)            :: size_spectrum
    real(kind(1d0))                :: rescaled_energy
    rescaled_energy=Energy*sqrt(pi)/size_spectrum
    if (rescaled_energy <= -1D0) approxcumulative=0d0
    if (rescaled_energy >= 1D0) approxcumulative=pi
    if ((-1d0<rescaled_energy).and.(rescaled_energy<1d0)) &
         approxcumulative=(asin(rescaled_energy)+rescaled_energy*sqrt(1-rescaled_energy**2)+pi/2)
    approxcumulative=(approxcumulative-pi/2)*size_spectrum/pi
  end function approxcumulative
end module GUE_generator
module GOE_generator
  !  use la_precision, only: wp => dp
  !  use f95_lapack, only : la_syev
  use math
  use mylinearalgebra
  implicit none
contains
  subroutine generate_GOE_matrix(matrix)
    ! Generates a matrix W with the following properties, see my thesis
    ! <W>_{ij}=0
    ! <W_{ij} W_{kl}>=\delta_{il}\delta_{jk}+\delta{ik}\delta_{jl}
    implicit none
    integer                         :: nsize,j1,j2
    real(kind(1d0)), intent(out)    :: matrix(:,:)
    nsize=size(matrix(1,:))
    if (size(matrix(:,1)).ne.nsize) then
       print*,"Error en generate_GOE_matrix, tamanos incompatibles ",size(matrix(:,1)),nsize
    end if
    do j1=1,nsize
       matrix(j1,j1)=real(RandomGaussian(sqrt(2d0),0d0),kind(1d0))
       do j2=j1+1,nsize
          matrix(j1,j2)=real(RandomGaussian(),kind(1d0))
          matrix(j2,j1)=matrix(j1,j2)
       end do
    end do
  end subroutine generate_GOE_matrix  
  subroutine generate_GOE_normal(matrix)
    ! El espectr de esta matriz tiene las siguientes cotas aproximadas
    ! minimo en -N/sqrt(pi) y maximo en N/sqrt(pi)
    integer                         :: nsize
    real(kind(1d0)), intent(out)    :: matrix(:,:)
    nsize=size(matrix(1,:))
    if (size(matrix(:,1)).ne.nsize) then
       print*,"Error en generate_GOE_matrix, tamanos incompatibles ",size(matrix(:,1)),nsize
       stop
    end if
    call generate_GOE_matrix(matrix)
    matrix=matrix*sqrt(nsize/(4*pi))
  end subroutine generate_GOE_normal
  subroutine generate_GOE_unfolded_normal_matrix(Hamiltonian)
    real(kind(1d0)), intent(out)           :: Hamiltonian(:,:)
    real(kind(1d0)), allocatable           :: Spectrum(:)
    integer                                :: sizesystem, j1
    sizesystem=size(Hamiltonian(:,1))
    allocate(Spectrum(sizesystem))
    call generate_GOE_normal(Hamiltonian)
    call LA_SYEV(Hamiltonian,spectrum,JOBZ='V')
    do j1=1,sizesystem
       Spectrum(j1)=approxcumulative(Spectrum(j1),sizesystem)
    end do
    Hamiltonian=Matmul(MatMul(Hamiltonian,diagonal_from_vector(Spectrum)),Transpose(Hamiltonian))
    deallocate(Spectrum)
  end subroutine generate_GOE_unfolded_normal_matrix
  subroutine generate_GOE_unfolded_normal_spectrum(spectrum)
    real(kind(1d0)), intent(out)           :: spectrum(:)
    real(kind(1d0)), allocatable           :: Hamiltonian(:,:)
    integer                                :: sizesystem, j1
    sizesystem=size(spectrum)
    allocate(Hamiltonian(sizesystem,sizesystem))
    call generate_GOE_normal(Hamiltonian)
    call LA_SYEV(Hamiltonian,spectrum)
    do j1=1,sizesystem
       Spectrum(j1)=approxcumulative(Spectrum(j1),sizesystem)
    end do
    deallocate(Hamiltonian)
  end subroutine generate_GOE_unfolded_normal_spectrum
  function approxcumulative(Energy,size_spectrum)
    ! se asume que entra un espectro con densidad de semicirculo, ordenado 
    ! con minimo en -N/sqrt(pi) y maximo en N/sqrt(pi)
    real(kind(1d0))                :: approxcumulative
    real(kind(1d0)), intent(in)    :: energy
    integer, intent(in)            :: size_spectrum
    real(kind(1d0))                :: rescaled_energy
    rescaled_energy=Energy*sqrt(pi)/size_spectrum
    if (rescaled_energy <= -1D0) approxcumulative=0d0
    if (rescaled_energy >= 1D0) approxcumulative=pi
    if ((-1d0<rescaled_energy).and.(rescaled_energy<1d0)) &
         approxcumulative=(asin(rescaled_energy)+rescaled_energy*sqrt(1-rescaled_energy**2)+pi/2)
    approxcumulative=(approxcumulative-pi/2)*size_spectrum/pi
  end function approxcumulative
end module GOE_generator
module GSE_generator
  !  use la_precision, only: wp => dp
  !  use f95_lapack, only : la_heev
  use GOE_generator
  use math
  !  use mylinearalgebra
  implicit none
contains
  subroutine generate_GSE_matrix(matrix)
    ! Los eigenvalores de esta matriz estan entre +- Sqrt[8 * nsize]
    implicit none
    integer                            :: nsize,nsize2
    complex(kind(1d0)), intent(out)    :: matrix(:,:)
    real(kind(1d0)), allocatable       :: Hs(:,:,:)
    nsize=size(matrix(1,:))
    if (mod(nsize,2).ne.0) then
       print*,"error en generate_GSE_matrix",nsize
       stop
    end if
    nsize2=nsize/2
    allocate(Hs(nsize2,nsize2,0:3))
    call generate_GOE_matrix(Hs(:,:,0))
    call generate_antisymmetric_matrix(Hs(:,:,1))
    call generate_antisymmetric_matrix(Hs(:,:,2))
    call generate_antisymmetric_matrix(Hs(:,:,3))
    matrix(1:nsize2,1:nsize2)=Hs(:,:,0)-I*Hs(:,:,3)
    matrix(nsize2+1:nsize,1:nsize2)=-Hs(:,:,2)-I*Hs(:,:,1)
    matrix(1:nsize2,nsize2+1:nsize)= Hs(:,:,2)-I*Hs(:,:,1)
    matrix(nsize2+1:nsize,nsize2+1:nsize)=Hs(:,:,0)+I*Hs(:,:,3)
    deallocate(Hs)
  end subroutine generate_GSE_matrix  
  
  subroutine generate_GSE_normal(matrix)
    !!CHECAR  ! El espectr de esta matriz tiene las siguientes cotas aproximadas
    ! minimo en -N/sqrt(pi) y maximo en N/sqrt(pi)
    integer                            :: nsize
    complex(kind(1d0)), intent(out)    :: matrix(:,:)
    nsize=size(matrix(1,:))
    if (size(matrix(:,1)).ne.nsize) then
       print*,"Error en generate_GSE_matrix, tamanos incompatibles ",size(matrix(:,1)),nsize
       stop
    end if
    call generate_GSE_matrix(matrix)
    matrix=matrix*sqrt(nsize/(8*pi))
  end subroutine generate_GSE_normal
  
  subroutine generate_GSE_unfolded_normal_spectrum(spectrum)
    real(kind(1d0)), intent(out)           :: spectrum(:)
    complex(kind(1d0)), allocatable        :: Hamiltonian(:,:)
    integer                                :: sizesystem, j1
    sizesystem=size(spectrum)
    allocate(Hamiltonian(sizesystem,sizesystem))
    call generate_GSE_normal(Hamiltonian)
    call LA_HEEV(Hamiltonian,spectrum)
    do j1=1,sizesystem
       Spectrum(j1)=approxcumulative(Spectrum(j1),sizesystem)
    end do
    deallocate(Hamiltonian)
  end subroutine generate_GSE_unfolded_normal_spectrum
  
  subroutine generate_antisymmetric_matrix(matrix)
    implicit none
    integer                         :: nsize,j1,j2
    real(kind(1d0)), intent(out)    :: matrix(:,:)
    nsize=size(matrix(1,:))
    do j1=1,nsize
       matrix(j1,j1)=0d0
       do j2=j1+1,nsize
          matrix(j1,j2)=real(RandomGaussian()); matrix(j2,j1)=-matrix(j1,j2)
       end do
    end do
  end subroutine generate_antisymmetric_matrix
end module GSE_generator
module circular_ensembles
  use math, only :random_permutation_generator,pi,I,arg
  implicit none
contains  
  subroutine generate_COE_spectrum(spectrum)
    use heap_sort_mod
    real(kind(1d0)),intent(out)     :: spectrum(:)
    integer                         :: n
    complex(kind(1d0)), allocatable :: matrix(:,:)
    real(kind(1d0)), allocatable    :: errors(:)
    n=size(spectrum)
    allocate(matrix(n,n),errors(n))
    call generate_COE_matrix(matrix)
    call diagonalize_unitary_symmetric(matrix,spectrum,errors)    
    call heap_sort(spectrum)
    deallocate(matrix,errors)
  end subroutine generate_COE_spectrum
  function generate_CUE_spectrum(n)
    use mylinearalgebra
    !use f95_lapack, only : la_geev
    use mkl95_lapack, only: geev            !  in program head
    use heap_sort_mod
    integer, intent(in)             :: n
    real(kind(1d0))                 :: generate_CUE_spectrum(n)
    complex(kind(1d0)), allocatable :: matrix(:,:),fullphases(:)
    real(kind(1d0)), allocatable    :: errors(:)
    integer                         :: error_messages1
    allocate(matrix(n,n),errors(n),fullphases(n))
    call generate_CUE_matrix(matrix)
    call GEEV(matrix,fullphases,INFO=error_messages1)
    !call LA_GEEV(matrix,fullphases,INFO=error_messages1)
    print*,"desviacino de uno",sum(abs(abs(fullphases)-1))
    generate_CUE_spectrum=arg(fullphases)+pi
    call heap_sort(generate_CUE_spectrum)
    deallocate(matrix,errors,fullphases)
  end function generate_CUE_spectrum
  subroutine generate_COE_matrix(matrix)
    use mkl95_blas, only: syrk            !  in program head
    complex(kind(1d0)), intent(out) :: matrix(:,:)
    complex(kind(1d0)),allocatable  :: cue_matrix(:,:)
    integer                         :: size_problem,j_1,j_2
    size_problem= size(matrix,1);
    if (size_problem.ne.size(matrix,2)) then; print*,"error generate_COE_matrix2";stop;endif
    !print*,"en la rutina generate_COE_matrix",size_problem
    allocate(cue_matrix(size_problem,size_problem))
    call generate_CUE_matrix(cue_matrix); 
    !print*,"la matriz CUE",sum(abs(cue_matrix)**2)
    matrix=0d0!call syrk(cue_matrix,matrix)
    call syrk(cue_matrix,matrix)
    do j_1=2,size_problem;do j_2=1,j_1-1;matrix(j_1,j_2)=matrix(j_2,j_1);enddo;enddo
    deallocate(cue_matrix)
  end subroutine generate_COE_matrix
  subroutine generate_CUE_matrix(matrix)
    use GUE_generator, only : generate_GUE_matrix,la_heev
    complex(kind(1d0)), intent(out) :: matrix(:,:)
    integer                         :: size_matrix,j_1
    real(kind(1d0)), allocatable    :: spectrum(:)
    complex(kind(1d0)), allocatable :: tmp_eigenvector(:)
    real(kind(1d0))                 :: random_phase
    integer, allocatable            :: rpg(:)
    size_matrix=size(matrix,1)
    !print*,"en la rutina generate_CUE_matrix"
    if (size_matrix.ne.size(matrix,2))then
       print*,"error en generate_CUE_matrix",size_matrix,size(matrix,2); stop;
    end if
    allocate(spectrum(size_matrix),rpg(size_matrix),tmp_eigenvector(size_matrix))
    call generate_GUE_matrix(matrix)       ! first get a H=GUE matrix
    call la_heev(matrix,spectrum,JOBZ='V') ! then diagonalize an keep the eigen vectors. 
    ! then multiply each eigenvector by a random phase
    do j_1=1,size_matrix; 
       call random_number(random_phase); matrix(:,j_1)=exp(-2*pi*I*random_phase)*matrix(:,j_1); 
    enddo
    rpg=random_permutation_generator(size_matrix)+1!random permute to forget about order of the eigenphases
    do j_1=1,size_matrix
       tmp_eigenvector=matrix(:,j_1);matrix(:,j_1)=matrix(:,rpg(j_1));matrix(:,rpg(j_1))=tmp_eigenvector
    end do
    deallocate(spectrum,rpg,tmp_eigenvector)
  end subroutine generate_CUE_matrix
  subroutine diagonalize_unitary_symmetric(U,phases,errors)
    use GOE_generator, only : la_syev
    use math, only : get_angulo_de_sincos_array2
    complex(kind(1d0)), intent(in)        :: U(:,:)
    real(kind(1d0)),intent(out)           :: phases(:)
    real(kind(1d0)),intent(out)           :: errors(:)
    integer                               :: ss,error_messages1,error_messages2
    real(kind(1d0)),allocatable           :: cosine_array(:),sine_array(:),tmpU(:,:)
    ss=size(U,1)
    if((ss.ne.size(U,2)).or.(ss.ne.size(phases)).or.(ss.ne.size(errors))) then
       print*,"error en diagonalize_unitary_symmetric",ss,size(U,2),size(phases),size(errors);stop
    end if
    allocate(cosine_array(ss),sine_array(ss),tmpU(ss,ss))
    tmpU=real(U); call LA_SYEV(tmpU,cosine_array,INFO=error_messages1)
    tmpU=aimag(U); call LA_SYEV(tmpU,sine_array,INFO=error_messages2)
    if((error_messages1.ne.0).or.(error_messages2.ne.0)) then
       print*,"algo malo.... muy malo...",error_messages1,error_messages2; stop
    end if
    call get_angulo_de_sincos_array2(sine_array,cosine_array,phases,errors)
    deallocate(cosine_array,tmpU,sine_array)
  end subroutine diagonalize_unitary_symmetric
end module circular_ensembles
