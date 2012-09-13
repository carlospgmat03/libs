module quantumwau
  !  use definitions
  !  use EvolutionOperators
  use bit_manipulation
  use math
  implicit none
contains
  function pauli(n)
    integer, intent(in)    :: n
    complex(kind(1d0))     :: pauli(0:1,0:1)
    pauli=0.
    if (n==0) then
      pauli(0,0)=1.
      pauli(1,1)=1.
    end if 
    if (n==1) then
      pauli(0,1)=1d0
      pauli(1,0)=1d0
    end if 
    if (n==2) then
      pauli(0,1)=(0d0,-1d0)
      pauli(1,0)=(0d0, 1d0)
    end if 
    if (n==3) then
      pauli(0,0)=1d0
      pauli(1,1)=-1d0
    end if 
  end function pauli
  subroutine sigma_x_pos(state,qubit)
    ! here they do the the sigma_x^{(j)} acting on a state
    ! IF where  sigma_y is bitflip
    ! use math, only : swap_values_complex
    complex(kind(1d0)), intent(inout) :: state(0:)
    integer,intent(in)                :: qubit 
    integer                           :: j,number_1
    do j=0,size(state)/2-1
      number_1=merge_two_integers(1,j,2**qubit)
      state(number_1)= -state(number_1)
    enddo
  end subroutine sigma_x_pos
  subroutine sigma_y_pos(state,qubit)
    ! here they do the the sigma_x^{(j)} acting on a state
    ! IF where  sigma_y is bitflip
    use math, only : swap_values_complex
    use bit_manipulation, only :merge_two_integers
    complex(kind(1d0)), intent(inout) :: state(0:)
    integer,intent(in)                :: qubit 
    integer                           :: j,qdos,number_1,number_2
    qdos=size(state)
    do j=0,qdos/2-1
      number_1=merge_two_integers(0,j,2**qubit)
      number_2=ibset(number_1,qubit)
      call swap_values_complex(state(number_1),state(number_2))
    enddo
  end subroutine sigma_y_pos
  subroutine sum_sigma_x(state)
    ! here they do the the \sum_j sigma_x^{(j)} acting on a state
    ! IF where  sigma_x is diagonal
    use math, only                 :  integerlogbase2
    use bit_manipulation, only     : bits_on_one
    complex(kind(1d0)), intent(inout) :: state(0:)
    integer                           :: qubits, j
    qubits=IntegerLogBase2(size(state))
    do j=0,size(state)-1
      state(j)=state(j)*(qubits-2*bits_on_one(j))
    enddo
  end subroutine sum_sigma_x
  subroutine sum_sigma_y(state)
    ! here they do the the \sum_j sigma_x^{(j)} acting on a state
    ! IF where  sigma_y is bit_flip
    use math, only                 :  integerlogbase2
    complex(kind(1d0)), intent(inout) :: state(0:)
    complex(kind(1d0)), allocatable   :: state_0(:)
    integer                           :: qubits, qubit,j
    complex(kind(1d0))                :: z
    allocate(state_0(0:size(state)-1))
    qubits=IntegerLogBase2(size(state))
    state_0=state;
    do j=0,size(state)-1
      z=0d0
      do qubit=0,qubits-1
        z=z+state_0(ieor(j,2**qubit))
      enddo 
      state(j)=z
    enddo
    deallocate(state_0)
  end subroutine sum_sigma_y
  subroutine Merge_States(state_in_1,state_in_2,nwhich, state_out)
    ! la opcion nwich me dice donde estan los qubits que se van a corresponder al state_in_1
    ! por ejemplo si 
    ! nwichi=14 =  0 0 1 1 1 0
    ! qubit 0  =>  5 4 3 2 1 0
    ! tonces se van a poner el state_in en los qubits 1, 2, y 3
    implicit none
    integer, intent(in)                                :: nwhich
    complex(kind(1d0)), intent(in)                     :: state_in_1(0:),state_in_2(0:)
    complex(kind(1d0)), intent(out)                    :: state_out(0:)
    integer                                            :: qin1, qin2,qout,j1, n1out,n2out
    qin1=IntegerLogBase2(size(state_in_1))
    qin2=IntegerLogBase2(size(state_in_2))
    qout=IntegerLogBase2(size(state_out))
    if ((qin1+qin2.ne.qout).or. (qin1.ne.bits_on_one(nwhich))) then
       print*,"error en merge states",qin1,qin2,qout,bits_on_one(nwhich),nwhich
       stop 1
    end if
    do j1=0,2**qout-1
       call  exdig(j1,qout,n1out,n2out,nwhich)
       state_out(j1)=state_in_1(n1out)*state_in_2(n2out)
       !       print*,"in mergestates"j1,n1out,n2out,&
       !           abs(state_in_1(n1out)),abs(state_in_2(n2out)),abs(state_out(j1))
    end do
  end subroutine Merge_States
  recursive subroutine state_preparation(state,option,internal_options)
    ! state_preparation options:
    !
    ! 1. Random gaussian variables, normalized
    ! 2. Base State (internal options 1 indica la base)
    ! 3. bell times random 
    ! 4. random times random 
    ! 5. infile
    ! 6. Random desintangled state
    ! 7. (Random desintangled) times (random)
    ! 8. (|0>+|1>)/sqrt(2) times (random)
    ! 9. (|0> ***OR*** |1>) times (random)
    ! 10. Bell times (base)
    ! 11. Bell in 0 and 1 times (random) times (random)
    ! 12. non random Bell times (random)
    ! 13. random times (base)
    ! 14. random real state with gaussian coefficients centered at zero
    ! 16. GHZ=(|0...0>+ |1...1>)/\sqrt{2}
    implicit none
    complex(kind(1d0)), intent(out) :: state(0:)
    complex(kind(1d0))              :: x1,x2
    complex(kind(1d0)), allocatable :: state_aux1(:), state_aux2(:), state_aux3(:)
    integer, intent(in)             :: option, internal_options(:)
    integer                         :: j1, nwhich, nb, bath1, sbath1,qubits
    integer,allocatable             :: previous_seed(:)
    real(kind(1d0)),allocatable     :: xhi(:)
    qubits=IntegerLogBase2(size(state))
    state=0d0
    allocate(xhi(2))
    if (option==1) then
       do j1=0,size(state)-1
          state(j1)=RandomGaussian()
       enddo
       state=state/norma(state)
    else if (option==2) then
       state=0d0
       state(internal_options(1))=1d0
    else if (option==3) then
       allocate (state_aux1(0:2**(qubits-2)-1),state_aux2(0:3))
       call state_preparation(state_aux1,1,(/0/))
       call random_bell(state_aux2)
       call Merge_States(state_aux1,state_aux2,2**qubits-1-internal_options(1), state)
       deallocate(state_aux1,state_aux2)
    else if (option==4) then
      nwhich=internal_options(1)
      nb=bits_on_one(nwhich)
      allocate (state_aux2(0:2**(qubits-nb)-1),state_aux1(0:2**nb-1))
      call state_preparation(state_aux1,1,(/0/))
      call state_preparation(state_aux2,1,(/0/))
      !       print*,"hola state preparation option 4", nwhich
      call Merge_States(state_aux1,state_aux2,nwhich, state)
      deallocate(state_aux1,state_aux2)
    else if (option==5) then
!     Export["~/investigacion/vicente/state.dat", 
!                Table[{Random[], Random[]}, {2^14}]]
       open(unit=98,file="state.dat",status="old")
       do j1=0,2**qubits-1
          read(98,*)x1,x2
          state(j1)=x1+I*x2
       end do
       state=state/norma(state)
       close(98)
    else if (option==6) then
      allocate(state_aux1(0:1))
      allocate(state_aux3(0:0))
      state_aux3=1d0
      do j1=1,IntegerLogBase2(size(state(:)))
        call state_preparation(state_aux1,1,(/0/))
        !call random_number(xhi)
        !call random_number(theta)
        !state_aux1(0)=exp(-2*pi*I*(xhi(1)))*sin(2*pi*theta)
        !state_aux1(1)=exp(-2*pi*I*(xhi(2)))*cos(2*pi*theta)
        allocate(state_aux2(0:size(state_aux3)-1))
        state_aux2=state_aux3
        deallocate(state_aux3)
        allocate(state_aux3(0:2*size(state_aux2)-1))
        call Merge_States(state_aux1,state_aux2,1, state_aux3)
        deallocate(state_aux2)
      enddo
      state=state_aux3
      state=state/norma(state)
      deallocate(state_aux1,state_aux3)
    else if (option==7) then
      nwhich=internal_options(1)
      nb=bits_on_one(nwhich)
      allocate (state_aux2(0:2**(qubits-nb)-1),state_aux1(0:2**nb-1))
      call state_preparation(state_aux1,6,(/0/))
      call state_preparation(state_aux2,1,(/0/))
      !print*,"hola state preparation option 4", nwhich
      call Merge_States(state_aux1,state_aux2,nwhich, state)
      state=state/norma(state)
      deallocate(state_aux1,state_aux2)
    else if (option==8) then
       nwhich=internal_options(1)
       allocate (state_aux2(0:2**(qubits-1)-1),state_aux1(0:1))
       state_aux1=1d0/sqrt(2d0)
       call state_preparation(state_aux2,1,(/0/))
       call Merge_States(state_aux1,state_aux2,nwhich, state)
       deallocate(state_aux1,state_aux2)
    else if (option==9) then
       nwhich=internal_options(1)
       allocate (state_aux2(0:2**(qubits-1)-1),state_aux1(0:1))
       if (nwhich==0) then
          state_aux1=(/1d0,0d0/)
       else if (nwhich==1) then
          state_aux1=(/0d0,1d0/)
       else
          print*,"opcion mala en state preparation, opcion 9"
          stop
       end if
       state_aux1=state_aux1/norma(state_aux1)
       call state_preparation(state_aux2,1,(/0/))
       call Merge_States(state_aux1,state_aux2,1, state)
       deallocate(state_aux1,state_aux2)
    else if (option==10) then
       allocate (state_aux1(0:2**(qubits-2)-1),state_aux2(0:3))
       call state_preparation(state_aux1,2,(/internal_options(2)/))
       call random_bell(state_aux2)
       call Merge_States(state_aux1,state_aux2,2**qubits-1-internal_options(1), state)
       deallocate(state_aux1,state_aux2)
    else if (option==11) then
      !   Lo primero es hacer que los dos aleatorios se junten. Para esto necesito rebir como
      !   parametro quien es el bano termico inicial
      !      internal_options(1)  Este parametro me dice quienes estan en el primer bano
      ! termico. 
      ! las condicinoes que tiene qeu tener este es que sea multiplo de 4,
      ! y que el numero de bits sea pues menor que q-2 y ya!
      bath1=internal_options(1)
      sbath1=bits_on_one(bath1)
      if ((mod(bath1,4).ne.0).or.(sbath1>=qubits-2)) then
        print*,"error en state preparation, opcion 11", bath1
        stop
      end if
      allocate (state_aux1(0:2**(sbath1)-1),state_aux2(0:2**(qubits-sbath1-2)-1),state_aux3(0:2**(qubits-2)-1))
      call state_preparation(state_aux1,1,(/0/))
      call state_preparation(state_aux2,1,(/0/))
      call Merge_States(state_aux1,state_aux2,bath1/4, state_aux3)
      deallocate(state_aux1,state_aux2)
      allocate (state_aux1(0:3))
      call random_bell(state_aux1)
      call Merge_States(state_aux1,state_aux3,3, state)
      deallocate(state_aux1,state_aux3)
    else if (option==12) then
       allocate (state_aux1(0:size(state)/4-1),state_aux2(0:3))
       call state_preparation(state_aux1,1,(/0/))
       state_aux2=0d0; state_aux2(0)=1/sqrt(2d0);state_aux2(3)=1/sqrt(2d0);
       call Merge_States(state_aux1,state_aux2,size(state)-1-internal_options(1), state)
       deallocate(state_aux1,state_aux2)
    else if (option==13) then
      nwhich=internal_options(1)
      nb=bits_on_one(nwhich)
      allocate (state_aux2(0:2**(qubits-nb)-1),state_aux1(0:2**nb-1))
      call state_preparation(state_aux1,1,(/0/))
      call state_preparation(state_aux2,2,(/internal_options(2)/))
      !print*,"hola state preparation option 4", nwhich
      call Merge_States(state_aux1,state_aux2,nwhich, state)
      deallocate(state_aux1,state_aux2)
    else if (option==14) then
      call state_preparation(state,1,(/0/));       state=abs(state)
    else if (option==15) then
      call state_preparation(state,1,(/0/));state=real(state); state=state/norma(state)
    else if (option==16) then
      state(0)=1/sqrt(2d0); state(size(state)-1)=state(0)
    else if (option==17) then
      nwhich=internal_options(1)
      nb=bits_on_one(nwhich)
      allocate (state_aux2(0:2**(qubits-nb)-1),state_aux1(0:2**nb-1))
      call state_preparation(state_aux1,16,(/0/))
      call state_preparation(state_aux2,1,(/0/))
      call Merge_States(state_aux1,state_aux2,nwhich, state)
      deallocate(state_aux1,state_aux2)
    else if (option==18) then
      call random_seed(size=nb); 
      if(nb.ne.size(internal_options)) stop "state_preparation 18"
      allocate(previous_seed(nb))
      call random_seed(get=previous_seed)
      call random_seed(put=internal_options)
      call state_preparation(state,1,(/0/))
      call random_seed(put=previous_seed)
      deallocate(previous_seed)
    else if (option==19) then
      call random_seed(size=nb); 
      allocate(previous_seed(nb))
      call random_seed(get=previous_seed)
      call random_seed(put=internal_options)
      call state_preparation(state,1,(/0/))
      call random_seed(put=previous_seed)
      deallocate(previous_seed)
    else
      print*,"Bad option in state preparation"
      stop
    endif
    deallocate(xhi)
  end subroutine state_preparation
  subroutine random_bell(state)
    implicit none
    complex(kind(1d0)), intent(out)            :: state(0:3)
    real(kind(1d0))                            :: phi(2), theta(2),global_phase,theta_rel
    !    real(kind(1d0))                            :: xhi(2), x(2)
    state=0d0
    call random_number(global_phase); global_phase=2*pi*global_phase
    call random_number(phi);          phi=2*pi*phi
    call random_number(theta);        theta=2*theta-1d0; theta(1)=acos(theta(1));theta(2)=acos(theta(2))
    theta_rel=(theta(1)-theta(2))/2
    state(0)=cos(theta_rel)
    state(1)=exp(I*phi(2))*sin(theta_rel)
    state(2)=exp(I*phi(1))*sin(theta_rel)
    state(3)=-exp(I*sum(phi))*state(0)
    state=exp(I*global_phase)*state/sqrt(2d0)
    !    state=0d0; state(0)=1/sqrt(2d0);state(3)=1/sqrt(2d0);
  end subroutine random_bell
  recursive subroutine state_preparation_non_qubits(state,option,internal_options)
    ! state_preparation options:
    !
    ! 1. Random gaussian variables, normalized
    implicit none
    complex(kind(1d0)), intent(out) :: state(0:)
    integer, intent(in)             :: option, internal_options(:)
    integer                         :: j1
    if(size(internal_options)<0) then
      write(*,*)"Atencion en state_preparation_non_qubits"
    endif
    state=0d0
    if (option==1) then
       do j1=0,size(state)-1
          state(j1)=RandomGaussian()
       enddo
       state=state/norma(state)
    else
      print*,"Bad option in state preparation"
      stop
    endif
  end subroutine state_preparation_non_qubits
  function Purity(rho)
    implicit none
    real (kind(1d0))  :: Purity
    complex(kind(1d0)), intent(in)::  rho(:,:)
    Purity=MatrixTrace(matmul(rho,rho))
  end function Purity
  function partly_entangled_pure2qubit(theta)
    complex(kind(1d0))                                 :: partly_entangled_pure2qubit(0:3)
    real(kind(1d0)),intent(in)                         :: theta
    complex(kind(1d0)),dimension(0:1)                  :: p1,p2
    complex(kind(1d0)),dimension(0:3)                  :: cc,oo
    call state_preparation(p1,1,(/0/))
    call state_preparation(p2,1,(/0/))
    cc=(/p1(0)*p2(0),p1(1)*p2(0),p1(0)*p2(1),p1(1)*p2(1)/)
    p1=conjg(p1(1:0:-1)*(/1d0,-1d0/));     p2=conjg(p2(1:0:-1)*(/1d0,-1d0/))
    oo=(/p1(0)*p2(0),p1(1)*p2(0),p1(0)*p2(1),p1(1)*p2(1)/)
    partly_entangled_pure2qubit=cos(theta)*cc+sin(theta)*oo
  end function partly_entangled_pure2qubit
  subroutine random_uniform(state)
    implicit none
    complex(kind(1d0)), intent(out)            :: state(0:3)
    real(kind(1d0))                            :: xhi(0:3), theta(0:2)
    call random_number(xhi)
    call random_number(theta)
    state(0)=exp(2*pi*I*xhi(0))*sqrt(1-theta(0)**2)
    state(1)=exp(2*pi*I*xhi(1))*theta(0)*           sqrt(1-theta(1)**2)
    state(2)=exp(2*pi*I*xhi(2))*theta(0)*           theta(1)*           sqrt(1-theta(2)**2)
    state(3)=exp(2*pi*I*xhi(3))*theta(0)*           theta(1)*           theta(2)
  end subroutine random_uniform
  subroutine safe_PartialTrace(statin,rho,nwhich)
    ! DANGER:: I am leaving qubits indicated by nwhich
    ! If the system is of 7 qubits, and mwhich is
    ! 0 0 1 0 0 1 1 = 19 = nwhich
    ! 6 5 4 3 2 1 0 => qubits             
    ! i trace out qubits 2, 3, 6 and 6, leaving qubits 0, 1 and 4.
    implicit none
    complex(kind(1d0))                         ::  statin(:)
    complex(kind(1d0))                         ::  rho(0:,0:)
    integer, intent(in)                        ::  nwhich
    integer                                    ::  err(3), n_final,j1,j2,qubits
    complex(kind(1d0)), pointer                ::  staout(:,:), st1(:), st2(:)
    qubits=IntegerLogBase2(size(statin))
    n_final=bits_on_one(nwhich)
    if ((size(rho,1).ne.2**n_final).or.(size(rho,2).ne.2**n_final)) then
      print*,"algun pedo en safe_PartialTrace"; stop
    endif
    allocate (staout(0:2**(qubits-n_final)-1,0:2**n_final-1), stat=err(2))
    allocate (st1(0:2**(qubits-n_final)-1),st2(0:2**(qubits-n_final)-1), stat=err(3))
    if (any(err /= 0)) then; print*,"error en PartialTrace"; stop ; endif
    call unmsta(statin,staout,nwhich,qubits)
    do j1=0,2**n_final-1
       do j2=0,2**n_final-1
          rho(j1,j2)=dot_product(staout(:,j2),staout(:,j1))
       enddo
    enddo
    deallocate(staout, st1, st2, stat=err(1))
    return
    contains
      subroutine unmsta(statin,staout,nwhich,qubits)
        implicit none
        integer, intent(in)             :: nwhich,qubits
        complex(kind(1d0)), intent(in)  :: statin(:)
        complex(kind(1d0)), intent(out) :: staout(:,:)
        integer                         :: j, ncol, nrow
        staout=0d0
        do j=0,size( statin)-1
          call exdig(j,qubits,ncol,nrow,nwhich)
          staout(nrow+1,ncol+1)=statin(j+1)
        enddo
        return
      end subroutine unmsta
  end subroutine safe_PartialTrace
  subroutine PartialTrace(statin,rho,nwhich)
    ! DANGER:: I am leaving qubits indicated by nwhich
    ! If the system is of 7 qubits, and mwhich is
    ! 0 0 1 0 0 1 1 = 19 = nwhich
    ! 6 5 4 3 2 1 0 => qubits             
    ! i trace out qubits 2, 3, 6 and 6, leaving qubits 0, 1 and 4.
    implicit none
    complex(kind(1d0))                         ::  statin(:)
    complex(kind(1d0)), pointer                ::  rho(:,:)
    integer, intent(in)                        ::  nwhich
    integer                                    ::  err(3), n_final,j1,j2,qubits
    complex(kind(1d0)), pointer                ::  staout(:,:), st1(:), st2(:)
    qubits=IntegerLogBase2(size(statin))
    n_final=bits_on_one(nwhich)
    allocate (rho(0:2**n_final-1,0:2**n_final-1), stat=err(1))
    allocate (staout(0:2**(qubits-n_final)-1,0:2**n_final-1), stat=err(2))
    allocate (st1(0:2**(qubits-n_final)-1),st2(0:2**(qubits-n_final)-1), stat=err(3))
    if (any(err /= 0)) then; print*,"error en PartialTrace"; stop ; endif
    allocate (st1(0:2**(qubits-n_final)-1),st2(0:2**(qubits-n_final)-1), stat=err(3))
    call unmsta(statin,staout,nwhich,qubits)
    do j1=0,2**n_final-1
       do j2=0,2**n_final-1
          rho(j1,j2)=dot_product(staout(:,j2),staout(:,j1))
       enddo
    enddo
    deallocate(staout, st1, st2, stat=err(1))
    return
    contains
      subroutine unmsta(statin,staout,nwhich,qubits)
        implicit none
        integer, intent(in)             :: nwhich,qubits
        complex(kind(1d0)), intent(in)  :: statin(:)
        complex(kind(1d0)), intent(out) :: staout(:,:)
        integer                         :: j, ncol, nrow
        staout=0d0
        do j=0,size( statin)-1
          call exdig(j,qubits,ncol,nrow,nwhich)
          staout(nrow+1,ncol+1)=statin(j+1)
        enddo
        return
      end subroutine unmsta
  end subroutine PartialTrace
end module  quantumwau
