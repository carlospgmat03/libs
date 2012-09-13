module quantum
  use definitions
  use EvolutionOperators
  use bit_manipulation
  use math
  implicit none
contains

  subroutine EnsembleBitFlip_Channelj(Ensemblestates,qubit_to_act,probability)
    complex(kind(1d0)), intent(inout)           :: ensemblestates(:,:)
    integer, intent(in)                         :: qubit_to_act
    real(kind(1d0)), intent(in)                 :: probability
    integer                                     :: j1
    do j1=1,size(ensemblestates(:,1))
       ensemblestates(j1,:)=BitFlip_Channelj(ensemblestates(j1,:),qubit_to_act,probability)
    end do
  end subroutine EnsembleBitFlip_Channelj

  function BitFlip_Channelj(state,qubit_to_act,probability)
    complex(kind(1d0)), intent(in)           :: state(0:)
    real(kind(1d0)), intent(in)              :: probability
    integer, intent(in)                      :: qubit_to_act
    complex(kind(1d0))                       :: BitFlip_Channelj(0:size(state)-1)
    real(kind(1d0))                          :: x
    if (probability<0d0.or.probability>1) then
       print*,"Error en BitFlip_Channel_individual_qubit, probability-",probability
    end if
    call random_number(x)
    if (x>probability) then
       BitFlip_Channelj=state
    else
       BitFlip_Channelj=BitFlip_individual_qubit(state,qubit_to_act)
    end if
  end function BitFlip_Channelj

  function PhaseFlip_Channelj(state,qubit_to_act,probability)
    complex(kind(1d0)), intent(in)           :: state(0:)
    real(kind(1d0)), intent(in)              :: probability
    integer, intent(in)                      :: qubit_to_act
    complex(kind(1d0))                       :: PhaseFlip_Channelj(0:size(state)-1)
    real(kind(1d0))                          :: x
    if (probability<0d0.or.probability>1) then
       print*,"Error en PhaseFlip_Channel_individual_qubit, probability-",probability
    end if
    call random_number(x)
    if (x>probability) then
       PhaseFlip_Channelj=state
    else
       PhaseFlip_Channelj=PhaseFlip_individual_qubit(state,qubit_to_act)
    end if
  end function PhaseFlip_Channelj

  function BitFlip_individual_qubit(state,qubit_to_act)
    implicit none
    complex(kind(1d0)), intent(in)           :: state(0:)
    complex(kind(1d0))                       :: BitFlip_individual_qubit(0:size(state)-1)
    integer, intent(in)                      :: qubit_to_act
    integer                                  :: j1,j2,power2j
    power2j=2**qubit_to_act
    do j1=0,size(state)-1
       j2=ieor(power2j,j1)
       BitFlip_individual_qubit(j2)=state(j1)
    end do
  end function BitFlip_individual_qubit

  function PhaseFlip_individual_qubit(state,qubit_to_act)
    implicit none
    complex(kind(1d0)), intent(in)           :: state(0:)
    complex(kind(1d0))                       :: PhaseFlip_individual_qubit(0:size(state)-1)
    integer, intent(in)                      :: qubit_to_act
    integer                                  :: j1
    PhaseFlip_individual_qubit=state
    do j1=0,size(state)-1
       if (btest(j1,qubit_to_act)) PhaseFlip_individual_qubit(j1)=-state(j1)
    end do
  end function PhaseFlip_individual_qubit

  subroutine MeasurmentFirstQubit(state,MeasurmentResult)
    complex(kind(1d0)), intent(inout)           :: state(0:)
    real (kind(1d0))                            :: x
    logical                                     :: MeasurmentResult
    call random_number(x)
    if (ProbablityFirstQubit(state,.true.)>x) then
       MeasurmentResult=.true.
       state=ProyectionFirstQubit(state,.true.)
    else
       MeasurmentResult=.false.
       state=ProyectionFirstQubit(state,.false.)
    end if
  end subroutine MeasurmentFirstQubit

  function ProyectionFirstQubit(state,updown)
    ! Hacer rutina que me mida la probabilidad de que quede en cero o en uno
    complex(kind(1d0)), intent(in)           :: state(0:)
    complex(kind(1d0))                       :: ProyectionFirstQubit(0:size(state)-1)
    logical, intent(in)                      :: updown
    real(kind(1d0))                          :: norm
    ProyectionFirstQubit=state
    If (updown) ProyectionFirstQubit(0::2)=0d0
    If (.not.updown) ProyectionFirstQubit(1::2)=0d0
    norm=norma(ProyectionFirstQubit)
    if (norm.eq.0d0) then
       print*,"Error en ProyectionFirstQubit",norm;  stop
    end if
    ProyectionFirstQubit=ProyectionFirstQubit/norma(ProyectionFirstQubit)
  end function ProyectionFirstQubit

  function ProbablityFirstQubit(state,updown)
    ! If updown is true, measures that the spin is in state "|1>"
    ! Beware, it can give negative probalilities due to numerical error
    real(kind(1d0))                          :: ProbablityFirstQubit
    complex(kind(1d0)), intent(in)           :: state(0:)
    logical, intent(in)                      :: updown
    ProbablityFirstQubit=sum(abs(state(0::2))**2)
    if (updown) ProbablityFirstQubit=1-ProbablityFirstQubit
  end function ProbablityFirstQubit

  subroutine print_state(state)
    implicit none
    complex(kind(1d0)), intent(in) :: state(0:)
    integer                        :: qubits,j1
    logical, allocatable           :: aux(:)
    qubits=IntegerLog(size(state),2)
    allocate(aux(0:qubits-1))
    do j1=0,size(state)-1
       aux=get_bits(j1,qubits)
      !aux=aux(qubits-1:0:1)
       write(*,*)j1,aux(qubits-1:0:-1),state(j1)
      !write(*,'(I4,2F7.3)')j1,state(j1),aux(qubits-1:0:-1)
      !write(*,'(I4,2X,qubitsL1,2F7.3)')j1,state(j1),aux(qubits-1:0:-1)
    end do
    deallocate(aux)
  end subroutine print_state

  function Purity(rho)
    implicit none
    real (kind(1d0))  :: Purity
    complex(kind(1d0)), intent(in)::  rho(:,:)
    Purity=MatrixTrace(matmul(rho,rho))
  end function Purity
  
  subroutine Merge_States(state_in_1,state_in_2,nwhich, state_out)
    implicit none
    integer, intent(in)                                :: nwhich
    complex(kind(1d0)), intent(in)                     :: state_in_1(0:),state_in_2(0:)
    complex(kind(1d0)), intent(out)                    :: state_out(0:)
    integer                                            :: qin1, qin2,qout,j1, n1out,n2out
    qin1=IntegerLog(size(state_in_1),2)
    qin2=IntegerLog(size(state_in_2),2)
    qout=IntegerLog(size(state_out),2)
    if ((qin1+qin2.ne.qout).or. (qin1.ne.bits_on_one(nwhich))) then
       print*,"error en merge states",qin1,qin2,qout,bits_on_one(nwhich),nwhich
       stop
    end if
    do j1=0,2**qout-1
       call  exdig(j1,qout,n1out,n2out,nwhich)
       state_out(j1)=state_in_1(n1out)*state_in_2(n2out)
       !       print*,"in mergestates"j1,n1out,n2out,&
       !           abs(state_in_1(n1out)),abs(state_in_2(n2out)),abs(state_out(j1))
    end do
  end subroutine Merge_States
  
  
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
    integer                                    ::  err(3), n_final,j1,j2
    complex(kind(1d0)), pointer                ::  staout(:,:), st1(:), st2(:)
    n_final=bits_on_one(nwhich)
    allocate (rho(0:2**n_final-1,0:2**n_final-1), stat=err(1))
    allocate (staout(0:2**(q-n_final)-1,0:2**n_final-1), stat=err(2))
    allocate (st1(0:2**(q-n_final)-1),st2(0:2**(q-n_final)-1), stat=err(3))
    if (any(err /= 0)) then; print*,"error en PartialTrace"; stop ; endif
    allocate (st1(0:2**(q-n_final)-1),st2(0:2**(q-n_final)-1), stat=err(3))
    call unmsta(statin,staout,nwhich)
    do j1=0,2**n_final-1
       do j2=0,2**n_final-1
          rho(j1,j2)=dot_product(staout(:,j2),staout(:,j1))
       enddo
    enddo
    deallocate(staout, st1, st2, stat=err(1))
    return
  contains
    subroutine unmsta(statin,staout,nwhich)
      implicit none
      integer, intent(in)  :: nwhich
      complex(kind(1d0)), intent(in)  :: statin(:)
      complex(kind(1d0)), intent(out) :: staout(:,:)
      integer            :: j, ncol, nrow
      staout=0d0
      do j=0,size( statin)-1
         call exdig(j,q,ncol,nrow,nwhich)
         staout(nrow+1,ncol+1)=statin(j+1)
      enddo
      return
    end subroutine unmsta
  end subroutine PartialTrace
  
 
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
    ! 14. random real
    implicit none
    complex(kind(1d0)), intent(out) :: state(0:)
    complex(kind(1d0))              :: x1,x2
    complex(kind(1d0)), allocatable :: state_aux1(:), state_aux2(:), state_aux3(:)
    integer, intent(in)             :: option, internal_options(:)
    integer                         :: j1, nwhich, nb, bath1, sbath1
    if (option==1) then
       do j1=0,size(state)-1
          state(j1)=RandomGaussian()
       enddo
       state=state/norma(state)
    else if (option==2) then
       state=0d0
       state(internal_options(1))=1d0
    else if (option==3) then
       allocate (state_aux1(0:2**(q-2)-1),state_aux2(0:3))
       call state_preparation(state_aux1,1,(/0/))
       call random_bell(state_aux2)
       call Merge_States(state_aux1,state_aux2,2**q-1-internal_options(1), state)
       deallocate(state_aux1,state_aux2)
    else if (option==4) then
      nwhich=internal_options(1)
      nb=bits_on_one(nwhich)
      allocate (state_aux2(0:2**(q-nb)-1),state_aux1(0:2**nb-1))
      call state_preparation(state_aux1,1,(/0/))
      call state_preparation(state_aux2,1,(/0/))
      !print*,"hola state preparation option 4", nwhich
      call Merge_States(state_aux1,state_aux2,nwhich, state)
      deallocate(state_aux1,state_aux2)
    else if (option==5) then
       open(unit=98,file="state.dat",status="old")
       do j1=0,2**q-1
          read(98,*)x1,x2
          state(j1)=x1+I*x2
       end do
       state=state/norma(state)
       close(98)
    else if (option==6) then
      allocate(state_aux1(0:1))
      allocate(state_aux3(0:0))
      state_aux3=1d0
      do j1=1,IntegerLog(size(state(:)),2)
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
       allocate (state_aux2(0:2**(q-nb)-1),state_aux1(0:2**nb-1))
       call state_preparation(state_aux1,6,(/0/))
       call state_preparation(state_aux2,1,(/0/))
       !print*,"hola state preparation option 4", nwhich
       call Merge_States(state_aux1,state_aux2,nwhich, state)
       state=state/norma(state)
       deallocate(state_aux1,state_aux2)
    else if (option==8) then
       nwhich=internal_options(1)
       allocate (state_aux2(0:2**(q-1)-1),state_aux1(0:1))
       state_aux1=1d0/sqrt(2d0)
       call state_preparation(state_aux2,1,(/0/))
       call Merge_States(state_aux1,state_aux2,nwhich, state)
       deallocate(state_aux1,state_aux2)
    else if (option==9) then
       nwhich=internal_options(1)
       allocate (state_aux2(0:2**(q-1)-1),state_aux1(0:1))
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
       allocate (state_aux1(0:2**(q-2)-1),state_aux2(0:3))
       call state_preparation(state_aux1,2,(/internal_options(2)/))
       call random_bell(state_aux2)
       call Merge_States(state_aux1,state_aux2,2**q-1-internal_options(1), state)
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
      if ((mod(bath1,4).ne.0).or.(sbath1>=q-2)) then
        print*,"error en state preparation, opcion 11", bath1
        stop
      end if
      allocate (state_aux1(0:2**(sbath1)-1),state_aux2(0:2**(q-sbath1-2)-1),state_aux3(0:2**(q-2)-1))
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
       allocate (state_aux2(0:2**(q-nb)-1),state_aux1(0:2**nb-1))
       call state_preparation(state_aux1,1,(/0/))
       call state_preparation(state_aux2,2,(/internal_options(2)/))
       !print*,"hola state preparation option 4", nwhich
       call Merge_States(state_aux1,state_aux2,nwhich, state)
       deallocate(state_aux1,state_aux2)
    else if (option==14) then
       call state_preparation(state,1,(/0/));       state=abs(state)
    else
       print*,"Bad option in state preparation"
       stop
    endif
    return
  end subroutine state_preparation
     
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
     
  subroutine random_bell(state)
    implicit none
    complex(kind(1d0)), intent(out)            :: state(0:3)
    real(kind(1d0))                            :: phi(2), theta(2),global_phase,theta_rel
    !real(kind(1d0))                            :: xhi(2), x(2)
    call random_number(global_phase); global_phase=2*pi*global_phase
    call random_number(phi);          phi=2*pi*phi
    call random_number(theta);        theta=2*theta-1d0;               where (theta==theta) theta=acos(theta)
    theta_rel=(theta(1)-theta(2))/2
    state(0)=cos(theta_rel)
    state(1)=exp(I*phi(2))*sin(theta_rel)
    state(2)=exp(I*phi(1))*sin(theta_rel)
    state(3)=-exp(I*sum(phi))*state(0)
    state=exp(I*global_phase)*state/sqrt(2d0)
    !state=0d0; state(0)=1/sqrt(2d0);state(3)=1/sqrt(2d0);
  end subroutine random_bell

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

  function Concurrence(rho)
    use diagslow
    !    use lin_eig_gen_int
    implicit none
    real (kind(1d0))                ::  concurrence
    complex(kind(1d0)),intent(in)   ::  rho(0:,0:)
    complex(kind(1d0)), allocatable ::  rho_tilde(:,:)
    real(kind(1d0)), allocatable    ::  eigen_values(:)
    integer                         ::  qubits, j1,j2, qubits2, err
    qubits=IntegerLog(size(rho(1,:)),2)
    qubits2=2**qubits
    if (mod(qubits,2).ne.0) then ;print*,"error en la concurrencia!", qubits; stop; end if
    allocate(rho_tilde(0:qubits2-1,0:qubits2-1))
    allocate(eigen_values(0:qubits2-1), stat=err); if (err/=0) print*,"mierda"
    do j1=0,2**qubits-1
       do j2=0,2**qubits-1
          rho_tilde(j1,j2)=conjg(rho(qubits2-j1-1,qubits2-j2-1))
          if (mod(bits_on_one(j1)+bits_on_one(j2),2)==1) rho_tilde(j1,j2)=-rho_tilde(j1,j2) 
       end do
    end do
    rho_tilde=matmul(rho,rho_tilde)
    !    call matrix_slow_diagonalize(rho_tilde, eigen_values) 
    call eigenvalues_FOR_CONCURRENCE(rho_tilde,eigen_values)
    if ( maxval(eigen_values)<0) then 
       Concurrence=0
    else
       Concurrence=2*sqrt(   (maxval(eigen_values) ))-&
            sum(sqrt(abs(eigen_values)))
    endif
    deallocate(rho_tilde,eigen_values, stat=err); if (err/=0) print*,"mierda"
  end function Concurrence
end module quantum
