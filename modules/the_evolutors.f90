module EvolutionOperatorsF90
  use math
  implicit none
contains
  subroutine Partial_Close_Ring_Env_All(state,Coupling_Positions,J_coupling,J_env,J_Closing,b_env,b_cen) 
    !use EvolutionOperatorsF90
    complex(kind(1d0)), intent(inout)   :: state(:)
    real(kind(1d0)), intent(in)         :: J_env,J_coupling,J_closing,b_env(:),b_cen(:)
    integer, intent(in)                 :: Coupling_Positions(0:)
    integer                             :: qubits_cen,qubits,Target_1,Target_2,which_qubit
    qubits=IntegerLogBase2(size(state))
    qubits_cen=size(Coupling_Positions)
    if (size(b_env).ne.3.or.size(b_cen).ne.3)   stop "pedo 1 en many_eye_Spectator"
    if (qubits_cen.gt.qubits)                   stop "pedo 2 en many_eye_Spectator"
    ! La evolucion libre de la cadena
    do which_qubit=qubits_cen,qubits-2; Target_1=which_qubit; Target_2=which_qubit+1; call Action_Ising_Interaction(state,Target_1,Target_2,J_env); enddo
    Target_1=qubits_cen; Target_2=qubits-1; call Action_Ising_Interaction(state,Target_1,Target_2,J_closing)
    ! La interaccion del qubit con la cadena
    do  which_qubit=0,qubits_cen-1;call Action_Ising_Interaction(state,which_qubit,Coupling_Positions(which_qubit)+qubits_cen,J_coupling); enddo
    ! Patada magnetica a la cadena
    do which_qubit=0,qubits_cen-1;      call Application_Of_Kick(state,b_cen,which_qubit); enddo
    do which_qubit=qubits_cen,qubits-1; call Application_Of_Kick(state,b_env,which_qubit); enddo
  end subroutine Partial_Close_Ring_Env_All
  subroutine Partial_Close_Ring_Env_Spectator(state,qubits_cen,coupled_qubit,Coupling_Position,J_coupling,J_env,J_Closing,b_env,b_cen) 
    !use EvolutionOperatorsF90
    complex(kind(1d0)), intent(inout)   :: state(:)
    real(kind(1d0)), intent(in)         :: J_env,J_coupling,J_closing,b_env(:),b_cen(:)
    integer, intent(in)                 :: qubits_cen,coupled_qubit,Coupling_Position
    integer                             :: qubits,Target_1,Target_2,which_qubit
    qubits=IntegerLogBase2(size(state))
    if (size(b_env).ne.3.or.size(b_cen).ne.3)   stop "pedo 1 en many_eye_Spectator"
    if (qubits_cen.gt.qubits)                   stop "pedo 2 en many_eye_Spectator"
    if (coupled_qubit.gt.qubits_cen)            stop "pedo 3 en many_eye_Spectator"
    if (Coupling_Position.gt.qubits-qubits_cen) stop "pedo 4 en many_eye_Spectator"
    ! La evolucion libre de la cadena
    do which_qubit=qubits_cen,qubits-2; Target_1=which_qubit; Target_2=which_qubit+1; call Action_Ising_Interaction(state,Target_1,Target_2,J_env); enddo
    Target_1=qubits_cen; Target_2=qubits-1; call Action_Ising_Interaction(state,Target_1,Target_2,J_closing)
    ! La interaccion del qubit con la cadena
    call Action_Ising_Interaction(state,Coupled_Qubit,Coupling_Position+qubits_cen,J_coupling)
    ! Patada magnetica a la cadena
    do which_qubit=0,qubits_cen-1;      call Application_Of_Kick(state,b_cen,which_qubit); enddo
    do which_qubit=qubits_cen,qubits-1; call Application_Of_Kick(state,b_env,which_qubit); enddo
  end subroutine Partial_Close_Ring_Env_Spectator
  subroutine many_eye_Spectator(state,J_Ising,J_coupling,uncoupled_qubits,CoupledQubit,Coupling_Position,b)
    !use EvolutionOperatorsF90
    complex(kind(1d0)), intent(inout)   :: state(:)
    real(kind(1d0)), intent(in)         :: J_Ising,J_coupling,b(:)
    integer, intent(in)                 :: uncoupled_qubits,Coupling_Position,CoupledQubit
    integer                             :: which_qubit,number_of_qubits,Target_1,Target_2
    number_of_qubits=IntegerLogBase2(size(state))
    if (size(b).ne.3) stop "pedo en many_eye_Spectator"
    ! La evolucion libre de la cadena
    do which_qubit=uncoupled_qubits,number_of_qubits-1
      Target_1=which_qubit; Target_2=which_qubit+1; If(Target_2.eq.number_of_qubits) Target_2=uncoupled_qubits
      !print*,"Durante la evolucion libre interaccionan los qubits",Target_1,Target_2
      call Action_Ising_Interaction(state,Target_1,Target_2,J_Ising)
    enddo
    ! La interaccion del qubit con la cadena
    call Action_Ising_Interaction(state,CoupledQubit,Coupling_Position,J_coupling)
    print*,"many_eye_Spectator",CoupledQubit,Coupling_Position,J_coupling
    ! Patada magnetica a la cadena
    do which_qubit=uncoupled_qubits,number_of_qubits-1
      call Application_Of_Kick(state,b,which_qubit)
    enddo
  end subroutine many_eye_Spectator
  subroutine many_eye_Ising_Chain(state,J_Ising,J_coupling,uncoupled_qubits,Coupling_Positions,b)
    !use EvolutionOperatorsF90
    complex(kind(1d0)), intent(inout)   :: state(:)
    real(kind(1d0)), intent(in)         :: J_Ising,J_coupling,b(:)
    integer, intent(in)                 :: uncoupled_qubits,Coupling_Positions(0:)

    integer                             :: which_qubit,number_of_qubits,Target_1,Target_2
    if (size(b).ne.3) stop "pedo en many_Ising_Chain"
    if (size(Coupling_Positions).ne.uncoupled_qubits) then;
      print*,"algun pedo yyy en many_eye_Ising_Chain",size(Coupling_Positions),uncoupled_qubits; stop
    endif

    number_of_qubits=IntegerLogBase2(size(state))
    ! La evolucion libre de la cadena
    do which_qubit=uncoupled_qubits,number_of_qubits-1
      Target_1=which_qubit; Target_2=which_qubit+1; If(Target_2.eq.number_of_qubits) Target_2=uncoupled_qubits
      !print*,"Durante la evolucion libre interaccionan los qubits",Target_1,Target_2
      call Action_Ising_Interaction(state,Target_1,Target_2,J_Ising)
    enddo
    ! La interaccion de cada qubit con la cadena
    do which_qubit=0,uncoupled_qubits-1
      Target_1=which_qubit;  Target_2=Coupling_Positions(which_qubit)
      if (Target_2<uncoupled_qubits.or.Target_2>number_of_qubits-1) then
        print*,"algun pedo xxx en many_eye_Ising_Chain",Target_2,uncoupled_qubits,number_of_qubits; stop
      endif
      !print*,"Durante el acoplamiento interaccionan los qubits",Target_1,Target_2
      call Action_Ising_Interaction(state,Target_1,Target_2,J_coupling)
    enddo
    ! Patada magnetica a la cadena
    do which_qubit=uncoupled_qubits,number_of_qubits-1
      call Application_Of_Kick(state,b,which_qubit)
    enddo

  end subroutine many_eye_Ising_Chain
  subroutine H_plus_Homogeneous(state,Magnetick_Kick,Ising_Strength)
    complex(kind(1d0)),intent(inout)  :: state(:)
    real(kind(1d0)),intent(in)        :: Magnetick_Kick(:),Ising_Strength
    complex(kind(1d0)),allocatable    :: state_p(:),state_m(:)
    allocate(state_p(size(state)),state_m(size(state)))
    !Hminusscalar=I*.5d0*(statePlus-stateminus)
    state_P=state;    state_M=state
    call EvolutionOperatorHomogeneous(state_p,Magnetick_Kick,Ising_Strength)
    call BackEvolutionOperatorHomogeneous(state_m,Magnetick_Kick,Ising_Strength)
    state=.5d0*(state_p+state_m)
    deallocate(state_p,state_m)
  end subroutine H_plus_Homogeneous
  subroutine H_minus_Homogeneous(state,Magnetick_Kick,Ising_Strength)
    complex(kind(1d0)),intent(inout)  :: state(:)
    real(kind(1d0)),intent(in)        :: Magnetick_Kick(:),Ising_Strength
    complex(kind(1d0)),allocatable    :: state_p(:),state_m(:)
    allocate(state_p(size(state)),state_m(size(state)))
    state_P=state;    state_M=state
    call EvolutionOperatorHomogeneous(state_p,Magnetick_Kick,Ising_Strength)
    call BackEvolutionOperatorHomogeneous(state_m,Magnetick_Kick,Ising_Strength)
    stop 'creo que aca teng un problema de signo, no se si va un signo menos adelante'
    state=I*.5d0*(state_p-state_m)
    deallocate(state_p,state_m)
  end subroutine H_minus_Homogeneous
  subroutine BackEvolutionOperatorHomogeneous(state,Magnetick_Kick,Ising_Strength)
    complex(kind(1d0)),intent(inout)  :: state(:)
    real(kind(1d0)),intent(in)        :: Magnetick_Kick(:),Ising_Strength
    call Application_Of_Kick_Homogeneous(state,-Magnetick_Kick)
    call Ising_Interaction_Homogeneous(state,-Ising_Strength)
  end subroutine BackEvolutionOperatorHomogeneous
  subroutine EvolutionOperatorHomogeneous(state,Magnetick_Kick,Ising_Strength)
    complex(kind(1d0)),intent(inout)  :: state(:)
    real(kind(1d0)),intent(in)        :: Magnetick_Kick(:),Ising_Strength
    call Ising_Interaction_Homogeneous(state,Ising_Strength)
    call Application_Of_Kick_Homogeneous(state,Magnetick_Kick)
  end subroutine EvolutionOperatorHomogeneous
  subroutine Ising_Interaction_Homogeneous(state,Ising_Strength)
    complex(kind(1d0)),intent(inout)  :: state(:)
    real(kind(1d0)),intent(in)        :: Ising_Strength
    integer                           :: Spin,qubits
    qubits=IntegerLogBase2(size(state))
    do Spin=0,qubits-2
       call Action_Ising_Interaction(state,Spin,Spin+1,Ising_Strength)
    end do
    call Action_Ising_Interaction(state,qubits-1,0,Ising_Strength)
  end subroutine Ising_Interaction_Homogeneous
  subroutine Application_Of_Kick_Homogeneous(state,Magnetick_Kick)
    complex(kind(1d0)),intent(inout)  :: state(:)
    real(kind(1d0)),intent(in)        :: Magnetick_Kick(:)
    integer                           :: qubits,Spin
    qubits=IntegerLogBase2(size(state))
    do Spin=0,qubits-1
       call Application_Of_Kick(state,Magnetick_Kick,Spin)
    end do
  end subroutine Application_Of_Kick_Homogeneous
  subroutine Application_Of_Kick(state,Magnetick_Kick,SpinToBeKicked)
    ! Esta rutina supone que la base es tal que Sigma_x = diag (1,-1) y de ahi ciclico (ver por ejemplo
    ! Chuan pagina... no se... pero las matrices de pauli que estan ahi son las buenas)
    ! es la implementacion de la operacion Exp[-I b.sigma^j]
    use bit_manipulation, only : AddOneZero
    complex(kind(1d0)),intent(inout)  :: state(0:)
    real(kind(1d0)),intent(in)        :: Magnetick_Kick(:)
    complex(kind(1d0))                :: Matrix_For_Kick(2,2)
    complex(kind(1d0))                :: z(2)
    integer,intent(in)                :: SpinToBeKicked
    integer                           :: UpperNumber,LowerNumber,qubits,j_left,j_right,left_digit
    if(2**(SpinToBeKicked)>=size(state)) then; print*,"error",SpinToBeKicked,size(state);  stop; end if
    if(norma(Magnetick_Kick).eq.0d0) return
    if(norma(Magnetick_Kick).lt.10d-7) stop "posible error en la runita Application_Of_Kick, mire bien los parametros"
    Matrix_For_Kick=Matrix_For_Magnetic_Kick(Magnetick_Kick)
    qubits=IntegerLogBase2(size(state))
    do j_left=0,2**(qubits-1-SpinToBeKicked)-1
       left_digit=ishft(j_left,SpinToBeKicked+1)
       do j_right=0,2**SpinToBeKicked-1
          LowerNumber=left_digit+j_right
          UpperNumber=IBset(LowerNumber,SpinToBeKicked)
          z(1)=state(LowerNumber); z(2)=state(UpperNumber)
          z=MatMul(Matrix_For_Kick,z)
          state(LowerNumber)=z(1); state(UpperNumber)=z(2)
       end do
    end do
  end subroutine Application_Of_Kick
  subroutine Action_Ising_Interaction(state,Target_Qubit_1,Target_Qubit_2,Ising_Strength)
    ! esta es la rutina de interaccion de ising en la direccion x, donde es diagonal 
    ! la matriz de Pauli correspondiente (ver doc de subroutine Application_Of_Kick).
    ! La entrada es un estado de varios qubits, dos enteros diferentes de 0 a el numero de qubits,
    ! y la fuerza de la interaccion. La salida es 
    ! Exp[-I J sigma_x^j sigma_x^j'] |state>
    ! donde j y j' son los indices de las posiciones relevantes: Target_Qubit_1,Target_Qubit_2
    complex(kind(1d0)),intent(inout)  :: state(0:)
    real(kind(1d0)),intent(in)        :: Ising_Strength
    integer,intent(in)                :: Target_Qubit_1,Target_Qubit_2
    integer                           :: J
    complex(kind(1d0))                :: exp_minus,exp_plus
    if((2**(Target_Qubit_1)>=size(state)).or.(2**(Target_Qubit_2)>=size(state))) then
       print*,"error",Target_Qubit_1,Target_Qubit_2,size(state);  stop
    end if
    if((Target_Qubit_1==Target_Qubit_2).or.(Target_Qubit_1<0.or.Target_Qubit_2<0)) then
       print*,"error",Target_Qubit_1,Target_Qubit_2;  stop
    end if
    exp_plus=exp(I*Ising_Strength);  exp_minus=conjg(exp_plus)
    do j=0,size(state)-1
       If(btest(j,Target_Qubit_1).eqv.btest(j,Target_Qubit_2)) then
          state(j)=exp_minus*state(j)
       else
          state(j)=exp_plus*state(j)
       end If
    end do
  end subroutine Action_Ising_Interaction
  function Matrix_For_Magnetic_Kick(b)
    complex(kind(1d0))         :: Matrix_For_Magnetic_Kick(2,2)
    real(kind(1d0)),intent(in) :: b(3)
    real(kind(1d0))            :: n(3),phi
    phi=norma(b);    n=b/phi
    Matrix_For_Magnetic_Kick(1,1)=cos(phi)-I*n(1)*sin(phi)
    Matrix_For_Magnetic_Kick(2,2)=conjg(Matrix_For_Magnetic_Kick(1,1))
    Matrix_For_Magnetic_Kick(1,2)=-(n(3)+I*n(2))*sin(phi)
    Matrix_For_Magnetic_Kick(2,1)=-conjg(Matrix_For_Magnetic_Kick(1,2))
  end function Matrix_For_Magnetic_Kick
  subroutine EO_SymmetricHomogeneous(state,Magnetick_Kick,Ising_Strength)
    complex(kind(1d0)),intent(inout)  :: state(0:)
    real(kind(1d0)),intent(in)        :: Magnetick_Kick(:),Ising_Strength
    call Application_Of_Kick_Homogeneous(state,Magnetick_Kick/2)
    call Ising_Interaction_Homogeneous(state,Ising_Strength)
    call Application_Of_Kick_Homogeneous(state,Magnetick_Kick/2)
  end subroutine EO_SymmetricHomogeneous
  subroutine BackwardsMultipleEvolution(J,H,state)
    implicit none    
    complex(kind(1d0)), intent(inout) :: state(:)
    real(kind(1d0)), intent(in)       :: J(:)
    real(kind(1d0)), intent(in)       :: H(:,:)
    integer                           :: j1
    if ((size(h(1,:)).ne.3).or.(size(h(:,1)).ne.size(j))) then
       print*," Error en MultipleEvolution",size(h(1,:)),size(h(:,1)),size(j)
       stop
    end if
    do j1=size(J),1,-1
       call Application_Of_Kick_Homogeneous(state,-H(j1,:))
       call Ising_Interaction_Homogeneous(state,-J(j1))
    end do 
   end subroutine BackwardsMultipleEvolution 
  subroutine MultipleEvolution(J,H,state)
    implicit none    
    complex(kind(1d0)), intent(inout) :: state(:)
    real(kind(1d0)), intent(in)       :: J(:)
    real(kind(1d0)), intent(in)       :: H(:,:)
    integer                           :: j1
    if ((size(h(1,:)).ne.3).or.(size(h(:,1)).ne.size(j))) then
       print*," Error en MultipleEvolution",size(h(1,:)),size(h(:,1)),size(j)
       stop
    end if
    do j1=1,size(J)
       call Ising_Interaction_Homogeneous(state,J(j1))
       !call IsingXgate(J(j1),state)
       call Application_Of_Kick_Homogeneous(state,H(j1,:))
       !call MagnetickKick(H(j1,:),state)
    end do 
   end subroutine MultipleEvolution 
end module EvolutionOperatorsF90
