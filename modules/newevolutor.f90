module EvolutionOperators
  use definitions
  use math
  use bit_manipulation
  implicit none
  interface EnsembleEvolutionOperator
     module procedure EnsembleEvolutionOperatorvJvB
  end interface
  interface EvolutionOperator
     module procedure EvolutionOperatorvJ
     module procedure EvolutionOperatorvJvB
     module procedure EvolutionOperatorscalar
     module procedure EvolutionOperatorMKVJVB
  end interface
  interface BackwardsEvolution
     module procedure BackwardsEvolutionscalar
     module procedure BackwardsEvolutionMKVJVB
  end interface
  interface IsingXGate
     module procedure IsingXGate_scalar
     module procedure IsingXGate_vector
  end interface
  interface MagnetickKick
     module procedure MagnetickKick_scalar
     module procedure MagnetickKick_vector
  end interface
  interface Hplus
     module procedure Hplusscalar
     module procedure HplusMKVJVB
  end interface
  interface Hminus
     module procedure Hminusscalar
     module procedure HminusMKVJVB
  end interface
contains
  function Hminusscalar(sc,state)
    implicit none
    type (system_car_scalar), intent(in) :: sc
    complex(kind(1d0)), intent(in)       :: state(:)
    complex(kind(1d0))                   :: Hminusscalar(size(state))
    complex(kind(1d0))                   :: statePlus(size(state))
    complex(kind(1d0))                   :: stateMinus(size(state))
    statePlus=state
    stateMinus=state
    call evolutionoperator(sc,StatePlus)
    call BackwardsEvolution(sc,StateMinus)
    Hminusscalar=I*.5d0*(statePlus-stateminus)
  end function Hminusscalar
  function HminusMKVJVB(sc,state)
    implicit none
    type (system_carvJvBvM), intent(in)  :: sc
    complex(kind(1d0)), intent(in)       :: state(:)
    complex(kind(1d0))                   :: HminusMKVJVB(size(state))
    complex(kind(1d0))                   :: statePlus(size(state))
    complex(kind(1d0))                   :: stateMinus(size(state))
    statePlus=state
    stateMinus=state
    call evolutionoperator(sc,StatePlus)
    call BackwardsEvolution(sc,StateMinus)
    HminusMKVJVB=I*.5d0*(statePlus-stateminus)
  end function HminusMKVJVB
  function HplusMKVJVB(sc,state)
    implicit none
    type (system_carvJvBvM), intent(in)  :: sc
    complex(kind(1d0)), intent(in)       :: state(:)
    complex(kind(1d0))                   :: HPlusMKVJVB(size(state))
    complex(kind(1d0))                   :: statePlus(size(state))
    complex(kind(1d0))                   :: stateMinus(size(state))
    statePlus=state;    stateMinus=state
    call evolutionoperator(sc,StatePlus)
    call BackwardsEvolution(sc,StateMinus)
    HplusMKVJVB=.5d0*(statePlus+stateminus)
  end function HplusMKVJVB
  function HplusMKVJVB_random_Phase(sc,state,random_phase)
    implicit none
    type (system_carvJvBvM), intent(in)  :: sc
    complex(kind(1d0)), intent(in)       :: state(:)
    complex(kind(1d0))                   :: HPlusMKVJVB_random_Phase(size(state))
    complex(kind(1d0))                   :: statePlus(size(state))
    complex(kind(1d0))                   :: stateMinus(size(state))
    real(kind(1d0))                      :: random_phase
    statePlus=state;    stateMinus=state
    call evolutionoperator(sc,StatePlus)   ; statePlus=exp(I*random_phase)*statePlus
    call BackwardsEvolution(sc,StateMinus) ; stateMinus=exp(-I*random_phase)*stateminus
    HplusMKVJVB_random_Phase=.5d0*(statePlus+stateminus)
  end function HplusMKVJVB_random_Phase
  function HminusMKVJVB_random_Phase(sc,state,random_phase)
    implicit none
    type (system_carvJvBvM), intent(in)  :: sc
    complex(kind(1d0)), intent(in)       :: state(:)
    complex(kind(1d0))                   :: HminusMKVJVB_random_Phase(size(state))
    complex(kind(1d0))                   :: statePlus(size(state))
    complex(kind(1d0))                   :: stateMinus(size(state))
    real(kind(1d0))                      :: random_phase
    statePlus=state
    stateMinus=state
    call evolutionoperator(sc,StatePlus)
    call BackwardsEvolution(sc,StateMinus)
    HminusMKVJVB_random_Phase=-I*.5d0*(exp(I*random_phase)*statePlus-exp(-I*random_phase)*stateminus)
!    HminusMKVJVB_random_Phase=I*.5d0*(exp(I*random_phase)*statePlus-exp(-I*random_phase)*stateminus)
  end function HminusMKVJVB_random_Phase
  function Hplusscalar(sc,state)
    implicit none
    type (system_car_scalar), intent(in) :: sc
    complex(kind(1d0)), intent(in)       :: state(:)
    complex(kind(1d0))                   :: HPlusscalar(size(state))
    complex(kind(1d0))                   :: statePlus(size(state))
    complex(kind(1d0))                   :: stateMinus(size(state))
    statePlus=state;    stateMinus=state
    call evolutionoperator(sc,StatePlus)
    call BackwardsEvolution(sc,StateMinus)
    Hplusscalar=.5d0*(statePlus+stateminus)
  end function Hplusscalar
  subroutine BackwardsEvolutionMKVJVB(sc,state)
    implicit none
    type (system_carvJvBvM), intent(in) :: sc
    type (system_carvJvB)               :: sc_aux
    complex(kind(1d0)), intent(inout)   :: state(0:)
    integer                             :: j
    sc_aux%qbits=sc%qbits
    allocate(sc_aux%J(sc_aux%qbits))
    allocate(sc_aux%h(3,sc_aux%qbits))
    do j=size(sc%j(:,1)),1,-1
       sc_aux%J=-sc%J(j,:);   sc_aux%h=-sc%h(j,:,:)
       call MagnetickKick(sc_aux%h,state)
       call IsingXgate(sc_aux%j,state)
    end do
    deallocate(sc_aux%J)
    deallocate(sc_aux%h)
  end subroutine BackwardsEvolutionMKVJVB
  subroutine EvolutionOperatorMKVJVB(sc,state)
    implicit none
    type (system_carvJvBvM), intent(in) :: sc
    type (system_carvJvB)               :: sc_aux
    complex(kind(1d0)), intent(inout)   :: state(0:)
    integer                             :: j
    sc_aux%qbits=sc%qbits
    allocate(sc_aux%J(sc_aux%qbits))
    allocate(sc_aux%h(3,sc_aux%qbits))
    do j=1,size(sc%j(:,1))
       sc_aux%J=sc%J(j,:);   sc_aux%h=sc%h(j,:,:)
       call EvolutionOperatorvJvB(sc_aux,state)
    end do
    deallocate(sc_aux%J)
    deallocate(sc_aux%h)
  end subroutine EvolutionOperatorMKVJVB
  subroutine BackwardsEvolutionScalar(sc,state)
    implicit none
    type (system_car_scalar), intent(in) :: sc
    type (system_car_scalar)             :: scaux
    complex(kind(1d0)), intent(inout) :: state(0:)
    scaux=sc; scaux%J=-scaux%J;scaux%h=-scaux%h
    call MagnetickKick(scaux%h,state)
    call IsingXgate(scaux%j,state)
  end subroutine BackwardsEvolutionScalar
  subroutine EvolutionOperatorRudi(sc,state)
    implicit none
    type (system_car_Rudi), intent(in) :: sc
!    real(kind(1d0))                   :: J1,J2,H1(3),H2(3)
    complex(kind(1d0)), intent(inout) :: state(0:)
    call IsingXgate(sc%j,state)
    call MagnetickKick( (/sc%h1,0d0,sc%h1/),state)
    call IsingXgate(sc%j,state)
    call MagnetickKick( (/sc%h2,sc%h2,0d0/),state)
  end subroutine EvolutionOperatorRudi
  subroutine EvolutionOperatorRudiGeneral(J1,J2,H1,H2,state)
    implicit none
!    type (system_car_scalar), intent(in) :: sc
    real(kind(1d0))                   :: J1,J2,H1(3),H2(3)
    complex(kind(1d0)), intent(inout) :: state(0:)
    call IsingXgate(J1,state)
    call MagnetickKick( H1,state)
    call IsingXgate(J2,state)
    call MagnetickKick( H2,state)
  end subroutine EvolutionOperatorRudiGeneral
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
       call IsingXgate(J(j1),state)
       call MagnetickKick(H(j1,:),state)
    end do
  end subroutine MultipleEvolution
  subroutine BackwardsMultipleEvolution(J,H,state)
    implicit none    
    complex(kind(1d0)), intent(inout) :: state(:)
    real(kind(1d0)), intent(in)       :: J(:)
    real(kind(1d0)), intent(in)       :: H(:,:)
    integer                           :: j1
    if ((size(h(1,:)).ne.3).or.(size(h(:,1)).ne.size(j))) then
       print*," Error en BackwardsMultipleEvolution",size(h(1,:)),size(h(:,1)),size(j)
       stop
    end if
    do j1=size(J),1,-1
       call MagnetickKick(-H(j1,:),state)
       call IsingXgate(-J(j1),state)
    end do
  end subroutine BackwardsMultipleEvolution
  subroutine EnsembleEvolutionOperatorvJvB(sc,EnsembleStates)
    implicit none
    type (system_carvJvB), intent(in) :: sc
    complex(kind(1d0)), intent(inout) :: EnsembleStates(:,:)
    integer                           :: j1,NumberOfMembersOfEnsemble
    NumberOfMembersOfEnsemble=size(EnsembleStates(:,1))
    do j1=1,NumberOfMembersOfEnsemble
       call EvolutionOperator(sc,EnsembleStates(j1,:))
    end do
  end subroutine EnsembleEvolutionOperatorvJvB
  subroutine EvolutionOperatorvJvB(sc,state)
    implicit none
    type (system_carvJvB), intent(in) :: sc
    complex(kind(1d0)), intent(inout) :: state(0:)
    call IsingXgate(sc%J,state)
    call MagnetickKick(sc%H,state)
  end subroutine EvolutionOperatorvJvB
  subroutine SymmetricEvolutionScalar(sc,state)
    implicit none
    type (system_car_scalar), intent(in) :: sc
    type (system_car_scalar)             :: scaux
    complex(kind(1d0)), intent(inout) :: state(0:)
    scaux=sc; scaux%J=-scaux%J;scaux%h=-scaux%h
    call IsingXgate(scaux%j/2,state)
    call MagnetickKick(scaux%h,state)
    call IsingXgate(scaux%j/2,state)
  end subroutine SymmetricEvolutionScalar
  subroutine EvolutionOperatorvJ(sc,state)
    implicit none
    type (system_carvJ), intent(in) :: sc
    complex(kind(1d0)), intent(inout) :: state(0:)
    call IsingXgate(sc%J,state)
    call MagnetickKick(sc%H,state)
  end subroutine EvolutionOperatorvJ
  subroutine EvolutionOperatorScalar(sc,state)
    implicit none
    type (system_car_scalar), intent(in) :: sc
    complex(kind(1d0)), intent(inout) :: state(0:)
    call IsingXgate(sc%J,state)
    call MagnetickKick(sc%H,state)
  end subroutine EvolutionOperatorScalar
  subroutine MagnetickKick_vector(hin,state)
    implicit none
    complex(kind(1d0)), intent(inout) :: state(0:)
    real(kind(1d0)), intent(in)       :: hin(:,0:)
    integer                           :: j
    integer                           :: qubits
    qubits=integerlog(size(state),2)
    if ( (size(hin(1,:)).ne.qubits).or. (size(hin(:,0)).ne.3)) then
       print*,"Something wrong in MagnetickKick_vector...",qubits,size(hin(1,:)),size(hin(:,0))
       stop
    end if
    do j=0,qubits-1
       call MagnetickKickj(hin(:,j),j,state)
    end do
  end subroutine MagnetickKick_vector
  subroutine MagnetickKick_scalar(hin,state)
    implicit none
    complex(kind(1d0)), intent(inout) :: state(0:)
    real(kind(1d0)), intent(in)       :: hin(3)
    integer                           :: qubits
    integer                           :: j
    qubits=integerlog(size(state),2)
    do j=0,qubits-1
       call MagnetickKickj(hin,j,state)
    end do
  end subroutine MagnetickKick_scalar
  subroutine MagnetickKickj(hin,j,state)
    implicit none
    integer, intent(in)               :: j
    real(kind(1d0)), intent(in)       :: hin(3)
    complex(kind(1d0)), intent(inout) :: state(0:)
    integer                           :: qubits
    complex(kind(1d0))                :: u(2,2),w
    real(kind(1d0))                   :: h
    integer                           :: jbit,k,ki,kf,kk,kx

    qubits=integerlog(size(state),2)
    h=norma(hin)


    if (h<10d0**(-15)) then
    else
       u(1,1)=cos(h)-I*sin(h)*hin(1)/h
       u(2,2)=conjg(u(1,1))
       u(2,1)=(hin(3)-I*hin(2))*sin(h)/h
       u(1,2)=-conjg(u(2,1))
       jbit=2**j
       do k = 0,2**(qubits-j-1) - 1
          ki = ishft(k,j+1)
          kf = ki + jbit - 1
          do kk = ki,kf
             kx = kk + jbit
             w         = u(1,1)*state(kk) + u(1,2)*state(kx)
             state(kx) = u(2,1)*state(kk) + u(2,2)*state(kx) 
             state(kk) = w
          end do
       end do
    end if
  end subroutine MagnetickKickj
  subroutine IsingXgate_scalar(alpha,state)
    implicit none
    complex(kind(1d0)), intent(inout) :: state(0:)
    real(kind(1d0)),intent(in)        :: alpha
    integer                           :: j
    integer                           :: qubits
    qubits=integerlog(size(state),2)
    do j=0,qubits-1
       call Xgatej(alpha,j,state)
    end do
  end subroutine IsingXgate_scalar
  subroutine IsingXgate_vector(alpha,state)
    implicit none
    complex(kind(1d0)), intent(inout) :: state(0:)
    real(kind(1d0)),intent(in)        :: alpha(0:)
    integer                           :: j
    integer                           :: qubits
    qubits=integerlog(size(state),2)
    if (size(alpha).ne.qubits) then
       print*,"Something wrong in IsingXgate_vector...",qubits,size(alpha)
       stop
    end if
    do j=0,qubits-1
       call Xgatej(alpha(j),j,state)
    end do
  end subroutine IsingXgate_vector
  subroutine Xgatej(alpha,j,state)
    implicit none
    complex(kind(1d0)),intent(inout) :: state(0:)
    real(kind(1d0)),intent(in)       ::  alpha
    integer,intent(in)               :: j
    integer                          ::  qubits
    complex(kind(1d0)) :: wp,wm
    integer            :: jbit01,jbit10,jbit11
    integer            :: lbit
    integer            :: k,ki,kf,kr,kk
!    print*,j
!    alpha=sc%J(j)
    qubits=integerlog(size(state),2)

    wp = cos(alpha)+I*sin(alpha)
    wm = conjg(wp)
    if ((j>=0) .and. (j <=(qubits-2))) then
       jbit01 = 2**j
       jbit10 = 2*jbit01
       jbit11 = 3*jbit01
       do k = 0,2**(qubits-j-2) - 1
          ki = ishft(k,j+2)
          kf = ki + jbit01 - 1            
          do kk = ki,kf
             state(kk) = wp*state(kk)
          end do
          do kk = ki + jbit01,kf + jbit01
             state(kk) = wm*state(kk)
          end do
          do kk = ki + jbit10,kf + jbit10
             state(kk) = wm*state(kk)
          end do
          do kk = ki + jbit11,kf + jbit11
             state(kk) = wp*state(kk)
          end do
       end do
    else
       lbit = 2**(qubits-1)
       kr = 2**(qubits-1) - 2
       do kk = 0,kr,2
          state(kk) = wp*state(kk)
       end do
       do kk = 1,kr + 1,2
          state(kk) = wm*state(kk)
       end do
       do kk = lbit,kr + lbit,2
          state(kk) = wm*state(kk)
       end do
       do kk = lbit + 1,kr + lbit + 1,2
          state(kk) = wp*state(kk)
       end do
    end if
  end subroutine Xgatej
  subroutine Ygate(beta,j,L,state)
  !implicit real*8 (a-h,o-z)      
  !parameter (maxL = 20)
  implicit none
  real (kind(1d0)), intent(in)     :: beta
  integer, intent(in)              :: j, L
  complex(kind(1d0)),intent(inout) :: state(:)
  real(kind(1d0))                  :: rstate(2,size(state))
  rstate(1,:)=real(state)
  rstate(2,:)=aimag(state)
  call Ygate_real(beta,j,L,rstate)
  state=rstate(1,:)+(1.,0.)*rstate(2,:)
  return
  end subroutine Ygate
  subroutine Ygate_real(beta,j,L,rstate)
  !implicit real*8 (a-h,o-z)      
  !parameter (maxL = 20)
  implicit none
  real (kind(1d0)), intent(in)   :: beta
  integer, intent(in)            :: j, L
  real(kind(1d0)), intent(inout) :: rstate(2,0:2**L-1)
  
  real (kind(1d0))               :: c,s, a, b
  integer                        :: k, kk, ki, kf, kx
  integer                        :: jbit01, jbit11, lbitp, lbitm, kr
  
  c = dcos(beta)
  s = dsin(beta)
  if (j .ge. 0 .and. j .le. L-2) then
     jbit01 = 2**j
     jbit11 = 3*jbit01
     do k = 0,2**(L-j-2) - 1
        ki = ishft(k,j+2)
        kf = ki + jbit01 - 1            
        do kk = ki,kf
           kx = kk + jbit11
           a = rstate(1,kk)
           b = rstate(2,kk)
           rstate(1,kk) = c*a + s*rstate(2,kx)
           rstate(2,kk) = c*b - s*rstate(1,kx)
           rstate(1,kx) = c*rstate(1,kx) + s*b
           rstate(2,kx) = c*rstate(2,kx) - s*a 
        end do
        do kk = ki + jbit01,kf + jbit01
           kx = kk + jbit01
           a = rstate(1,kk)
           b = rstate(2,kk)
           rstate(1,kk) = c*a + s*rstate(2,kx)
           rstate(2,kk) = c*b - s*rstate(1,kx)
           rstate(1,kx) = c*rstate(1,kx) + s*b
           rstate(2,kx) = c*rstate(2,kx) - s*a 
        end do
     end do
  else
     lbitp = 2**(L-1) + 1
     lbitm = 2**(L-1) - 1
     kr = 2**(L-1) - 2
     do kk = 0,kr,2
        kx = kk + lbitp
        a = rstate(1,kk)
        b = rstate(2,kk)
        rstate(1,kk) = c*a + s*rstate(2,kx)
        rstate(2,kk) = c*b - s*rstate(1,kx)
        rstate(1,kx) = c*rstate(1,kx) + s*b
        rstate(2,kx) = c*rstate(2,kx) - s*a 
     end do
     do kk = 1,kr + 1,2
        kx = kk + lbitm
        a = rstate(1,kk)
        b = rstate(2,kk)
        rstate(1,kk) = c*a + s*rstate(2,kx)
        rstate(2,kk) = c*b - s*rstate(1,kx)
        rstate(1,kx) = c*rstate(1,kx) + s*b
        rstate(2,kx) = c*rstate(2,kx) - s*a 
     end do 
  end if
  return
  end subroutine Ygate_real
  subroutine IsingGeneralX(Jising,k,kp,state)
    implicit none
    complex(kind(1d0)), intent(inout) :: state(0:)
    real(kind(1d0)), intent(in)       :: Jising
    integer,intent(in)                :: k,kp
    complex(kind(1d0))                :: eij,eijc
    integer                           :: j
    eij  = exp(I*Jising);            eijc = conjg(eij)
    do j=0,size(state)-1
       if (btest(j,k).neqv.btest(j,kp)) then
          state(j)=state(j)*eijc
       else
          state(j)=state(j)*eij
       end if
    end do
  end subroutine IsingGeneralX
  subroutine IsingGeneralY(Jising,k,kp,state)
    implicit none
    complex(kind(1d0)), intent(inout) :: state(0:)
    real(kind(1d0)), intent(in)       :: Jising
    integer,intent(in)                :: k,kp
    complex(kind(1d0))                :: R(2,2),x(2)
    integer                           :: j,zz,zo,oz,oo
    R(1,1)=cos(Jising);    R(1,2)=I*sin(Jising)
    R(2,1)=R(1,2);         R(2,2)=R(1,1)
    do j=0,size(state)/4-1
       zz=AddZeros(j,k,kp);     zo=Ibset(zz,kp)
       oz=Ibset(zz,k);          oo=Ibset(zo,k)
       x(1)=state(zz);    x(2)=state(oo)
       x=matmul(R,x)
       state(zz)=x(1);    state(oo)=x(2)
       x(1)=state(zo);    x(2)=state(oz)
       x=matmul(R,x)
       state(zo)=x(1);    state(oz)=x(2)
    end do
  end subroutine IsingGeneralY
  subroutine IsingGeneralZ(Jising,k,kp,state)
    implicit none
    complex(kind(1d0)), intent(inout) :: state(0:)
    real(kind(1d0)), intent(in)       :: Jising
    integer,intent(in)                :: k,kp
    complex(kind(1d0))                :: Rp(2,2),Rm(2,2),x(2)
    integer                           :: j,zz,zo,oz,oo
    Rp(1,1)=cos(Jising);     Rp(1,2)=I*sin(Jising)
    Rp(2,1)=Rp(1,2);         Rp(2,2)=Rp(1,1)
    Rm(1,1)=Rp(1,1);         Rm(1,2)=-Rp(1,2)
    Rp(2,1)=Rm(1,2);         Rp(2,2)=Rp(1,1)
    do j=0,size(state)/4-1
       zz=AddZeros(j,k,kp);     zo=Ibset(zz,kp)
       oz=Ibset(zz,k);          oo=Ibset(zo,k)
       x(1)=state(zz);    x(2)=state(oo)
       x=matmul(Rp,x)
       state(zz)=x(1);    state(oo)=x(2)
       x(1)=state(zo);    x(2)=state(oz)
       x=matmul(Rm,x)
       state(zo)=x(1);    state(oz)=x(2)
    end do
  end subroutine IsingGeneralZ
! Ok, tonces el operador de evolucion. El qubit especial sera el primer
! qubit. Por lo tanto tengo que tener un campo magentico para el primer
! qubit y otro para el resto. 
  subroutine ModelHeiner(Jbath,Jinteraction,phi,MagneticFieldBath, &
     MagneticFieldQubit, state)
    implicit none
    real(kind(1d0)), intent(in)  :: Jbath,Jinteraction,phi
    real(kind(1d0)), intent(in)  :: MagneticFieldBath(3), MagneticFieldQubit(3)
    complex(kind(1d0)), intent(inout) :: state(:)
    
    call FreeEvolutionRingButOne(MagneticFieldBath, state)
    call MagnetickKickj( MagneticFieldQubit,0,state)
    call InternalInteractionRingButOneX(Jbath, state)
    call InteractionRingWithFirstQubitX(Jinteraction, state)
    call FreeEvolutionRingButOne(MagneticFieldBath, state)
    call MagnetickKickj( MagneticFieldQubit,0,state)
    call InternalInteractionRingButOneX(Jbath, state)

    call InteractionRingWithFirstQubitX(Jinteraction*cos(phi), state)
    call InteractionRingWithFirstQubitY(Jinteraction*sin(phi), state)


  end subroutine ModelHeiner
  subroutine InternalInteractionRingButOneX(J, state)
    implicit none
    real(kind(1d0)), intent(in)       :: J
!    real(kind(1d0)), intent(in)       :: MagneticFieldFirstQubit(3)
    complex(kind(1d0)), intent(inout) :: state(0:)
    integer                           :: i, qubits
    qubits=integerlog(size(state),2)
    do i=1,qubits-2
      call Xgatej( J, i, state)
    enddo
    call IsingGeneralX(J,1,qubits-1,state)
  end subroutine InternalInteractionRingButOneX
  subroutine InteractionRingWithFirstQubitY(J, state)
    implicit none
    real(kind(1d0)), intent(in)       :: J
!    real(kind(1d0)), intent(in)       :: MagneticFieldFirstQubit(3)
    complex(kind(1d0)), intent(inout) :: state(0:)
    integer                           :: i, qubits
    qubits=integerlog(size(state),2)
    do i=1,qubits-1
      call IsingGeneralY(J,0,i,state)
    enddo
  end subroutine InteractionRingWithFirstQubitY
  subroutine InteractionRingWithFirstQubitX(J, state)
    implicit none
    real(kind(1d0)), intent(in)       :: J
!    real(kind(1d0)), intent(in)       :: MagneticFieldFirstQubit(3)
    complex(kind(1d0)), intent(inout) :: state(0:)
    integer                           :: i, qubits
    qubits=integerlog(size(state),2)
    do i=1,qubits-1
      call IsingGeneralX(J,0,i,state)
    enddo
  end subroutine InteractionRingWithFirstQubitX
  subroutine FreeEvolutionRingButOne(Magneticfield, state)
    implicit none
    real(kind(1d0)), intent(in)       :: Magneticfield(3)
    complex(kind(1d0)), intent(inout) :: state(0:)
    integer                           :: i, qubits
    qubits=integerlog(size(state),2)
    !call MagnetickKickj( MagneticFieldFirstQubit,0,state)
    do i=1,qubits-1
      call MagnetickKickj( MagneticField,i,state)
    enddo
  end subroutine FreeEvolutionRingButOne
end module EvolutionOperators
