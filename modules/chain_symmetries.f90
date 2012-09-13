module QubitChainSymmetries
  use math
  implicit none
  type base_generator_real
    integer     :: base_number
    logical     :: degenerate
    logical     :: plus_sign
  end type base_generator_real
contains
  function size_all_k(qubits)
    integer, intent(in) :: qubits
    integer             :: size_all_k(0:qubits-1)
    integer             :: j,j1,k
    logical,allocatable :: consider(:)
    allocate(consider(0:2**qubits-1))
    consider=.true.
    size_all_k=0
    do j=0,2**qubits-1
       if (consider(j)) then
          do j1=0,qubits-1; consider(T_number(j,j1,qubits))=.false.; end do
          do k=0,qubits-1;  If (.not.ProyectToZero(k,SizeCycle(j,qubits))) size_all_k(k)=size_all_k(k)+1; enddo
       end if
    end do
    deallocate(consider)
  end function size_all_k
  function size_k_space(k,qubits)
    integer             :: size_k_space
    integer, intent(in) :: k,qubits
    integer             :: j,j1
    logical,allocatable :: consider(:)
    allocate(consider(0:2**qubits-1))
    consider=.true.
    size_k_space=0
    do j=0,2**qubits-1
       if (consider(j)) then
          do j1=0,qubits-1; 
             consider(T_number(j,j1,qubits))=.false.; 
          end do
          If (.not.ProyectToZero(k,SizeCycle(j,qubits))) size_k_space=size_k_space+1; 
       end if
    end do
    deallocate(consider)
  end function size_k_space
  function SizeCycle(n,qubits)
    implicit none
    integer                :: SizeCycle
    integer,intent(in)     :: n,qubits
    integer                :: ntmp,nmod
    nmod=mod(n,2**qubits);     ntmp=T_number(nmod,1,qubits);    SizeCycle=1
    do while (ntmp.ne.nmod)
       ntmp=T_number(ntmp,1,qubits);       SizeCycle=SizeCycle+1
    end do
    SizeCycle=qubits/SizeCycle
  end function SizeCycle
  subroutine Get_cycles_n(n,n_rotated,qubits)
    integer, intent(in)       :: n,qubits
    integer, intent(out)      :: n_rotated(:)
    integer                   :: j,sc
    sc= qubits/SizeCycle(n,qubits)
    if (size(n_rotated).ne. SC)then
       print*,"algun pedo en Get_cycles_n",size(n_rotated),SizeCycle(n,qubits),n,qubits; stop
    end if
    do J=0,sc-1; n_rotated(J+1)=T_number(n,j,qubits); end do
  end subroutine Get_cycles_n
  function ProyectToZero(k,SCycle)
    logical             :: ProyectToZero
    integer, intent(in) :: k,SCycle
    ProyectToZero=(modulo(k,SCycle).ne.0)
  end function ProyectToZero
  function P_Number(number,qubits)
    integer               :: P_number,j
    integer, intent(in)   :: number,qubits
    P_Number=0
    do j=0,qubits-1; if(btest(number,j)) P_number=ibset(P_Number,qubits-j-1); end do
  end function P_Number
  function P_state(state)
    implicit none
    complex(kind(1d0)), intent(in)    :: state(0:)
    complex(kind(1d0))                :: P_state(0:size(state)-1)
    integer                           :: qubits, j1
    qubits=IntegerLogBase2(size(state))
    do j1=0,2**qubits-1
       P_state(   P_number(j1,qubits)         )=state(j1)
    enddo
  end function P_state
  function T_number(number,power,qubits)
    implicit none
    integer               :: T_number
    integer, intent(in)   :: number,power,qubits
    T_number=ishftc(number,power,qubits)
  end function T_number
  function anti_unitary_keep_k(state)
    implicit none
    complex(kind(1d0)), intent(in)    :: state(0:)
    complex(kind(1d0))                :: anti_unitary_keep_k(0:size(state)-1)
    anti_unitary_keep_k=P_state(conjg(state))
  end function anti_unitary_keep_k
  function Proyector_k(state,k)
    implicit none
    complex(kind(1d0)), intent(in)    :: state(0:)
    complex(kind(1d0))                :: Proyector_k(0:size(state)-1)
    integer, intent(in)               :: k
    integer                           :: qubits, j1
    qubits=IntegerLogBase2(size(state))
    Proyector_k=0d0
    do j1=0,qubits-1
       Proyector_k=Proyector_k + exp(-2*I*pi*j1*k/real(qubits,kind(1d0)))*T_state(state,j1)/qubits
    enddo
  end function Proyector_k
  function T_state(state,power)
    implicit none
    complex(kind(1d0)), intent(in)    :: state(0:)
    complex(kind(1d0))                :: T_state(0:size(state)-1)
    integer, intent(in)               :: power
    integer                           :: qubits, j1
    qubits=IntegerLogBase2(size(state))
    do j1=0,2**qubits-1
       T_state(   T_number(j1,power,qubits)         )=state(j1)
    enddo
  end function T_state
  subroutine safe_T_state(state_in,power,state_out)
    implicit none
    complex(kind(1d0)), intent(in)    :: state_in(0:)
    complex(kind(1d0)), intent(out)   :: state_out(0:)
    integer, intent(in)               :: power
    integer                           :: qubits, j1
    qubits=IntegerLogBase2(size(state_in))
    do j1=0,2**qubits-1
       state_out(   T_number(j1,power,qubits)         )=state_in(j1)
    enddo
  end subroutine safe_T_state
  subroutine safe_Proyector_k(state,k)
    implicit none
    complex(kind(1d0)), intent(inout) :: state(0:)
    complex(kind(1d0)),allocatable    :: Proyector_k(:),state_tmp(:)
    integer, intent(in)               :: k
    integer                           :: qubits, j1
    qubits=IntegerLogBase2(size(state))
    allocate(Proyector_k(0:2**qubits-1),state_tmp(0:2**qubits-1))
    Proyector_k=0d0
    do j1=0,qubits-1
       call safe_T_state(state,j1,state_tmp)
       Proyector_k=Proyector_k + exp(-2*I*pi*j1*k/real(qubits,kind(1d0)))*state_tmp/qubits
    enddo
    state=Proyector_k
    deallocate(Proyector_k,state_tmp)
  end subroutine safe_Proyector_k
  subroutine encode_a_good_base(k,qubits,base_generators)
    integer, intent(in)                      :: k,qubits
    type(base_generator_real), intent(out)   :: base_generators(:)
    logical, allocatable                     :: consider(:)
    integer                                  :: n,size_cycle,j_1,counter_dim,degeneration_period(1)
    integer,allocatable                      :: n_rotated_s(:),n_reflected_rotated_s(:)
    allocate(consider(0:2**qubits-1)); consider=.true.; counter_dim=0
    do n=0,2**qubits-1
      If (.not.consider(n)) cycle
      size_cycle=qubits/SizeCycle(n,qubits); allocate(n_rotated_s(size_cycle)); allocate(n_reflected_rotated_s(size_cycle))
      !call Get_cycles_n(n,n_rotated_s,qubits);call Get_cycles_n(P_number(n,qubits),n_reflected_rotated_s,qubits)
      call Get_cycles_n(n,n_rotated_s,qubits);
      do j_1=1,size_cycle; n_reflected_rotated_s(j_1)=P_number(n_rotated_s(j_1),qubits); enddo
      If (ProyectToZero(k,SizeCycle(n,qubits))) then
        do j_1=1,size_cycle; 
          consider(n_rotated_s(j_1))=.false.;
          consider(n_reflected_rotated_s(j_1))=.false.;
          end do
        deallocate(n_rotated_s,n_reflected_rotated_s); cycle
      end If
      if (any(n==n_reflected_rotated_s))then
        counter_dim=counter_dim+1
        base_generators(counter_dim)%degenerate  = .true.
        base_generators(counter_dim)%base_number = n
        degeneration_period=minloc(abs(n_reflected_rotated_s-n))
        !print*,n,"es degenerado y tiene periodo",degeneration_period(1),qubits,"=qubits"
      !!! Aca voy en la traduccion
        if (mod(2*(degeneration_period(1)-1)*k,qubits)==0.and.mod((degeneration_period(1)-1)*k,qubits).ne.0) then
          !if (mod(2*degeneration_period(1)*k,qubits)==0.and.mod(degeneration_period(1)*k,qubits).ne.0) then
          !if (mod(qubits,2)==0.and.mod(degeneration_period(1)*k,qubits/2)==0) then
          base_generators(counter_dim)%plus_sign=.false.
          else; 
           base_generators(counter_dim)%plus_sign=.true.
        end if
        do j_1=2,size_cycle;consider(n_rotated_s(j_1))=.false.;enddo
      else
        counter_dim=counter_dim+2
        base_generators(counter_dim-1:counter_dim)%degenerate  = .false.
        base_generators(counter_dim-1:counter_dim)%base_number = (/n,n/)
        base_generators(counter_dim-1:counter_dim)%plus_sign = (/.true.,.false./)
        do j_1=2,size_cycle;consider(n_rotated_s(j_1))=.false.;enddo
        do j_1=1,size_cycle;consider(n_reflected_rotated_s(j_1))=.false.;enddo
      end if
      deallocate(n_rotated_s,n_reflected_rotated_s)
      !stop
    end do
    if(counter_dim.ne.size(base_generators)) then
      print*,"eeeeeeror",k,qubits,counter_dim,size(base_generators);    stop
    end if
    deallocate(consider)
  end subroutine encode_a_good_base
  function decode_base_element(base_generator,qubits,k)
    type(base_generator_real), intent(in)   :: base_generator
    integer, intent(in)                     :: qubits,k
    !complex(kind(1d0))                      :: decode_base_element(0:)
    complex(kind(1d0))                      :: decode_base_element(0:2**qubits-1)
    integer, allocatable                    :: n_r(:),pn_r(:),positions_pn_r(:)
    complex(kind(1d0)), allocatable         :: exponentials(:)
    integer                                 :: size_cycle,j_1
    real(kind(1d0))                         :: x
    decode_base_element=0d0; 
    size_cycle=qubits/SizeCycle(base_generator%base_number,qubits)
    allocate(n_r(size_cycle),pn_r(size_cycle),positions_pn_r(size_cycle),exponentials(size_cycle)); 
    call Get_cycles_n(base_generator%base_number,n_r,qubits)
    do J_1=1,size_cycle;
       exponentials(j_1)=exp(-(2*pi*I*(j_1-1)*k)/qubits)
       pn_r(j_1)=P_number(n_r(j_1),qubits)
    enddo
    if (base_generator%degenerate) then ! si son degenerados
       do J_1=1,size_cycle; positions_pn_r(j_1:j_1)=minloc(abs(n_r-pn_r(j_1))); enddo;
       if (base_generator%plus_sign) then
          x=1/(abs(exponentials(1)+conjg(exponentials(positions_pn_r(1))))*sqrt(real(size_cycle,kind(1d0))))
          do J_1=1,size_cycle;decode_base_element(n_r(j_1))=x*(exponentials(j_1)+conjg(exponentials(positions_pn_r(j_1))));enddo
       else
          x=1/(abs(exponentials(1)-conjg(exponentials(positions_pn_r(1))))*sqrt(real(size_cycle,kind(1d0))))
          do J_1=1,size_cycle;decode_base_element(n_r(j_1))=x*I*(exponentials(j_1)-conjg(exponentials(positions_pn_r(j_1))));enddo
       end if
    else ! si no son degenerados
       x=1/sqrt(2d0*size_cycle)
       if (base_generator%plus_sign) then
          do J_1=1,size_cycle
             decode_base_element(n_r(j_1)) = exponentials(j_1)*x; decode_base_element(pn_r(j_1)) = conjg(exponentials(j_1))*x
          enddo
       else
          do J_1=1,size_cycle
             decode_base_element(n_r(j_1)) = I*exponentials(j_1)*x;decode_base_element(pn_r(j_1)) = -I*conjg(exponentials(j_1))*x
          enddo
       end if
    end if
    deallocate(n_r,exponentials,positions_pn_r,pn_r)
  end function decode_base_element
  subroutine safe_decode_base_element(base_generator,qubits,k,decode_base_element)
    type(base_generator_real), intent(in)   :: base_generator
    integer, intent(in)                     :: qubits,k
    complex(kind(1d0)),intent(out)          :: decode_base_element(0:)
    integer, allocatable                    :: n_r(:),pn_r(:),positions_pn_r(:)
    complex(kind(1d0)), allocatable         :: exponentials(:)
    integer                                 :: size_cycle,j_1
    real(kind(1d0))                         :: x
    decode_base_element=0d0; 
    size_cycle=qubits/SizeCycle(base_generator%base_number,qubits)
    allocate(n_r(size_cycle),pn_r(size_cycle),positions_pn_r(size_cycle),exponentials(size_cycle)); 
    call Get_cycles_n(base_generator%base_number,n_r,qubits)
    do J_1=1,size_cycle;
       exponentials(j_1)=exp(-(2*pi*I*(j_1-1)*k)/qubits)
       pn_r(j_1)=P_number(n_r(j_1),qubits)
    enddo
    if (base_generator%degenerate) then ! si son degenerados
       do J_1=1,size_cycle; positions_pn_r(j_1:j_1)=minloc(abs(n_r-pn_r(j_1))); enddo;
       if (base_generator%plus_sign) then
          x=1/(abs(exponentials(1)+conjg(exponentials(positions_pn_r(1))))*sqrt(real(size_cycle,kind(1d0))))
          do J_1=1,size_cycle;decode_base_element(n_r(j_1))=x*(exponentials(j_1)+conjg(exponentials(positions_pn_r(j_1))));enddo
       else
          x=1/(abs(exponentials(1)-conjg(exponentials(positions_pn_r(1))))*sqrt(real(size_cycle,kind(1d0))))
          do J_1=1,size_cycle;decode_base_element(n_r(j_1))=x*I*(exponentials(j_1)-conjg(exponentials(positions_pn_r(j_1))));enddo
       end if
    else ! si no son degenerados
       x=1/sqrt(2d0*size_cycle)
       if (base_generator%plus_sign) then
          do J_1=1,size_cycle
             decode_base_element(n_r(j_1)) = exponentials(j_1)*x; decode_base_element(pn_r(j_1)) = conjg(exponentials(j_1))*x
          enddo
       else
          do J_1=1,size_cycle
             decode_base_element(n_r(j_1)) = I*exponentials(j_1)*x;decode_base_element(pn_r(j_1)) = -I*conjg(exponentials(j_1))*x
          enddo
       end if
    end if
    deallocate(n_r,exponentials,positions_pn_r,pn_r)
  end subroutine safe_decode_base_element
  function decode_base_element_o(base_generator,qubits,k)
    type(base_generator_real), intent(in)   :: base_generator
    integer, intent(in)                     :: qubits,k
    complex(kind(1d0))                      :: decode_base_element_o(0:2**qubits-1)
    integer, allocatable                    :: n_r(:)
    integer                                 :: size_cycle,j_1
    decode_base_element_o=0d0; size_cycle=qubits/SizeCycle(base_generator%base_number,qubits)
    allocate(n_r(size_cycle)); call Get_cycles_n(base_generator%base_number,n_r,qubits)
    do J_1=1,size_cycle
       decode_base_element_o(n_r(j_1))=exp(-(2*pi*I*(j_1-1)*k)/qubits)
       if (.not.base_generator%plus_sign)decode_base_element_o(n_r(j_1))=I*decode_base_element_o(n_r(j_1))
    end do
    do J_1=1,size_cycle
       if (base_generator%plus_sign) then
          decode_base_element_o(P_number(n_r(j_1),qubits)) =&
               decode_base_element_o(P_number(n_r(j_1),qubits))+exp((2*pi*I*(j_1-1)*k)/qubits)
       else
          decode_base_element_o(P_number(n_r(j_1),qubits)) =&
               decode_base_element_o(P_number(n_r(j_1),qubits))-I*exp((2*pi*I*(j_1-1)*k)/qubits)

       end if
    end do
    decode_base_element_o=decode_base_element_o/norma(decode_base_element_o)
    deallocate(n_r)
  end function decode_base_element_o
  function maximum_relevant_k(q)
    integer            :: maximum_relevant_k
    integer,intent(in) :: q
    if (mod(q,2)==0) maximum_relevant_k=q/2-1
    if (mod(q,2)==1) maximum_relevant_k=(q-1)/2
  end function maximum_relevant_k
  subroutine random_k_state(state,k)
     complex(kind(1d0)),intent(inout)  :: state(:)
     integer, intent(in)               :: k
     integer                           :: j
     do j=1,size(state)
        state(j)=RandomGaussian()
     enddo
     call safe_Proyector_k(state,k)
     state=state/norma(state)
  end subroutine random_k_state

end module QubitChainSymmetries
