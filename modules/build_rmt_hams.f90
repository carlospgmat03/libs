module build_rmt_hams
  use GUE_generator,   only : generate_GUE_unfolded_normal_matrix, generate_GUE_matrix
  use math,            only : IntegerLogBase2
  use mylinearalgebra, only : Add_operator_some_qubits
  implicit none
  interface  build_gue_spectator_ham
     module procedure build_gue_spectator_ham_first
     module procedure build_gue_spectator_ham_n
  end interface
contains 
  subroutine build_gue_spectator_ham_first(Hamiltonian,delta)
    complex(kind(1d0)), intent(out)      :: Hamiltonian(:,:)
    real(kind(1d0)), intent(in)          :: delta
    integer                              :: size_h,qubits
     complex(kind(1d0)), allocatable     :: Henv(:,:),H_1env(:,:)
    Hamiltonian=0d0
    size_h=size(Hamiltonian,1)
    qubits=IntegerLogBase2(size_h)
    if (size_h.ne.size(Hamiltonian,2)) stop "error en build_gue_spectator_ham"
    if (qubits.le.2) stop "error en build_gue_spectator_ham 2"
    allocate(Henv(2**(qubits-2),2**(qubits-2)),H_1env(2**(qubits-1),2**(qubits-1)))
    call generate_GUE_unfolded_normal_matrix(Henv); call Add_operator_some_qubits(Hamiltonian,Henv,3)
    call generate_GUE_matrix(H_1env); H_1env=delta*H_1env; call Add_operator_some_qubits(Hamiltonian,H_1env,1)
    deallocate(Henv,H_1env)
  end subroutine build_gue_spectator_ham_first
! 
!   subroutine build_gue_spectator_ham_n(Hamiltonian,delta,qubit)
!     complex(kind(1d0)), intent(out)      :: Hamiltonian(:,:)
!     real(kind(1d0)), intent(in)          :: delta
!     integer, intent(in)                  :: qubit
!     integer                              :: size_h,qubits
!      complex(kind(1d0)), allocatable     :: Henv(:,:),H_1env(:,:)
!     Hamiltonian=0d0
!     size_h=size(Hamiltonian,1)
!     qubits=IntegerLogBase2(size_h)
!     if (size_h.ne.size(Hamiltonian,2)) stop "error en build_gue_spectator_ham"
!     if (qubits.le.2) stop "error en build_gue_spectator_ham 2"
!     allocate(Henv(2**(qubits-2),2**(qubits-2)),H_1env(2**(qubits-1),2**(qubits-1)))
!     call generate_GUE_unfolded_normal_matrix(Henv); call Add_operator_some_qubits(Hamiltonian,Henv,3)
!     call generate_GUE_matrix(H_1env); H_1env=delta*H_1env; call Add_operator_some_qubits(Hamiltonian,H_1env,1)
!     deallocate(Henv,H_1env)
!   end subroutine build_gue_spectator_ham_n
end module
