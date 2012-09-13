module varios
  !  use definitions
  !  use EvolutionOperators
  !use bit_manipulation
  implicit none
contains
  subroutine prepare_IC(qubits_cen,state)
    use ask
    use math, only : integerlogbase2
    use quantumwau,only :Merge_States,state_preparation
    integer, intent(in)             :: qubits_cen
    complex(kind(1d0)), intent(out) :: state(:)
    integer                         :: qubits_env,qubits,Type_env,Type_cen,n_o_env,n_o_cen,j1
    integer, allocatable            :: options_cen(:),options_env(:)
    complex(kind(1d0)), allocatable :: state_cen(:),state_env(:)

    Type_cen=16; n_o_cen=0; Type_env   =18 ; n_o_env=2;
    print*,"Inserte las primeras opciones para la condicion inicial:",&
        "Type_central_system,np_cen,Type_environment,np_env"
    if (debe_preguntar) read(*,*)Type_cen,n_o_cen,Type_env,n_o_env
    qubits=IntegerLogBase2(size(state)); qubits_env=qubits-qubits_cen
    allocate(options_cen(n_o_cen),options_env(n_o_env)); options_cen=3;options_env=0;
    if (n_o_env==2) options_env=(/-498234,923842398/)
    allocate(state_cen(2**qubits_cen),state_env(2**qubits_env))
    print*,"Y ahora las ",n_o_cen,"opciones para el estado central";
    if (debe_preguntar) read(*,*)(options_cen(j1),j1=1,n_o_cen)
    print*,"Y ahora las ",n_o_env,"opciones para el estado environment";
    if (debe_preguntar) read(*,*)(options_env(j1),j1=1,n_o_env)
    call state_preparation(state_env, Type_env,options_env);
    call state_preparation(state_cen, Type_cen,options_cen); 
    call Merge_States(state_cen,state_env,2**qubits_cen-1, state);
    deallocate(state_cen,state_env,options_cen,options_env)
  end subroutine prepare_IC
end module varios
