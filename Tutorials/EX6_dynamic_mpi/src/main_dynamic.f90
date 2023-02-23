!
!

!>
module myapp

use pf_mod_mpi
use pf_mod_zndarray
use pf_mod_dynres
use mpidynres
use, intrinsic :: iso_c_binding

contains
  subroutine entry() bind(c)
    implicit none

    integer ::  ierror
    type(c_ptr) ::  session


    !> Initialize MPI Session
    call mpi_session_init(MPI_INFO_NULL, MPI_ERRHANDLER_NULL, session, ierror)
    if (ierror /= 0) &
        stop "ERROR: Can't initialize MPI."

    !> Call the solver
    call run(session)

    !> Close mpi
    call mpi_session_finalize(session)
  end subroutine entry


  !>  This subroutine implements pfasst to solve the equation
  subroutine run(session)
    use pfasst  !< This module has include statements for the main pfasst routines
    use my_sweeper  !< Local module for sweeper and function evaluations
    use my_level    !< Local module for the levels
    use hooks   !< Local module for diagnostics and i/o
    use probin  !< Local module reading/parsing problem parameters


    implicit none

    type(c_ptr), intent(in), target ::  session

    !>  Local variables
    type(pf_pfasst_t) :: pf       !<  the main pfasst structure
    type(pf_comm_t)   :: comm     !<  the communicator
    type(pf_dynres_t)   :: dynres     !<  the dynamic resource handler
    type(pf_zndarray_t):: y_0      !<  the initial condition
    character(256)    :: pf_fname   !<  file name for input of PFASST parameters
    logical           :: finished

    integer           :: l   !  loop variable over levels
    integer           :: ierror

    type(pf_zndarray_t):: qend      !<  the solution
    complex(pfdp),         pointer :: sol(:)




    !> Read problem parameters
    call probin_init(pf_fname)

    !>  Create the pfasst structure
    call pf_pfasst_create_dynamic(pf, dynres=dynres, fname=pf_fname)

    if (pf%use_dynamic_mpi) then
        !> Set up dynres handler
        call pf_dynres_create(dynres, session)

        !> Create communicator by either joining existing run or by mpi://WORLD
        call pf_dynres_create_comm(dynres, comm)
    else
        !>  Set up communicator
        call pf_mpi_create(comm, MPI_COMM_WORLD)
    end if


    !> Loop over levels and set some level specific parameters
    do l = 1, pf%nlevels
       !>  Allocate the user specific level object
       allocate(my_level_t::pf%levels(l)%ulevel)

       !>  Allocate the user specific data factory
       allocate(pf_zndarray_factory_t::pf%levels(l)%ulevel%factory)

       !>  Add the sweeper to the level
       allocate(my_sweeper_t::pf%levels(l)%ulevel%sweeper)
       ! my_sweeper_t::pf%levels(l)%ulevel%sweeper.implicit = .TRUE. ! ONLY EXPLICIT PART

       !>  Set the size of the data on this level
       ! buflen must be set explicitly as complex numbers need double the bufsize
       call pf_level_set_size(pf,l,[1], 2)
    end do

    !>  Set up some pfasst stuff
    call pf_pfasst_setup_dynamic(pf, comm)

    !> Add some hooks for output
    call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)

    if (.not. pf%dynres%is_dynamic_start) then
       !>  Output the run options
       call pf_print_options(pf,un_opt=6)
       !>  Output local parameters
       call print_loc_options(pf,un_opt=6)
    end if


    !>  Allocate initial consdition
    call zndarray_build(y_0, [ 1 ])

    !> Set the initial condition
    call y_0%setval(1.0_pfdp)

    !>  Allocate solution
    call zndarray_build(qend, [ 0 ])

    !> Do the PFASST stepping
    if (pf%use_dynamic_mpi) then
       call pf_pfasst_run(pf, y_0, dt, 0.0_pfdp, nsteps, qend, join_existing=pf%dynres%is_dynamic_start)
    else
       call pf_pfasst_run(pf, y_0, dt, 0.0_pfdp, nsteps, qend)
    end if

    finished = .TRUE.
    if (pf%use_dynamic_mpi) then
       finished = .not. pf%dynres%needs_shutdown
    end if

    if (finished .and. pf%state%step == (pf%state%nsteps - 1)) then
       sol  => get_array1d(qend)
       print *, "Finished!"
       print *, sol
       !call qend%pack(sol)
       !print *,sol
    end if

    if (pf%use_dynamic_mpi) then
        !> destroy dynres handler
        call pf_dynres_destroy(dynres)
    end if

    !>  Deallocate initial condition and final solution
    call zndarray_destroy(y_0)

    !>  Deallocate pfasst structure
    call pf_pfasst_destroy(pf)

  end subroutine run











  subroutine seed_rng()
    integer :: n
    integer,allocatable :: seed(:)

    call random_seed(size=n)
    allocate(seed(n))
    seed = 123456789    ! putting arbitrary seed to all elements
    call random_seed(put=seed)
    deallocate(seed)
  end subroutine seed_rng


  character(len=64) function itoa(n) result(s)
      integer, intent(in) :: n
      write (s, *) n
      s = adjustl(s)
  end function itoa


end module myapp


program main
  use myapp, only: entry,itoa

  use pf_mod_mpi
  use mpidynres_sim

  implicit none

  type(MPIDYNRES_SIM_CONFIG) :: config
  integer :: ierror, world_size
  procedure(mpidynres_main_func), pointer :: main_func => entry

  call MPI_INIT(ierror)
  !call MPIDYNRES_SIM_GET_DEFAULT_CONFIG(config)
  config%base_communicator = MPI_COMM_WORLD
  call MPI_INFO_CREATE(config%manager_config, ierror)
  ! call MPI_INFO_SET(config%manager_config, "manager_initial_number_random", "yes")
  call MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, ierror)
  call MPI_INFO_SET(config%manager_config, "manager_initial_number", itoa(world_size - 1), ierror)

  call MPIDYNRES_SIM_START(config, main_func)

  call MPI_INFO_FREE(config%manager_config, ierror)
  call MPI_FINALIZE(ierror)
end program main
