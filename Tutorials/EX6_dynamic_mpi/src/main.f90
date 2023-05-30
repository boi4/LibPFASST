!
!

!>
program main
  use pf_mod_mpi
  use pf_mod_zndarray


  integer ::  ierr
  integer :: session

  !> Initialize MPI Session
  call mpi_session_init(MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, session, ierr)
  if (ierr /= 0) &
      stop "ERROR: Can't initialize MPI."

  !> Call the  solver
  call run_pfasst(session)

  !> Close mpi
  call mpi_session_finalize(session, ierr)

contains

  !>  This subroutine implements pfasst to solve the advection diffusion equation
  subroutine run_pfasst(session)
    use pfasst  !< This module has include statements for the main pfasst routines
    use my_sweeper  !< Local module for sweeper and function evaluations
    use my_level    !< Local module for the levels
    use hooks   !< Local module for diagnostics and i/o
    use probin  !< Local module reading/parsing problem parameters

    implicit none

    !> argument
    integer, intent(in) :: session

    !>  Local variables
    type(pf_pfasst_t)              :: pf       !<  the main pfasst structure
    type(pf_dynprocs_t)            :: dynprocs     !<  instead of a pf_comm_t object, we use a pf_dynprocs_t object!
    type(pf_zndarray_t)            :: y_0      !<  the initial condition
    character(256)                 :: pf_fname   !<  file name for input of PFASST parameters
    logical                        :: is_dynamic
    logical                        :: premature_exit

    integer                        ::  l   !  loop variable over levels

    type(pf_zndarray_t)            :: qend      !<  the solution
    complex(pfdp),         pointer :: sol(:)




    !> Read problem parameters
    call probin_init(pf_fname)

    ! determine if dynamic
    call pf_dynprocs_check_dynamic(session, is_dynamic)

    !>  Set up communicator
    call pf_dynprocs_create(dynprocs, session, "mpi://WORLD")

    !>  Create the pfasst structure
    call pf_pfasst_create_dynamic(pf, dynprocs, fname=pf_fname)


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
    call pf_pfasst_setup(pf)

    !> Add some hooks for output
    call pf_add_hook(pf, -1, PF_POST_SWEEP, echo_error)

    !> Add hook for resizing logic
    call pf_add_hook(pf, -1, PF_PRE_POT_RESIZE, resize_decider)

    !>  Output the run options
    call pf_print_options(pf,un_opt=6)

    !>  Output local parameters
    call print_loc_options(pf,un_opt=6)
    
    !>  Allocate initial consdition
    call zndarray_build(y_0, [ 1 ])

    !> Set the initial condition
    call y_0%setval(1.0_pfdp)

    !>  Allocate solution
    call zndarray_build(qend, [ 0 ])

    !> Do the PFASST stepping
    call pf_pfasst_run(pf, y_0, dt, 0.0_pfdp, nsteps, qend, join_existing=is_dynamic, premature_exit=premature_exit)

    if (.not. premature_exit) then
       sol  => get_array1d(qend)
       print *, sol
       !call qend%pack(sol)
       !print *,sol

       !>  Wait for everyone to be done
       call mpi_barrier(pf%comm%comm, ierr)

    end if

    !>  Deallocate initial condition and final solution
    call zndarray_destroy(y_0)

    !>  Deallocate pfasst structure
    call pf_pfasst_destroy(pf)

    !> free PFASST communicator (already)
    call mpi_comm_disconnect(pf%comm%comm, ierr)
  end subroutine run_pfasst

end program
