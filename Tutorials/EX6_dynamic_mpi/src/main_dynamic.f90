module myapp

use pf_mod_mpi
use mpi
use pf_mod_zndarray
use pf_mod_dynres
use, intrinsic :: iso_c_binding

contains
  subroutine entry() bind(c)
    implicit none

    integer :: ierror
    integer :: world_rank
    logical :: finished
    logical :: initial_start
    integer :: size

    !> Initialize MPI
    call mpi_init(ierror)
    if (ierror /= 0) &
        stop "ERROR: Can't initialize MPI."

    call mpi_comm_rank(MPI_COMM_WORLD, world_rank, ierror)


    initial_start = .TRUE.
    do while (.not. finished)
      call run(initial_start, finished)
      initial_start = .FALSE.

      do while (.TRUE.)
         ! receive size of new communicator
         call mpi_bcast(size, 1, MPI_INT, 0, MPI_COMM_WORLD, ierror)
         if (world_rank < size) then
            exit
         else
            ! dummy
            call mpi_bcast(size, 1, MPI_INT, 0, MPI_COMM_WORLD, ierror)
         end if
      end do
    end do

    !> Close mpi
    call mpi_finalize(ierror)
  end subroutine entry


  !>  This subroutine implements pfasst to solve the advection diffusion equation
  subroutine run(initial_start, finished)
    use pfasst  !< This module has include statements for the main pfasst routines
    use my_sweeper  !< Local module for sweeper and function evaluations
    use my_level    !< Local module for the levels
    use hooks   !< Local module for diagnostics and i/o
    use probin  !< Local module reading/parsing problem parameters


    implicit none

    logical,intent(in)  :: initial_start
    logical,intent(out) :: finished

    !>  Local variables
    type(pf_pfasst_t)  :: pf       !<  the main pfasst structure
    type(pf_comm_t)    :: comm     !<  the communicator
    type(pf_dynres_t)  :: dynres     !<  the dynamic resource handler
    type(pf_zndarray_t):: y_0      !<  the initial condition
    character(256)     :: pf_fname   !<  file name for input of PFASST parameters

    integer           :: l   !  loop variable over levels
    integer           :: ierror

    type(pf_zndarray_t):: qend      !<  the solution
    complex(pfdp),         pointer :: sol(:)




    !> Read problem parameters
    call probin_init(pf_fname)

    !>  Create the pfasst structure
    call pf_pfasst_create_dynamic(pf, dynres=dynres, fname=pf_fname)

    !> Set up dynres handler
    call pf_dynres_create(dynres)
    dynres%is_dynamic_start = .not. initial_start

    !> Create communicator by either joining existing run or by mpi://WORLD
    call pf_dynres_create_comm(dynres, comm)


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

    !>  Output the run options
    if (.not. pf%dynres%is_dynamic_start) then
       call pf_print_options(pf,un_opt=6)
    end if

    !>  Output local parameters
    call print_loc_options(pf,un_opt=6)

    !>  Allocate initial consdition
    call zndarray_build(y_0, [ 1 ])

    !> Set the initial condition
    call y_0%setval(1.0_pfdp)

    !>  Allocate solution
    call zndarray_build(qend, [ 0 ])

    call pf_pfasst_run(pf, y_0, dt, 0.0_pfdp, nsteps, qend, join_existing=pf%dynres%is_dynamic_start)

    finished = .not. pf%dynres%needs_shutdown

    if (finished .and. pf%rank == (pf%comm%nproc - 1)) then
       sol  => get_array1d(qend)
       print *, "Finished!"
       print *, sol
       !call qend%pack(sol)
       !print *,sol
    end if

    !> destroy dynres handler
    call pf_dynres_destroy(dynres)

    !>  Deallocate initial condition and final solution
    call zndarray_destroy(y_0)

    !>  Deallocate pfasst structure
    call pf_pfasst_destroy(pf)

  end subroutine run


end module myapp


program main
  use myapp, only: entry

  call entry()
end program main
