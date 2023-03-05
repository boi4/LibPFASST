!!  Implementation of Dynamic MPI support
!
! This file is part of LIBPFASST.
!
! TODO: add comments
module pf_mod_dynres
  use pf_mod_pfasst
  use mpi, only: mpi_info_get, mpi_info_free, MPI_ERRHANDLER_NULL, MPI_INFO_NULL, MPI_INT, &
                 mpi_allreduce, MPI_SUM, MPI_CHAR, MPI_COMM_WORLD, MPI_GROUP_NULL

  implicit none
contains


  ! create new dynres object using initialized MPI session
  subroutine pf_dynres_create(this)
    type(pf_dynres_t), intent(out)   :: this

    this%needs_shutdown = .FALSE.
  end subroutine pf_dynres_create

  ! destroy dynres object
  ! Empty for now, as no dynamic allocation
  subroutine pf_dynres_destroy(this)
    type(pf_dynres_t), intent(out)   :: this
  end subroutine pf_dynres_destroy


  ! create new communicator from mpi_comm_world with ranks 0...size-1
  subroutine get_comm_size(size, comm)
    integer, intent(in) :: size
    integer, intent(out):: comm

    integer :: tmp
    integer :: ierror
    integer :: world_size
    integer :: mpicomm
    integer :: group,newgroup
    integer, dimension(3,1) :: range

    ! construct mpi group
    call mpi_comm_group(MPI_COMM_WORLD, group, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi comm group fail, error=',ierror)
    ! call mpi_group_size(group, tmp, ierror)
    ! if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi group size fail, error=',ierror)
    ! print *, "Group comm world size", tmp

    call mpi_comm_size(MPI_COMM_WORLD, world_size, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi comm size fail, error=',ierror)

    if (size < world_size) then
       range(1,1) = size
       range(2,1) = world_size-1
       range(3,1) = 1 ! stride
       !print *, range
       call mpi_group_range_excl(group, 1, range, newgroup, ierror)
       if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi group range excl fail, error=',ierror)
    else
       newgroup = group
    end if

    ! sanity check
    call mpi_group_size(newgroup, tmp, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi group size fail, error=',ierror)
    if (tmp /= size) call pf_stop(__FILE__,__LINE__,'something went wrong',ierror)

    ! construct new communicator
    call mpi_comm_create_group(MPI_COMM_WORLD, newgroup, 1337, mpicomm, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi comm create fail, error=',ierror)
    call mpi_comm_size(mpicomm, tmp, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi comm size fail, error=',ierror)
    if (tmp /= size) call pf_stop(__FILE__,__LINE__,'something went wrong',ierror)

    if (group /= MPI_GROUP_NULL) then
       call mpi_group_free(group, ierror)
       if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierror)
    end if
    if (newgroup /= MPI_GROUP_NULL) then
       call mpi_group_free(newgroup, ierror)
       if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierror)
    end if

    ! if (comm /= MPI_COMM_WORLD) then
    !    call mpi_comm_free(comm, ierror)
    !    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi comm free fail, error=',ierror)
    ! end if
    print *, "get_comm_size finished with a comm of size ", size
    comm = mpicomm
    call mpi_barrier(comm, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi barrier fail, error=',ierror)
  end subroutine get_comm_size

  ! join an existing pfasst run and create communication
  subroutine pf_dynres_create_comm(this, comm)
    type(pf_dynres_t), intent(inout) :: this
    type(pf_comm_t)  , intent(out)   :: comm

    integer :: ierror
    integer :: mpicomm
    integer :: size

    if (this%is_dynamic_start) then
       print *, "dynamic start!!!!!!!!!!!!!! create comm"
       ! receive size of new communicator
       call mpi_bcast(size, 1, MPI_INT, 0, MPI_COMM_WORLD, ierror)

       ! create communicator of that size
       call get_comm_size(size, mpicomm)
    else
        ! start with mpi_comm_world
        mpicomm = MPI_COMM_WORLD
    end if

    ! create pfasst communicator
    call pf_mpi_create(comm, mpicomm)

    ! BARRIER A
    call mpi_barrier(comm%comm, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi barrier fail, error=',ierror)

    print *, "Passed Barrier A"
  end subroutine pf_dynres_create_comm


  ! submit resource change preference to resource manager
  subroutine pf_dynres_submit_preference(pf)
    type(pf_pfasst_t), intent(inout)   :: pf
    print *, "TODO"
  end subroutine pf_dynres_submit_preference

  ! check for resource change & optionally resize
  subroutine pf_dynres_resize(pf, q0, q0len)
    type(pf_pfasst_t), intent(inout) :: pf
    integer          , intent(in)    :: q0len
    real(pfdp)       , intent(inout) :: q0(q0len)

    real    :: r
    integer :: old_size, new_size, world_size
    integer :: ierror
    integer :: mpicomm

    call mpi_comm_size(pf%comm%comm, old_size, ierror)

    if (pf%rank == 0) then
       call mpi_comm_size(MPI_COMM_WORLD, world_size, ierror)
       call random_number(r)
       new_size = 1+FLOOR(r * world_size)
    end if

    ! broadcast to world (including idle ranks)
    call mpi_bcast(new_size, 1, MPI_INT, 0, MPI_COMM_WORLD, ierror)

    ! broadcast again to world (including idle ranks for pf_dynres_create_comm)
    call mpi_bcast(new_size, 1, MPI_INT, 0, MPI_COMM_WORLD, ierror)

    print *, "resize, new size: ", new_size

    if (new_size > old_size) then
       ! new ranks will be started
       ! create new communicator
       call get_comm_size(new_size, mpicomm)

       ! update pf communicator
       call pf_mpi_destroy(pf%comm)
       call pf_mpi_create(pf%comm, mpicomm)
       call pf_mpi_setup(pf%comm, pf, ierror)

       ! BARRIER A
       call mpi_barrier(pf%comm%comm, ierror)
       print *, "Passed Barrier A"

       ! sync the state
       print *, "Sync state"
       if (pf%rank == 0) pf%dynres%sync_root = .TRUE.
       call pf_dynres_sync_state(pf, q0, q0len)
    end if
    if (new_size < old_size) then
       if (pf%rank >= new_size) then
          pf%dynres%needs_shutdown = .TRUE.
       else
          ! create new communicator
          call get_comm_size(new_size, mpicomm)

          ! update pf communicator
          call pf_mpi_destroy(pf%comm)
          call pf_mpi_create(pf%comm, mpicomm)
          call pf_mpi_setup(pf%comm, pf, ierror)
       end if
    end if
  end subroutine pf_dynres_resize



  ! Receive information on current pfasst run
  subroutine pf_dynres_join_run(pf, q0, q0len)
    type(pf_pfasst_t), intent(inout) :: pf
    integer          , intent(in)    :: q0len
    real(pfdp)       , intent(out)   :: q0(q0len)

    pf%dynres%sync_root = .FALSE. ! we are only receiving
    print *, "Sync state (join_run)"
    call pf_dynres_sync_state(pf, q0, q0len)
  end subroutine pf_dynres_join_run


  subroutine pf_dynres_sync_state(pf, q0, q0len)
    type(pf_pfasst_t), intent(inout) :: pf
    integer          , intent(in)    :: q0len
    real(pfdp)       , intent(inout) :: q0(q0len)

    integer               :: root, rank
    integer               :: ierror
    integer, dimension(2) :: buf

    ! determine root via pf%dynres%sync_root
    rank = 0
    if (pf%dynres%sync_root) rank = pf%rank
    call mpi_barrier(pf%comm%comm, ierror)
    call mpi_allreduce(rank, root, 1, MPI_INT, MPI_SUM, pf%comm%comm, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi allreduce fail, error=',ierror)

    ! share run state
    buf(1) = pf%state%steps_done
    buf(2) = pf%state%pfblock
    call mpi_bcast(buf, 2, MPI_INT, root, pf%comm%comm, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierror)
    pf%state%steps_done = buf(1)
    pf%state%pfblock = buf(2)
    print *, "Shared/received steps_done: ", pf%state%steps_done, ", pfblock: ", pf%state%pfblock

    ! share initial condition of next block
    call mpi_bcast(q0, q0len, myMPI_Datatype, root, pf%comm%comm, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierror)

  end subroutine pf_dynres_sync_state







  !> ====================================================
  !> Alternative pfasst setup routines
  !> Don't do anything magical, but we don't want to break compatibility
  !> ====================================================

  !> Same as pf_pfasst_create, but without setting up comm and results dir and with adding dynres
  subroutine pf_pfasst_create_dynamic(pf, nlevels, dynres, fname, nocmd)
    use pf_mod_hooks, only: PF_MAX_HOOK


    type(pf_pfasst_t), intent(inout)                   :: pf        !! Main pfasst object
    integer,           intent(in   ), optional         :: nlevels   !! number of pfasst levels
    type(pf_dynres_t), intent(in   ), target, optional :: dynres    !! dynres object to add to pfasst
    character(len=*),  intent(in   ), optional         :: fname     !! Input file for pfasst parameters
    logical,           intent(in   ), optional         :: nocmd     !! Determines if command line variables are to be read

    logical :: read_cmd              !! Local version of nocmd
    integer :: ierr                  !! Record system call error
    integer :: l                     !!  Loop variable for levels
    integer :: system                !!  For opening directory
    if (present(nlevels)) pf%nlevels = nlevels


    !> gather some input from a file and command line
    read_cmd = .true.
    if (present(nocmd)) then
         if (nocmd) read_cmd = .false.
    end if
    if (present(fname)) then      !!  fname  present,  read inputs from a file (and maybe command line)
       call pf_read_opts(pf, read_cmd, fname)
    else                           !!  fname not present, only call read_opts if we want command line read
       if (read_cmd) call pf_read_opts(pf, read_cmd)
    end if

    ! add dynres to pf if present
    if (present(dynres)) then
       allocate(pf%dynres,stat=ierr)
       if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate error dynres")

       pf%dynres => dynres
    end if

    !>  allocate level pointers
    allocate(pf%levels(pf%nlevels),stat=ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate error",pf%nlevels)
    !>  loop over levels to set parameters
    do l = 1, pf%nlevels
       pf%levels(l)%index = l
       pf%levels(l)%nsweeps = pf%nsweeps(l)
       pf%levels(l)%nsweeps_pred = pf%nsweeps_pred(l)
       pf%levels(l)%nnodes = pf%nnodes(l)
       pf%levels(l)%Finterp = pf%Finterp
    end do

    !>  allocate hooks
    allocate(pf%hooks(pf%nlevels, PF_MAX_HOOK, PF_MAX_HOOKS),stat=ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate error hooks")
    allocate(pf%nhooks(pf%nlevels, PF_MAX_HOOK),stat=ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate error nhooks")
    pf%nhooks = 0

    !>  allocate status
    allocate(pf%state,stat=ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate error state")
    pf%state%pstatus = 0
    pf%state%status  = 0

    ! Create the output directory if it is not there
    ierr= system('mkdir -p dat')
    if (ierr .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make directory dat")

  end subroutine pf_pfasst_create_dynamic




  !> Setup both the PFASST object and the comm object
  subroutine pf_pfasst_setup_dynamic(pf, comm)
    type(pf_pfasst_t), intent(inout), target :: pf   !!  Main pfasst structure
    type(pf_comm_t),   intent(inout), target :: comm      !! Communicator

    class(pf_level_t), pointer :: f_lev, c_lev  !!  Pointers to level structures for brevity
    integer                   :: l                      !!  Level loop index
    integer                   :: ierr                   !!  error flag
    character(len=5)          :: dirname     ! used to append output directory


    !>  set communicator
    pf%comm => comm

    !>  Set up the mpi communicator buffers and rank
    call pf_mpi_setup(pf%comm, pf,ierr)
    if (ierr /=0 )  call pf_stop(__FILE__,__LINE__,"ERROR: mpi_setup failed")

    if (pf%rank < 0) then
       call pf_stop(__FILE__,__LINE__,&
            "Invalid PF rank: did you call setup correctly?")
    end if


    !>  loop over levels to set parameters
    do l = 1, pf%nlevels
       call pf_level_setup(pf, l)
    end do
    !>  set default finest level
    pf%state%finest_level=pf%nlevels
    !>  Loop over levels setting interpolation and restriction matrices (in time)
    do l = pf%nlevels, 2, -1
       f_lev => pf%levels(l); c_lev => pf%levels(l-1)
       allocate(f_lev%tmat(f_lev%nnodes,c_lev%nnodes),stat=ierr)
       if (ierr /= 0)  call pf_stop(__FILE__,__LINE__,"allocate fail",f_lev%nnodes)

       allocate(f_lev%rmat(c_lev%nnodes,f_lev%nnodes),stat=ierr)
       if (ierr /= 0) call pf_stop(__FILE__,__LINE__,"allocate fail",f_lev%nnodes)


       ! with the RK stepper, no need to interpolate and restrict in time
       ! we only copy the first node and last node betweem levels
       if (pf%use_rk_stepper .eqv. .true.) then
          f_lev%tmat = 0.0_pfdp
          f_lev%rmat = 0.0_pfdp

          f_lev%tmat(1,1) = 1.0_pfdp
          f_lev%tmat(f_lev%nnodes,c_lev%nnodes) = 1.0_pfdp

          f_lev%rmat(1,1) = 1.0_pfdp
          f_lev%rmat(c_lev%nnodes,f_lev%nnodes) = 1.0_pfdp
       else         ! else compute the interpolation matrix
          call pf_time_interpolation_matrix(f_lev%nodes, f_lev%nnodes, c_lev%nodes, c_lev%nnodes, f_lev%tmat)
          call pf_time_interpolation_matrix(c_lev%nodes, c_lev%nnodes, f_lev%nodes, f_lev%nnodes, f_lev%rmat)
       endif
    end do


    !  Stick the number of processors on the end of the output directory
    write (dirname, "(A1,I0.4)") 'P',pf%comm%nproc
    pf%outdir       = trim(pf%outdir)//trim(dirname)
    ierr= system('mkdir -p dat/' // trim(pf%outdir))
    if (ierr .ne. 0) call pf_stop(__FILE__,__LINE__, "Cannot make base directory")

  end subroutine pf_pfasst_setup_dynamic

end module pf_mod_dynres
