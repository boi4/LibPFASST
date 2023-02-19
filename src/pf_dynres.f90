!!  Implementation of Dynamic MPI support
!
! This file is part of LIBPFASST.
!
module pf_mod_dynres
  use pf_mod_pfasst
  use mpidynres ! for mpi sessions and resource change api
  use mpi, only: mpi_info_get, mpi_info_free, MPI_ERRHANDLER_NULL, MPI_INFO_NULL, MPI_INT, &
                 mpi_allreduce, MPI_SUM, MPI_CHAR

  implicit none
contains


  ! create new dynres object using initialized MPI session
  subroutine pf_dynres_create(this, session)
    type(pf_dynres_t), intent(out)   :: this
    !integer,           intent(in)    :: session
    type(c_ptr),           intent(in)    :: session

    character(len=100)               :: buffer
    integer                          :: psets ! mpi info object
    integer                          :: ierror
    logical                          :: contains_key


    this%session = session
    this%needs_shutdown = .FALSE.

    ! determine if dynamic start or not
    !
    ! JAN: TODO overhead could be reduced by storing psets in dynres object
    call mpi_session_get_psets(session, MPI_INFO_NULL, psets)

    call mpi_info_get(psets, "mpi://WORLD", 100, buffer, contains_key, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi info get fail, error=',ierror)

    if (contains_key) then
       this%is_dynamic_start = .FALSE.
    else
       this%is_dynamic_start = .TRUE.
    end if

    call mpi_info_free(psets, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierror)
  end subroutine pf_dynres_create

  ! destroy dynres object
  ! Empty for now, as no dynamic allocation
  ! MPI Session must be finalized by application!
  ! TODO
  subroutine pf_dynres_destroy(this)
    type(pf_dynres_t), intent(out)   :: this
  end subroutine pf_dynres_destroy

  ! join an existing pfasst run and create communication
  subroutine pf_dynres_create_comm(this, comm)
    type(pf_dynres_t), intent(inout) :: this
    type(pf_comm_t)  , intent(out)   :: comm

    integer :: main_mpi_comm
    integer :: main_mpi_group
    integer :: info
    integer :: ierror
    integer :: size
    integer :: psets
    logical :: contains_key

    if (this%is_dynamic_start) then
       call mpi_session_get_info(this%session, info, ierror)
       if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi session get info fail, error=',ierror)

       ! TODO: update length
       call mpi_info_get(info, "pfasst://main_pset", 4096, this%main_pset, contains_key, ierror)
       if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi session get info fail, error=',ierror)
       if (.NOT. contains_key) call pf_stop(__FILE__,__LINE__,'fatal: main_pset not found in mpi session info, error=',ierror)

       call mpi_info_free(info, ierror)
       if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierror)
    else
       this%main_pset = "mpi://WORLD"
    end if

    call mpi_group_from_session_pset(this%session, this%main_pset, main_mpi_group, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi group from pset fail, error=',ierror)
    call mpi_comm_create_from_group(main_mpi_group, "", MPI_INFO_NULL, MPI_ERRHANDLER_NULL, main_mpi_comm, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi comm create from group fail, error=',ierror)
    call mpi_group_free(main_mpi_group, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierror)

    call pf_mpi_create(comm, main_mpi_comm)

    call mpi_barrier(comm%comm, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi barrier fail, error=',ierror)
  end subroutine pf_dynres_create_comm




  ! JAN: TODO
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

    ! check for resource change
    call pf_dynres_check_rc(pf)
    if (pf%dynres%rc_type /= MPIDYNRES_RC_NONE) then
       ! check if process needs to shutdown and set pf%dynres%needs_shutdown if yes
       call pf_dynres_check_shutdown(pf)

       ! determine who is responsible for syncing to new processes
       ! this needs to be done before rank information gets lost
       pf%dynres%sync_root = (pf%rank == 0)

       ! apply resource change: mpi related attributes of pf get updated
       call pf_dynres_apply_rc(pf)

       if (pf%dynres%rc_type == MPIDYNRES_RC_ADD) then
          call pf_dynres_sync_state(pf, q0, q0len)
       end if
    end if
  end subroutine pf_dynres_resize



  ! Receive information on current pfasst run
  subroutine pf_dynres_join_run(pf, q0, q0len)
    type(pf_pfasst_t), intent(inout) :: pf
    integer          , intent(in)    :: q0len
    real(pfdp)       , intent(out)   :: q0(q0len)

    pf%dynres%sync_root = .FALSE. ! we are only receiving
    call pf_dynres_sync_state(pf, q0, q0len)
  end subroutine pf_dynres_join_run




  subroutine pf_dynres_check_rc(pf)
    type(pf_pfasst_t), intent(inout) :: pf

    integer :: info ! resource change information, ignored
    integer :: ierror

    if (pf%rank == 0) then
       call mpidynres_rc_get(pf%dynres%session, pf%dynres%rc_type, &
            pf%dynres%delta_pset, pf%dynres%rc_tag, info, ierror)
       if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpidynres rc get fail, error=',ierror)
       if (info /= MPI_INFO_NULL) then
          call MPI_Info_free(info, ierror)
          if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierror)
       end if
    end if

    call mpi_bcast(pf%dynres%rc_type, 1, MPI_INT, 0, pf%comm%comm, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierror)

    if (pf%dynres%rc_type /= MPIDYNRES_RC_NONE) then
       ! TODO: adjust length
       call mpi_bcast(pf%dynres%delta_pset, 4096, MPI_CHAR, 0, pf%comm%comm, ierror)
       if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierror)
    end if
  end subroutine pf_dynres_check_rc


  subroutine pf_dynres_check_shutdown(pf)
    type(pf_pfasst_t), intent(inout) :: pf
    integer                          :: psets
    character(len=4096)              :: buffer
    integer                          :: ierror
    logical                          :: contains_key

    integer :: l

    if (pf%dynres%rc_type == MPIDYNRES_RC_SUB) then
       ! get psets that contain this process
       call mpi_session_get_psets(pf%dynres%session, MPI_INFO_NULL, psets)

       ! check if delta pset contains this process
       call mpi_info_get(psets, pf%dynres%delta_pset, 4096, buffer, contains_key, ierror)
       if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi info get fail, error=',ierror)

       if (contains_key) then
          print *, "Need to shutdown" ! TODO: remove
          pf%dynres%needs_shutdown = .TRUE.
       end if
    end if
  end subroutine pf_dynres_check_shutdown


  subroutine pf_dynres_apply_rc(pf)
    type(pf_pfasst_t), intent(inout) :: pf
    !character(len=MPI_AX_PSET_NAME_LEN) :: new_main_pset
    character(len=4096) :: new_main_pset

    integer :: main_mpi_comm
    integer :: main_mpi_group
    integer :: info
    integer :: ierror

    if (pf%dynres%rc_type == MPIDYNRES_RC_ADD) then
       if (pf%rank == 0) then
          ! do union set operation
          call mpidynres_pset_create_op(pf%dynres%session, MPI_INFO_NULL, &
               pf%dynres%main_pset, pf%dynres%delta_pset, &
               MPIDYNRES_PSET_UNION, new_main_pset)
       end if
    else if (pf%dynres%rc_type == MPIDYNRES_RC_SUB) then
       if (pf%rank == 0) then
          ! do difference set operation
          call mpidynres_pset_create_op(pf%dynres%session, MPI_INFO_NULL, &
               pf%dynres%main_pset, pf%dynres%delta_pset, &
               MPIDYNRES_PSET_DIFFERENCE, new_main_pset)
       end if
    end if

    ! broadcast new main pset
    call mpi_bcast(new_main_pset, MPI_MAX_PSET_NAME_LEN, MPI_CHAR, 0, pf%comm%comm, ierror)
    if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierror)

    pf%dynres%main_pset = new_main_pset

    ! accept resource change
    if (pf%rank == 0) then
      call mpi_info_create(info, ierror)
      if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierror)
      call mpi_info_set(info, "pfasst://main_pset", pf%dynres%main_pset, ierror)
      if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi info set fail, error=',ierror)

      call mpidynres_rc_accept(pf%dynres%session, pf%dynres%rc_tag, info, ierror)
      if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi rc accept fail, error=',ierror)

      call mpi_info_free(info, ierror)
      if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierror)
    end if

    if (.NOT. pf%dynres%needs_shutdown) then
       ! create communicator
       call mpi_group_from_session_pset(pf%dynres%session, pf%dynres%main_pset, main_mpi_group, ierror)
       if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi group from pset fail, error=',ierror)
       call mpi_comm_create_from_group(main_mpi_group, "", MPI_INFO_NULL, MPI_ERRHANDLER_NULL, main_mpi_comm, ierror)
       if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi comm create from group fail, error=',ierror)
       call mpi_group_free(main_mpi_group, ierror)
       if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierror)

       ! TODO: optimization potential
       ! create new pfasst communicator
       call pf_mpi_destroy(pf%comm)
       call pf_mpi_create(pf%comm, main_mpi_comm)

       ! do some setup
       call pf_mpi_setup(pf%comm, pf,ierror)
       if (ierror /=0 )  call pf_stop(__FILE__,__LINE__,"ERROR: mpi_setup failed")

       call mpi_barrier(pf%comm%comm, ierror)
       if (ierror /=0) call pf_stop(__FILE__,__LINE__,'mpi barrier fail, error=',ierror)

       pf%state%proc = pf%rank+1
    end if
  end subroutine pf_dynres_apply_rc


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
