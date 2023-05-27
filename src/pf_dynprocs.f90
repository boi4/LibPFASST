!!  Implementation of Dynamic MPI support
!
! This file is part of LIBPFASST.
!
module pf_mod_dynprocs
  use pf_mod_pfasst
  use pf_mod_hooks
  use mpi

  implicit none
contains


  !> =========================================================
  !> Main functions to be called by the user
  !> =========================================================

  ! Constructor for new dynprocs object using initialized MPI session
  ! when running multiple pfasst instances in parallel, global_set and horizontal_set must be provided
  subroutine pf_dynprocs_create(this, session, main_pset, global_pset, horizontal_pset)
    type(pf_dynprocs_t), intent(out)          :: this
    integer            , intent(in)           :: session
    character(len=*)   , intent(in)           :: main_pset
    character(len=*)   , intent(in), optional :: global_pset
    character(len=*)   , intent(in), optional :: horizontal_pset

    character(len=100)               :: buffer
    integer                          :: info
    integer                          :: ierr
    logical                          :: contains_key
    integer                          :: global_mpi_group


    ! allocate strings
    allocate(character(len=MPI_MAX_PSET_NAME_LEN) :: this%main_pset)
    allocate(character(len=MPI_MAX_PSET_NAME_LEN) :: this%delta_pset)
    allocate(character(len=MPI_MAX_PSET_NAME_LEN) :: this%global_pset)
    allocate(character(len=MPI_MAX_PSET_NAME_LEN) :: this%horizontal_pset)

    ! setup variables
    this%session = session
    this%needs_shutdown = .FALSE.
    this%main_pset = main_pset

    if (present(global_pset)) then
       if (.not. present(horizontal_pset)) call pf_stop(__FILE__,__LINE__,'fatal: horizontal_pset must be present if global_pset is present')
       this%global_used = .true.
       this%global_pset = global_pset
       this%horizontal_pset = horizontal_pset
    end if

    ! Check if dynamic start
    call MPI_Session_get_pset_info(session, "mpi://WORLD", info, ierr)
    call MPI_Info_get(info, "mpi_dyn", 100, buffer, contains_key, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info get fail, error=',ierr)
    call MPI_Info_free(info, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)
    this%is_dynamic_start = (contains_key .and. trim(buffer) == "True")
  end subroutine pf_dynprocs_create


  ! destroy dynprocs object
  ! MPI Session must be finalized by application!
  subroutine pf_dynprocs_destroy(this)
    type(pf_dynprocs_t), intent(out)   :: this

    integer :: ierr
    if (allocated(this%main_pset)) deallocate(this%main_pset)
    if (allocated(this%delta_pset)) deallocate(this%delta_pset)
    if (allocated(this%global_pset)) then
       deallocate(this%global_pset)
       ! call mpi_comm_free(this%global_comm, ierr)
    end if
    if (allocated(this%horizontal_pset)) then
       deallocate(this%horizontal_pset)
       ! call mpi_comm_free(this%horizontal_comm, ierr)
    end if
  end subroutine pf_dynprocs_destroy



  !> like pf_pfasst_create but with a dynproc structure instead of a communicator
  subroutine pf_pfasst_create_dynamic(pf, dynprocs, nlevels, fname, nocmd)
    use pf_mod_hooks, only: PF_MAX_HOOK

    type(pf_pfasst_t),   intent(inout)        :: pf        !! Main pfasst object
    type(pf_dynprocs_t), intent(in), target   :: dynprocs  !! dynprocs object to add to pfasst
    integer,             intent(in), optional :: nlevels   !! number of pfasst levels
    character(len=*),    intent(in), optional :: fname     !! orInput file for pfasst parameters
    logical,             intent(in), optional :: nocmd     !! Determines if command line variables are to be read

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

    ! add dynprocs to pf
    pf%dynprocs => dynprocs
    pf%use_dynprocs = .true.

    ! create pfasst communicator
    ! and start communication if joining run
    allocate(pf%comm)
    call pf_dynprocs_create_comm(pf%dynprocs, pf%comm)

    ! and setup it up
    call pf_mpi_setup(pf%comm, pf,ierr)
    if (ierr /=0 )  call pf_stop(__FILE__,__LINE__,"ERROR: mpi_setup failed")

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




  !> =========================================================
  !> Generic helper functions
  !> =========================================================

  !> helper to create an MPI communicator from the given process set
  !> can be used by applications too as it's not pfasst specific
  subroutine pf_dynprocs_comm_from_pset(session, pset, comm)
     integer,          intent(in)  :: session
     character(len=*), intent(in)  :: pset
     integer,          intent(out) :: comm

     integer :: mgroup
     integer :: ierr

    ! create communicator
    call mpi_group_from_session_pset(session, pset, mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group from pset fail, error=',ierr)
    call mpi_comm_create_from_group(mgroup, pset, MPI_INFO_NULL, MPI_ERRHANDLER_NULL, comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm create from group fail, error=',ierr)
    call mpi_group_free(mgroup, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi group free fail, error=',ierr)
  end subroutine pf_dynprocs_comm_from_pset


  !> helper to check if process set pset contains the calling process
  !> can be used by applications too as it's not pfasst specific
  subroutine pf_dynprocs_pset_contains_me(session, pset, contains_me)
     integer,          intent(in)  :: session
     character(len=*), intent(in)  :: pset
     logical,          intent(out) :: contains_me

     integer :: info
     integer :: ierr
     logical :: contains_key
     character(len=20) :: boolean_string

     call mpi_session_get_pset_info(session, pset, info, ierr)
     if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get pset info fail, error=',ierr)
     call mpi_info_get(info, "mpi_included", 20, boolean_string, contains_key, ierr)
     if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info get fail, error=',ierr)
     call mpi_info_free(info, ierr)
     if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)

     contains_me = (trim(boolean_string) == "True")
  end subroutine pf_dynprocs_pset_contains_me


  !> helper to check if process was started dynamically
  !> can be used by applications too as it's not pfasst specific
  subroutine pf_dynprocs_check_dynamic(session, is_dynamic)
    integer, intent(in) :: session
    logical, intent(out) :: is_dynamic
    integer info,ierr
    character(len=20) :: boolean_string
    logical :: contains_key

   ! Get the info from our mpi://WORLD pset
   call mpi_session_get_pset_info(session, "mpi://WORLD", info, ierr)

   ! get value for the 'mpi_dyn' key -> if true, this process was added dynamically
   call mpi_info_get(info, "mpi_dyn", 20, boolean_string, contains_key, ierr)
   if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info get fail, error=',ierr)

   is_dynamic = (contains_key .and. trim(boolean_string) == "True")
  end subroutine pf_dynprocs_check_dynamic


  !> helper to get the name of a pset op
  !> useful for debugging
  subroutine pf_dynprocs_psetop2str(psetop, str)
    integer, intent(in) :: psetop
    character(len=*), intent(out) :: str

    character(len=23) :: myStrings(11) = [ &
      "MPI_PSETOP_NULL        ", &
      "MPI_PSETOP_ADD         ", &
      "MPI_PSETOP_SUB         ", &
      "MPI_PSETOP_REPLACE     ", &
      "MPI_PSETOP_MALLEABLE   ", &
      "MPI_PSETOP_GROW        ", &
      "MPI_PSETOP_SHRINK      ", &
      "MPI_PSETOP_UNION       ", &
      "MPI_PSETOP_DIFFERENCE  ", &
      "MPI_PSETOP_INTERSECTION", &
      "MPI_PSETOP_SPLIT       " &
    ]

    str = mystrings(psetop+1)
  end subroutine pf_dynprocs_psetop2str



  !> =========================================================
  !> Functions below should usually not be called by the user
  !> =========================================================


  !> Create a PFASST communicator (pf_comm) from a dynprocs object
  !> Should only be called by pf_dynprocs_create_dynamic
  !> comm must be allocated by the user
  !> this routine also sets up properly all necessary communication
  !> related to pf_dynprocs_t
  !> Many MPI collectives are called together from here and from pf_dynprocs_handle_grow_global
  subroutine pf_dynprocs_create_comm(this, comm)
    type(pf_dynprocs_t), intent(inout) :: this
    type(pf_comm_t) , intent(out)   :: comm

    integer :: main_mpi_comm
    integer :: main_mpi_group
    integer :: info
    integer :: ierr
    integer :: size
    integer :: psets
    integer :: i
    integer :: other_instance_identifer
    integer :: local_rank
    integer :: status(MPI_STATUS_SIZE)
    integer :: tmp
    logical :: contains_key
    logical :: contains_me
    character(len=4096) :: keys(1)
    character(len=MPI_MAX_PSET_NAME_LEN) :: tmp_pset

    ! create global communicator if it exists
    if (this%global_used) then
       ! if dynamic, first get actual global_pset name from session
       ! here we expect that the union operation with the new procs was already done
       if (this%is_dynamic_start) then
          keys(1) = "pfasst://global_pset"
          call mpi_session_get_pset_data(this%session, "mpi://WORLD", "mpi://WORLD", keys, 1, 1, info, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get info fail, error=',ierr)
          ! we need to use tmp_pset, as open mpi ignores valuelen
          call mpi_info_get(info, "pfasst://global_pset", MPI_MAX_PSET_NAME_LEN, tmp_pset, contains_key, ierr)
          this%global_pset = tmp_pset


          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get info fail, error=',ierr)
          if (.NOT. contains_key) call pf_stop(__FILE__,__LINE__,'fatal: global_pset not found in mpi session info')
          call mpi_info_free(info, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)
       end if

       ! create communicator
       call pf_dynprocs_comm_from_pset(this%session, this%global_pset, this%global_comm)

       ! get global size and rank
       call mpi_comm_size(this%global_comm, this%global_size, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm size fail, error=',ierr)
       call mpi_comm_rank(this%global_comm, this%global_rank, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm rank fail, error=',ierr)

       ! do the same with horizontal pset
       call pf_dynprocs_comm_from_pset(this%session, this%horizontal_pset, this%horizontal_comm)

       ! get global size and rank
       call mpi_comm_size(this%horizontal_comm, this%horizontal_size, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm size fail, error=',ierr)
       call mpi_comm_rank(this%horizontal_comm, this%horizontal_rank, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm rank fail, error=',ierr)
    end if




    ! if we are dynamic we need to find out main pset name
    if (this%is_dynamic_start) then
       if (this%global_used) then
          ! create main communicator from main pset (currently only containing the new processes)
          call pf_dynprocs_comm_from_pset(this%session, this%main_pset, main_mpi_comm)

          call mpi_comm_rank(main_mpi_comm, local_rank, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm rank fail, error=',ierr)

          if (local_rank == 0) then
             ! we need to match with right existing pfasst instance
             ! global rank 0 will match them by checking their horizontal rank
             ! send main pset and horizontal rank to rank 0 of global communicator
             call mpi_send(this%horizontal_rank, 1, MPI_INTEGER, 0, 1338, this%global_comm, ierr)
             if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi send fail, error=',ierr)
             ! also send main pset name to global rank 0
             call mpi_send(this%main_pset, MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, 0, 1339, this%global_comm, ierr)
             if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi send fail, error=',ierr)
             ! print *, "Sending main_pset to global comm: ", trim(this%main_pset)

             !> global rank 0 creates a new main pset with the new processes and the old ones

             ! receive union pset in main_pset
             call mpi_recv(this%main_pset, MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, 0, 1342, this%global_comm, status, ierr)
             if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi recv fail, error=',ierr)
             ! print *, "Received new main_pset from global comm: ", trim(this%main_pset)
          end if
          ! broadcast main pset to all other new processes in this new libpfasst split
          call mpi_bcast(this%main_pset, MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, 0, main_mpi_comm, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierr)

          ! ! destroy old communicator as it will get replaced by the new one
          call mpi_comm_free(main_mpi_comm, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm free fail, error=',ierr)
       else
          ! if no dynamic global communicator, we need to get main pset name from session
          keys(1) = "pfasst://main_pset"
          call mpi_session_get_pset_data(this%session, "mpi://WORLD", "mpi://WORLD", keys, 1, 1, info, ierr)
          call mpi_info_get(info, "pfasst://main_pset", MPI_MAX_PSET_NAME_LEN, tmp_pset, contains_key, ierr)
          this%main_pset = tmp_pset
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get info fail, error=',ierr)
          if (.NOT. contains_key) call pf_stop(__FILE__,__LINE__,'fatal: main_pset not found in mpi session info')

          call mpi_info_free(info, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)
       end if
    end if

    !> here, a process should have this%main_pset set to a value that contains all processes that will be used
    !> when running multiple pfasst instances in parallel, this%global_pset and related stuff is set as well

    ! print *, "Create main communicator from main pset: ", trim(this%main_pset)
    ! WARNING: this call is important, as it syncs new psets here
    call pf_dynprocs_pset_contains_me(this%session, this%main_pset, contains_me)
    ! print *, "contains_me=", contains_me
    !
    ! create main communicator from main pset
    call pf_dynprocs_comm_from_pset(this%session, this%main_pset, main_mpi_comm)

    ! create rank
    call mpi_comm_rank(main_mpi_comm, local_rank, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm rank fail, error=',ierr)

    ! create pfasst communicator
    call pf_mpi_create(comm, main_mpi_comm)

    ! barrier on local pfasst instance
    call mpi_barrier(comm%comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi barrier fail, error=',ierr)
  end subroutine pf_dynprocs_create_comm



  !> helper function that does union
  !> on the last num_steps_to_shrink horizontal psets
  !> union_pset is only set on first pfasst instance (horizonal_rank == 0)
  !> if not global_used, then union_pset
  subroutine pf_dynprocs_get_shrink_union(pf, num_steps_to_shrink, union_pset)
    type(pf_pfasst_t),                    intent(inout) :: pf
    integer,                              intent(in)    :: num_steps_to_shrink
    character(len=*),                     intent(inout) :: union_pset

    character(len=MPI_MAX_PSET_NAME_LEN), allocatable  :: input_psets(:)
    character(len=MPI_MAX_PSET_NAME_LEN)               :: output_psets(2)
    character(len=30) :: splitstr
    character(len=30) :: tmpstr
    integer :: i
    integer :: ierr
    integer :: info
    integer :: noutput
    integer :: op

    if (pf%debug) print *, "pf_dynprocs_get_shrink_union"
    if (pf%debug) print *, "num_steps_to_shrink:", num_steps_to_shrink

    if (pf%dynprocs%global_used) then
       if (pf%dynprocs%horizontal_rank == 0) then
          ! get the names of the horizontal psets to remove
          allocate(input_psets(num_steps_to_shrink))
          do i=1,num_steps_to_shrink
              if (pf%rank == pf%comm%nproc - i) then
                  input_psets(i) = pf%dynprocs%horizontal_pset
              end if
              call mpi_bcast(input_psets(i), MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, pf%comm%nproc - i, pf%comm%comm, ierr)
          end do

          if (pf%rank == 0) then
              ! do a union psetop
              op = MPI_PSETOP_UNION
              call mpi_info_create(info, ierr)
              if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)

              noutput = 1
              call mpi_session_dyn_v2a_psetop(pf%dynprocs%session, op, input_psets, num_steps_to_shrink, union_pset, noutput, info, ierr)
              if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session dyn v2a psetop fail, error=',ierr)

              call mpi_info_free(info, ierr)
              if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)

              if (pf%debug) print *, "Union pset: ", trim(union_pset)
          end if

          ! broadcast union pset
          call mpi_bcast(union_pset, MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, 0, pf%comm%comm, ierr)

          deallocate(input_psets)
       end if
    else
       if (pf%rank == 0) then
          ! do a split operation on main_pset
          op = MPI_PSETOP_SPLIT
          call mpi_info_create(info, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)

          splitstr = ""
          write(tmpstr,'(I0)') (pf%comm%nproc - num_steps_to_shrink)
          splitstr = trim(splitstr)//trim(tmpstr)//","

          write(tmpstr,'(I0)') num_steps_to_shrink
          splitstr = trim(splitstr)//trim(tmpstr)

          call mpi_info_set(info, "mpi_part_sizes", splitstr, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info set fail, error=',ierr)

          noutput = 2
          print *, "calling mpi_session_dyn_v2a_psetop with part_sizes=", splitstr
          call mpi_session_dyn_v2a_psetop(pf%dynprocs%session, op, pf%dynprocs%main_pset, 1, output_psets, noutput, info, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session dyn v2a psetop fail, error=',ierr)

          call mpi_info_free(info, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)

          union_pset = output_psets(2)
       end if
    end if
  end subroutine pf_dynprocs_get_shrink_union


  !> update main pset after a shrink when running multiple pfasst instances in parallel
  !> only call this if global_used is true
  !> if we are running multiple pfasst instances in parallel,
  !> we need to remove the new delta_pset from every one of the parallel pfasst instances
  !> This should be called before psetop_finalize
  subroutine pf_dynprocs_handle_shrink_global(pf)
    type(pf_pfasst_t), intent(inout) :: pf

    character(len=MPI_MAX_PSET_NAME_LEN), allocatable  :: input_psets(:)
    integer :: delta_size
    integer :: num_steps_to_shrink
    integer :: ierr
    integer :: info
    integer :: noutput
    integer :: i
    integer :: op
    integer :: indicator
    logical :: contains_key
    logical :: contains_me
    character(len=20) buf
    character(len=MPI_MAX_PSET_NAME_LEN)  :: union_pset

    ! mpi_size key is wrongly uninitialized (as of May 2023)
    ! once it's fixed, this approach will be (probably) faster than the allreduce
    !call mpi_session_get_pset_info(pf%dynprocs%session, pf%dynprocs%delta_pset, info, ierr);
    !if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session get pset info fail, error=',ierr)

    !call mpi_info_get(info, "mpi_size", 20, buf, contains_key, ierr);
    !if (ierr /=0 .or. .not. contains_key) call pf_stop(__FILE__,__LINE__,'mpi info get fail, error=',ierr)

    !call mpi_info_free(info, ierr)
    !if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)
    !
    !read(buf,*) delta_size

    ! determine size of delta pset
    if (pf%dynprocs%needs_shutdown) then
       indicator = 1
    else
        indicator = 0
    end if
    call mpi_allreduce(indicator, delta_size, 1, MPI_INTEGER, MPI_SUM, pf%dynprocs%global_comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi allreduce fail, error=',ierr)

    ! convert to integer
    num_steps_to_shrink = delta_size / pf%dynprocs%horizontal_size
    if (MOD(delta_size, pf%dynprocs%horizontal_size) /= 0) then
       call pf_stop(__FILE__,__LINE__,'fatal: delta pset size not divisible by horizontal size')
    end if

    if (pf%debug) print *, "Will shrink by ", num_steps_to_shrink

    ! get the union pset of the last num_to_shrink timesteps
    call pf_dynprocs_get_shrink_union(pf, num_steps_to_shrink, union_pset)

    ! do diff operations across the time psets (using horizontal comm at time 0 ranks)
    if (pf%rank == 0) then
        !bcast union pset across horizontal comm
        call mpi_bcast(union_pset, MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, 0, pf%dynprocs%horizontal_comm, ierr)

        ! do diff operations across the time psets
        op = MPI_PSETOP_DIFFERENCE
        call mpi_info_create(info, ierr)
        if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)

        allocate(input_psets(2))
        input_psets(1) = pf%dynprocs%main_pset
        input_psets(2) = union_pset

        noutput = 1
        call mpi_session_dyn_v2a_psetop(pf%dynprocs%session, op, input_psets, 2, pf%dynprocs%main_pset, noutput, info, ierr)
        if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session dyn v2a psetop fail, error=',ierr)

        call mpi_info_free(info, ierr)
        if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)

        if (pf%debug) print *, "New main pset: ", trim(pf%dynprocs%main_pset)
    end if

    ! broadcast new shrunk main_pset across pf%comm%comm
    call mpi_bcast(pf%dynprocs%main_pset, MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, 0, pf%comm%comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierr)

    ! check if new main pset contains me
    call pf_dynprocs_pset_contains_me(pf%dynprocs%session, pf%dynprocs%main_pset, contains_me)
    if ((contains_me .and. pf%dynprocs%needs_shutdown) .or. (.not. contains_me .and. .not. pf%dynprocs%needs_shutdown)) then
       call pf_stop(__FILE__,__LINE__,'fatal: shrink inconsistency')
    end if
  end subroutine pf_dynprocs_handle_shrink_global



  !> update main pset after a grow when running multiple pfasst instances in parallel
  !> if we are running multiple pfasst instances in parallel,
  !> we need to split the new delta_pset into one pset per instance
  !> this is happening on the user's side
  !> then, each split pset needs to find it's correct instance pset
  !> see pf_dynprocs_create_comm for the new processes side of this
  !> on the old processes side (here), the global leader is
  !> receiving the instance id (=horizontal rank) and the split pset from the split pset leader
  !> then here it does a union and sends the union back and also to the old instance leader
  !> Some MPI collectives are called together with collectives in pf_dynprocs_create_comm
  !> This should be called after psetop_finalize
  subroutine pf_dynprocs_handle_grow_global(pf)
    type(pf_pfasst_t), intent(inout) :: pf

    integer :: main_mpi_comm
    integer :: main_mpi_group
    integer :: global_mpi_group
    integer :: info
    integer :: ierr
    integer :: noutput
    integer :: op, orig_op
    integer :: status(MPI_STATUS_SIZE)
    integer :: split_instance_identifer
    integer :: split_instance_leader
    integer :: i
    integer :: j
    integer :: tmp
    logical :: contains_me
    character(len=40)                    :: key_name
    character(len=MPI_MAX_PSET_NAME_LEN) :: base_pset
    character(len=MPI_MAX_PSET_NAME_LEN) :: input_psets(2)
    character(len=MPI_MAX_PSET_NAME_LEN) :: output_psets(1)

    if (pf%dynprocs%global_rank == 0) then
       ! receive split pset name and instance id from each new split proc
       do i = 1,pf%dynprocs%horizontal_size

          ! step 1: receive a split identifier (=horizontal rank) from a new split pset
          call mpi_recv(split_instance_identifer, 1, MPI_INTEGER, MPI_ANY_SOURCE, 1338, pf%dynprocs%global_comm, status, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi recv fail, error=',ierr)
          ! extract global rank from status
          split_instance_leader = status(MPI_SOURCE)

          ! step 2: receive additionally the split pset name from NEW process leader
          call mpi_recv(input_psets(2), MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, split_instance_leader, 1339, pf%dynprocs%global_comm, status, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi recv fail, error=',ierr)
          if (pf%debug) print *, "============================================================================"
          if (pf%debug) print *, "Received message from global rank ", split_instance_leader, " with split pset ", trim(input_psets(2)), " and instance id ", split_instance_identifer

          ! step 3: receive additionally the main pset from OLD process leader
          if (split_instance_identifer == pf%dynprocs%horizontal_rank .and. pf%rank == 0) then
             input_psets(1) = pf%dynprocs%main_pset
          else
             if (pf%debug) print *, "Waiting to receive pset name from leader ", split_instance_identifer
             ! receive existing pset name from OLD process leader
             call mpi_recv(input_psets(1), MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, split_instance_identifer, 1340, pf%dynprocs%horizontal_comm, status, ierr)
             if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi recv fail, error=',ierr)
             if (pf%debug) print *, "Received pset ", trim(input_psets(1)), " from ", status(MPI_SOURCE)
          end if

          ! step 4: do union operation
          op = MPI_PSETOP_UNION
          orig_op = op
          noutput = 1
          call mpi_info_create(info, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)

          if (pf%debug) print *, "About to do psetop"
          if (pf%debug) print *, "Unioning ", trim(input_psets(1)), " with ", trim(input_psets(2))
          call MPI_Session_dyn_v2a_psetop(pf%dynprocs%session, op, input_psets, 2, output_psets, noutput, info, ierr)
          if (pf%debug) print *, "-> Unioned to ", trim(output_psets(1))
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session dyn psetop fail, error=',ierr)
          if (op /= orig_op) call pf_stop(__FILE__,__LINE__,'mpi session dyn psetop returned wrong op')

          call mpi_info_free(info, ierr)
          if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)

          if (pf%debug) print *, "Did union op ", trim(output_psets(1))

          ! step 5: send union set back to OLD local leaders
          if (split_instance_identifer == pf%dynprocs%horizontal_rank .and. pf%rank == 0) then
             pf%dynprocs%main_pset = output_psets(1)
          else
             if (pf%debug) print *, "Sending union set to local leader ",split_instance_identifer
             call mpi_send(output_psets(1), MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, split_instance_identifer, 1341, pf%dynprocs%horizontal_comm, ierr)
             if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi send fail, error=',ierr)
             if (pf%debug) print *, "Done"
          end if

          ! step 6: send union set back to NEW split leader
          call mpi_send(output_psets(1), MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, split_instance_leader, 1342, pf%dynprocs%global_comm, ierr)
       end do
    end if ! end of global leader part

    if (pf%rank == 0 .and. pf%dynprocs%global_rank /= 0) then
       ! send local main pset name to global leader
       if (pf%debug) print *, "Sending main pset to global leader (= horizontal rank 0): ", pf%dynprocs%main_pset
       call mpi_send(pf%dynprocs%main_pset, MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, 0, 1340, pf%dynprocs%horizontal_comm, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi send fail, error=',ierr)
       if (pf%debug) print *, "Waiting to receive union set from global leader"
       ! receive union set from global leader
       call mpi_recv(pf%dynprocs%main_pset, MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, 0, 1341, pf%dynprocs%horizontal_comm, status, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi recv fail, error=',ierr)
       if (pf%debug) print *, "Done receiving union set from global leader", pf%dynprocs%main_pset
    end if

    ! broadcast union set to other local pfasst instance procs
    call mpi_bcast(pf%dynprocs%main_pset, MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, 0, pf%comm%comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierr)

    ! WARNING: this call is important, as it syncs new psets here
    call pf_dynprocs_pset_contains_me(pf%dynprocs%session, pf%dynprocs%main_pset, contains_me)
    ! print *, "contains_me=", contains_me
  end subroutine pf_dynprocs_handle_grow_global



  !> Apply a pending resource change
  !> If new processess appear, create communication with them
  subroutine pf_dynprocs_apply_rc(pf)
    type(pf_pfasst_t), intent(inout) :: pf

    character(len=MPI_MAX_PSET_NAME_LEN) :: base_pset
    character(len=30) :: key_name
    integer :: main_mpi_comm
    integer :: info
    integer :: ierr
    integer :: noutput
    logical :: am_leader

    if (pf%debug) print *, "Apply RC called"

    ! check who is the leader for resource change psetops
    if (pf%dynprocs%global_used) then
       am_leader = (pf%dynprocs%global_rank == 0)
       key_name = "pfasst://global_pset"
       base_pset = pf%dynprocs%global_pset
    else
       am_leader = (pf%rank == 0)
       key_name = "pfasst://main_pset"
       base_pset = pf%dynprocs%main_pset
    end if


    ! if we are growing, we need to give the new processes some kick-off information
    ! to establish communication
    if (am_leader .and. pf%dynprocs%rc_op == MPI_PSETOP_GROW) then
       ! Publish the name of the new main PSet on the delta Pset
       call mpi_info_create(info, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)

       call mpi_info_set(info, key_name, base_pset, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info set fail, error=',ierr)
       call mpi_session_set_pset_data(pf%dynprocs%session, pf%dynprocs%delta_pset, info, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi set pset data fail, error=',ierr)
       call mpi_info_free(info, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)
    end if

    if (pf%dynprocs%rc_op == MPI_PSETOP_SHRINK .and. pf%dynprocs%global_used) then
       ! update main psets
       call pf_dynprocs_handle_shrink_global(pf)
    end if

    ! finalize resource change
    if (am_leader) then
       if (pf%debug) print *, "About to finalize psetop"
       call MPI_Session_dyn_finalize_psetop(pf%dynprocs%session, base_pset, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session dyn finalize psetop fail, error=',ierr)
       if (pf%debug) print *, "Done finalizing psetop"
    end if

    ! warning: new processes are started now if growing

    if (pf%dynprocs%global_used .and. .not. pf%dynprocs%needs_shutdown) then
       ! create new global communicator
       ! this will hang until new processes have called the same in pf_dynprocs_create_comm
       call mpi_comm_free(pf%dynprocs%global_comm, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm free fail, error=',ierr)
       call pf_dynprocs_comm_from_pset(pf%dynprocs%session, pf%dynprocs%global_pset, pf%dynprocs%global_comm)

       ! get global size and rank (although these should stay the same)
       call mpi_comm_size(pf%dynprocs%global_comm, pf%dynprocs%global_size, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm size fail, error=',ierr)
       call mpi_comm_rank(pf%dynprocs%global_comm, pf%dynprocs%global_rank, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi comm rank fail, error=',ierr)
    end if

    if (pf%dynprocs%rc_op == MPI_PSETOP_GROW .and. pf%dynprocs%global_used) then
       ! update main pset separately
       call pf_dynprocs_handle_grow_global(pf)
    end if


    if (.not. pf%dynprocs%needs_shutdown) then
       if (pf%debug) print *, "Update main communicator from ", trim(pf%dynprocs%main_pset)
       ! create new communicator from main pset
       call pf_dynprocs_comm_from_pset(pf%dynprocs%session, pf%dynprocs%main_pset, main_mpi_comm)

       ! create new pfasst communicator
       call pf_mpi_destroy(pf%comm)
       call pf_mpi_create(pf%comm, main_mpi_comm)

       ! do some setup
       call pf_mpi_setup(pf%comm, pf, ierr)
       if (ierr /=0 )  call pf_stop(__FILE__,__LINE__,"ERROR: mpi_setup failed")

       call mpi_barrier(pf%comm%comm, ierr)
       if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi barrier fail, error=',ierr)
    end if

  end subroutine pf_dynprocs_apply_rc






  !> This function is called directly from the block loop
  !> It checks if a resource change is necessary and if yes, tries to apply it
  !> Then it checks for pending resource changes, which could either be the result of a previous
  !> LibPFASST suggestion or a suggestion from the runtime
  !> If there are pending resource changes, it applies them by updating communicators,
  !> starting new processes and flagging processes for shutdown
  subroutine pf_dynprocs_resize(pf, q0, q0len)
    type(pf_pfasst_t), intent(inout) :: pf
    integer          , intent(in)    :: q0len
    real(pfdp)       , intent(inout) :: q0(q0len)

    integer :: ierr

    if (pf%debug) print *, "Resize called"

    call call_hooks(pf, 1, PF_PRE_POT_RESIZE)

    call pf_dynprocs_suggest_rc(pf)

    ! check for resource change
    ! set rc_op, delta_pset and main/global pset
    call pf_dynprocs_check_rc(pf)
    if (pf%dynprocs%rc_op /= MPI_PSETOP_NULL) then
       call call_hooks(pf, 1, PF_PRE_RESIZE)
       ! check if process needs to shutdown and set pf%dynprocs%needs_shutdown if yes
       call pf_dynprocs_check_shutdown(pf)

       ! apply resource change
       call pf_dynprocs_apply_rc(pf)

       if (pf%dynprocs%rc_op == MPI_PSETOP_GROW) then
          call pf_dynprocs_sync_state(pf, q0, q0len)
       end if
       call call_hooks(pf, 1, PF_POST_RESIZE)
    end if

    call call_hooks(pf, 1, PF_POST_POT_RESIZE)
  end subroutine pf_dynprocs_resize



  !> if resize_delta is set, we create a psetop to grow/shrink
  !> create a GROW/SHRINK psetop to grow by pf%dynprocs%resize_delta timesteps
  !> can later query_psetop and finalize_psetop
  !> pf%dynprocs%resize_delta is reset to 0
  subroutine pf_dynprocs_suggest_rc(pf)
    type(pf_pfasst_t), intent(inout) :: pf

    character(len=30) :: str
    character(len=MPI_MAX_PSET_NAME_LEN) :: base_pset
    character(len=MPI_MAX_PSET_NAME_LEN) :: input_psets(2)
    character(len=MPI_MAX_PSET_NAME_LEN) :: output_psets(2)
    character(len=MPI_MAX_PSET_NAME_LEN) :: union_pset
    integer :: op
    integer :: ierr
    integer :: info
    integer :: ninput
    integer :: noutput
    integer :: num_proc_delta
    integer :: union_size
    logical :: am_leader

    ! broadcast resize_delta from 0
    ! only needs to be valid at horizontal_rank == 0
    call mpi_bcast(pf%dynprocs%resize_delta, 1, MPI_INTEGER, 0, pf%comm%comm, ierr)
    if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierr)


    if (pf%dynprocs%global_used) then
       am_leader = (pf%dynprocs%global_rank == 0)
       base_pset = pf%dynprocs%global_pset
       num_proc_delta = pf%dynprocs%horizontal_size * pf%dynprocs%resize_delta
    else
       am_leader = (pf%rank == 0)
       input_psets(1) = pf%dynprocs%main_pset
       base_pset = pf%dynprocs%main_pset
       num_proc_delta = pf%dynprocs%resize_delta
    end if


    if (num_proc_delta < 0 .and. (.not. pf%dynprocs%global_used .or. pf%dynprocs%horizontal_rank == 0)) then
       if (pf%comm%nproc <= 1) then
          print *, 'WARNING: cannot shrink below 1 time process'
          num_proc_delta = 0
       else
          if (am_leader) then
            ! need to create a pset of processes to shutdown
            call pf_dynprocs_get_shrink_union(pf, -pf%dynprocs%resize_delta, union_pset)
            if (pf%debug) print *, "Creating psetop shrink union pset: ", trim(union_pset)
          end if
       end if
    end if

    if (am_leader .and. num_proc_delta /= 0) then
       call mpi_info_create(info, ierr)
       if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)

       ! call psetop to create resource change
       if (num_proc_delta > 0) then
          input_psets(1) = base_pset
          ninput = 1
          op = MPI_PSETOP_GROW
          write (str, "(I4)") num_proc_delta
          call mpi_info_set(info, "mpi_num_procs_add", str, ierr)
          if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi info set fail, error=',ierr)

          if (pf%debug) print *, "Set op to GROW by ", num_proc_delta
       else
          ! only shrink by the union pset
          input_psets(1) = base_pset
          input_psets(2) = union_pset
          ninput = 2
          ! This can be replaced with a single shrink in the future
          op = MPI_PSETOP_SHRINK
          write (str, "(I4)") -num_proc_delta
          call mpi_info_set(info, "mpi_num_procs_sub", str, ierr)
          if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi info set fail, error=',ierr)

          if (pf%debug) print *, "Set op to SHRINK by ", -num_proc_delta, " using union pset ", trim(union_pset)
       end if

       noutput = 2
       if (pf%debug) print *, "Calling psetop"
       if (pf%debug) print *, "input psets(1): ", trim(input_psets(1))
       if (pf%debug .and. ninput == 2) print *, "input psets(2): ", trim(input_psets(2))
       call mpi_session_dyn_v2a_psetop(pf%dynprocs%session, op, input_psets, ninput, output_psets, noutput, info, ierr)
       if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi session dyn v2a psetop fail, error=',ierr)

       if (op == MPI_PSETOP_NULL) then
          if (pf%debug) print *, "Psetop rejected"
       end if

       call mpi_info_free(info, ierr)
       if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)
    end if

    ! reset resize_delta
    pf%dynprocs%resize_delta = 0
  end subroutine pf_dynprocs_suggest_rc


  ! Receive information on current pfasst run
  ! A hook could be added here in the future
  subroutine pf_dynprocs_join_run(pf, q0, q0len)
    type(pf_pfasst_t), intent(inout) :: pf
    integer          , intent(in)    :: q0len
    real(pfdp)       , intent(out)   :: q0(q0len)

    ! print *, "Join Run called"

    call pf_dynprocs_sync_state(pf, q0, q0len)
  end subroutine pf_dynprocs_join_run



  !> Check if there is a pending resource change psetop
  !> If yes, update delta_pset and global/main_pset but do not establish communication yet
  subroutine pf_dynprocs_check_rc(pf)
    type(pf_pfasst_t), intent(inout) :: pf

    character(len=MPI_MAX_PSET_NAME_LEN) :: pset_name
    character(len=MPI_MAX_PSET_NAME_LEN) :: output_psets(2)
    character(len=MPI_MAX_PSET_NAME_LEN) :: input_psets(2)
    integer :: ierr
    integer :: noutput
    integer :: op
    integer :: info
    integer :: comm
    logical :: am_leader

    if (pf%dynprocs%global_used) then
       pset_name = pf%dynprocs%global_pset
       comm = pf%dynprocs%global_comm
       am_leader = (pf%dynprocs%global_rank == 0)
    else
       pset_name = pf%dynprocs%main_pset
       comm = pf%comm%comm
       am_leader = (pf%rank == 0)
    end if

    noutput = 2
    if (pf%debug) print *, "Calling mpi session dyn v2a query psetop with pset_name: ", trim(pset_name)
    call mpi_session_dyn_v2a_query_psetop(pf%dynprocs%session, pset_name, pset_name, &
                                          pf%dynprocs%rc_op, output_psets, noutput, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi session dyn v2a query psetop fail, error=',ierr)
    if (pf%debug) print *, "got psetop: ", pf%dynprocs%rc_op

    if (pf%dynprocs%rc_op /= MPI_PSETOP_NULL) then
       pf%dynprocs%delta_pset = output_psets(1)

       ! in case of ADD or SUB, we do GROW/SHRINK manually
       if (pf%dynprocs%rc_op == MPI_PSETOP_ADD .or. pf%dynprocs%rc_op == MPI_PSETOP_SUB) then
          if (am_leader) then
             ! setup arguments
             call mpi_info_create(info, ierr)
             if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi info create fail, error=',ierr)

             input_psets(1) = pset_name
             input_psets(2) = pf%dynprocs%delta_pset

             if (pf%dynprocs%rc_op == MPI_PSETOP_ADD) then
                op = MPI_PSETOP_UNION
             else
                op = MPI_PSETOP_DIFFERENCE
             end if

             noutput = 1
             call mpi_session_dyn_v2a_psetop(pf%dynprocs%session, op, input_psets, 2, output_psets, noutput, info, ierr)
             if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi session dyn v2a psetop fail, error=',ierr)


             call mpi_info_free(info, ierr)
             if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi info free fail, error=',ierr)
          end if

          ! bcast union/diff result
          call mpi_bcast(output_psets(1), MPI_MAX_PSET_NAME_LEN, MPI_CHARACTER, 0, comm, ierr)
          if (ierr /= 0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierr)

          ! make output_psets consisten with a grow/shrink request
          output_psets(2) = output_psets(1)
          output_psets(1) = pf%dynprocs%delta_pset

          ! update op to grow/shrink, so apply_rc can handle accordingly
          if (pf%dynprocs%rc_op == MPI_PSETOP_ADD) then
             pf%dynprocs%rc_op = MPI_PSETOP_GROW
          else
             pf%dynprocs%rc_op = MPI_PSETOP_SHRINK
          end if
       end if

       if (pf%dynprocs%global_used) then
          pf%dynprocs%global_pset = output_psets(2)
       else
          pf%dynprocs%main_pset = output_psets(2)
       end if
    end if
  end subroutine pf_dynprocs_check_rc


  !> check if we are part of delta pset
  !> and set needs_shutdown accordingly
  subroutine pf_dynprocs_check_shutdown(pf)
    type(pf_pfasst_t), intent(inout) :: pf
    integer                          :: psets
    character(len=20)                :: boolean_string
    integer                          :: ierr
    logical                          :: contains_key
    integer                          :: info

    integer :: l

    if (pf%debug) print *, "Check Shutdown called"

    if (pf%dynprocs%rc_op == MPI_PSETOP_SHRINK) then
       call pf_dynprocs_pset_contains_me(pf%dynprocs%session, pf%dynprocs%delta_pset, pf%dynprocs%needs_shutdown)
    end if
  end subroutine pf_dynprocs_check_shutdown


  !> sync state from process w. rank 0 to other ranks
  !> syncs all relevant infos in between blocks to join a new block
  subroutine pf_dynprocs_sync_state(pf, q0, q0len)
    type(pf_pfasst_t), intent(inout) :: pf
    integer          , intent(in)    :: q0len
    real(pfdp)       , intent(inout) :: q0(q0len)

    integer               :: rank
    integer               :: ierr
    integer, dimension(2) :: buf

    if (pf%debug) print *, "Sync state called"

    call call_hooks(pf, 1, PF_PRE_SYNC)

    ! share run state
    buf(1) = pf%state%steps_done
    buf(2) = pf%state%pfblock
    call mpi_bcast(buf, 2, MPI_INT, 0, pf%comm%comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierr)
    pf%state%steps_done = buf(1)
    pf%state%pfblock = buf(2)
    if (pf%debug) print *, "Shared/received steps_done: ", pf%state%steps_done, ", pfblock: ", pf%state%pfblock

    ! share initial condition of next block
    call mpi_bcast(q0, q0len, myMPI_Datatype, 0, pf%comm%comm, ierr)
    if (ierr /=0) call pf_stop(__FILE__,__LINE__,'mpi bcast fail, error=',ierr)

    call call_hooks(pf, 1, PF_POST_SYNC)
  end subroutine pf_dynprocs_sync_state



end module pf_mod_dynprocs
