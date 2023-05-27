!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pfasst
  use pf_mod_zndarray
  implicit none
contains

  !>  Resize libpfasst randomly
  subroutine resize_decider(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    !integer :: max_timesteps = 8
    integer :: max_timesteps = 32
    integer :: cur_timesteps
    integer :: new_timesteps
    real :: u
    character(len=20) :: node_size_s
    integer :: node_size
    integer :: status
    integer :: len


    ! we only can set resize_delta at the process that calls the psetop
    if (pf%rank == 0 .and. ((.not. pf%dynprocs%global_used) .or. pf%dynprocs%horizontal_rank == 0)) then
        ! currently, we can only grow and shrink in the same granularity as procs per node
        len = 20
        call get_environment_variable("OMPI_COMM_WORLD_LOCAL_SIZE", node_size_s, len, status)
        if (status == 0) then
            read(node_size_s,*) node_size
        else
           print *, "Could not get environment variable OMPI_COMM_WORLD_LOCAL_SIZE"
           ! set default to 1
           node_size = 1
        end if

        cur_timesteps = pf%comm%nproc
        ! get random number between 1 and max_timesteps/node_size
        ! and subtract cur_timesteps from it
        call random_number(u)
        new_timesteps = 1 + floor(u * (max_timesteps/node_size +1 - 1))
        new_timesteps = new_timesteps*node_size
        print *, "Trying to resize to ", new_timesteps, " parallel timesteps"
        pf%dynprocs%resize_delta = new_timesteps - cur_timesteps
        print *, "Set resize_delta to ", pf%dynprocs%resize_delta
    end if
  end subroutine resize_decider



  !>  Output the error and residual in the solution
  subroutine echo_error(pf, level_index)
    type(pf_pfasst_t), intent(inout) :: pf
    integer, intent(in) :: level_index

    ! real(pfdp) :: yexact(pf%levels(level_index)%lev_shape(1))
    ! real(pfdp) :: maxerr
    complex(pfdp), pointer :: y_end(:)
    real(pfdp) ::   time !,resid
    integer ::   step,rank,iter
    time=pf%state%t0+pf%state%dt
    step=pf%state%step+1
    rank=pf%rank
    iter=pf%state%iter
    !resid=pf%levels(level_index)%residual

    !>  compute the error at last end point
    y_end => get_array1d(pf%levels(level_index)%qend)
    ! call exact(time, yexact)
    ! maxerr = maxval(abs(y_end-yexact))

     !print '("time: ", f10.4," step: ", i7.7," rank: ",i3.3," iter: ",i4.3," level: ",i2.2," real: ",f20.10," complex: ",f20.20)', &
          !time, step, rank, iter, level_index, real(y_end), aimag(y_end)

    ! call pf_set_error(pf,level_index,maxerr)

  end subroutine echo_error
  


end module hooks
