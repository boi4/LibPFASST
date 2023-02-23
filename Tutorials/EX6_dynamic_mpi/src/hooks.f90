!
! This file is part of LIBPFASST.
!
!>  User defined routines that can be called from inside libpfasst using hooks
module hooks
  use pfasst
  use pf_mod_zndarray
  implicit none
contains

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
