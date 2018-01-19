!
! Copyright (C) 2012, 2013 Matthew Emmett and Michael Minion.
!
! This file is part of LIBPFASST.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!
!> Module with random subroutines that don't seem to fit in other modules
module pf_mod_utils
  use pf_mod_dtype
  use pf_mod_timer
  implicit none
contains



  !
  ! Compute full residual
  !
  ! During the process of computing the residual we compute the '0 to
  ! node' integral and store it in I.  This is used later when doing
  ! restriction (see restrict_time_space_fas).
  !
  subroutine pf_residual(pf, lev, dt)
    type(pf_pfasst_t), intent(inout) :: pf
    class(pf_level_t),  intent(inout) :: lev
    real(pfdp),        intent(in)    :: dt

    real(pfdp) :: norms(lev%nnodes-1)
    integer :: m

!    if (pf%nlevels == 1 .and. pf%abs_res_tol == 0 .and. pf%rel_res_tol == 0) return
!   I think we often want the residual for diagnostics.  Maybe need flag to turn this off
!   for efficiency?

    call start_timer(pf, TRESIDUAL)

    call lev%ulevel%sweeper%residual(lev, dt)

    ! compute max residual norm
    do m = 1, lev%nnodes-1
       norms(m) = lev%R(m)%norm()
!       print *, 'norm(', m, ') = ', norms(m)
    end do
!    lev%residual = maxval(abs(norms))
    lev%residual = norms(lev%nnodes-1)

    call end_timer(pf, TRESIDUAL)

  end subroutine pf_residual

  !
  ! Generic residual
  !
  subroutine pf_generic_residual(this, lev, dt)
    class(pf_sweeper_t), intent(in)  :: this
    class(pf_level_t),  intent(inout) :: lev
    real(pfdp),        intent(in)    :: dt

    integer :: m

    call lev%ulevel%sweeper%integrate(lev, lev%Q, lev%F, dt, lev%I)
!    do m = 2, lev%nnodes-1
!       call lev%I(M)%axpy(1.0_pfdp, lev%I(m-1))
!    end do

    ! add tau (which is 'node to node')
    if (allocated(lev%tauQ)) then
       do m = 1, lev%nnodes-1
          call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m))
       end do
    end if

    ! subtract out Q
    do m = 1, lev%nnodes-1
       call lev%R(m)%copy(lev%I(m))
       call lev%R(m)%axpy(1.0_pfdp, lev%Q(1))
       call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1))
    end do

  end subroutine pf_generic_residual

  !
  ! Generic evaluate all
  !
  subroutine pf_generic_evaluate_all(this, lev, t)
    class(pf_sweeper_t), intent(in)  :: this
    class(pf_level_t),  intent(inout) :: lev
    real(pfdp),        intent(in)    :: t(:)

    integer :: m
    do m = 1, lev%nnodes
       call lev%ulevel%sweeper%evaluate(lev, t(m), m)
    end do
  end subroutine pf_generic_evaluate_all

  subroutine myLUq(Q,Qtil,Nnodes,fillq)
    real(pfdp),       intent(in)    :: Q(Nnodes-1,Nnodes)
    real(pfdp),     intent(inout)   :: Qtil(Nnodes-1,Nnodes)
    integer,        intent (in)     :: Nnodes
    integer,        intent (in)     :: fillq

    ! Return the Qtil=U^T where U is the LU decomposition of Q without pivoting
    ! if fillq is positive, then the first row of Qtil is filled to make
    ! the matrix consistent

    integer :: i,j,N
    real(pfdp) :: c
    real(pfdp)  :: U(Nnodes-1,Nnodes-1)
    real(pfdp)  :: L(Nnodes-1,Nnodes-1)
    L = 0.0_pfdp
    U = 0.0_pfdp
    N = Nnodes-1
    U=transpose(Q(1:Nnodes-1,2:Nnodes))

    do i = 1,N
       if (U(i,i) /= 0.0) then
          do j=i+1,N
             c = U(j,i)/U(i,i)
             U(j,i:N)=U(j,i:N)-c*U(i,i:N)
             L(j,i)=c
          end do
       end if
       L(i,i) = 1.0_pfdp
    end do

    !  Check
    print *,'LU error',matmul(L,U)-transpose(Q(1:Nnodes-1,2:Nnodes))

    Qtil = 0.0_pfdp
    Qtil(1:Nnodes-1,2:Nnodes)=transpose(U)
    !  Now scale the columns of U to match the sum of A
    if (fillq .eq. 1) then
       do j=1,Nnodes-1
          Qtil(j,1)=sum(Q(j,1:Nnodes))-sum(U(j,1:Nnodes-1))
       end do
    end if

  end subroutine myLUq


end module pf_mod_utils
