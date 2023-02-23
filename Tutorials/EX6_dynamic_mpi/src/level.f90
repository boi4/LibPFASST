!
! This file is part of LIBPFASST.
!
module my_level
  use pf_mod_dtype
  use pf_mod_fftpackage
  implicit none

  !>  extend the generic level type by defining transfer operators
  type, extends(pf_user_level_t) :: my_level_t
   contains
     procedure :: restrict
     procedure :: interpolate
  end type my_level_t

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  These are the transfer functions that must be  provided for the level
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine interpolate(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    use my_sweeper, only: my_sweeper_t

    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t                  !  Equation time
    integer, intent(in), optional :: flags                 !  Optional flags (not used here)

    print *, "This should never be called!"
  end subroutine interpolate

  !>  Restrict from fine level to coarse
  subroutine restrict(this, f_lev, c_lev, f_vec, c_vec, t, flags)
    class(my_level_t), intent(inout) :: this
    class(pf_level_t), intent(inout)      :: f_lev, c_lev  !  fine and coarse levels
    class(pf_encap_t),   intent(inout)    :: f_vec, c_vec  !  fine and coarse vectors
    real(pfdp),        intent(in   ) :: t      !<  time of solution
    integer, intent(in), optional :: flags

    print *, "This should never be called!"
  end subroutine restrict

end module my_level

