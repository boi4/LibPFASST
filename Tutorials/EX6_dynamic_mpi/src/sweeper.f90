!
! This file is part of LIBPFASST.
!


!> Sweeper and RHS routines
module my_sweeper
  use pf_mod_dtype
  use pf_mod_zndarray
  use pf_mod_imex_sweeper
  use pf_mod_solutions

  use, intrinsic :: iso_c_binding

  implicit none

  !>  extend the imex sweeper type with stuff we need to compute rhs
  type, extends(pf_imex_sweeper_t) :: my_sweeper_t

   contains

     procedure :: f_eval
     procedure :: f_comp
     procedure :: initialize

  end type my_sweeper_t

contains


  !>  Routine to initialize sweeper (bypasses imex sweeper initialize)
  subroutine initialize(this, pf, level_index)
    class(my_sweeper_t), intent(inout) :: this
    type(pf_pfasst_t),   intent(inout),target :: pf
    integer,             intent(in)    :: level_index

    !  Call the imex sweeper initialize
    call this%imex_initialize(pf, level_index)

    this%explicit=.TRUE.
    this%implicit=.FALSE.
  end subroutine initialize

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! These routines must be provided for the sweeper
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Evaluate the explicit function at y, t.
  subroutine f_eval(this, y, t, level_index, f, piece)
    use probin, only:  lam
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(in   ) :: y
    class(pf_encap_t),   intent(inout) :: f
    real(pfdp),          intent(in   ) :: t
    integer,             intent(in   ) :: level_index
    integer,             intent(in   ) :: piece  !  Which piece to solve for

    complex(pfdp),       pointer :: yvec(:), fvec(:)

    !  Grab the arrays from the encap
    yvec  => get_array1d(y)
    fvec => get_array1d(f)

    select case (piece)
    case (1)  ! Explicit piece
       fvec = lam * yvec
    case DEFAULT
       print *,'Bad case for piece in f_eval ', piece
       call exit(0)
    end select
  end subroutine f_eval

  ! Dummy function
  subroutine f_comp(this, y, t, dtq, rhs, level_index, f,piece)
    use probin, only:  lam
    class(my_sweeper_t), intent(inout) :: this
    class(pf_encap_t),   intent(inout) :: y       !  The solution we seek
    real(pfdp),          intent(in   ) :: t       !  Equation time of implicit solve
    real(pfdp),          intent(in   ) :: dtq     !  The 
    class(pf_encap_t),   intent(in   ) :: rhs     !  The right hand side of the solve
    integer,             intent(in   ) :: level_index !  Which level this is
    class(pf_encap_t),   intent(inout) :: f       !  The function value
    integer,             intent(in   ) :: piece   !  Designates which piece to solve for (here implicit)

    print *,"This should never be called"
  end subroutine f_comp


end module my_sweeper
