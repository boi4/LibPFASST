!
! This file is part of LIBPFASST.
!
!


!> Sweeper and RHS routines
module my_sweeper
  use pf_mod_dtype
  use pf_mod_zndarray
  use pf_mod_imex_sweeper
  use pf_mod_fftpackage
  use pf_mod_solutions

  use, intrinsic :: iso_c_binding

  implicit none

  !>  extend the imex sweeper type with stuff we need to compute rhs
  type, extends(pf_imex_sweeper_t) :: my_sweeper_t
     ! integer ::     nx   !  Grid size

     !>  FFT and Spectral derivatives
     ! type(pf_fft_t), pointer :: fft_tool
     ! complex(pfdp), allocatable :: opE(:) ! Explicit spectral operator
     ! complex(pfdp), allocatable :: opI(:) ! Implicit spectral operator
     !

   contains

     procedure :: f_eval    !  Computes the advection and diffusion terms
     procedure :: f_comp    !  Does implicit solves
     procedure :: initialize  !  Bypasses base sweeper initialize
     ! procedure :: destroy     !  Bypasses base sweeper destroy

  end type my_sweeper_t

contains

  ! !>  Helper function to return sweeper pointer
  ! function as_my_sweeper(sweeper) result(r)
  !   class(pf_sweeper_t), intent(inout), target :: sweeper
  !   class(my_sweeper_t), pointer :: r
  !   select type(sweeper)
  !   type is (my_sweeper_t)
  !      r => sweeper
  !   class default
  !      stop
  !   end select
  ! end function as_my_sweeper


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

!  subroutine usleep (useconds)  bind ( C, name="usleep" )
!      integer(c_int32_t), value :: useconds
!  end subroutine usleep


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

!    integer                            :: nfact
!    integer                            :: n

    complex(pfdp),       pointer :: yvec(:), fvec(:)


!    do n = 1, 10000
!      nfact = nfact * n
!    end do

    !  Grab the arrays from the encap
    yvec  => get_array1d(y)
    fvec => get_array1d(f)


    select case (piece)
    case (1)  ! Explicit piece
       fvec = lam * yvec
       ! print *, lam
       ! print *, "*"
       ! print *, yvec
       ! print *, "="
       ! print *, fvec
    ! case (2)  ! Implicit piece
    !    print *, "This should never be called!"
    case DEFAULT
       print *,'Bad case for piece in f_eval ', piece
       call exit(0)
    end select


    !print *,"f_eval, time:", t, "yvec:", yvec

  end subroutine f_eval

  ! Solve for y and return f2 also
  !   y-dtq*f(y,t) = rhs
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

    complex(pfdp),         pointer :: yvec(:), rhsvec(:), fvec(:)

    print *,"This should never be called"

    !  Grab the arrays from the encaps
    yvec  => get_array1d(y)
    rhsvec => get_array1d(rhs)
    fvec => get_array1d(f)

    ! if (imex_stat .eq. 0)  then
    !    print *,'We should not be calling fcomp for fully explicit'
    !    yvec=rhsvec
    !    fvec=0.0_pfdp
    !    return
    ! endif

    ! ! Grab the fft workspace
    ! fft => this%fft_tool

    ! if (piece == 2) then
    !    ! Apply the inverse operator with the FFT convolution
    !    call fft%conv(rhsvec,1.0_pfdp/(1.0_pfdp - dtq*this%opI),yvec)

    !    !  The function is easy to derive
    !    fvec = (yvec - rhsvec) / dtq
    ! else
    !    print *,'Bad piece in f_comp ',piece
    !    call exit(0)
    ! end if
  end subroutine f_comp



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>  Here are some extra routines which are problem dependent  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Routine to set initial condition.
  ! subroutine initial_sol(y_0)
  !   type(pf_ndarray_t), intent(inout) :: y_0
  !   call exact(0.0_pfdp, y_0%flatarray)
  ! end subroutine initial_sol

  !> Routine to return the exact solution
  ! subroutine exact(t, yex)
  !   use probin, only: nu, v,kfreq,Lx,ic_type
  !   real(pfdp), intent(in)  :: t
  !   real(pfdp), intent(out) :: yex(:)


  !   !  Call exact solution from Libpfasst for ad problem
  !   if (ic_type .eq. 1) then
  !      call exact_ad_cos(t,yex,nu,v,kfreq,Lx)  ! Cosine wave
  !   else
  !      call exact_ad_exp(t,yex,nu,v,Lx)   !  Exponential
  !   endif
  ! end subroutine exact


end module my_sweeper
