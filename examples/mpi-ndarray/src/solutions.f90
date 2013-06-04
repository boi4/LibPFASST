module solutions
  use pf_mod_dtype
  use probin
  use encap
  implicit none

contains

  ! Set initial condition.
  subroutine initial(q0)
    class(array1d), intent(inout) :: q0

    select case(problem)
    case (PROB_AD)
       call exact(0.0_pfdp, size(q0%array), q0%array)
    case (PROB_HEAT)
       call gaussian(q0%array)
    case (PROB_VB)
       call gaussian(q0%array)
    case (PROB_KS)
       call sinusoidal(q0%array)
    case default
       stop "ERROR: Unknown problem type."
    end select
  end subroutine initial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine gaussian(q)
    real(pfdp), intent(inout) :: q(:)
    integer :: nvars, ii, i
    double precision :: x

    nvars = size(q)

    q = 0.0d0
    do ii = -3, 3
       do i = 1, nvars
          x = Lx * (dble(i-nvars/2-1)/dble(nvars) + ii)
          q(i) = q(i) + dexp(-x**2/sigma)
       end do
    end do
  end subroutine gaussian

  subroutine sinusoidal(q)
    real(pfdp), intent(inout) :: q(:)
    integer :: nvars, ii, i
    double precision :: x

    nvars = size(q)

    q = 0.0d0
    do ii = -3, 3
       do i = 1, nvars
          x = Lx * (dble(i-nvars/2-1)/dble(nvars) + ii)
          q(i) = q(i) + 0.1*sin(3*two_pi*x/Lx) &
                      + 0.2*sin(4*two_pi*x/Lx) &
                      + 0.3*sin(7*two_pi*x/Lx)
       end do
    end do
  end subroutine sinusoidal

  subroutine exact(t, nvars, yex)
    real(pfdp), intent(in)  :: t
    integer,    intent(in)  :: nvars
    real(pfdp), intent(out) :: yex(nvars)

    integer :: i, ii
    real(pfdp) :: x

    yex = 0.0_pfdp

    
    if (nu > 0) then

       do ii = -3, 3
          do i = 1, nvars
             x = dble(i-nvars/2-1)/dble(nvars) + ii - t*v
             yex(i) = yex(i) + 1/(4.0_pfdp*pi*nu*(t+t0))**(0.5)*dexp(-x**2/(4.0_pfdp*nu*(t+t0)))
             !yex(i) = yex(i) + dcos(2.0_pfdp*pi*x)*dexp(-4.0_pfdp*pi*pi*nu*t)
          end do
       end do

    else

       do ii = -3, 3
          do i = 1, nvars
             x = dble(i-nvars/2-1)/dble(nvars) + ii - t*v
             yex(i) = yex(i) + 1.0/(4.0*pi*t0)**(0.5)*dexp(-x**2/(4.0*t0))
          end do
       end do

    end if

  end subroutine exact


end module solutions