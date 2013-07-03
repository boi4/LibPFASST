!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! RHS routines for advection/diffusion example.

module spec_ad
  use pf_mod_dtype
  use pf_mod_ndarray
  implicit none
  include 'fftw3.f03'

  real(pfdp), parameter :: &
       Lx     = 1.0_pfdp, &        ! domain size
       v      = 1.0_pfdp, &        ! velocity
       nu     = 0.02_pfdp, &       ! viscosity
       t00    = 0.15_pfdp           ! initial time for exact solution

  real(pfdp), parameter :: pi = 3.141592653589793_pfdp
  real(pfdp), parameter :: two_pi = 6.2831853071795862_pfdp

  type :: ad_work_t
     type(c_ptr) :: ffft, ifft
     complex(pfdp), pointer :: wk(:)              ! work space
     complex(pfdp), pointer :: ddx(:), lap(:)     ! operators
  end type ad_work_t

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine feval_create_workspace(ctx, nvars)
    type(c_ptr), intent(out) :: ctx
    integer,     intent(in)  :: nvars

    type(ad_work_t), pointer :: work

    integer     :: i
    type(c_ptr) :: wk
    real(pfdp)  :: kx

    allocate(work)

    ! create in-place, complex fft plans
    wk = fftw_alloc_complex(int(nvars, c_size_t))
    call c_f_pointer(wk, work%wk, [nvars])

    work%ffft = fftw_plan_dft_1d(nvars, &
         work%wk, work%wk, FFTW_FORWARD, FFTW_ESTIMATE)
    work%ifft = fftw_plan_dft_1d(nvars, &
         work%wk, work%wk, FFTW_BACKWARD, FFTW_ESTIMATE)

    ! create operators
    allocate(work%ddx(nvars))
    allocate(work%lap(nvars))
    do i = 1, nvars
       if (i <= nvars/2+1) then
          kx = two_pi / Lx * dble(i-1)
       else
          kx = two_pi / Lx * dble(-nvars + i - 1)
       end if

       work%ddx(i) = (0.0_pfdp, 1.0_pfdp) * kx

       if (kx**2 < 1e-13) then
          work%lap(i) = 0.0_pfdp
       else
          work%lap(i) = -kx**2
       end if
    end do

    ctx = c_loc(work)
  end subroutine feval_create_workspace

  subroutine feval_destroy_workspace(ctx)
    type(c_ptr), intent(in) :: ctx
    type(ad_work_t), pointer :: work

    call c_f_pointer(ctx, work)

    deallocate(work%wk)
    deallocate(work%ddx)
    deallocate(work%lap)
    call fftw_destroy_plan(work%ffft)
    call fftw_destroy_plan(work%ifft)
    deallocate(work)
  end subroutine feval_destroy_workspace

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set initial condition.
  subroutine initial(q0)
    type(ndarray), intent(inout) :: q0
    call exact(0.0_pfdp, q0%flatarray)
  end subroutine initial

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine exact(t, yex)
    real(pfdp), intent(in)  :: t
    real(pfdp), intent(out) :: yex(:)

    integer    :: nvars, i, ii, nbox
    real(pfdp) :: tol, x

    nvars = size(yex)
    yex   = 0

    if (nu .gt. 0) then
       do i = 1, nvars
          x = Lx*dble(i-nvars/2-1)/dble(nvars) - t*v
          yex(i) = yex(i) + dcos(2.0_pfdp*pi*x)*dexp(-4.0_pfdp*pi*pi*nu*t)
       end do
    else

       ! decide how many images so that contribution is neglible
       tol  = 1e-16
       nbox = 1 + ceiling( sqrt( -(4.0*t00)*log((4.0*pi*(t00))**(0.5)*tol) ))

       do ii = -nbox, nbox
          do i = 1, nvars
             x = Lx*dble(i-nvars/2-1)/dble(nvars) + ii*Lx - t*v
             yex(i) = yex(i) + 1.0/(4.0*pi*t00)**(0.5)*dexp(-x**2/(4.0*t00))
          end do
       end do
    end if


  end subroutine exact

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the explicit function at y, t.
  subroutine eval_f1(yptr, t, level, ctx, f1ptr)
    type(c_ptr), intent(in), value :: yptr, f1ptr, ctx
    real(pfdp),  intent(in)        :: t
    integer,     intent(in)        :: level

    type(ad_work_t), pointer :: work
    real(pfdp),      pointer :: y(:), f1(:)
    complex(pfdp),   pointer :: wk(:)

    call c_f_pointer(ctx, work)

    y  => array1(yptr)
    f1 => array1(f1ptr)
    wk => work%wk

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)
    wk = -v * work%ddx * wk / size(wk)
    call fftw_execute_dft(work%ifft, wk, wk)

    f1 = real(wk)

  end subroutine eval_f1

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Evaluate the implicit function at y, t.
  subroutine eval_f2(yptr, t, level, ctx, f2ptr)
    type(c_ptr), intent(in), value :: yptr, f2ptr, ctx
    real(pfdp),  intent(in)        :: t
    integer,     intent(in)        :: level

    type(ad_work_t), pointer :: work
    real(pfdp),      pointer :: y(:), f2(:)
    complex(pfdp),   pointer :: wk(:)

    call c_f_pointer(ctx, work)

    y  => array1(yptr)
    f2 => array1(f2ptr)
    wk => work%wk

    wk = y
    call fftw_execute_dft(work%ffft, wk, wk)
    wk = nu * work%lap * wk / size(wk)
    call fftw_execute_dft(work%ifft, wk, wk)

    f2 = real(wk)
  end subroutine eval_f2

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solve for y and return f2 also.
  subroutine comp_f2(yptr, t, dt, rhsptr, level, ctx, f2ptr)
    type(c_ptr), intent(in), value :: yptr, rhsptr, f2ptr, ctx
    real(pfdp),  intent(in)        :: t, dt
    integer,     intent(in)        :: level

    type(ad_work_t), pointer :: work
    real(pfdp),      pointer :: y(:), rhs(:), f2(:)
    complex(pfdp),   pointer :: wk(:)

    call c_f_pointer(ctx, work)

    y  => array1(yptr)
    rhs => array1(rhsptr)
    f2 => array1(f2ptr)
    wk => work%wk

    wk = rhs
    call fftw_execute_dft(work%ffft, wk, wk)
    wk = wk / (1.0_pfdp - nu*dt*work%lap) / size(wk)
    call fftw_execute_dft(work%ifft, wk, wk)

    y  = real(wk)
    f2 = (y - rhs) / dt
  end subroutine comp_f2

  subroutine interpolate(qFp, qGp, levelF, ctxF, levelG, ctxG)
    type(c_ptr), intent(in), value :: qFp, qGp, ctxF, ctxG
    integer,     intent(in)        :: levelF, levelG

    type(ad_work_t), pointer :: workF, workG
    real(pfdp),      pointer :: qF(:), qG(:)
    complex(kind=8), pointer :: wkF(:), wkG(:)

    integer :: nvarF, nvarG, xrat

    call c_f_pointer(ctxF, workF)
    call c_f_pointer(ctxG, workG)

    qF => array1(qFp)
    qG => array1(qGp)

    nvarF = size(qF)
    nvarG = size(qG)
    xrat  = nvarF / nvarG

    if (xrat == 1) then
       qF = qG
       return
    endif

    wkF => workF%wk
    wkG => workG%wk

    wkG = qG
    call fftw_execute_dft(workG%ffft, wkG, wkG)
    wkG = wkG / nvarG

    wkF = 0.0d0
    wkF(1:nvarG/2) = wkG(1:nvarG/2)
    wkF(nvarF-nvarG/2+2:nvarF) = wkG(nvarG/2+2:nvarG)

    call fftw_execute_dft(workF%ifft, wkF, wkF)

    qF = real(wkF)
  end subroutine interpolate

  subroutine restrict(qFp, qGp, levelF, ctxF, levelG, ctxG)
    type(c_ptr), intent(in), value :: qFp, qGp, ctxF, ctxG
    integer,     intent(in)        :: levelF, levelG

    real(pfdp), pointer :: qF(:), qG(:)

    integer :: nvarF, nvarG, xrat

    qF => array1(qFp)
    qG => array1(qGp)

    nvarF = size(qF)
    nvarG = size(qG)
    xrat  = nvarF / nvarG

    qG = qF(::xrat)
  end subroutine restrict

end module spec_ad