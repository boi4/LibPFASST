!
! Copyright (C) 2016 Matthew Emmett, Michael Minion, and Sebastian Goetschel
!
! This file is part of LIBPFASST_OC.
!
! LIBPFASST is free software: you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! LIBPFASST is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with LIBPFASST.  If not, see <http://www.gnu.org/licenses/>.
!

module pf_mod_imexQ_oc
  use pf_mod_dtype
  use pf_mod_utils
  
  implicit none

  type, extends(pf_sweeper_t), abstract :: pf_imexQ_oc_t
     real(pfdp), allocatable :: QtilE(:,:)   !>  Approximate explcit quadrature rule
     real(pfdp), allocatable :: QtilI(:,:)   !>  Approximate implcit quadrature rule
     real(pfdp), allocatable :: dtsdc(:)     !>  SDC step sizes
     real(pfdp), allocatable :: QdiffE(:,:)  !>  qmat-QtilE
     real(pfdp), allocatable :: QdiffI(:,:)  !>  qmat-QtilI

     logical                 :: explicit = .true. !>  Is there an explicit piece
     logical                 :: implicit = .true. !>  Is there an implicit piece
     logical                 :: use_LUq = .true.  !>  Use the LU trick

     class(pf_encap_t), allocatable :: rhs   !> holds rhs for implicit solve
     
  contains
     procedure(pf_f_eval_p), deferred :: f_eval   !>  RHS function evaluations
     procedure(pf_f_comp_p), deferred :: f_comp   !>  Implicit solver
     !>  Set the generic functions
     procedure :: sweep        => imexQ_oc_sweep
     procedure :: initialize   => imexQ_oc_initialize
     procedure :: evaluate     => imexQ_oc_evaluate
     procedure :: integrate    => imexQ_oc_integrate
     procedure :: residual     => imexQ_oc_residual
     procedure :: evaluate_all => imexQ_oc_evaluate_all
     procedure :: spreadq0     => imexQ_oc_spreadq0
     procedure :: destroy      => imexQ_oc_destroy
     procedure :: imexQ_oc_destroy
  end type pf_imexQ_oc_t
    
  interface
     !>  This is the interface for the routine to compute the RHS function values
     subroutine pf_f_eval_p(this,y, t, level_index, f, piece, flags, idx, step)
       !>  Evaluae f_piece(y), where piece is one or two 
       import pf_imexQ_oc_t, pf_encap_t, pfdp
       class(pf_imexQ_oc_t),  intent(inout) :: this
       class(pf_encap_t),  intent(in   ) :: y        !>  Argument for evaluation
       real(pfdp),         intent(in   ) :: t        !>  Time at evaluation
       integer,            intent(in   ) :: level_index     !>  Level index
       class(pf_encap_t),  intent(inout) :: f        !>  RHS function value
       integer,            intent(in   ) :: piece           !>  Which piece to evaluate
       
       integer, intent(in)              :: flags
       integer, intent(in), optional    :: idx       ! index of quadrature node
       integer, intent(in), optional    :: step   ! time step for sequential version
       
     end subroutine pf_f_eval_p
     subroutine pf_f_comp_p(this, y, t, dtq, rhs, level_index, f, piece, flags)
       !>  Solve the equation y - dtq*f_2(y) =rhs
       import pf_imexQ_oc_t, pf_encap_t, pfdp
       class(pf_imexQ_oc_t),  intent(inout) :: this
       class(pf_encap_t), intent(inout) :: y      !>  Solution of implicit solve 
       real(pfdp),        intent(in   ) :: t      !>  Time of solve
       real(pfdp),        intent(in   ) :: dtq    !>  dt*quadrature weight
       class(pf_encap_t), intent(in   ) :: rhs    !>  RHS for solve
       integer,    intent(in   ) :: level_index   !>  Level index
       class(pf_encap_t), intent(inout) :: f      !>  f_2 of solution y
       integer,    intent(in   ) :: piece         !>  Which piece to evaluate
       
       integer,            intent(in)    :: flags
     end subroutine pf_f_comp_p
  end interface
contains

  ! Perform on SDC sweep on level Lev and set qend appropriately.
  subroutine imexQ_oc_sweep(this, pf, level_index, t0, dt, nsweeps, flags)
    use pf_mod_timer
    use pf_mod_hooks

    class(pf_imexQ_oc_t), intent(inout) :: this
    type(pf_pfasst_t), intent(inout),target :: pf  !<  PFASST structure
    real(pfdp),        intent(in   ) :: t0         !>  Time at beginning of time step
    real(pfdp),        intent(in   ) :: dt     !>  time step size
    integer,             intent(in)    :: level_index  !>  which level this is
    integer,             intent(in)    :: nsweeps      !>  number of sweeps to do
    integer, optional, intent(in   ) :: flags    

    class(pf_level_t), pointer :: lev    !>  points to current level
    
    ! indicate if sweep on both (0, default; might skip y or p if tolerance satisfied), just y (1), just p (2)
    
    integer     :: k, m, n,Nnodes, which
    real(pfdp)  :: t,tend
    real(pfdp)  :: dtq
    type(c_ptr) :: rhs

    logical     :: sweep_y, sweep_p
!     real(pfdp)  :: norms_y(Lev%nnodes-1), norms_p(Lev%nnodes-1)
    integer     ::step
    
    lev => pf%levels(level_index)   !<  Assign level pointer
    
    step = pf%state%step+1
!     print *, 'sweep on step', step


    which = 0
    if (present(flags)) which = flags
    
    Nnodes = lev%nnodes
    tend = t0+dt


    call start_timer(pf, TLEVEL+lev%index-1)

    if (which .eq. 1) then
        sweep_y = .true.
        sweep_p = .false.
    else if (which .eq. 2) then
        sweep_y = .false.
        sweep_p = .true.
    else
       sweep_y = .true.
       sweep_p = .true.
    end if

!     if( sweep_p .and. pf%rank == 0)  print *, "sweep on p with which = ", which

    
  do k = 1,nsweeps   !>  Loop over sweeps
    call call_hooks(pf, level_index, PF_PRE_SWEEP)    

    ! compute integrals from previous iteration and add fas correction
    do m = 1, Nnodes-1
       call lev%I(m)%setval(0.0_pfdp, 1)
       call lev%I(m)%setval(0.0_pfdp, 2)
!        call Lev%encap%setval(Lev%S(m), 0.0_pfdp, 1)
!        call Lev%encap%setval(Lev%S(m), 0.0_pfdp, 2)

       if( sweep_y ) then
          !  Forward in y
          if (this%explicit) then
             do n = 1, Nnodes
                call lev%I(m)%axpy(dt*this%QdiffE(m,n), lev%F(n,1),1)
!                 call Lev%encap%axpy(Lev%S(m), dt*imexQ_oc%QdiffE(m,n), Lev%F(n,1),1)
             end do
          end if
          if (this%implicit) then
             do n = 1, lev%nnodes
                call lev%I(m)%axpy(dt*this%QdiffI(m,n), lev%F(n,2), 1)
!                 call Lev%encap%axpy(Lev%S(m), dt*imexQ_oc%QdiffI(m,n), Lev%F(n,2),1)
             end do
          end if
          if (allocated(lev%tauQ)) then
            call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m), 1)
!              call Lev%encap%axpy(Lev%S(m), 1.0_pfdp, Lev%tauQ(m),1)
          end if
        end if

       if( sweep_p ) then
          !  Backward in p, note S(m) goes backward now
!2          do n =  1,Nnodes
!2             call Lev%encap%axpy(Lev%S(Nnodes-m), dt*imexQ_oc%QdiffE(m,n), Lev%F(Nnodes+1-n,1),2)
!2             call Lev%encap%axpy(Lev%S(Nnodes-m), dt*imexQ_oc%QdiffI(m,n), Lev%F(Nnodes+1-n,2),2)
!2          end do
          if (this%explicit) then
            do n =  Nnodes,1,-1
              call lev%I(m)%axpy(dt*this%QdiffE(Nnodes-m,Nnodes+1-n), lev%F(n,1), 2)
              !call Lev%encap%axpy(Lev%S(m), dt*imexQ_oc%QdiffE(Nnodes-m,Nnodes+1-n), Lev%F(n,1),2)
            end do
          end if
          if (this%implicit) then
            do n =  Nnodes,1,-1
              call lev%I(m)%axpy(dt*this%QdiffI(Nnodes-m,Nnodes+1-n), lev%F(n,2), 2)
!               call Lev%encap%axpy(Lev%S(m), dt*imexQ_oc%QdiffI(Nnodes-m,Nnodes+1-n), Lev%F(n,2),2)
            end do
          end if
          if (allocated(lev%tauQ)) then
             call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m), 2)
          end if
       end if
    end do

    ! Reload the newest initial values
    ! Recompute first function values
    if (sweep_y) then
      if (k .eq. 1) then
        call lev%Q(1)%copy(lev%q0, 1)
        if (this%explicit) &
          call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,1), 1, 1, 1, step)
        if (this%implicit) &
          call this%f_eval(lev%Q(1), t0, lev%index, lev%F(1,2), 2, 1, 1, step)
      end if
!        call Lev%encap%copy(Lev%Q(1), Lev%q0, 1)
!        call imexQ_oc%f1eval(Lev%Q(1), t0, Lev%level, Lev%ctx, Lev%F(1,1), 1, 1, step)
!        call imexQ_oc%f2eval(Lev%Q(1), t0, Lev%level, Lev%ctx, Lev%F(1,2), 1)
    end if
    !else 
    if( sweep_p ) then
      if (k .eq. 1) then
        call lev%Q(Nnodes)%copy(lev%qend, 2)
        if (this%explicit) &
          call this%f_eval(lev%Q(Nnodes), tend, lev%index, lev%F(Nnodes,1), 1, 2, Nnodes, step)
        if (this%implicit) &
          call this%f_eval(lev%Q(Nnodes), tend, lev%index, lev%F(Nnodes,2), 2, 2, Nnodes, step)
      end if
!        call Lev%encap%copy(Lev%Q(Nnodes), Lev%qend, 2)
!        call imexQ_oc%f1eval(Lev%Q(Nnodes), tend, Lev%level, Lev%ctx, Lev%F(Nnodes,1), 2, Nnodes, step)
!        call imexQ_oc%f2eval(Lev%Q(Nnodes), tend, Lev%level, Lev%ctx, Lev%F(Nnodes,2), 2)
    end if

    !  Make some space
!     call Lev%encap%create(rhs, Lev%level, SDC_KIND_SOL_FEVAL, Lev%nvars, Lev%shape, Lev%ctx)

    if (sweep_y) then
      !  Forward sweep on y
      t = t0
      do m = 1, Nnodes-1
         t = t + dt*this%dtsdc(m)  !  forward running time

         !  Form rhs with all explicit terms
         call this%rhs%setval(0.0_pfdp, 1)
         do n = 1, m
            if (this%explicit) &
              call this%rhs%axpy(dt*this%QtilE(m,n), lev%F(n,1), 1)  
            if (this%implicit) &
              call this%rhs%axpy(dt*this%QtilI(m,n), lev%F(n,2), 1)  
         end do
         
         call this%rhs%axpy(1.0_pfdp, lev%I(m), 1)
         call this%rhs%axpy(1.0_pfdp, lev%Q(1), 1)
       
        ! Do implicit solve  
         if (this%implicit) then
            call this%f_comp(lev%Q(m+1), t, dt*this%QtilI(m,m+1), this%rhs, lev%index, lev%F(m+1,2), 2, 1)
         else
            call lev%Q(m+1)%copy(this%rhs)            
         end if
          !  Compute explicit piece on new value
         if (this%explicit) &
            call this%f_eval(lev%Q(m+1), t, lev%index, lev%F(m+1,1), 1, 1, m+1, step)
      end do
      !  Reset last values
      call lev%qend%copy(lev%Q(Nnodes), 1)
      call pf_residual(pf, lev, dt, 1)
    end if

   if (sweep_p) then
   !  Backward  sweep on p
      t = tend
      do m =  Nnodes-1,1,-1
         t = t - dt*this%dtsdc(m)

         ! Do the dirk parts
         call this%rhs%setval(0.0_pfdp, 2)
         do n = Nnodes, m+1,-1
            if (this%explicit) &
              call this%rhs%axpy(dt*this%QtilE(Nnodes-m,Nnodes-n+1), lev%F(n,1), 2)  
            if (this%implicit) &
              call this%rhs%axpy(dt*this%QtilI(Nnodes-m,Nnodes-n+1), lev%F(n,2), 2)
         end do
       
         call this%rhs%axpy(1.0_pfdp, lev%I(m), 2)
         call this%rhs%axpy(1.0_pfdp, lev%Q(Nnodes), 2)

         !  Do implicit solve  
         if (this%implicit) then
           call this%f_comp(lev%Q(m), t, dt*this%QtilI(Nnodes-m,Nnodes-m+1), this%rhs, lev%index, lev%F(m,2), 2, 2)
         else
            call lev%Q(m)%copy(this%rhs)
         end if
         if (this%explicit) &
           call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1), 1, 2, m, step) 
      end do
      ! reset first value
      call lev%q0%copy(lev%Q(1), 2)
      call pf_residual(pf, lev, dt, 2)
    end if

    ! done
    call call_hooks(pf, level_index, PF_POST_SWEEP)

   end do !nsweeps
   
   call end_timer(pf, TLEVEL+lev%index-1)
  end subroutine imexQ_oc_sweep

  
  ! Evaluate function values
  subroutine imexQ_oc_evaluate(this, lev, t, m, flags, step)
    class(pf_imexQ_oc_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: t
    integer,           intent(in   ) :: m
    integer, intent(in), optional   :: flags, step
    integer                         :: which, mystep

    which = 0
    if (present(flags)) which = flags
    
    mystep = 1
    if(present(step)) then 
      mystep = step
    else
      print *, "step not present in evaluate", which
      stop
    end if
!     print *, 'evaluate step = ', mystep
    
    if (this%explicit) &
      call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,1), 1, which, m, mystep)
    if (this%implicit) &
      call this%f_eval(lev%Q(m), t, lev%index, lev%F(m,2), 2, which, m, mystep)
  end subroutine imexQ_oc_evaluate

  
  subroutine imexQ_oc_evaluate_all(this, lev, t, flags, step)
    !> Evaluate all function values
    class(pf_imexQ_oc_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    real(pfdp),        intent(in   ) :: t(:)
    integer, intent(in), optional   :: flags, step
!     call pf_generic_evaluate_all(this, lev, t, flags, step)
    integer :: m
    do m = 1, lev%nnodes
      call this%evaluate(lev, t(m), m, flags, step)
    end do
  end subroutine imexQ_oc_evaluate_all
  
  
  ! Initialize matrices
  subroutine imexQ_oc_initialize(this, lev)
    class(pf_imexQ_oc_t), intent(inout) :: this
    class(pf_level_t),    intent(inout) :: lev  
    
    integer    :: m, n, Nnodes
    
    this%npieces = 2

    Nnodes = lev%nnodes
    allocate(this%QdiffE(Nnodes-1,Nnodes))  !  S-FE
    allocate(this%QdiffI(Nnodes-1,Nnodes))  !  S-BE
    allocate(this%QtilE(Nnodes-1,Nnodes))  !  S-FE
    allocate(this%QtilI(Nnodes-1,Nnodes))  !  S-BE
    
    this%QtilE = 0.0_pfdp
    this%QtilI = 0.0_pfdp

    this%dtsdc = lev%nodes(2:Nnodes) - lev%nodes(1:Nnodes-1)

    ! Implicit matrix
    if (this%use_LUq) then 
       ! Get the LU
       call myLUq(lev%qmat,lev%LUmat,lev%nnodes,0)
       this%QtilI = lev%LUmat
!       print *,'LU',lev%LUmat
!       print *,'BE',lev%qmatBE
    else 
       this%QtilI =  lev%qmatBE
    end if

    ! Explicit matrix
    this%QtilE =  lev%qmatFE

    this%QdiffE = lev%qmat-this%QtilE
    this%QdiffI = lev%qmat-this%QtilI

    !>  Make space for rhs
    call lev%ulevel%factory%create_single(this%rhs, lev%index, lev%shape)
    
  end subroutine imexQ_oc_initialize

  
  
  ! Compute SDC integral
  subroutine imexQ_oc_integrate(this, lev, qSDC, fSDC, dt, fintSDC, flags)
    class(pf_imexQ_oc_t), intent(inout) :: this
    class(pf_level_t), intent(in   ) :: lev
    class(pf_encap_t), intent(in   ) :: qSDC(:), fSDC(:, :) !qSDC unused?
    real(pfdp),       intent(in)     :: dt
    class(pf_encap_t), intent(inout) :: fintSDC(:)
    integer, intent(in), optional    :: flags

    integer :: n, m, p,Nnodes, which
    Nnodes=lev%nnodes
    
    which = 0
    if(present(flags)) then 
      which = flags
    else
      print *, "flags not present in integrate", which
      stop
    end if
    
    do n = 1, Nnodes-1
       !  Forward in y
       if( (which .eq. 0) .or. (which .eq. 1) ) then
          call fintSDC(n)%setval(0.0_pfdp, 1)
          do m = 1, Nnodes
!               do p = 1, npieces
!                 call Lev%encap%axpy(fintSDC(n), dt*Lev%qmat(n,m), fSDC(m,p),1)
!               end do
            if (this%explicit) &
               call fintSDC(n)%axpy(dt*lev%qmat(n,m), fSDC(m,1), 1)
            if (this%implicit) &
               call fintSDC(n)%axpy(dt*lev%qmat(n,m), fSDC(m,2), 1)
          end do
       end if
       
       !  Backward in p
       if( (which .eq. 0) .or. (which .eq. 2) ) then
          call fintSDC(Nnodes-n)%setval(0.0_pfdp, 2)
          do m = 1, Nnodes
!               do p = 1, npieces
!                 call Lev%encap%axpy(fintSDC(Nnodes-n), dt*Lev%qmat(n,m), fSDC(Nnodes+1-m,p),2)
!               end do
            if (this%explicit) &
               call fintSDC(Nnodes-n)%axpy(dt*lev%qmat(n,m), fSDC(Nnodes+1-m,1), 2)
            if (this%implicit) &
               call fintSDC(Nnodes-n)%axpy(dt*lev%qmat(n,m), fSDC(Nnodes+1-m,2), 2)
          end do
       end if
    end do
  end subroutine imexQ_oc_integrate
  
  
  
  subroutine imexQ_oc_residual(this, lev, dt, flags)
    class(pf_imexQ_oc_t),  intent(inout) :: this
    class(pf_level_t), intent(inout)  :: lev
    real(pfdp),        intent(in)     :: dt
    integer, intent(in), optional     :: flags
    
    integer :: m, n, which

    which = 0
    if(present(flags)) then
      which = flags
    else
      print *, "flags not present in residual", which
      stop
    end if
    
    call this%integrate(lev, lev%Q, lev%F, dt, lev%I, which)

    ! add tau (which is 'node to node')
    if (allocated(lev%tauQ)) then
       do m = 1, lev%nnodes-1
          call lev%I(m)%axpy(1.0_pfdp, lev%tauQ(m), which)
       end do
    end if

    ! subtract out Q
    do m = 1, lev%nnodes-1
       if( (which .eq. 0) .or. (which .eq. 1) ) then
          call lev%R(m)%copy(lev%I(m), 1)
          call lev%R(m)%axpy(1.0_pfdp, lev%Q(1), 1)
          call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m+1), 1)
       end if
       if( (which .eq. 0) .or. (which .eq. 2) ) then
          call lev%R(m)%copy(lev%I(m), 2)
          call lev%R(m)%axpy(1.0_pfdp,  lev%Q(lev%nnodes), 2)
          call lev%R(m)%axpy(-1.0_pfdp, lev%Q(m), 2)
       end if
    end do

  end subroutine imexQ_oc_residual

  subroutine imexQ_oc_spreadq0(this, lev, t0, flags, step)
    class(pf_imexQ_oc_t),  intent(inout) :: this
    class(pf_level_t),     intent(inout) :: lev
    real(pfdp),            intent(in   ) :: t0
    integer,     optional, intent(in)    :: flags, step

    integer :: m, p, which, mystep
    which = 3
    if(present(flags)) which = flags
    mystep = 1
    if(present(step))  then
      mystep = step !needed for sequential version
    else
      print *, "step not present in spreadq0", which
      stop
    end if
    
    select case(which)
      case(1)
        !  Stick initial condition into first node slot
        call lev%Q(1)%copy(lev%q0, 1)
        !  Evaluate F at first spot
        call lev%ulevel%sweeper%evaluate(lev, t0, 1, 1, mystep)
        ! Spread F and solution to all nodes
        do m = 2, lev%nnodes
          call lev%Q(m)%copy(lev%Q(1), 1)
          do p = 1, lev%ulevel%sweeper%npieces
            call lev%F(m,p)%copy(lev%F(1,p), 1)
          end do
        end do
      case(2)
        !  Stick terminal condition into last node slot
        call lev%Q(lev%nnodes)%copy(lev%qend, 2)
        !  Evaluate F at first spot
        call lev%ulevel%sweeper%evaluate(lev, t0, lev%nnodes, 2, mystep)
        ! Spread F and solution to all nodes
        do m = lev%nnodes-1, 1, -1
          call lev%Q(m)%copy(lev%Q(lev%nnodes), 2)
          do p = 1, lev%ulevel%sweeper%npieces
            call lev%F(m,p)%copy(lev%F(lev%nnodes,p), 2)
          end do
        end do
      case default
        print *, "imexQ_oc_spreadq0 needs a flag", which
        stop
    end select
    
  end subroutine imexQ_oc_spreadq0
  
  subroutine imexQ_oc_destroy(this, lev)
    !>  deallocate
    class(pf_imexQ_oc_t),  intent(inout) :: this
    class(pf_level_t), intent(inout) :: lev
    
    deallocate(this%QdiffE)
    deallocate(this%QdiffI)
    deallocate(this%QtilE)
    deallocate(this%QtilI)
    deallocate(this%dtsdc)

    call lev%ulevel%factory%destroy_single(this%rhs, lev%index, lev%shape)
  end subroutine imexQ_oc_destroy

end module pf_mod_imexQ_oc
