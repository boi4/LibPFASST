!
! This file is part of LIBPFASST.
!
!>  Module for reading parameters for the problem
module probin
  use pf_mod_dtype


  character(len=64), save :: problem_type

  complex(pfdp), save :: lam

  real(pfdp), save :: dt     ! time step
  real(pfdp), save :: Tfin   ! Final time

  integer, save :: nsteps          ! number of time steps

  character(len=128), save :: pfasst_nml
  character(len=64), save :: output ! directory name for output

  namelist /params/  lam, dt, Tfin, nsteps, pfasst_nml, output

contains

  subroutine probin_init(pf_fname)
    character(len=*), intent(inout) :: pf_fname
    integer :: i   !  loop variable
    integer :: un  !  file read unit
    character(len=128) :: arg  !  command line argument
    character(128)    :: probin_fname   !<  file name for input parameters
    character(len=256) :: istring           ! stores command line argument
    character(len=1024) :: message          ! use for I/O error messages
    complex(pfdp)      ::   lam
    integer :: ios,iostat
    
    !> Set the name of the input file
    probin_fname = "probin.nml" ! default file name - can be overwritten on the command line
    if (command_argument_count() >= 1) &
         call get_command_argument(1, value=probin_fname)
    
    !> set defaults
    nsteps  = -1
    dt      = 0.01_pfdp
    Tfin    = 0.0_pfdp
    lam    = cmplx(1.0, 0.0, pfdp)
    pfasst_nml=probin_fname

    !>  Read in stuff from input file
    un = 9
    write(*,*) 'opening file ',TRIM(probin_fname), '  for input'
    open(unit=un, file = probin_fname, status = 'old', action = 'read')
    read(unit=un, nml = params)
    close(unit=un)
          
    !>  Read the command line
    i = 0
    do
       call get_command_argument(i, arg)
       if (LEN_TRIM(arg) == 0) EXIT
       if (i > 0) then
          istring="&PARAMS "//TRIM(arg)//" /"    
          read(istring,nml=params,iostat=ios,iomsg=message) ! internal read of NAMELIST
       end if
       i = i+1
    end do

    !  Reset dt if Tfin is set
    if (Tfin .gt. 0.0) dt = Tfin/dble(nsteps)

    !  Return the name of the file from which to read PFASST parameters
    pf_fname=pfasst_nml

  end subroutine probin_init

  subroutine print_loc_options(pf, un_opt)
    type(pf_pfasst_t), intent(inout)           :: pf   
    integer,           intent(in   ), optional :: un_opt
    integer :: un = 6

    if (pf%rank /= 0) return
    if (present(un_opt)) un = un_opt

    !  Print out the local parameters
    write(un,*) '=================================================='
    write(un,*) ' '
    write(un,*) 'Local Variables'
    write(un,*) '----------------'
    write(un,*) 'nsteps: ', nsteps, '! Number of steps'
    write(un,*) 'Dt:     ', Dt, '! Time step size'
    write(un,*) 'Tfin:   ', Tfin,   '! Final time of run'
    write(un,*) 'lam:   ', lam, '! linear constant'

    write(un,*) 'PFASST parameters read from input file ', pfasst_nml
    write(un,*) '=================================================='
  end subroutine print_loc_options
  

end module probin
