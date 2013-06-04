!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

program fpfasst
  use pf_mod_dtype
  use pf_mod_pfasst
  use pf_mod_parallel
  use pf_mod_imex
  use pf_mod_implicit
  use pf_mod_mpi
  use pf_mod_comm_mpi

  use encap
  use feval
  use hooks
  use probin
  use solutions
  use transfer

  implicit none

  type(pf_pfasst_t) :: pf
  type(pf_comm_t)   :: comm
  type(array1d)     :: q1

  type(pf_imex_t),     target :: imex1
  type(pf_implicit_t), target :: implicit1

  integer        :: ierror, l
  character(256) :: probin_fname


  !
  ! read options
  !

  if (command_argument_count() == 1) then
     call get_command_argument(1, value=probin_fname)
  else
     probin_fname = "probin.nml"
  end if
  call probin_init(probin_fname)


  !
  ! initialize mpi
  !

  call mpi_init(ierror)
  if (ierror .ne. 0) &
       stop "ERROR: Can't initialize MPI."


  !
  ! initialize pfasst
  !

  call pf_mpi_create(comm, MPI_COMM_WORLD)
  call pf_imex_create(imex1, f1eval1, f2eval1, f2comp1)
  call pf_implicit_create(implicit1, f2eval1, f2comp1)
  call pf_pfasst_create(pf, comm, nlevs)

  pf%niters = niters
  pf%qtype  = SDC_GAUSS_LOBATTO ! + SDC_PROPER_NODES

  pf%echo_timings = .true.

  if (nlevs > 1) then
     pf%levels(1)%nsweeps = 2
  end if

  do l = 1, nlevs
     pf%levels(l)%nvars  = nvars(l)
     pf%levels(l)%nnodes = nnodes(l)

     call feval_create_workspace(pf%levels(l)%ctx, pf%levels(l)%nvars)

     pf%levels(l)%interpolate => interpolate
     pf%levels(l)%restrict    => restrict
     pf%levels(l)%new         => encap_new
     if (problem == PROB_HEAT) then
        pf%levels(l)%sweeper  => implicit1
     else
        pf%levels(l)%sweeper  => imex1
     end if
  end do

  call pf_mpi_setup(comm, pf)
  call pf_pfasst_setup(pf)

  if (pf%rank == 0) then
     print *, 'nvars: ', pf%levels(:)%nvars
     print *, 'nnodes:', pf%levels(:)%nnodes
  end if


  !
  ! run
  !
  allocate(q1%array(pf%levels(nlevs)%nvars))
  call initial(q1)

  if (problem == PROB_AD) then
     call pf_add_hook(pf, nlevs, PF_POST_ITERATION, echo_error_hook)
  end if

  if (len_trim(output) > 0) then
     call dump_mkdir(output, len_trim(output))
     call pf_add_hook(pf, nlevs, PF_POST_SWEEP, dump_hook)
  end if

  if (nsteps < comm%nproc) then
     nsteps = comm%nproc
  end if

  call pf_pfasst_run(pf, q1, dt, 0.0_pfdp, nsteps=nsteps)


  !
  ! cleanup
  !
  deallocate(q1%array)

  do l = 1, nlevs
     call feval_destroy_workspace(pf%levels(l)%ctx)
  end do

  call pf_pfasst_destroy(pf)
  call pf_mpi_destroy(comm)
  call mpi_finalize(ierror)
  call fftw_cleanup()

end program fpfasst