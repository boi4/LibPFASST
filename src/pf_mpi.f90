!
! Copyright (c) 2012, Matthew Emmett and Michael Minion.  All rights reserved.
!

! This module implements MPI communications.

module pf_mod_mpi
  include "mpif.h"
end module pf_mod_mpi

module pf_mod_comm
  use encap
  use pf_mod_mpi
  use pf_mod_dtype
  use pf_mod_timer

  implicit none

  interface create
     module procedure pf_comm_create
  end interface create

  interface setup
     module procedure pf_comm_setup
  end interface setup

  interface destroy
     module procedure pf_comm_destroy
  end interface destroy

contains

  ! Create an MPI based PFASST communicator
  subroutine pf_comm_create(pf_comm, mpi_comm)
    type(pf_comm_t), intent(out) :: pf_comm
    integer,         intent(in)  :: mpi_comm

    integer :: ierror

    pf_comm%comm = mpi_comm
    call mpi_comm_size(mpi_comm, pf_comm%nproc, ierror)
  end subroutine pf_comm_create

  ! Setup
  subroutine pf_comm_setup(pf_comm, pf)
    use pf_mod_mpi, only: MPI_REQUEST_NULL

    type(pf_comm_t),   intent(inout) :: pf_comm
    type(pf_pfasst_t), intent(inout) :: pf

    integer :: ierror

    call mpi_comm_rank(pf_comm%comm, pf%rank, ierror)

    allocate(pf_comm%recvreq(pf%nlevels))
    allocate(pf_comm%sendreq(pf%nlevels))

    pf_comm%sendreq = MPI_REQUEST_NULL
  end subroutine pf_comm_setup

  ! Destroy
  subroutine pf_comm_destroy(pf_comm)
    type(pf_comm_t), intent(inout) :: pf_comm

    deallocate(pf_comm%recvreq)
    deallocate(pf_comm%sendreq)
  end subroutine pf_comm_destroy

  ! Post
  subroutine post(pf, level, tag)
    use pf_mod_mpi, only: MPI_REAL8

    type(pf_pfasst_t), intent(in)    :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag

    integer :: ierror

    if (pf%comm%nproc > 1 .and. pf%rank > 0) then
       call mpi_irecv(level%recv, level%nvars, MPI_REAL8, &
            pf%rank-1, tag, pf%comm%comm, pf%comm%recvreq(level%level), ierror)
    end if
  end subroutine post

  ! Receive
  subroutine recv(pf, level, tag, blocking)
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking

    integer :: ierror, stat(MPI_STATUS_SIZE)

    call start_timer(pf, TRECEIVE + level%level - 1)

    if (pf%rank > 0) then
       if (blocking) then
          call mpi_recv(level%recv, level%nvars, MPI_REAL8, &
               pf%rank-1, tag, pf%comm%comm, stat, ierror)
       else
          call mpi_wait(pf%comm%recvreq(level%level), stat, ierror)
       end if

       if (ierror .ne. 0) then
          print *, 'WARNING: MPI ERROR DURING RECEIVE', ierror
       end if

       level%q0 = level%recv
    end if

    call end_timer(pf, TRECEIVE + level%level - 1)
  end subroutine recv

  ! Send
  subroutine send(pf, level, tag, blocking)
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(inout) :: pf
    type(pf_level_t),  intent(inout) :: level
    integer,           intent(in)    :: tag
    logical,           intent(in)    :: blocking

    integer :: ierror, stat(MPI_STATUS_SIZE)

    call start_timer(pf, TSEND + level%level - 1)

    if (pf%rank < pf%comm%nproc-1) then

       if (blocking) then
          call pack(level%send, level%qend)
          call mpi_send(level%send, level%nvars, MPI_REAL8, &
               pf%rank+1, tag, pf%comm%comm, stat, ierror)
       else
          call mpi_wait(pf%comm%sendreq(level%level), stat, ierror)
          call pack(level%send, level%qend)
          call mpi_isend(level%send, level%nvars, MPI_REAL8, &
               pf%rank+1, tag, pf%comm%comm, pf%comm%sendreq(level%level), ierror)
       end if
    end if

    call end_timer(pf, TSEND + level%level - 1)
  end subroutine send

  ! Send
  subroutine wait(pf, level)
    use pf_mod_mpi, only: MPI_REAL8, MPI_STATUS_SIZE

    type(pf_pfasst_t), intent(in) :: pf
    integer,           intent(in) :: level

    integer :: ierror, stat(MPI_STATUS_SIZE)

    call mpi_wait(pf%comm%sendreq(level), stat, ierror)
  end subroutine wait


  ! Broadcast
  subroutine broadcast(pf, y, nvar, root)
    use pf_mod_mpi, only: MPI_REAL8

    type(pf_pfasst_t), intent(inout) :: pf
    real(pfdp)  ,      intent(in)    :: y(nvar)
    integer,           intent(in)    :: nvar, root

    integer :: ierror

    call start_timer(pf, TSEND)

    call mpi_bcast(y, nvar, MPI_REAL8, root, pf%comm%comm, ierror)

    call end_timer(pf, TSEND)
  end subroutine broadcast

end module pf_mod_comm
