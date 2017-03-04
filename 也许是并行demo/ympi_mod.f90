!###############################################################################
!-------------------------------------------------------------------------------

  module ympi_mod

!-------------------------------------------------------------------------------
! ******************************************************************************
!-------------------------------------------------------------------------------
!                                                Copyright (C) 2005 Xunqiang Yin
!                                                MODULE NAME : ympi_mod
!                                                Current VERSION : 2011/10/19
!
! --- USAGE : Tools for parallel computation using MPI.
! --- DEPEND: netcdf_mod developed by Xunqiang Yin.
!
! --- NOTE for describing of subroutine / function :
!  A. The parameters bracketed with [], means optional parameter.
!  B. The describe for the parameters of subroutine / function, started with:
!   * It means input prameter;
!   # It means output prameter;
!   @ It means input and output prameter(it will be changed inside).
!
!-------------------------------------------------------------------------------
! ******************************************************************************
! ***                       INTERFACE DESCRIBE                               ***
! ******************************************************************************
!-------------------------------------------------------------------------------
!
!  1. subroutine ympi_init
!     --- Initialize this module.
!
!-------------------------------------------------------------------------------
!
!  2. subroutine ympi_final
!     --- Finalize this module.
!
!-------------------------------------------------------------------------------
!
!  3. subroutine init_pebox(filename, xname, yname, mname, n, myid)
!     * character(len=*) :: filename = the name of grid file.
!     * character(len=*) :: xname = the name of x-direction.
!     * character(len=*) :: yname = the name of y-direction.
!     * character(len=*) :: mname = the name of mask arrary.
!     * integer          :: n = number of all PEs.
!     * integer          :: myid = index of the current PE.
!
!-------------------------------------------------------------------------------
!
!  4. subroutine opennc(ncid, filename, action)
!     # integer :: ncid = unit number of opened netcdf file.
!     * character(len=*) :: filename = the name of netcdf file.
!     * character :: action = 'create', 'define', 'write' or 'read'.
!
!-------------------------------------------------------------------------------
!
!  5. subroutine closenc(ncid)
!     * integer :: ncid = unit number of opened netcdf file.
!
!-------------------------------------------------------------------------------
!
!  6. subroutine get_mpimin(x)
!     @ integer/real/double precision :: x = the value need to get min.
!
!-------------------------------------------------------------------------------
!
!  7. subroutine bcast_to_all(x)
!     @ integer/real/double precision :: x = the value need to be broadcasted.
!
!-------------------------------------------------------------------------------
!
!  8. subroutine updatev(ee, ixl, iyl[, kl[, jl]], halo[, flag])
!     @ real(4/8) :: ee = the array need to be updated, could be 2d,3d & 4d.
!     * integer :: ixl, iyl, kl, jl = is for the shape of ee.
!     * integer :: halo = the size of halo in spatial partition.
!     * flag = the flag for 4d direction of 4d array.
!       --- If it is 4d & flag is given, the direction of dimensions is 
!           (ixl, iyl, kl, jl), otherwise, it is (kl, jl, ixl, iyl)
!       --- For 3d array, it is defined as ee(ixl, iyl, kl)
!
!-------------------------------------------------------------------------------
!
! --- Useful variables:
!
!     lon, lat, im, jm, halosize, flagxcycle, ixoverlay
!     myid, npe, mypebox, loc2, loc3, loc4, is, ie, js, je
!
!-------------------------------------------------------------------------------
!
!                                                --- Xunqiang Yin, 2011/10/19
!                                                 E-Mail: XunqiangYin@gmail.com
!
!-------------------------------------------------------------------------------
! ******************************************************************************
!-------------------------------------------------------------------------------

  use netcdf_mod

  implicit none

!-------------------------------------------------------------------------------

  public ympi_init, ympi_final, init_pebox
  public get_mpimin, bcast_to_all, updatev, opennc, closenc, runbyturn
  public lon, lat, im, jm, halosize, flagxcycle, ixoverlay
  public myid, npe, mypebox, loc2, loc3, loc4, is, ie, js, je
!  private

!-------------------------------------------------------------------------------

  include 'mpif.h'
  integer :: myid, npe, mpi_comm_ympi
  integer :: loc2(2), loc3(3), loc4(4), is, ie, js, je
  integer, allocatable :: qrs(:), qrr(:), qls(:), qlr(:), tgl(:), tgr(:)
  integer, allocatable :: statr(:, :), statl(:, :)
  integer :: im, jm, flagxcycle, halosize, ixoverlay
  real(4), allocatable :: lon(:), lat(:), mask(:, :)

!-------------------------------------------------------------------------------

  type pebox_type
    integer :: myid, i1, i2, j1, j2   ! 1,2,3,4,5
    integer :: up, dn, nr, nl, np     ! 7, 6, 8, 9, 10
    integer :: ei1, ei2, ej1, ej2     ! 11, 12, 13, 14
    integer :: halo, isize, jsize     ! 15, 16
    integer, pointer :: left(:, :)    ! (nl, 5): id, rj1, rj2, sj1, sj2
    integer, pointer :: right(:, :)   ! (nr, 5): id, rj1, rj2, sj1, sj2
  end type pebox_type
  type(pebox_type) :: mypebox
  
!-------------------------------------------------------------------------------

  interface updatev
    module procedure update2d_sgl, update3d_sgl, update4d_sgl, update4d_sgle, &
                     update2d_dbl, update3d_dbl, update4d_dbl, update4d_dble
  end interface updatev

  interface get_mpimin
    module procedure get_mpimini, get_mpiminr, get_mpimind
  end interface get_mpimin

  interface bcast_to_all
    module procedure bcast_to_alli, bcast_to_allr, bcast_to_alld
  end interface bcast_to_all
  
!-------------------------------------------------------------------------------

  contains

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: runbyturn

  subroutine runbyturn(flag)
    integer, intent(in) :: flag
    integer :: ierr, i = 0
    integer :: stat(mpi_status_size)
    if(flag == 0 .and. myid /= 0)then
      ! recv message from myid-1, tag=229
      call mpi_recv(i, 1, mpi_integer, myid-1, 229, mpi_comm_ympi, stat, ierr)
    endif
    if(flag == 1 .and. myid < npe-1)then
      ! send message to myid+1, tag=229
      call mpi_send(i, 1, mpi_integer, myid+1, 229, mpi_comm_ympi, ierr)
    endif
  end subroutine runbyturn

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: opennc

  subroutine opennc(ncid, filename, action)
    character(len=*), intent(in) :: filename, action
    integer, intent(in) :: ncid
    integer :: flag, ierr
    integer :: stat(mpi_status_size)
    if(myid /= 0)then
      ! recv message from myid-1, tag=229
      call mpi_recv(flag, 1, mpi_integer, myid-1, 229, mpi_comm_ympi, stat, ierr)
    endif
    call open_nc(ncid, filename, action)
  end subroutine opennc

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: closenc

  subroutine closenc(ncid)
    integer, intent(in) :: ncid
    integer :: flag=0, ierr
    integer :: stat(mpi_status_size)
    call close_nc(ncid)
    if(myid < npe-1)then
      ! send message to myid+1, tag=229
      call mpi_send(flag, 1, mpi_integer, myid+1, 229, mpi_comm_ympi, ierr           )
    endif
  end subroutine closenc

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: ympi_init

  subroutine ympi_init
    integer :: ierr
    character(len=10) :: str
    ! --- Initial mpi, get myid & npe
    call mpi_init(ierr)
    mpi_comm_ympi = mpi_comm_world
    call mpi_comm_rank(mpi_comm_ympi, myid, ierr)
    call mpi_comm_size(mpi_comm_ympi, npe, ierr)
  end subroutine ympi_init

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: ympi_final

  subroutine ympi_final
    integer :: rc
    close(6)
    call mpi_finalize(rc); stop
  end subroutine ympi_final

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: get_mpimini

  subroutine get_mpimini(x)
    integer, intent(inout) :: x
    integer :: y
    integer :: ierr
    y = x
    call mpi_reduce(y, x, 1, mpi_integer, mpi_min, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(x, 1, mpi_integer, 0, mpi_comm_ympi, ierr)
  end subroutine get_mpimini
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: get_mpiminr

  subroutine get_mpiminr(x)
    real(4), intent(inout) :: x
    real(4) :: y
    integer :: ierr
    y = x
    call mpi_reduce(y, x, 1, mpi_real, mpi_min, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(x, 1, mpi_real, 0, mpi_comm_ympi, ierr)
  end subroutine get_mpiminr
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: get_mpimind

  subroutine get_mpimind(x)
    real(8), intent(inout) :: x
    real(8) :: y
    integer :: ierr
    y = x
    call mpi_reduce(y, x, 1, mpi_double_precision, mpi_min, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(x, 1, mpi_double_precision, 0, mpi_comm_ympi, ierr)
  end subroutine get_mpimind
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: bcast_to_alli

  subroutine bcast_to_alli(ii)
    integer, intent(inout) :: ii
    integer :: ierr
    call mpi_bcast(ii, 1, mpi_integer, 0, mpi_comm_ympi, ierr)
  end subroutine bcast_to_alli
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: bcast_to_allr

  subroutine bcast_to_allr(ii)
    real(4), intent(inout) :: ii
    integer :: ierr
    call mpi_bcast(ii, 1, mpi_real, 0, mpi_comm_ympi, ierr)
  end subroutine bcast_to_allr
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: bcast_to_alld

  subroutine bcast_to_alld(ii)
    real(8), intent(inout) :: ii
    integer :: ierr
    call mpi_bcast(ii, 1, mpi_double_precision, 0, mpi_comm_ympi, ierr)
  end subroutine bcast_to_alld
  
!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update2d_sgl

  subroutine update2d_sgl(ee, ixl, iyl, halo)
    integer, intent(in) :: ixl, iyl, halo
    real(4), intent(inout) :: ee(ixl, iyl)
    integer :: ncount, i, ierr
    integer :: stat(mpi_status_size)
    real(4), dimension(ixl, halo) :: dnee, dnee1, upee, upee1
    real(4), dimension(halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    call mpi_barrier(mpi_comm_ympi, ierr)
    dnee1(:, :) = ee(:, halosize+1:halosize*2)
    upee1(:, :) = ee(:, iyl-halosize*2+1:iyl-halosize)
    ncount = ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, mpi_real, mypebox%up, 100, dnee , ncount, mpi_real, mypebox%dn, 100, mpi_comm_ympi, stat, ierr)
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, mpi_real, mypebox%dn, 200, upee , ncount, mpi_real, mypebox%up, 200, mpi_comm_ympi, stat, ierr)
    if(mypebox%dn /= MPI_PROC_NULL)ee(:, 1:halosize) = dnee
    if(mypebox%up /= MPI_PROC_NULL)ee(:, iyl-halosize+1:iyl) = upee
    lfee1(:, :) = ee(halosize+1:halosize*2, :)
    rtee1(:, :) = ee(ixl-halosize*2+1:ixl-halosize, :)
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, mypebox%left(i, 4)), ncount, mpi_real, mypebox%left(i, 1), 101, mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, mypebox%left(i, 4)), ncount, mpi_real, mypebox%left(i, 1), 201, mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr, 1 !myrightno
      ! --- Send to right.
      ncount = halosize*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, mypebox%right(i, 4)), ncount, mpi_real, mypebox%right(i, 1), 201, mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv( rtee(1, mypebox%right(i, 4)), ncount, mpi_real, mypebox%right(i, 1), 101, mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      ee(ixl-halosize+1:ixl, :) = rtee
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      ee(1:halosize, :) = lfee
    endif
  end subroutine update2d_sgl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update3d_sgl

  subroutine update3d_sgl(ee, ixl, iyl, kl, halo)
    integer, intent(in) :: ixl, iyl, kl, halo
    real(4), intent(inout) :: ee(ixl, iyl, kl)
    integer :: ncount, i, ierr, j, k
    integer :: stat(mpi_status_size)
    real(4), dimension(kl, ixl, halo) :: dnee, dnee1, upee, upee1
    real(4), dimension(kl, halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    do k = 1, kl
    do i = 1, ixl
      dnee1(k, i, :) = ee(i, halosize+1:halosize*2, k)
      upee1(k, i, :) = ee(i, iyl-halosize*2+1:iyl-halosize, k)
    enddo
    enddo
    ncount = kl*ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, mpi_real, mypebox%up, 100, dnee , ncount, mpi_real, mypebox%dn, 100, mpi_comm_ympi, stat, ierr                )
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, mpi_real, mypebox%dn, 200, upee , ncount, mpi_real, mypebox%up, 200, mpi_comm_ympi, stat, ierr                )
    if(mypebox%dn /= MPI_PROC_NULL .and. mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do i = 1, ixl
        ee(i, 1:halosize        , k) = dnee(k, i, :)
        ee(i, iyl-halosize+1:iyl, k) = upee(k, i, :)
      enddo
      enddo
    elseif(mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do i = 1, ixl
        ee(i, iyl-halosize+1:iyl, k) = upee(k, i, :)
      enddo
      enddo
    elseif(mypebox%dn /= MPI_PROC_NULL)then
      do k = 1, kl
      do i = 1, ixl
        ee(i, 1:halosize        , k) = dnee(k, i, :)
      enddo
      enddo
    endif
    do k = 1, kl
    do i = 1, iyl
      lfee1(k, :, i) = ee(halosize+1:halosize*2        , i, k)
      rtee1(k, :, i) = ee(ixl-halosize*2+1:ixl-halosize, i, k)
    enddo
    enddo
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*kl*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, 1, mypebox%left(i, 4)),ncount, mpi_real, mypebox%left(i, 1), 101,mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, 1, mypebox%left(i, 4)),ncount, mpi_real, mypebox%left(i, 1), 201,mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr !myrightno
      ! --- Send to right.
      ncount = halosize*kl*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, 1, mypebox%right(i, 4)),ncount, mpi_real, mypebox%right(i, 1), 201,mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv(rtee(1, 1, mypebox%right(i, 4)),ncount, mpi_real, mypebox%right(i, 1), 101,mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      do k = 1, kl
      do i = 1, iyl
        ee(ixl-halosize+1:ixl, i, k) = rtee(k, :, i)
      enddo
      enddo
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      do k = 1, kl
      do i = 1, iyl
        ee(1:halosize, i, k) = lfee(k, :, i)
      enddo
      enddo
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
    endif
  end subroutine update3d_sgl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update4d_sgle

  subroutine update4d_sgle(ee, ixl, iyl, kl, jl, halo)
    integer, intent(in) :: ixl, iyl, kl, jl, halo
    real(4), intent(inout) :: ee(kl, jl, ixl, iyl)
    integer :: ncount, i, ierr
    integer :: stat(mpi_status_size)
    real(4), dimension(kl, jl, ixl, halo) :: dnee, dnee1, upee, upee1
    real(4), dimension(kl, jl, halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    dnee1(:, :, :, :) = ee(:, :, :, halosize+1:halosize*2)
    upee1(:, :, :, :) = ee(:, :, :, iyl-halosize*2+1:iyl-halosize)
    ncount = kl*jl*ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, mpi_real, mypebox%up, 100, dnee , ncount, mpi_real, mypebox%dn, 100, mpi_comm_ympi, stat, ierr                )
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, mpi_real, mypebox%dn, 200, upee , ncount, mpi_real, mypebox%up, 200, mpi_comm_ympi, stat, ierr                )
    if(mypebox%dn /= MPI_PROC_NULL)ee(:, :, :, 1:halosize) = dnee
    if(mypebox%up /= MPI_PROC_NULL)ee(:, :, :, iyl-halosize+1:iyl) = upee
    lfee1(:, :, :, :) = ee(:, :, halosize+1:halosize*2, :)
    rtee1(:, :, :, :) = ee(:, :, ixl-halosize*2+1:ixl-halosize, :)
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*kl*jl*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, 1, 1, mypebox%left(i, 4)),ncount, mpi_real, mypebox%left(i, 1), 101,mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, 1, 1, mypebox%left(i, 4)),ncount, mpi_real, mypebox%left(i, 1), 201,mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr !myrightno
      ! --- Send to right.
      ncount = halosize*kl*jl*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, 1, 1, mypebox%right(i, 4)),ncount, mpi_real, mypebox%right(i, 1), 201,mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv(rtee(1, 1, 1, mypebox%right(i, 4)),ncount, mpi_real, mypebox%right(i, 1), 101,mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      ee(:, :, ixl-halosize+1:ixl, :) = rtee
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      ee(:, :, 1:halosize, :) = lfee
    endif
  end subroutine update4d_sgle

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update4d_sgl

  subroutine update4d_sgl(ee, ixl, iyl, kl, jl, halo, flag)
    integer, intent(in) :: ixl, iyl, kl, jl, halo, flag
    real(4), intent(inout) :: ee(ixl, iyl, kl, jl)
    integer :: ncount, i, ierr, j, k
    integer :: stat(mpi_status_size)
    real(4), dimension(kl, jl, ixl, halo) :: dnee, dnee1, upee, upee1
    real(4), dimension(kl, jl, halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    do k = 1, kl
    do j = 1, jl
    do i = 1, ixl
      dnee1(k, j, i, :) = ee(i, halosize+1:halosize*2, k, j)
      upee1(k, j, i, :) = ee(i, iyl-halosize*2+1:iyl-halosize, k, j)
    enddo
    enddo
    enddo
    ncount = kl*jl*ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, mpi_real, mypebox%up, 100, dnee , ncount, mpi_real, mypebox%dn, 100, mpi_comm_ympi, stat, ierr                )
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, mpi_real, mypebox%dn, 200, upee , ncount, mpi_real, mypebox%up, 200, mpi_comm_ympi, stat, ierr                )
    if(mypebox%dn /= MPI_PROC_NULL .and. mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do j = 1, jl
      do i = 1, ixl
        ee(i, 1:halosize        , k, j) = dnee(k, j, i, :)
        ee(i, iyl-halosize+1:iyl, k, j) = upee(k, j, i, :)
      enddo
      enddo
      enddo
    elseif(mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do j = 1, jl
      do i = 1, ixl
        ee(i, iyl-halosize+1:iyl, k, j) = upee(k, j, i, :)
      enddo
      enddo
      enddo
    elseif(mypebox%dn /= MPI_PROC_NULL)then
      do k = 1, kl
      do j = 1, jl
      do i = 1, ixl
        ee(i, 1:halosize        , k, j) = dnee(k, j, i, :)
      enddo
      enddo
      enddo
    endif
    do k = 1, kl
    do j = 1, jl
    do i = 1, iyl
      lfee1(k, j, :, i) = ee(halosize+1:halosize*2        , i, k, j)
      rtee1(k, j, :, i) = ee(iyl-halosize*2+1:iyl-halosize, i, k, j)
    enddo
    enddo
    enddo
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*kl*jl*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, 1, 1, mypebox%left(i, 4)),ncount, mpi_real, mypebox%left(i, 1), 101,mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, 1, 1, mypebox%left(i, 4)),ncount, mpi_real, mypebox%left(i, 1), 201,mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr !myrightno
      ! --- Send to right.
      ncount = halosize*kl*jl*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, 1, 1, mypebox%right(i, 4)),ncount, mpi_real, mypebox%right(i, 1), 201,mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv(rtee(1, 1, 1, mypebox%right(i, 4)),ncount, mpi_real, mypebox%right(i, 1), 101,mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      do k = 1, kl
      do j = 1, jl
      do i = 1, iyl
        ee(ixl-halosize+1:ixl, i, k, j) = rtee(k, j, :, i)
      enddo
      enddo
      enddo
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      do k = 1, kl
      do j = 1, jl
      do i = 1, iyl
        ee(1:halosize, i, k, j) = lfee(k, j, :, i)
      enddo
      enddo
      enddo
    endif
  end subroutine update4d_sgl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update2d_dbl

  subroutine update2d_dbl(ee, ixl, iyl, halo)
    integer, intent(in) :: ixl, iyl, halo
    real(8), intent(inout) :: ee(ixl, iyl)
    integer :: ncount, i, ierr
    integer :: stat(mpi_status_size)
    real(8), dimension(ixl, halo) :: dnee, dnee1, upee, upee1
    real(8), dimension(halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    call mpi_barrier(mpi_comm_ympi, ierr)
    dnee1(:, :) = ee(:, halosize+1:halosize*2)
    upee1(:, :) = ee(:, iyl-halosize*2+1:iyl-halosize)
    ncount = ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, MPI_DOUBLE_PRECISION, mypebox%up, 100, dnee , ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 100, mpi_comm_ympi, stat, ierr)
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 200, upee , ncount, MPI_DOUBLE_PRECISION, mypebox%up, 200, mpi_comm_ympi, stat, ierr)
    if(mypebox%dn /= MPI_PROC_NULL)ee(:, 1:halosize) = dnee
    if(mypebox%up /= MPI_PROC_NULL)ee(:, iyl-halosize+1:iyl) = upee
    lfee1(:, :) = ee(halosize+1:halosize*2, :)
    rtee1(:, :) = ee(ixl-halosize*2+1:ixl-halosize, :)
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, mypebox%left(i, 4)), ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 101, mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, mypebox%left(i, 4)), ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 201, mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr, 1 !myrightno
      ! --- Send to right.
      ncount = halosize*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, mypebox%right(i, 4)), ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 201, mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv( rtee(1, mypebox%right(i, 4)), ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 101, mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      ee(ixl-halosize+1:ixl, :) = rtee
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      ee(1:halosize, :) = lfee
    endif
  end subroutine update2d_dbl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update3d_dbl

  subroutine update3d_dbl(ee, ixl, iyl, kl, halo)
    integer, intent(in) :: ixl, iyl, kl, halo
    real(8), intent(inout) :: ee(ixl, iyl, kl)
    integer :: ncount, i, ierr, j, k
    integer :: stat(mpi_status_size)
    real(8), dimension(kl, ixl, halo) :: dnee, dnee1, upee, upee1
    real(8), dimension(kl, halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    do k = 1, kl
    do i = 1, ixl
      dnee1(k, i, :) = ee(i, halosize+1:halosize*2, k)
      upee1(k, i, :) = ee(i, iyl-halosize*2+1:iyl-halosize, k)
    enddo
    enddo
    ncount = kl*ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, MPI_DOUBLE_PRECISION, mypebox%up, 100, dnee , ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 100, mpi_comm_ympi, stat, ierr                )
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 200, upee , ncount, MPI_DOUBLE_PRECISION, mypebox%up, 200, mpi_comm_ympi, stat, ierr                )
    if(mypebox%dn /= MPI_PROC_NULL .and. mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do i = 1, ixl
        ee(i, 1:halosize        , k) = dnee(k, i, :)
        ee(i, iyl-halosize+1:iyl, k) = upee(k, i, :)
      enddo
      enddo
    elseif(mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do i = 1, ixl
        ee(i, iyl-halosize+1:iyl, k) = upee(k, i, :)
      enddo
      enddo
    elseif(mypebox%dn /= MPI_PROC_NULL)then
      do k = 1, kl
      do i = 1, ixl
        ee(i, 1:halosize        , k) = dnee(k, i, :)
      enddo
      enddo
    endif
    do k = 1, kl
    do i = 1, iyl
      lfee1(k, :, i) = ee(halosize+1:halosize*2        , i, k)
      rtee1(k, :, i) = ee(ixl-halosize*2+1:ixl-halosize, i, k)
    enddo
    enddo
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*kl*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, 1, mypebox%left(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 101,mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, 1, mypebox%left(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 201,mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr !myrightno
      ! --- Send to right.
      ncount = halosize*kl*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, 1, mypebox%right(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 201,mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv(rtee(1, 1, mypebox%right(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 101,mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      do k = 1, kl
      do i = 1, iyl
        ee(ixl-halosize+1:ixl, i, k) = rtee(k, :, i)
      enddo
      enddo
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      do k = 1, kl
      do i = 1, iyl
        ee(1:halosize, i, k) = lfee(k, :, i)
      enddo
      enddo
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
    endif
  end subroutine update3d_dbl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update4d_dble

  subroutine update4d_dble(ee, ixl, iyl, kl, jl, halo)
    integer, intent(in) :: ixl, iyl, kl, jl, halo
    real(8), intent(inout) :: ee(kl, jl, ixl, iyl)
    integer :: ncount, i, ierr
    integer :: stat(mpi_status_size)
    real(8), dimension(kl, jl, ixl, halo) :: dnee, dnee1, upee, upee1
    real(8), dimension(kl, jl, halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    dnee1(:, :, :, :) = ee(:, :, :, halosize+1:halosize*2)
    upee1(:, :, :, :) = ee(:, :, :, iyl-halosize*2+1:iyl-halosize)
    ncount = kl*jl*ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, MPI_DOUBLE_PRECISION, mypebox%up, 100, dnee , ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 100, mpi_comm_ympi, stat, ierr                )
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 200, upee , ncount, MPI_DOUBLE_PRECISION, mypebox%up, 200, mpi_comm_ympi, stat, ierr                )
    if(mypebox%dn /= MPI_PROC_NULL)ee(:, :, :, 1:halosize) = dnee
    if(mypebox%up /= MPI_PROC_NULL)ee(:, :, :, iyl-halosize+1:iyl) = upee
    lfee1(:, :, :, :) = ee(:, :, halosize+1:halosize*2, :)
    rtee1(:, :, :, :) = ee(:, :, ixl-halosize*2+1:ixl-halosize, :)
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*kl*jl*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, 1, 1, mypebox%left(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 101,mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, 1, 1, mypebox%left(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 201,mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr !myrightno
      ! --- Send to right.
      ncount = halosize*kl*jl*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, 1, 1, mypebox%right(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 201,mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv(rtee(1, 1, 1, mypebox%right(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 101,mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      ee(:, :, ixl-halosize+1:ixl, :) = rtee
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      ee(:, :, 1:halosize, :) = lfee
    endif
  end subroutine update4d_dble

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: update4d_dbl

  subroutine update4d_dbl(ee, ixl, iyl, kl, jl, halo, flag)
    integer, intent(in) :: ixl, iyl, kl, jl, halo, flag
    real(8), intent(inout) :: ee(ixl, iyl, kl, jl)
    integer :: ncount, i, ierr, j, k
    integer :: stat(mpi_status_size)
    real(8), dimension(kl, jl, ixl, halo) :: dnee, dnee1, upee, upee1
    real(8), dimension(kl, jl, halo, iyl) :: lfee, lfee1, rtee, rtee1
    if(halosize /= halo)then
      write(6, *)'ERROR: The halo size is incorrect.';stop
    endif
    do k = 1, kl
    do j = 1, jl
    do i = 1, ixl
      dnee1(k, j, i, :) = ee(i, halosize+1:halosize*2, k, j)
      upee1(k, j, i, :) = ee(i, iyl-halosize*2+1:iyl-halosize, k, j)
    enddo
    enddo
    enddo
    ncount = kl*jl*ixl*halosize
    ! --- From down to up.
    call mpi_sendrecv(upee1, ncount, MPI_DOUBLE_PRECISION, mypebox%up, 100, dnee , ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 100, mpi_comm_ympi, stat, ierr                )
    ! --- From up to down.
    call mpi_sendrecv(dnee1, ncount, MPI_DOUBLE_PRECISION, mypebox%dn, 200, upee , ncount, MPI_DOUBLE_PRECISION, mypebox%up, 200, mpi_comm_ympi, stat, ierr                )
    if(mypebox%dn /= MPI_PROC_NULL .and. mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do j = 1, jl
      do i = 1, ixl
        ee(i, 1:halosize        , k, j) = dnee(k, j, i, :)
        ee(i, iyl-halosize+1:iyl, k, j) = upee(k, j, i, :)
      enddo
      enddo
      enddo
    elseif(mypebox%up /= MPI_PROC_NULL)then
      do k = 1, kl
      do j = 1, jl
      do i = 1, ixl
        ee(i, iyl-halosize+1:iyl, k, j) = upee(k, j, i, :)
      enddo
      enddo
      enddo
    elseif(mypebox%dn /= MPI_PROC_NULL)then
      do k = 1, kl
      do j = 1, jl
      do i = 1, ixl
        ee(i, 1:halosize        , k, j) = dnee(k, j, i, :)
      enddo
      enddo
      enddo
    endif
    do k = 1, kl
    do j = 1, jl
    do i = 1, iyl
      lfee1(k, j, :, i) = ee(halosize+1:halosize*2        , i, k, j)
      rtee1(k, j, :, i) = ee(iyl-halosize*2+1:iyl-halosize, i, k, j)
    enddo
    enddo
    enddo
    do i = 1, mypebox%nl !myleftno
      ! --- Send to left.
      ncount = halosize*kl*jl*(mypebox%left(i, 5)-mypebox%left(i, 4)+1)
      call mpi_isend(lfee1(1, 1, 1, mypebox%left(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 101,mpi_comm_ympi, qls(i), ierr                         )
      ! --- Receive from left.
      call mpi_irecv(lfee(1, 1, 1, mypebox%left(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%left(i, 1), 201,mpi_comm_ympi, qlr(i), ierr                         )
    enddo
    do i = 1, mypebox%nr !myrightno
      ! --- Send to right.
      ncount = halosize*kl*jl*(mypebox%right(i, 5)-mypebox%right(i, 4)+1)
      call mpi_isend(rtee1(1, 1, 1, mypebox%right(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 201,mpi_comm_ympi, qrs(i), ierr                         )
      ! --- Receive from right.
      call mpi_irecv(rtee(1, 1, 1, mypebox%right(i, 4)),ncount, MPI_DOUBLE_PRECISION, mypebox%right(i, 1), 101,mpi_comm_ympi, qrr(i), ierr                      )
    enddo
    if(mypebox%nr > 0)then
      call mpi_waitall(mypebox%nr, qrs, statr(:, 1:mypebox%nr), ierr)
      call mpi_waitall(mypebox%nr, qrr, statr(:, 1:mypebox%nr), ierr)
      do k = 1, kl
      do j = 1, jl
      do i = 1, iyl
        ee(ixl-halosize+1:ixl, i, k, j) = rtee(k, j, :, i)
      enddo
      enddo
      enddo
    endif
    if(mypebox%nl > 0)then
      call mpi_waitall(mypebox%nl, qls, statl(:, 1:mypebox%nl), ierr)
      call mpi_waitall(mypebox%nl, qlr, statl(:, 1:mypebox%nl), ierr)
      do k = 1, kl
      do j = 1, jl
      do i = 1, iyl
        ee(1:halosize, i, k, j) = lfee(k, j, :, i)
      enddo
      enddo
      enddo
    endif
  end subroutine update4d_dbl

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: set_pebox !input parameters is needed before this sub.

  subroutine init_pebox(filename, xname, yname, mname, n, myid)
    integer, intent(in) :: n, myid
    character(len=*), intent(in) :: filename, xname, yname, mname
    integer :: ncid, i, ierr, im1
    integer, allocatable :: pebox(:, :), peleft(:, :, :), peright(:, :, :)
    allocate(pebox(17, n), peleft(n, n, 5), peright(n, n, 5))
    if(myid == 0)then
      ! --- Input lon, lat & mask.
      call open_nc(ncid, filename, 'r')
      im = get_dimension_len(ncid, xname)
      jm = get_dimension_len(ncid, yname)
      allocate(lon(im), lat(jm), mask(im, jm))
      call readnc(ncid, xname, lon)
      call readnc(ncid, yname, lat)
      call readnc(ncid, mname, mask)
      call close_nc(ncid)
      ! --- Check xcycle & ixoverlay for global model.
      flagxcycle = 1 ! --- For regional model.
      !if((lon(1) + 360. - lon(im)) < (lon(2) - lon(1)))flagxcycle = 0
      if((2*lon(1) + 360. - lon(im) - lon(2)) < 0)flagxcycle = 0 ! --- global
      ixoverlay = 0
      if(flagxcycle == 0)then
        do i = 1, im
          if(lon(i) > lon(im) - 360)then
            exit
          else
            ixoverlay = i
          endif
        enddo
      endif
      im1 = im - ixoverlay
      ! --- Set pebox for all processors.
      call set_pebox(im1, jm, n, flagxcycle, halosize, lon(1:im-1), lat, mask(1:im1, :), pebox, peleft, peright)
    endif
  ! --- Bcast pebox, peleft, peright & some parameters.
    call mpi_bcast(pebox  , 17*n, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(peleft , 5*n*n, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(peright, 5*n*n, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(im, 1, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(jm, 1, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(flagxcycle, 1, mpi_integer, 0, mpi_comm_ympi, ierr)
    call mpi_bcast(ixoverlay, 1, mpi_integer, 0, mpi_comm_ympi, ierr)
    mypebox%myid  = pebox(1, myid+1)
    mypebox%i1    = pebox(2, myid+1)
    mypebox%i2    = pebox(3, myid+1)
    mypebox%j1    = pebox(4, myid+1)
    mypebox%j2    = pebox(5, myid+1)
    mypebox%dn    = pebox(6, myid+1)
    mypebox%up    = pebox(7, myid+1)
    mypebox%nr    = pebox(8, myid+1)
    mypebox%nl    = pebox(9, myid+1)
    mypebox%np    = pebox(10, myid+1)
    mypebox%ei1   = pebox(11, myid+1)
    mypebox%ei2   = pebox(12, myid+1)
    mypebox%ej1   = pebox(13, myid+1)
    mypebox%ej2   = pebox(14, myid+1)
    mypebox%halo  = pebox(15, myid+1)
    mypebox%isize = pebox(16, myid+1)
    mypebox%jsize = pebox(17, myid+1)
    allocate(mypebox%left(mypebox%nl, 5), mypebox%right(mypebox%nr, 5))
    mypebox%left(:, :) = peleft(myid+1, 1:mypebox%nl, 1:5)
    mypebox%right(:, :) = peright(myid+1, 1:mypebox%nr, 1:5)
    loc2 = [mypebox%i1, mypebox%j1]
    loc3 = [mypebox%i1, mypebox%j1, 1]
    loc4 = [1, 1, mypebox%i1, mypebox%j1]
    is = mypebox%i1 - mypebox%ei1 + 1
    ie = mypebox%i2 - mypebox%ei1 + 1
    js = mypebox%j1 - mypebox%ej1 + 1
    je = mypebox%j2 - mypebox%ej1 + 1
    if(mypebox%nr>0)then
      allocate(qrs(mypebox%nr), qrr(mypebox%nr), tgr(mypebox%nr), statr(mpi_status_size, mypebox%nr))
      qrs = 0; qrr = 0;
    endif
    if(mypebox%nl>0)then
      allocate(qls(mypebox%nl), qlr(mypebox%nl), tgl(mypebox%nl), statl(mpi_status_size, mypebox%nl))
      qlr = 0; qls = 0
    endif
    if(myid == 0)then
      write(6, *)'myid  = ', mypebox%myid      
      write(6, *)'i1    = ', mypebox%i1        
      write(6, *)'i2    = ', mypebox%i2        
      write(6, *)'j1    = ', mypebox%j1        
      write(6, *)'j2    = ', mypebox%j2        
      write(6, *)'dn    = ', mypebox%dn        
      write(6, *)'up    = ', mypebox%up        
      write(6, *)'nr    = ', mypebox%nr        
      write(6, *)'nl    = ', mypebox%nl        
      write(6, *)'np    = ', mypebox%np        
      write(6, *)'ei1   = ', mypebox%ei1       
      write(6, *)'ei2   = ', mypebox%ei2       
      write(6, *)'ej1   = ', mypebox%ej1       
      write(6, *)'ej2   = ', mypebox%ej2       
      write(6, *)'halo  = ', mypebox%halo 
      write(6, *)'isize = ', mypebox%isize     
      write(6, *)'jsize = ', mypebox%jsize     
      write(6, "('left  (id, rj1, rj2, sj1, sj2) : ', 5i5)")mypebox%left(:, :)
      write(6, "('right (id, rj1, rj2, sj1, sj2) : ', 5i5)")mypebox%right(:, :)
!      print *,"CISTIME=",istime
    endif
    deallocate(pebox, peleft, peright)
  end subroutine init_pebox

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: set_pebox

  subroutine set_pebox(im, jm, n, flagxcycle, halosize, lon, lat, mask, pebox, peleft, peright)
    integer, intent(in) :: im, jm, n, flagxcycle, halosize
    real(4), intent(in) :: lon(im), lat(jm), mask(im, jm)
    integer, intent(out) :: pebox(17, n), peleft(n, n, 5), peright(n, n, 5)
    real(4) :: avenpx, avenpy, a, aa, avepoint
    integer, allocatable :: npx(:), npxstart(:), npxend(:)
    integer :: i1, i2, i, n1, n2, ii, j1, j2, j, jj, nn
    integer :: nr, nl, nx, ny, nxmore, ny1
    !  nx = 1; ny = n; aa = im*jm
    !  do i = 1, n/2
    !    j = n / i
    !    if(i * j /= n)cycle
    !    a = abs(i-j)
    !    if(aa > a)then
    !      aa = a; nx = i; ny = j
    !    endif
    !  enddo
    !
    !  if(im > jm .and. nx < ny)then
    !    i = nx; nx = ny; ny = i
    !  endif
    !  
    !  if(im < jm .and. nx > ny)then
    !    i = nx; nx = ny; ny = i
    !  endif
    nx = sqrt(n * im / float(jm))
    ny = n / nx
    nxmore = n - nx * ny
    !write(6, *)'nx, ny, nxmore', nx, ny, nxmore, im, jm
    allocate(npx(im), npxstart(nx), npxend(nx))
  ! --- Along x-direction.
    avepoint = sum(mask) / float(n)
    npxstart = 1; npxend = im
    do ii = 1, nx
      if(ii <= nxmore)then
        avenpx = avepoint * (ny + 1)
      else
        avenpx = avepoint * ny
      endif
  
      if(ii /= 1)npxstart(ii) = npxend(ii-1) + 1
      i1 = npxstart(ii)
  
      do i = i1+1, im
        n1 = sum(mask(i1:i-1, :))
        n2 = sum(mask(i1:i, :))
        if(n1 <= avenpx .and. n2 > avenpx)then
          npxend(ii) = i;  exit
        endif
      enddo
    enddo
  ! --- Along y-direction.
    nn = 0
    do ii = 1, nx
      i1 = npxstart(ii)
      i2 = npxend(ii)
      if(ii <= nxmore)then
        avenpy = sum(mask(i1:i2, :)) / (ny + 1); ny1 = ny + 1
      else
        avenpy = sum(mask(i1:i2, :)) / ny; ny1 = ny
      endif
      do jj = 1, ny1
        if(jj == 1)then
          j1 = 1
        else
          j1 = j2 + 1
        endif
        do j = j1+1, jm
          j2 = j
          n1 = sum(mask(i1:i2, j1:j2-1))
          n2 = sum(mask(i1:i2, j1:j2))
          if(avenpy >= n1 .and. avenpy < n2)exit
        enddo
        nn = nn + 1
        pebox(1:5, nn) = [nn-1, i1, i2, j1, j2]
      enddo
    enddo
  ! --- Ids of up & down.
    pebox(6, :) = MPI_PROC_NULL
    pebox(7, :) = MPI_PROC_NULL
    do nn = 1, n
      i1 = pebox(2, nn)
      i2 = pebox(3, nn)
      j1 = pebox(4, nn)
      j2 = pebox(5, nn)
      do i = 1, n
        if(i1 == pebox(2, i) .and. i2 == pebox(3, i) &
                             .and. j1 == pebox(5, i) + 1)then
          pebox(6, nn) = pebox(1, i)  ! my down
        endif
        if(i1 == pebox(2, i) .and. i2 == pebox(3, i) &
                             .and. j2 == pebox(4, i) - 1)then
          pebox(7, nn) = pebox(1, i)  ! my up
        endif
      enddo
    enddo
  ! --- halosize
    do nn = 1, n
      i1 = pebox(2, nn)
      i2 = pebox(3, nn)
      j1 = pebox(4, nn)
      j2 = pebox(5, nn)
      if(flagxcycle == 1)then
        pebox(11, nn) = max(i1 - halosize, 1)
        pebox(12, nn) = min(i2 + halosize, im)
      else
        pebox(11, nn) = i1 - halosize
        pebox(12, nn) = i2 + halosize
      endif
      pebox(13, nn) = max(1, j1 - halosize)
      pebox(14, nn) = min(jm, j2 + halosize)
      pebox(15, nn) = halosize
      pebox(16, nn) = pebox(12, nn) - pebox(11, nn) + 1 ! isize
      pebox(17, nn) = pebox(14, nn) - pebox(13, nn) + 1 ! jsize
    enddo
  ! --- Ids, js and je of right & left.
    peright = MPI_PROC_NULL; peleft = MPI_PROC_NULL; nl = 0; nr = 0
    do nn = 1, n
      i1 = pebox(2, nn)
      i2 = pebox(3, nn)
      j1 = pebox(4, nn)
      j2 = pebox(5, nn)
      nr = 0 ! ---- My right.
      do i = 1, n
        if(flagxcycle == 1)then
          if(i2 /= pebox(2, i) - 1)cycle
        else
          if(i2 /= pebox(2, i) - 1 .and. i2 /= im)cycle
          if(i2 == im .and. pebox(2, i) /= 1)cycle
        endif
        if(max(j2, pebox(5, i)) - min(j1, pebox(4, i)) &
            <= j2 + pebox(5, i) - j1 - pebox(4, i))then
          nr = nr + 1
          peright(nn, nr, 1) = pebox(1, i)
          peright(nn, nr, 2) = max(pebox(4, i), j1) !- j1 + 1
          peright(nn, nr, 3) = min(pebox(5, i), j2) !- j1 + 1
          peright(nn, nr, 4) = max(1 , peright(nn, nr, 2) - halosize) &
                             - pebox(13, nn) + 1
          peright(nn, nr, 5) = min(jm, peright(nn, nr, 3) + halosize) &
                             - pebox(13, nn) + 1
        endif
      enddo
      nl = 0 ! ---- My left.
      do i = 1, n
        if(flagxcycle == 1)then
          if(i1 /= pebox(3, i) + 1)cycle
        else
          if(i1 /= pebox(3, i) + 1 .and. i1 /= 1)cycle
          if(i1 == 1 .and. pebox(3, i) /= im)cycle
        endif
        if(max(j2, pebox(5, i)) - min(j1, pebox(4, i)) &
            <= j2 + pebox(5, i) - j1 - pebox(4, i))then
          nl = nl + 1
          peleft(nn, nl, 1) = pebox(1, i)
          peleft(nn, nl, 2) = max(pebox(4, i), j1) !- j1 + 1
          peleft(nn, nl, 3) = min(pebox(5, i), j2) !- j1 + 1
          peleft(nn, nl, 4) = max(1 , peleft(nn, nl, 2) - halosize) &
                            - pebox(13, nn) + 1
          peleft(nn, nl, 5) = min(jm, peleft(nn, nl, 3) + halosize) &
                            - pebox(13, nn) + 1
        endif
      enddo
      pebox(8, nn) = nr
      pebox(9, nn) = nl
      pebox(10, nn) = sum(mask(i1:i2, j1:j2))
    enddo
    do i = 1, n
      write(11, '(17(i5,x))')pebox(:, i)
      write(12, '(17(i5,x))')peleft(i, :, 1)
      write(13, '(17(i5,x))')peright(i, :, 1)
    enddo
  end subroutine set_pebox

!-------------------------------------------------------------------------------

  end module ympi_mod

!-------------------------------------------------------------------------------
!###############################################################################
