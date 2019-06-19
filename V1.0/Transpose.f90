module Transposeclass
  use Timerclass
  use Pgrid2dclass

contains

  subroutine trans2dnew_TRANSP(ny,nx,nysizex,nsizex,xin,tempmtr,np,stab,&
       rtab,comm)
    !use Timerclass
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ny,nx,nysizex,nsizex,np,comm
    double complex, dimension(ny,nsizex), intent(in) :: xin
    integer, dimension(0:np-1), intent(in) :: stab,rtab
    double complex, dimension(nx,nysizex), intent(out) :: tempmtr
    double complex, dimension(ny*nsizex) :: sendbuf
    double complex, dimension(nx*nysizex) :: recvbuf
    integer, dimension(0:np-1) :: sendcount,recvcount,senddisp, &
         recvdisp,rtabdisp,stabdisp 
    integer :: ierr,i,j,i0,j0,k
    double precision :: t0 !, t1
    integer :: ntmp

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    ntmp = nsizex

    if(np.ne.1) then

       sendcount = stab*ntmp
       senddisp(0) = 0
       do i = 1, np-1
          senddisp(i) = senddisp(i-1) + sendcount(i-1)
       enddo

       stabdisp = senddisp/ntmp

       do i = 0, np-1 
          i0 = senddisp(i)
          j0 = stabdisp(i)
          do j = 1, nsizex
             do k = 1, stab(i)
                sendbuf(k+(j-1)*stab(i)+i0) = xin(k+j0,j)
             enddo
          enddo
       enddo

       ntmp = nysizex
       recvcount = rtab*ntmp
       recvdisp(0) = 0
       do i = 1, np-1
          recvdisp(i) = recvdisp(i-1) + recvcount(i-1)
       enddo

       rtabdisp = recvdisp/ntmp

       ! do all-to-all scatter.
       call MPI_ALLTOALLV(sendbuf,sendcount,senddisp,MPI_DOUBLE_COMPLEX,&
            recvbuf,recvcount,recvdisp,MPI_DOUBLE_COMPLEX,comm,ierr)
       !        call MPI_BARRIER(comm,ierr)
       !        t_init = t_init + elapsedtime_Timer(t1)

       ! reshape the 1d receiving buffer into 2d array (nx*nysizex).
       do i = 0, np-1
          i0 = rtabdisp(i)
          j0 = recvdisp(i)
          do j = 1, nysizex
             do k = 1, rtab(i)
                tempmtr(i0+k,j) = recvbuf(j0+(k-1)*nysizex+j)
             enddo
          enddo
       enddo

    else

       do j = 1, nsizex
          do i = 1, ny
             tempmtr(j,i) = xin(i,j)
          enddo
       enddo

    endif

    !call MPI_BARRIER(mpicommwd,ierr)
    t_transp = t_transp + elapsedtime_Timer(t0)

  end subroutine trans2dnew_TRANSP



  ! Subroutine Trans2d: get the transpose of 2D double precision array with
  ! uneven distribution.
  subroutine trans2d2old_TRANSP(ny,nx,nysizex,nsizex,xin,xout,np,stab,&
       rtab,comm,myid)
    use Timerclass
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ny,nx,nysizex,nsizex,np,comm,myid
    double complex, dimension(ny,nsizex), intent(in) :: xin
    integer, dimension(0:np-1), intent(in) :: stab,rtab
    double complex, dimension(nx,nysizex), intent(out) :: xout
    double complex, allocatable, dimension(:,:) :: s
    double complex, allocatable, dimension(:,:) :: t
    integer, dimension(0:np-1) :: senddisp, &
         recvdisp 
    integer :: ierr,i,j,k,kstrt,ks,ir,msid,sendnum,recnum,joff,koff
    integer status(MPI_STATUS_SIZE)
    double precision :: t0

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    if(np.ne.1) then

       senddisp(0) = 0
       do i = 1, np-1
          senddisp(i) = senddisp(i-1) + stab(i-1)
       enddo

       recvdisp(0) = 0
       do i = 1, np-1
          recvdisp(i) = recvdisp(i-1) + rtab(i-1)
       enddo

       kstrt = myid + 1
       ks = myid -1
       do i = 1, np
          ir = i - (1+ks)
          if(ir.lt.1) ir = ir + np
          ! post receive
          allocate(t(nysizex,rtab(ir-1)))
          recnum = rtab(ir-1)*nysizex
          if(ir.ne.kstrt) call MPI_IRECV(t,recnum,MPI_DOUBLE_COMPLEX,&
               ir-1,i,comm,msid,ierr)
          ! send data
          allocate(s(stab(ir-1),nsizex))
          joff = senddisp(ir-1)
          do k = 1, nsizex
             do j = 1, stab(ir-1)
                s(j,k) = xin(j+joff,k)
             enddo
          enddo
          ! copy data to oneself directly.
          if(ir.eq.kstrt) then
             do k = 1, nsizex
                do j = 1, stab(ir-1)
                   t(j,k) = s(j,k)
                enddo
             enddo
          else
             sendnum = stab(ir-1)*nsizex
             call MPI_SEND(s,sendnum,MPI_DOUBLE_COMPLEX,ir-1,i,comm,ierr)
          endif
          !receive data
          if(ir.ne.kstrt) call MPI_WAIT(msid,status,ierr)
          koff = recvdisp(ir-1)
          do k = 1, rtab(ir-1)
             do j = 1, nysizex
                xout(k+koff,j) = t(j,k)
             enddo
          enddo
          deallocate(s)
          deallocate(t)
       enddo

    else

       do j = 1, nsizex
          do i = 1, ny
             xout(j,i) = xin(i,j)
          enddo
       enddo

    endif

    !call MPI_BARRIER(mpicommwd,ierr)
    t_transp = t_transp + elapsedtime_Timer(t0)

  end subroutine trans2d2old_TRANSP



  subroutine trans2d2_TRANSP(ny,nx,nysizex,nsizex,xin,xout,np,stab,&
       rtab,comm,myid)
    use Timerclass
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ny,nx,nysizex,nsizex,np,comm,myid
    double complex, dimension(ny,nsizex), intent(in) :: xin
    integer, dimension(0:np-1), intent(in) :: stab,rtab
    double complex, dimension(nx,nysizex), intent(out) :: xout
    double complex, allocatable, dimension(:,:) :: s
    double complex, allocatable, dimension(:,:) :: t
    integer, dimension(0:np-1) :: senddisp, &
         recvdisp 
    integer :: ierr,i,j,k,kstrt,ks,ir,msid,sendnum,recnum,joff,koff,nhalf
    integer status(MPI_STATUS_SIZE)
    double precision :: t0

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    if(np.ne.1) then

       senddisp(0) = 0
       do i = 1, np-1
          senddisp(i) = senddisp(i-1) + stab(i-1)
       enddo

       recvdisp(0) = 0
       do i = 1, np-1
          recvdisp(i) = recvdisp(i-1) + rtab(i-1)
       enddo

       if(myid.lt.np) then
          nhalf = 0
       else
          nhalf = np
       endif

       kstrt = myid + 1
       ks = myid -1
       do i = 1, np
          ir = i - (1+ks) + nhalf
          if(ir.lt.1)  then
             ir = ir + np + nhalf 
          else
             ir = ir + nhalf
          endif
          ! post receive
          allocate(t(nysizex,rtab(ir-1-nhalf)))
          recnum = rtab(ir-1-nhalf)*nysizex
          if(ir.ne.kstrt) call MPI_IRECV(t,recnum,MPI_DOUBLE_COMPLEX,&
               ir-1,i,comm,msid,ierr)
          ! send data
          allocate(s(stab(ir-1-nhalf),nsizex))
          joff = senddisp(ir-1-nhalf)
          do k = 1, nsizex
             do j = 1, stab(ir-1-nhalf)
                s(j,k) = xin(j+joff,k)
             enddo
          enddo
          ! copy data to oneself directly.
          if(ir.eq.kstrt) then
             do k = 1, nsizex
                do j = 1, stab(ir-1-nhalf)
                   t(j,k) = s(j,k)
                enddo
             enddo
          else
             sendnum = stab(ir-1-nhalf)*nsizex
             call MPI_SEND(s,sendnum,MPI_DOUBLE_COMPLEX,ir-1,i,comm,ierr)
          endif
          !receive data
          if(ir.ne.kstrt) call MPI_WAIT(msid,status,ierr)
          koff = recvdisp(ir-1-nhalf)
          do k = 1, rtab(ir-1-nhalf)
             do j = 1, nysizex
                xout(k+koff,j) = t(j,k)
             enddo
          enddo
          deallocate(s)
          deallocate(t)
       enddo

    else

       do j = 1, nsizex
          do i = 1, ny
             xout(j,i) = xin(i,j)
          enddo
       enddo

    endif

    !call MPI_BARRIER(mpicommwd,ierr)
    t_transp = t_transp + elapsedtime_Timer(t0)

  end subroutine trans2d2_TRANSP



  subroutine trans2d2real_TRANSP(ny,nx,nysizex,nsizex,xin,xout,np,stab,&
       rtab,comm,myid)
    use Timerclass
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ny,nx,nysizex,nsizex,np,comm,myid
    double precision, dimension(ny,nsizex), intent(in) :: xin
    integer, dimension(0:np-1), intent(in) :: stab,rtab
    double precision, dimension(nx,nysizex), intent(out) :: xout
    double precision, allocatable, dimension(:,:) :: s
    double precision, allocatable, dimension(:,:) :: t
    integer, dimension(0:np-1) :: senddisp, &
         recvdisp 
    integer :: ierr,i,j,k,kstrt,ks,ir,msid,sendnum,recnum,joff,koff,nhalf
    integer status(MPI_STATUS_SIZE)
    double precision :: t0

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    if(np.ne.1) then

       senddisp(0) = 0
       do i = 1, np-1
          senddisp(i) = senddisp(i-1) + stab(i-1)
       enddo

       recvdisp(0) = 0
       do i = 1, np-1
          recvdisp(i) = recvdisp(i-1) + rtab(i-1)
       enddo

       if(myid.lt.np) then
          nhalf = 0
       else
          nhalf = np
       endif

       kstrt = myid + 1
       ks = myid -1
       do i = 1, np
          ir = i - (1+ks) + nhalf
          if(ir.lt.1)  then
             ir = ir + np + nhalf 
          else
             ir = ir + nhalf
          endif
          ! post receive
          allocate(t(nysizex,rtab(ir-1-nhalf)))
          recnum = rtab(ir-1-nhalf)*nysizex
          if(ir.ne.kstrt) call MPI_IRECV(t,recnum,MPI_DOUBLE_PRECISION,&
               ir-1,i,comm,msid,ierr)
          ! send data
          allocate(s(stab(ir-1-nhalf),nsizex))
          joff = senddisp(ir-1-nhalf)
          do k = 1, nsizex
             do j = 1, stab(ir-1-nhalf)
                s(j,k) = xin(j+joff,k)
             enddo
          enddo
          ! copy data to oneself directly.
          if(ir.eq.kstrt) then
             do k = 1, nsizex
                do j = 1, stab(ir-1-nhalf)
                   t(j,k) = s(j,k)
                enddo
             enddo
          else
             sendnum = stab(ir-1-nhalf)*nsizex
             call MPI_SEND(s,sendnum,MPI_DOUBLE_PRECISION,ir-1,i,comm,ierr)
          endif
          !receive data
          if(ir.ne.kstrt) call MPI_WAIT(msid,status,ierr)
          koff = recvdisp(ir-1-nhalf)
          do k = 1, rtab(ir-1-nhalf)
             do j = 1, nysizex
                xout(k+koff,j) = t(j,k)
             enddo
          enddo
          deallocate(s)
          deallocate(t)
       enddo

    else

       do j = 1, nsizex
          do i = 1, ny
             xout(j,i) = xin(i,j)
          enddo
       enddo

    endif

    !call MPI_BARRIER(mpicommwd,ierr)
    t_transp = t_transp + elapsedtime_Timer(t0)

  end subroutine trans2d2real_TRANSP



  ! Subroutine Trans2d: get the transpose of 2D double precision array with
  ! uneven distribution using 1 group processors..
  subroutine trans2d21G_TRANSP(ny,nx,nysizex,nsizex,xin,xout,np,stab,&
       rtab,comm,myid)
    use Timerclass
    implicit none
    include 'mpif.h'
    integer, intent(in) :: ny,nx,nysizex,nsizex,np,comm,myid
    double complex, dimension(ny,nsizex), intent(in) :: xin
    integer, dimension(0:np-1), intent(in) :: stab,rtab
    double complex, dimension(nx,nysizex), intent(out) :: xout
    double complex, allocatable, dimension(:,:) :: s
    double complex, allocatable, dimension(:,:) :: t
    integer, dimension(0:np-1) :: senddisp, &
         recvdisp 
    integer :: ierr,i,j,k,kstrt,ks,ir,msid,sendnum,recnum,joff,koff
    integer status(MPI_STATUS_SIZE)
    double precision :: t0

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    if(np.ne.1) then

       senddisp(0) = 0
       do i = 1, np-1
          senddisp(i) = senddisp(i-1) + stab(i-1)
       enddo

       recvdisp(0) = 0
       do i = 1, np-1
          recvdisp(i) = recvdisp(i-1) + rtab(i-1)
       enddo

       kstrt = myid + 1
       ks = myid -1
       do i = 1, np
          ir = i - (1+ks)
          if(ir.lt.1) ir = ir + np
          ! post receive
          allocate(t(nysizex,rtab(ir-1)))
          recnum = rtab(ir-1)*nysizex
          if(ir.ne.kstrt) call MPI_IRECV(t,recnum,MPI_DOUBLE_COMPLEX,&
               ir-1,i,comm,msid,ierr)
          ! send data
          allocate(s(stab(ir-1),nsizex))
          joff = senddisp(ir-1)
          do k = 1, nsizex
             do j = 1, stab(ir-1)
                s(j,k) = xin(j+joff,k)
             enddo
          enddo
          ! copy data to oneself directly.
          if(ir.eq.kstrt) then
             do k = 1, nsizex
                do j = 1, stab(ir-1)
                   t(j,k) = s(j,k)
                enddo
             enddo
          else
             sendnum = stab(ir-1)*nsizex
             call MPI_SEND(s,sendnum,MPI_DOUBLE_COMPLEX,ir-1,i,comm,ierr)
          endif
          !receive data
          if(ir.ne.kstrt) call MPI_WAIT(msid,status,ierr)
          koff = recvdisp(ir-1)
          do k = 1, rtab(ir-1)
             do j = 1, nysizex
                xout(k+koff,j) = t(j,k)
             enddo
          enddo
          deallocate(s)
          deallocate(t)
       enddo

    else

       do j = 1, nsizex
          do i = 1, ny
             xout(j,i) = xin(i,j)
          enddo
       enddo

    endif

    !call MPI_BARRIER(mpicommwd,ierr)
    t_transp = t_transp + elapsedtime_Timer(t0)

  end subroutine trans2d21G_TRANSP


end module Transposeclass
