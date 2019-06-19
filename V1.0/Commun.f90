module Communclass
  use Pgrid2dclass
  use Timerclass
  include 'mpif.h'
contains
  ! In all the following subroutine 1-4, we have assumed that the
  ! source is 0 or np and total processor is 2np
  subroutine reduce(tmplc,tmpgl,ndata,datatype,op,source,np,myid,comm)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: ndata,datatype,op,comm,source,myid,np
    double precision, dimension(ndata), intent(in) :: tmplc
    double precision, dimension(ndata), intent(out) :: tmpgl
    double precision, dimension(ndata) :: recvbuf1
    integer :: i,j,ierr
    integer status(MPI_STATUS_SIZE)

    if(myid.eq.source) then
       tmpgl = tmplc
       do i = 1,np-1
          call MPI_RECV(recvbuf1,ndata,MPI_DOUBLE_PRECISION,i,&
               i,comm,status,ierr)
          if(op.eq.1) then
             tmpgl = tmpgl + recvbuf1
          else if(op.eq.2) then
             do j = 1, ndata
                if(recvbuf1(j).gt.tmpgl(j)) then
                   tmpgl(j) = recvbuf1(j)
                endif
             enddo
          else if(op.eq.3) then
             do j = 1, ndata
                if(recvbuf1(j).lt.tmpgl(j)) then
                   tmpgl(j) = recvbuf1(j)
                endif
             enddo
          endif
       enddo
    else
       if(myid.lt.np) then
          call MPI_SEND(tmplc,ndata,MPI_DOUBLE_PRECISION,source,&
               myid,comm,ierr)
       endif
    endif

  end subroutine reduce



  subroutine reduce2(tmplc,tmpgl,ndata,datatype,op,source,np,myid,comm)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: ndata,datatype,op,comm,source,np,myid
    double precision, dimension(ndata), intent(in) :: tmplc
    double precision, dimension(ndata), intent(out) :: tmpgl
    double precision, dimension(ndata) :: recvbuf1
    integer :: i,j,ierr
    integer status(MPI_STATUS_SIZE)

    if(myid.eq.source) then

       if(op.eq.1) then
          tmpgl = 0.0
       else if(op.eq.2) then
          tmpgl = -1.0e20
       else if(op.eq.3) then
          tmpgl = 1.0e20
       endif

       do i = np,2*np - 1
          call MPI_RECV(recvbuf1,ndata,MPI_DOUBLE_PRECISION,i,&
               i,comm,status,ierr)
          if(op.eq.1) then
             tmpgl = tmpgl + recvbuf1
          else if(op.eq.2) then
             do j = 1, ndata
                if(recvbuf1(j).gt.tmpgl(j)) then
                   tmpgl(j) = recvbuf1(j)
                endif
             enddo
          else if(op.eq.3) then
             do j = 1, ndata
                if(recvbuf1(j).lt.tmpgl(j)) then
                   tmpgl(j) = recvbuf1(j)
                endif
             enddo
          endif
       enddo
    else
       if(myid.ge.np) then
          call MPI_SEND(tmplc,ndata,MPI_DOUBLE_PRECISION,source,&
               myid,comm,ierr)
       endif
    endif

  end subroutine reduce2



  subroutine reduce3(tmplc,tmpgl,ndata,datatype,op,source,np,myid,comm)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: ndata,datatype,op,comm,source,np,myid
    double precision, dimension(ndata), intent(in) :: tmplc
    double precision, dimension(ndata), intent(out) :: tmpgl
    double precision, dimension(ndata) :: recvbuf1
    integer :: i,j,ierr
    integer status(MPI_STATUS_SIZE)

    if(myid.eq.source) then
       tmpgl = tmplc
       do i = np+1,np+np-1
          call MPI_RECV(recvbuf1,ndata,MPI_DOUBLE_PRECISION,i,&
               i,comm,status,ierr)
          if(op.eq.1) then
             tmpgl = tmpgl + recvbuf1
          else if(op.eq.2) then
             do j = 1, ndata
                if(recvbuf1(j).gt.tmpgl(j)) then
                   tmpgl(j) = recvbuf1(j)
                endif
             enddo
          else if(op.eq.3) then
             do j = 1, ndata
                if(recvbuf1(j).lt.tmpgl(j)) then
                   tmpgl(j) = recvbuf1(j)
                endif
             enddo
          endif
       enddo
    else
       if(myid.ge.np) then
          call MPI_SEND(tmplc,ndata,MPI_DOUBLE_PRECISION,source,&
               myid,comm,ierr)
       endif
    endif

  end subroutine reduce3



  subroutine reduce4(tmplc,tmpgl,ndata,datatype,op,source1,source2,np,&
       myid,comm)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: ndata,datatype,op,comm,source1,source2,np,myid
    double precision, dimension(ndata), intent(in) :: tmplc
    double precision, dimension(ndata), intent(out) :: tmpgl
    double precision, dimension(ndata) :: recvbuf1
    integer :: i,j,ierr
    integer status(MPI_STATUS_SIZE)

    if(myid.eq.source1) then
       tmpgl = tmplc
       do i = 1,np-1
          call MPI_RECV(recvbuf1,ndata,MPI_DOUBLE_PRECISION,i,&
               i,comm,status,ierr)
          if(op.eq.1) then
             tmpgl = tmpgl + recvbuf1
          else if(op.eq.2) then
             do j = 1, ndata
                if(recvbuf1(j).gt.tmpgl(j)) then
                   tmpgl(j) = recvbuf1(j)
                endif
             enddo
          else if(op.eq.3) then
             do j = 1, ndata
                if(recvbuf1(j).lt.tmpgl(j)) then
                   tmpgl(j) = recvbuf1(j)
                endif
             enddo
          endif
       enddo
    else if(myid.eq.source2) then
       tmpgl = tmplc
       do i = np+1,np+np-1
          call MPI_RECV(recvbuf1,ndata,MPI_DOUBLE_PRECISION,i,&
               i,comm,status,ierr)
          if(op.eq.1) then
             tmpgl = tmpgl + recvbuf1
          else if(op.eq.2) then
             do j = 1, ndata
                if(recvbuf1(j).gt.tmpgl(j)) then
                   tmpgl(j) = recvbuf1(j)
                endif
             enddo
          else if(op.eq.3) then
             do j = 1, ndata
                if(recvbuf1(j).lt.tmpgl(j)) then
                   tmpgl(j) = recvbuf1(j)
                endif
             enddo
          endif
       enddo
    else if(myid.lt.np) then
       call MPI_SEND(tmplc,ndata,MPI_DOUBLE_PRECISION,source1,&
            myid,comm,ierr)
    else if(myid.ge.np) then
       call MPI_SEND(tmplc,ndata,MPI_DOUBLE_PRECISION,source2,&
            myid,comm,ierr)
    endif

  end subroutine reduce4



  ! sum up the contributions of rho from all other processors
  subroutine allreduce(rhotmp,rholc,ndata,npy,myidy,commcol,op) 
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: ndata,npy,commcol,myidy,op
    double precision, dimension(ndata), intent(in) :: rhotmp
    double precision, dimension(ndata), intent(out) :: rholc
    double precision, dimension(ndata) :: sendbuf1
    double precision, dimension(ndata) :: recvbuf1
    integer :: i,ierr,msid,k,mydes,mysou,nsend
    integer status(MPI_STATUS_SIZE)

    call MPI_BARRIER(commcol,ierr)

    !Here, npy is half of total Npy since half PE is assigned to one beam.
    if(myidy.lt.npy) then

       do k = 1, ndata
          rholc(k) = rhotmp(k)
       enddo

       do i = 1, npy-1
          !find the target processor
          if((myidy-i).lt.0) then
             mydes = myidy - i + npy
          else
             mydes = myidy - i
          endif
          !find the source processor
          if((myidy+i).ge.npy) then
             mysou = myidy + i - npy
          else
             mysou = myidy + i
          endif

          do k = 1, ndata
             !Nytot is uniformly distributed along npy.
             sendbuf1(k) = rhotmp(k)
          enddo

          nsend = ndata

          call MPI_IRECV(recvbuf1,nsend,MPI_DOUBLE_PRECISION,mysou,&
               i,commcol,msid,ierr)
          call MPI_SEND(sendbuf1,nsend,MPI_DOUBLE_PRECISION,mydes,&
               i,commcol,ierr)
          call MPI_WAIT(msid,status,ierr)

          if(op.eq.1) then ! summation
             do k = 1, ndata
                rholc(k) = rholc(k) + recvbuf1(k)
             enddo
          else if(op.eq.2) then ! maximum
             do k = 1, ndata
                if(recvbuf1(k).gt.rholc(k)) then
                   rholc(k) = recvbuf1(k)
                endif
             enddo
          else if(op.eq.3) then ! minimum
             do k = 1, ndata
                if(recvbuf1(k).lt.rholc(k)) then
                   rholc(k) = recvbuf1(k)
                endif
             enddo
          else
             print*,"wrong operation!"
             stop
          endif
       enddo

    else
       do k = 1, ndata
          rholc(k) = rhotmp(k)
       enddo
       do i = 1, npy-1
          !find the target processor
          if((myidy-i-npy).lt.0) then
             mydes = myidy - i + npy
          else
             mydes = myidy - i
          endif
          !find the source processor
          if((myidy+i-npy).ge.npy) then
             mysou = myidy + i - npy
          else
             mysou = myidy + i
          endif

          do k = 1, ndata
             sendbuf1(k) = rhotmp(k)
          enddo

          nsend = ndata

          call MPI_IRECV(recvbuf1,nsend,MPI_DOUBLE_PRECISION,mysou,&
               i,commcol,msid,ierr)
          call MPI_SEND(sendbuf1,nsend,MPI_DOUBLE_PRECISION,mydes,&
               i,commcol,ierr)
          call MPI_WAIT(msid,status,ierr)

          if(op.eq.1) then ! summation
             do k = 1, ndata
                rholc(k) = rholc(k) + recvbuf1(k)
             enddo
          else if(op.eq.2) then ! maximum
             do k = 1, ndata
                if(recvbuf1(k).gt.rholc(k)) then
                   rholc(k) = recvbuf1(k)
                endif
             enddo
          else if(op.eq.3) then ! minimum
             do k = 1, ndata
                if(recvbuf1(k).lt.rholc(k)) then
                   rholc(k) = recvbuf1(k)
                endif
             enddo
          else
             print*,"wrong operation!"
             stop
          endif

       enddo
    endif

    call MPI_BARRIER(commcol,ierr)

  end subroutine allreduce



  !particle manager for 1d decomposition
  subroutine ptsmv1d_ptclmger(Ptsl,Nptlocal,pdim,lcrange,myidy,npy,&
       commcol,maxpt)
    implicit none
    !include 'mpif.h'
    integer, intent(inout) :: Nptlocal
    integer, intent(in) :: pdim,npy,myidy,commcol,maxpt
    !double precision, pointer, dimension(:,:) :: Ptsl
    double precision, dimension(6,maxpt) :: Ptsl
    double precision, dimension(:),intent(in) :: lcrange
    integer, parameter :: nptmv = 100000
    double precision, dimension(6,3*nptmv) :: up,down,recvup,recvdown
    double precision, dimension(6,6*nptmv) :: recv
    double precision, dimension(6,Nptlocal) :: temp1
    integer :: iup,idown
    integer :: jup,jdown
    integer :: myup,mydown
    integer :: i,j,numpts
    integer :: nout,ii,nmv,totnmv
    integer :: msid,ierr,ikeep
    integer status(MPI_STATUS_SIZE) 
    integer statarry(MPI_STATUS_SIZE,4), req(4)
    double precision :: t0

    call starttime_Timer(t0)

    if(myidy.ne.npy-1) then
       myup = myidy + 1
    else
       myup = MPI_PROC_NULL
    endif
    if(myidy.ne.0) then
       mydown = myidy -1
    else
       mydown = MPI_PROC_NULL
    endif

    ikeep = 0
    iup = 0
    idown = 0
    do i = 1, Nptlocal
       if(Ptsl(3,i).gt.lcrange(4)) then
          if(myidy.ne.(npy-1)) then
             iup = iup + 1
             up(:,iup) = Ptsl(:,i)
          else
             ikeep = ikeep + 1
             temp1(:,ikeep) = Ptsl(:,i)
          endif
       else if(Ptsl(3,i).le.lcrange(3)) then
          if(myidy.ne.0) then
             idown = idown + 1
             down(:,idown) = Ptsl(:,i)
          else
             ikeep = ikeep + 1
             temp1(:,ikeep) = Ptsl(:,i)
          endif
       else
          ikeep = ikeep + 1
          temp1(:,ikeep) = Ptsl(:,i)
       endif
    enddo

    ikeep = 0
    nout = iup + idown
    do
       jup = 0
       jdown = 0

       call MPI_IRECV(jdown,1,MPI_INTEGER,myup,0,commcol,req(1),&
            ierr)
       call MPI_IRECV(jup,1,MPI_INTEGER,mydown,0,commcol,req(2),&
            ierr)
       call MPI_ISEND(idown,1,MPI_INTEGER,mydown,0,commcol,req(3),&
            ierr)
       call MPI_ISEND(iup,1,MPI_INTEGER,myup,0,commcol,req(4),&
            ierr)
       call MPI_WAITALL(4,req,statarry,ierr) 

       !if(jdown+jup.gt.maxpt) then
       !  stop
       !endif
       !if(jup.gt.maxpt) then
       !  stop
       !endif

       !send outgoing particles to down neibhoring processor.
       jdown = 6*jdown
       idown = 6*idown
       call MPI_IRECV(recvup(1,1),jdown,MPI_DOUBLE_PRECISION,myup,&
            0,commcol,msid,ierr)
       call MPI_SEND(down(1,1),idown,MPI_DOUBLE_PRECISION,mydown,&
            0,commcol,ierr)
       call MPI_WAIT(msid,status,ierr)
       idown = idown/6
       jdown = jdown/6

       !send outgoing particles to up neibhoring processor.
       jup = 6*jup
       iup = 6*iup
       call MPI_IRECV(recvdown(1,1),jup,MPI_DOUBLE_PRECISION,mydown,&
            0,commcol,msid,ierr)
       call MPI_SEND(up(1,1),iup,MPI_DOUBLE_PRECISION,myup,&
            0,commcol,ierr)
       call MPI_WAIT(msid,status,ierr)
       iup = iup/6
       jup = jup/6

       iup = 0
       do i = 1, jup
          if(recvdown(3,i).gt.lcrange(4)) then
             if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recvdown(:,i)
             else
                ikeep = ikeep + 1
                !if(ikeep.gt.maxpt) then
                !  stop
                !endif
                recv(:,ikeep) = recvdown(:,i)
             endif
          else 
             ikeep = ikeep + 1
             !if(ikeep.gt.maxpt) then
             !    stop
             !endif
             recv(:,ikeep) = recvdown(:,i)
          endif
       enddo

       idown = 0
       do i = 1, jdown
          if(recvup(3,i).lt.lcrange(3)) then
             if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recvup(:,i)
             else
                ikeep = ikeep + 1
                !if(ikeep.gt.maxpt) then
                !  stop
                !endif
                recv(:,ikeep) = recvup(:,i)
             endif
          else 
             ikeep = ikeep + 1
             !if(ikeep.gt.maxpt) then
             !    stop
             !endif
             recv(:,ikeep) = recvup(:,i)
          endif
       enddo

       nmv = idown+iup
       call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
            commcol,ierr)
       if(totnmv.eq.0) then
          exit
       endif
    enddo

    numpts = Nptlocal-nout

    !recopy the remaining local particles back to Ptsl which has 
    !a new size now.
    Nptlocal = numpts+ikeep
    if(Nptlocal.gt.maxpt) then
       print*,"overflow in ptsmv, Nptlocal: ",Nptlocal
       stop
    endif
    !deallocate(Ptsl)
    !allocate(Ptsl(6,Nptlocal))
    do i = 1, numpts
       do j = 1, 6
          Ptsl(j,i) = temp1(j,i)
       enddo
    enddo

    do i = 1, ikeep
       ii = i + numpts
       do j = 1, 6
          Ptsl(j,ii) = recv(j,i)
       enddo
    enddo

    t_ptsmv = t_ptsmv + elapsedtime_Timer(t0)

  end subroutine ptsmv1d_ptclmger



  ! move particles from one processor to 4 neighboring processors.
  subroutine ptsmv1dold_ptclmger(Ptsl,Nptlocal,pdim,lcrange,myidy,npy,&
       commcol)
    implicit none
    !include 'mpif.h'
    integer, intent(inout) :: Nptlocal
    integer, intent(in) :: pdim,npy,myidy,commcol
    double precision, pointer, dimension(:,:) :: Ptsl
    double precision, dimension(:),intent(in) :: lcrange
    integer, parameter :: nptmv = 10000
    double precision, dimension(6,3*nptmv) :: up,down,recvup,recvdown
    double precision, dimension(6,6*nptmv) :: recv
    double precision, dimension(6,Nptlocal) :: temp1
    integer :: iup,idown
    integer :: jup,jdown
    integer :: myup,mydown
    integer :: i,j,numpts
    integer :: nout,ii,nmv,totnmv
    integer :: msid,ierr,ikeep
    integer status(MPI_STATUS_SIZE) 
    integer statarry(MPI_STATUS_SIZE,4), req(4)
    double precision :: t0

    call starttime_Timer(t0)

    if(myidy.ne.npy-1) then
       myup = myidy + 1
    else
       myup = MPI_PROC_NULL
    endif
    if(myidy.ne.0) then
       mydown = myidy -1
    else
       mydown = MPI_PROC_NULL
    endif

    ikeep = 0
    iup = 0
    idown = 0
    do i = 1, Nptlocal
       if(Ptsl(3,i).gt.lcrange(4)) then
          if(myidy.ne.(npy-1)) then
             iup = iup + 1
             up(:,iup) = Ptsl(:,i)
          else
             ikeep = ikeep + 1
             temp1(:,ikeep) = Ptsl(:,i)
          endif
       else if(Ptsl(3,i).le.lcrange(3)) then
          if(myidy.ne.0) then
             idown = idown + 1
             down(:,idown) = Ptsl(:,i)
          else
             ikeep = ikeep + 1
             temp1(:,ikeep) = Ptsl(:,i)
          endif
       else
          ikeep = ikeep + 1
          temp1(:,ikeep) = Ptsl(:,i)
       endif
    enddo

    ikeep = 0
    nout = iup + idown
    do
       jup = 0
       jdown = 0

       call MPI_IRECV(jdown,1,MPI_INTEGER,myup,0,commcol,req(1),&
            ierr)
       call MPI_IRECV(jup,1,MPI_INTEGER,mydown,0,commcol,req(2),&
            ierr)
       call MPI_ISEND(idown,1,MPI_INTEGER,mydown,0,commcol,req(3),&
            ierr)
       call MPI_ISEND(iup,1,MPI_INTEGER,myup,0,commcol,req(4),&
            ierr)
       call MPI_WAITALL(4,req,statarry,ierr) 

       !send outgoing particles to down neibhoring processor.
       jdown = 6*jdown
       idown = 6*idown
       call MPI_IRECV(recvup(1,1),jdown,MPI_DOUBLE_PRECISION,myup,&
            0,commcol,msid,ierr)
       call MPI_SEND(down(1,1),idown,MPI_DOUBLE_PRECISION,mydown,&
            0,commcol,ierr)
       call MPI_WAIT(msid,status,ierr)
       idown = idown/6
       jdown = jdown/6

       !send outgoing particles to up neibhoring processor.
       jup = 6*jup
       iup = 6*iup
       call MPI_IRECV(recvdown(1,1),jup,MPI_DOUBLE_PRECISION,mydown,&
            0,commcol,msid,ierr)
       call MPI_SEND(up(1,1),iup,MPI_DOUBLE_PRECISION,myup,&
            0,commcol,ierr)
       call MPI_WAIT(msid,status,ierr)
       iup = iup/6
       jup = jup/6

       iup = 0
       do i = 1, jup
          if(recvdown(3,i).gt.lcrange(4)) then
             if(myidy.ne.(npy-1)) then
                iup = iup + 1
                up(:,iup) = recvdown(:,i)
             else
                ikeep = ikeep + 1
                recv(:,ikeep) = recvdown(:,i)
             endif
          else 
             ikeep = ikeep + 1
             recv(:,ikeep) = recvdown(:,i)
          endif
       enddo

       idown = 0
       do i = 1, jdown
          if(recvup(3,i).lt.lcrange(3)) then
             if(myidy.ne.0) then
                idown = idown + 1
                down(:,idown) = recvup(:,i)
             else
                ikeep = ikeep + 1
                recv(:,ikeep) = recvup(:,i)
             endif
          else 
             ikeep = ikeep + 1
             recv(:,ikeep) = recvup(:,i)
          endif
       enddo

       nmv = idown+iup
       call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
            commcol,ierr)
       if(totnmv.eq.0) then
          exit
       endif
    enddo

    numpts = Nptlocal-nout

    !recopy the remaining local particles back to Ptsl which has 
    !a new size now.
    Nptlocal = numpts+ikeep
    deallocate(Ptsl)
    allocate(Ptsl(6,Nptlocal))
    do i = 1, numpts
       do j = 1, 6
          Ptsl(j,i) = temp1(j,i)
       enddo
    enddo

    do i = 1, ikeep
       ii = i + numpts
       do j = 1, 6
          Ptsl(j,ii) = recv(j,i)
       enddo
    enddo

    t_ptsmv = t_ptsmv + elapsedtime_Timer(t0)

  end subroutine ptsmv1dold_ptclmger



  ! move particles from one processor to 4 neighboring processors.
  subroutine ptsmv1dold2_ptclmger(Ptsl,Nptlocal,grid,pdim,lcrange)
    implicit none
    !include 'mpif.h'
    type (Pgrid2d), intent(in) :: grid
    integer, intent(inout) :: Nptlocal
    integer, intent(in) :: pdim
    double precision, pointer, dimension(:,:) :: Ptsl
    double precision, dimension(:),intent(in) :: lcrange
    integer, parameter :: nptmv = 10000
    double precision, dimension(6,3*nptmv) :: up,down
    double precision, allocatable, dimension(:,:) :: temp1,recv
    integer :: myid,myidx,myidy,totnp,npy,npx, &
         comm2d,commcol,commrow
    integer :: ileft,iright,iup,idown
    integer :: jleft,jright,jup,jdown
    integer :: myleft,myright,myup,mydown
    integer :: i,j,numpts,ic
    logical, dimension(Nptlocal) :: msk
    logical, allocatable, dimension(:) :: mmsk
    integer :: numbuf,nmv,nmv0,nout,ii,totnmv,ij
    integer :: msid,ierr
    integer status(MPI_STATUS_SIZE) 
    integer statarry(MPI_STATUS_SIZE,4), req(4)
    integer :: flag,Nptlocal0,nst,nout0,iileft,iiright,iiup,iidown
    integer :: totflag
    double precision :: t0

    call starttime_Timer(t0)

    call getsize_Pgrid2d(grid,totnp,npy,npx)
    call getpost_Pgrid2d(grid,myid,myidy,myidx)
    call getcomm_Pgrid2d(grid,comm2d,commcol,commrow)
    if(myidx.ne.(npx-1)) then
       myright = myidx + 1
    else
       myright = MPI_PROC_NULL
    endif
    if(myidx.ne.0) then
       myleft = myidx - 1
    else
       myleft = MPI_PROC_NULL
    endif

    if(myidy.ne.npy-1) then
       myup = myidy + 1
    else
       myup = MPI_PROC_NULL
    endif
    if(myidy.ne.0) then
       mydown = myidy -1
    else
       mydown = MPI_PROC_NULL
    endif

    !        call MPI_BARRIER(comm2d,ierr)

    flag = 0
    Nptlocal0 = Nptlocal
    nout0 = 0

    do 
       ileft = 0
       iright = 0
       iup = 0
       idown = 0
       iileft = 0
       iiright = 0
       iiup = 0
       iidown = 0
       do i = 1, Nptlocal0 - nout0
          msk(i) = .true.
          if(Ptsl(3,i).gt.lcrange(4)) then
             if(myidy.ne.(npy-1)) then
                iiup = iiup + 1
                if(iiup.le.nptmv) then
                   iup = iup + 1
                   up(:,iup) = Ptsl(:,i)
                   msk(i) = .false.
                endif
             endif
          else if(Ptsl(3,i).le.lcrange(3)) then
             if(myidy.ne.0) then
                iidown = iidown + 1
                if(iidown.le.nptmv) then
                   idown = idown + 1
                   down(:,idown) = Ptsl(:,i)
                   msk(i) = .false.
                endif
             endif
          else
          endif
       enddo

       if((iileft.gt.nptmv).or.(iiright.gt.nptmv).or.(iiup.gt.nptmv) &
            .or.(iidown.gt.nptmv)) then
          flag = 1
       else
          flag = 0
       endif

       nmv0 = 0
       nout = ileft+iright+iup+idown
       allocate(recv(6,nmv0))
       allocate(temp1(6,nmv0))
       ij = 0
       call MPI_BARRIER(comm2d,ierr)
       jleft = 0
       jright = 0

       do
          ij = ij + 1

          !        call MPI_BARRIER(commrow,ierr)
          !        if(myid.eq.0) then
          !        endif

          jup = 0
          jdown = 0

          call MPI_IRECV(jdown,1,MPI_INTEGER,myup,0,commcol,req(1),&
               ierr)
          call MPI_IRECV(jup,1,MPI_INTEGER,mydown,0,commcol,req(2),&
               ierr)
          call MPI_ISEND(idown,1,MPI_INTEGER,mydown,0,commcol,req(3),&
               ierr)
          call MPI_ISEND(iup,1,MPI_INTEGER,myup,0,commcol,req(4),&
               ierr)
          call MPI_WAITALL(4,req,statarry,ierr) 

          numbuf = jup+jdown 

          !        call MPI_BARRIER(commcol,ierr)
          !        if(myid.eq.0) then
          !        endif

          deallocate(recv)
          allocate(recv(6,numbuf+nmv0))
          do i = 1, nmv0
             recv(:,i) = temp1(:,i)
          enddo
          deallocate(temp1)

          !        call MPI_BARRIER(commrow,ierr)
          !        if(myid.eq.0) then
          !        endif

          nst = nmv0 + 1
          !send outgoing particles to down neibhoring processor.
          jdown = 6*jdown
          idown = 6*idown
          call MPI_IRECV(recv(1,nst),jdown,MPI_DOUBLE_PRECISION,myup,&
               0,commcol,msid,ierr)
          call MPI_SEND(down(1,1),idown,MPI_DOUBLE_PRECISION,mydown,&
               0,commcol,ierr)
          call MPI_WAIT(msid,status,ierr)
          idown = idown/6
          jdown = jdown/6
          nmv0 = nmv0 + jdown

          nst = nmv0 + 1
          !send outgoing particles to up neibhoring processor.
          jup = 6*jup
          iup = 6*iup
          call MPI_IRECV(recv(1,nst),jup,MPI_DOUBLE_PRECISION,mydown,&
               0,commcol,msid,ierr)
          call MPI_SEND(up(1,1),iup,MPI_DOUBLE_PRECISION,myup,&
               0,commcol,ierr)
          call MPI_WAIT(msid,status,ierr)
          iup = iup/6
          jup = jup/6
          nmv0 = nmv0 + jup

          !        call MPI_BARRIER(commcol,ierr)
          !        if(myid.eq.0) then
          !        endif

          allocate(mmsk(numbuf))
          ileft = 0
          iright = 0
          iup = 0
          idown = 0
          nmv0 = nmv0 - numbuf
          do i = 1, numbuf
             mmsk(i) = .true.
             ii = i+nmv0
             if(recv(3,ii).gt.lcrange(4)) then
                if(myidy.ne.(npy-1)) then
                   iup = iup + 1
                   up(:,iup) = recv(:,ii)
                   mmsk(i) = .false.
                endif
             else if(recv(3,ii).le.lcrange(3)) then
                if(myidy.ne.0) then
                   idown = idown + 1
                   down(:,idown) = recv(:,ii)
                   mmsk(i) = .false.
                endif
             else
             endif
          enddo
          nmv = ileft+iright+idown+iup
          call MPI_ALLREDUCE(nmv,totnmv,1,MPI_INTEGER,MPI_SUM, &
               comm2d,ierr)
          if(totnmv.eq.0) then
             nmv0 = nmv0 + numbuf - nmv
             deallocate(mmsk)
             exit
          endif

          ic = 0
          allocate(temp1(6,nmv0+numbuf-nmv))
          do i = 1, nmv0
             temp1(:,i) = recv(:,i)
          enddo
          do i = 1, numbuf
             ii = i + nmv0
             if(mmsk(i)) then
                ic = ic + 1
                temp1(:,ic+nmv0) = recv(:,ii)
             endif
          enddo
          nmv0 = nmv0 + numbuf - nmv
          deallocate(mmsk)

          !        call MPI_BARRIER(comm2d,ierr)

       enddo

       !copy the remaining local particles into a temporary array.
       numpts = Nptlocal-nout
       allocate(temp1(6,numpts))
       ic = 0
       do i = 1, Nptlocal0-nout0
          if(msk(i)) then
             ic = ic + 1
             do j = 1, 6
                temp1(j,ic) = Ptsl(j,i)
             enddo
          endif
       enddo
       do i = Nptlocal0-nout0+1, Nptlocal
          ii = i-nout
          do j = 1, 6
             temp1(j,ii) = Ptsl(j,i)
          enddo
       enddo

       !        call MPI_BARRIER(comm2d,ierr)
       !recopy the remaining local particles back to Ptsl which has 
       !a new size now.
       Nptlocal = numpts+nmv0 
       deallocate(Ptsl)
       allocate(Ptsl(6,Nptlocal))
       do i = 1, numpts
          do j = 1, 6
             Ptsl(j,i) = temp1(j,i)
          enddo
       enddo
       deallocate(temp1)
       do i = 1, nmv0
          ii = i + numpts
          do j = 1, 6
             Ptsl(j,ii) = recv(j,i)
          enddo
       enddo

       deallocate(recv)

       nout0 = nout0 + nout

       call MPI_ALLREDUCE(flag,totflag,1,MPI_INTEGER,MPI_SUM, &
            comm2d,ierr)
       if(totflag.eq.0) exit

    enddo

    t_ptsmv = t_ptsmv + elapsedtime_Timer(t0)

  end subroutine ptsmv1dold2_ptclmger



  !--------------------------------------------------------------------------
  ! neighboring grid communication for the 3D open boundary conditions
  ! sum up the contributions of rho from all other processors
  subroutine guardsum2dold(rhotmp,rholc,innx,inny,nytot,npy,commcol) 
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innx,inny,nytot,npy,commcol
    double precision, dimension(innx,nytot), intent(in) :: rhotmp
    double precision, dimension(innx,inny), intent(out) :: rholc
    double precision, dimension(innx*nytot) :: sendbuf1
    double precision, dimension(innx*inny) :: recvbuf1
    integer :: i,j,ii,ierr
    integer, dimension(npy) :: recvcounts
    double precision :: t0

    call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    !        if(npy.gt.1) then

    do i = 1, npy
       recvcounts(i) = innx*inny
    enddo

    do j = 1, nytot
       do i = 1, innx
          ii = (j-1)*innx + i
          sendbuf1(ii) = rhotmp(i,j)
       enddo
    enddo

    call MPI_REDUCE_SCATTER(sendbuf1,recvbuf1,recvcounts,&
         MPI_DOUBLE_PRECISION,MPI_SUM,commcol,ierr)

    do j = 1, inny
       do i = 1, innx
          ii = (j-1)*innx + i
          rholc(i,j) = recvbuf1(ii)
       enddo
    enddo

    !        else
    !          rholc = rhotmp  
    !        endif

    call MPI_BARRIER(mpicommwd,ierr)
    t_guardsum = t_guardsum + elapsedtime_Timer(t0)

  end subroutine guardsum2dold



  ! exchange grid information between neighboring guard cell
  ! to calculate E from phi. 
  subroutine guardexch2dold(phiin,phiout,innx,inny,nytot,commcol)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innx, inny, nytot,commcol
    double precision, dimension(innx,inny), intent(in) :: phiin
    double precision, dimension(innx,nytot), intent(out) :: phiout
    double precision, dimension(innx*inny) :: sendbuf1 
    double precision, dimension(innx*nytot) :: recvbuf1
    integer :: i,j,ii,ierr,nsend
    double precision :: t0

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    nsend = innx*inny

    do j = 1, inny
       do i = 1, innx
          ii = (j-1)*innx+i
          sendbuf1(ii) = phiin(i,j)
       enddo
    enddo

    call MPI_ALLGATHER(sendbuf1,nsend,MPI_DOUBLE_PRECISION,recvbuf1,&
         nsend,MPI_DOUBLE_PRECISION,commcol,ierr)

    do j = 1, nytot
       do i = 1, innx
          ii = (j-1)*innx+i
          phiout(i,j) = recvbuf1(ii)
       enddo
    enddo

    !call MPI_BARRIER(mpicommwd,ierr)
    t_guardexch = t_guardexch + elapsedtime_Timer(t0)

  end subroutine guardexch2dold



  ! sum up the contributions of rho from all other processors
  subroutine guardsum2d(rhotmp,rholc,innx,inny,nytot,npy,myidy,commcol) 
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innx,inny,npy,commcol,myidy,nytot
    double precision, dimension(innx,nytot), intent(in) :: rhotmp
    double precision, dimension(innx,inny), intent(out) :: rholc
    double precision, dimension(innx*inny) :: sendbuf1
    double precision, dimension(innx*inny) :: recvbuf1
    integer :: i,j,ii,ierr,msid,k,kk,mydes,mysou,nsend
    integer status(MPI_STATUS_SIZE)
    double precision :: t0

    call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    if(mod(nytot,npy).ne.0) then
       print*,"nytot has to be an integer multiple of npy!!"
       stop
    endif

    !Here, npy is half of total Npy since half PE is assigned to one beam.
    if(myidy.lt.npy) then

       do k = 1, inny
          do j = 1, innx
             !Nytot is uniformly distributed along npy.
             kk = myidy*inny + k
             rholc(j,k) = rhotmp(j,kk)
          enddo
       enddo

       do i = 1, npy-1
          !find the target processor
          if((myidy-i).lt.0) then
             mydes = myidy - i + npy
          else
             mydes = myidy - i
          endif
          !find the source processor
          if((myidy+i).ge.npy) then
             mysou = myidy + i - npy
          else
             mysou = myidy + i
          endif

          do k = 1, inny
             do j = 1, innx
                ii = (k-1)*innx + j
                !Nytot is uniformly distributed along npy.
                kk = mydes*inny + k
                sendbuf1(ii) = rhotmp(j,kk)
             enddo
          enddo

          nsend = inny*innx

          call MPI_IRECV(recvbuf1,nsend,MPI_DOUBLE_PRECISION,mysou,&
               i,commcol,msid,ierr)
          call MPI_SEND(sendbuf1,nsend,MPI_DOUBLE_PRECISION,mydes,&
               i,commcol,ierr)
          call MPI_WAIT(msid,status,ierr)

          do k = 1, inny
             do j = 1, innx
                ii = (k-1)*innx + j
                rholc(j,k) = rholc(j,k) + recvbuf1(ii)
             enddo
          enddo
       enddo

    else
       do k = 1, inny
          do j = 1, innx
             !Nytot is uniformly distributed along npy.
             kk = (myidy-npy)*inny + k
             rholc(j,k) = rhotmp(j,kk)
          enddo
       enddo
       do i = 1, npy-1
          !find the target processor
          if((myidy-i-npy).lt.0) then
             mydes = myidy - i + npy
          else
             mydes = myidy - i
          endif
          !find the source processor
          if((myidy+i-npy).ge.npy) then
             mysou = myidy + i - npy
          else
             mysou = myidy + i
          endif

          do k = 1, inny
             do j = 1, innx
                ii = (k-1)*innx + j
                !Nytot is uniformly distributed along npy.
                kk = (mydes-npy)*inny + k
                sendbuf1(ii) = rhotmp(j,kk)
             enddo
          enddo

          nsend = inny*innx

          call MPI_IRECV(recvbuf1,nsend,MPI_DOUBLE_PRECISION,mysou,&
               i,commcol,msid,ierr)
          call MPI_SEND(sendbuf1,nsend,MPI_DOUBLE_PRECISION,mydes,&
               i,commcol,ierr)
          call MPI_WAIT(msid,status,ierr)

          do k = 1, inny
             do j = 1, innx
                ii = (k-1)*innx + j
                rholc(j,k) = rholc(j,k) + recvbuf1(ii)
             enddo
          enddo
       enddo

    endif

    call MPI_BARRIER(mpicommwd,ierr)
    t_guardsum = t_guardsum + elapsedtime_Timer(t0)

  end subroutine guardsum2d



  ! exchange grid information between neighboring guard cell
  ! to calculate E from phi. 
  subroutine guardexch2d(phiin,phiout,innx,inny,nytot,npy,myidy,commcol)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innx,inny,nytot,commcol,myidy,npy
    double precision, dimension(innx,inny), intent(in) :: phiin
    double precision, dimension(innx,nytot), intent(out) :: phiout
    double precision, dimension(innx*inny) :: sendbuf1 
    !double precision, dimension(innx*nytot) :: recvbuf1
    double precision, dimension(innx*inny) :: recvbuf1
    integer :: i,j,ii,ierr,nsend,k,kk,msid,mysou,mydes
    integer status(MPI_STATUS_SIZE)
    double precision :: t0

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    if(mod(nytot,npy).ne.0) then
       print*,"nytot has to be an integer multiple of npy!!"
       stop
    endif

    !Here, npy is half of total Npy since half PE is assigned to one beam.
    if(myidy.lt.npy) then

       do i = 1, npy
          !find the target processor
          if((myidy-i).lt.0) then
             mydes = myidy - i + npy + npy
          else
             mydes = myidy - i + npy
          endif
          !find the source processor
          if((myidy+i).ge.npy) then
             mysou = myidy + i - npy + npy
          else
             mysou = myidy + i + npy
          endif

          do k = 1, inny
             do j = 1, innx
                ii = (k-1)*innx + j
                sendbuf1(ii) = phiin(j,k)
             enddo
          enddo

          nsend = inny*innx

          call MPI_IRECV(recvbuf1,nsend,MPI_DOUBLE_PRECISION,mysou,&
               i,commcol,msid,ierr)
          call MPI_SEND(sendbuf1,nsend,MPI_DOUBLE_PRECISION,mydes,&
               i,commcol,ierr)
          call MPI_WAIT(msid,status,ierr)

          do k = 1, inny
             do j = 1, innx
                ii = (k-1)*innx + j
                !Nytot is uniformly distributed along npy.
                kk = (mysou-npy)*inny + k
                phiout(j,kk) = recvbuf1(ii)
             enddo
          enddo
       enddo

    else
       do i = 1, npy
          !find the target processor
          if((myidy-i-npy).lt.0) then
             mydes = myidy - i + npy - npy
          else
             mydes = myidy - i - npy
          endif
          !find the source processor
          if((myidy+i-npy).ge.npy) then
             mysou = myidy + i - npy - npy
          else
             mysou = myidy + i - npy
          endif

          do k = 1, inny
             do j = 1, innx
                ii = (k-1)*innx + j
                sendbuf1(ii) = phiin(j,k)
             enddo
          enddo

          nsend = inny*innx

          call MPI_IRECV(recvbuf1,nsend,MPI_DOUBLE_PRECISION,mysou,&
               i,commcol,msid,ierr)
          call MPI_SEND(sendbuf1,nsend,MPI_DOUBLE_PRECISION,mydes,&
               i,commcol,ierr)
          call MPI_WAIT(msid,status,ierr)

          do k = 1, inny
             do j = 1, innx
                ii = (k-1)*innx + j
                !Nytot is uniformly distributed along npy.
                kk = mysou*inny + k
                phiout(j,kk) = recvbuf1(ii)
             enddo
          enddo
       enddo

    endif

    !call MPI_BARRIER(mpicommwd,ierr)
    t_guardexch = t_guardexch + elapsedtime_Timer(t0)

  end subroutine guardexch2d



  ! sum up the contributions of rho from all other processors along the row
  subroutine guardsum2drow(rhotmp,innx,inny,myjend,jend,npx,myidx,commrow) 
    implicit none
    include 'mpif.h'
    integer, intent(in) :: innx,inny,npx,commrow,myidx,myjend,jend
    double precision, dimension(innx,inny,jend), intent(inout) :: rhotmp
    double precision, dimension(innx,inny,myjend) :: rholc
    double precision, dimension(innx*inny*myjend) :: sendbuf1
    double precision, dimension(innx*inny*myjend) :: recvbuf1
    integer :: i,j,ii,ierr,msid,k,mydes,mysou,nsend,j1,mydesjend
    integer :: jj
    integer status(MPI_STATUS_SIZE)
    double precision :: t0

    call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    mydesjend = myjend
    do j1 = 1, myjend
       jj = (j1-1)*npx + myidx + 1
       do k = 1, inny
          do j = 1, innx
             rholc(j,k,j1) = rhotmp(j,k,jj)
          enddo
       enddo
    enddo

    do i = 1, npx-1
       !find the target processor
       if((myidx-i).lt.0) then
          mydes = myidx - i + npx
       else
          mydes = myidx - i
       endif
       !find the source processor
       if((myidx+i).ge.npx) then
          mysou = myidx + i - npx
       else
          mysou = myidx + i
       endif

       do j1 = 1, mydesjend
          jj = (j1-1)*npx + mydes + 1
          do k = 1, inny
             do j = 1, innx
                ii = (j1-1)*innx*inny + (k-1)*innx + j
                sendbuf1(ii) = rhotmp(j,k,jj)
             enddo
          enddo
       enddo

       nsend = inny*innx*mydesjend

       call MPI_IRECV(recvbuf1,nsend,MPI_DOUBLE_PRECISION,mysou,&
            i,commrow,msid,ierr)
       call MPI_SEND(sendbuf1,nsend,MPI_DOUBLE_PRECISION,mydes,&
            i,commrow,ierr)
       call MPI_WAIT(msid,status,ierr)

       do j1 = 1, myjend
          do k = 1, inny
             do j = 1, innx
                ii = (j1-1)*innx*inny + (k-1)*innx + j
                rholc(j,k,j1) = rholc(j,k,j1) + recvbuf1(ii)
             enddo
          enddo
       enddo

    enddo

    do j1 = 1, myjend
       do k = 1, inny
          do j = 1, innx
             rhotmp(j,k,j1) = rholc(j,k,j1)
          enddo
       enddo
    enddo

    call MPI_BARRIER(mpicommwd,ierr)
    t_guardsum = t_guardsum + elapsedtime_Timer(t0)

  end subroutine guardsum2drow



  ! collect local potential along row processors. 
  subroutine guardexch2drow(phiin,innx,inny,jend,myjend,npx,myidx,commrow)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innx,inny,jend,myjend,commrow,myidx,npx
    double precision, dimension(innx,inny,jend), intent(inout) :: phiin
    double precision, dimension(innx,inny,myjend) :: phiout
    double precision, dimension(innx*inny*myjend) :: sendbuf1 
    double precision, dimension(innx*inny*myjend) :: recvbuf1
    integer :: i,j,ierr,nsend,k,msid,mysou,mydes,jj,j1,ii
    integer status(MPI_STATUS_SIZE)
    double precision :: t0

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    do j1 = 1, myjend
       do k = 1, inny
          do j = 1, innx
             phiout(j,k,j1) = phiin(j,k,j1)
          enddo
       enddo
    enddo

    do j1 = 1, myjend
       jj = (j1-1)*npx + myidx + 1
       do k = 1, inny
          do j = 1, innx
             phiin(j,k,jj) = phiout(j,k,j1)
          enddo
       enddo
    enddo

    do i = 1, npx-1
       !find the target processor
       if((myidx-i).lt.0) then
          mydes = myidx - i + npx 
       else
          mydes = myidx - i 
       endif
       !find the source processor
       if((myidx+i).ge.npx) then
          mysou = myidx + i - npx 
       else
          mysou = myidx + i 
       endif

       do j1 = 1, myjend
          do k = 1, inny
             do j = 1, innx
                ii = (j1-1)*innx*inny + (k-1)*innx + j
                sendbuf1(ii) = phiout(j,k,j1)
             enddo
          enddo
       enddo

       nsend = inny*innx*myjend

       call MPI_IRECV(recvbuf1,nsend,MPI_DOUBLE_PRECISION,mysou,&
            i,commrow,msid,ierr)
       call MPI_SEND(sendbuf1,nsend,MPI_DOUBLE_PRECISION,mydes,&
            i,commrow,ierr)
       call MPI_WAIT(msid,status,ierr)

       do j1 = 1, myjend
          do k = 1, inny
             do j = 1, innx
                ii = (j1-1)*innx*inny + (k-1)*innx + j
                jj = (j1-1)*npx + mysou + 1
                phiin(j,k,jj) = recvbuf1(ii)
             enddo
          enddo
       enddo
    enddo

    !call MPI_BARRIER(mpicommwd,ierr)
    t_guardexch = t_guardexch + elapsedtime_Timer(t0)

  end subroutine guardexch2drow



  ! exchange the information for interpolation 
  ! from neighboring guard cells.
  subroutine boundint2d(x1,x2,innx,inny,myidy,npy,commcol) 
    implicit none
    include 'mpif.h'
    integer, intent(in) :: innx, inny,myidy,npy,commcol
    double precision, dimension(innx,inny), intent(inout) :: x1,x2
    double precision, dimension(innx,2) :: sendbuf1, recvbuf1
    !        double precision, allocatable, dimension(:,:) :: sendbuf1, recvbuf1
    integer :: bottom,top,msid,ierr,i
    integer status(MPI_STATUS_SIZE)
    double precision :: t0
    integer :: nsend

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    top = myidy + 1
    bottom = myidy - 1

    ! This could be modified to improve speed
    if(npy.gt.1) then

       nsend = innx*2
       !allocate(recvbuf1(innx,2))
       !allocate(sendbuf1(innx,2))

       if(myidy.ne.0) then
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,1,commcol,&
               msid,ierr)
       endif
       if(myidy.ne.(npy-1)) then
          do i = 1, innx
             sendbuf1(i,1) = x1(i,inny-1)
          enddo
          do i = 1, innx
             sendbuf1(i,2) = x2(i,inny-1)
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,1,&
               commcol,ierr)
       endif
       if(myidy.ne.0) then
          call MPI_WAIT(msid,status,ierr)
          do i = 1, innx
             x1(i,1) = recvbuf1(i,1)
          enddo
          do i = 1, innx
             x2(i,1) = recvbuf1(i,2)
          enddo
       endif

       if(myidy.ne.(npy-1)) then
          call MPI_IRECV(recvbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,top,2,commcol,&
               msid,ierr)
       endif
       if(myidy.ne.0) then
          do i = 1, innx
             sendbuf1(i,1) = x1(i,2)
          enddo
          do i = 1, innx
             sendbuf1(i,2) = x2(i,2)
          enddo
          call MPI_SEND(sendbuf1(1,1),nsend,MPI_DOUBLE_PRECISION,bottom,2,commcol,&
               ierr)
       endif
       if(myidy.ne.(npy-1)) then
          call MPI_WAIT(msid,status,ierr)
          do i = 1, innx
             x1(i,inny) = recvbuf1(i,1)
          enddo
          do i = 1, innx
             x2(i,inny) = recvbuf1(i,2)
          enddo
       endif

       !deallocate(recvbuf1)
       !deallocate(sendbuf1)

    endif

    t_boundint = t_boundint + elapsedtime_Timer(t0)

  end subroutine boundint2d



  !--------------------------------------------------------------------------
  ! neighboring grid communication for the 2D open boundary conditions
  ! using 1 group processors.
  ! sum up the contributions of rho from all other processors
  subroutine guardsum2d1G(rhotmp,rholc,innx,inny,nytot,npy,commcol) 
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innx,inny,nytot,npy,commcol
    double precision, dimension(innx,nytot), intent(in) :: rhotmp
    double precision, dimension(innx,inny), intent(out) :: rholc
    double precision, dimension(innx*nytot) :: sendbuf1
    double precision, dimension(innx*inny) :: recvbuf1
    integer :: i,j,ii,ierr
    integer, dimension(npy) :: recvcounts
    double precision :: t0

    call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    !        if(npy.gt.1) then

    do i = 1, npy
       recvcounts(i) = innx*inny
    enddo

    do j = 1, nytot
       do i = 1, innx
          ii = (j-1)*innx + i
          sendbuf1(ii) = rhotmp(i,j)
       enddo
    enddo

    call MPI_REDUCE_SCATTER(sendbuf1,recvbuf1,recvcounts,&
         MPI_DOUBLE_PRECISION,MPI_SUM,commcol,ierr)

    do j = 1, inny
       do i = 1, innx
          ii = (j-1)*innx + i
          rholc(i,j) = recvbuf1(ii)
       enddo
    enddo

    !        else
    !          rholc = rhotmp  
    !        endif

    call MPI_BARRIER(mpicommwd,ierr)
    t_guardsum = t_guardsum + elapsedtime_Timer(t0)

  end subroutine guardsum2d1G



  ! exchange grid information between neighboring guard cell
  ! to calculate E from phi. 
  subroutine guardexch2d1G(phiin,phiout,innx,inny,nytot,commcol)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innx, inny, nytot,commcol
    double precision, dimension(innx,inny), intent(in) :: phiin
    double precision, dimension(innx,nytot), intent(out) :: phiout
    double precision, dimension(innx*inny) :: sendbuf1 
    double precision, dimension(innx*nytot) :: recvbuf1
    integer :: i,j,ii,ierr,nsend
    double precision :: t0

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    nsend = innx*inny

    do j = 1, inny
       do i = 1, innx
          ii = (j-1)*innx+i
          sendbuf1(ii) = phiin(i,j)
       enddo
    enddo

    call MPI_ALLGATHER(sendbuf1,nsend,MPI_DOUBLE_PRECISION,recvbuf1,&
         nsend,MPI_DOUBLE_PRECISION,commcol,ierr)

    do j = 1, nytot
       do i = 1, innx
          ii = (j-1)*innx+i
          phiout(i,j) = recvbuf1(ii)
       enddo
    enddo

    !call MPI_BARRIER(mpicommwd,ierr)
    t_guardexch = t_guardexch + elapsedtime_Timer(t0)

  end subroutine guardexch2d1G


end module Communclass
