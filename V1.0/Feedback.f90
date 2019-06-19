!-----------------------------------------------------------------------------
! (c) Copyright, 2011 by the Regents of the University of California.
! Author: Stefan Paret
! Description: This module defines a feedback system following closely
!              the setup in the LHC, with 2 beam position monitors and 
!              one kicker per beam and plane.
!-----------------------------------------------------------------------------


module Feedbackclass

  use Utilityclass
  use Distributionclass
  use FFTclass

  save

  integer,parameter,private :: fbAvg = 7  !\\ number of turns to average
  integer,parameter,private :: nInput = 10  !\\ number of input parameters per pick-up

  type fbType
     double precision,private,dimension(fbAvg) :: Hilbert7,Hilbert9  ! Hilbert filter
     double precision,private,dimension(:,:),allocatable :: off7Hist,off9Hist
     double precision,private,dimension(2,2) :: R7,R9
     double precision,private,dimension(2) :: R7err,R9err
     integer,private :: delay,start,hIter,histLen  ! start flag,marker for interval most recent offsets
  endtype fbType

  type fbSystem
     type(fbType),private :: hor,ver
  endtype fbSystem


contains

  !// initialize FB for one beam and plane
  subroutine init_fb(fbs,betaIP,alphaIP,gains,myidy,npyhalf,Qx1,Qy1,Qx2,Qy2)
    implicit none
    include 'mpif.h'

    double precision,intent(in),dimension(2) :: betaIP,alphaIP
    double precision,intent(in),dimension(8) :: gains
    integer,intent(in) :: myidy,npyhalf
    type(fbSystem),intent(inout) :: fbs
    double precision,intent(in) :: Qx1,Qy1,Qx2,Qy2

    double precision,dimension(nInput) :: mpiArr
    integer :: ii,myid,ierr,fbState

    fbState = 0  !// FB off
    do ii = 1,8  !// switch FB on only if one gain is not 0.
       if(gains(ii)/=0) fbState = 1
    end do
    call MPI_COMM_RANK(mpicommwd,myid,ierr)

    if(myid==0 .and. fbState==1) open(42,file="feedback.in",status="old")
    do ii = 1,4
       if(fbState==1) then
          if(myid==0) then
             read(42,*) mpiArr(1:2)
             read(42,*) mpiArr(3:4)
             read(42,*) mpiArr(5)
             read(42,*) mpiArr(6:7)
             read(42,*) mpiArr(8:9)
             read(42,*) mpiArr(10)
          endif
          call MPI_BCAST(mpiArr,nInput,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
       endif

       if(myidy<npyhalf) then  !// associate parameters with proper beam and plane
          if(ii==1) call fb_initInternal(fbs%hor,fbState,mpiArr,betaIP(1),alphaIP(1),gains(1:2),Qx1)
          if(ii==2) call fb_initInternal(fbs%ver,fbState,mpiArr,betaIP(2),alphaIP(2),gains(3:4),Qy1)
       else
          if(ii==3) call fb_initInternal(fbs%hor,fbState,mpiArr,betaIP(1),alphaIP(1),gains(5:6),Qx2)
          if(ii==4) call fb_initInternal(fbs%ver,fbState,mpiArr,betaIP(2),alphaIP(2),gains(7:8),Qy2)
       end if
    end do
    if(myid==0 .and. fbState==1) close(42)

  end subroutine init_fb



  subroutine fb_initInternal(fbObj,fbState,mpiArr,betaIP,alphaIP,gains,tune)
    implicit none
    double precision,parameter :: pi = acos(-1.d0)

    integer,intent(in) :: fbState
    double precision,intent(in),dimension(nInput) :: mpiArr
    double precision,intent(in),dimension(2) :: gains
    double precision,intent(in) :: betaIP,alphaIP,tune
    type(fbType),intent(out) :: fbObj

    double precision :: phik7,phik9,phiH7,phiH9,muIP7,muIP9,a_0,theta,beta7,beta9,err7,err9

    if(fbState==1) then
       phik7 = mpiArr(1)/180.*pi  !// phases from pick-ups to kicker
       phik9 = mpiArr(2)/180.*pi
       phiH7 = mpiArr(3)/180.*pi  !// phases for Hilbert filter
       phiH9 = mpiArr(4)/180.*pi
       fbObj%delay = int(mpiArr(5))
       err7 = mpiArr(6)  !// rms measuring error
       err9 = mpiArr(7)
       beta7 = mpiArr(8)  !// beta function at pick-up
       beta9 = mpiArr(9)
       muIP7 = mpiArr(10)/180.*pi  !// phase advance from IP to pick-up 7
       fbObj%start = fbObj%delay

       fbObj%histLen = fbAvg+fbObj%delay+1
       allocate(fbObj%off7Hist(fbObj%histLen,2))
       allocate(fbObj%off9Hist(fbObj%histLen,2))
       fbObj%off7Hist = 0.
       fbObj%off9Hist = 0.

       muIP9 = muIP7 - phik9 + phik7  ! phase advance from IP to pick-up 9

       !// fill Hilbert filters
       fbObj%Hilbert7(6) = 0.
       fbObj%Hilbert7(5) = 2./pi*sin(phiH7)
       fbObj%Hilbert7(4) = cos(phiH7)
       fbObj%Hilbert7(7) = fbObj%Hilbert7(5)/3.  ! associated with most recent datum
       fbObj%Hilbert7(3) = -fbObj%Hilbert7(5)
       fbObj%Hilbert7(2) = 0.
       fbObj%Hilbert7(1) = -fbObj%Hilbert7(7)  ! associated with least recent datum
       fbObj%Hilbert9(6) = 0.
       fbObj%Hilbert9(5) = 2./pi*sin(phiH9)
       fbObj%Hilbert9(4) = cos(phiH9)
       fbObj%Hilbert9(7) = fbObj%Hilbert9(5)/3.
       fbObj%Hilbert9(3) = -fbObj%Hilbert9(5)
       fbObj%Hilbert9(2) = 0.
       fbObj%Hilbert9(1) = -fbObj%Hilbert9(7)

       !// Coefficients for transformation of offset and error to effective change of coordinates at previoud IP
       theta = 2.d0*pi*tune
       a_0 = 1./sqrt(2*(1-cos(theta))*(cos(phiH7)**2 + 8./(9.*pi**2)*sin(phiH7)**2* &
            (10.-3.*cos(2*theta)-6*cos(4*theta)-cos(6*theta))))  !// amplitude scaling for Hilbert notch filter (see Zhabitsky)
       fbObj%R7(1,1) = gains(1)*a_0*(cos(muIP7)+alphaIP*sin(muIP7))*sin(muIP7+phik7)
       fbObj%R7(1,2) = gains(1)*a_0*betaIP*sin(muIP7)*sin(muIP7+phik7)
       fbObj%R7(2,1) = -gains(1)*a_0*(cos(muIP7)+alphaIP*sin(muIP7))*(cos(muIP7+phik7) + &
            alphaIP*sin(muIP7+phik7))/betaIP
       fbObj%R7(2,2) = -gains(1)*a_0*sin(muIP7)*(cos(muIP7+phik7)+alphaIP*sin(muIP7+phik7))
       a_0 = 1./sqrt(2*(1-cos(theta))*(cos(phiH7)**2 + 8./(9.*pi**2)*sin(phiH7)**2* &
            (10.-3.*cos(2*theta)-6*cos(4*theta)-cos(6*theta))))
       fbObj%R9(1,1) = gains(2)*a_0*(cos(muIP9)+alphaIP*sin(muIP9))*sin(muIP9+phik9)
       fbObj%R9(1,2) = gains(2)*a_0*betaIP*sin(muIP9)*sin(muIP9+phik9)
       fbObj%R9(2,1) = -gains(2)*a_0*(cos(muIP9)+alphaIP*sin(muIP9))*(cos(muIP9+phik9) + &
            alphaIP*sin(muIP9+phik9))/betaIP
       fbObj%R9(2,2) = -gains(2)*a_0*sin(muIP9)*(cos(muIP9+phik9)+alphaIP*sin(muIP9+phik9))
       !// transfer matrix elements for pick-up errors
       fbObj%R7err(1) = -gains(1)*a_0*sqrt(betaIP/beta7)*sin(phik7)*err7
       fbObj%R7err(2) = gains(1)*a_0/sqrt(beta7*betaIP)*(cos(phik7)+alphaIP*sin(phik7))*err7
       fbObj%R9err(1) = -gains(1)*a_0*sqrt(betaIP/beta9)*sin(phik9)*err9
       fbObj%R9err(2) = gains(1)*a_0/sqrt(beta9*betaIP)*(cos(phik9)+alphaIP*sin(phik9))*err9
    else  !// set everything (except start) to 0 if FB is off
       fbObj%Hilbert7 = 0.
       fbObj%Hilbert9 = 0.
       fbObj%delay = 0
       fbObj%start = -1
    endif
    fbObj%hIter = iToIter(fbObj,1)  !// Initialize iterator

  end subroutine fb_initInternal



  !// Apply feedback kick
  subroutine fb_kick(fbs,Bpts,Nplocal,Np,myidy,npyhalf,errSeed)
    implicit none
    include 'mpif.h'

    integer,intent(in) :: Nplocal,Np,myidy,npyhalf
    integer,intent(inout) :: errSeed
    type(fbSystem),intent(inout) :: fbs
    double precision,pointer,dimension(:,:),intent(inout) :: Bpts

    type(fbType) :: fbObj
    double precision,dimension(6) :: off1,off2,off
    double precision,dimension(4) :: normErrs
    double precision,dimension(2,2) :: tempErr
    double precision,dimension(2) :: offCorr,myErr

    integer :: i,j,ind,myid,ierr


    call MPI_COMM_RANK(mpicommwd,myid,ierr)
    call findcentroids_Utility(Bpts,Nplocal,Np,myidy,npyhalf,off1,off2)

    do j = 1,2  !// process first horizontal plane,then vertical
       !// Provide random numbers for pick-up noise
       if(myid==0) call normVecnew(errSeed,tempErr,2)
       normErrs(1:2) = tempErr(1,1:2)
       normErrs(3:4) = tempErr(2,1:2)
       call MPI_BCAST(normErrs,4,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)

       !// Initialize local variables depending on beam
       if(j==1) then  !// horizontal plane
          fbObj = fbs%hor  !// horizontal plane
          if(myidy<npyhalf) then
             myErr = normErrs(1:2)
             off = off1  !// beam 1
          else
             myErr = normErrs(3:4)
             off = off2  !// beam 2
          endif
          ind = 1  !// point to x
       else  !// vertical plane
          fbObj = fbs%ver  !// vertical plane
          if(myidy<npyhalf) then
             myErr = normErrs(1:2)
             off = off1
          else
             myErr = normErrs(3:4)
             off = off2
          endif
          ind = 3  !// point to y
       endif

       if(fbObj%start>=0) then
          !// Insert current offset (transformed back to IP) into FB memory
          fbObj%off7Hist(iToIter(fbObj,fbObj%hIter-2),:) = &
               matmul(fbObj%R7,off(ind:ind+1)) + fbObj%R7err*myErr(1)
          fbObj%off9Hist(iToIter(fbObj,fbObj%hIter-2),:) = &
               matmul(fbObj%R9,off(ind:ind+1)) + fbObj%R9err*myErr(2)

          !// Determine correction (transformed back to IP)
          offCorr = 0.
          do i=1,fbAvg
             offCorr = offCorr + fbObj%Hilbert7(i) * &
                  (fbObj%off7Hist(iToIter(fbObj,i+fbObj%hIter-1),:) - &
                  fbObj%off7Hist(iToIter(fbObj,i+fbObj%hIter-2),:))
          enddo
          do i=1,fbAvg
             offCorr = offCorr + fbObj%Hilbert9(i) * &
                  (fbObj%off9Hist(iToIter(fbObj,i+fbObj%hIter-1),:) - &
                  fbObj%off9Hist(iToIter(fbObj,i+fbObj%hIter-2),:))
          enddo

          !// Update iterator
          fbObj%hIter = iToIter(fbObj,fbObj%hIter+1)

          !// Apply correction
          if(fbObj%start==0) then  !// only after countdown finished
             do i = 1,Nplocal
                Bpts(ind:ind+1,i) = Bpts(ind:ind+1,i) + offCorr
             enddo
          else  !// count down to beginning of damping
             fbObj%start = fbObj%start - 1
          endif

          !// Update persistent feedback object
          if(j==1) then !// horizontal plane
             fbs%hor = fbObj
          else  !// vertical plane
             fbs%ver = fbObj
          endif
       endif
    end do

  end subroutine fb_kick



  subroutine destroy_fb(fbs)
    type(fbSystem),intent(inout) :: fbs

    if(fbs%hor%start>=0) then
       deallocate(fbs%hor%off7Hist)
       deallocate(fbs%hor%off9Hist)
    endif
    if(fbs%ver%start>=0) then
       deallocate(fbs%ver%off7Hist)
       deallocate(fbs%ver%off9Hist)
    endif

  end subroutine destroy_fb



  integer function iToIter(fbObj,i)  !// set iterator for offset history
    type(fbType),intent(in) :: fbObj
    integer,intent(in) :: i

    iToIter = i
    do
       if(iToIter>0) exit
       iToIter = iToIter + fbObj%histLen
    enddo
    iToIter = mod(iToIter-1,fbObj%histLen) + 1

  end function iToIter

      !find the centroid of each slice, apply the band-pass filter, and damping
      subroutine feedbackslice(pts,innp,nslice,npyhalf,myidy,frqlowx,frqhighx,frqlowy,frqhighy,&
                               gainx,gainy,close2g) 
      implicit none
      include 'mpif.h'
      integer, intent(in) :: innp,nslice,npyhalf,myidy
      real*8, pointer, dimension(:,:) :: pts
      real*8, intent(in) :: frqlowx,frqhighx,frqlowy,frqhighy,gainx,gainy
      real*8, intent(in), dimension(2) :: close2g
      real*8, dimension(nslice) :: xx,yy,count
      real*8, dimension(nslice,6) :: tmparry,tmpglb
      integer :: i,iz,iz1,nsend,ierr,my_rank
      real*8 :: zmin,zmax,hz,zz,zavg,zmingl,zmaxgl
      integer :: ksign,ilow,ihigh,ittp
      double complex, dimension(nslice/2+1,1) :: ctmpx,ctmpx2
      real*8, dimension(nslice,1) :: tmpct
      real*8 :: xs,ys,deltafrq,ab
      real*8 :: scale

      zmin = 1.e12
      zmax = -1.e12
      do i = 1, innp
        if(zmin.ge.pts(5,i)) then
          zmin = pts(5,i)
        endif
        if(zmax.le.pts(5,i)) then
          zmax = pts(5,i)
        endif
      enddo
      call MPI_ALLREDUCE(zmin,zmingl,1,MPI_DOUBLE_PRECISION,&
                            MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(zmax,zmaxgl,1,MPI_DOUBLE_PRECISION,&
                            MPI_MAX,MPI_COMM_WORLD,ierr)

      !avoid index overflow
      zmin = zmingl - (zmaxgl-zmingl)*0.5e-5
      zmax = zmaxgl + (zmaxgl-zmingl)*0.5e-5
      hz = (zmax-zmin)/(nslice-1)

      deltafrq = 1.0d0/(hz*nslice)
      !print*,"zmin, zmax and hz: ",zmin*scxlt,zmax*scxlt,hz*scxlt
      count = 0.0
      xx = 0.0
      yy = 0.0

      !calculate the mean of each slice by linear deposition
      do i = 1, innp
        iz = (pts(5,i)-zmin)/hz + 1
        if(iz.lt.1) iz = 1
        if(iz.ge.nslice) iz = nslice-1
        iz1 = iz + 1
        ab = ((zmin-pts(5,i))+iz*hz)/hz
        count(iz) = count(iz) + ab
        count(iz1) = count(iz1) + 1.0d0 - ab
        xx(iz) = xx(iz) + pts(1,i)*ab
        xx(iz1) = xx(iz1) + pts(1,i)*(1.0d0-ab) 
        yy(iz) = yy(iz) + pts(3,i)*ab
        yy(iz1) = yy(iz1) + pts(3,i)*(1.0d0-ab)
      enddo 
  
      if(myidy.lt.npyhalf) then
        tmparry(:,1) = count
        tmparry(:,2) = xx
        tmparry(:,3) = yy
        tmparry(:,4) = 0.0d0
        tmparry(:,5) = 0.0d0
        tmparry(:,6) = 0.0d0
      else
        tmparry(:,1) = 0.0d0
        tmparry(:,2) = 0.0d0
        tmparry(:,3) = 0.0d0
        tmparry(:,4) = count
        tmparry(:,5) = xx
        tmparry(:,6) = yy
      endif

      nsend = nslice*6

      tmpglb = 0.0d0
      call MPI_ALLREDUCE(tmparry,tmpglb,nsend,MPI_DOUBLE_PRECISION,&
                            MPI_SUM,MPI_COMM_WORLD,ierr)
      if(myidy.lt.npyhalf) then
        count = tmpglb(:,1)
        xx = tmpglb(:,2)
        yy = tmpglb(:,3)
      else
        count = tmpglb(:,4)
        xx = tmpglb(:,5)
        yy = tmpglb(:,6)
      endif

      !print*,"sum count: ",sum(count)
      do i = 1, nslice
        if(count(i).gt.0.0d0) then
          xx(i) = xx(i)/count(i)
          yy(i) = yy(i)/count(i)
        else
          xx(i) = 0.0d0
          yy(i) = 0.0d0
        endif
      enddo

      !go to freq. domain
      ksign = 1
      scale = 1.0d0
      !print*,"sumxx1:",sum(xx)
      !if(myidy.eq.npyhalf) then
      !do i = 1, nslice
      !  write(1,*)i,xx(i)
      !enddo
      !endif
      tmpct(:,1) = xx(:)
      call fftrclocal1_FFT(ksign,scale,tmpct,nslice,1,ctmpx)

      !apply the band-pass filter
      ilow =  frqlowx/deltafrq + 1
      ittp = frqhighx/deltafrq + 1
      ihigh =  min(ittp,nslice/2+1)

      !print*,"frqlowhxx: ",frqhighx,deltafrq,ilow,ihigh,gainx
      ctmpx2 = dcmplx(0.0d0,0.0d0)
      ctmpx2(ilow:ihigh,1) = ctmpx(ilow:ihigh,1)

      !print*,"ctmpx: ",sum(real(ctmpx(:,1))),sum(aimag(ctmpx(:,1)))
      !print*,"ctmpx2: ",sum(real(ctmpx2(:,1))),sum(aimag(ctmpx2(:,1)))
      !back to time domain
      !if(myidy.eq.npyhalf) then
      !do i = 1, nslice/2+1
      !  write(3,123)i*1.0d0,real(ctmpx(i,1)),aimag(ctmpx(i,1)),&
      !              real(ctmpx2(i,1)),aimag(ctmpx2(i,1))
      !enddo
      !endif
!123   format(5(1x,e16.8))

      tmpct = 0.0d0
      ksign = -1
      scale = 1.0d0/nslice
      call fftcrlocal1_FFT(ksign,scale,ctmpx2,nslice,1,tmpct)
      xx(:) = tmpct(:,1)
      !print*,"sumxx2:",sum(xx)
      !if(myidy.eq.npyhalf) then
      !do i = 1, nslice
      !  write(2,*)i,xx(i)
      !enddo
      !endif

      tmpct(:,1) = yy(:)
      !go to freq. domain
      ksign = 1
      scale = 1.0d0
      call fftrclocal1_FFT(ksign,scale,tmpct,nslice,1,ctmpx)

      !apply the band-pass filter
      ilow =  frqlowy/deltafrq + 1
      !ihigh =  min(frqhighy/deltafrq+1,nslice/2+1)
      ittp = frqhighy/deltafrq + 1
      ihigh =  min(ittp,nslice/2+1)
      !print*,"frqlowhyy: ",frqhighy,deltafrq,ilow,ihigh,gainy
      ctmpx2 = dcmplx(0.0d0,0.0d0)
      ctmpx2(ilow:ihigh,1) = ctmpx(ilow:ihigh,1)

      !back to time domain
      tmpct = 0.0d0
      ksign = -1
      scale = 1.0d0/nslice
      call fftcrlocal1_FFT(ksign,scale,ctmpx2,nslice,1,tmpct)
      yy(:) = tmpct(:,1)

      if(myidy.lt.npyhalf) then
        xx = 0.0d0
        yy = 0.0d0
      endif

      !apply feedback
      do i = 1, innp
        iz = (pts(5,i)-zmin)/hz + 1
        if(iz.lt.1) iz = 1
        if(iz.ge.nslice) iz = nslice-1
        iz1 = iz + 1
        ab = ((zmin-pts(5,i))+iz*hz)/hz
        xs = xx(iz)*(1.0d0-ab) + xx(iz1)*ab
        pts(1,i) = pts(1,i) - gainx*xs + close2g(1)
        ys = yy(iz)*(1.0d0-ab) + yy(iz1)*ab
        pts(3,i) = pts(3,i) - gainy*ys + close2g(2)
      enddo


      !call flush(1)
      !call flush(2)
      !call flush(3)

      end subroutine feedbackslice

end module feedbackclass
