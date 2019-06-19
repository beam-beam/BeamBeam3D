!----------------------------------------------------------------
! (c) Copyright, 2011 by the Regents of the University of California.
! Distributionclass: Initial distribution of charged beam bunch class in 
!                    Beam module of APPLICATION layer.
! Version: 2.6
! Author: Ji Qiang 
! Description: This class defines initial distributions for the charged 
!              particle beam bunch information in the accelerator.
! Comments:
!----------------------------------------------------------------

module Distributionclass
  use Pgrid2dclass
  use Timerclass
  use Utilityclass


contains
  ! sample the particles with intial distribution.
  subroutine sample_Dist(Pts1,distparam,nparam,flagdist,grid,&
       Npt,Nptlocal,ngroup,ibunch)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nparam,Npt,Nptlocal,ngroup,ibunch
    double precision, dimension(nparam) :: distparam
    double precision, pointer, dimension(:,:) :: Pts1
    type (Pgrid2d), intent(in) :: grid
    integer, intent(in) :: flagdist
    integer :: myid, myidx, myidy,seedsize,i,npyhalf,npy,npx,&
         totnp,nptot
    !integer seedarray(1)
    integer, allocatable, dimension(:) :: seedarray
    real*8 rancheck
    integer :: iseed,meanpts20


    call getpost_Pgrid2d(grid,myid,myidy,myidx)
    call getsize_Pgrid2d(grid,totnp,npy,npx)
    npyhalf = npy/ngroup
    meanpts20 = (Npt/totnp)*20
    !        seedarray(1)=(100001+myid)*(myid+7)
    !        call random_seed(put=seedarray(1:1))
    !        write(6,*)'seedarray=',seedarray

    call random_seed(SIZE=seedsize)
    allocate(seedarray(seedsize))
    do i = 1, seedsize
       !seedarray(i) = (1000+5*myid)*(myid+7)+i-1+300*(ibunch-1)
       !          seedarray(i) = (1000+50*myid)*(10*myid+7)+i-1+300*(ibunch-1)
       !          seedarray(i) = (1000+5*myid)*(myid+7)+i-1
       !seedarray(i) = (1000)*(7)+i-1
       seedarray(i) = int(10.0d0 + myid*1.0d0*meanpts20+i*1.0d0*myid +300.0d0*(ibunch-1))
    enddo
    call random_seed(PUT=seedarray)
    call random_number(rancheck)
    !warm up the random number generator
    !do i = 1, 3000
    do i = 1, 3000*(myid+1)
       call random_number(rancheck)
    enddo
    !        write(6,*)'myid,rancheck=',seedarray,myid,rancheck

    if(flagdist.eq.1) then
       print*,"not available yet!"
    else if(flagdist.eq.2) then
       iseed = -((1000+5*myid)*(myid+7)+i-1+300*(ibunch-1))
       call Gauss3_Dist(Pts1,nparam,distparam,grid,0,Npt,Nptlocal,ngroup,iseed)
       !call Gauss3Sob_Dist(Pts1,nparam,distparam,grid,0,Npt,Nptlocal,ngroup,iseed)
       !call Gauss3test_Dist(Pts1,nparam,distparam,grid,0,Npt,Nptlocal,&
       !ngroup)
    else if(flagdist.eq.3) then
       call Waterbag_Dist(Pts1,nparam,distparam,grid,0,Npt,Nptlocal,ngroup)
    else if(flagdist.eq.4) then
       print*,"not available yet!"
    else if(flagdist.eq.5) then
       call KV3d_Dist(Pts1,nparam,distparam,grid,Npt,Nptlocal,ngroup)
    else if(flagdist.eq.6) then
       print*,"not available yet!"
    else if(flagdist.eq.7) then
       print*,"not available yet!"
    else if(flagdist.eq.8) then
       print*,"not available yet!"
    else if(flagdist.eq.9) then
       print*,"not available yet!"
    else if(flagdist.eq.10) then
       print*,"not available yet!"
    else if(flagdist.eq.11) then
       if(ngroup.eq.2) then
          call read2g_Dist(Pts1,nparam,distparam,grid,0,Npt,Nptlocal,ngroup)
          !else if(ngroup.eq.1) then
       else
          print*,"wrong with ngroup: ",ngroup
       endif
    else if(flagdist.eq.12) then
       call regen_Dist(Pts1,myidx,myidy,npx,npyhalf,nptot)
       if(nptot.ne.Npt) then
          print*,"please check the total number of particles in both files"
          stop
       endif
    else
       print*,"Initial distribution not available!!"
       stop
    endif

    deallocate(seedarray)

  end subroutine sample_Dist



  subroutine Gauss3_Dist(Pts1,nparam,distparam,grid,flagalloc,&
       Npt,Nptlocal,ngroup,iseed)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: nparam,flagalloc,Npt,Nptlocal,ngroup
    double precision, dimension(nparam) :: distparam
    type (Pgrid2d), intent(in) :: grid
    integer, intent(inout) :: iseed
    double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
         sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
    double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
    double precision :: sig1,sig2,sig3,sig4,sig5,sig6
    double precision :: sq12,sq34,sq56
    double precision, allocatable, dimension(:,:) :: x1,x2,x3 
    integer :: totnp,npy,npx
    integer :: avgpts
    integer :: myid,myidx,myidy,i,j,k,intvsamp
    !        integer seedarray(1)
    double precision :: t0,x11
    double precision :: xrms,pxrms,yrms,pyrms,xavg,x2avg,pxavg,px2avg,&
         yavg,y2avg,pyavg,py2avg
    double precision, dimension(16) :: tmplc,tmpgl
    integer :: npyhalf,ierr

    call starttime_Timer(t0)

    sigx = distparam(1)
    sigpx = distparam(2)
    muxpx = distparam(3)
    xscale = distparam(4)
    pxscale = distparam(5)
    xmu1 = distparam(6)
    xmu2 = distparam(7)
    sigy = distparam(8)
    sigpy = distparam(9)
    muypy = distparam(10)
    yscale = distparam(11)
    pyscale = distparam(12)
    xmu3 = distparam(13)
    xmu4 = distparam(14)
    sigz = distparam(15)
    sigpz = distparam(16)
    muzpz = distparam(17)
    zscale = distparam(18)
    pzscale = distparam(19)
    xmu5 = distparam(20)
    xmu6 = distparam(21)

    call getsize_Pgrid2d(grid,totnp,npy,npx)

    call getpost_Pgrid2d(grid,myid,myidy,myidx)
    !        seedarray(1)=(1001+myid)*(myid+7)
    !        write(6,*)'seedarray=',seedarray
    !        call random_seed(PUT=seedarray(1:1))
    call random_number(x11)

    avgpts = Npt/(npx*npy/ngroup)

    sig1 = sigx*xscale
    sig2 = sigpx*pxscale
    sig3 = sigy*yscale
    sig4 = sigpy*pyscale
    sig5 = sigz*zscale
    sig6 = sigpz*pzscale

    sq12=sqrt(1.-muxpx*muxpx)
    sq34=sqrt(1.-muypy*muypy)
    sq56=sqrt(1.-muzpz*muzpz)

    ! initial allocate 'avgpts' particles on each processor.
    !if(flagalloc.eq.1) then
    !  Pts1 = 0.0
    !else
    !  allocate(Pts1(6,avgpts))
    !  Pts1 = 0.0
    !endif

    !        allocate(x1(2,avgpts))
    !        allocate(x2(2,avgpts))
    !        allocate(x3(2,avgpts))
    !        call normVec(x1,avgpts)
    !        call normVec(x2,avgpts)
    !        call normVec(x3,avgpts)

    !intvsamp = 10
    intvsamp = avgpts
    allocate(x1(2,intvsamp))
    allocate(x2(2,intvsamp))
    allocate(x3(2,intvsamp))

    xavg = 0.0
    x2avg = 0.0
    pxavg = 0.0
    px2avg = 0.0
    yavg = 0.0
    y2avg = 0.0
    pyavg = 0.0
    py2avg = 0.0
    do j = 1, avgpts/intvsamp
       !call normVecnew(iseed,x1,intvsamp)
       !call normVecnew(iseed,x2,intvsamp)
       !call normVecnew(iseed,x3,intvsamp)
       call normVec(x1,intvsamp)
       call normVec(x2,intvsamp)
       call normVec(x3,intvsamp)
       do k = 1, intvsamp
          !x-px:
          !            call normdv(x1)
          !           Correct Gaussian distribution.
          i = (j-1)*intvsamp + k
          !Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
          !Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
          !miguel's sampling. here, muxpx->alphax,sigpx->betax
          Pts1(1,i) = xmu1 + sigx*x1(1,k)*xscale
          xavg = xavg + Pts1(1,i)
          x2avg = x2avg + Pts1(1,i)*Pts1(1,i)
          Pts1(2,i) = xmu2 + sigx  &
               *(-muxpx*x1(1,k)+x1(2,k))/sigpx*pxscale
          pxavg = pxavg + Pts1(2,i)
          px2avg = px2avg + Pts1(2,i)*Pts1(2,i)
          !y-py
          !            call normdv(x1)
          !           Correct Gaussian distribution.
          !Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
          !Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))
          !miguel's sampling. here, muxpx->alphax,sigpx->betax
          Pts1(3,i) = xmu3 + sigy*x2(1,k)*yscale
          yavg = yavg + Pts1(3,i)
          y2avg = y2avg + Pts1(3,i)*Pts1(3,i)
          Pts1(4,i) = xmu4 + sigy  &
               *(-muypy*x2(1,k)+x2(2,k))/sigpy*pyscale
          pyavg = pyavg + Pts1(4,i)
          py2avg = py2avg + Pts1(4,i)*Pts1(4,i)
          !z-pz
          !            call normdv(x1)
          !           Correct Gaussian distribution.
          !Pts1(5,i) = xmu5 + sig5*x3(1,k)/sq56
          !Pts1(6,i) = xmu6 + sig6*(-muzpz*x3(1,k)/sq56+x3(2,k))
          !miguel's sampling. here, muxpx->alphax,sigpx->betax
          Pts1(5,i) = xmu5 + x3(1,k)*sig5
          Pts1(6,i) = xmu6 + x3(2,k)*sig6
       enddo
    enddo
    npyhalf = npy/ngroup
    tmplc = 0.0d0
    if(myidy.lt.npyhalf) then
       tmplc(1) = xavg
       tmplc(2) = x2avg
       tmplc(3) = pxavg
       tmplc(4) = px2avg
       tmplc(5) = yavg
       tmplc(6) = y2avg
       tmplc(7) = pyavg
       tmplc(8) = py2avg
    else
       tmplc(9) = xavg
       tmplc(10) = x2avg
       tmplc(11) = pxavg
       tmplc(12) = px2avg
       tmplc(13) = yavg
       tmplc(14) = y2avg
       tmplc(15) = pyavg
       tmplc(16) = py2avg
    endif

    call MPI_ALLREDUCE(tmplc,tmpgl,16,MPI_DOUBLE_PRECISION,&
         MPI_SUM,mpicommwd,ierr)

    do i = 1, 16
       tmpgl(i) = tmpgl(i)/Npt
    enddo

    if(myidy.lt.npyhalf) then
       xrms = sqrt(tmpgl(2)-tmpgl(1)**2)
       pxrms = sqrt(tmpgl(4)-tmpgl(3)**2)
       yrms = sqrt(tmpgl(6)-tmpgl(5)**2)
       pyrms = sqrt(tmpgl(8)-tmpgl(7)**2)
    else
       xrms = sqrt(tmpgl(10)-tmpgl(9)**2)
       pxrms = sqrt(tmpgl(12)-tmpgl(11)**2)
       yrms = sqrt(tmpgl(14)-tmpgl(13)**2)
       pyrms = sqrt(tmpgl(16)-tmpgl(15)**2)
    endif

    !        do j = 1, avgpts
    !          Pts1(1,j) = Pts1(1,j)*sigx/xrms
    !          Pts1(2,j) = Pts1(2,j)*(sigx/sigpx)/pxrms
    !          Pts1(3,j) = Pts1(3,j)*sigy/yrms
    !          Pts1(4,j) = Pts1(4,j)*(sigy/sigpy)/pyrms
    !        enddo

    goto 1000

!    Pts1(1,avgpts+1)= xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+1)= xmu2 + 0.0
!    Pts1(3,avgpts+1) = xmu3 + 5.0*sigy
!    Pts1(4,avgpts+1) = xmu4 + 0.0
!    Pts1(5,avgpts+1) = 0.0
!    Pts1(6,avgpts+1) = 0.0
!    Pts1(1,avgpts+2) = xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+2) = xmu2 + 0.0
!    Pts1(3,avgpts+2) = xmu3 - 5.0*sigy
!    Pts1(4,avgpts+2) = xmu4 + 0.0
!    Pts1(5,avgpts+2) = 0.0
!    Pts1(6,avgpts+2) = 0.0
!    Pts1(1,avgpts+3) = xmu1 + 5.0*sigx 
!    Pts1(2,avgpts+3) = xmu2 + 0.0
!    Pts1(3,avgpts+3) = xmu3 - 0.0*sigy
!    Pts1(4,avgpts+3) = xmu4 + 0.0
!    Pts1(5,avgpts+3) = 0.0
!    Pts1(6,avgpts+3) = 0.0
!    Pts1(1,avgpts+4) = xmu1 - 5.0*sigx 
!    Pts1(2,avgpts+4) = xmu2 + 0.0
!    Pts1(3,avgpts+4) = xmu3 + 0.0*sigy
!    Pts1(4,avgpts+4) = xmu4 + 0.0
!    Pts1(5,avgpts+4) = 0.0
!    Pts1(6,avgpts+4) = 0.0
!    Pts1(1,avgpts+5) = xmu1 + 5.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+5) = xmu2 + 0.0
!    Pts1(3,avgpts+5) = xmu3 + 5.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+5) = xmu4 + 0.0
!    Pts1(5,avgpts+5) = 0.0
!    Pts1(6,avgpts+5) = 0.0
!    Pts1(1,avgpts+6) = xmu1 - 5.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+6) = xmu2 + 0.0
!    Pts1(3,avgpts+6) = xmu3 - 5.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+6) = xmu4 + 0.0
!    Pts1(5,avgpts+6) = 0.0
!    Pts1(6,avgpts+6) = 0.0
!
!    Pts1(1,avgpts+7)= xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+7)= xmu2 + 0.0
!    Pts1(3,avgpts+7) = xmu3 + 4.0*sigy
!    Pts1(4,avgpts+7) = xmu4 + 0.0
!    Pts1(5,avgpts+7) = 0.0
!    Pts1(6,avgpts+7) = 0.0
!    Pts1(1,avgpts+8) = xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+8) = xmu2 + 0.0
!    Pts1(3,avgpts+8) = xmu3 - 4.0*sigy
!    Pts1(4,avgpts+8) = xmu4 + 0.0
!    Pts1(5,avgpts+8) = 0.0
!    Pts1(6,avgpts+8) = 0.0
!    Pts1(1,avgpts+9) = xmu1 + 4.0*sigx 
!    Pts1(2,avgpts+9) = xmu2 + 0.0
!    Pts1(3,avgpts+9) = xmu3 - 0.0*sigy
!    Pts1(4,avgpts+9) = xmu4 + 0.0
!    Pts1(5,avgpts+9) = 0.0
!    Pts1(6,avgpts+9) = 0.0
!    Pts1(1,avgpts+10) = xmu1 - 4.0*sigx 
!    Pts1(2,avgpts+10) = xmu2 + 0.0
!    Pts1(3,avgpts+10) = xmu3 + 0.0*sigy
!    Pts1(4,avgpts+10) = xmu4 + 0.0
!    Pts1(5,avgpts+10) = 0.0
!    Pts1(6,avgpts+10) = 0.0
!    Pts1(1,avgpts+11) = xmu1 - 4.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+11) = xmu2 + 0.0
!    Pts1(3,avgpts+11) = xmu3 + 4.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+11) = xmu4 + 0.0
!    Pts1(5,avgpts+11) = 0.0
!    Pts1(6,avgpts+11) = 0.0
!    Pts1(1,avgpts+12) = xmu1 - 4.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+12) = xmu2 + 0.0
!    Pts1(3,avgpts+12) = xmu3 - 4.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+12) = xmu4 + 0.0
!    Pts1(5,avgpts+12) = 0.0
!    Pts1(6,avgpts+12) = 0.0
!    Pts1(1,avgpts+13) = xmu1 + 4.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+13) = xmu2 + 0.0
!    Pts1(3,avgpts+13) = xmu3 + 4.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+13) = xmu4 + 0.0
!    Pts1(5,avgpts+13) = 0.0
!    Pts1(6,avgpts+13) = 0.0
!    Pts1(1,avgpts+14) = xmu1 + 4.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+14) = xmu2 + 0.0
!    Pts1(3,avgpts+14) = xmu3 - 4.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+14) = xmu4 + 0.0
!    Pts1(5,avgpts+14) = 0.0
!    Pts1(6,avgpts+14) = 0.0

!    Pts1(1,avgpts+15)= xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+15)= xmu2 + 0.0
!    Pts1(3,avgpts+15) = xmu3 + 3.0*sigy
!    Pts1(4,avgpts+15) = xmu4 + 0.0
!    Pts1(5,avgpts+15) = 0.0
!    Pts1(6,avgpts+15) = 0.0
!    Pts1(1,avgpts+16) = xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+16) = xmu2 + 0.0
!    Pts1(3,avgpts+16) = xmu3 - 3.0*sigy
!    Pts1(4,avgpts+16) = xmu4 + 0.0
!    Pts1(5,avgpts+16) = 0.0
!    Pts1(6,avgpts+16) = 0.0
!    Pts1(1,avgpts+17) = xmu1 + 3.0*sigx 
!    Pts1(2,avgpts+17) = xmu2 + 0.0
!    Pts1(3,avgpts+17) = xmu3 - 0.0*sigy
!    Pts1(4,avgpts+17) = xmu4 + 0.0
!    Pts1(5,avgpts+17) = 0.0
!    Pts1(6,avgpts+17) = 0.0
!    Pts1(1,avgpts+18) = xmu1 - 3.0*sigx 
!    Pts1(2,avgpts+18) = xmu2 + 0.0
!    Pts1(3,avgpts+18) = xmu3 + 0.0*sigy
!    Pts1(4,avgpts+18) = xmu4 + 0.0
!    Pts1(5,avgpts+18) = 0.0
!    Pts1(6,avgpts+18) = 0.0
!    Pts1(1,avgpts+19) = xmu1 - 3.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+19) = xmu2 + 0.0
!    Pts1(3,avgpts+19) = xmu3 + 3.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+19) = xmu4 + 0.0
!    Pts1(5,avgpts+19) = 0.0
!    Pts1(6,avgpts+19) = 0.0
!    Pts1(1,avgpts+20) = xmu1 - 3.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+20) = xmu2 + 0.0
!    Pts1(3,avgpts+20) = xmu3 - 3.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+20) = xmu4 + 0.0
!    Pts1(5,avgpts+20) = 0.0
!    Pts1(6,avgpts+20) = 0.0
!    Pts1(1,avgpts+21) = xmu1 + 3.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+21) = xmu2 + 0.0
!    Pts1(3,avgpts+21) = xmu3 + 3.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+21) = xmu4 + 0.0
!    Pts1(5,avgpts+21) = 0.0
!    Pts1(6,avgpts+21) = 0.0
!    Pts1(1,avgpts+22) = xmu1 + 3.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+22) = xmu2 + 0.0
!    Pts1(3,avgpts+22) = xmu3 - 3.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+22) = xmu4 + 0.0
!    Pts1(5,avgpts+22) = 0.0
!    Pts1(6,avgpts+22) = 0.0
!
!    Pts1(1,avgpts+23) = xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+23) = xmu2 + 0.0
!    Pts1(3,avgpts+23) = xmu3 + 2.0*sigy
!    Pts1(4,avgpts+23) = xmu4 + 0.0
!    Pts1(5,avgpts+23) = 0.0
!    Pts1(6,avgpts+23) = 0.0
!    Pts1(1,avgpts+24) = xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+24) = xmu2 + 0.0
!    Pts1(3,avgpts+24) = xmu3 - 2.0*sigy
!    Pts1(4,avgpts+24) = xmu4 + 0.0
!    Pts1(5,avgpts+24) = 0.0
!    Pts1(6,avgpts+24) = 0.0
!    Pts1(1,avgpts+25) = xmu1 + 2.0*sigx 
!    Pts1(2,avgpts+25) = xmu2 + 0.0
!    Pts1(3,avgpts+25) = xmu3 - 0.0*sigy
!    Pts1(4,avgpts+25) = xmu4 + 0.0
!    Pts1(5,avgpts+25) = 0.0
!    Pts1(6,avgpts+25) = 0.0
!    Pts1(1,avgpts+26) = xmu1 - 2.0*sigx 
!    Pts1(2,avgpts+26) = xmu2 + 0.0
!    Pts1(3,avgpts+26) = xmu3 + 0.0*sigy
!    Pts1(4,avgpts+26) = xmu4 + 0.0
!    Pts1(5,avgpts+26) = 0.0
!    Pts1(6,avgpts+26) = 0.0
!    Pts1(1,avgpts+27) = xmu1 - 2.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+27) = xmu2 + 0.0
!    Pts1(3,avgpts+27) = xmu3 + 2.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+27) = xmu4 + 0.0
!    Pts1(5,avgpts+27) = 0.0
!    Pts1(6,avgpts+27) = 0.0
!    Pts1(1,avgpts+28) = xmu1 - 2.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+28) = xmu2 + 0.0
!    Pts1(3,avgpts+28) = xmu3 - 2.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+28) = xmu4 + 0.0
!    Pts1(5,avgpts+28) = 0.0
!    Pts1(6,avgpts+28) = 0.0
!    Pts1(1,avgpts+29) = xmu1 + 2.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+29) = xmu2 + 0.0
!    Pts1(3,avgpts+29) = xmu3 + 2.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+29) = xmu4 + 0.0
!    Pts1(5,avgpts+29) = 0.0
!    Pts1(6,avgpts+29) = 0.0
!    Pts1(1,avgpts+30) = xmu1 + 2.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+30) = xmu2 + 0.0
!    Pts1(3,avgpts+30) = xmu3 - 2.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+30) = xmu4 + 0.0
!    Pts1(5,avgpts+30) = 0.0
!    Pts1(6,avgpts+30) = 0.0

!    Pts1(1,avgpts+31) = xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+31) = xmu2 + 0.0
!    Pts1(3,avgpts+31) = xmu3 + 1.5*sigy
!    Pts1(4,avgpts+31) = xmu4 + 0.0
!!!    Pts1(5,avgpts+31) = 0.0
!    Pts1(6,avgpts+31) = 0.0
!    Pts1(1,avgpts+32) = xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+32) = xmu2 + 0.0
!    Pts1(3,avgpts+32) = xmu3 - 1.5*sigy
!    Pts1(4,avgpts+32) = xmu4 + 0.0
!    Pts1(5,avgpts+32) = 0.0
!    Pts1(6,avgpts+32) = 0.0
!    Pts1(1,avgpts+33) = xmu1 + 1.5*sigx 
!    Pts1(2,avgpts+33) = xmu2 + 0.0
!    Pts1(3,avgpts+33) = xmu3 - 0.0*sigy
!    Pts1(4,avgpts+33) = xmu4 + 0.0
!    Pts1(5,avgpts+33) = 0.0
!    Pts1(6,avgpts+33) = 0.0
!    Pts1(1,avgpts+34) = xmu1 - 1.5*sigx 
!    Pts1(2,avgpts+34) = xmu2 + 0.0
!    Pts1(3,avgpts+34) = xmu3 + 0.0*sigy
!    Pts1(4,avgpts+34) = xmu4 + 0.0
!    Pts1(5,avgpts+34) = 0.0
!    Pts1(6,avgpts+34) = 0.0
!    Pts1(1,avgpts+35) = xmu1 - 1.5*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+35) = xmu2 + 0.0
!    Pts1(3,avgpts+35) = xmu3 + 1.5*sigy/sqrt(2.0)
!    Pts1(4,avgpts+35) = xmu4 + 0.0
!    Pts1(5,avgpts+35) = 0.0
!    Pts1(6,avgpts+35) = 0.0
!    Pts1(1,avgpts+36) = xmu1 - 1.5*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+36) = xmu2 + 0.0
!    Pts1(3,avgpts+36) = xmu3 - 1.5*sigy/sqrt(2.0)
!    Pts1(4,avgpts+36) = xmu4 + 0.0
!    Pts1(5,avgpts+36) = 0.0
!    Pts1(6,avgpts+36) = 0.0
!    Pts1(1,avgpts+37) = xmu1 + 1.5*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+37) = xmu2 + 0.0
!    Pts1(3,avgpts+37) = xmu3 + 1.5*sigy/sqrt(2.0)
!    Pts1(4,avgpts+37) = xmu4 + 0.0
!    Pts1(5,avgpts+37) = 0.0
!    Pts1(6,avgpts+37) = 0.0
!    Pts1(1,avgpts+38) = xmu1 + 1.5*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+38) = xmu2 + 0.0
!    Pts1(3,avgpts+38) = xmu3 - 1.5*sigy/sqrt(2.0)
!    Pts1(4,avgpts+38) = xmu4 + 0.0
!    Pts1(5,avgpts+38) = 0.0
!    Pts1(6,avgpts+38) = 0.0
!
!    Pts1(1,avgpts+39) = xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+39) = xmu2 + 0.0
!    Pts1(3,avgpts+39) = xmu3 + 1.0*sigy
!    Pts1(4,avgpts+39) = xmu4 + 0.0
!    Pts1(5,avgpts+39) = 0.0
!    Pts1(6,avgpts+39) = 0.0
!    Pts1(1,avgpts+40) = xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+40) = xmu2 + 0.0
!    Pts1(3,avgpts+40) = xmu3 - 1.0*sigy
!    Pts1(4,avgpts+40) = xmu4 + 0.0
!    Pts1(5,avgpts+40) = 0.0
!    Pts1(6,avgpts+40) = 0.0
!    Pts1(1,avgpts+41) = xmu1 + 1.0*sigx 
!    Pts1(2,avgpts+41) = xmu2 + 0.0
!    Pts1(3,avgpts+41) = xmu3 - 0.0*sigy
!    Pts1(4,avgpts+41) = xmu4 + 0.0
!    Pts1(5,avgpts+41) = 0.0
!    Pts1(6,avgpts+41) = 0.0
!    Pts1(1,avgpts+42) = xmu1 - 1.0*sigx 
!    Pts1(2,avgpts+42) = xmu2 + 0.0
!    Pts1(3,avgpts+42) = xmu3 + 0.0*sigy
!!    Pts1(4,avgpts+42) = xmu4 + 0.0
!    Pts1(5,avgpts+42) = 0.0
!    Pts1(6,avgpts+42) = 0.0
!    Pts1(1,avgpts+43) = xmu1 - 1.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+43) = xmu2 + 0.0
!    Pts1(3,avgpts+43) = xmu3 + 1.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+43) = xmu4 + 0.0
!    Pts1(5,avgpts+43) = 0.0
!    Pts1(6,avgpts+43) = 0.0
!    Pts1(1,avgpts+44) = xmu1 - 1.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+44) = xmu2 + 0.0
!    Pts1(3,avgpts+44) = xmu3 - 1.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+44) = xmu4 + 0.0
!    Pts1(5,avgpts+44) = 0.0
!    Pts1(6,avgpts+44) = 0.0
!    Pts1(1,avgpts+45) = xmu1 + 1.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+45) = xmu2 + 0.0
!    Pts1(3,avgpts+45) = xmu3 + 1.0*sigy/sqrt(2.0)
!!    Pts1(4,avgpts+45) = xmu4 + 0.0
!    Pts1(5,avgpts+45) = 0.0
!    Pts1(6,avgpts+45) = 0.0
!    Pts1(1,avgpts+46) = xmu1 + 1.0*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+46) = xmu2 + 0.0
!    Pts1(3,avgpts+46) = xmu3 - 1.0*sigy/sqrt(2.0)
!    Pts1(4,avgpts+46) = xmu4 + 0.0
!    Pts1(5,avgpts+46) = 0.0
!    Pts1(6,avgpts+46) = 0.0

!    Pts1(1,avgpts+47) = xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+47) = xmu2 + 0.0
!    Pts1(3,avgpts+47) = xmu3 + 0.5*sigy
!    Pts1(4,avgpts+47) = xmu4 + 0.0
!!    Pts1(5,avgpts+47) = 0.0
!    Pts1(6,avgpts+47) = 0.0
!    Pts1(1,avgpts+48) = xmu1 + 0.0*sigx 
!    Pts1(2,avgpts+48) = xmu2 + 0.0
!    Pts1(3,avgpts+48) = xmu3 - 0.5*sigy
!    Pts1(4,avgpts+48) = xmu4 + 0.0
!    Pts1(5,avgpts+48) = 0.0
!    Pts1(6,avgpts+48) = 0.0
!    Pts1(1,avgpts+49) = xmu1 + 0.5*sigx 
!    Pts1(2,avgpts+49) = xmu2 + 0.0
!    Pts1(3,avgpts+49) = xmu3 - 0.0*sigy
!    Pts1(4,avgpts+49) = xmu4 + 0.0
!    Pts1(5,avgpts+49) = 0.0
!    Pts1(6,avgpts+49) = 0.0
!    Pts1(1,avgpts+50) = xmu1 - 0.5*sigx 
!    Pts1(2,avgpts+50) = xmu2 + 0.0
!    Pts1(3,avgpts+50) = xmu3 + 0.0*sigy
!    Pts1(4,avgpts+50) = xmu4 + 0.0
!    Pts1(5,avgpts+50) = 0.0
!    Pts1(6,avgpts+50) = 0.0
!    Pts1(1,avgpts+51) = xmu1 - 0.5*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+51) = xmu2 + 0.0
!    Pts1(3,avgpts+51) = xmu3 + 0.5*sigy/sqrt(2.0)
!    Pts1(4,avgpts+51) = xmu4 + 0.0
!    Pts1(5,avgpts+51) = 0.0
!    Pts1(6,avgpts+51) = 0.0
!!    Pts1(1,avgpts+52) = xmu1 - 0.5*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+52) = xmu2 + 0.0
!    Pts1(3,avgpts+52) = xmu3 - 0.5*sigy/sqrt(2.0)
!    Pts1(4,avgpts+52) = xmu4 + 0.0
!    Pts1(5,avgpts+52) = 0.0
!    Pts1(6,avgpts+52) = 0.0
!    Pts1(1,avgpts+53) = xmu1 + 0.5*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+53) = xmu2 + 0.0
!    Pts1(3,avgpts+53) = xmu3 + 0.5*sigy/sqrt(2.0)
!    Pts1(4,avgpts+53) = xmu4 + 0.0
!    Pts1(5,avgpts+53) = 0.0
!    Pts1(6,avgpts+53) = 0.0
!    Pts1(1,avgpts+54) = xmu1 + 0.5*sigx/sqrt(2.0) 
!    Pts1(2,avgpts+54) = xmu2 + 0.0
!    Pts1(3,avgpts+54) = xmu3 - 0.5*sigy/sqrt(2.0)
!    Pts1(4,avgpts+54) = xmu4 + 0.0
!    Pts1(5,avgpts+54) = 0.0
!    Pts1(6,avgpts+54) = 0.0

1000 continue

    deallocate(x1)
    deallocate(x2)
    deallocate(x3)

    !Nptlocal = avgpts

    t_kvdist = t_kvdist + elapsedtime_Timer(t0)

  end subroutine Gauss3_Dist



  !Note: particle transfer is done from IP6 to LR-BB.
  subroutine Gauss3Sob_Dist(Pts1,nparam,distparam,grid,flagalloc,&
       Npt,Nptlocal,ngroup,iseed)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: nparam,flagalloc,Npt,Nptlocal,ngroup
    double precision, dimension(nparam) :: distparam
    type (Pgrid2d), intent(in) :: grid
    integer, intent(inout) :: iseed
    double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
         sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
    double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
    double precision :: sig1,sig2,sig3,sig4,sig5,sig6
    double precision :: sq12,sq34,sq56
    double precision, allocatable, dimension(:,:) :: x1,x2,x3 
    integer :: totnp,npy,npx
    integer :: avgpts
    integer :: myid,myidx,myidy,i,j,k,intvsamp
    !        integer seedarray(1)
    double precision :: t0,x11
    double precision :: xrms,pxrms,yrms,pyrms,xavg,x2avg,pxavg,px2avg,&
         yavg,y2avg,pyavg,py2avg
    double precision, dimension(16) :: tmplc,tmpgl
    integer :: npyhalf,ierr,ist,ilow,ihigh
    real*8 :: val1,val2,val3,val4,val5,val6,tmpval1,tmpval2,tmpval3,&
         tmpval4,tmpval5,tmpval6
    real*8, dimension(6,6) :: tmp

    call starttime_Timer(t0)

    sigx = distparam(1)
    sigpx = distparam(2)
    muxpx = distparam(3)
    xscale = distparam(4)
    pxscale = distparam(5)
    xmu1 = distparam(6)
    xmu2 = distparam(7)
    sigy = distparam(8)
    sigpy = distparam(9)
    muypy = distparam(10)
    yscale = distparam(11)
    pyscale = distparam(12)
    xmu3 = distparam(13)
    xmu4 = distparam(14)
    sigz = distparam(15)
    sigpz = distparam(16)
    muzpz = distparam(17)
    zscale = distparam(18)
    pzscale = distparam(19)
    xmu5 = distparam(20)
    xmu6 = distparam(21)
    !for LHC offset study only
    !no offset because offset is done in AccSimulator.f90
    !xmu1 = 0.0
    !xmu2 = 0.0
    !xmu3 = 0.0
    !xmu4 = 0.0
    !xmu5 = 0.0
    !xmu6 = 0.0

    call getsize_Pgrid2d(grid,totnp,npy,npx)

    call getpost_Pgrid2d(grid,myid,myidy,myidx)
    !        seedarray(1)=(1001+myid)*(myid+7)
    !        write(6,*)'seedarray=',seedarray
    !        call random_seed(PUT=seedarray(1:1))
    call random_number(x11)

    avgpts = Npt/(npx*npy/ngroup)

    if(myidy.lt.npy/2) then
       ilow = myid*avgpts
       ihigh = (myid+1)*avgpts
       ist = 0
       !          ist = 303
    else
       ilow = (myid-totnp/2)*avgpts
       ihigh = (myid-totnp/2+1)*avgpts
       ist = 3030
       !ist = 0 !make two beams the same initial distribution.
    endif

    sig1 = sigx*xscale
    sig2 = sigpx*pxscale
    sig3 = sigy*yscale
    sig4 = sigpy*pyscale
    sig5 = sigz*zscale
    sig6 = sigpz*pzscale

    sq12=sqrt(1.-muxpx*muxpx)
    sq34=sqrt(1.-muypy*muypy)
    sq56=sqrt(1.-muzpz*muzpz)

    ! initial allocate 'avgpts' particles on each processor.
    !if(flagalloc.eq.1) then
    !  Pts1 = 0.0
    !else
    !  allocate(Pts1(6,avgpts))
    !  Pts1 = 0.0
    !endif

    !        allocate(x1(2,avgpts))
    !        allocate(x2(2,avgpts))
    !        allocate(x3(2,avgpts))
    !        call normVec(x1,avgpts)
    !        call normVec(x2,avgpts)
    !        call normVec(x3,avgpts)

    !intvsamp = 10
    intvsamp = avgpts
    allocate(x1(2,Npt))
    allocate(x2(2,Npt))
    allocate(x3(2,Npt))

    xavg = 0.0
    x2avg = 0.0
    pxavg = 0.0
    px2avg = 0.0
    yavg = 0.0
    y2avg = 0.0
    pyavg = 0.0
    py2avg = 0.0
    !        call MPI_BARRIER(mpicommwd,ierr)
    do j = 1, 1
       !          call normVecnew(iseed,x1,intvsamp)
       !          call normVecnew(iseed,x2,intvsamp)
       !          call normVecnew(iseed,x3,intvsamp)
       call normVecSob(x1,x2,x3,Npt,ist)
       do k = 1, Npt
          !x-px:
          !            call normdv(x1)
          !           Correct Gaussian distribution.
          if(k.gt.ilow .and. k.le.ihigh) then
             i = k-ilow
             !Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
             !Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
             !miguel's sampling. here, muxpx->alphax,sigpx->betax
             Pts1(1,i) = xmu1 + sigx*x1(1,k)*xscale
             xavg = xavg + Pts1(1,i)
             x2avg = x2avg + Pts1(1,i)*Pts1(1,i)
             Pts1(2,i) = xmu2 + sigx  &
                  *(-muxpx*x1(1,k)+x1(2,k))/sigpx*pxscale
             pxavg = pxavg + Pts1(2,i)
             px2avg = px2avg + Pts1(2,i)*Pts1(2,i)
             !y-py
             !              call normdv(x1)
             !             Correct Gaussian distribution.
             !Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
             !Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))
             !miguel's sampling. here, muxpx->alphax,sigpx->betax
             Pts1(3,i) = xmu3 + sigy*x2(1,k)*yscale
             yavg = yavg + Pts1(3,i)
             y2avg = y2avg + Pts1(3,i)*Pts1(3,i)
             Pts1(4,i) = xmu4 + sigy  &
                  *(-muypy*x2(1,k)+x2(2,k))/sigpy*pyscale
             pyavg = pyavg + Pts1(4,i)
             py2avg = py2avg + Pts1(4,i)*Pts1(4,i)
             !z-pz
             !              call normdv(x1)
             !             Correct Gaussian distribution.
             !Pts1(5,i) = xmu5 + sig5*x3(1,k)/sq56
             !Pts1(6,i) = xmu6 + sig6*(-muzpz*x3(1,k)/sq56+x3(2,k))
             !miguel's sampling. here, muxpx->alphax,sigpx->betax
             Pts1(5,i) = xmu5 + x3(1,k)*sig5
             Pts1(6,i) = xmu6 + x3(2,k)*sig6
          endif
       enddo
    enddo

    npyhalf = npy/ngroup
    tmplc = 0.0d0
    if(myidy.lt.npyhalf) then
       tmplc(1) = xavg
       tmplc(2) = x2avg
       tmplc(3) = pxavg
       tmplc(4) = px2avg
       tmplc(5) = yavg
       tmplc(6) = y2avg
       tmplc(7) = pyavg
       tmplc(8) = py2avg
    else
       tmplc(9) = xavg
       tmplc(10) = x2avg
       tmplc(11) = pxavg
       tmplc(12) = px2avg
       tmplc(13) = yavg
       tmplc(14) = y2avg
       tmplc(15) = pyavg
       tmplc(16) = py2avg
    endif

    !        call MPI_BARRIER(mpicommwd,ierr)

    call MPI_ALLREDUCE(tmplc,tmpgl,16,MPI_DOUBLE_PRECISION,&
         MPI_SUM,mpicommwd,ierr)
    do i = 1, 16
       tmpgl(i) = tmpgl(i)/Npt
    enddo

    if(myidy.lt.npyhalf) then
       xrms = sqrt(tmpgl(2)-tmpgl(1)**2)
       pxrms = sqrt(tmpgl(4)-tmpgl(3)**2)
       yrms = sqrt(tmpgl(6)-tmpgl(5)**2)
       pyrms = sqrt(tmpgl(8)-tmpgl(7)**2)
    else
       xrms = sqrt(tmpgl(10)-tmpgl(9)**2)
       pxrms = sqrt(tmpgl(12)-tmpgl(11)**2)
       yrms = sqrt(tmpgl(14)-tmpgl(13)**2)
       pyrms = sqrt(tmpgl(16)-tmpgl(15)**2)
    endif

    do j = 1, avgpts
       Pts1(1,j) = Pts1(1,j)*sigx/xrms
       Pts1(2,j) = Pts1(2,j)*(sigx/sigpx)*sqrt(1+muxpx**2)/pxrms
       Pts1(3,j) = Pts1(3,j)*sigy/yrms
       Pts1(4,j) = Pts1(4,j)*(sigy/sigpy)*sqrt(1+muypy**2)/pyrms
    enddo
    !        Pts1(1,avgpts) = 0.0
    !        Pts1(2,avgpts) = 0.0
    !        Pts1(3,avgpts) = xmu3 + 5.0*sigy
    !        Pts1(4,avgpts) = 0.0
    !        Pts1(5,avgpts) = 0.0
    !        Pts1(6,avgpts) = 0.0

    goto 1000

    !sigx = 0.0
    !sigy = 0.0

    Pts1(1,avgpts+1)= xmu1 + 0.0*sigx 
    Pts1(2,avgpts+1)= xmu2 + 0.0
    Pts1(3,avgpts+1) = xmu3 + 5.0*sigy
    Pts1(4,avgpts+1) = xmu4 + 0.0
    Pts1(5,avgpts+1) = 0.0
    Pts1(6,avgpts+1) = 0.0
    Pts1(1,avgpts+2) = xmu1 + 0.0*sigx 
    Pts1(2,avgpts+2) = xmu2 + 0.0
    Pts1(3,avgpts+2) = xmu3 - 5.0*sigy
    Pts1(4,avgpts+2) = xmu4 + 0.0
    Pts1(5,avgpts+2) = 0.0
    Pts1(6,avgpts+2) = 0.0
    Pts1(1,avgpts+3) = xmu1 + 5.0*sigx 
    Pts1(2,avgpts+3) = xmu2 + 0.0
    Pts1(3,avgpts+3) = xmu3 - 0.0*sigy
    Pts1(4,avgpts+3) = xmu4 + 0.0
    Pts1(5,avgpts+3) = 0.0
    Pts1(6,avgpts+3) = 0.0
    Pts1(1,avgpts+4) = xmu1 - 5.0*sigx 
    Pts1(2,avgpts+4) = xmu2 + 0.0
    Pts1(3,avgpts+4) = xmu3 + 0.0*sigy
    Pts1(4,avgpts+4) = xmu4 + 0.0
    Pts1(5,avgpts+4) = 0.0
    Pts1(6,avgpts+4) = 0.0
    Pts1(1,avgpts+5) = xmu1 + 5.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+5) = xmu2 + 0.0
    Pts1(3,avgpts+5) = xmu3 + 5.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+5) = xmu4 + 0.0
    Pts1(5,avgpts+5) = 0.0
    Pts1(6,avgpts+5) = 0.0
    Pts1(1,avgpts+6) = xmu1 - 5.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+6) = xmu2 + 0.0
    Pts1(3,avgpts+6) = xmu3 - 5.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+6) = xmu4 + 0.0
    Pts1(5,avgpts+6) = 0.0
    Pts1(6,avgpts+6) = 0.0

    Pts1(1,avgpts+7)= xmu1 + 0.0*sigx 
    Pts1(2,avgpts+7)= xmu2 + 0.0
    Pts1(3,avgpts+7) = xmu3 + 4.0*sigy
    Pts1(4,avgpts+7) = xmu4 + 0.0
    Pts1(5,avgpts+7) = 0.0
    Pts1(6,avgpts+7) = 0.0
    Pts1(1,avgpts+8) = xmu1 + 0.0*sigx 
    Pts1(2,avgpts+8) = xmu2 + 0.0
    Pts1(3,avgpts+8) = xmu3 - 4.0*sigy
    Pts1(4,avgpts+8) = xmu4 + 0.0
    Pts1(5,avgpts+8) = 0.0
    Pts1(6,avgpts+8) = 0.0
    Pts1(1,avgpts+9) = xmu1 + 4.0*sigx 
    Pts1(2,avgpts+9) = xmu2 + 0.0
    Pts1(3,avgpts+9) = xmu3 - 0.0*sigy
    Pts1(4,avgpts+9) = xmu4 + 0.0
    Pts1(5,avgpts+9) = 0.0
    Pts1(6,avgpts+9) = 0.0
    Pts1(1,avgpts+10) = xmu1 - 4.0*sigx 
    Pts1(2,avgpts+10) = xmu2 + 0.0
    Pts1(3,avgpts+10) = xmu3 + 0.0*sigy
    Pts1(4,avgpts+10) = xmu4 + 0.0
    Pts1(5,avgpts+10) = 0.0
    Pts1(6,avgpts+10) = 0.0
    Pts1(1,avgpts+11) = xmu1 - 4.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+11) = xmu2 + 0.0
    Pts1(3,avgpts+11) = xmu3 + 4.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+11) = xmu4 + 0.0
    Pts1(5,avgpts+11) = 0.0
    Pts1(6,avgpts+11) = 0.0
    Pts1(1,avgpts+12) = xmu1 - 4.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+12) = xmu2 + 0.0
    Pts1(3,avgpts+12) = xmu3 - 4.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+12) = xmu4 + 0.0
    Pts1(5,avgpts+12) = 0.0
    Pts1(6,avgpts+12) = 0.0
    Pts1(1,avgpts+13) = xmu1 + 4.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+13) = xmu2 + 0.0
    Pts1(3,avgpts+13) = xmu3 + 4.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+13) = xmu4 + 0.0
    Pts1(5,avgpts+13) = 0.0
    Pts1(6,avgpts+13) = 0.0
    Pts1(1,avgpts+14) = xmu1 + 4.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+14) = xmu2 + 0.0
    Pts1(3,avgpts+14) = xmu3 - 4.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+14) = xmu4 + 0.0
    Pts1(5,avgpts+14) = 0.0
    Pts1(6,avgpts+14) = 0.0

    Pts1(1,avgpts+15)= xmu1 + 0.0*sigx 
    Pts1(2,avgpts+15)= xmu2 + 0.0
    Pts1(3,avgpts+15) = xmu3 + 3.0*sigy
    Pts1(4,avgpts+15) = xmu4 + 0.0
    Pts1(5,avgpts+15) = 0.0
    Pts1(6,avgpts+15) = 0.0
    Pts1(1,avgpts+16) = xmu1 + 0.0*sigx 
    Pts1(2,avgpts+16) = xmu2 + 0.0
    Pts1(3,avgpts+16) = xmu3 - 3.0*sigy
    Pts1(4,avgpts+16) = xmu4 + 0.0
    Pts1(5,avgpts+16) = 0.0
    Pts1(6,avgpts+16) = 0.0
    Pts1(1,avgpts+17) = xmu1 + 3.0*sigx 
    Pts1(2,avgpts+17) = xmu2 + 0.0
    Pts1(3,avgpts+17) = xmu3 - 0.0*sigy
    Pts1(4,avgpts+17) = xmu4 + 0.0
    Pts1(5,avgpts+17) = 0.0
    Pts1(6,avgpts+17) = 0.0
    Pts1(1,avgpts+18) = xmu1 - 3.0*sigx 
    Pts1(2,avgpts+18) = xmu2 + 0.0
    Pts1(3,avgpts+18) = xmu3 + 0.0*sigy
    Pts1(4,avgpts+18) = xmu4 + 0.0
    Pts1(5,avgpts+18) = 0.0
    Pts1(6,avgpts+18) = 0.0
    Pts1(1,avgpts+19) = xmu1 - 3.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+19) = xmu2 + 0.0
    Pts1(3,avgpts+19) = xmu3 + 3.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+19) = xmu4 + 0.0
    Pts1(5,avgpts+19) = 0.0
    Pts1(6,avgpts+19) = 0.0
    Pts1(1,avgpts+20) = xmu1 - 3.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+20) = xmu2 + 0.0
    Pts1(3,avgpts+20) = xmu3 - 3.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+20) = xmu4 + 0.0
    Pts1(5,avgpts+20) = 0.0
    Pts1(6,avgpts+20) = 0.0
    Pts1(1,avgpts+21) = xmu1 + 3.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+21) = xmu2 + 0.0
    Pts1(3,avgpts+21) = xmu3 + 3.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+21) = xmu4 + 0.0
    Pts1(5,avgpts+21) = 0.0
    Pts1(6,avgpts+21) = 0.0
    Pts1(1,avgpts+22) = xmu1 + 3.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+22) = xmu2 + 0.0
    Pts1(3,avgpts+22) = xmu3 - 3.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+22) = xmu4 + 0.0
    Pts1(5,avgpts+22) = 0.0
    Pts1(6,avgpts+22) = 0.0

    Pts1(1,avgpts+23) = xmu1 + 0.0*sigx 
    Pts1(2,avgpts+23) = xmu2 + 0.0
    Pts1(3,avgpts+23) = xmu3 + 2.0*sigy
    Pts1(4,avgpts+23) = xmu4 + 0.0
    Pts1(5,avgpts+23) = 0.0
    Pts1(6,avgpts+23) = 0.0
    Pts1(1,avgpts+24) = xmu1 + 0.0*sigx 
    Pts1(2,avgpts+24) = xmu2 + 0.0
    Pts1(3,avgpts+24) = xmu3 - 2.0*sigy
    Pts1(4,avgpts+24) = xmu4 + 0.0
    Pts1(5,avgpts+24) = 0.0
    Pts1(6,avgpts+24) = 0.0
    Pts1(1,avgpts+25) = xmu1 + 2.0*sigx 
    Pts1(2,avgpts+25) = xmu2 + 0.0
    Pts1(3,avgpts+25) = xmu3 - 0.0*sigy
    Pts1(4,avgpts+25) = xmu4 + 0.0
    Pts1(5,avgpts+25) = 0.0
    Pts1(6,avgpts+25) = 0.0
    Pts1(1,avgpts+26) = xmu1 - 2.0*sigx 
    Pts1(2,avgpts+26) = xmu2 + 0.0
    Pts1(3,avgpts+26) = xmu3 + 0.0*sigy
    Pts1(4,avgpts+26) = xmu4 + 0.0
    Pts1(5,avgpts+26) = 0.0
    Pts1(6,avgpts+26) = 0.0
    Pts1(1,avgpts+27) = xmu1 - 2.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+27) = xmu2 + 0.0
    Pts1(3,avgpts+27) = xmu3 + 2.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+27) = xmu4 + 0.0
    Pts1(5,avgpts+27) = 0.0
    Pts1(6,avgpts+27) = 0.0
    Pts1(1,avgpts+28) = xmu1 - 2.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+28) = xmu2 + 0.0
    Pts1(3,avgpts+28) = xmu3 - 2.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+28) = xmu4 + 0.0
    Pts1(5,avgpts+28) = 0.0
    Pts1(6,avgpts+28) = 0.0
    Pts1(1,avgpts+29) = xmu1 + 2.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+29) = xmu2 + 0.0
    Pts1(3,avgpts+29) = xmu3 + 2.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+29) = xmu4 + 0.0
    Pts1(5,avgpts+29) = 0.0
    Pts1(6,avgpts+29) = 0.0
    Pts1(1,avgpts+30) = xmu1 + 2.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+30) = xmu2 + 0.0
    Pts1(3,avgpts+30) = xmu3 - 2.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+30) = xmu4 + 0.0
    Pts1(5,avgpts+30) = 0.0
    Pts1(6,avgpts+30) = 0.0

    Pts1(1,avgpts+31) = xmu1 + 0.0*sigx 
    Pts1(2,avgpts+31) = xmu2 + 0.0
    Pts1(3,avgpts+31) = xmu3 + 1.5*sigy
    Pts1(4,avgpts+31) = xmu4 + 0.0
    Pts1(5,avgpts+31) = 0.0
    Pts1(6,avgpts+31) = 0.0
    Pts1(1,avgpts+32) = xmu1 + 0.0*sigx 
    Pts1(2,avgpts+32) = xmu2 + 0.0
    Pts1(3,avgpts+32) = xmu3 - 1.5*sigy
    Pts1(4,avgpts+32) = xmu4 + 0.0
    Pts1(5,avgpts+32) = 0.0
    Pts1(6,avgpts+32) = 0.0
    Pts1(1,avgpts+33) = xmu1 + 1.5*sigx 
    Pts1(2,avgpts+33) = xmu2 + 0.0
    Pts1(3,avgpts+33) = xmu3 - 0.0*sigy
    Pts1(4,avgpts+33) = xmu4 + 0.0
    Pts1(5,avgpts+33) = 0.0
    Pts1(6,avgpts+33) = 0.0
    Pts1(1,avgpts+34) = xmu1 - 1.5*sigx 
    Pts1(2,avgpts+34) = xmu2 + 0.0
    Pts1(3,avgpts+34) = xmu3 + 0.0*sigy
    Pts1(4,avgpts+34) = xmu4 + 0.0
    Pts1(5,avgpts+34) = 0.0
    Pts1(6,avgpts+34) = 0.0
    Pts1(1,avgpts+35) = xmu1 - 1.5*sigx/sqrt(2.0) 
    Pts1(2,avgpts+35) = xmu2 + 0.0
    Pts1(3,avgpts+35) = xmu3 + 1.5*sigy/sqrt(2.0)
    Pts1(4,avgpts+35) = xmu4 + 0.0
    Pts1(5,avgpts+35) = 0.0
    Pts1(6,avgpts+35) = 0.0
    Pts1(1,avgpts+36) = xmu1 - 1.5*sigx/sqrt(2.0) 
    Pts1(2,avgpts+36) = xmu2 + 0.0
    Pts1(3,avgpts+36) = xmu3 - 1.5*sigy/sqrt(2.0)
    Pts1(4,avgpts+36) = xmu4 + 0.0
    Pts1(5,avgpts+36) = 0.0
    Pts1(6,avgpts+36) = 0.0
    Pts1(1,avgpts+37) = xmu1 + 1.5*sigx/sqrt(2.0) 
    Pts1(2,avgpts+37) = xmu2 + 0.0
    Pts1(3,avgpts+37) = xmu3 + 1.5*sigy/sqrt(2.0)
    Pts1(4,avgpts+37) = xmu4 + 0.0
    Pts1(5,avgpts+37) = 0.0
    Pts1(6,avgpts+37) = 0.0
    Pts1(1,avgpts+38) = xmu1 + 1.5*sigx/sqrt(2.0) 
    Pts1(2,avgpts+38) = xmu2 + 0.0
    Pts1(3,avgpts+38) = xmu3 - 1.5*sigy/sqrt(2.0)
    Pts1(4,avgpts+38) = xmu4 + 0.0
    Pts1(5,avgpts+38) = 0.0
    Pts1(6,avgpts+38) = 0.0

    Pts1(1,avgpts+39) = xmu1 + 0.0*sigx 
    Pts1(2,avgpts+39) = xmu2 + 0.0
    Pts1(3,avgpts+39) = xmu3 + 1.0*sigy
    Pts1(4,avgpts+39) = xmu4 + 0.0
    Pts1(5,avgpts+39) = 0.0
    Pts1(6,avgpts+39) = 0.0
    Pts1(1,avgpts+40) = xmu1 + 0.0*sigx 
    Pts1(2,avgpts+40) = xmu2 + 0.0
    Pts1(3,avgpts+40) = xmu3 - 1.0*sigy
    Pts1(4,avgpts+40) = xmu4 + 0.0
    Pts1(5,avgpts+40) = 0.0
    Pts1(6,avgpts+40) = 0.0
    Pts1(1,avgpts+41) = xmu1 + 1.0*sigx 
    Pts1(2,avgpts+41) = xmu2 + 0.0
    Pts1(3,avgpts+41) = xmu3 - 0.0*sigy
    Pts1(4,avgpts+41) = xmu4 + 0.0
    Pts1(5,avgpts+41) = 0.0
    Pts1(6,avgpts+41) = 0.0
    Pts1(1,avgpts+42) = xmu1 - 1.0*sigx 
    Pts1(2,avgpts+42) = xmu2 + 0.0
    Pts1(3,avgpts+42) = xmu3 + 0.0*sigy
    Pts1(4,avgpts+42) = xmu4 + 0.0
    Pts1(5,avgpts+42) = 0.0
    Pts1(6,avgpts+42) = 0.0
    Pts1(1,avgpts+43) = xmu1 - 1.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+43) = xmu2 + 0.0
    Pts1(3,avgpts+43) = xmu3 + 1.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+43) = xmu4 + 0.0
    Pts1(5,avgpts+43) = 0.0
    Pts1(6,avgpts+43) = 0.0
    Pts1(1,avgpts+44) = xmu1 - 1.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+44) = xmu2 + 0.0
    Pts1(3,avgpts+44) = xmu3 - 1.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+44) = xmu4 + 0.0
    Pts1(5,avgpts+44) = 0.0
    Pts1(6,avgpts+44) = 0.0
    Pts1(1,avgpts+45) = xmu1 + 1.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+45) = xmu2 + 0.0
    Pts1(3,avgpts+45) = xmu3 + 1.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+45) = xmu4 + 0.0
    Pts1(5,avgpts+45) = 0.0
    Pts1(6,avgpts+45) = 0.0
    Pts1(1,avgpts+46) = xmu1 + 1.0*sigx/sqrt(2.0) 
    Pts1(2,avgpts+46) = xmu2 + 0.0
    Pts1(3,avgpts+46) = xmu3 - 1.0*sigy/sqrt(2.0)
    Pts1(4,avgpts+46) = xmu4 + 0.0
    Pts1(5,avgpts+46) = 0.0
    Pts1(6,avgpts+46) = 0.0

    Pts1(1,avgpts+47) = xmu1 + 0.0*sigx 
    Pts1(2,avgpts+47) = xmu2 + 0.0
    Pts1(3,avgpts+47) = xmu3 + 0.5*sigy
    Pts1(4,avgpts+47) = xmu4 + 0.0
    Pts1(5,avgpts+47) = 0.0
    Pts1(6,avgpts+47) = 0.0
    Pts1(1,avgpts+48) = xmu1 + 0.0*sigx 
    Pts1(2,avgpts+48) = xmu2 + 0.0
    Pts1(3,avgpts+48) = xmu3 - 0.5*sigy
    Pts1(4,avgpts+48) = xmu4 + 0.0
    Pts1(5,avgpts+48) = 0.0
    Pts1(6,avgpts+48) = 0.0
    Pts1(1,avgpts+49) = xmu1 + 0.5*sigx 
    Pts1(2,avgpts+49) = xmu2 + 0.0
    Pts1(3,avgpts+49) = xmu3 - 0.0*sigy
    Pts1(4,avgpts+49) = xmu4 + 0.0
    Pts1(5,avgpts+49) = 0.0
    Pts1(6,avgpts+49) = 0.0
    Pts1(1,avgpts+50) = xmu1 - 0.5*sigx 
    Pts1(2,avgpts+50) = xmu2 + 0.0
    Pts1(3,avgpts+50) = xmu3 + 0.0*sigy
    Pts1(4,avgpts+50) = xmu4 + 0.0
    Pts1(5,avgpts+50) = 0.0
    Pts1(6,avgpts+50) = 0.0
    Pts1(1,avgpts+51) = xmu1 - 0.5*sigx/sqrt(2.0) 
    Pts1(2,avgpts+51) = xmu2 + 0.0
    Pts1(3,avgpts+51) = xmu3 + 0.5*sigy/sqrt(2.0)
    Pts1(4,avgpts+51) = xmu4 + 0.0
    Pts1(5,avgpts+51) = 0.0
    Pts1(6,avgpts+51) = 0.0
    Pts1(1,avgpts+52) = xmu1 - 0.5*sigx/sqrt(2.0) 
    Pts1(2,avgpts+52) = xmu2 + 0.0
    Pts1(3,avgpts+52) = xmu3 - 0.5*sigy/sqrt(2.0)
    Pts1(4,avgpts+52) = xmu4 + 0.0
    Pts1(5,avgpts+52) = 0.0
    Pts1(6,avgpts+52) = 0.0
    Pts1(1,avgpts+53) = xmu1 + 0.5*sigx/sqrt(2.0) 
    Pts1(2,avgpts+53) = xmu2 + 0.0
    Pts1(3,avgpts+53) = xmu3 + 0.5*sigy/sqrt(2.0)
    Pts1(4,avgpts+53) = xmu4 + 0.0
    Pts1(5,avgpts+53) = 0.0
    Pts1(6,avgpts+53) = 0.0
    Pts1(1,avgpts+54) = xmu1 + 0.5*sigx/sqrt(2.0) 
    Pts1(2,avgpts+54) = xmu2 + 0.0
    Pts1(3,avgpts+54) = xmu3 - 0.5*sigy/sqrt(2.0)
    Pts1(4,avgpts+54) = xmu4 + 0.0
    Pts1(5,avgpts+54) = 0.0
    Pts1(6,avgpts+54) = 0.0

1000 continue

    deallocate(x1)
    deallocate(x2)
    deallocate(x3)

    !Nptlocal = avgpts

    goto 999 !neglecting transfer as a test

    !Nptlocal = avgpts
    !The following is due to the fact that long-range bb is NOT at IP
    !However, we generate our initial distribution at IP.
    if(myidy.lt.npyhalf) then
       open(1,file='map1',status='old')
       do i = 1, 6
          read(1,*)tmp(i,1:6)
       enddo
       close(1)
       !do i = 1, 6
       !  write(11,*)tmp(i,1:6)
       !enddo
       do i = 1, avgpts
          val1 =  Pts1(1,i) - xmu1
          val2 =  Pts1(2,i) - xmu2
          val3 =  Pts1(3,i) - xmu3
          val4 =  Pts1(4,i) - xmu4
          val5 =  Pts1(5,i) - xmu5
          val6 =  Pts1(6,i) - xmu6
          tmpval1 = val1*tmp(1,1)+val2*tmp(1,2)+val3*tmp(1,3)+val4*tmp(1,4)&
               +val5*tmp(1,5)+val6*tmp(1,6)
          tmpval2 = val1*tmp(2,1)+val2*tmp(2,2)+val3*tmp(2,3)+val4*tmp(2,4)&
               +val5*tmp(2,5)+val6*tmp(2,6)
          tmpval3 = val1*tmp(3,1)+val2*tmp(3,2)+val3*tmp(3,3)+val4*tmp(3,4)&
               +val5*tmp(3,5)+val6*tmp(3,6)
          tmpval4 = val1*tmp(4,1)+val2*tmp(4,2)+val3*tmp(4,3)+val4*tmp(4,4)&
               +val5*tmp(4,5)+val6*tmp(4,6)
          tmpval5 = val1*tmp(5,1)+val2*tmp(5,2)+val3*tmp(5,3)+val4*tmp(5,4)&
               +val5*tmp(5,5)+val6*tmp(5,6)
          tmpval6 = val1*tmp(6,1)+val2*tmp(6,2)+val3*tmp(6,3)+val4*tmp(6,4)&
               +val5*tmp(6,5)+val6*tmp(6,6)
          Pts1(1,i) = tmpval1 + xmu1
          Pts1(2,i) = tmpval2 + xmu2
          Pts1(3,i) = tmpval3 + xmu3
          Pts1(4,i) = tmpval4 + xmu4
          Pts1(5,i) = tmpval5 + xmu5
          Pts1(6,i) = tmpval6 + xmu6
       enddo
    else
       open(2,file='map2',status='old')
       do i = 1, 6
          read(2,*)tmp(i,1:6)
       enddo
       close(2)
       !do i = 1, 6
       !  write(12,*)tmp(i,1:6)
       !enddo
       do i = 1, avgpts
          val1 =  Pts1(1,i) - xmu1
          val2 =  Pts1(2,i) - xmu2
          val3 =  Pts1(3,i) - xmu3
          val4 =  Pts1(4,i) - xmu4
          val5 =  Pts1(5,i) - xmu5
          val6 =  Pts1(6,i) - xmu6
          tmpval1 = val1*tmp(1,1)+val2*tmp(1,2)+val3*tmp(1,3)+val4*tmp(1,4)&
               +val5*tmp(1,5)+val6*tmp(1,6)
          tmpval2 = val1*tmp(2,1)+val2*tmp(2,2)+val3*tmp(2,3)+val4*tmp(2,4)&
               +val5*tmp(2,5)+val6*tmp(2,6)
          tmpval3 = val1*tmp(3,1)+val2*tmp(3,2)+val3*tmp(3,3)+val4*tmp(3,4)&
               +val5*tmp(3,5)+val6*tmp(3,6)
          tmpval4 = val1*tmp(4,1)+val2*tmp(4,2)+val3*tmp(4,3)+val4*tmp(4,4)&
               +val5*tmp(4,5)+val6*tmp(4,6)
          tmpval5 = val1*tmp(5,1)+val2*tmp(5,2)+val3*tmp(5,3)+val4*tmp(5,4)&
               +val5*tmp(5,5)+val6*tmp(5,6)
          tmpval6 = val1*tmp(6,1)+val2*tmp(6,2)+val3*tmp(6,3)+val4*tmp(6,4)&
               +val5*tmp(6,5)+val6*tmp(6,6)
          Pts1(1,i) = tmpval1 + xmu1
          Pts1(2,i) = tmpval2 + xmu2
          Pts1(3,i) = tmpval3 + xmu3
          Pts1(4,i) = tmpval4 + xmu4
          Pts1(5,i) = tmpval5 + xmu5
          Pts1(6,i) = tmpval6 + xmu6
       enddo
    endif

999 continue

    t_kvdist = t_kvdist + elapsedtime_Timer(t0)

  end subroutine Gauss3Sob_Dist



  subroutine Gauss3test_Dist(Pts1,nparam,distparam,grid,flagalloc,&
       Npt,Nptlocal,ngroup)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: nparam,flagalloc,Npt,Nptlocal,ngroup
    double precision, dimension(nparam) :: distparam
    type (Pgrid2d), intent(in) :: grid
    double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
         sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
    double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
    double precision :: sig1,sig2,sig3,sig4,sig5,sig6
    double precision :: sq12,sq34,sq56
    double precision, allocatable, dimension(:,:) :: x1,x2,x3 
    integer :: totnp,npy,npx
    integer :: avgpts
    integer :: myid,myidx,myidy,i,j,k,intvsamp,jlow,jhigh
    !        integer seedarray(1)
    double precision :: t0,x11

    call starttime_Timer(t0)

    sigx = distparam(1)
    sigpx = distparam(2)
    muxpx = distparam(3)
    xscale = distparam(4)
    pxscale = distparam(5)
    xmu1 = distparam(6)
    xmu2 = distparam(7)
    sigy = distparam(8)
    sigpy = distparam(9)
    muypy = distparam(10)
    yscale = distparam(11)
    pyscale = distparam(12)
    xmu3 = distparam(13)
    xmu4 = distparam(14)
    sigz = distparam(15)
    sigpz = distparam(16)
    muzpz = distparam(17)
    zscale = distparam(18)
    pzscale = distparam(19)
    xmu5 = distparam(20)
    xmu6 = distparam(21)

    call getsize_Pgrid2d(grid,totnp,npy,npx)

    call getpost_Pgrid2d(grid,myid,myidy,myidx)
    !        seedarray(1)=(1001+myid)*(myid+7)
    !        write(6,*)'seedarray=',seedarray
    !        call random_seed(PUT=seedarray(1:1))
    call random_number(x11)

    avgpts = Npt/(npx*npy)

    sig1 = sigx*xscale
    sig2 = sigpx*pxscale
    sig3 = sigy*yscale
    sig4 = sigpy*pyscale
    sig5 = sigz*zscale
    sig6 = sigpz*pzscale

    sq12=sqrt(1.-muxpx*muxpx)
    sq34=sqrt(1.-muypy*muypy)
    sq56=sqrt(1.-muzpz*muzpz)

    ! initial allocate 'avgpts' particles on each processor.
    !if(flagalloc.eq.1) then
    !  Pts1 = 0.0
    !else
    !  allocate(Pts1(6,avgpts))
    !  Pts1 = 0.0
    !endif

    !        allocate(x1(2,avgpts))
    !        allocate(x2(2,avgpts))
    !        allocate(x3(2,avgpts))
    !        call normVec(x1,avgpts)
    !        call normVec(x2,avgpts)
    !        call normVec(x3,avgpts)

    intvsamp = Npt
    !        intvsamp = 10
    allocate(x1(2,intvsamp))
    allocate(x2(2,intvsamp))
    allocate(x3(2,intvsamp))

    call normVec(x1,intvsamp)
    call normVec(x2,intvsamp)
    call normVec(x3,intvsamp)

    jlow = myid*avgpts + 1
    jhigh = (myid+1)*avgpts

    do j = 1, Npt
       if( (j.ge.jlow).and.(j.le.jhigh)) then
          i = j - jlow + 1
          k = j
          !x-px:
          !            call normdv(x1)
          !           Correct Gaussian distribution.
          !i = (j-1)*intvsamp + k
          !Pts1(1,i) = xmu1 + sig1*x1(1,k)/sq12
          !Pts1(2,i) = xmu2 + sig2*(-muxpx*x1(1,k)/sq12+x1(2,k))
          !miguel's sampling. here, muxpx->alphax,sigpx->betax
          Pts1(1,i) = xmu1 + sigx*x1(1,k)*xscale
          Pts1(2,i) = xmu2 + sigx  &
               *(-muxpx*x1(1,k)+x1(2,k))/sigpx*pxscale
          !y-py
          !            call normdv(x1)
          !           Correct Gaussian distribution.
          !Pts1(3,i) = xmu3 + sig3*x2(1,k)/sq34
          !Pts1(4,i) = xmu4 + sig4*(-muypy*x2(1,k)/sq34+x2(2,k))
          !miguel's sampling. here, muxpx->alphax,sigpx->betax
          Pts1(3,i) = xmu3 + sigy*x2(1,k)*yscale
          Pts1(4,i) = xmu4 + sigy  &
               *(-muypy*x2(1,k)+x2(2,k))/sigpy*pyscale
          !z-pz
          !            call normdv(x1)
          !           Correct Gaussian distribution.
          !Pts1(5,i) = xmu5 + sig5*x3(1,k)/sq56
          !Pts1(6,i) = xmu6 + sig6*(-muzpz*x3(1,k)/sq56+x3(2,k))
          !miguel's sampling. here, muxpx->alphax,sigpx->betax
          Pts1(5,i) = xmu5 + x3(1,k)*sig5
          Pts1(6,i) = xmu6 + x3(2,k)*sig6
       endif
    enddo

    deallocate(x1)
    deallocate(x2)
    deallocate(x3)

    !Nptlocal = avgpts

    t_kvdist = t_kvdist + elapsedtime_Timer(t0)

  end subroutine Gauss3test_Dist



  ! sample the particles with intial distribution 
  ! using rejection method. 
  subroutine Waterbag_Dist(Pts1,nparam,distparam,grid,flagalloc,&
       Npt,Nptlocal,ngroup)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: nparam,flagalloc,Npt,Nptlocal,ngroup
    double precision, dimension(nparam) :: distparam
    type (Pgrid2d), intent(in) :: grid
    double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
         sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
    double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
    double precision :: sig1,sig2,sig3,sig4,sig5,sig6
    double precision :: rootx,rooty,rootz,r1,r2,x1,x2
    double precision :: r3,r4,r5,r6,x3,x4,x5,x6
    integer :: totnp,npy,npx
    integer :: avgpts,numpts,isamz,isamy
    integer :: myid,myidx,myidy,iran,intvsamp
    !        integer seedarray(2)
    double precision :: t0,x11
    double precision, allocatable, dimension(:) :: ranum6

    call starttime_Timer(t0)

    sigx = distparam(1)
    sigpx = distparam(2)
    muxpx = distparam(3)
    xscale = distparam(4)
    pxscale = distparam(5)
    xmu1 = distparam(6)
    xmu2 = distparam(7)
    sigy = distparam(8)
    sigpy = distparam(9)
    muypy = distparam(10)
    yscale = distparam(11)
    pyscale = distparam(12)
    xmu3 = distparam(13)
    xmu4 = distparam(14)
    sigz = distparam(15)
    sigpz = distparam(16)
    muzpz = distparam(17)
    zscale = distparam(18)
    pzscale = distparam(19)
    xmu5 = distparam(20)
    xmu6 = distparam(21)

    call getsize_Pgrid2d(grid,totnp,npy,npx)

    call getpost_Pgrid2d(grid,myid,myidy,myidx)
    !        seedarray(1)=(1001+myid)*(myid+7)
    !        seedarray(2)=(101+2*myid)*(myid+4)
    !        write(6,*)'seedarray=',seedarray
    !        call random_seed(PUT=seedarray)
    call random_number(x11)

    avgpts = Npt/(npx*npy)

    sig1 = sigx*xscale
    sig2 = sigpx*pxscale
    sig3 = sigy*yscale
    sig4 = sigpy*pyscale
    sig5 = sigz*zscale
    sig6 = sigpz*pzscale

    rootx=sqrt(1.-muxpx*muxpx)
    rooty=sqrt(1.-muypy*muypy)
    rootz=sqrt(1.-muzpz*muzpz)

    ! initial allocate 'avgpts' particles on each processor.
    !if(flagalloc.eq.1) then
    !  Pts1 = 0.0
    !else
    !  allocate(Pts1(6,avgpts))
    !  Pts1 = 0.0
    !endif
    numpts = 0
    isamz = 0
    isamy = 0
    !intvsamp = avgpts
    intvsamp = 8
    allocate(ranum6(6*intvsamp))

    do 
       ! rejection sample.
10     continue 
       isamz = isamz + 1
       if(mod(isamz-1,intvsamp).eq.0) then
          call random_number(ranum6)
       endif
       iran = 6*mod(isamz-1,intvsamp)
       r1 = 2.0*ranum6(iran+1)-1.0
       r2 = 2.0*ranum6(iran+2)-1.0
       r3 = 2.0*ranum6(iran+3)-1.0
       r4 = 2.0*ranum6(iran+4)-1.0
       r5 = 2.0*ranum6(iran+5)-1.0
       r6 = 2.0*ranum6(iran+6)-1.0
       if(r1**2+r2**2+r3**2+r4**2+r5**2+r6**2.gt.1.0) goto 10
       isamy = isamy + 1
       numpts = numpts + 1
       if(numpts.gt.avgpts) exit
       !x-px:
       x1 = r1*sqrt(8.0)
       x2 = r2*sqrt(8.0)
       !Correct transformation.
       Pts1(1,numpts) = xmu1 + sig1*x1/rootx
       Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
       !Rob's transformation
       !Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
       !Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
       !y-py:
       x3 = r3*sqrt(8.0)
       x4 = r4*sqrt(8.0)
       !correct transformation
       Pts1(3,numpts) = xmu3 + sig3*x3/rooty
       Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
       !Rob's transformation
       !Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
       !Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
       !t-pt:
       x5 = r5*sqrt(8.0)
       x6 = r6*sqrt(8.0)
       !correct transformation
       Pts1(5,numpts) = xmu5 + sig5*x5/rootz
       Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
       !Rob's transformation
       !Pts1(5,numpts) = (xmu5 + sig5*x5)*zscale
       !Pts1(6,numpts) = (xmu6 + sig6*(muzpz*x5+rootz*x6))/zscale
    enddo

    deallocate(ranum6)

    !Nptlocal = avgpts

    t_kvdist = t_kvdist + elapsedtime_Timer(t0)

  end subroutine Waterbag_Dist



  subroutine KV3d_Dist(Pts1,nparam,distparam,grid,Npt,Nptlocal,ngroup)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: nparam,Npt,Nptlocal,ngroup
    double precision, dimension(nparam) :: distparam
    type (Pgrid2d), intent(in) :: grid
    double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
         sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
    double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
    double precision :: sig1,sig2,sig3,sig4,sig5,sig6
    double precision :: rootx,rooty,rootz,r1,r2,x1,x2
    double precision :: r3,r4,r5,r6,x3,x4,x5,x6
    integer :: totnp,npy,npx
    integer :: avgpts,numpts
    integer :: myid,myidx,myidy
    !        integer seedarray(1)
    double precision :: t0,x11,twopi

    call starttime_Timer(t0)

    sigx = distparam(1)
    sigpx = distparam(2)
    muxpx = distparam(3)
    xscale = distparam(4)
    pxscale = distparam(5)
    xmu1 = distparam(6)
    xmu2 = distparam(7)
    sigy = distparam(8)
    sigpy = distparam(9)
    muypy = distparam(10)
    yscale = distparam(11)
    pyscale = distparam(12)
    xmu3 = distparam(13)
    xmu4 = distparam(14)
    sigz = distparam(15)
    sigpz = distparam(16)
    muzpz = distparam(17)
    zscale = distparam(18)
    pzscale = distparam(19)
    xmu5 = distparam(20)
    xmu6 = distparam(21)

    call getsize_Pgrid2d(grid,totnp,npy,npx)

    call getpost_Pgrid2d(grid,myid,myidy,myidx)
    !        seedarray(1)=(1001+myid)*(myid+7)
    !        write(6,*)'seedarray=',seedarray
    !        call random_seed(PUT=seedarray(1:1))
    call random_number(x11)

    avgpts = Npt/(npx*npy)

    sig1 = sigx*xscale
    sig2 = sigpx*pxscale
    sig3 = sigy*yscale
    sig4 = sigpy*pyscale
    sig5 = sigz*zscale
    sig6 = sigpz*pzscale

    rootx=sqrt(1.-muxpx*muxpx)
    rooty=sqrt(1.-muypy*muypy)
    rootz=sqrt(1.-muzpz*muzpz)

    ! initial allocate 'avgpts' particles on each processor.
    !allocate(Pts1(6,avgpts))
    Pts1 = 0.0
    twopi = 4*asin(1.0d0)

    do numpts = 1, avgpts
       call random_number(r1)
       call random_number(r2)
       call random_number(r3)
       r4 = sqrt(r1)
       r5 = sqrt(1.0-r1)
       r2 = r2*twopi
       r3 = r3*twopi
       x1 = 2*r4*cos(r2)
       x2 = 2*r4*sin(r2)
       x3 = 2*r5*cos(r3)
       x4 = 2*r5*sin(r3)
       !x-px:
       !Correct transformation.
       Pts1(1,numpts) = xmu1 + sig1*x1/rootx
       Pts1(2,numpts) = xmu2 + sig2*(-muxpx*x1/rootx+x2)
       !Rob's transformation.
       !Pts1(1,numpts) = (xmu1 + sig1*x1)*xscale
       !Pts1(2,numpts) = (xmu2 + sig2*(muxpx*x1+rootx*x2))/xscale
       !y-py:
       !correct transformation
       Pts1(3,numpts) = xmu3 + sig3*x3/rooty
       Pts1(4,numpts) = xmu4 + sig4*(-muypy*x3/rooty+x4)
       !Rob's transformation
       !Pts1(3,numpts) = (xmu3 + sig3*x3)*yscale
       !Pts1(4,numpts) = (xmu4 + sig4*(muypy*x3+rooty*x4))/yscale
       !t-pt:
       call random_number(r5)
       r5 = 2*r5 - 1.0
       call random_number(r6)
       r6 = 2*r6 - 1.0
       x5 = r5*sqrt(3.0)
       x6 = r6*sqrt(3.0)
       !correct transformation
       Pts1(5,numpts) = xmu5 + sig5*x5/rootz
       Pts1(6,numpts) = xmu6 + sig6*(-muzpz*x5/rootz+x6)
       !Rob's transformation
       !Pts1(5,numpts) = (xmu5 + sig5*x5)*zscale
       !Pts1(6,numpts) = (xmu6 + sig6*(muzpz*x5+rootz*x6))/zscale
    enddo

    !Nptlocal = avgpts

    t_kvdist = t_kvdist + elapsedtime_Timer(t0)

  end subroutine KV3d_Dist



  subroutine normdv(y)
    implicit none
    include 'mpif.h'
    double precision, dimension(2), intent(out) :: y
    double precision :: twopi,x1,x2,epsilon

    epsilon = 1.0e-18

    twopi = 4.0*asin(1.0d0)
    call random_number(x2)
10  call random_number(x1)
    !        x1 = 0.5
    !10      x2 = 0.6
    if(x1.eq.0.0) goto 10
    !        if(x1.eq.0.0) x1 = epsilon
    y(1) = sqrt(-2.0*log(x1))*cos(twopi*x2)
    y(2) = sqrt(-2.0*log(x1))*sin(twopi*x2)

  end subroutine normdv



  subroutine normVec(y,num)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: num
    double precision, dimension(2,num), intent(out) :: y
    double precision :: twopi,epsilon
    double precision, dimension(num) :: x1,x2
    integer :: i

    epsilon = 1.0d-18
    twopi = 4.0*asin(1.0d0)

    call random_number(x2)
    call random_number(x1)
    do i = 1, num
       if(x1(i).eq.0.0d0) x1(i) = epsilon
       y(1,i) = sqrt(-2.0*log(x1(i)))*cos(twopi*x2(i))
       y(2,i) = sqrt(-2.0*log(x1(i)))*sin(twopi*x2(i))
    enddo

  end subroutine normVec



  subroutine normVecnew(iseed,y,num)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: num
    integer, intent(inout) :: iseed
    double precision, dimension(2,num), intent(out) :: y
    double precision :: twopi,epsilon
    double precision, dimension(num) :: x1,x2
    integer :: i

    epsilon = 1.0e-18

    twopi = 4.0*asin(1.0d0)
    call ran2(iseed,x2,num)
    call ran2(iseed,x1,num)
    do i = 1, num
       if(x1(i).eq.0.0) x1(i) = epsilon
       y(1,i) = sqrt(-2.0*log(x1(i)))*cos(twopi*x2(i))
       y(2,i) = sqrt(-2.0*log(x1(i)))*sin(twopi*x2(i))
    enddo

  end subroutine normVecnew



  subroutine regen_Dist(Pts1,myidx,myidy,npx,npyhalf,nptot)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: myidy,npyhalf,myidx,npx
    integer, intent(out) :: nptot
    integer :: i,j,jlow,jhigh,avgpts
    double precision, dimension(6) :: tmptcl
    double precision :: sum1,sum2,sum3,sum4

    open(unit=12,file='partcl.data1',status='old')
    open(unit=13,file='partcl.data2',status='old')

    sum1 = 0.0
    sum2 = 0.0
    sum3 = 0.0
    sum4 = 0.0

    if(myidy.lt.npyhalf) then
       read(12,*)nptot
       avgpts = nptot/npyhalf/npx
       jlow = myidx*npyhalf*avgpts + myidy*avgpts + 1
       jhigh = myidx*npyhalf*avgpts + (myidy+1)*avgpts
       do j = 1, nptot
          read(12,*)tmptcl(1:6)
          sum1 = sum1 + tmptcl(1)
          sum2 = sum2 + tmptcl(3)
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
             i = j - jlow + 1
             Pts1(1:6,i) = tmptcl(1:6)
          endif
       enddo
    else
       read(13,*)nptot
       avgpts = nptot/npyhalf/npx
       jlow = myidx*npyhalf*avgpts+(myidy-npyhalf)*avgpts + 1
       jhigh = myidx*npyhalf*avgpts+(myidy-npyhalf+1)*avgpts
       do j = 1, nptot
          read(13,*)tmptcl(1:6)
          sum3 = sum3 + tmptcl(1)
          sum4 = sum4 + tmptcl(3)
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
             i = j - jlow + 1
             Pts1(1:6,i) = tmptcl(1:6)
          endif
       enddo
    endif

    close(12)
    close(13)

  end subroutine regen_Dist



  subroutine normVecSob(y1,y2,y3,num,iseed)
    implicit none
    include 'mpif.h'
    integer :: num
    double precision, dimension(2,num) :: y1,y2,y3
    double precision :: twopi,epsilon
    integer :: i,ndim,iseed
    real*8, dimension(6) :: xtmp 

    epsilon = 1.0e-18
    ndim = 6
    twopi = 4.0*asin(1.0d0)

    do i = 1, num
       call sobseq(ndim,xtmp,iseed)
       if(xtmp(1).eq.0.0) xtmp(1) = epsilon
       if(xtmp(3).eq.0.0) xtmp(3) = epsilon
       if(xtmp(5).eq.0.0) xtmp(5) = epsilon
       y1(1,i) = sqrt(-2.0*log(xtmp(1)))*cos(twopi*xtmp(2))
       y1(2,i) = sqrt(-2.0*log(xtmp(1)))*sin(twopi*xtmp(2))
       y2(1,i) = sqrt(-2.0*log(xtmp(3)))*cos(twopi*xtmp(4))
       y2(2,i) = sqrt(-2.0*log(xtmp(3)))*sin(twopi*xtmp(4))
       y3(1,i) = sqrt(-2.0*log(xtmp(5)))*cos(twopi*xtmp(6))
       y3(2,i) = sqrt(-2.0*log(xtmp(5)))*sin(twopi*xtmp(6))
    enddo

  end subroutine normVecSob



  SUBROUTINE sobseq(n,x,iseed)
    INTEGER n,MAXBIT,MAXDIM,iseed
    REAL*8 x(*)
    PARAMETER (MAXBIT=30,MAXDIM=6)
    INTEGER i,im,in,ipp,j,k,l,ip(MAXDIM),iu(MAXDIM,MAXBIT),&
         iv(MAXBIT*MAXDIM),ix(MAXDIM),mdeg(MAXDIM)
    REAL*8 fac
    SAVE ip,mdeg,ix,iv,in,fac
    EQUIVALENCE (iv,iu)
    DATA ip /0,1,1,2,1,4/, mdeg /1,2,3,3,4,4/, ix /6*0/
    DATA iv /6*1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9,156*0/
    if (n.lt.0) then
       do 14 k=1,MAXDIM
          do 11 j=1,mdeg(k)
             iu(k,j)=iu(k,j)*2**(MAXBIT-j)
11           continue
          do 13 j=mdeg(k)+1,MAXBIT
             ipp=ip(k)
             i=iu(k,j-mdeg(k))
             i=ieor(i,i/2**mdeg(k))
             do 12 l=mdeg(k)-1,1,-1
                if(iand(ipp,1).ne.0)i=ieor(i,iu(k,j-l))
                ipp=ipp/2
12              continue
             iu(k,j)=i
13           continue
14        continue
       fac=1./2.**MAXBIT
       in=iseed
    else
       im=in
       do 15 j=1,MAXBIT
          if(iand(im,1).eq.0)goto 1
          im=im/2
15        continue
       pause 'MAXBIT too small in sobseq'
1      im=(j-1)*MAXDIM
       do 16 k=1,min(n,MAXDIM)
          ix(k)=ieor(ix(k),iv(im+k))
          x(k)=ix(k)*fac
16     continue
       in=in+1
    endif
    return
  END subroutine sobseq



  subroutine read2g_Dist(Pts1,nparam,distparam,grid,flagalloc,&
       Npt,Nptlocal,ngroup)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: nparam,flagalloc,Npt,Nptlocal,ngroup
    double precision, dimension(nparam) :: distparam
    type (Pgrid2d), intent(in) :: grid
    double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
         sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
    double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
    integer :: totnp,npy,npx,nproc
    integer :: avgpts
    integer :: myid,myidx,myidy,i,j,myidtmp
    !        integer seedarray(1)
    double precision :: t0,x11
    double precision, dimension(6) :: tmptcl
    integer :: npyhalf
    integer :: nptot,nleft,jlow,jhigh

    call starttime_Timer(t0)

    sigx = distparam(1)
    sigpx = distparam(2)
    muxpx = distparam(3)
    xscale = distparam(4)
    pxscale = distparam(5)
    xmu1 = distparam(6)
    xmu2 = distparam(7)
    sigy = distparam(8)
    sigpy = distparam(9)
    muypy = distparam(10)
    yscale = distparam(11)
    pyscale = distparam(12)
    xmu3 = distparam(13)
    xmu4 = distparam(14)
    sigz = distparam(15)
    sigpz = distparam(16)
    muzpz = distparam(17)
    zscale = distparam(18)
    pzscale = distparam(19)
    xmu5 = distparam(20)
    xmu6 = distparam(21)

    call getsize_Pgrid2d(grid,totnp,npy,npx)
    call getpost_Pgrid2d(grid,myid,myidy,myidx)

    !        seedarray(1)=(1001+myid)*(myid+7)
    !        write(6,*)'seedarray=',seedarray
    !        call random_seed(PUT=seedarray(1:1))
    call random_number(x11)

    avgpts = Npt/(npx*npy/ngroup)

    npyhalf = npy/ngroup
    nproc = npx*npy/ngroup

    open(unit=11,file='partcl1.data',status='old')
    open(unit=12,file='partcl2.data',status='old')

    if(myidy.lt.npyhalf) then
       read(11,*)nptot
       if(nptot.ne.Npt) print*,"input particle # is not correct!"
       avgpts = nptot/nproc
       nleft = nptot - avgpts*nproc
       if(myid.lt.nleft) then
          avgpts = avgpts+1
          jlow = myid*avgpts + 1
          jhigh = (myid+1)*avgpts
       else
          jlow = myid*avgpts + 1 + nleft
          jhigh = (myid+1)*avgpts + nleft
       endif
       allocate(Pts1(6,avgpts))
       Pts1 = 0.0
       do j = 1, nptot
          read(11,*)tmptcl(1:6)
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
             i = j - jlow + 1
             Pts1(1:6,i) = tmptcl(1:6)
          endif
       enddo
    else
       read(12,*)nptot
       if(nptot.ne.Npt) print*,"input particle # is not correct!"
       avgpts = nptot/nproc
       nleft = nptot - avgpts*nproc
       myidtmp = myid-nproc
       if(myidtmp.lt.nleft) then
          avgpts = avgpts+1
          jlow = myidtmp*avgpts + 1
          jhigh = (myidtmp+1)*avgpts
       else
          jlow = myidtmp*avgpts + 1 + nleft
          jhigh = (myidtmp+1)*avgpts + nleft
       endif
       allocate(Pts1(6,avgpts))
       Pts1 = 0.0
       do j = 1, nptot
          read(12,*)tmptcl(1:6)
          if( (j.ge.jlow).and.(j.le.jhigh) ) then
             i = j - jlow + 1
             Pts1(1:6,i) = tmptcl(1:6)
          endif
       enddo
    endif

    !Nptlocal = avgpts

    t_kvdist = t_kvdist + elapsedtime_Timer(t0)
    close(11)
    close(12)

  end subroutine read2g_Dist


  
  subroutine read1g_Dist(Pts1,nparam,distparam,grid,&
       Npt,Nptlocal,ngroup,filename)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: nparam,Npt,Nptlocal,ngroup
    double precision, dimension(nparam) :: distparam
    type (Pgrid2d), intent(in) :: grid
    character*12 :: filename
    double precision  :: sigx,sigpx,muxpx,xscale,sigy,&
         sigpy,muypy, yscale,sigz,sigpz,muzpz,zscale,pxscale,pyscale,pzscale
    double precision :: xmu1,xmu2,xmu3,xmu4,xmu5,xmu6
    integer :: totnp,npy,npx,nproc
    integer :: avgpts
    integer :: myid,myidx,myidy,i,j
    !        integer seedarray(1)
    double precision :: t0,x11
    double precision, dimension(6) :: tmptcl
    integer :: npyhalf
    integer :: nptot,nleft,jlow,jhigh

    call starttime_Timer(t0)

    sigx = distparam(1)
    sigpx = distparam(2)
    muxpx = distparam(3)
    xscale = distparam(4)
    pxscale = distparam(5)
    xmu1 = distparam(6)
    xmu2 = distparam(7)
    sigy = distparam(8)
    sigpy = distparam(9)
    muypy = distparam(10)
    yscale = distparam(11)
    pyscale = distparam(12)
    xmu3 = distparam(13)
    xmu4 = distparam(14)
    sigz = distparam(15)
    sigpz = distparam(16)
    muzpz = distparam(17)
    zscale = distparam(18)
    pzscale = distparam(19)
    xmu5 = distparam(20)
    xmu6 = distparam(21)

    call getsize_Pgrid2d(grid,totnp,npy,npx)
    call getpost_Pgrid2d(grid,myid,myidy,myidx)

    !        seedarray(1)=(1001+myid)*(myid+7)
    !        write(6,*)'seedarray=',seedarray
    !        call random_seed(PUT=seedarray(1:1))
    call random_number(x11)

    avgpts = Npt/(npx*npy/ngroup)

    npyhalf = npy/ngroup
    nproc = npx*npy/ngroup


    open(unit=11,file=filename,status='old')

    read(11,*)nptot
    avgpts = nptot/nproc
    nleft = nptot - avgpts*nproc
    if(myid.lt.nleft) then
       avgpts = avgpts+1
       jlow = myid*avgpts + 1
       jhigh = (myid+1)*avgpts
    else
       jlow = myid*avgpts + 1 + nleft
       jhigh = (myid+1)*avgpts + nleft
    endif
    allocate(Pts1(6,avgpts))
    Pts1 = 0.0
    do j = 1, nptot
       read(11,*)tmptcl(1:6)
       if( (j.ge.jlow).and.(j.le.jhigh) ) then
          i = j - jlow + 1
          Pts1(1:6,i) = tmptcl(1:6)
       endif
    enddo

    !Nptlocal = avgpts

    t_kvdist = t_kvdist + elapsedtime_Timer(t0)
    close(11)

  end subroutine read1g_Dist


end module Distributionclass
