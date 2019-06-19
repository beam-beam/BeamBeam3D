!----------------------------------------------------------------
! (c) Copyright, 2011 by the Regents of the University of California.
! Outputclass: Output class in I/O module of CONTROL layer. 
! Version: 2.6
! Author: Ji Qiang, LBNL
! Description: This class defines functions to print out the charged
!              particle beam information in the accelerator.
! Comments:
!----------------------------------------------------------------

module Outputclass
  use Timerclass
  use Pgrid2dclass
  use Communclass
  use utilityclass

contains

  ! calculate <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
  ! <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance using 2 group PEs.
  subroutine diagnostic2G_Output(z,Pts1,innp,nptot,&
       myidx,myidy,npx,npy,commrow,commcol,comm2d)
    implicit none
    include 'mpif.h'
    double precision, intent(in) :: z
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot,myidx,myidy,npx,npy,commrow,&
         commcol,comm2d
    double precision:: den1,sqsum1,sqsum2,sqsum3,sqsum4,&
         epsx2,epsy2
    double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
         xpxfac,ypyfac
    double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
    double precision:: xpxlocal,ypylocal,zpzlocal
    double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
    double precision:: sqsum5local,sqsum6local
    double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
         pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
         z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
         z03,z04,pz03,pz04
    double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
         sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
    integer :: i,my_rank,ierr,j
    double precision, dimension(6) :: localmax
    double precision, dimension(14) :: localmax2, glmax
    double precision, dimension(54) :: tmplc,tmpgl
    double precision :: t0,lcrmax
    !        double precision :: alphax,alphay,alphaz
    integer :: npctmin,npctmax,npyhalf

    call starttime_Timer(t0)

    npctmin = 1
    npctmax = 1

    call MPI_COMM_RANK(mpicommwd,my_rank,ierr)
    !        call MPI_COMM_SIZE(mpicommwd,nproc,ierr)
    npyhalf = npy/2 

    den1 = 1.0/dble(nptot)
    x0lc = 0.0
    px0lc = 0.0
    y0lc = 0.0
    py0lc = 0.0
    z0lc = 0.0
    pz0lc = 0.0
    sqsum1local = 0.0
    sqsum2local = 0.0
    sqsum3local = 0.0
    sqsum4local = 0.0
    sqsum5local = 0.0
    sqsum6local = 0.0
    xpxlocal = 0.0
    ypylocal = 0.0
    zpzlocal = 0.0
    x0lc3 = 0.0
    x0lc4 = 0.0
    px0lc3 = 0.0
    px0lc4 = 0.0
    y0lc3 = 0.0
    y0lc4 = 0.0
    py0lc3 = 0.0
    py0lc4 = 0.0
    z0lc3 = 0.0
    z0lc4 = 0.0
    pz0lc3 = 0.0
    pz0lc4 = 0.0

    ! for cache optimization.
    if(innp.ne.0) then
       do i = 1, 6
          localmax(i) = abs(Pts1(i,1))
       enddo
       lcrmax = Pts1(1,1)**2+Pts1(3,1)**2
    else
       do i = 1, 6
          localmax(i) = 0.0
       enddo
       lcrmax = 0.0
    endif
    do i = 1, innp
       x0lc = x0lc + Pts1(1,i)
       sqsum1local = sqsum1local + Pts1(1,i)*Pts1(1,i)
       x0lc3 = x0lc3 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)
       x0lc4 = x0lc4 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)*&
            Pts1(1,i)
       xpxlocal = xpxlocal + Pts1(1,i)*Pts1(2,i)
       px0lc = px0lc + Pts1(2,i)
       sqsum2local = sqsum2local + Pts1(2,i)*Pts1(2,i)
       px0lc3 = px0lc3 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)
       px0lc4 = px0lc4 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)*&
            Pts1(2,i)
       y0lc = y0lc + Pts1(3,i)
       sqsum3local = sqsum3local + Pts1(3,i)*Pts1(3,i)
       y0lc3 = y0lc3 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)
       y0lc4 = y0lc4 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)*&
            Pts1(3,i)
       ypylocal = ypylocal + Pts1(3,i)*Pts1(4,i)
       py0lc = py0lc + Pts1(4,i)
       sqsum4local = sqsum4local + Pts1(4,i)*Pts1(4,i)
       py0lc3 = py0lc3 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)
       py0lc4 = py0lc4 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)*&
            Pts1(4,i)
       z0lc = z0lc + Pts1(5,i)
       sqsum5local = sqsum5local + Pts1(5,i)*Pts1(5,i)
       z0lc3 = z0lc3 + abs(Pts1(5,i)*Pts1(5,i)*Pts1(5,i))
       z0lc4 = z0lc4 + Pts1(5,i)*Pts1(5,i)*Pts1(5,i)*&
            Pts1(5,i)

       zpzlocal = zpzlocal + Pts1(5,i)*Pts1(6,i)
       pz0lc = pz0lc + Pts1(6,i)
       sqsum6local = sqsum6local + Pts1(6,i)*Pts1(6,i)
       pz0lc3 = pz0lc3 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)
       pz0lc4 = pz0lc4 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)*&
            Pts1(6,i)
       do j = 1, 6
          if(localmax(j).lt.abs(Pts1(j,i))) then
             localmax(j) = abs(Pts1(j,i))
          endif
       enddo
       if(lcrmax.lt.(Pts1(1,i)**2+Pts1(3,i)**2)) then
          lcrmax = Pts1(1,i)**2 + Pts1(3,i)**2
       endif
    enddo

    localmax2 = 0.0
    tmplc = 0.0
    if(myidy.lt.npyhalf) then
       tmplc(1) = x0lc
       tmplc(2) = px0lc
       tmplc(3) = y0lc
       tmplc(4) = py0lc
       tmplc(5) = z0lc
       tmplc(6) = pz0lc
       tmplc(7) = sqsum1local
       tmplc(8) = sqsum2local
       tmplc(9) = sqsum3local
       tmplc(10) = sqsum4local
       tmplc(11) = sqsum5local
       tmplc(12) = sqsum6local
       tmplc(13) = xpxlocal
       tmplc(14) = ypylocal
       tmplc(15) = zpzlocal
       tmplc(16) = x0lc3
       tmplc(17) = x0lc4
       tmplc(18) = px0lc3
       tmplc(19) = px0lc4
       tmplc(20) = y0lc3
       tmplc(21) = y0lc4
       tmplc(22) = py0lc3
       tmplc(23) = py0lc4
       tmplc(24) = z0lc3
       tmplc(25) = z0lc4
       tmplc(26) = pz0lc3
       tmplc(27) = pz0lc4
       localmax2(1) = localmax(1)
       localmax2(2) = localmax(2)
       localmax2(3) = localmax(3)
       localmax2(4) = localmax(4)
       localmax2(5) = localmax(5)
       localmax2(6) = localmax(6)
       localmax2(7) = lcrmax
    else
       tmplc(28) = x0lc
       tmplc(29) = px0lc
       tmplc(30) = y0lc
       tmplc(31) = py0lc
       tmplc(32) = z0lc
       tmplc(33) = pz0lc
       tmplc(34) = sqsum1local
       tmplc(35) = sqsum2local
       tmplc(36) = sqsum3local
       tmplc(37) = sqsum4local
       tmplc(38) = sqsum5local
       tmplc(39) = sqsum6local
       tmplc(40) = xpxlocal
       tmplc(41) = ypylocal
       tmplc(42) = zpzlocal
       tmplc(43) = x0lc3
       tmplc(44) = x0lc4
       tmplc(45) = px0lc3
       tmplc(46) = px0lc4
       tmplc(47) = y0lc3
       tmplc(48) = y0lc4
       tmplc(49) = py0lc3
       tmplc(50) = py0lc4
       tmplc(51) = z0lc3
       tmplc(52) = z0lc4
       tmplc(53) = pz0lc3
       tmplc(54) = pz0lc4
       localmax2(8) = localmax(1)
       localmax2(9) = localmax(2)
       localmax2(10) = localmax(3)
       localmax2(11) = localmax(4)
       localmax2(12) = localmax(5)
       localmax2(13) = localmax(6)
       localmax2(14) = lcrmax
    endif

    call MPI_REDUCE(tmplc,tmpgl,54,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpicommwd,ierr)
    call MPI_REDUCE(localmax2,glmax,14,MPI_DOUBLE_PRECISION,MPI_MAX,0,mpicommwd,ierr)

    if(my_rank.eq.0) then
       !//output beam 1
       x0 = tmpgl(1)*den1
       px0 = tmpgl(2)*den1
       y0 = tmpgl(3)*den1
       py0 = tmpgl(4)*den1
       z0 = tmpgl(5)*den1
       pz0 = tmpgl(6)*den1
       sqx = tmpgl(7)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(8)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(9)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(10)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(11)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(12)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(13)*den1 - x0*px0
       ypy = tmpgl(14)*den1 - y0*py0
       zpz = tmpgl(15)*den1 - z0*pz0
       cubx = tmpgl(16)*den1
       fthx = tmpgl(17)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(18)*den1
       fthpx = tmpgl(19)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(20)*den1
       fthy = tmpgl(21)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(22)*den1
       fthpy = tmpgl(23)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(24)*den1
       fthz = tmpgl(25)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(26)*den1
       fthpz = tmpgl(27)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0
       ypyfac = 0.0
       zpzfac = 0.0
       if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
       if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
       if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
       write(24,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
       write(25,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
       write(26,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
       write(27,102)z,glmax(1),glmax(2),glmax(3),glmax(4),glmax(5),&
            glmax(6),sqrt(glmax(7))
       write(28,101)z,npctmin,npctmax,nptot
       write(29,100)z,x03,px03,y03,py03,z03,pz03
       write(30,100)z,x04,px04,y04,py04,z04,pz04

       call flush(24)
       call flush(25)
       call flush(26)
       call flush(27)
       call flush(28)
       call flush(29)
       call flush(30)

       !//output beam 2
       x0 = tmpgl(28)*den1
       px0 = tmpgl(29)*den1
       y0 = tmpgl(30)*den1
       py0 = tmpgl(31)*den1
       z0 = tmpgl(32)*den1
       pz0 = tmpgl(33)*den1
       sqx = tmpgl(34)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(35)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(36)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(37)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(38)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(39)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(40)*den1 - x0*px0
       ypy = tmpgl(41)*den1 - y0*py0
       zpz = tmpgl(42)*den1 - z0*pz0
       cubx = tmpgl(43)*den1
       fthx = tmpgl(44)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(45)*den1
       fthpx = tmpgl(46)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(47)*den1
       fthy = tmpgl(48)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(49)*den1
       fthpy = tmpgl(50)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(51)*den1
       fthz = tmpgl(52)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(53)*den1
       fthpz = tmpgl(54)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0
       ypyfac = 0.0
       zpzfac = 0.0
       if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
       if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
       if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
       write(34,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
       write(35,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
       write(36,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
       write(37,102)z,glmax(8),glmax(9),glmax(10),glmax(11),glmax(12),&
            glmax(13),sqrt(glmax(14))
       write(38,101)z,npctmin,npctmax,nptot
       write(39,100)z,x03,px03,y03,py03,z03,pz03
       write(40,100)z,x04,px04,y04,py04,z04,pz04

       call flush(34)
       call flush(35)
       call flush(36)
       call flush(37)
       call flush(38)
       call flush(39)
       call flush(40)
    endif

100 format(7(1x,es13.6))
101 format(1x,es13.6,3I10)
102 format(8(1x,es13.6))

    t_diag = t_diag + elapsedtime_Timer(t0)

  end subroutine diagnostic2G_Output



  !Multiple bunch output
  ! calculate <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
  ! <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance using 2 group PEs.
  subroutine diagnostic2GMBold_Output(z,Pts1,innp,nptot,&
       myidx,myidy,npx,npy,commrow,commcol,comm2d,bid)
    implicit none
    include 'mpif.h'
    double precision, intent(in) :: z
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot,myidx,myidy,npx,npy,commrow,&
         commcol,comm2d,bid
    double precision:: den1,sqsum1,sqsum2,sqsum3,sqsum4,&
         epsx2,epsy2
    double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
         xpxfac,ypyfac
    double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
    double precision:: xpxlocal,ypylocal,zpzlocal
    double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
    double precision:: sqsum5local,sqsum6local
    double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
         pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
         z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
         z03,z04,pz03,pz04
    double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
         sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
    integer :: i,my_rank,ierr,j
    double precision, dimension(6) :: localmax
    double precision, dimension(14) :: localmax2, glmax
    double precision, dimension(54) :: tmplc,tmpgl
    double precision :: t0,lcrmax
    !        double precision :: alphax,alphay,alphaz
    integer :: npctmin,npctmax,npyhalf,nfile

    call starttime_Timer(t0)

    npctmin = 1
    npctmax = 1

    call MPI_COMM_RANK(mpicommwd,my_rank,ierr)
    !        call MPI_COMM_SIZE(mpicommwd,nproc,ierr)
    npyhalf = npy/2 

    den1 = 1.0/dble(nptot)
    x0lc = 0.0
    px0lc = 0.0
    y0lc = 0.0
    py0lc = 0.0
    z0lc = 0.0
    pz0lc = 0.0
    sqsum1local = 0.0
    sqsum2local = 0.0
    sqsum3local = 0.0
    sqsum4local = 0.0
    sqsum5local = 0.0
    sqsum6local = 0.0
    xpxlocal = 0.0
    ypylocal = 0.0
    zpzlocal = 0.0
    x0lc3 = 0.0
    x0lc4 = 0.0
    px0lc3 = 0.0
    px0lc4 = 0.0
    y0lc3 = 0.0
    y0lc4 = 0.0
    py0lc3 = 0.0
    py0lc4 = 0.0
    z0lc3 = 0.0
    z0lc4 = 0.0
    pz0lc3 = 0.0
    pz0lc4 = 0.0

    ! for cache optimization.
    if(innp.ne.0) then
       do i = 1, 6
          localmax(i) = abs(Pts1(i,1))
       enddo
       lcrmax = Pts1(1,1)**2+Pts1(3,1)**2
    else
       do i = 1, 6
          localmax(i) = 0.0
       enddo
       lcrmax = 0.0
    endif
    do i = 1, innp
       x0lc = x0lc + Pts1(1,i)
       sqsum1local = sqsum1local + Pts1(1,i)*Pts1(1,i)
       x0lc3 = x0lc3 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)
       x0lc4 = x0lc4 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)*&
            Pts1(1,i)
       xpxlocal = xpxlocal + Pts1(1,i)*Pts1(2,i)
       px0lc = px0lc + Pts1(2,i)
       sqsum2local = sqsum2local + Pts1(2,i)*Pts1(2,i)
       px0lc3 = px0lc3 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)
       px0lc4 = px0lc4 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)*&
            Pts1(2,i)
       y0lc = y0lc + Pts1(3,i)
       sqsum3local = sqsum3local + Pts1(3,i)*Pts1(3,i)
       y0lc3 = y0lc3 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)
       y0lc4 = y0lc4 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)*&
            Pts1(3,i)
       ypylocal = ypylocal + Pts1(3,i)*Pts1(4,i)
       py0lc = py0lc + Pts1(4,i)
       sqsum4local = sqsum4local + Pts1(4,i)*Pts1(4,i)
       py0lc3 = py0lc3 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)
       py0lc4 = py0lc4 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)*&
            Pts1(4,i)
       z0lc = z0lc + Pts1(5,i)
       sqsum5local = sqsum5local + Pts1(5,i)*Pts1(5,i)
       z0lc3 = z0lc3 + abs(Pts1(5,i)*Pts1(5,i)*Pts1(5,i))
       z0lc4 = z0lc4 + Pts1(5,i)*Pts1(5,i)*Pts1(5,i)*&
            Pts1(5,i)

       zpzlocal = zpzlocal + Pts1(5,i)*Pts1(6,i)
       pz0lc = pz0lc + Pts1(6,i)
       sqsum6local = sqsum6local + Pts1(6,i)*Pts1(6,i)
       pz0lc3 = pz0lc3 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)
       pz0lc4 = pz0lc4 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)*&
            Pts1(6,i)
       do j = 1, 6
          if(localmax(j).lt.abs(Pts1(j,i))) then
             localmax(j) = abs(Pts1(j,i))
          endif
       enddo
       if(lcrmax.lt.(Pts1(1,i)**2+Pts1(3,i)**2)) then
          lcrmax = Pts1(1,i)**2 + Pts1(3,i)**2
       endif
    enddo

    localmax2 = 0.0
    tmplc = 0.0
    if(myidy.lt.npyhalf) then
       tmplc(1) = x0lc
       tmplc(2) = px0lc
       tmplc(3) = y0lc
       tmplc(4) = py0lc
       tmplc(5) = z0lc
       tmplc(6) = pz0lc
       tmplc(7) = sqsum1local
       tmplc(8) = sqsum2local
       tmplc(9) = sqsum3local
       tmplc(10) = sqsum4local
       tmplc(11) = sqsum5local
       tmplc(12) = sqsum6local
       tmplc(13) = xpxlocal
       tmplc(14) = ypylocal
       tmplc(15) = zpzlocal
       tmplc(16) = x0lc3
       tmplc(17) = x0lc4
       tmplc(18) = px0lc3
       tmplc(19) = px0lc4
       tmplc(20) = y0lc3
       tmplc(21) = y0lc4
       tmplc(22) = py0lc3
       tmplc(23) = py0lc4
       tmplc(24) = z0lc3
       tmplc(25) = z0lc4
       tmplc(26) = pz0lc3
       tmplc(27) = pz0lc4
       localmax2(1) = localmax(1)
       localmax2(2) = localmax(2)
       localmax2(3) = localmax(3)
       localmax2(4) = localmax(4)
       localmax2(5) = localmax(5)
       localmax2(6) = localmax(6)
       localmax2(7) = lcrmax
    else
       tmplc(28) = x0lc
       tmplc(29) = px0lc
       tmplc(30) = y0lc
       tmplc(31) = py0lc
       tmplc(32) = z0lc
       tmplc(33) = pz0lc
       tmplc(34) = sqsum1local
       tmplc(35) = sqsum2local
       tmplc(36) = sqsum3local
       tmplc(37) = sqsum4local
       tmplc(38) = sqsum5local
       tmplc(39) = sqsum6local
       tmplc(40) = xpxlocal
       tmplc(41) = ypylocal
       tmplc(42) = zpzlocal
       tmplc(43) = x0lc3
       tmplc(44) = x0lc4
       tmplc(45) = px0lc3
       tmplc(46) = px0lc4
       tmplc(47) = y0lc3
       tmplc(48) = y0lc4
       tmplc(49) = py0lc3
       tmplc(50) = py0lc4
       tmplc(51) = z0lc3
       tmplc(52) = z0lc4
       tmplc(53) = pz0lc3
       tmplc(54) = pz0lc4
       localmax2(8) = localmax(1)
       localmax2(9) = localmax(2)
       localmax2(10) = localmax(3)
       localmax2(11) = localmax(4)
       localmax2(12) = localmax(5)
       localmax2(13) = localmax(6)
       localmax2(14) = lcrmax
    endif

    call MPI_REDUCE(tmplc,tmpgl,54,MPI_DOUBLE_PRECISION,&
         MPI_SUM,0,mpicommwd,ierr)
    call MPI_REDUCE(localmax2,glmax,14,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
         mpicommwd,ierr)

    if(my_rank.eq.0) then
       !//output beam 1
       x0 = tmpgl(1)*den1
       px0 = tmpgl(2)*den1
       y0 = tmpgl(3)*den1
       py0 = tmpgl(4)*den1
       z0 = tmpgl(5)*den1
       pz0 = tmpgl(6)*den1
       sqx = tmpgl(7)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(8)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(9)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(10)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(11)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(12)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(13)*den1 - x0*px0
       ypy = tmpgl(14)*den1 - y0*py0
       zpz = tmpgl(15)*den1 - z0*pz0
       cubx = tmpgl(16)*den1
       fthx = tmpgl(17)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(18)*den1
       fthpx = tmpgl(19)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(20)*den1
       fthy = tmpgl(21)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(22)*den1
       fthpy = tmpgl(23)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(24)*den1
       fthz = tmpgl(25)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(26)*den1
       fthpz = tmpgl(27)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0
       ypyfac = 0.0
       zpzfac = 0.0
       if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
       if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
       if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
       nfile = 4+bid*20
       write(nfile,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
       call flush(nfile)
       nfile = 5+bid*20
       write(nfile,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
       call flush(nfile)
       nfile = 6+bid*20
       write(nfile,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
       call flush(nfile)
       nfile = 7+bid*20
       write(nfile,102)z,glmax(1),glmax(2),glmax(3),glmax(4),glmax(5),&
            glmax(6),sqrt(glmax(7))
       call flush(nfile)
       nfile = 8+bid*20
       write(nfile,101)z,npctmin,npctmax,nptot
       call flush(nfile)
       nfile = 9+bid*20
       write(nfile,100)z,x03,px03,y03,py03,z03,pz03
       call flush(nfile)
       nfile = 10+bid*20
       write(nfile,100)z,x04,px04,y04,py04,z04,pz04
       call flush(nfile)

       !//output beam 2
       x0 = tmpgl(28)*den1
       px0 = tmpgl(29)*den1
       y0 = tmpgl(30)*den1
       py0 = tmpgl(31)*den1
       z0 = tmpgl(32)*den1
       pz0 = tmpgl(33)*den1
       sqx = tmpgl(34)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(35)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(36)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(37)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(38)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(39)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(40)*den1 - x0*px0
       ypy = tmpgl(41)*den1 - y0*py0
       zpz = tmpgl(42)*den1 - z0*pz0
       cubx = tmpgl(43)*den1
       fthx = tmpgl(44)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(45)*den1
       fthpx = tmpgl(46)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(47)*den1
       fthy = tmpgl(48)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(49)*den1
       fthpy = tmpgl(50)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(51)*den1
       fthz = tmpgl(52)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(53)*den1
       fthpz = tmpgl(54)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0
       ypyfac = 0.0
       zpzfac = 0.0
       if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
       if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
       if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
       nfile = 14+bid*20
       write(nfile,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
       call flush(nfile)
       nfile = 15+bid*20
       write(nfile,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
       call flush(nfile)
       nfile = 16+bid*20
       write(nfile,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
       call flush(nfile)
       nfile = 17+bid*20
       write(nfile,102)z,glmax(8),glmax(9),glmax(10),glmax(11),glmax(12),&
            glmax(13),sqrt(glmax(14))
       call flush(nfile)
       nfile = 18+bid*20
       write(nfile,101)z,npctmin,npctmax,nptot
       call flush(nfile)
       nfile = 19+bid*20
       write(nfile,100)z,x03,px03,y03,py03,z03,pz03
       call flush(nfile)
       nfile = 20+bid*20
       write(nfile,100)z,x04,px04,y04,py04,z04,pz04
       call flush(nfile)

    endif

100 format(7(1x,es13.6))
101 format(1x,es13.6,3I10)
102 format(8(1x,es13.6))

    t_diag = t_diag + elapsedtime_Timer(t0)

  end subroutine diagnostic2GMBold_Output

  subroutine diagnostic2GMBeic_Output(z,Pts1,innp,nptot,&
       myidx,myidy,npx,npy,commrow,commcol,comm2d,bid,alphamx,betamx,&
                           alphamy,betamy)
    implicit none
    include 'mpif.h'
    double precision, intent(in) :: z
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot,myidx,myidy,npx,npy,commrow,&
         commcol,comm2d,bid
    double precision:: den1,sqsum1,sqsum2,sqsum3,sqsum4,&
         epsx2,epsy2,gambet,betamx,alphamx,betamy,alphamy
    double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
         xpxfac,ypyfac
    double precision :: sqsum1local,sqsum2local,sqsum3local,sqsum4local
    double precision :: xpxlocal,ypylocal,zpzlocal
    double precision :: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
    double precision :: sqsum5local,sqsum6local
    double precision :: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
         pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
         z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
         z03,z04,pz03,pz04
    double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
         sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
    integer :: i,my_rank,ierr,j
    double precision, dimension(6) :: localmax
    double precision, dimension(14) :: localmax2, glmax
    double precision, dimension(58) :: tmplc,tmpgl
    double precision :: t0,lcrmax
    !        double precision :: alphax,alphay,alphaz
    integer :: npctmin,npctmax,npyhalf,nfile
    integer, parameter :: nbin = 400
    integer, dimension(2*nbin) :: tmpbin,glbin
    double precision :: f90,f95,f99,ex90,ex95,ex99,ex902,ex952,&
         ex992,ey90,ey95,ey99,ey902,ey952,ey992
    double precision,allocatable,dimension(:) :: epsiontmp
    double precision, dimension(2) :: tmpepslc,tmpep,Ealpha,Ebeta,&
         Egamma
    double precision :: xtmp,pxtmp,ytmp,pytmp,tmp1,tmp2,epsmylc,hyeps,&
         epsmxlc,hxeps,epsmx,epsmy
    integer :: nii,iitmp
    !effective  emittance
    real*8 :: effemtxlc,effemtx,effemtylc,effemty,gammamx,gammamy

    call starttime_Timer(t0)

    npctmin = 1
    npctmax = 1

    !machine Twiss parameters
    gammamx = (1.0d0+alphamx**2)/betamx
    gammamy = (1.0d0+alphamy**2)/betamy

    call MPI_COMM_RANK(mpicommwd,my_rank,ierr)
    !        call MPI_COMM_SIZE(mpicommwd,nproc,ierr)
    npyhalf = npy/2 

    den1 = 1.0/dble(nptot)
    x0lc = 0.0d0
    px0lc = 0.0d0
    y0lc = 0.0d0
    py0lc = 0.0d0
    z0lc = 0.0d0
    pz0lc = 0.0d0
    sqsum1local = 0.0d0
    sqsum2local = 0.0d0
    sqsum3local = 0.0d0
    sqsum4local = 0.0d0
    sqsum5local = 0.0d0
    sqsum6local = 0.0d0
    xpxlocal = 0.0d0
    ypylocal = 0.0d0
    zpzlocal = 0.0d0
    x0lc3 = 0.0d0
    x0lc4 = 0.0d0
    px0lc3 = 0.0d0
    px0lc4 = 0.0d0
    y0lc3 = 0.0d0
    y0lc4 = 0.0d0
    py0lc3 = 0.0d0
    py0lc4 = 0.0d0
    z0lc3 = 0.0d0
    z0lc4 = 0.0d0
    pz0lc3 = 0.0d0
    pz0lc4 = 0.0d0

    ! for cache optimization.
    if(innp.ne.0) then
       do i = 1, 6
          localmax(i) = abs(Pts1(i,1))
       enddo
       lcrmax = Pts1(1,1)**2+Pts1(3,1)**2
    else
       do i = 1, 6
          localmax(i) = 0.0d0
       enddo
       lcrmax = 0.0d0
    endif

    effemtxlc = 0.0d0
    effemtylc = 0.0d0
    do i = 1, innp
       x0lc = x0lc + Pts1(1,i)
       sqsum1local = sqsum1local + Pts1(1,i)*Pts1(1,i)
       x0lc3 = x0lc3 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)
       x0lc4 = x0lc4 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)*&
            Pts1(1,i)
       xpxlocal = xpxlocal + Pts1(1,i)*Pts1(2,i)
       px0lc = px0lc + Pts1(2,i)
       sqsum2local = sqsum2local + Pts1(2,i)*Pts1(2,i)
       px0lc3 = px0lc3 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)
       px0lc4 = px0lc4 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)*&
            Pts1(2,i)
       effemtxlc = effemtxlc + gammamx*Pts1(1,i)**2+ &
         2*alphamx*Pts1(1,i)*Pts1(2,i)+betamx*Pts1(2,i)**2
       y0lc = y0lc + Pts1(3,i)
       sqsum3local = sqsum3local + Pts1(3,i)*Pts1(3,i)
       y0lc3 = y0lc3 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)
       y0lc4 = y0lc4 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)*&
            Pts1(3,i)
       ypylocal = ypylocal + Pts1(3,i)*Pts1(4,i)
       py0lc = py0lc + Pts1(4,i)
       sqsum4local = sqsum4local + Pts1(4,i)*Pts1(4,i)
       py0lc3 = py0lc3 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)
       py0lc4 = py0lc4 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)*&
            Pts1(4,i)
       effemtylc = effemtylc + gammamy*Pts1(3,i)**2+ &
         2*alphamy*Pts1(3,i)*Pts1(4,i)+betamy*Pts1(4,i)**2
       z0lc = z0lc + Pts1(5,i)
       sqsum5local = sqsum5local + Pts1(5,i)*Pts1(5,i)
       z0lc3 = z0lc3 + abs(Pts1(5,i)*Pts1(5,i)*Pts1(5,i))
       z0lc4 = z0lc4 + Pts1(5,i)*Pts1(5,i)*Pts1(5,i)*&
            Pts1(5,i)

       zpzlocal = zpzlocal + Pts1(5,i)*Pts1(6,i)
       pz0lc = pz0lc + Pts1(6,i)
       sqsum6local = sqsum6local + Pts1(6,i)*Pts1(6,i)
       pz0lc3 = pz0lc3 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)
       pz0lc4 = pz0lc4 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)*&
            Pts1(6,i)
       do j = 1, 6
          if(localmax(j).lt.abs(Pts1(j,i))) then
             localmax(j) = abs(Pts1(j,i))
          endif
       enddo
       if(lcrmax.lt.(Pts1(1,i)**2+Pts1(3,i)**2)) then
          lcrmax = Pts1(1,i)**2 + Pts1(3,i)**2
       endif
    enddo

    localmax2 = 0.0d0
    tmplc = 0.0d0
    if(myidy.lt.npyhalf) then
       tmplc(1) = x0lc
       tmplc(2) = px0lc
       tmplc(3) = y0lc
       tmplc(4) = py0lc
       tmplc(5) = z0lc
       tmplc(6) = pz0lc
       tmplc(7) = sqsum1local
       tmplc(8) = sqsum2local
       tmplc(9) = sqsum3local
       tmplc(10) = sqsum4local
       tmplc(11) = sqsum5local
       tmplc(12) = sqsum6local
       tmplc(13) = xpxlocal
       tmplc(14) = ypylocal
       tmplc(15) = zpzlocal
       tmplc(16) = x0lc3
       tmplc(17) = x0lc4
       tmplc(18) = px0lc3
       tmplc(19) = px0lc4
       tmplc(20) = y0lc3
       tmplc(21) = y0lc4
       tmplc(22) = py0lc3
       tmplc(23) = py0lc4
       tmplc(24) = z0lc3
       tmplc(25) = z0lc4
       tmplc(26) = pz0lc3
       tmplc(27) = pz0lc4
       tmplc(55) = effemtxlc
       tmplc(56) = effemtylc
       localmax2(1) = localmax(1)
       localmax2(2) = localmax(2)
       localmax2(3) = localmax(3)
       localmax2(4) = localmax(4)
       localmax2(5) = localmax(5)
       localmax2(6) = localmax(6)
       localmax2(7) = lcrmax
    else
       tmplc(28) = x0lc
       tmplc(29) = px0lc
       tmplc(30) = y0lc
       tmplc(31) = py0lc
       tmplc(32) = z0lc
       tmplc(33) = pz0lc
       tmplc(34) = sqsum1local
       tmplc(35) = sqsum2local
       tmplc(36) = sqsum3local
       tmplc(37) = sqsum4local
       tmplc(38) = sqsum5local
       tmplc(39) = sqsum6local
       tmplc(40) = xpxlocal
       tmplc(41) = ypylocal
       tmplc(42) = zpzlocal
       tmplc(43) = x0lc3
       tmplc(44) = x0lc4
       tmplc(45) = px0lc3
       tmplc(46) = px0lc4
       tmplc(47) = y0lc3
       tmplc(48) = y0lc4
       tmplc(49) = py0lc3
       tmplc(50) = py0lc4
       tmplc(51) = z0lc3
       tmplc(52) = z0lc4
       tmplc(53) = pz0lc3
       tmplc(54) = pz0lc4
       localmax2(8) = localmax(1)
       localmax2(9) = localmax(2)
       localmax2(10) = localmax(3)
       localmax2(11) = localmax(4)
       localmax2(12) = localmax(5)
       localmax2(13) = localmax(6)
       localmax2(14) = lcrmax
       tmplc(57) = effemtxlc
       tmplc(58) = effemtylc
    endif

    !call MPI_REDUCE(tmplc,tmpgl,54,MPI_DOUBLE_PRECISION,&
    !                MPI_SUM,0,mpicommwd,ierr)
    !call MPI_ALLREDUCE(tmplc,tmpgl,54,MPI_DOUBLE_PRECISION,&
    !                MPI_SUM,mpicommwd,ierr)
    call MPI_ALLREDUCE(tmplc,tmpgl,58,MPI_DOUBLE_PRECISION,&
         MPI_SUM,comm2d,ierr)
    !call MPI_REDUCE(localmax2,glmax,14,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
    !                mpicommwd,ierr)
    !call MPI_ALLREDUCE(localmax2,glmax,14,MPI_DOUBLE_PRECISION,MPI_MAX,&
    !                mpicommwd,ierr)
    call MPI_ALLREDUCE(localmax2,glmax,14,MPI_DOUBLE_PRECISION,MPI_MAX,&
         comm2d,ierr)

    if(myidy.lt.npyhalf) then
       !//output beam 1
       x0 = tmpgl(1)*den1
       px0 = tmpgl(2)*den1
       y0 = tmpgl(3)*den1
       py0 = tmpgl(4)*den1
       z0 = tmpgl(5)*den1
       pz0 = tmpgl(6)*den1
       sqx = tmpgl(7)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(8)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(9)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(10)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(11)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(12)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(13)*den1 - x0*px0
       ypy = tmpgl(14)*den1 - y0*py0
       zpz = tmpgl(15)*den1 - z0*pz0
       cubx = tmpgl(16)*den1
       fthx = tmpgl(17)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0d0/3.0d0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(18)*den1
       fthpx = tmpgl(19)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0d0/3.0d0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(20)*den1
       fthy = tmpgl(21)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0d0/3.0d0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(22)*den1
       fthpy = tmpgl(23)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0d0/3.0d0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(24)*den1
       fthz = tmpgl(25)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0d0/3.0d0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(26)*den1
       fthpz = tmpgl(27)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0d0/3.0d0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0d0
       ypyfac = 0.0d0
       zpzfac = 0.0d0
       if(xrms.ne.0.0d0 .and. pxrms.ne.0.0d0)xpxfac=1.0d0/(xrms*pxrms)
       if(yrms.ne.0.0d0 .and. pyrms.ne.0.0d0)ypyfac=1.0d0/(yrms*pyrms)
       if(zrms.ne.0.0d0 .and. pzrms.ne.0.0d0)zpzfac=1.0d0/(zrms*pzrms)
       effemtx = tmpgl(55)*den1
       effemty = tmpgl(56)*den1

       if(myidx.eq.0 .and. myidy.eq.0) then
          nfile = 4+bid*20
          write(nfile,102)z,x0,xrms,px0,pxrms,-xpx/epx,epx,effemtx/2
          call flush(nfile)
          nfile = 5+bid*20
          write(nfile,102)z,y0,yrms,py0,pyrms,-ypy/epy,epy,effemty/2
          call flush(nfile)
          nfile = 6+bid*20
          write(nfile,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
          call flush(nfile)
          nfile = 7+bid*20
          write(nfile,102)z,glmax(1),glmax(2),glmax(3),glmax(4),glmax(5),&
               glmax(6),sqrt(glmax(7))
          call flush(nfile)
          nfile = 8+bid*20
          write(nfile,101)z,npctmin,npctmax,nptot
          call flush(nfile)
          nfile = 9+bid*20
          write(nfile,100)z,x03,px03,y03,py03,z03,pz03
          call flush(nfile)
          nfile = 10+bid*20
          write(nfile,100)z,x04,px04,y04,py04,z04,pz04
          call flush(nfile)
       endif
    else
       !//output beam 2
       x0 = tmpgl(28)*den1
       px0 = tmpgl(29)*den1
       y0 = tmpgl(30)*den1
       py0 = tmpgl(31)*den1
       z0 = tmpgl(32)*den1
       pz0 = tmpgl(33)*den1
       sqx = tmpgl(34)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(35)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(36)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(37)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(38)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(39)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(40)*den1 - x0*px0
       ypy = tmpgl(41)*den1 - y0*py0
       zpz = tmpgl(42)*den1 - z0*pz0
       cubx = tmpgl(43)*den1
       fthx = tmpgl(44)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0d0/3.0d0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(45)*den1
       fthpx = tmpgl(46)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0d0/3.0d0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(47)*den1
       fthy = tmpgl(48)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0d0/3.0d0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(49)*den1
       fthpy = tmpgl(50)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0d0/3.0d0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(51)*den1
       fthz = tmpgl(52)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0d0/3.0d0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(53)*den1
       fthpz = tmpgl(54)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0d0/3.0d0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0d0
       ypyfac = 0.0d0
       zpzfac = 0.0d0
       if(xrms.ne.0.0d0 .and. pxrms.ne.0.0d0)xpxfac=1.0d0/(xrms*pxrms)
       if(yrms.ne.0.0d0 .and. pyrms.ne.0.0d0)ypyfac=1.0d0/(yrms*pyrms)
       if(zrms.ne.0.0d0 .and. pzrms.ne.0.0d0)zpzfac=1.0d0/(zrms*pzrms)
       effemtx = tmpgl(57)*den1
       effemty = tmpgl(58)*den1

       if(myidx.eq.0 .and. myidy.eq.npyhalf) then
          nfile = 14+bid*20
          write(nfile,102)z,x0,xrms,px0,pxrms,-xpx/epx,epx,effemtx/2
          call flush(nfile)
          nfile = 15+bid*20
          write(nfile,102)z,y0,yrms,py0,pyrms,-ypy/epy,epy,effemty/2
          call flush(nfile)
          nfile = 16+bid*20
          write(nfile,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
          call flush(nfile)
          nfile = 17+bid*20
          write(nfile,102)z,glmax(8),glmax(9),glmax(10),glmax(11),glmax(12),&
               glmax(13),sqrt(glmax(14))
          call flush(nfile)
          nfile = 18+bid*20
          write(nfile,101)z,npctmin,npctmax,nptot
          call flush(nfile)
          nfile = 19+bid*20
          write(nfile,100)z,x03,px03,y03,py03,z03,pz03
          call flush(nfile)
          nfile = 20+bid*20
          write(nfile,100)z,x04,px04,y04,py04,z04,pz04
          call flush(nfile)
       endif
    endif

100 format(7(1x,es13.6))
101 format(1x,es13.6,3I10)
102 format(8(1x,es13.6))

    t_diag = t_diag + elapsedtime_Timer(t0)

  end subroutine diagnostic2GMBeic_Output

  subroutine diagnostic2GMB_Output(z,Pts1,innp,nptot,&
       myidx,myidy,npx,npy,commrow,commcol,comm2d,bid,gambet)
    implicit none
    include 'mpif.h'
    double precision, intent(in) :: z
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot,myidx,myidy,npx,npy,commrow,&
         commcol,comm2d,bid
    double precision:: den1,sqsum1,sqsum2,sqsum3,sqsum4,&
         epsx2,epsy2,gambet
    double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
         xpxfac,ypyfac
    double precision :: sqsum1local,sqsum2local,sqsum3local,sqsum4local
    double precision :: xpxlocal,ypylocal,zpzlocal
    double precision :: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
    double precision :: sqsum5local,sqsum6local
    double precision :: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
         pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
         z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
         z03,z04,pz03,pz04
    double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
         sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
    integer :: i,my_rank,ierr,j
    double precision, dimension(6) :: localmax
    double precision, dimension(14) :: localmax2, glmax
    double precision, dimension(54) :: tmplc,tmpgl
    double precision :: t0,lcrmax
    !        double precision :: alphax,alphay,alphaz
    integer :: npctmin,npctmax,npyhalf,nfile
    integer, parameter :: nbin = 400
    integer, dimension(2*nbin) :: tmpbin,glbin
    double precision :: f90,f95,f99,ex90,ex95,ex99,ex902,ex952,&
         ex992,ey90,ey95,ey99,ey902,ey952,ey992
    double precision,allocatable,dimension(:) :: epsiontmp
    double precision, dimension(2) :: tmpepslc,tmpep,Ealpha,Ebeta,&
         Egamma
    double precision :: xtmp,pxtmp,ytmp,pytmp,tmp1,tmp2,epsmylc,hyeps,&
         epsmxlc,hxeps,epsmx,epsmy
    integer :: nii,iitmp

    call starttime_Timer(t0)

    npctmin = 1
    npctmax = 1

    call MPI_COMM_RANK(mpicommwd,my_rank,ierr)
    !        call MPI_COMM_SIZE(mpicommwd,nproc,ierr)
    npyhalf = npy/2 

    den1 = 1.0/dble(nptot)
    x0lc = 0.0d0
    px0lc = 0.0d0
    y0lc = 0.0d0
    py0lc = 0.0d0
    z0lc = 0.0d0
    pz0lc = 0.0d0
    sqsum1local = 0.0d0
    sqsum2local = 0.0d0
    sqsum3local = 0.0d0
    sqsum4local = 0.0d0
    sqsum5local = 0.0d0
    sqsum6local = 0.0d0
    xpxlocal = 0.0d0
    ypylocal = 0.0d0
    zpzlocal = 0.0d0
    x0lc3 = 0.0d0
    x0lc4 = 0.0d0
    px0lc3 = 0.0d0
    px0lc4 = 0.0d0
    y0lc3 = 0.0d0
    y0lc4 = 0.0d0
    py0lc3 = 0.0d0
    py0lc4 = 0.0d0
    z0lc3 = 0.0d0
    z0lc4 = 0.0d0
    pz0lc3 = 0.0d0
    pz0lc4 = 0.0d0

    ! for cache optimization.
    if(innp.ne.0) then
       do i = 1, 6
          localmax(i) = abs(Pts1(i,1))
       enddo
       lcrmax = Pts1(1,1)**2+Pts1(3,1)**2
    else
       do i = 1, 6
          localmax(i) = 0.0d0
       enddo
       lcrmax = 0.0d0
    endif
    do i = 1, innp
       x0lc = x0lc + Pts1(1,i)
       sqsum1local = sqsum1local + Pts1(1,i)*Pts1(1,i)
       x0lc3 = x0lc3 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)
       x0lc4 = x0lc4 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)*&
            Pts1(1,i)
       xpxlocal = xpxlocal + Pts1(1,i)*Pts1(2,i)
       px0lc = px0lc + Pts1(2,i)
       sqsum2local = sqsum2local + Pts1(2,i)*Pts1(2,i)
       px0lc3 = px0lc3 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)
       px0lc4 = px0lc4 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)*&
            Pts1(2,i)
       y0lc = y0lc + Pts1(3,i)
       sqsum3local = sqsum3local + Pts1(3,i)*Pts1(3,i)
       y0lc3 = y0lc3 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)
       y0lc4 = y0lc4 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)*&
            Pts1(3,i)
       ypylocal = ypylocal + Pts1(3,i)*Pts1(4,i)
       py0lc = py0lc + Pts1(4,i)
       sqsum4local = sqsum4local + Pts1(4,i)*Pts1(4,i)
       py0lc3 = py0lc3 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)
       py0lc4 = py0lc4 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)*&
            Pts1(4,i)
       z0lc = z0lc + Pts1(5,i)
       sqsum5local = sqsum5local + Pts1(5,i)*Pts1(5,i)
       z0lc3 = z0lc3 + abs(Pts1(5,i)*Pts1(5,i)*Pts1(5,i))
       z0lc4 = z0lc4 + Pts1(5,i)*Pts1(5,i)*Pts1(5,i)*&
            Pts1(5,i)

       zpzlocal = zpzlocal + Pts1(5,i)*Pts1(6,i)
       pz0lc = pz0lc + Pts1(6,i)
       sqsum6local = sqsum6local + Pts1(6,i)*Pts1(6,i)
       pz0lc3 = pz0lc3 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)
       pz0lc4 = pz0lc4 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)*&
            Pts1(6,i)
       do j = 1, 6
          if(localmax(j).lt.abs(Pts1(j,i))) then
             localmax(j) = abs(Pts1(j,i))
          endif
       enddo
       if(lcrmax.lt.(Pts1(1,i)**2+Pts1(3,i)**2)) then
          lcrmax = Pts1(1,i)**2 + Pts1(3,i)**2
       endif
    enddo

    localmax2 = 0.0d0
    tmplc = 0.0d0
    if(myidy.lt.npyhalf) then
       tmplc(1) = x0lc
       tmplc(2) = px0lc
       tmplc(3) = y0lc
       tmplc(4) = py0lc
       tmplc(5) = z0lc
       tmplc(6) = pz0lc
       tmplc(7) = sqsum1local
       tmplc(8) = sqsum2local
       tmplc(9) = sqsum3local
       tmplc(10) = sqsum4local
       tmplc(11) = sqsum5local
       tmplc(12) = sqsum6local
       tmplc(13) = xpxlocal
       tmplc(14) = ypylocal
       tmplc(15) = zpzlocal
       tmplc(16) = x0lc3
       tmplc(17) = x0lc4
       tmplc(18) = px0lc3
       tmplc(19) = px0lc4
       tmplc(20) = y0lc3
       tmplc(21) = y0lc4
       tmplc(22) = py0lc3
       tmplc(23) = py0lc4
       tmplc(24) = z0lc3
       tmplc(25) = z0lc4
       tmplc(26) = pz0lc3
       tmplc(27) = pz0lc4
       localmax2(1) = localmax(1)
       localmax2(2) = localmax(2)
       localmax2(3) = localmax(3)
       localmax2(4) = localmax(4)
       localmax2(5) = localmax(5)
       localmax2(6) = localmax(6)
       localmax2(7) = lcrmax
    else
       tmplc(28) = x0lc
       tmplc(29) = px0lc
       tmplc(30) = y0lc
       tmplc(31) = py0lc
       tmplc(32) = z0lc
       tmplc(33) = pz0lc
       tmplc(34) = sqsum1local
       tmplc(35) = sqsum2local
       tmplc(36) = sqsum3local
       tmplc(37) = sqsum4local
       tmplc(38) = sqsum5local
       tmplc(39) = sqsum6local
       tmplc(40) = xpxlocal
       tmplc(41) = ypylocal
       tmplc(42) = zpzlocal
       tmplc(43) = x0lc3
       tmplc(44) = x0lc4
       tmplc(45) = px0lc3
       tmplc(46) = px0lc4
       tmplc(47) = y0lc3
       tmplc(48) = y0lc4
       tmplc(49) = py0lc3
       tmplc(50) = py0lc4
       tmplc(51) = z0lc3
       tmplc(52) = z0lc4
       tmplc(53) = pz0lc3
       tmplc(54) = pz0lc4
       localmax2(8) = localmax(1)
       localmax2(9) = localmax(2)
       localmax2(10) = localmax(3)
       localmax2(11) = localmax(4)
       localmax2(12) = localmax(5)
       localmax2(13) = localmax(6)
       localmax2(14) = lcrmax
    endif

    !call MPI_REDUCE(tmplc,tmpgl,54,MPI_DOUBLE_PRECISION,&
    !                MPI_SUM,0,mpicommwd,ierr)
    !call MPI_ALLREDUCE(tmplc,tmpgl,54,MPI_DOUBLE_PRECISION,&
    !                MPI_SUM,mpicommwd,ierr)
    call MPI_ALLREDUCE(tmplc,tmpgl,54,MPI_DOUBLE_PRECISION,&
         MPI_SUM,comm2d,ierr)
    !call MPI_REDUCE(localmax2,glmax,14,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
    !                mpicommwd,ierr)
    !call MPI_ALLREDUCE(localmax2,glmax,14,MPI_DOUBLE_PRECISION,MPI_MAX,&
    !                mpicommwd,ierr)
    call MPI_ALLREDUCE(localmax2,glmax,14,MPI_DOUBLE_PRECISION,MPI_MAX,&
         comm2d,ierr)

    if(myidy.lt.npyhalf) then
       !//output beam 1
       x0 = tmpgl(1)*den1
       px0 = tmpgl(2)*den1
       y0 = tmpgl(3)*den1
       py0 = tmpgl(4)*den1
       z0 = tmpgl(5)*den1
       pz0 = tmpgl(6)*den1
       sqx = tmpgl(7)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(8)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(9)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(10)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(11)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(12)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(13)*den1 - x0*px0
       ypy = tmpgl(14)*den1 - y0*py0
       zpz = tmpgl(15)*den1 - z0*pz0
       cubx = tmpgl(16)*den1
       fthx = tmpgl(17)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0d0/3.0d0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(18)*den1
       fthpx = tmpgl(19)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0d0/3.0d0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(20)*den1
       fthy = tmpgl(21)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0d0/3.0d0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(22)*den1
       fthpy = tmpgl(23)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0d0/3.0d0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(24)*den1
       fthz = tmpgl(25)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0d0/3.0d0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(26)*den1
       fthpz = tmpgl(27)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0d0/3.0d0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0d0
       ypyfac = 0.0d0
       zpzfac = 0.0d0
       if(xrms.ne.0.0d0 .and. pxrms.ne.0.0d0)xpxfac=1.0d0/(xrms*pxrms)
       if(yrms.ne.0.0d0 .and. pyrms.ne.0.0d0)ypyfac=1.0d0/(yrms*pyrms)
       if(zrms.ne.0.0d0 .and. pzrms.ne.0.0d0)zpzfac=1.0d0/(zrms*pzrms)

       if(myidx.eq.0 .and. myidy.eq.0) then
          nfile = 4+bid*20
          write(nfile,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
          call flush(nfile)
          nfile = 5+bid*20
          write(nfile,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
          call flush(nfile)
          nfile = 6+bid*20
          write(nfile,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
          call flush(nfile)
          nfile = 7+bid*20
          write(nfile,102)z,glmax(1),glmax(2),glmax(3),glmax(4),glmax(5),&
               glmax(6),sqrt(glmax(7))
          call flush(nfile)
          nfile = 8+bid*20
          write(nfile,101)z,npctmin,npctmax,nptot
          call flush(nfile)
          nfile = 9+bid*20
          write(nfile,100)z,x03,px03,y03,py03,z03,pz03
          call flush(nfile)
          nfile = 10+bid*20
          write(nfile,100)z,x04,px04,y04,py04,z04,pz04
          call flush(nfile)
       endif
    else
       !//output beam 2
       x0 = tmpgl(28)*den1
       px0 = tmpgl(29)*den1
       y0 = tmpgl(30)*den1
       py0 = tmpgl(31)*den1
       z0 = tmpgl(32)*den1
       pz0 = tmpgl(33)*den1
       sqx = tmpgl(34)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(35)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(36)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(37)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(38)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(39)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(40)*den1 - x0*px0
       ypy = tmpgl(41)*den1 - y0*py0
       zpz = tmpgl(42)*den1 - z0*pz0
       cubx = tmpgl(43)*den1
       fthx = tmpgl(44)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0d0/3.0d0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(45)*den1
       fthpx = tmpgl(46)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0d0/3.0d0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(47)*den1
       fthy = tmpgl(48)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0d0/3.0d0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(49)*den1
       fthpy = tmpgl(50)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0d0/3.0d0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(51)*den1
       fthz = tmpgl(52)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0d0/3.0d0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(53)*den1
       fthpz = tmpgl(54)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0d0/3.0d0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0d0
       ypyfac = 0.0d0
       zpzfac = 0.0d0
       if(xrms.ne.0.0d0 .and. pxrms.ne.0.0d0)xpxfac=1.0d0/(xrms*pxrms)
       if(yrms.ne.0.0d0 .and. pyrms.ne.0.0d0)ypyfac=1.0d0/(yrms*pyrms)
       if(zrms.ne.0.0d0 .and. pzrms.ne.0.0d0)zpzfac=1.0d0/(zrms*pzrms)

       if(myidx.eq.0 .and. myidy.eq.npyhalf) then
          nfile = 14+bid*20
          write(nfile,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
          call flush(nfile)
          nfile = 15+bid*20
          write(nfile,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
          call flush(nfile)
          nfile = 16+bid*20
          write(nfile,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
          call flush(nfile)
          nfile = 17+bid*20
          write(nfile,102)z,glmax(8),glmax(9),glmax(10),glmax(11),glmax(12),&
               glmax(13),sqrt(glmax(14))
          call flush(nfile)
          nfile = 18+bid*20
          write(nfile,101)z,npctmin,npctmax,nptot
          call flush(nfile)
          nfile = 19+bid*20
          write(nfile,100)z,x03,px03,y03,py03,z03,pz03
          call flush(nfile)
          nfile = 20+bid*20
          write(nfile,100)z,x04,px04,y04,py04,z04,pz04
          call flush(nfile)
       endif
    endif


    !goto 999

    Ealpha(1) = -xpx/epx
    Ealpha(2) = -ypy/epy
    Ebeta(1) = xrms*xrms*gambet/epx
    Ebeta(2) = yrms*yrms*gambet/epy
    Egamma(:) = (1.0+Ealpha(:)*Ealpha(:))/Ebeta(:)

    allocate(epsiontmp(innp))
    epsmxlc = -1.0e10
    do i = 1, innp
       xtmp = Pts1(1,i) - x0
       tmp1 = Pts1(2,i)
       pxtmp = (tmp1 - px0)/gambet
       epsiontmp(i)=Egamma(1)*xtmp*xtmp+2*Ealpha(1)*xtmp*pxtmp+&
            Ebeta(1)*pxtmp*pxtmp
       if(epsmxlc.le.epsiontmp(i)) epsmxlc = epsiontmp(i)
    enddo

    tmpepslc = 0.0d0
    if(myidy.lt.npyhalf) then
       tmpepslc(1) = epsmxlc
    else
       tmpepslc(2) = epsmxlc
    endif

    call MPI_ALLREDUCE(tmpepslc,tmpep,2,MPI_DOUBLE_PRECISION,&
         MPI_MAX,mpicommwd,ierr)
    if(myidy.lt.npyhalf) then
       epsmx = tmpep(1) 
    else
       epsmx = tmpep(2) 
    endif

    hxeps = epsmx*1.0001d0/nbin
    tmpbin = 0


    if(myidy.lt.npyhalf) then
       do i = 1, innp
          iitmp = int(epsiontmp(i)/hxeps)
          nii = iitmp+1
          tmpbin(nii) = tmpbin(nii) + 1
       enddo
    else
       do i = 1, innp
          iitmp = int(epsiontmp(i)/hxeps)
          nii = iitmp+1
          tmpbin(nii+nbin) = tmpbin(nii+nbin) + 1
       enddo
    endif

    if(myidy.eq.0) then
       do i = 1, 2*nbin
       enddo
    endif

    call MPI_REDUCE(tmpbin,glbin,2*nbin,MPI_INTEGER,&
         MPI_SUM,0,mpicommwd,ierr)

    !        if(my_rank.eq.0) then
    !          do i = 1, 2*nbin
    !          enddo
    !        endif

    f90 = 0.999d0*nptot
    f95 = 0.9999d0*nptot
    f99 = 0.99999d0*nptot
    if(my_rank.eq.0) then

       hxeps = tmpep(1)*1.0001d0/nbin
       do i = 2, nbin
          glbin(i) = glbin(i) + glbin(i-1)
       enddo
       do i = 1, nbin
          if(glbin(i).gt.f90) then
             ex90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                  hxeps*(i-1)
             exit
          endif
       enddo
       do i = 1, nbin
          if(glbin(i).gt.f95) then
             ex95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                  hxeps*(i-1)
             exit
          endif
       enddo
       do i =1, nbin
          if(glbin(i).gt.f99) then
             ex99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                  hxeps*(i-1)
             exit
          endif
       enddo
       !          ex90 = ex90*gambet
       !          ex95 = ex95*gambet
       !          ex99 = ex99*gambet

       hxeps = tmpep(2)*1.0001/nbin
       do i = nbin+2, 2*nbin
          glbin(i) = glbin(i) + glbin(i-1)
       enddo
       do i = nbin+1, 2*nbin
          if(glbin(i).gt.f90) then
             ex902 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                  hxeps*(i-nbin-1)
             exit
          endif
       enddo
       do i = nbin+1, 2*nbin
          if(glbin(i).gt.f95) then
             ex952 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                  hxeps*(i-nbin-1)
             exit
          endif
       enddo
       do i =nbin+1, 2*nbin
          if(glbin(i).gt.f99) then
             ex992 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                  hxeps*(i-nbin-1)
             exit
          endif
       enddo
       !          ex902 = ex902*gambet
       !          ex952 = ex952*gambet
       !          ex992 = ex992*gambet

       nfile = 11+bid*20
       write(nfile,100)z,ex90,ex95,ex99,ex902,ex952,ex992
       call flush(nfile)
    endif

    epsmylc = -1.0e10
    do i = 1, innp
       ytmp = Pts1(3,i) - y0
       tmp2 = Pts1(4,i)
       pytmp = (tmp2 - py0)/gambet
       epsiontmp(i)=Egamma(2)*ytmp*ytmp+2*Ealpha(2)*ytmp*pytmp+&
            Ebeta(2)*pytmp*pytmp
       if(epsmylc.le.epsiontmp(i)) epsmylc = epsiontmp(i)
    enddo

    tmpepslc = 0.0d0
    if(myidy.lt.npyhalf) then
       tmpepslc(1) = epsmylc
    else
       tmpepslc(2) = epsmylc
    endif
    call MPI_ALLREDUCE(tmpepslc,tmpep,2,MPI_DOUBLE_PRECISION,&
         MPI_MAX,mpicommwd,ierr)
    if(myidy.lt.npyhalf) then
       epsmy = tmpep(1)
    else
       epsmy = tmpep(2)
    endif

    hyeps = epsmy*1.0001d0/nbin
    tmpbin = 0
    if(myidy.lt.npyhalf) then
       do i = 1, innp
          iitmp = int(epsiontmp(i)/hyeps)
          nii = iitmp+1
          tmpbin(nii) = tmpbin(nii) + 1
       enddo
    else
       do i = 1, innp
          iitmp = int(epsiontmp(i)/hyeps)
          nii = iitmp+1
          tmpbin(nii+nbin) = tmpbin(nii+nbin) + 1
       enddo
    endif
    call MPI_REDUCE(tmpbin,glbin,2*nbin,MPI_INTEGER,&
         MPI_SUM,0,mpicommwd,ierr)

    if(my_rank.eq.0) then
       hyeps = tmpep(1)*1.0001d0/nbin
       do i = 2, nbin
          glbin(i) = glbin(i) + glbin(i-1)
       enddo
       do i = 1, nbin
          if(glbin(i).gt.f90) then
             ey90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                  hyeps*(i-1)
             exit
          endif
       enddo
       do i = 1, nbin
          if(glbin(i).gt.f95) then
             ey95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                  hyeps*(i-1)
             exit
          endif
       enddo
       do i = 1, nbin
          if(glbin(i).gt.f99) then
             ey99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                  hyeps*(i-1)
             exit
          endif
       enddo
       !          ey90 = ey90*gambet
       !          ey95 = ey95*gambet
       !          ey99 = ey99*gambet

       hyeps = tmpep(2)*1.0001d0/nbin
       do i = nbin+2, 2*nbin
          glbin(i) = glbin(i) + glbin(i-1)
       enddo
       do i = nbin+1, 2*nbin
          if(glbin(i).gt.f90) then
             ey902 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                  hyeps*(i-nbin-1)
             exit
          endif
       enddo
       do i = nbin+1, 2*nbin
          if(glbin(i).gt.f95) then
             ey952 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                  hyeps*(i-nbin-1)
             exit
          endif
       enddo
       do i =nbin+1, 2*nbin
          if(glbin(i).gt.f99) then
             ey992 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                  hyeps*(i-nbin-1)
             exit
          endif
       enddo
       !          ey902 = ey902*gambet
       !          ey952 = ey952*gambet
       !          ey992 = ey992*gambet

       nfile = 21+bid*20
       write(nfile,100)z,ey90,ey95,ey99,ey902,ey952,ey992
       call flush(nfile)
    endif
    deallocate(epsiontmp)

!999 continue

100 format(7(1x,es13.6))
101 format(1x,es13.6,3I10)
102 format(8(1x,es13.6))

    t_diag = t_diag + elapsedtime_Timer(t0)

  end subroutine diagnostic2GMB_Output



  ! calculate <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
  ! <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance using 2 group PEs.
  subroutine diagnostic2Gbak_Output(z,Pts1,innp,nptot,center1,center2,&
       sigma,myidx,myidy,npx,npy,commrow,commcol,comm2d)
    implicit none
    include 'mpif.h'
    double precision, intent(in) :: z
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot,myidx,myidy,npx,npy,commrow,&
         commcol,comm2d
    double precision, dimension(3), intent(out) :: center1,center2
    double precision, dimension(6), intent(out) :: sigma
    double precision:: den1,sqsum1,sqsum2,sqsum3,sqsum4,&
         epsx2,epsy2
    double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
         xpxfac,ypyfac
    double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
    double precision:: xpxlocal,ypylocal,zpzlocal
    double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
    double precision:: sqsum5local,sqsum6local
    double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
         pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
         z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
         z03,z04,pz03,pz04
    double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
         sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
    integer :: i,my_rank,ierr,j
    double precision, dimension(6) :: localmax
    double precision, dimension(14) :: localmax2, glmax
    double precision, dimension(54) :: tmplc,tmpgl
    double precision, dimension(6) :: sigma1,sigma2
    double precision :: t0,lcrmax
    !        double precision :: alphax,alphay,alphaz
    integer :: npctmin,npctmax,npyhalf

    call starttime_Timer(t0)

    npctmin = 1
    npctmax = 1

    call MPI_COMM_RANK(mpicommwd,my_rank,ierr)
    !        call MPI_COMM_SIZE(mpicommwd,nproc,ierr)
    npyhalf = npy/2 

    den1 = 1.0/dble(nptot)
    x0lc = 0.0
    px0lc = 0.0
    y0lc = 0.0
    py0lc = 0.0
    z0lc = 0.0
    pz0lc = 0.0
    sqsum1local = 0.0
    sqsum2local = 0.0
    sqsum3local = 0.0
    sqsum4local = 0.0
    sqsum5local = 0.0
    sqsum6local = 0.0
    xpxlocal = 0.0
    ypylocal = 0.0
    zpzlocal = 0.0
    x0lc3 = 0.0
    x0lc4 = 0.0
    px0lc3 = 0.0
    px0lc4 = 0.0
    y0lc3 = 0.0
    y0lc4 = 0.0
    py0lc3 = 0.0
    py0lc4 = 0.0
    z0lc3 = 0.0
    z0lc4 = 0.0
    pz0lc3 = 0.0
    pz0lc4 = 0.0

    ! for cache optimization.
    if(innp.ne.0) then
       do i = 1, 6
          localmax(i) = abs(Pts1(i,1))
       enddo
       lcrmax = Pts1(1,1)**2+Pts1(3,1)**2
    else
       do i = 1, 6
          localmax(i) = 0.0
       enddo
       lcrmax = 0.0
    endif
    do i = 1, innp
       x0lc = x0lc + Pts1(1,i)
       sqsum1local = sqsum1local + Pts1(1,i)*Pts1(1,i)
       x0lc3 = x0lc3 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)
       x0lc4 = x0lc4 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)*&
            Pts1(1,i)
       xpxlocal = xpxlocal + Pts1(1,i)*Pts1(2,i)
       px0lc = px0lc + Pts1(2,i)
       sqsum2local = sqsum2local + Pts1(2,i)*Pts1(2,i)
       px0lc3 = px0lc3 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)
       px0lc4 = px0lc4 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)*&
            Pts1(2,i)
       y0lc = y0lc + Pts1(3,i)
       sqsum3local = sqsum3local + Pts1(3,i)*Pts1(3,i)
       y0lc3 = y0lc3 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)
       y0lc4 = y0lc4 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)*&
            Pts1(3,i)
       ypylocal = ypylocal + Pts1(3,i)*Pts1(4,i)
       py0lc = py0lc + Pts1(4,i)
       sqsum4local = sqsum4local + Pts1(4,i)*Pts1(4,i)
       py0lc3 = py0lc3 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)
       py0lc4 = py0lc4 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)*&
            Pts1(4,i)
       z0lc = z0lc + Pts1(5,i)
       sqsum5local = sqsum5local + Pts1(5,i)*Pts1(5,i)
       z0lc3 = z0lc3 + abs(Pts1(5,i)*Pts1(5,i)*Pts1(5,i))
       z0lc4 = z0lc4 + Pts1(5,i)*Pts1(5,i)*Pts1(5,i)*&
            Pts1(5,i)

       zpzlocal = zpzlocal + Pts1(5,i)*Pts1(6,i)
       pz0lc = pz0lc + Pts1(6,i)
       sqsum6local = sqsum6local + Pts1(6,i)*Pts1(6,i)
       pz0lc3 = pz0lc3 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)
       pz0lc4 = pz0lc4 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)*&
            Pts1(6,i)
       do j = 1, 6
          if(localmax(j).lt.abs(Pts1(j,i))) then
             localmax(j) = abs(Pts1(j,i))
          endif
       enddo
       if(lcrmax.lt.(Pts1(1,i)**2+Pts1(3,i)**2)) then
          lcrmax = Pts1(1,i)**2 + Pts1(3,i)**2
       endif
    enddo

    localmax2 = 0.0
    tmplc = 0.0
    if(myidy.lt.npyhalf) then
       tmplc(1) = x0lc
       tmplc(2) = px0lc
       tmplc(3) = y0lc
       tmplc(4) = py0lc
       tmplc(5) = z0lc
       tmplc(6) = pz0lc
       tmplc(7) = sqsum1local
       tmplc(8) = sqsum2local
       tmplc(9) = sqsum3local
       tmplc(10) = sqsum4local
       tmplc(11) = sqsum5local
       tmplc(12) = sqsum6local
       tmplc(13) = xpxlocal
       tmplc(14) = ypylocal
       tmplc(15) = zpzlocal
       tmplc(16) = x0lc3
       tmplc(17) = x0lc4
       tmplc(18) = px0lc3
       tmplc(19) = px0lc4
       tmplc(20) = y0lc3
       tmplc(21) = y0lc4
       tmplc(22) = py0lc3
       tmplc(23) = py0lc4
       tmplc(24) = z0lc3
       tmplc(25) = z0lc4
       tmplc(26) = pz0lc3
       tmplc(27) = pz0lc4
       localmax2(1) = localmax(1)
       localmax2(2) = localmax(2)
       localmax2(3) = localmax(3)
       localmax2(4) = localmax(4)
       localmax2(5) = localmax(5)
       localmax2(6) = localmax(6)
       localmax2(7) = lcrmax
    else
       tmplc(28) = x0lc
       tmplc(29) = px0lc
       tmplc(30) = y0lc
       tmplc(31) = py0lc
       tmplc(32) = z0lc
       tmplc(33) = pz0lc
       tmplc(34) = sqsum1local
       tmplc(35) = sqsum2local
       tmplc(36) = sqsum3local
       tmplc(37) = sqsum4local
       tmplc(38) = sqsum5local
       tmplc(39) = sqsum6local
       tmplc(40) = xpxlocal
       tmplc(41) = ypylocal
       tmplc(42) = zpzlocal
       tmplc(43) = x0lc3
       tmplc(44) = x0lc4
       tmplc(45) = px0lc3
       tmplc(46) = px0lc4
       tmplc(47) = y0lc3
       tmplc(48) = y0lc4
       tmplc(49) = py0lc3
       tmplc(50) = py0lc4
       tmplc(51) = z0lc3
       tmplc(52) = z0lc4
       tmplc(53) = pz0lc3
       tmplc(54) = pz0lc4
       localmax2(8) = localmax(1)
       localmax2(9) = localmax(2)
       localmax2(10) = localmax(3)
       localmax2(11) = localmax(4)
       localmax2(12) = localmax(5)
       localmax2(13) = localmax(6)
       localmax2(14) = lcrmax
    endif

    call MPI_REDUCE(tmplc,tmpgl,54,MPI_DOUBLE_PRECISION,&
         MPI_SUM,0,mpicommwd,ierr)
    call MPI_REDUCE(localmax2,glmax,14,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
         mpicommwd,ierr)

    if(my_rank.eq.0) then
       !//output beam 1
       x0 = tmpgl(1)*den1
       px0 = tmpgl(2)*den1
       y0 = tmpgl(3)*den1
       py0 = tmpgl(4)*den1
       z0 = tmpgl(5)*den1
       pz0 = tmpgl(6)*den1
       sqx = tmpgl(7)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(8)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(9)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(10)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(11)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(12)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(13)*den1 - x0*px0
       ypy = tmpgl(14)*den1 - y0*py0
       zpz = tmpgl(15)*den1 - z0*pz0
       cubx = tmpgl(16)*den1
       fthx = tmpgl(17)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(18)*den1
       fthpx = tmpgl(19)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(20)*den1
       fthy = tmpgl(21)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(22)*den1
       fthpy = tmpgl(23)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(24)*den1
       fthz = tmpgl(25)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(26)*den1
       fthpz = tmpgl(27)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0
       ypyfac = 0.0
       zpzfac = 0.0
       if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
       if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
       if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
       write(24,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
       write(25,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
       write(26,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
       write(27,102)z,glmax(1),glmax(2),glmax(3),glmax(4),glmax(5),&
            glmax(6),sqrt(glmax(7))
       write(28,101)z,npctmin,npctmax,nptot
       write(29,100)z,x03,px03,y03,py03,z03,pz03
       write(30,100)z,x04,px04,y04,py04,z04,pz04

       call flush(24)
       call flush(25)
       call flush(26)
       call flush(27)
       call flush(28)
       call flush(29)
       call flush(30)
       center1(1) = x0
       center1(2) = y0
       center1(3) = z0
       sigma1(1) = xrms
       sigma1(2) = pxrms
       sigma1(3) = yrms
       sigma1(4) = pyrms
       sigma1(5) = zrms
       sigma1(6) = pzrms

       !//output beam 2
       x0 = tmpgl(28)*den1
       px0 = tmpgl(29)*den1
       y0 = tmpgl(30)*den1
       py0 = tmpgl(31)*den1
       z0 = tmpgl(32)*den1
       pz0 = tmpgl(33)*den1
       sqx = tmpgl(34)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(35)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(36)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(37)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(38)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(39)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(40)*den1 - x0*px0
       ypy = tmpgl(41)*den1 - y0*py0
       zpz = tmpgl(42)*den1 - z0*pz0
       cubx = tmpgl(43)*den1
       fthx = tmpgl(44)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(45)*den1
       fthpx = tmpgl(46)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(47)*den1
       fthy = tmpgl(48)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(49)*den1
       fthpy = tmpgl(50)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(51)*den1
       fthz = tmpgl(52)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(53)*den1
       fthpz = tmpgl(54)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0
       ypyfac = 0.0
       zpzfac = 0.0
       if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
       if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
       if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
       write(34,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
       write(35,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
       write(36,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
       write(37,102)z,glmax(8),glmax(9),glmax(10),glmax(11),glmax(12),&
            glmax(13),sqrt(glmax(14))
       write(38,101)z,npctmin,npctmax,nptot
       write(39,100)z,x03,px03,y03,py03,z03,pz03
       write(40,100)z,x04,px04,y04,py04,z04,pz04

       call flush(34)
       call flush(35)
       call flush(36)
       call flush(37)
       call flush(38)
       call flush(39)
       call flush(40)
       center2(1) = x0
       center2(2) = y0
       center2(3) = z0
       sigma2(1) = xrms
       sigma2(2) = pxrms
       sigma2(3) = yrms
       sigma2(4) = pyrms
       sigma2(5) = zrms
       sigma2(6) = pzrms
    endif

    call MPI_BCAST(center1,3,MPI_DOUBLE_PRECISION,0,mpicommwd,&
         ierr)
    call MPI_BCAST(sigma1,6,MPI_DOUBLE_PRECISION,0,mpicommwd,&
         ierr)
    call MPI_BCAST(center2,3,MPI_DOUBLE_PRECISION,0,mpicommwd,&
         ierr)
    call MPI_BCAST(sigma2,6,MPI_DOUBLE_PRECISION,0,mpicommwd,&
         ierr)

    if(myidy.lt.npyhalf) then
       sigma = sigma1
    else
       sigma = sigma2
    endif

100 format(7(1x,es13.6))
101 format(1x,es13.6,3I10)
102 format(8(1x,es13.6))

    t_diag = t_diag + elapsedtime_Timer(t0)

  end subroutine diagnostic2Gbak_Output



  subroutine diagnostic2Gold_Output(z,Pts1,innp,nptot,center1,center2,&
       sigma,myidx,myidy,npx,npy,commrow,commcol,comm2d)
    implicit none
    include 'mpif.h'
    double precision, intent(in) :: z
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot,myidx,myidy,npx,npy,commrow,&
         commcol,comm2d
    double precision, dimension(3), intent(out) :: center1,center2
    double precision, dimension(6), intent(out) :: sigma
    double precision:: den1,sqsum1,sqsum2,sqsum3,sqsum4,&
         epsx2,epsy2
    double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
         xpxfac,ypyfac
    double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
    double precision:: xpxlocal,ypylocal,zpzlocal
    double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
    double precision:: sqsum5local,sqsum6local
    double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
         pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
         z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
         z03,z04,pz03,pz04
    double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
         sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
    integer :: i,ierr,j
    double precision, dimension(7) :: localmax, glmax
    double precision, dimension(27) :: tmplc,tmpgl
    double precision, dimension(6) :: sigma1,sigma2
    double precision :: t0,lcrmax,glrmax
    integer :: npctmin,npctmax,npyhalf
    !        double precision :: alphax,alphay,alphaz
    !integer :: my_rank

    call starttime_Timer(t0)

    npctmin = 1
    npctmax = 1

    !        call MPI_COMM_RANK(mpicommwd,my_rank,ierr)
    !        call MPI_COMM_SIZE(mpicommwd,nproc,ierr)
    npyhalf = npy/2 

    den1 = 1.0/dble(nptot)
    x0lc = 0.0
    px0lc = 0.0
    y0lc = 0.0
    py0lc = 0.0
    z0lc = 0.0
    pz0lc = 0.0
    sqsum1local = 0.0
    sqsum2local = 0.0
    sqsum3local = 0.0
    sqsum4local = 0.0
    sqsum5local = 0.0
    sqsum6local = 0.0
    xpxlocal = 0.0
    ypylocal = 0.0
    zpzlocal = 0.0
    x0lc3 = 0.0
    x0lc4 = 0.0
    px0lc3 = 0.0
    px0lc4 = 0.0
    y0lc3 = 0.0
    y0lc4 = 0.0
    py0lc3 = 0.0
    py0lc4 = 0.0
    z0lc3 = 0.0
    z0lc4 = 0.0
    pz0lc3 = 0.0
    pz0lc4 = 0.0

    ! for cache optimization.
    if(innp.ne.0) then
       do i = 1, 6
          localmax(i) = abs(Pts1(i,1))
       enddo
       lcrmax = Pts1(1,1)**2+Pts1(3,1)**2
    else
       do i = 1, 6
          localmax(i) = 0.0
       enddo
       lcrmax = 0.0
    endif
    do i = 1, innp
       x0lc = x0lc + Pts1(1,i)
       sqsum1local = sqsum1local + Pts1(1,i)*Pts1(1,i)
       x0lc3 = x0lc3 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)
       x0lc4 = x0lc4 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)*&
            Pts1(1,i)
       xpxlocal = xpxlocal + Pts1(1,i)*Pts1(2,i)
       px0lc = px0lc + Pts1(2,i)
       sqsum2local = sqsum2local + Pts1(2,i)*Pts1(2,i)
       px0lc3 = px0lc3 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)
       px0lc4 = px0lc4 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)*&
            Pts1(2,i)
       y0lc = y0lc + Pts1(3,i)
       sqsum3local = sqsum3local + Pts1(3,i)*Pts1(3,i)
       y0lc3 = y0lc3 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)
       y0lc4 = y0lc4 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)*&
            Pts1(3,i)
       ypylocal = ypylocal + Pts1(3,i)*Pts1(4,i)
       py0lc = py0lc + Pts1(4,i)
       sqsum4local = sqsum4local + Pts1(4,i)*Pts1(4,i)
       py0lc3 = py0lc3 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)
       py0lc4 = py0lc4 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)*&
            Pts1(4,i)
       z0lc = z0lc + Pts1(5,i)
       sqsum5local = sqsum5local + Pts1(5,i)*Pts1(5,i)
       z0lc3 = z0lc3 + abs(Pts1(5,i)*Pts1(5,i)*Pts1(5,i))
       z0lc4 = z0lc4 + Pts1(5,i)*Pts1(5,i)*Pts1(5,i)*&
            Pts1(5,i)

       zpzlocal = zpzlocal + Pts1(5,i)*Pts1(6,i)
       pz0lc = pz0lc + Pts1(6,i)
       sqsum6local = sqsum6local + Pts1(6,i)*Pts1(6,i)
       pz0lc3 = pz0lc3 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)
       pz0lc4 = pz0lc4 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)*&
            Pts1(6,i)
       do j = 1, 6
          if(localmax(j).lt.abs(Pts1(j,i))) then
             localmax(j) = abs(Pts1(j,i))
          endif
       enddo
       if(lcrmax.lt.(Pts1(1,i)**2+Pts1(3,i)**2)) then
          lcrmax = Pts1(1,i)**2 + Pts1(3,i)**2
       endif
    enddo

    tmplc(1) = x0lc
    tmplc(2) = px0lc
    tmplc(3) = y0lc
    tmplc(4) = py0lc
    tmplc(5) = z0lc
    tmplc(6) = pz0lc
    tmplc(7) = sqsum1local
    tmplc(8) = sqsum2local
    tmplc(9) = sqsum3local
    tmplc(10) = sqsum4local
    tmplc(11) = sqsum5local
    tmplc(12) = sqsum6local
    tmplc(13) = xpxlocal
    tmplc(14) = ypylocal
    tmplc(15) = zpzlocal
    tmplc(16) = x0lc3
    tmplc(17) = x0lc4
    tmplc(18) = px0lc3
    tmplc(19) = px0lc4
    tmplc(20) = y0lc3
    tmplc(21) = y0lc4
    tmplc(22) = py0lc3
    tmplc(23) = py0lc4
    tmplc(24) = z0lc3
    tmplc(25) = z0lc4
    tmplc(26) = pz0lc3
    tmplc(27) = pz0lc4

    call MPI_ALLREDUCE(tmplc,tmpgl,27,MPI_DOUBLE_PRECISION,&
         MPI_SUM,commrow,ierr)
    tmplc = tmpgl
    call reduce4(tmplc,tmpgl,27,1,1,0,npyhalf,npyhalf,myidy,commcol)
    localmax(7) = lcrmax
    call MPI_ALLREDUCE(localmax,glmax,7,MPI_DOUBLE_PRECISION,&
         MPI_MAX,commrow,ierr)
    call reduce4(localmax,glmax,7,1,2,0,npyhalf,npyhalf,myidy,commcol)
    glrmax = glmax(7)

    if((myidx.eq.0).and.(myidy.eq.0)) then
       x0 = tmpgl(1)*den1
       px0 = tmpgl(2)*den1
       y0 = tmpgl(3)*den1
       py0 = tmpgl(4)*den1
       z0 = tmpgl(5)*den1
       pz0 = tmpgl(6)*den1
       sqx = tmpgl(7)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(8)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(9)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(10)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(11)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(12)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(13)*den1 - x0*px0
       ypy = tmpgl(14)*den1 - y0*py0
       zpz = tmpgl(15)*den1 - z0*pz0
       cubx = tmpgl(16)*den1
       fthx = tmpgl(17)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(18)*den1
       fthpx = tmpgl(19)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(20)*den1
       fthy = tmpgl(21)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(22)*den1
       fthpy = tmpgl(23)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(24)*den1
       fthz = tmpgl(25)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(26)*den1
       fthpz = tmpgl(27)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0
       ypyfac = 0.0
       zpzfac = 0.0
       if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
       if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
       if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
       write(24,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
       write(25,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
       write(26,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
       write(27,102)z,glmax(1),glmax(2),glmax(3),glmax(4),glmax(5),&
            glmax(6),sqrt(glrmax)
       write(28,101)z,npctmin,npctmax,nptot
       write(29,100)z,x03,px03,y03,py03,z03,pz03
       write(30,100)z,x04,px04,y04,py04,z04,pz04

       call flush(24)
       call flush(25)
       call flush(26)
       call flush(27)
       call flush(28)
       call flush(29)
       call flush(30)
       center1(1) = x0
       center1(2) = y0
       center1(3) = z0
       sigma1(1) = xrms
       sigma1(2) = pxrms
       sigma1(3) = yrms
       sigma1(4) = pyrms
       sigma1(5) = zrms
       sigma1(6) = pzrms
    endif

    call MPI_BCAST(center1,3,MPI_DOUBLE_PRECISION,0,comm2d,&
         ierr)
    call MPI_BCAST(sigma1,6,MPI_DOUBLE_PRECISION,0,comm2d,&
         ierr)

    if((myidx.eq.0).and.(myidy.eq.npyhalf)) then
       x0 = tmpgl(1)*den1
       px0 = tmpgl(2)*den1
       y0 = tmpgl(3)*den1
       py0 = tmpgl(4)*den1
       z0 = tmpgl(5)*den1
       pz0 = tmpgl(6)*den1
       sqx = tmpgl(7)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(8)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(9)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(10)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(11)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(12)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(13)*den1 - x0*px0
       ypy = tmpgl(14)*den1 - y0*py0
       zpz = tmpgl(15)*den1 - z0*pz0
       cubx = tmpgl(16)*den1
       fthx = tmpgl(17)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(18)*den1
       fthpx = tmpgl(19)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(20)*den1
       fthy = tmpgl(21)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(22)*den1
       fthpy = tmpgl(23)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(24)*den1
       fthz = tmpgl(25)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(26)*den1
       fthpz = tmpgl(27)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0
       ypyfac = 0.0
       zpzfac = 0.0
       if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
       if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
       if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
       write(34,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
       write(35,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
       write(36,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
       write(37,102)z,glmax(1),glmax(2),glmax(3),glmax(4),glmax(5),&
            glmax(6),sqrt(glrmax)
       write(38,101)z,npctmin,npctmax,nptot
       write(39,100)z,x03,px03,y03,py03,z03,pz03
       write(40,100)z,x04,px04,y04,py04,z04,pz04

       call flush(34)
       call flush(35)
       call flush(36)
       call flush(37)
       call flush(38)
       call flush(39)
       call flush(40)
       center2(1) = x0
       center2(2) = y0
       center2(3) = z0
       sigma2(1) = xrms
       sigma2(2) = pxrms
       sigma2(3) = yrms
       sigma2(4) = pyrms
       sigma2(5) = zrms
       sigma2(6) = pzrms
    endif

    call MPI_BCAST(center2,3,MPI_DOUBLE_PRECISION,npx*npyhalf,comm2d,&
         ierr)
    call MPI_BCAST(sigma2,6,MPI_DOUBLE_PRECISION,npx*npyhalf,comm2d,&
         ierr)

    if(myidy.lt.npyhalf) then
       sigma = sigma1
    else
       sigma = sigma2
    endif

100 format(7(1x,es13.6))
101 format(1x,es13.6,3I10)
102 format(8(1x,es13.6))

    t_diag = t_diag + elapsedtime_Timer(t0)

  end subroutine diagnostic2Gold_Output



  ! calculate <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
  ! <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance.
  ! using 1 Group processors for beam 1.
  subroutine diagnostic1G1_Output(z,Pts1,innp,nptot,center)
    implicit none
    include 'mpif.h'
    double precision, intent(in) :: z
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot
    double precision, dimension(3), intent(out) :: center
    double precision:: den1,sqsum1,sqsum2,sqsum3,sqsum4,&
         epsx2,epsy2
    double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
         xpxfac,ypyfac
    double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
    double precision:: xpxlocal,ypylocal,zpzlocal
    double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
    double precision:: sqsum5local,sqsum6local
    double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
         pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
         z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
         z03,z04,pz03,pz04
    double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
         sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
    integer :: i,my_rank,ierr,j
    double precision, dimension(6) :: localmax, glmax
    double precision, dimension(27) :: tmplc,tmpgl
    double precision :: t0,lcrmax,glrmax
    !        double precision :: alphax,alphay,alphaz
    integer :: npctmin,npctmax

    call starttime_Timer(t0)

    call MPI_COMM_RANK(mpicommwd,my_rank,ierr)

    den1 = 1.0/dble(nptot)
    x0lc = 0.0
    px0lc = 0.0
    y0lc = 0.0
    py0lc = 0.0
    z0lc = 0.0
    pz0lc = 0.0
    sqsum1local = 0.0
    sqsum2local = 0.0
    sqsum3local = 0.0
    sqsum4local = 0.0
    sqsum5local = 0.0
    sqsum6local = 0.0
    xpxlocal = 0.0
    ypylocal = 0.0
    zpzlocal = 0.0
    x0lc3 = 0.0
    x0lc4 = 0.0
    px0lc3 = 0.0
    px0lc4 = 0.0
    y0lc3 = 0.0
    y0lc4 = 0.0
    py0lc3 = 0.0
    py0lc4 = 0.0
    z0lc3 = 0.0
    z0lc4 = 0.0
    pz0lc3 = 0.0
    pz0lc4 = 0.0

    ! for cache optimization.
    if(innp.ne.0) then
       do i = 1, 6
          localmax(i) = abs(Pts1(i,1))
       enddo
       lcrmax = Pts1(1,1)**2+Pts1(3,1)**2
    else
       do i = 1, 6
          localmax(i) = 0.0
       enddo
       lcrmax = 0.0
    endif
    do i = 1, innp
       x0lc = x0lc + Pts1(1,i)
       sqsum1local = sqsum1local + Pts1(1,i)*Pts1(1,i)
       x0lc3 = x0lc3 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)
       x0lc4 = x0lc4 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)*&
            Pts1(1,i)
       xpxlocal = xpxlocal + Pts1(1,i)*Pts1(2,i)
       px0lc = px0lc + Pts1(2,i)
       sqsum2local = sqsum2local + Pts1(2,i)*Pts1(2,i)
       px0lc3 = px0lc3 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)
       px0lc4 = px0lc4 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)*&
            Pts1(2,i)
       y0lc = y0lc + Pts1(3,i)
       sqsum3local = sqsum3local + Pts1(3,i)*Pts1(3,i)
       y0lc3 = y0lc3 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)
       y0lc4 = y0lc4 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)*&
            Pts1(3,i)
       ypylocal = ypylocal + Pts1(3,i)*Pts1(4,i)
       py0lc = py0lc + Pts1(4,i)
       sqsum4local = sqsum4local + Pts1(4,i)*Pts1(4,i)
       py0lc3 = py0lc3 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)
       py0lc4 = py0lc4 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)*&
            Pts1(4,i)
       z0lc = z0lc + Pts1(5,i)
       sqsum5local = sqsum5local + Pts1(5,i)*Pts1(5,i)
       z0lc3 = z0lc3 + abs(Pts1(5,i)*Pts1(5,i)*Pts1(5,i))
       z0lc4 = z0lc4 + Pts1(5,i)*Pts1(5,i)*Pts1(5,i)*&
            Pts1(5,i)

       zpzlocal = zpzlocal + Pts1(5,i)*Pts1(6,i)
       pz0lc = pz0lc + Pts1(6,i)
       sqsum6local = sqsum6local + Pts1(6,i)*Pts1(6,i)
       pz0lc3 = pz0lc3 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)
       pz0lc4 = pz0lc4 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)*&
            Pts1(6,i)
       do j = 1, 6
          if(localmax(j).lt.abs(Pts1(j,i))) then
             localmax(j) = abs(Pts1(j,i))
          endif
       enddo
       if(lcrmax.lt.(Pts1(1,i)**2+Pts1(3,i)**2)) then
          lcrmax = Pts1(1,i)**2 + Pts1(3,i)**2
       endif
    enddo

    tmplc(1) = x0lc
    tmplc(2) = px0lc
    tmplc(3) = y0lc
    tmplc(4) = py0lc
    tmplc(5) = z0lc
    tmplc(6) = pz0lc
    tmplc(7) = sqsum1local
    tmplc(8) = sqsum2local
    tmplc(9) = sqsum3local
    tmplc(10) = sqsum4local
    tmplc(11) = sqsum5local
    tmplc(12) = sqsum6local
    tmplc(13) = xpxlocal
    tmplc(14) = ypylocal
    tmplc(15) = zpzlocal
    tmplc(16) = x0lc3
    tmplc(17) = x0lc4
    tmplc(18) = px0lc3
    tmplc(19) = px0lc4
    tmplc(20) = y0lc3
    tmplc(21) = y0lc4
    tmplc(22) = py0lc3
    tmplc(23) = py0lc4
    tmplc(24) = z0lc3
    tmplc(25) = z0lc4
    tmplc(26) = pz0lc3
    tmplc(27) = pz0lc4

    call MPI_REDUCE(tmplc,tmpgl,27,MPI_DOUBLE_PRECISION,&
         MPI_SUM,0,mpicommwd,ierr)
    call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
         mpicommwd,ierr)
    call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
         mpicommwd,ierr)
    call MPI_REDUCE(innp,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
         mpicommwd,ierr)
    call MPI_REDUCE(innp,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
         mpicommwd,ierr)

    if(my_rank.eq.0) then
       x0 = tmpgl(1)*den1
       px0 = tmpgl(2)*den1
       y0 = tmpgl(3)*den1
       py0 = tmpgl(4)*den1
       z0 = tmpgl(5)*den1
       pz0 = tmpgl(6)*den1
       sqx = tmpgl(7)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(8)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(9)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(10)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(11)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(12)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(13)*den1 - x0*px0
       ypy = tmpgl(14)*den1 - y0*py0
       zpz = tmpgl(15)*den1 - z0*pz0
       cubx = tmpgl(16)*den1
       fthx = tmpgl(17)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(18)*den1
       fthpx = tmpgl(19)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(20)*den1
       fthy = tmpgl(21)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(22)*den1
       fthpy = tmpgl(23)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(24)*den1
       fthz = tmpgl(25)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(26)*den1
       fthpz = tmpgl(27)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0
       ypyfac = 0.0
       zpzfac = 0.0
       if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
       if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
       if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
       write(24,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
       write(25,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
       write(26,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
       write(27,102)z,glmax(1),glmax(2),glmax(3),glmax(4),glmax(5),&
            glmax(6),sqrt(glrmax)
       write(28,101)z,npctmin,npctmax,nptot
       write(29,100)z,x03,px03,y03,py03,z03,pz03
       write(30,100)z,x04,px04,y04,py04,z04,pz04

       call flush(24)
       call flush(25)
       call flush(26)
       call flush(27)
       call flush(28)
       call flush(29)
       call flush(30)
       center(1) = x0
       center(2) = y0
       center(3) = z0
    endif

    call MPI_BCAST(center,3,MPI_DOUBLE_PRECISION,0,mpicommwd,&
         ierr)

100 format(7(1x,es13.6))
101 format(1x,es13.6,3I10)
102 format(8(1x,es13.6))

    t_diag = t_diag + elapsedtime_Timer(t0)

  end subroutine diagnostic1G1_Output



  ! calculate <x^2>,<xp>,<px^2>,x emittance, <y^2>,<ypy>,
  ! <py^2> and y emittance, <z^2>,<zp>,<pz^2>,z emittance.
  ! for beam 2
  subroutine diagnostic2G1_Output(z,Pts1,innp,nptot,center)
    implicit none
    include 'mpif.h'
    double precision, intent(in) :: z
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot
    double precision, dimension(3), intent(out) :: center
    double precision:: den1,sqsum1,sqsum2,sqsum3,sqsum4,&
         epsx2,epsy2
    double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
         xpxfac,ypyfac
    double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
    double precision:: xpxlocal,ypylocal,zpzlocal
    double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
    double precision:: sqsum5local,sqsum6local
    double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
         pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
         z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
         z03,z04,pz03,pz04
    double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
         sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
    integer :: i,my_rank,ierr,j
    double precision, dimension(6) :: localmax, glmax
    double precision, dimension(27) :: tmplc,tmpgl
    double precision :: t0,lcrmax,glrmax
    !        double precision :: alphax,alphay,alphaz
    integer :: npctmin,npctmax

    call starttime_Timer(t0)

    call MPI_COMM_RANK(mpicommwd,my_rank,ierr)

    den1 = 1.0/dble(nptot)
    x0lc = 0.0
    px0lc = 0.0
    y0lc = 0.0
    py0lc = 0.0
    z0lc = 0.0
    pz0lc = 0.0
    sqsum1local = 0.0
    sqsum2local = 0.0
    sqsum3local = 0.0
    sqsum4local = 0.0
    sqsum5local = 0.0
    sqsum6local = 0.0
    xpxlocal = 0.0
    ypylocal = 0.0
    zpzlocal = 0.0
    x0lc3 = 0.0
    x0lc4 = 0.0
    px0lc3 = 0.0
    px0lc4 = 0.0
    y0lc3 = 0.0
    y0lc4 = 0.0
    py0lc3 = 0.0
    py0lc4 = 0.0
    z0lc3 = 0.0
    z0lc4 = 0.0
    pz0lc3 = 0.0
    pz0lc4 = 0.0

    ! for cache optimization.
    if(innp.ne.0) then
       do i = 1, 6
          localmax(i) = abs(Pts1(i,1))
       enddo
       lcrmax = Pts1(1,1)**2+Pts1(3,1)**2
    else
       do i = 1, 6
          localmax(i) = 0.0
       enddo
       lcrmax = 0.0
    endif
    do i = 1, innp
       x0lc = x0lc + Pts1(1,i)
       sqsum1local = sqsum1local + Pts1(1,i)*Pts1(1,i)
       x0lc3 = x0lc3 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)
       x0lc4 = x0lc4 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)*&
            Pts1(1,i)
       xpxlocal = xpxlocal + Pts1(1,i)*Pts1(2,i)
       px0lc = px0lc + Pts1(2,i)
       sqsum2local = sqsum2local + Pts1(2,i)*Pts1(2,i)
       px0lc3 = px0lc3 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)
       px0lc4 = px0lc4 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)*&
            Pts1(2,i)
       y0lc = y0lc + Pts1(3,i)
       sqsum3local = sqsum3local + Pts1(3,i)*Pts1(3,i)
       y0lc3 = y0lc3 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)
       y0lc4 = y0lc4 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)*&
            Pts1(3,i)
       ypylocal = ypylocal + Pts1(3,i)*Pts1(4,i)
       py0lc = py0lc + Pts1(4,i)
       sqsum4local = sqsum4local + Pts1(4,i)*Pts1(4,i)
       py0lc3 = py0lc3 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)
       py0lc4 = py0lc4 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)*&
            Pts1(4,i)
       z0lc = z0lc + Pts1(5,i)
       sqsum5local = sqsum5local + Pts1(5,i)*Pts1(5,i)
       z0lc3 = z0lc3 + abs(Pts1(5,i)*Pts1(5,i)*Pts1(5,i))
       z0lc4 = z0lc4 + Pts1(5,i)*Pts1(5,i)*Pts1(5,i)*&
            Pts1(5,i)

       zpzlocal = zpzlocal + Pts1(5,i)*Pts1(6,i)
       pz0lc = pz0lc + Pts1(6,i)
       sqsum6local = sqsum6local + Pts1(6,i)*Pts1(6,i)
       pz0lc3 = pz0lc3 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)
       pz0lc4 = pz0lc4 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)*&
            Pts1(6,i)
       do j = 1, 6
          if(localmax(j).lt.abs(Pts1(j,i))) then
             localmax(j) = abs(Pts1(j,i))
          endif
       enddo
       if(lcrmax.lt.(Pts1(1,i)**2+Pts1(3,i)**2)) then
          lcrmax = Pts1(1,i)**2 + Pts1(3,i)**2
       endif
    enddo

    tmplc(1) = x0lc
    tmplc(2) = px0lc
    tmplc(3) = y0lc
    tmplc(4) = py0lc
    tmplc(5) = z0lc
    tmplc(6) = pz0lc
    tmplc(7) = sqsum1local
    tmplc(8) = sqsum2local
    tmplc(9) = sqsum3local
    tmplc(10) = sqsum4local
    tmplc(11) = sqsum5local
    tmplc(12) = sqsum6local
    tmplc(13) = xpxlocal
    tmplc(14) = ypylocal
    tmplc(15) = zpzlocal
    tmplc(16) = x0lc3
    tmplc(17) = x0lc4
    tmplc(18) = px0lc3
    tmplc(19) = px0lc4
    tmplc(20) = y0lc3
    tmplc(21) = y0lc4
    tmplc(22) = py0lc3
    tmplc(23) = py0lc4
    tmplc(24) = z0lc3
    tmplc(25) = z0lc4
    tmplc(26) = pz0lc3
    tmplc(27) = pz0lc4

    call MPI_REDUCE(tmplc,tmpgl,27,MPI_DOUBLE_PRECISION,&
         MPI_SUM,0,mpicommwd,ierr)
    call MPI_REDUCE(localmax,glmax,6,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
         mpicommwd,ierr)
    call MPI_REDUCE(lcrmax,glrmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
         mpicommwd,ierr)
    call MPI_REDUCE(innp,npctmin,1,MPI_INTEGER,MPI_MIN,0,&
         mpicommwd,ierr)
    call MPI_REDUCE(innp,npctmax,1,MPI_INTEGER,MPI_MAX,0,&
         mpicommwd,ierr)

    if(my_rank.eq.0) then
       x0 = tmpgl(1)*den1
       px0 = tmpgl(2)*den1
       y0 = tmpgl(3)*den1
       py0 = tmpgl(4)*den1
       z0 = tmpgl(5)*den1
       pz0 = tmpgl(6)*den1
       sqx = tmpgl(7)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(8)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(9)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(10)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(11)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(12)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(13)*den1 - x0*px0
       ypy = tmpgl(14)*den1 - y0*py0
       zpz = tmpgl(15)*den1 - z0*pz0
       cubx = tmpgl(16)*den1
       fthx = tmpgl(17)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(18)*den1
       fthpx = tmpgl(19)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(20)*den1
       fthy = tmpgl(21)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(22)*den1
       fthpy = tmpgl(23)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(24)*den1
       fthz = tmpgl(25)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(26)*den1
       fthpz = tmpgl(27)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0
       ypyfac = 0.0
       zpzfac = 0.0
       if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
       if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
       if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)
       write(34,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
       write(35,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
       write(36,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
       write(37,102)z,glmax(1),glmax(2),glmax(3),glmax(4),glmax(5),&
            glmax(6),sqrt(glrmax)
       write(38,101)z,npctmin,npctmax,nptot
       write(39,100)z,x03,px03,y03,py03,z03,pz03
       write(40,100)z,x04,px04,y04,py04,z04,pz04

       call flush(34)
       call flush(35)
       call flush(36)
       call flush(37)
       call flush(38)
       call flush(39)
       call flush(40)
       center(1) = x0
       center(2) = y0
       center(3) = z0
    endif

100 format(7(1x,es13.6))
101 format(1x,es13.6,3I10)
102 format(8(1x,es13.6))

    call MPI_BCAST(center,3,MPI_DOUBLE_PRECISION,0,mpicommwd,&
         ierr)

    t_diag = t_diag + elapsedtime_Timer(t0)

  end subroutine diagnostic2G1_Output



  !//Terminate MPI
  subroutine end_Output(time)
    implicit none
    include 'mpif.h'
    double precision, intent(inout) :: time
    double precision :: endtime, mtime
    integer :: my_rank,ierr

    call MPI_COMM_RANK(mpicommwd,my_rank,ierr)
    endtime = MPI_WTIME()
    time = endtime - time
    call MPI_REDUCE(time,mtime,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
         mpicommwd,ierr)

    if(my_rank.eq.0) then
       print*,"time: ",mtime
    endif

    !//for measurement of memory
    !        call system_stats()

    call MPI_Finalize(ierr)

  end subroutine end_Output



  !// calculate the 2d luminosity of the beam using 1 group processor.
  subroutine luminosity_Output(Pts1,Pts2,nptlc1,nptlc2,nxlum,nylum,it,&
       bcurr1,bcurr2,Np1,Np2)
    implicit none
    include 'mpif.h'
    integer :: nptlc1,nptlc2,nxlum,nylum,it,Np1,Np2
    double precision :: bcurr1,bcurr2
    double precision, pointer, dimension(:,:) :: Pts1,Pts2
    double precision, dimension(4) :: range1,range2
    double precision, dimension(2) :: localrange1,localrange2,temp1,temp2
    double precision :: xmin,xmax,ymin,ymax,epsx,epsy,hx,hy,hxi,hyi
    double precision :: ab,cd,lum
    double precision, dimension(nxlum,nylum) :: rho1,rho2,rholc
    integer :: i,j,k,ix,jx,ix1,jx1,ngrid,myid,ierr

    call MPI_COMM_RANK(mpicommwd,myid,ierr)

    range1(1) = Pts1(1,1)
    range1(2) = Pts1(1,1)
    range1(3) = Pts1(3,1)
    range1(4) = Pts1(3,1)
    do k = 1, nptlc1
       if(range1(1).gt.Pts1(1,k)) then
          range1(1) = Pts1(1,k)
       else if(range1(2).le.Pts1(1,k)) then
          range1(2) = Pts1(1,k)
       else
       endif
       if(range1(3).gt.Pts1(3,k)) then
          range1(3) = Pts1(3,k)
       else if(range1(4).le.Pts1(3,k)) then
          range1(4) = Pts1(3,k)
       else
       endif
    enddo

    range2(1) = Pts2(1,1)
    range2(2) = Pts2(1,1)
    range2(3) = Pts2(3,1)
    range2(4) = Pts2(3,1)
    do k = 1, nptlc2
       if(range2(1).gt.Pts2(1,k)) then
          range2(1) = Pts2(1,k)
       else if(range2(2).le.Pts2(1,k)) then
          range2(2) = Pts2(1,k)
       else
       endif
       if(range2(3).gt.Pts2(3,k)) then
          range2(3) = Pts2(3,k)
       else if(range2(4).le.Pts2(3,k)) then
          range2(4) = Pts2(3,k)
       else
       endif
    enddo

    do i = 1, 2
       if(range1(2*i-1).lt.range2(2*i-1)) then
          localrange1(i) = range1(2*i-1)
       else
          localrange1(i) = range2(2*i-1)
       endif
       if(range1(2*i).gt.range2(2*i)) then
          localrange2(i) = range1(2*i)
       else
          localrange2(i) = range2(2*i)
       endif
    enddo

    call MPI_ALLREDUCE(localrange1,temp1,2,MPI_DOUBLE_PRECISION,&
         MPI_MIN,mpicommwd,ierr)
    call MPI_ALLREDUCE(localrange2,temp2,2,MPI_DOUBLE_PRECISION,&
         MPI_MAX,mpicommwd,ierr)

    epsx = 1.0/(nxlum-3.0)
    epsy = 1.0/(nylum-3.0)

    xmin = temp1(1)-epsx*(temp2(1)-temp1(1))
    xmax = temp2(1)+epsx*(temp2(1)-temp1(1))
    ymin = temp1(2)-epsy*(temp2(2)-temp1(2))
    ymax = temp2(2)+epsy*(temp2(2)-temp1(2))

    hx = (xmax - xmin)/(nxlum - 1.0)
    hy = (ymax - ymin)/(nylum - 1.0)
    hxi = 1.0/hx
    hyi = 1.0/hy

    rholc = 0.0
    do i = 1, nptlc1
       ix=int((Pts1(1,i)-xmin)*hxi) + 1
       ab=((xmin-Pts1(1,i))+ix*hx)*hxi
       jx=int((Pts1(3,i)-ymin)*hyi) + 1
       cd=((ymin-Pts1(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.inny)) &
       !then
       !   innp
       !   stop
       !endif
       ! (i,j):
       rholc(ix,jx) = rholc(ix,jx) + ab*cd
       ! (i,j+1):
       rholc(ix,jx1) = rholc(ix,jx1) + ab*(1.0-cd)
       ! (i+1,j):
       rholc(ix1,jx) = rholc(ix1,jx)+(1.0-ab)*cd
       ! (i+1,j+1):
       rholc(ix1,jx1) = rholc(ix1,jx1)+(1.0-ab)*(1.0-cd)
    enddo
    ngrid = nxlum*nylum

    call MPI_REDUCE(rholc,rho1,ngrid,MPI_DOUBLE_PRECISION,&
         MPI_SUM,0,mpicommwd,ierr)

    rholc = 0.0
    do i = 1, nptlc2
       ix=int((Pts2(1,i)-xmin)*hxi) + 1
       ab=((xmin-Pts2(1,i))+ix*hx)*hxi
       jx=int((Pts2(3,i)-ymin)*hyi) + 1
       cd=((ymin-Pts2(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.inny)) &
       !then
       !   innp
       !   stop
       !endif
       ! (i,j):
       rholc(ix,jx) = rholc(ix,jx) + ab*cd
       ! (i,j+1):
       rholc(ix,jx1) = rholc(ix,jx1) + ab*(1.0-cd)
       ! (i+1,j):
       rholc(ix1,jx) = rholc(ix1,jx)+(1.0-ab)*cd
       ! (i+1,j+1):
       rholc(ix1,jx1) = rholc(ix1,jx1)+(1.0-ab)*(1.0-cd)
    enddo

    call MPI_REDUCE(rholc,rho2,ngrid,MPI_DOUBLE_PRECISION,&
         MPI_SUM,0,mpicommwd,ierr)

    if(myid.eq.0) then
       lum = 0.0
       do j = 1, nylum
          do i = 1, nxlum
             lum = lum + rho1(i,j)*rho2(i,j)
          enddo
       enddo
       open(1,file="luminosity.data",status="unknown",position="append")
       write(1,1000)it*1.0,&
            lum*bcurr1/Np1*bcurr2/Np2/(hx*hy)
       close(1)
    endif
1000 format(2(1x,es14.7))

  end subroutine luminosity_Output



  !//caculate the 2d luminosity using 2 group processors.
  subroutine luminosity2G_Output(Pts,nptlc,nxlum,nylum,it,&
       bcurr,Np,myidy,halfnpy)
    implicit none
    include 'mpif.h'
    integer :: nptlc,nxlum,nylum,it,Np,myidy,halfnpy
    double precision :: bcurr
    double precision, pointer, dimension(:,:) :: Pts
    double precision, dimension(2) :: range1,range2
    double precision, dimension(2) :: temp1,temp2
    double precision :: xmin,xmax,ymin,ymax,epsx,epsy,hx,hy,hxi,hyi
    double precision :: ab,cd,lum
    double precision, dimension(nxlum,nylum) :: rho1,rho2,rholc,rholc1,rholc2
    integer :: i,j,k,ix,jx,ix1,jx1,ngrid,myid,ierr

    call MPI_COMM_RANK(mpicommwd,myid,ierr)

    range1(1) = Pts(1,1) !minimum x
    range1(2) = Pts(3,1) !minimum y
    range2(1) = Pts(1,1) !maximum x
    range2(2) = Pts(3,1) !maximum y
    do k = 1, nptlc
       if(range1(1).gt.Pts(1,k)) then
          range1(1) = Pts(1,k)
       else if(range2(1).le.Pts(1,k)) then
          range2(1) = Pts(1,k)
       endif
       if(range1(2).gt.Pts(3,k)) then
          range1(2) = Pts(3,k)
       else if(range2(2).le.Pts(3,k)) then
          range2(2) = Pts(3,k)
       endif
    enddo

    call MPI_ALLREDUCE(range1,temp1,2,MPI_DOUBLE_PRECISION,&
         MPI_MIN,mpicommwd,ierr)
    call MPI_ALLREDUCE(range2,temp2,2,MPI_DOUBLE_PRECISION,&
         MPI_MAX,mpicommwd,ierr)

    epsx = 1.0/(nxlum-3.0)
    epsy = 1.0/(nylum-3.0)

    xmin = temp1(1)-epsx*(temp2(1)-temp1(1))
    xmax = temp2(1)+epsx*(temp2(1)-temp1(1))
    ymin = temp1(2)-epsy*(temp2(2)-temp1(2))
    ymax = temp2(2)+epsy*(temp2(2)-temp1(2))

    hx = (xmax - xmin)/(nxlum - 1.0)
    hy = (ymax - ymin)/(nylum - 1.0)
    hxi = 1.0/hx
    hyi = 1.0/hy

    !//get charge density on the grid using a CIC deposition scheme.
    rholc = 0.0
    do i = 1, nptlc
       ix=int((Pts(1,i)-xmin)*hxi) + 1
       ab=((xmin-Pts(1,i))+ix*hx)*hxi
       jx=int((Pts(3,i)-ymin)*hyi) + 1
       cd=((ymin-Pts(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.inny)) &
       !then
       !   innp
       !   stop
       !endif
       ! (i,j):
       rholc(ix,jx) = rholc(ix,jx) + ab*cd
       ! (i,j+1):
       rholc(ix,jx1) = rholc(ix,jx1) + ab*(1.0-cd)
       ! (i+1,j):
       rholc(ix1,jx) = rholc(ix1,jx)+(1.0-ab)*cd
       ! (i+1,j+1):
       rholc(ix1,jx1) = rholc(ix1,jx1)+(1.0-ab)*(1.0-cd)
    enddo

    !//separate out two beam density
    if(myidy.lt.halfnpy) then
       rholc1 = rholc*Bcurr/Np
       rholc2 = 0.0
    else
       rholc1 = 0.0
       rholc2 = rholc*Bcurr/Np
    endif

    ngrid = nxlum*nylum

    call MPI_REDUCE(rholc1,rho1,ngrid,MPI_DOUBLE_PRECISION,&
         MPI_SUM,0,mpicommwd,ierr)
    call MPI_REDUCE(rholc2,rho2,ngrid,MPI_DOUBLE_PRECISION,&
         MPI_SUM,0,mpicommwd,ierr)

    !//calculate the luminosity from density overlaping.
    if(myid.eq.0) then
       lum = 0.0
       do j = 1, nylum
          do i = 1, nxlum
             lum = lum + rho1(i,j)*rho2(i,j)
          enddo
       enddo
       open(1,file="luminosity2d.data",status="unknown",position="append")
       write(1,1000)it*1.0,lum/(hx*hy)
       !            lum*bcurr1/Np1*bcurr2/Np2/(hx*hy)
       close(1)
       !                          (xmax-xmin),(ymax-ymin)
    endif
1000 format(2(1x,es14.7))

  end subroutine luminosity2G_Output



  !//caculate the 2d luminosity using 2 group processors.
  subroutine luminosity2G3d_Output(Pts,nptlc,nxlum,nylum,&
       bcurr,wt,Np,myidy,halfnpy,lum)
    implicit none
    include 'mpif.h'
    integer :: nptlc,nxlum,nylum,Np,myidy,halfnpy
    double precision :: bcurr,wt
    double precision, dimension(:,:) :: Pts
    double precision, dimension(2) :: range1,range2
    double precision, dimension(2) :: temp1,temp2
    double precision :: xmin,xmax,ymin,ymax,epsx,epsy,hx,hy,hxi,hyi
    double precision :: ab,cd,lum
    double precision, dimension(nxlum,nylum) :: rho1,rho2,rholc,rholc1,rholc2
    integer :: i,j,k,ix,jx,ix1,jx1,ngrid,myid,ierr

    call MPI_COMM_RANK(mpicommwd,myid,ierr)

    range1(1) = Pts(1,1) !minimum x
    range1(2) = Pts(3,1) !minimum y
    range2(1) = Pts(1,1) !maximum x
    range2(2) = Pts(3,1) !maximum y
    do k = 1, nptlc
       if(range1(1).gt.Pts(1,k)) then
          range1(1) = Pts(1,k)
       else if(range2(1).le.Pts(1,k)) then
          range2(1) = Pts(1,k)
       else
       endif
       if(range1(2).gt.Pts(3,k)) then
          range1(2) = Pts(3,k)
       else if(range2(2).le.Pts(3,k)) then
          range2(2) = Pts(3,k)
       else
       endif
    enddo

    !//get the range of the global domain containing both beams.
    call MPI_ALLREDUCE(range1,temp1,2,MPI_DOUBLE_PRECISION,&
         MPI_MIN,mpicommwd,ierr)
    call MPI_ALLREDUCE(range2,temp2,2,MPI_DOUBLE_PRECISION,&
         MPI_MAX,mpicommwd,ierr)

    epsx = 1.0/(nxlum-3.0)
    epsy = 1.0/(nylum-3.0)

    xmin = temp1(1)-epsx*(temp2(1)-temp1(1))
    xmax = temp2(1)+epsx*(temp2(1)-temp1(1))
    ymin = temp1(2)-epsy*(temp2(2)-temp1(2))
    ymax = temp2(2)+epsy*(temp2(2)-temp1(2))

    hx = (xmax - xmin)/(nxlum - 1.0)
    hy = (ymax - ymin)/(nylum - 1.0)
    hxi = 1.0/hx
    hyi = 1.0/hy

    !//get charge density on the grid using a CIC deposition scheme.
    rholc = 0.0
    do i = 1, nptlc
       ix=int((Pts(1,i)-xmin)*hxi) + 1
       ab=((xmin-Pts(1,i))+ix*hx)*hxi
       jx=int((Pts(3,i)-ymin)*hyi) + 1
       cd=((ymin-Pts(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.inny)) &
       !then
       !   innp
       !   stop
       !endif
       ! (i,j):
       rholc(ix,jx) = rholc(ix,jx) + ab*cd
       ! (i,j+1):
       rholc(ix,jx1) = rholc(ix,jx1) + ab*(1.0-cd)
       ! (i+1,j):
       rholc(ix1,jx) = rholc(ix1,jx)+(1.0-ab)*cd
       ! (i+1,j+1):
       rholc(ix1,jx1) = rholc(ix1,jx1)+(1.0-ab)*(1.0-cd)
    enddo

    !//separate out two beam density
    if(myidy.lt.halfnpy) then
       rholc1 = rholc*bcurr*wt/Np
       rholc2 = 0.0
    else
       rholc1 = 0.0
       rholc2 = rholc*bcurr*wt/Np
    endif

    ngrid = nxlum*nylum

    call MPI_ALLREDUCE(rholc1,rho1,ngrid,MPI_DOUBLE_PRECISION,&
         MPI_SUM,mpicommwd,ierr)
    call MPI_ALLREDUCE(rholc2,rho2,ngrid,MPI_DOUBLE_PRECISION,&
         MPI_SUM,mpicommwd,ierr)

    !//calculate the luminosity from density overlaping.
    lum = 0.0
    do j = 1, nylum
       do i = 1, nxlum
          lum = lum + rho1(i,j)*rho2(i,j)/(hx*hy)
       enddo
    enddo

  end subroutine luminosity2G3d_Output



  !//caculate the 3d luminosity using 1 group processors.
  subroutine luminosity1G3d_Output(Pts1,nptlc1,Pts2,nptlc2,nxlum,nylum,&
       bcurr1,bcurr2,wt1,wt2,Np1,Np2,lum)
    implicit none
    include 'mpif.h'
    integer :: nptlc1,nptlc2,nxlum,nylum,Np1,Np2
    double precision :: bcurr1,bcurr2,wt1,wt2
    double precision, dimension(:,:) :: Pts1,Pts2
    double precision, dimension(2) :: range1,range2
    double precision, dimension(2) :: temp1,temp2
    double precision :: xmin,xmax,ymin,ymax,epsx,epsy,hx,hy,hxi,hyi
    double precision :: ab,cd,lum
    double precision, dimension(nxlum,nylum) :: rho1,rho2,rholc1,rholc2
    integer :: i,j,k,ix,jx,ix1,jx1,ngrid,myid,ierr

    call MPI_COMM_RANK(mpicommwd,myid,ierr)

    range1(1) = Pts1(1,1) !minimum x
    range1(2) = Pts1(3,1) !minimum y
    range2(1) = Pts1(1,1) !maximum x
    range2(2) = Pts1(3,1) !maximum y
    do k = 1, nptlc1
       if(range1(1).gt.Pts1(1,k)) then
          range1(1) = Pts1(1,k)
       else if(range2(1).le.Pts1(1,k)) then
          range2(1) = Pts1(1,k)
       else
       endif
       if(range1(2).gt.Pts1(3,k)) then
          range1(2) = Pts1(3,k)
       else if(range2(2).le.Pts1(3,k)) then
          range2(2) = Pts1(3,k)
       else
       endif
    enddo
    do k = 1, nptlc2
       if(range1(1).gt.Pts2(1,k)) then
          range1(1) = Pts2(1,k)
       else if(range2(1).le.Pts2(1,k)) then
          range2(1) = Pts2(1,k)
       else
       endif
       if(range1(2).gt.Pts2(3,k)) then
          range1(2) = Pts2(3,k)
       else if(range2(2).le.Pts2(3,k)) then
          range2(2) = Pts2(3,k)
       else
       endif
    enddo

    !//get the range of the global domain containing both beams.
    call MPI_ALLREDUCE(range1,temp1,2,MPI_DOUBLE_PRECISION,&
         MPI_MIN,mpicommwd,ierr)
    call MPI_ALLREDUCE(range2,temp2,2,MPI_DOUBLE_PRECISION,&
         MPI_MAX,mpicommwd,ierr)

    epsx = 1.0/(nxlum-3.0)
    epsy = 1.0/(nylum-3.0)

    xmin = temp1(1)-epsx*(temp2(1)-temp1(1))
    xmax = temp2(1)+epsx*(temp2(1)-temp1(1))
    ymin = temp1(2)-epsy*(temp2(2)-temp1(2))
    ymax = temp2(2)+epsy*(temp2(2)-temp1(2))

    hx = (xmax - xmin)/(nxlum - 1.0)
    hy = (ymax - ymin)/(nylum - 1.0)
    hxi = 1.0/hx
    hyi = 1.0/hy

    !//get charge density on the grid using a CIC deposition scheme.
    rholc1 = 0.0
    do i = 1, nptlc1
       ix=int((Pts1(1,i)-xmin)*hxi) + 1
       ab=((xmin-Pts1(1,i))+ix*hx)*hxi
       jx=int((Pts1(3,i)-ymin)*hyi) + 1
       cd=((ymin-Pts1(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       ! (i,j):
       rholc1(ix,jx) = rholc1(ix,jx) + ab*cd
       ! (i,j+1):
       rholc1(ix,jx1) = rholc1(ix,jx1) + ab*(1.0-cd)
       ! (i+1,j):
       rholc1(ix1,jx) = rholc1(ix1,jx)+(1.0-ab)*cd
       ! (i+1,j+1):
       rholc1(ix1,jx1) = rholc1(ix1,jx1)+(1.0-ab)*(1.0-cd)
    enddo
    rholc1 = rholc1*bcurr1*wt1/Np1
    rholc2 = 0.0
    do i = 1, nptlc2
       ix=int((Pts2(1,i)-xmin)*hxi) + 1
       ab=((xmin-Pts2(1,i))+ix*hx)*hxi
       jx=int((Pts2(3,i)-ymin)*hyi) + 1
       cd=((ymin-Pts2(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       ! (i,j):
       rholc2(ix,jx) = rholc2(ix,jx) + ab*cd
       ! (i,j+1):
       rholc2(ix,jx1) = rholc2(ix,jx1) + ab*(1.0-cd)
       ! (i+1,j):
       rholc2(ix1,jx) = rholc2(ix1,jx)+(1.0-ab)*cd
       ! (i+1,j+1):
       rholc2(ix1,jx1) = rholc2(ix1,jx1)+(1.0-ab)*(1.0-cd)
    enddo
    rholc2 = rholc2*bcurr2*wt2/Np2

    ngrid = nxlum*nylum
    call MPI_ALLREDUCE(rholc1,rho1,ngrid,MPI_DOUBLE_PRECISION,&
         MPI_SUM,mpicommwd,ierr)
    call MPI_ALLREDUCE(rholc2,rho2,ngrid,MPI_DOUBLE_PRECISION,&
         MPI_SUM,mpicommwd,ierr)

    !//calculate the luminosity from density overlaping.
    lum = 0.0
    do j = 1, nylum
       do i = 1, nxlum
          lum = lum + rho1(i,j)*rho2(i,j)/(hx*hy)
       enddo
    enddo

  end subroutine luminosity1G3d_Output



  subroutine phase_Output(nfile,Pts1,Nptlocal,nsamp)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nfile,Nptlocal,nsamp
    double precision, pointer, dimension(:,:) :: Pts1
    integer :: np,my_rank,ierr
    integer status(MPI_STATUS_SIZE)
    integer :: i,j,sixnpt,mnpt,nstride
    integer, allocatable, dimension(:) :: nptlist
    double precision, allocatable,dimension(:,:) :: recvbuf

    call MPI_COMM_RANK(mpicommwd,my_rank,ierr)
    call MPI_COMM_SIZE(mpicommwd,np,ierr)
    call MPI_ALLREDUCE(Nptlocal,mnpt,1,MPI_INTEGER,MPI_MAX,&
         mpicommwd,ierr)

    !nstride = 1
    nstride = nsamp
    allocate(nptlist(0:np-1))
    nptlist = 0
    allocate(recvbuf(6,mnpt))
    sixnpt = 6*Nptlocal

    call MPI_GATHER(Nptlocal,1,MPI_INTEGER,nptlist,1,&
         MPI_INTEGER,0,mpicommwd,ierr)

    nptlist = 6*nptlist

    if(my_rank.eq.0) then
       open(nfile,status='unknown')
       do i = 1, Nptlocal, nstride
          write(nfile,100)Pts1(1,i),Pts1(2,i),Pts1(3,i),&
               Pts1(4,i),Pts1(5,i),Pts1(6,i)
       enddo
       do i = 1, np-1
          call MPI_RECV(recvbuf(1,1),nptlist(i),MPI_DOUBLE_PRECISION,&
               i,1,mpicommwd,status,ierr) 

          do j = 1, nptlist(i)/6, nstride
             write(nfile,100)recvbuf(1,j),recvbuf(2,j),recvbuf(3,j),&
                  recvbuf(4,j),recvbuf(5,j),recvbuf(6,j)
          enddo
       enddo
       close(nfile)
    else
       call MPI_SEND(Pts1(1,1),sixnpt,MPI_DOUBLE_PRECISION,0,1,&
            mpicommwd,ierr)
    endif

100 format(6(1x,es14.7))

    deallocate(nptlist)
    deallocate(recvbuf)

  end subroutine phase_Output



  !readin the particle distribution for restart purpose
  subroutine phaseinBB_Output(Bpts,myid,iturnend,idbunch,Nplocal,&
       close2g,iticend,iseedend)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: myid,idbunch
    integer, intent(inout) :: Nplocal
    integer, intent(out) :: iturnend,iticend
    integer*8, intent(out) :: iseedend
    double precision, dimension(:,:) :: Bpts
    double precision, dimension(2) :: close2g
    integer :: i,j,k,l,m,n
    character*8 name1
    character*9 name2
    character*10 name3
    character*11 name4

    name1 = 'ph0000sx'
    name2 = 'ph0000sxx'
    name3 = 'ph0000sxxx'
    name4 = 'ph0000sxxxx'

    if(idbunch < 10) then
       if(myid.lt.10) then
          name1(6:6) = char(idbunch+48)
          name1(8:8) = char(myid+48)
          open(9,file=name1,status="unknown",form="unformatted")
       else if((myid.ge.10).and.(myid.lt.100)) then
          name2(6:6) = char(idbunch+48)
          i = myid/10
          j = myid - 10*i
          name2(8:8) = char(i+48)
          name2(9:9) = char(j+48)
          open(9,file=name2,status="unknown",form="unformatted")
       else if((myid.ge.100).and.(myid.lt.1000)) then
          name3(6:6) = char(idbunch+48)
          i = myid/100
          j = myid - 100*i
          k = j/10
          l = j - 10*k
          name3(8:8) = char(i+48)
          name3(9:9) = char(k+48)
          name3(10:10) = char(l+48)
          open(9,file=name3,status="unknown",form="unformatted")
       else
          name4(6:6) = char(idbunch+48)
          i = myid/1000
          j = myid - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - 10*m
          name4(8:8) = char(i+48)
          name4(9:9) = char(k+48)
          name4(10:10) = char(m+48)
          name4(11:11) = char(n+48)
          open(9,file=name4,status="unknown",form="unformatted")
       endif
    else if((idbunch.ge.10).and.(idbunch.lt.100)) then
       if(myid.lt.10) then
          i = idbunch/10
          j = idbunch - 10*i
          name1(5:5) = char(i+48)
          name1(6:6) = char(j+48)
          name1(8:8) = char(myid+48)
          open(9,file=name1,status="unknown",form="unformatted")
       else if((myid.ge.10).and.(myid.lt.100)) then
          i = idbunch/10
          j = idbunch - 10*i
          name2(5:5) = char(i+48)
          name2(6:6) = char(j+48)
          i = myid/10
          j = myid - 10*i
          name2(8:8) = char(i+48)
          name2(9:9) = char(j+48)
          open(9,file=name2,status="unknown",form="unformatted")
       else if((myid.ge.100).and.(myid.lt.1000)) then
          i = idbunch/10
          j = idbunch - 10*i
          name3(5:5) = char(i+48)
          name3(6:6) = char(j+48)
          i = myid/100
          j = myid - 100*i
          k = j/10
          l = j - 10*k
          name3(8:8) = char(i+48)
          name3(9:9) = char(k+48)
          name3(10:10) = char(l+48)
          open(9,file=name3,status="unknown",form="unformatted")
       else
          i = idbunch/10
          j = idbunch - 10*i
          name4(5:5) = char(i+48)
          name4(6:6) = char(j+48)
          i = myid/1000
          j = myid - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - 10*m
          name4(8:8) = char(i+48)
          name4(9:9) = char(k+48)
          name4(10:10) = char(m+48)
          name4(11:11) = char(n+48)
          open(9,file=name4,status="unknown",form="unformatted")
       endif
    else if((idbunch.ge.100).and.(idbunch.lt.1000)) then
       if(myid.lt.10) then
          i = idbunch/100
          j = idbunch - 100*i
          k = j/10
          l = j - 10*k
          name1(4:4) = char(i+48)
          name1(5:5) = char(k+48)
          name1(6:6) = char(l+48)
          name1(8:8) = char(myid+48)
          open(9,file=name1,status="unknown",form="unformatted")
       else if((myid.ge.10).and.(myid.lt.100)) then
          i = idbunch/100
          j = idbunch - 100*i
          k = j/10
          l = j - 10*k
          name2(4:4) = char(i+48)
          name2(5:5) = char(k+48)
          name2(6:6) = char(l+48)
          i = myid/10
          j = myid - 10*i
          name2(8:8) = char(i+48)
          name2(9:9) = char(j+48)
          open(9,file=name2,status="unknown",form="unformatted")
       else if((myid.ge.100).and.(myid.lt.1000)) then
          i = idbunch/100
          j = idbunch - 100*i
          k = j/10
          l = j - 10*k
          name3(4:4) = char(i+48)
          name3(5:5) = char(k+48)
          name3(6:6) = char(l+48)
          i = myid/100
          j = myid - 100*i
          k = j/10
          l = j - 10*k
          name3(8:8) = char(i+48)
          name3(9:9) = char(k+48)
          name3(10:10) = char(l+48)
          open(9,file=name3,status="unknown",form="unformatted")
       else
          i = idbunch/100
          j = idbunch - 100*i
          k = j/10
          l = j - 10*k
          name4(4:4) = char(i+48)
          name4(5:5) = char(k+48)
          name4(6:6) = char(l+48)
          i = myid/1000
          j = myid - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - 10*m
          name4(8:8) = char(i+48)
          name4(9:9) = char(k+48)
          name4(10:10) = char(m+48)
          name4(11:11) = char(n+48)
          open(9,file=name4,status="unknown",form="unformatted")
       endif
    else
       if(myid.lt.10) then
          i = idbunch/1000
          j = idbunch - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - m*10
          name1(3:3) = char(i+48)
          name1(4:4) = char(k+48)
          name1(5:5) = char(m+48)
          name1(6:6) = char(n+48)
          name1(8:8) = char(myid+48)
          open(9,file=name1,status="unknown",form="unformatted")
       else if((myid.ge.10).and.(myid.lt.100)) then
          i = idbunch/1000
          j = idbunch - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - m*10
          name2(3:3) = char(i+48)
          name2(4:4) = char(k+48)
          name2(5:5) = char(m+48)
          name2(6:6) = char(n+48)
          i = myid/10
          j = myid - 10*i
          name2(8:8) = char(i+48)
          name2(9:9) = char(j+48)
          open(9,file=name2,status="unknown",form="unformatted")
       else if((myid.ge.100).and.(myid.lt.1000)) then
          i = idbunch/1000
          j = idbunch - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - m*10
          name3(3:3) = char(i+48)
          name3(4:4) = char(k+48)
          name3(5:5) = char(m+48)
          name3(6:6) = char(n+48)
          i = myid/100
          j = myid - 100*i
          k = j/10
          l = j - 10*k
          name3(8:8) = char(i+48)
          name3(9:9) = char(k+48)
          name3(10:10) = char(l+48)
          open(9,file=name3,status="unknown",form="unformatted")
       else
          i = idbunch/1000
          j = idbunch - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - m*10
          name4(3:3) = char(i+48)
          name4(4:4) = char(k+48)
          name4(5:5) = char(m+48)
          name4(6:6) = char(n+48)
          i = myid/1000
          j = myid - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - 10*m
          name4(8:8) = char(i+48)
          name4(9:9) = char(k+48)
          name4(10:10) = char(m+48)
          name4(11:11) = char(n+48)
          open(9,file=name4,status="unknown",form="unformatted")
       endif
    endif

    read(9)Nplocal,iturnend,iticend,iseedend,close2g(1),close2g(2)
    !        allocate(Bpts(6,Nplocal))
    read(9)Bpts(1:6,1:Nplocal)

    close(9)

  end subroutine phaseinBB_Output



  !output the particle information for restart function
  subroutine phaseoutBB_Output(Bpts,myid,iturn,idbunch,Nplocal,&
       close2g,itic,iseedinit)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: myid,iturn,idbunch,Nplocal,itic
    integer*8, intent(in) :: iseedinit
    double precision, dimension(:,:) :: Bpts
    double precision, dimension(2) :: close2g
    integer :: i,j,k,l,m,n
    character*8 name1
    character*9 name2
    character*10 name3
    character*11 name4

    name1 = 'ph0000sx'
    name2 = 'ph0000sxx'
    name3 = 'ph0000sxxx'
    name4 = 'ph0000sxxxx'

    if(idbunch < 10) then
       if(myid.lt.10) then
          name1(6:6) = char(idbunch+48)
          name1(8:8) = char(myid+48)
          open(9,file=name1,status="unknown",form="unformatted")
       else if((myid.ge.10).and.(myid.lt.100)) then
          name2(6:6) = char(idbunch+48)
          i = myid/10
          j = myid - 10*i
          name2(8:8) = char(i+48)
          name2(9:9) = char(j+48)
          open(9,file=name2,status="unknown",form="unformatted")
       else if((myid.ge.100).and.(myid.lt.1000)) then
          name3(6:6) = char(idbunch+48)
          i = myid/100
          j = myid - 100*i
          k = j/10
          l = j - 10*k
          name3(8:8) = char(i+48)
          name3(9:9) = char(k+48)
          name3(10:10) = char(l+48)
          open(9,file=name3,status="unknown",form="unformatted")
       else
          name4(6:6) = char(idbunch+48)
          i = myid/1000
          j = myid - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - 10*m
          name4(8:8) = char(i+48)
          name4(9:9) = char(k+48)
          name4(10:10) = char(m+48)
          name4(11:11) = char(n+48)
          open(9,file=name4,status="unknown",form="unformatted")
       endif
    else if((idbunch.ge.10).and.(idbunch.lt.100)) then
       if(myid.lt.10) then
          i = idbunch/10
          j = idbunch - 10*i
          name1(5:5) = char(i+48)
          name1(6:6) = char(j+48)
          name1(8:8) = char(myid+48)
          open(9,file=name1,status="unknown",form="unformatted")
       else if((myid.ge.10).and.(myid.lt.100)) then
          i = idbunch/10
          j = idbunch - 10*i
          name2(5:5) = char(i+48)
          name2(6:6) = char(j+48)
          i = myid/10
          j = myid - 10*i
          name2(8:8) = char(i+48)
          name2(9:9) = char(j+48)
          open(9,file=name2,status="unknown",form="unformatted")
       else if((myid.ge.100).and.(myid.lt.1000)) then
          i = idbunch/10
          j = idbunch - 10*i
          name3(5:5) = char(i+48)
          name3(6:6) = char(j+48)
          i = myid/100
          j = myid - 100*i
          k = j/10
          l = j - 10*k
          name3(8:8) = char(i+48)
          name3(9:9) = char(k+48)
          name3(10:10) = char(l+48)
          open(9,file=name3,status="unknown",form="unformatted")
       else
          i = idbunch/10
          j = idbunch - 10*i
          name4(5:5) = char(i+48)
          name4(6:6) = char(j+48)
          i = myid/1000
          j = myid - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - 10*m
          name4(8:8) = char(i+48)
          name4(9:9) = char(k+48)
          name4(10:10) = char(m+48)
          name4(11:11) = char(n+48)
          open(9,file=name4,status="unknown",form="unformatted")
       endif
    else if((idbunch.ge.100).and.(idbunch.lt.1000)) then
       if(myid.lt.10) then
          i = idbunch/100
          j = idbunch - 100*i
          k = j/10
          l = j - 10*k
          name1(4:4) = char(i+48)
          name1(5:5) = char(k+48)
          name1(6:6) = char(l+48)
          name1(8:8) = char(myid+48)
          open(9,file=name1,status="unknown",form="unformatted")
       else if((myid.ge.10).and.(myid.lt.100)) then
          i = idbunch/100
          j = idbunch - 100*i
          k = j/10
          l = j - 10*k
          name2(4:4) = char(i+48)
          name2(5:5) = char(k+48)
          name2(6:6) = char(l+48)
          i = myid/10
          j = myid - 10*i
          name2(8:8) = char(i+48)
          name2(9:9) = char(j+48)
          open(9,file=name2,status="unknown",form="unformatted")
       else if((myid.ge.100).and.(myid.lt.1000)) then
          i = idbunch/100
          j = idbunch - 100*i
          k = j/10
          l = j - 10*k
          name3(4:4) = char(i+48)
          name3(5:5) = char(k+48)
          name3(6:6) = char(l+48)
          i = myid/100
          j = myid - 100*i
          k = j/10
          l = j - 10*k
          name3(8:8) = char(i+48)
          name3(9:9) = char(k+48)
          name3(10:10) = char(l+48)
          open(9,file=name3,status="unknown",form="unformatted")
       else
          i = idbunch/100
          j = idbunch - 100*i
          k = j/10
          l = j - 10*k
          name4(4:4) = char(i+48)
          name4(5:5) = char(k+48)
          name4(6:6) = char(l+48)
          i = myid/1000
          j = myid - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - 10*m
          name4(8:8) = char(i+48)
          name4(9:9) = char(k+48)
          name4(10:10) = char(m+48)
          name4(11:11) = char(n+48)
          open(9,file=name4,status="unknown",form="unformatted")
       endif
    else
       if(myid.lt.10) then
          i = idbunch/1000
          j = idbunch - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - m*10
          name1(3:3) = char(i+48)
          name1(4:4) = char(k+48)
          name1(5:5) = char(m+48)
          name1(6:6) = char(n+48)
          name1(8:8) = char(myid+48)
          open(9,file=name1,status="unknown",form="unformatted")
       else if((myid.ge.10).and.(myid.lt.100)) then
          i = idbunch/1000
          j = idbunch - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - m*10
          name2(3:3) = char(i+48)
          name2(4:4) = char(k+48)
          name2(5:5) = char(m+48)
          name2(6:6) = char(n+48)
          i = myid/10
          j = myid - 10*i
          name2(8:8) = char(i+48)
          name2(9:9) = char(j+48)
          open(9,file=name2,status="unknown",form="unformatted")
       else if((myid.ge.100).and.(myid.lt.1000)) then
          i = idbunch/1000
          j = idbunch - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - m*10
          name3(3:3) = char(i+48)
          name3(4:4) = char(k+48)
          name3(5:5) = char(m+48)
          name3(6:6) = char(n+48)
          i = myid/100
          j = myid - 100*i
          k = j/10
          l = j - 10*k
          name3(8:8) = char(i+48)
          name3(9:9) = char(k+48)
          name3(10:10) = char(l+48)
          open(9,file=name3,status="unknown",form="unformatted")
       else
          i = idbunch/1000
          j = idbunch - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - m*10
          name4(3:3) = char(i+48)
          name4(4:4) = char(k+48)
          name4(5:5) = char(m+48)
          name4(6:6) = char(n+48)
          i = myid/1000
          j = myid - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - 10*m
          name4(8:8) = char(i+48)
          name4(9:9) = char(k+48)
          name4(10:10) = char(m+48)
          name4(11:11) = char(n+48)
          open(9,file=name4,status="unknown",form="unformatted")
       endif
    endif

    write(9)Nplocal,iturn,itic,iseedinit,close2g(1),close2g(2)
    write(9)Bpts(1:6,1:Nplocal)

    close(9)

  end subroutine phaseoutBB_Output



  !output the particle information for restart function
  subroutine ptclout_Output(Bpts,myid,iturn,idbunch,Nplocal,npOut,nturns,outRate,sampleInd)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: myid,iturn,idbunch,Nplocal,npOut,nturns,outRate
    double precision, dimension(:,:) :: Bpts
    double precision, allocatable, dimension(:,:) :: tmppts
    integer :: i,j,k,l,m,n
    character*8 name1
    character*9 name2
    character*10 name3
    character*11 name4
    integer, allocatable, dimension(:), intent(in) :: sampleInd

    name1 = 'pt0000sx'
    name2 = 'pt0000sxx'
    name3 = 'pt0000sxxx'
    name4 = 'pt0000sxxxx'

    if(idbunch < 10) then
       if(myid.lt.10) then
          name1(6:6) = char(idbunch+48)
          name1(8:8) = char(myid+48)
          open(9,file=name1,status="unknown",position="append",form="unformatted")
       else if((myid.ge.10).and.(myid.lt.100)) then
          name2(6:6) = char(idbunch+48)
          i = myid/10
          j = myid - 10*i
          name2(8:8) = char(i+48)
          name2(9:9) = char(j+48)
          open(9,file=name2,status="unknown",position="append",form="unformatted")
       else if((myid.ge.100).and.(myid.lt.1000)) then
          name3(6:6) = char(idbunch+48)
          i = myid/100
          j = myid - 100*i
          k = j/10
          l = j - 10*k
          name3(8:8) = char(i+48)
          name3(9:9) = char(k+48)
          name3(10:10) = char(l+48)
          open(9,file=name3,status="unknown",position="append",form="unformatted")
       else
          name4(6:6) = char(idbunch+48)
          i = myid/1000
          j = myid - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - 10*m
          name4(8:8) = char(i+48)
          name4(9:9) = char(k+48)
          name4(10:10) = char(m+48)
          name4(11:11) = char(n+48)
          open(9,file=name4,status="unknown",position="append",form="unformatted")
       endif
    else if((idbunch.ge.10).and.(idbunch.lt.100)) then
       if(myid.lt.10) then
          i = idbunch/10
          j = idbunch - 10*i
          name1(5:5) = char(i+48)
          name1(6:6) = char(j+48)
          name1(8:8) = char(myid+48)
          open(9,file=name1,status="unknown",position="append",form="unformatted")
       else if((myid.ge.10).and.(myid.lt.100)) then
          i = idbunch/10
          j = idbunch - 10*i
          name2(5:5) = char(i+48)
          name2(6:6) = char(j+48)
          i = myid/10
          j = myid - 10*i
          name2(8:8) = char(i+48)
          name2(9:9) = char(j+48)
          open(9,file=name2,status="unknown",position="append",form="unformatted")
       else if((myid.ge.100).and.(myid.lt.1000)) then
          i = idbunch/10
          j = idbunch - 10*i
          name3(5:5) = char(i+48)
          name3(6:6) = char(j+48)
          i = myid/100
          j = myid - 100*i
          k = j/10
          l = j - 10*k
          name3(8:8) = char(i+48)
          name3(9:9) = char(k+48)
          name3(10:10) = char(l+48)
          open(9,file=name3,status="unknown",position="append",form="unformatted")
       else
          i = idbunch/10
          j = idbunch - 10*i
          name4(5:5) = char(i+48)
          name4(6:6) = char(j+48)
          i = myid/1000
          j = myid - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - 10*m
          name4(8:8) = char(i+48)
          name4(9:9) = char(k+48)
          name4(10:10) = char(m+48)
          name4(11:11) = char(n+48)
          open(9,file=name4,status="unknown",position="append",form="unformatted")
       endif
    else if((idbunch.ge.100).and.(idbunch.lt.1000)) then
       if(myid.lt.10) then
          i = idbunch/100
          j = idbunch - 100*i
          k = j/10
          l = j - 10*k
          name1(4:4) = char(i+48)
          name1(5:5) = char(k+48)
          name1(6:6) = char(l+48)
          name1(8:8) = char(myid+48)
          open(9,file=name1,status="unknown",position="append",form="unformatted")
       else if((myid.ge.10).and.(myid.lt.100)) then
          i = idbunch/100
          j = idbunch - 100*i
          k = j/10
          l = j - 10*k
          name2(4:4) = char(i+48)
          name2(5:5) = char(k+48)
          name2(6:6) = char(l+48)
          i = myid/10
          j = myid - 10*i
          name2(8:8) = char(i+48)
          name2(9:9) = char(j+48)
          open(9,file=name2,status="unknown",position="append",form="unformatted")
       else if((myid.ge.100).and.(myid.lt.1000)) then
          i = idbunch/100
          j = idbunch - 100*i
          k = j/10
          l = j - 10*k
          name3(4:4) = char(i+48)
          name3(5:5) = char(k+48)
          name3(6:6) = char(l+48)
          i = myid/100
          j = myid - 100*i
          k = j/10
          l = j - 10*k
          name3(8:8) = char(i+48)
          name3(9:9) = char(k+48)
          name3(10:10) = char(l+48)
          open(9,file=name3,status="unknown",position="append",form="unformatted")
       else
          i = idbunch/100
          j = idbunch - 100*i
          k = j/10
          l = j - 10*k
          name4(4:4) = char(i+48)
          name4(5:5) = char(k+48)
          name4(6:6) = char(l+48)
          i = myid/1000
          j = myid - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - 10*m
          name4(8:8) = char(i+48)
          name4(9:9) = char(k+48)
          name4(10:10) = char(m+48)
          name4(11:11) = char(n+48)
          open(9,file=name4,status="unknown",position="append",form="unformatted")
       endif
    else
       if(myid.lt.10) then
          i = idbunch/1000
          j = idbunch - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - m*10
          name1(3:3) = char(i+48)
          name1(4:4) = char(k+48)
          name1(5:5) = char(m+48)
          name1(6:6) = char(n+48)
          name1(8:8) = char(myid+48)
          open(9,file=name1,status="unknown",position="append",form="unformatted")
       else if((myid.ge.10).and.(myid.lt.100)) then
          i = idbunch/1000
          j = idbunch - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - m*10
          name2(3:3) = char(i+48)
          name2(4:4) = char(k+48)
          name2(5:5) = char(m+48)
          name2(6:6) = char(n+48)
          i = myid/10
          j = myid - 10*i
          name2(8:8) = char(i+48)
          name2(9:9) = char(j+48)
          open(9,file=name2,status="unknown",position="append",form="unformatted")
       else if((myid.ge.100).and.(myid.lt.1000)) then
          i = idbunch/1000
          j = idbunch - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - m*10
          name3(3:3) = char(i+48)
          name3(4:4) = char(k+48)
          name3(5:5) = char(m+48)
          name3(6:6) = char(n+48)
          i = myid/100
          j = myid - 100*i
          k = j/10
          l = j - 10*k
          name3(8:8) = char(i+48)
          name3(9:9) = char(k+48)
          name3(10:10) = char(l+48)
          open(9,file=name3,status="unknown",position="append",form="unformatted")
       else
          i = idbunch/1000
          j = idbunch - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - m*10
          name4(3:3) = char(i+48)
          name4(4:4) = char(k+48)
          name4(5:5) = char(m+48)
          name4(6:6) = char(n+48)
          i = myid/1000
          j = myid - 1000*i
          k = j/100
          l = j - 100*k
          m = l/10
          n = l - 10*m
          name4(8:8) = char(i+48)
          name4(9:9) = char(k+48)
          name4(10:10) = char(m+48)
          name4(11:11) = char(n+48)
          open(9,file=name4,status="unknown",position="append",form="unformatted")
       endif
    endif

    allocate(tmppts(6,npOut))
    do i = 1, npOut
       tmppts(:,i) = Bpts(:,sampleInd(i))
    enddo

    !print*,"sampleInd: ",sum(sampleInd),npOut

    if(iturn==0) then
       write(9)nturns,outRate,npOut
    endif
    write(9)iturn
    write(9)tmppts

    close(9)
    deallocate(tmppts)

  end subroutine ptclout_Output



  subroutine diagnostic2GMBg_Output(z,Pts1,innp,nptot,&
       myidx,myidy,npx,npy,commrow,commcol,comm2d,bid,&
       gambet,mygroupid)
    implicit none
    include 'mpif.h'
    double precision, intent(in) :: z
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot,myidx,myidy,npx,npy,commrow,&
         commcol,comm2d,bid,mygroupid
    double precision:: den1,sqsum1,sqsum2,sqsum3,sqsum4,&
         epsx2,epsy2,gambet
    double precision:: xpx,ypy,epx,epy,xrms,pxrms,yrms,pyrms,&
         xpxfac,ypyfac
    double precision:: sqsum1local,sqsum2local,sqsum3local,sqsum4local
    double precision:: xpxlocal,ypylocal,zpzlocal
    double precision:: sqsum5,sqsum6,epsz2,zpz,epz,zrms,pzrms,zpzfac
    double precision:: sqsum5local,sqsum6local
    double precision:: x0lc,x0,px0lc,px0,y0lc,y0,py0lc,py0,z0lc,z0,&
         pz0lc,pz0,x0lc3,x0lc4,px0lc3,px0lc4,y0lc3,y0lc4,py0lc3,py0lc4,&
         z0lc3,z0lc4,pz0lc3,pz0lc4,x03,x04,px03,px04,y03,y04,py03,py04,&
         z03,z04,pz03,pz04
    double precision ::sqx,cubx,fthx,sqpx,cubpx,fthpx,sqy,cuby,fthy,&
         sqpy,cubpy,fthpy,sqz,cubz,fthz,sqpz,cubpz,fthpz
    integer :: i,my_rank,ierr,j
    double precision, dimension(6) :: localmax
    double precision, dimension(14) :: localmax2, glmax
    double precision, dimension(54) :: tmplc,tmpgl
    double precision :: t0,lcrmax
    !        double precision :: alphax,alphay,alphaz
    integer :: npctmin,npctmax,npyhalf,nfile
    integer, parameter :: nbin = 400
    integer, dimension(2*nbin) :: tmpbin,glbin
    double precision :: f90,f95,f99,ex90,ex95,ex99,ex902,ex952,&
         ex992,ey90,ey95,ey99,ey902,ey952,ey992
    double precision,allocatable,dimension(:) :: epsiontmp
    double precision, dimension(2) :: tmpepslc,tmpep,Ealpha,Ebeta,&
         Egamma
    double precision :: xtmp,pxtmp,ytmp,pytmp,tmp1,tmp2,epsmylc,hyeps,&
         epsmxlc,hxeps,epsmx,epsmy
    integer :: nii,iitmp

    call starttime_Timer(t0)

    npctmin = 1
    npctmax = 1

    call MPI_COMM_RANK(mpicommwd,my_rank,ierr)
    !        call MPI_COMM_SIZE(mpicommwd,nproc,ierr)
    npyhalf = npy/2 

    den1 = 1.0/dble(nptot)
    x0lc = 0.0
    px0lc = 0.0
    y0lc = 0.0
    py0lc = 0.0
    z0lc = 0.0
    pz0lc = 0.0
    sqsum1local = 0.0
    sqsum2local = 0.0
    sqsum3local = 0.0
    sqsum4local = 0.0
    sqsum5local = 0.0
    sqsum6local = 0.0
    xpxlocal = 0.0
    ypylocal = 0.0
    zpzlocal = 0.0
    x0lc3 = 0.0
    x0lc4 = 0.0
    px0lc3 = 0.0
    px0lc4 = 0.0
    y0lc3 = 0.0
    y0lc4 = 0.0
    py0lc3 = 0.0
    py0lc4 = 0.0
    z0lc3 = 0.0
    z0lc4 = 0.0
    pz0lc3 = 0.0
    pz0lc4 = 0.0

    ! for cache optimization.
    if(innp.ne.0) then
       do i = 1, 6
          localmax(i) = abs(Pts1(i,1))
       enddo
       lcrmax = Pts1(1,1)**2+Pts1(3,1)**2
    else
       do i = 1, 6
          localmax(i) = 0.0
       enddo
       lcrmax = 0.0
    endif
    do i = 1, innp
       x0lc = x0lc + Pts1(1,i)
       sqsum1local = sqsum1local + Pts1(1,i)*Pts1(1,i)
       x0lc3 = x0lc3 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)
       x0lc4 = x0lc4 + Pts1(1,i)*Pts1(1,i)*Pts1(1,i)*&
            Pts1(1,i)
       xpxlocal = xpxlocal + Pts1(1,i)*Pts1(2,i)
       px0lc = px0lc + Pts1(2,i)
       sqsum2local = sqsum2local + Pts1(2,i)*Pts1(2,i)
       px0lc3 = px0lc3 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)
       px0lc4 = px0lc4 + Pts1(2,i)*Pts1(2,i)*Pts1(2,i)*&
            Pts1(2,i)
       y0lc = y0lc + Pts1(3,i)
       sqsum3local = sqsum3local + Pts1(3,i)*Pts1(3,i)
       y0lc3 = y0lc3 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)
       y0lc4 = y0lc4 + Pts1(3,i)*Pts1(3,i)*Pts1(3,i)*&
            Pts1(3,i)
       ypylocal = ypylocal + Pts1(3,i)*Pts1(4,i)
       py0lc = py0lc + Pts1(4,i)
       sqsum4local = sqsum4local + Pts1(4,i)*Pts1(4,i)
       py0lc3 = py0lc3 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)
       py0lc4 = py0lc4 + Pts1(4,i)*Pts1(4,i)*Pts1(4,i)*&
            Pts1(4,i)
       z0lc = z0lc + Pts1(5,i)
       sqsum5local = sqsum5local + Pts1(5,i)*Pts1(5,i)
       z0lc3 = z0lc3 + abs(Pts1(5,i)*Pts1(5,i)*Pts1(5,i))
       z0lc4 = z0lc4 + Pts1(5,i)*Pts1(5,i)*Pts1(5,i)*&
            Pts1(5,i)

       zpzlocal = zpzlocal + Pts1(5,i)*Pts1(6,i)
       pz0lc = pz0lc + Pts1(6,i)
       sqsum6local = sqsum6local + Pts1(6,i)*Pts1(6,i)
       pz0lc3 = pz0lc3 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)
       pz0lc4 = pz0lc4 + Pts1(6,i)*Pts1(6,i)*Pts1(6,i)*&
            Pts1(6,i)
       do j = 1, 6
          if(localmax(j).lt.abs(Pts1(j,i))) then
             localmax(j) = abs(Pts1(j,i))
          endif
       enddo
       if(lcrmax.lt.(Pts1(1,i)**2+Pts1(3,i)**2)) then
          lcrmax = Pts1(1,i)**2 + Pts1(3,i)**2
       endif
    enddo

    localmax2 = 0.0
    tmplc = 0.0
    if(myidy.lt.npyhalf) then
       tmplc(1) = x0lc
       tmplc(2) = px0lc
       tmplc(3) = y0lc
       tmplc(4) = py0lc
       tmplc(5) = z0lc
       tmplc(6) = pz0lc
       tmplc(7) = sqsum1local
       tmplc(8) = sqsum2local
       tmplc(9) = sqsum3local
       tmplc(10) = sqsum4local
       tmplc(11) = sqsum5local
       tmplc(12) = sqsum6local
       tmplc(13) = xpxlocal
       tmplc(14) = ypylocal
       tmplc(15) = zpzlocal
       tmplc(16) = x0lc3
       tmplc(17) = x0lc4
       tmplc(18) = px0lc3
       tmplc(19) = px0lc4
       tmplc(20) = y0lc3
       tmplc(21) = y0lc4
       tmplc(22) = py0lc3
       tmplc(23) = py0lc4
       tmplc(24) = z0lc3
       tmplc(25) = z0lc4
       tmplc(26) = pz0lc3
       tmplc(27) = pz0lc4
       localmax2(1) = localmax(1)
       localmax2(2) = localmax(2)
       localmax2(3) = localmax(3)
       localmax2(4) = localmax(4)
       localmax2(5) = localmax(5)
       localmax2(6) = localmax(6)
       localmax2(7) = lcrmax
    else
       tmplc(28) = x0lc
       tmplc(29) = px0lc
       tmplc(30) = y0lc
       tmplc(31) = py0lc
       tmplc(32) = z0lc
       tmplc(33) = pz0lc
       tmplc(34) = sqsum1local
       tmplc(35) = sqsum2local
       tmplc(36) = sqsum3local
       tmplc(37) = sqsum4local
       tmplc(38) = sqsum5local
       tmplc(39) = sqsum6local
       tmplc(40) = xpxlocal
       tmplc(41) = ypylocal
       tmplc(42) = zpzlocal
       tmplc(43) = x0lc3
       tmplc(44) = x0lc4
       tmplc(45) = px0lc3
       tmplc(46) = px0lc4
       tmplc(47) = y0lc3
       tmplc(48) = y0lc4
       tmplc(49) = py0lc3
       tmplc(50) = py0lc4
       tmplc(51) = z0lc3
       tmplc(52) = z0lc4
       tmplc(53) = pz0lc3
       tmplc(54) = pz0lc4
       localmax2(8) = localmax(1)
       localmax2(9) = localmax(2)
       localmax2(10) = localmax(3)
       localmax2(11) = localmax(4)
       localmax2(12) = localmax(5)
       localmax2(13) = localmax(6)
       localmax2(14) = lcrmax
    endif

    !call MPI_REDUCE(tmplc,tmpgl,54,MPI_DOUBLE_PRECISION,&
    !                MPI_SUM,0,mpicommwd,ierr)
    !call MPI_ALLREDUCE(tmplc,tmpgl,54,MPI_DOUBLE_PRECISION,&
    !                MPI_SUM,mpicommwd,ierr)
    call MPI_ALLREDUCE(tmplc,tmpgl,54,MPI_DOUBLE_PRECISION,&
         MPI_SUM,comm2d,ierr)
    !call MPI_REDUCE(localmax2,glmax,14,MPI_DOUBLE_PRECISION,MPI_MAX,0,&
    !                mpicommwd,ierr)
    !call MPI_ALLREDUCE(localmax2,glmax,14,MPI_DOUBLE_PRECISION,MPI_MAX,&
    !                mpicommwd,ierr)
    call MPI_ALLREDUCE(localmax2,glmax,14,MPI_DOUBLE_PRECISION,MPI_MAX,&
         comm2d,ierr)

    if(myidy.lt.npyhalf) then
       !//output beam 1
       x0 = tmpgl(1)*den1
       px0 = tmpgl(2)*den1
       y0 = tmpgl(3)*den1
       py0 = tmpgl(4)*den1
       z0 = tmpgl(5)*den1
       pz0 = tmpgl(6)*den1
       sqx = tmpgl(7)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(8)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(9)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(10)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(11)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(12)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(13)*den1 - x0*px0
       ypy = tmpgl(14)*den1 - y0*py0
       zpz = tmpgl(15)*den1 - z0*pz0
       cubx = tmpgl(16)*den1
       fthx = tmpgl(17)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(18)*den1
       fthpx = tmpgl(19)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(20)*den1
       fthy = tmpgl(21)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(22)*den1
       fthpy = tmpgl(23)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(24)*den1
       fthz = tmpgl(25)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(26)*den1
       fthpz = tmpgl(27)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0
       ypyfac = 0.0
       zpzfac = 0.0
       if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
       if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
       if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)

       if(myidx.eq.0 .and. myidy.eq.0) then
          nfile = 504+bid*50 + mygroupid*200
          write(nfile,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
          call flush(nfile)
          nfile = 505+bid*50 + mygroupid*200
          write(nfile,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
          call flush(nfile)
          !          nfile = 6+bid*40 + mygroupid*200
          !          write(nfile,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
          !          call flush(nfile)
          !          nfile = 7+bid*40 + mygroupid*200
          !          write(nfile,102)z,glmax(1),glmax(2),glmax(3),glmax(4),glmax(5),&
          !                       glmax(6),sqrt(glmax(7))
          !          call flush(nfile)
          !          nfile = 8+bid*40
          !          write(nfile,101)z,npctmin,npctmax,nptot
          !          call flush(nfile)
          !          nfile = 9+bid*40 + mygroupid*200
          !          write(nfile,100)z,x03,px03,y03,py03,z03,pz03
          !          call flush(nfile)
          !          nfile = 10+bid*40 + mygroupid*200
          !          write(nfile,100)z,x04,px04,y04,py04,z04,pz04
          !          call flush(nfile)
       endif
    else
       !//output beam 2
       x0 = tmpgl(28)*den1
       px0 = tmpgl(29)*den1
       y0 = tmpgl(30)*den1
       py0 = tmpgl(31)*den1
       z0 = tmpgl(32)*den1
       pz0 = tmpgl(33)*den1
       sqx = tmpgl(34)*den1
       sqsum1 = sqx - x0*x0
       sqpx = tmpgl(35)*den1
       sqsum2 = sqpx - px0*px0
       sqy = tmpgl(36)*den1
       sqsum3 = sqy - y0*y0
       sqpy = tmpgl(37)*den1
       sqsum4 = sqpy - py0*py0
       sqz = tmpgl(38)*den1
       sqsum5 = sqz - z0*z0
       sqpz = tmpgl(39)*den1
       sqsum6 = sqpz - pz0*pz0
       xpx = tmpgl(40)*den1 - x0*px0
       ypy = tmpgl(41)*den1 - y0*py0
       zpz = tmpgl(42)*den1 - z0*pz0
       cubx = tmpgl(43)*den1
       fthx = tmpgl(44)*den1
       x03 = (abs(cubx-3*sqx*x0+2*x0*x0*x0))**(1.0/3.0)
       x04 = sqrt(sqrt(abs(fthx-4*cubx*x0+6*sqx*x0*x0-3*x0*x0*x0*x0)))
       cubpx = tmpgl(45)*den1
       fthpx = tmpgl(46)*den1
       px03 = (abs(cubpx-3*sqpx*px0+2*px0*px0*px0))**(1.0/3.0)
       px04 = sqrt(sqrt(abs(fthpx-4*cubpx*px0+6*sqpx*px0*px0-&
            3*px0*px0*px0*px0)))
       cuby = tmpgl(47)*den1
       fthy = tmpgl(48)*den1
       y03 = (abs(cuby-3*sqy*y0+2*y0*y0*y0))**(1.0/3.0)
       y04 = sqrt(sqrt(abs(fthy-4*cuby*y0+6*sqy*y0*y0-3*y0*y0*y0*y0)))
       cubpy = tmpgl(49)*den1
       fthpy = tmpgl(50)*den1
       py03 = (abs(cubpy-3*sqpy*py0+2*py0*py0*py0))**(1.0/3.0)
       py04 = sqrt(sqrt(abs(fthpy-4*cubpy*py0+6*sqpy*py0*py0-&
            3*py0*py0*py0*py0)))
       cubz = tmpgl(51)*den1
       fthz = tmpgl(52)*den1
       z03 = (abs(cubz-3*sqz*z0+2*z0*z0*z0))**(1.0/3.0)
       z04 = sqrt(sqrt(abs(fthz-4*cubz*z0+6*sqz*z0*z0-3*z0*z0*z0*z0)))
       cubpz = tmpgl(53)*den1
       fthpz = tmpgl(54)*den1
       pz03 = (abs(cubpz-3*sqpz*pz0+2*pz0*pz0*pz0))**(1.0/3.0)
       pz04 = sqrt(sqrt(abs(fthpz-4*cubpz*pz0+6*sqpz*pz0*pz0-&
            3*pz0*pz0*pz0*pz0)))
       epsx2 = (sqsum1*sqsum2-xpx*xpx)
       epsy2 = (sqsum3*sqsum4-ypy*ypy)
       epsz2 = (sqsum5*sqsum6-zpz*zpz)
       epx = sqrt(max(epsx2,0.0d0))
       epy = sqrt(max(epsy2,0.0d0))
       epz = sqrt(max(epsz2,0.0d0))
       xrms = sqrt(sqsum1)
       pxrms = sqrt(sqsum2)
       yrms = sqrt(sqsum3)
       pyrms = sqrt(sqsum4)
       zrms = sqrt(sqsum5)
       pzrms = sqrt(sqsum6)
       xpxfac = 0.0
       ypyfac = 0.0
       zpzfac = 0.0
       if(xrms.ne.0.0 .and. pxrms.ne.0.0)xpxfac=1.0/(xrms*pxrms)
       if(yrms.ne.0.0 .and. pyrms.ne.0.0)ypyfac=1.0/(yrms*pyrms)
       if(zrms.ne.0.0 .and. pzrms.ne.0.0)zpzfac=1.0/(zrms*pzrms)

       if(myidx.eq.0 .and. myidy.eq.npyhalf) then
          nfile = 514+bid*50 + mygroupid*200
          write(nfile,100)z,x0,xrms,px0,pxrms,-xpx/epx,epx
          call flush(nfile)
          nfile = 515+bid*50 + mygroupid*200
          write(nfile,100)z,y0,yrms,py0,pyrms,-ypy/epy,epy
          call flush(nfile)
          !          nfile = 16+bid*40 + mygroupid*200
          !          write(nfile,100)z,z0,zrms,pz0,pzrms,-zpz/epz,epz
          !          call flush(nfile)
          !          nfile = 17+bid*40 + mygroupid*200
          !          write(nfile,102)z,glmax(8),glmax(9),glmax(10),glmax(11),glmax(12),&
          !                       glmax(13),sqrt(glmax(14))
          !          call flush(nfile)
          !          nfile = 18+bid*40
          !          write(nfile,101)z,npctmin,npctmax,nptot
          !          call flush(nfile)
          !          nfile = 19+bid*40 + mygroupid*200
          !          write(nfile,100)z,x03,px03,y03,py03,z03,pz03
          !          call flush(nfile)
          !          nfile = 20+bid*40 + mygroupid*200
          !          write(nfile,100)z,x04,px04,y04,py04,z04,pz04
          !          call flush(nfile)
       endif
    endif


    goto 999

    Ealpha(1) = -xpx/epx
    Ealpha(2) = -ypy/epy
    Ebeta(1) = xrms*xrms*gambet/epx
    Ebeta(2) = yrms*yrms*gambet/epy
    Egamma(:) = (1.0+Ealpha(:)*Ealpha(:))/Ebeta(:)

    allocate(epsiontmp(innp))
    epsmxlc = -1.0e10
    do i = 1, innp
       xtmp = Pts1(1,i) - x0
       tmp1 = Pts1(2,i)
       pxtmp = (tmp1 - px0)/gambet
       epsiontmp(i)=Egamma(1)*xtmp*xtmp+2*Ealpha(1)*xtmp*pxtmp+&
            Ebeta(1)*pxtmp*pxtmp
       if(epsmxlc.le.epsiontmp(i)) epsmxlc = epsiontmp(i)
    enddo

    tmpepslc = 0.0d0
    if(myidy.lt.npyhalf) then
       tmpepslc(1) = epsmxlc
    else
       tmpepslc(2) = epsmxlc
    endif

    call MPI_ALLREDUCE(tmpepslc,tmpep,2,MPI_DOUBLE_PRECISION,&
         MPI_MAX,mpicommwd,ierr)
    if(myidy.lt.npyhalf) then
       epsmx = tmpep(1) 
    else
       epsmx = tmpep(2) 
    endif

    hxeps = epsmx*1.0001/nbin
    tmpbin = 0


    if(myidy.lt.npyhalf) then
       do i = 1, innp
          iitmp = int(epsiontmp(i)/hxeps)
          nii = iitmp+1
          tmpbin(nii) = tmpbin(nii) + 1
       enddo
    else
       do i = 1, innp
          iitmp = int(epsiontmp(i)/hxeps)
          nii = iitmp+1
          tmpbin(nii+nbin) = tmpbin(nii+nbin) + 1
       enddo
    endif

    if(myidy.eq.0) then
       do i = 1, 2*nbin
       enddo
    endif

    call MPI_REDUCE(tmpbin,glbin,2*nbin,MPI_INTEGER,&
         MPI_SUM,0,mpicommwd,ierr)

    !        if(my_rank.eq.0) then
    !          do i = 1, 2*nbin
    !          enddo
    !        endif

    f90 = 0.999*nptot
    f95 = 0.9999*nptot
    f99 = 0.99999*nptot
    if(my_rank.eq.0) then

       hxeps = tmpep(1)*1.0001/nbin
       do i = 2, nbin
          glbin(i) = glbin(i) + glbin(i-1)
       enddo
       do i = 1, nbin
          if(glbin(i).gt.f90) then
             ex90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                  hxeps*(i-1)
             exit
          endif
       enddo
       do i = 1, nbin
          if(glbin(i).gt.f95) then
             ex95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                  hxeps*(i-1)
             exit
          endif
       enddo
       do i =1, nbin
          if(glbin(i).gt.f99) then
             ex99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                  hxeps*(i-1)
             exit
          endif
       enddo
       !          ex90 = ex90*gambet
       !          ex95 = ex95*gambet
       !          ex99 = ex99*gambet

       hxeps = tmpep(2)*1.0001/nbin
       do i = nbin+2, 2*nbin
          glbin(i) = glbin(i) + glbin(i-1)
       enddo
       do i = nbin+1, 2*nbin
          if(glbin(i).gt.f90) then
             ex902 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                  hxeps*(i-nbin-1)
             exit
          endif
       enddo
       do i = nbin+1, 2*nbin
          if(glbin(i).gt.f95) then
             ex952 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                  hxeps*(i-nbin-1)
             exit
          endif
       enddo
       do i =nbin+1, 2*nbin
          if(glbin(i).gt.f99) then
             ex992 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hxeps+&
                  hxeps*(i-nbin-1)
             exit
          endif
       enddo
       !          ex902 = ex902*gambet
       !          ex952 = ex952*gambet
       !          ex992 = ex992*gambet

       nfile = 11+bid*20
       write(nfile,100)z,ex90,ex95,ex99,ex902,ex952,ex992
       call flush(nfile)
    endif

    epsmylc = -1.0e10
    do i = 1, innp
       ytmp = Pts1(3,i) - y0
       tmp2 = Pts1(4,i)
       pytmp = (tmp2 - py0)/gambet
       epsiontmp(i)=Egamma(2)*ytmp*ytmp+2*Ealpha(2)*ytmp*pytmp+&
            Ebeta(2)*pytmp*pytmp
       if(epsmylc.le.epsiontmp(i)) epsmylc = epsiontmp(i)
    enddo

    tmpepslc = 0.0d0
    if(myidy.lt.npyhalf) then
       tmpepslc(1) = epsmylc
    else
       tmpepslc(2) = epsmylc
    endif
    call MPI_ALLREDUCE(tmpepslc,tmpep,2,MPI_DOUBLE_PRECISION,&
         MPI_MAX,mpicommwd,ierr)
    if(myidy.lt.npyhalf) then
       epsmy = tmpep(1)
    else
       epsmy = tmpep(2)
    endif

    hyeps = epsmy*1.0001/nbin
    tmpbin = 0
    if(myidy.lt.npyhalf) then
       do i = 1, innp
          iitmp = int(epsiontmp(i)/hyeps)
          nii = iitmp+1
          tmpbin(nii) = tmpbin(nii) + 1
       enddo
    else
       do i = 1, innp
          iitmp = int(epsiontmp(i)/hyeps)
          nii = iitmp+1
          tmpbin(nii+nbin) = tmpbin(nii+nbin) + 1
       enddo
    endif
    call MPI_REDUCE(tmpbin,glbin,2*nbin,MPI_INTEGER,&
         MPI_SUM,0,mpicommwd,ierr)

    if(my_rank.eq.0) then
       hyeps = tmpep(1)*1.0001/nbin
       do i = 2, nbin
          glbin(i) = glbin(i) + glbin(i-1)
       enddo
       do i = 1, nbin
          if(glbin(i).gt.f90) then
             ey90 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                  hyeps*(i-1)
             exit
          endif
       enddo
       do i = 1, nbin
          if(glbin(i).gt.f95) then
             ey95 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                  hyeps*(i-1)
             exit
          endif
       enddo
       do i = 1, nbin
          if(glbin(i).gt.f99) then
             ey99 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                  hyeps*(i-1)
             exit
          endif
       enddo
       !          ey90 = ey90*gambet
       !          ey95 = ey95*gambet
       !          ey99 = ey99*gambet

       hyeps = tmpep(2)*1.0001/nbin
       do i = nbin+2, 2*nbin
          glbin(i) = glbin(i) + glbin(i-1)
       enddo
       do i = nbin+1, 2*nbin
          if(glbin(i).gt.f90) then
             ey902 = ((f90 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                  hyeps*(i-nbin-1)
             exit
          endif
       enddo
       do i = nbin+1, 2*nbin
          if(glbin(i).gt.f95) then
             ey952 = ((f95 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                  hyeps*(i-nbin-1)
             exit
          endif
       enddo
       do i =nbin+1, 2*nbin
          if(glbin(i).gt.f99) then
             ey992 = ((f99 - glbin(i-1))/(glbin(i)-glbin(i-1)))*hyeps+&
                  hyeps*(i-nbin-1)
             exit
          endif
       enddo
       !          ey902 = ey902*gambet
       !          ey952 = ey952*gambet
       !          ey992 = ey992*gambet

       nfile = 21+bid*20
       write(nfile,100)z,ey90,ey95,ey99,ey902,ey952,ey992
       call flush(nfile)
    endif
    deallocate(epsiontmp)

999 continue

100 format(7(1x,es13.6))
!101 format(1x,es13.6,3I10)
!102 format(8(1x,es13.6))

    t_diag = t_diag + elapsedtime_Timer(t0)

  end subroutine diagnostic2GMBg_Output



  !//caculate the 2d luminosity using 2 group processors.
  subroutine luminosity2Gg_Output(Pts,nptlc,nxlum,nylum,it,&
       bcurr,Np,myidy,halfnpy,mygroupid)
    implicit none
    include 'mpif.h'
    integer :: nptlc,nxlum,nylum,it,Np,myidy,halfnpy,mygroupid
    double precision :: bcurr
    double precision, pointer, dimension(:,:) :: Pts
    double precision, dimension(2) :: range1,range2
    double precision, dimension(2) :: temp1,temp2
    double precision :: xmin,xmax,ymin,ymax,epsx,epsy,hx,hy,hxi,hyi
    double precision :: ab,cd,lum
    double precision, dimension(nxlum,nylum) :: rho1,rho2,rholc,rholc1,rholc2
    integer :: i,j,k,ix,jx,ix1,jx1,ngrid,myid,ierr,nfilelum

    call MPI_COMM_RANK(mpicommwd,myid,ierr)

    range1(1) = Pts(1,1) !minimum x
    range1(2) = Pts(3,1) !minimum y
    range2(1) = Pts(1,1) !maximum x
    range2(2) = Pts(3,1) !maximum y
    do k = 1, nptlc
       if(range1(1).gt.Pts(1,k)) then
          range1(1) = Pts(1,k)
       else if(range2(1).le.Pts(1,k)) then
          range2(1) = Pts(1,k)
       else
       endif
       if(range1(2).gt.Pts(3,k)) then
          range1(2) = Pts(3,k)
       else if(range2(2).le.Pts(3,k)) then
          range2(2) = Pts(3,k)
       else
       endif
    enddo

    call MPI_ALLREDUCE(range1,temp1,2,MPI_DOUBLE_PRECISION,&
         MPI_MIN,mpicommwd,ierr)
    call MPI_ALLREDUCE(range2,temp2,2,MPI_DOUBLE_PRECISION,&
         MPI_MAX,mpicommwd,ierr)

    epsx = 1.0/(nxlum-3.0)
    epsy = 1.0/(nylum-3.0)

    xmin = temp1(1)-epsx*(temp2(1)-temp1(1))
    xmax = temp2(1)+epsx*(temp2(1)-temp1(1))
    ymin = temp1(2)-epsy*(temp2(2)-temp1(2))
    ymax = temp2(2)+epsy*(temp2(2)-temp1(2))

    hx = (xmax - xmin)/(nxlum - 1.0)
    hy = (ymax - ymin)/(nylum - 1.0)
    hxi = 1.0/hx
    hyi = 1.0/hy

    !//get charge density on the grid using a CIC deposition scheme.
    rholc = 0.0
    do i = 1, nptlc
       ix=int((Pts(1,i)-xmin)*hxi) + 1
       ab=((xmin-Pts(1,i))+ix*hx)*hxi
       jx=int((Pts(3,i)-ymin)*hyi) + 1
       cd=((ymin-Pts(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.inny)) &
       !then
       !   innp
       !   stop
       !endif
       ! (i,j):
       rholc(ix,jx) = rholc(ix,jx) + ab*cd
       ! (i,j+1):
       rholc(ix,jx1) = rholc(ix,jx1) + ab*(1.0-cd)
       ! (i+1,j):
       rholc(ix1,jx) = rholc(ix1,jx)+(1.0-ab)*cd
       ! (i+1,j+1):
       rholc(ix1,jx1) = rholc(ix1,jx1)+(1.0-ab)*(1.0-cd)
    enddo

    !//separate out two beam density
    if(myidy.lt.halfnpy) then
       rholc1 = rholc*Bcurr/Np
       rholc2 = 0.0
    else
       rholc1 = 0.0
       rholc2 = rholc*Bcurr/Np
    endif

    ngrid = nxlum*nylum

    call MPI_REDUCE(rholc1,rho1,ngrid,MPI_DOUBLE_PRECISION,&
         MPI_SUM,0,mpicommwd,ierr)
    call MPI_REDUCE(rholc2,rho2,ngrid,MPI_DOUBLE_PRECISION,&
         MPI_SUM,0,mpicommwd,ierr)

    !//calculate the luminosity from density overlaping.
    if(myid.eq.0) then
       lum = 0.0
       do j = 1, nylum
          do i = 1, nxlum
             lum = lum + rho1(i,j)*rho2(i,j)
          enddo
       enddo
       !          open(1,file="luminosity2d.data",status="unknown",position="append")
       !write(1,1000)it*1.0,lum/(hx*hy)
       nfilelum = 14+mygroupid
       write(nfilelum,1000)it*1.0,lum/(hx*hy)
       call flush(nfilelum)
       !            lum*bcurr1/Np1*bcurr2/Np2/(hx*hy)
       !          close(1)
       !                          (xmax-xmin),(ymax-ymin)
    endif
1000 format(2(1x,es14.7))

  end subroutine luminosity2Gg_Output



  !//caculate the 2d luminosity using 2 group processors.
  subroutine luminosity2Gg2_Output(Pts,nptlc,nxlum,nylum,it,&
       bcurr,Np,myidy,halfnpy,mygroupid,lum2,comm2dlc)
    implicit none
    include 'mpif.h'
    integer :: nptlc,nxlum,nylum,it,Np,myidy,halfnpy,mygroupid,comm2dlc
    double precision :: bcurr
    double precision, pointer, dimension(:,:) :: Pts
    double precision, dimension(2) :: range1,range2
    double precision, dimension(2) :: temp1,temp2
    real*8 :: lum2
    double precision :: xmin,xmax,ymin,ymax,epsx,epsy,hx,hy,hxi,hyi
    double precision :: ab,cd,lum
    double precision, dimension(nxlum,nylum) :: rho1,rho2,rholc,rholc1,rholc2
    integer :: i,j,k,ix,jx,ix1,jx1,ngrid,myid,ierr
    !integer :: nfilelum

    !call MPI_COMM_RANK(mpicommwd,myid,ierr)
    call MPI_COMM_RANK(comm2dlc,myid,ierr)

    range1(1) = Pts(1,1) !minimum x
    range1(2) = Pts(3,1) !minimum y
    range2(1) = Pts(1,1) !maximum x
    range2(2) = Pts(3,1) !maximum y
    do k = 1, nptlc
       if(range1(1).gt.Pts(1,k)) then
          range1(1) = Pts(1,k)
       else if(range2(1).le.Pts(1,k)) then
          range2(1) = Pts(1,k)
       else
       endif
       if(range1(2).gt.Pts(3,k)) then
          range1(2) = Pts(3,k)
       else if(range2(2).le.Pts(3,k)) then
          range2(2) = Pts(3,k)
       else
       endif
    enddo

    !call MPI_ALLREDUCE(range1,temp1,2,MPI_DOUBLE_PRECISION,&
    !                 MPI_MIN,mpicommwd,ierr)
    !call MPI_ALLREDUCE(range2,temp2,2,MPI_DOUBLE_PRECISION,&
    !                 MPI_MAX,mpicommwd,ierr)
    call MPI_ALLREDUCE(range1,temp1,2,MPI_DOUBLE_PRECISION,&
         MPI_MIN,comm2dlc,ierr)
    call MPI_ALLREDUCE(range2,temp2,2,MPI_DOUBLE_PRECISION,&
         MPI_MAX,comm2dlc,ierr)

    epsx = 1.0/(nxlum-3.0)
    epsy = 1.0/(nylum-3.0)

    xmin = temp1(1)-epsx*(temp2(1)-temp1(1))
    xmax = temp2(1)+epsx*(temp2(1)-temp1(1))
    ymin = temp1(2)-epsy*(temp2(2)-temp1(2))
    ymax = temp2(2)+epsy*(temp2(2)-temp1(2))

    hx = (xmax - xmin)/(nxlum - 1.0)
    hy = (ymax - ymin)/(nylum - 1.0)
    hxi = 1.0/hx
    hyi = 1.0/hy

    !//get charge density on the grid using a CIC deposition scheme.
    rholc = 0.0
    do i = 1, nptlc
       ix=int((Pts(1,i)-xmin)*hxi) + 1
       ab=((xmin-Pts(1,i))+ix*hx)*hxi
       jx=int((Pts(3,i)-ymin)*hyi) + 1
       cd=((ymin-Pts(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.inny)) &
       !then
       !   innp
       !   stop
       !endif
       ! (i,j):
       rholc(ix,jx) = rholc(ix,jx) + ab*cd
       ! (i,j+1):
       rholc(ix,jx1) = rholc(ix,jx1) + ab*(1.0-cd)
       ! (i+1,j):
       rholc(ix1,jx) = rholc(ix1,jx)+(1.0-ab)*cd
       ! (i+1,j+1):
       rholc(ix1,jx1) = rholc(ix1,jx1)+(1.0-ab)*(1.0-cd)
    enddo

    !//separate out two beam density
    if(myidy.lt.halfnpy) then
       rholc1 = rholc*Bcurr/Np
       rholc2 = 0.0
    else
       rholc1 = 0.0
       rholc2 = rholc*Bcurr/Np
    endif

    ngrid = nxlum*nylum

    !call MPI_REDUCE(rholc1,rho1,ngrid,MPI_DOUBLE_PRECISION,&
    !                MPI_SUM,0,mpicommwd,ierr)
    !call MPI_REDUCE(rholc2,rho2,ngrid,MPI_DOUBLE_PRECISION,&
    !                MPI_SUM,0,mpicommwd,ierr)
    call MPI_REDUCE(rholc1,rho1,ngrid,MPI_DOUBLE_PRECISION,&
         MPI_SUM,0,comm2dlc,ierr)
    call MPI_REDUCE(rholc2,rho2,ngrid,MPI_DOUBLE_PRECISION,&
         MPI_SUM,0,comm2dlc,ierr)

    lum2 = 0.0d0
    !//calculate the luminosity from density overlaping.
    if(myid.eq.0) then
       lum = 0.0
       do j = 1, nylum
          do i = 1, nxlum
             lum = lum + rho1(i,j)*rho2(i,j)
          enddo
       enddo
       lum2 = lum
       !          open(1,file="luminosity2d.data",status="unknown",position="append")
       !write(1,1000)it*1.0,lum/(hx*hy)
       !          nfilelum = 10+mygroupid
       !          write(nfilelum,1000)it*1.0,lum/(hx*hy)
       !          call flush(nfilelum)
       !            lum*bcurr1/Np1*bcurr2/Np2/(hx*hy)
       !          close(1)
       !                          (xmax-xmin),(ymax-ymin)
    endif
!1000 format(2(1x,es14.7))

  end subroutine luminosity2Gg2_Output



  subroutine xySeries_Output(Bpts,myid,iturn,idbunch,Nplocal,npOut,repRate,lenSeries,nTurns,sampleInd)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: myid,iturn,idbunch,Nplocal,npOut,repRate,lenSeries,nTurns
    double precision, dimension(:,:) :: Bpts
    double precision, allocatable, dimension(:) :: X, Y
    integer :: ii,sampleSize,unit,ierr
    character(len=4) :: pid
    character(len=10) :: filename
    integer, dimension(Nplocal) :: ipt
    integer, dimension(:), pointer, intent(inout) :: sampleInd

    if(mod(iturn,repRate)<lenSeries)then
       if(npOut<Nplocal)then
          sampleSize=npOut
       else
          sampleSize=Nplocal
       endif

       write(pid,'(I4.4)') myid
       filename = 'xy_'//pid
       open(newunit=unit,file=filename,status="unknown",position="append",form="unformatted",&
            iostat=ierr)
       if(ierr==0) then
          allocate(X(npOut))
          allocate(Y(npOut))
          X(:) = Bpts(1,sampleInd(:))
          Y(:) = Bpts(3,sampleInd(:))
          
          if(iturn==0) then
             write(unit)lenSeries,npOut,(nTurns-lenSeries)/repRate+1
          endif
          write(unit)iturn,idbunch
          write(unit)X
          write(unit)Y

          close(unit,iostat=ierr)
          if(ierr/=0) print*,"Warning: Failed to close ",filename,". Expect trouble."
          deallocate(X)
          deallocate(Y)
       else
          print*,"Warning: ",filename," could not be opened. Skipped outputing xy data."
       end if
    end if
  end subroutine xyseries_Output


end module Outputclass
