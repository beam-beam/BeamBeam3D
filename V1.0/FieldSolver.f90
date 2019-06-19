module FieldSolverclass
  use CompDomclass
  use Pgrid2dclass
  use FFTclass
  use Timerclass
  use CompDomclass
  use Transposeclass
  !include 'mpif.h'

contains

  !----------------------------------------------------------------------
  ! update potential (solving Possion's equation) with 2D isolated 
  ! boundary conditions using 2 group PEs and inputed Green function.
  subroutine fieldsolver2d(source,nxlc,nylc,nprocrow,&
       nproccol,nytot,nxpylc2,xpystable,pytable,hx,hy,myidy,commcol,grn)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: nxlc,nylc,nproccol,nytot,&
         myidy,commcol,nprocrow,nxpylc2
    double precision,dimension(nxlc,nylc),intent(inout) :: source
    integer, dimension(0:nproccol-1),intent(in) :: xpystable,pytable
    double complex, intent (in), &
         dimension (2*nytot,nxpylc2) :: grn
    double precision,intent(in) :: hx,hy
    double complex, dimension(2*nytot,nxpylc2) :: rho2out
    integer :: n1,n2
    integer :: i,j !,ierr
    double precision :: t0,scalex,scaley

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    n1 = 2*nxlc
    n2 = 2*nytot

    scalex = 1.0d0
    scaley = 1.0d0

    call fft2d_FFT(n1,n2,nylc,nxpylc2,1,scalex,scaley,source,&
         xpystable,pytable,nproccol,commcol,myidy,rho2out)


    ! multiply transformed charge density and transformed Green 
    ! function:
    do j = 1, nxpylc2
       do i = 1, n2
          rho2out(i,j) = rho2out(i,j)*grn(i,j)
       enddo
    enddo

    ! inverse FFT:
    scalex = 1.0d0/dble(n1)
    scaley = 1.0d0/dble(n2)
    call invfft2d_FFT(n2,n1,nxpylc2,nylc,-1,scalex,scaley,rho2out,&
         pytable,xpystable,nproccol,commcol,myidy,source)


    !call MPI_BARRIER(mpicommwd,ierr)
    t_field = t_field + elapsedtime_Timer(t0)

  end subroutine fieldsolver2d


  !----------------------------------------------------------------------
  ! update potential (solving Possion's equation) with 2D isolated 
  ! boundary conditions for 1 slice model.
  subroutine fieldsolver2d1slc(source,nxlc,nylc,nprocrow,&
       nproccol,nylcr,nytot,LocalTable,hx,hy,myidy,commcol,shift)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: nxlc,nylc,nproccol,nylcr,nytot,&
         myidy,commcol,nprocrow
    double precision,dimension(nxlc,nylc),intent(inout) :: source
    integer,dimension(2,0:nprocrow-1,0:nproccol-1),intent(in)::LocalTable
    double precision, dimension(2), intent(in) :: shift
    double precision,intent(in) :: hx,hy
    integer, dimension(0:nproccol-1) :: xpystable,pytable
    integer :: nxpylc2
    integer :: nsxy1,nsxy2
    integer :: i,inxglb,inyglb,innx,inny
    integer :: nphalf,npbc !,ierr
    double precision :: t0

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    nphalf = nproccol
    if(myidy.lt.nphalf) then
       npbc = 0
    else
       npbc = nphalf
    endif
    do i = 0, nphalf-1
       pytable(i) = LocalTable(2,0,i)
    enddo

    inxglb = nxlc 
    inyglb = nytot
    inny = nylcr
    innx = nxlc

    ! +1 is from the real to complex fft.
    nsxy1 = (inxglb+1)/nphalf
    nsxy2 = (inxglb+1) - nphalf*nsxy1
    do i = 0, nphalf-1
       if(i.le.(nsxy2-1)) then
          xpystable(i) = nsxy1+1
       else
          xpystable(i) = nsxy1
       endif
    enddo

    nxpylc2 = xpystable(myidy-npbc)

    ! Open boundary conditions!
    call openBC2D(innx,inny,source,hx,hy,&
         nxpylc2,myidy,nphalf,commcol,pytable,xpystable,inxglb,inyglb,&
         shift)

    !call MPI_BARRIER(mpicommwd,ierr)
    t_field = t_field + elapsedtime_Timer(t0)

  end subroutine fieldsolver2d1slc


  ! Solving Poisson's equation with open BCs.
  subroutine openBC2D(innx,inny,rho,hx,hy,&
       nxpylc2,myidy,npy,commcol,pytable,xpystable,inxglb,inyglb,&
       shift)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innx,inny,inxglb,inyglb
    integer, intent(in) :: nxpylc2
    integer, intent(in) :: myidy,npy,commcol
    integer, dimension(0:npy-1), intent(in) :: pytable,xpystable
    double precision, dimension(innx,inny), intent(inout) :: rho
    double precision, intent(in) :: hx, hy
    double precision, dimension(2),intent(in) :: shift
    double precision :: scalex,scaley
    integer :: i,j,n1,n2,nylc22
    double complex, dimension(2*inyglb,nxpylc2) :: rho2out,grn
    integer :: ginny

    n1 = 2*inxglb
    n2 = 2*inyglb
    nylc22 = nxpylc2
    scalex = 1.0d0
    scaley = 1.0d0

    call fft2d_FFT(n1,n2,inny,nylc22,1,scalex,scaley,rho,&
         xpystable,pytable,npy,commcol,myidy,rho2out)

    ! compute FFT of the Green function on the grid:
    ! here the +1 is from the unsymmetry of green function
    ! on double-sized grid.
    if(myidy.eq.(npy-1)) then
       ginny = inny + 1
    else
       ginny = inny
    endif
    !call greenf2d(inxglb,inyglb,ginny,nylc22, &
    !       hx,hy,myidy,npy,commcol,xpystable,pytable,grn)
    !call greenf2d(inxglb,inyglb,nylc22,myidy,npy,xpystable,hx,hy,grn,&
    !              shift)
    call greenf2d(inxglb,inyglb,nylc22,hx,hy,myidy,npy,commcol,&
         xpystable,grn,shift)

    ! multiply transformed charge density and transformed Green 
    ! function:
    do j = 1, nylc22
       do i = 1, n2
          rho2out(i,j) = rho2out(i,j)*grn(i,j)
       enddo
    enddo

    ! inverse FFT:
    scalex = 1.0d0/dble(n1)
    scaley = 1.0d0/dble(n2)
    call invfft2d_FFT(n2,n1,nylc22,inny,-1,scalex,scaley,rho2out,&
         pytable,xpystable,npy,commcol,myidy,rho)

  end subroutine openBC2D

  ! parallel integrated green function phi for extended array. (2 group PEs)
  subroutine greenf2d(nx,ny,nsizexy,&
       hx,hy,myidy,npy,commcol,ystable,grnout,shift)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: nx,ny,nsizexy,commcol
    integer, intent(in) :: myidy,npy
    integer, dimension(0:npy-1),intent(in) :: ystable
    double precision, intent(in) :: hx, hy
    double precision, dimension(2), intent(in) :: shift
    double complex, intent (out), &
         dimension (2*ny,nsizexy) :: grnout
    integer :: i,j,iii,jjj,n1,n2
    integer :: nblocky,ksign,nx1,ny1
    !        double precision, dimension (nx,ny) :: grn
    double precision :: scalex,scaley
    double precision :: t0
    !        double precision, dimension(2*nx,2*ny) :: tmp1
    double precision, dimension(2*nx,2*ny/npy) :: tmp2
    double complex, dimension(nx+1,2*ny/npy) :: tmp10
    double precision :: xcent,ycent,rr1
    integer :: flag0,inny,npbc,jgb
    integer, dimension(0:npy-1) :: yrtable
    double precision :: xtmp1,ytmp1,fgrn1,fgrn2,fgrn3,fgrn4,xoy,&
         hxhalf,hyhalf,hxi,hyi

    call starttime_Timer(t0)

    if(myidy.lt.npy) then
       npbc = 0
    else
       npbc = npy
    endif
    inny = 2*ny/npy
    nblocky = (myidy-npbc)*inny 
    n1 = 2*nx
    n2 = 2*ny
    nx1 = n1/2 + 1
    ny1 = n2/2 + 1

    flag0 = 0
    xcent = shift(1)
    ycent = shift(2)
    !        xcent = 0.0
    !        ycent = 0.0
    hxhalf = hx/2
    hyhalf = hy/2
    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    do j = 1, inny
       jgb = j + nblocky
       if(jgb.lt.ny1) then
          do i = 1, nx
             jjj = jgb - 1
             iii = i - 1
             xtmp1 = xcent+hx*iii+hxhalf
             ytmp1 = ycent+hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn1 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent+hx*iii-hxhalf
             ytmp1 = ycent+hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn2 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent+hx*iii+hxhalf
             ytmp1 = ycent+hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn3 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent+hx*iii-hxhalf
             ytmp1 = ycent+hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn4 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             tmp2(i,j) =-(fgrn1-fgrn2-fgrn3+fgrn4)*hxi*hyi/2
          enddo
          do i = nx1,n1
             jjj = jgb - 1
             iii = n1-i+1
             xtmp1 = xcent-hx*iii+hxhalf
             ytmp1 = ycent+hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn1 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent-hx*iii-hxhalf
             ytmp1 = ycent+hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn2 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent-hx*iii+hxhalf
             ytmp1 = ycent+hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn3 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent-hx*iii-hxhalf
             ytmp1 = ycent+hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn4 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             tmp2(i,j) =-(fgrn1-fgrn2-fgrn3+fgrn4)*hxi*hyi/2
          enddo
       else
          do i = 1, nx
             jjj = n2-jgb+1
             iii = i - 1
             xtmp1 = xcent+hx*iii+hxhalf
             ytmp1 = ycent-hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn1 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent+hx*iii-hxhalf
             ytmp1 = ycent-hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn2 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent+hx*iii+hxhalf
             ytmp1 = ycent-hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn3 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent+hx*iii-hxhalf
             ytmp1 = ycent-hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn4 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             tmp2(i,j) =-(fgrn1-fgrn2-fgrn3+fgrn4)*hxi*hyi/2
          enddo
          do i = nx1,n1
             jjj = n2-jgb+1
             iii = n1-i+1
             xtmp1 = xcent-hx*iii+hxhalf
             ytmp1 = ycent-hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn1 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent-hx*iii-hxhalf
             ytmp1 = ycent-hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn2 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent-hx*iii+hxhalf
             ytmp1 = ycent-hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn3 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent-hx*iii-hxhalf
             ytmp1 = ycent-hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn4 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             tmp2(i,j) =-(fgrn1-fgrn2-fgrn3+fgrn4)*hxi*hyi/2
          enddo
       endif
    enddo

    scalex = 1.0d0
    scaley = 1.0d0
    ksign = 1

    ! FFTs along x dimensions:
    call fftrclocal_FFT(ksign,scalex,tmp2,n1,inny,tmp10)

    ! FFTs along y dimensions:
    yrtable = n2/npy
    call trans2d2_TRANSP(nx1,n2,nsizexy,inny,tmp10,grnout,npy,&
         ystable,yrtable,commcol,myidy)

    call fftlocal_FFT(ksign,scaley,grnout,n2,nsizexy)

    t_greenf = t_greenf + elapsedtime_Timer(t0)

  end subroutine greenf2d


  ! FFT for 2D open boundary conditions. 
  ! The original computational domain is doubled in each dimension
  ! to apply the FFT for the new domain.
  ! 2_D FFT.
  subroutine fft2d_FFT(nx,ny,nsizey,nsizexy,ksign,scalex,scaley,x,&
       ystable,yrtable,nproccol,commcol,myidy,xout)
    implicit none
    !include 'mpif.h'
    integer,intent(in) :: nx,ny,nsizey,nsizexy
    integer,intent(in) :: nproccol,commcol
    integer,intent(in) :: ksign,myidy
    double precision, intent(in) :: scalex,scaley
    integer,dimension(0:nproccol-1),intent(in) :: ystable,yrtable
    double precision, dimension(nx/2,nsizey), intent(in) &
         :: x
    double complex, dimension(ny,nsizexy), intent(out) &
         :: xout
    double precision, dimension(nx,nsizey) :: tmp1
    !        double complex, dimension(nx/2+1,nsizey) :: tmp10
    double complex, dimension(nx/2+1,nsizey) :: x0
    !        double complex, dimension(ny,nsizexy) :: tmp2
    double complex, dimension(ny/2,nsizexy) :: x1
    integer :: i,j
    double precision :: t0
    integer :: nxx !,ierr

    call starttime_Timer(t0)

    nxx = nx/2 + 1
    !FFTs along x dimensions: could be a lot of cache miss.
    do j = 1, nsizey 
       do i = 1, nx/2
          tmp1(i,j) = x(i,j)
       enddo
       do i = nx/2+1, nx
          tmp1(i,j) = 0.0
       enddo
    enddo

    ! FFTs along x dimensions:
    call fftrclocal_FFT(ksign,scalex,tmp1,nx,nsizey,x0)

    !        do j = 1, nsizey 
    !          do i = 1, nx/2+1
    !            x0(i,j) = tmp10(i,j)
    !          enddo
    !        enddo

    ! FFTs along y dimensions:
    !        call MPI_BARRIER(commcol,ierr)
    ! yrtable needs to be changed.
    ! transpose between the x and y dimension.
    !        call trans2dnew_TRANSP(nxx,ny/2,nsizexy,nsizey,x0,x1,nproccol,&
    !                     ystable,yrtable,commcol)
    call trans2d2_TRANSP(nxx,ny/2,nsizexy,nsizey,x0,x1,nproccol,&
         ystable,yrtable,commcol,myidy)

    do j = 1, nsizexy 
       do i = 1, ny/2
          xout(i,j) = x1(i,j)
       enddo
       do i = ny/2+1,ny
          xout(i,j) = (0.0,0.0)
       enddo
    enddo

    call fftlocal_FFT(ksign,scaley,xout,ny,nsizexy)

    !        do j = 1, nsizexy 
    !          do i = 1, ny
    !            xout(i,j) = tmp2(i,j) 
    !          enddo
    !        enddo

    t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

    return
  end subroutine fft2d_FFT



  ! 2_D inverse FFT for open BCs.
  subroutine invfft2d_FFT(ny,nx,nsizexy,nsizey,ksign,scalex,scaley,x,&
       ystable,yrtable,nproccol,commcol,&
       myidy,xout)
    implicit none
    !include 'mpif.h'
    integer,intent(in) :: nx,ny,nsizey,nsizexy
    integer,intent(in) :: nproccol,commcol
    integer,intent(in) :: ksign,myidy
    double precision, intent(in) :: scalex,scaley
    integer,dimension(0:nproccol-1),intent(in) :: ystable,yrtable
    double complex, dimension(ny,nsizexy), intent(inout) &
         :: x
    double precision, dimension(nx/2,nsizey), intent(out) &
         :: xout
    !        double complex, dimension(ny,nsizexy) :: tmp2
    !        double complex, dimension(nx/2+1,nsizey) :: tmp3
    double precision, dimension(nx,nsizey) :: tmp30
    double complex, dimension(ny/2,nsizexy) :: x0
    double complex, dimension(nx/2+1,nsizey) :: x1
    integer :: i,j
    double precision :: t0
    integer :: nxx !,ierr
    !double complex, allocatable, dimension(:,:) :: x0
    !double complex, allocatable, dimension(:,:) :: x1

    call starttime_Timer(t0)

    nxx = nx/2 + 1

    !        do j = 1, nsizexy 
    !          do i = 1, ny
    !            tmp2(i,j) = x(i,j)
    !          enddo
    !        enddo

    call fftlocal_FFT(ksign,scaley,x,ny,nsizexy)

    do j = 1, nsizexy 
       do i = 1, ny/2
          x0(i,j) = x(i,j) 
          !            x0(i,j,k) = tmp2(i,j)*scaley 
       enddo
    enddo

    !        call trans2dnew_TRANSP(ny/2,nxx,nsizey,nsizexy,x0,x1,nproccol,&
    !                     ystable,yrtable,commcol)
    call trans2d2_TRANSP(ny/2,nxx,nsizey,nsizexy,x0,x1,nproccol,&
         ystable,yrtable,commcol,myidy)
    !        call MPI_BARRIER(commcol,ierr)

    !        do j = 1, nsizey
    !          do i = 1, nxx
    !            tmp3(i,j) = x1(i,j) 
    !          enddo
    !        enddo
    call fftcrlocal_FFT(ksign,scalex,x1,nx,nsizey,tmp30)

    do j = 1, nsizey
       do i = 1, nx/2
          xout(i,j) = tmp30(i,j)
          !            xout(i,j) = tmp30(i,j)*scalex*2
       enddo
    enddo

    t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

    return

  end subroutine invfft2d_FFT


  !----------------------------------------------------------------------------
  ! update potential (solving Possion's equation) with 2D isolated 
  ! boundary conditions using 1 group PEs.
  ! update potential (solving Possion's equation) with 3D isolated 
  ! boundary conditions for strong beam only.
  subroutine fieldsolver2d1Gwk(source,nxlc,nylc,nprocrow,&
       nproccol,nytot,nxpylc2,xpystable,pytable,hx,hy,myidy,commcol,grn)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: nxlc,nylc,nproccol,nytot,&
         myidy,commcol,nprocrow,nxpylc2
    double precision,dimension(nxlc,nylc),intent(inout) :: source
    integer, dimension(0:nproccol-1),intent(in) :: xpystable,pytable
    double complex, intent (in), &
         dimension (2*nytot,nxpylc2) :: grn
    double precision,intent(in) :: hx,hy
    double complex, dimension(2*nytot,nxpylc2) :: rho2out
    integer :: n1,n2
    integer :: i,j !,ierr
    double precision :: t0,scalex,scaley

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    n1 = 2*nxlc
    n2 = 2*nytot

    scalex = 1.0d0
    scaley = 1.0d0

    call fft2d1G_FFT(n1,n2,nylc,nxpylc2,1,scalex,scaley,source,&
         xpystable,pytable,nproccol,commcol,myidy,rho2out)


    ! multiply transformed charge density and transformed Green 
    ! function:
    do j = 1, nxpylc2
       do i = 1, n2
          rho2out(i,j) = rho2out(i,j)*grn(i,j)
       enddo
    enddo

    ! inverse FFT:
    scalex = 1.0d0/dble(n1)
    scaley = 1.0d0/dble(n2)
    call invfft2d1G_FFT(n2,n1,nxpylc2,nylc,-1,scalex,scaley,rho2out,&
         pytable,xpystable,nproccol,commcol,myidy,source)


    !call MPI_BARRIER(mpicommwd,ierr)
    t_field = t_field + elapsedtime_Timer(t0)

  end subroutine fieldsolver2d1Gwk



  subroutine fieldsolver2d1G(source,nxlc,nylc,nprocrow,&
       nproccol,nylcr,nytot,LocalTable,hx,hy,myidy,commcol,shift)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: nxlc,nylc,nproccol,nylcr,nytot,&
         myidy,commcol,nprocrow
    double precision,dimension(nxlc,nylc),intent(inout) :: source
    integer,dimension(2,0:nprocrow-1,0:nproccol-1),intent(in)::LocalTable
    double precision, dimension(2), intent(in) :: shift
    double precision,intent(in) :: hx,hy
    integer, dimension(0:nproccol-1) :: xpystable,pytable
    integer :: nxpylc2,nsxy1,nsxy2
    integer :: i,inxglb,inyglb,innx,inny
    !integer :: ierr
    double precision :: t0

    !call MPI_BARRIER(mpicommwd,ierr)
    call starttime_Timer(t0)

    do i = 0, nproccol-1
       pytable(i) = LocalTable(2,0,i)
    enddo

    inxglb = nxlc 
    inyglb = nytot
    inny = nylcr
    innx = nxlc

    ! +1 is from the real to complex fft.
    nsxy1 = (inxglb+1)/nproccol
    nsxy2 = (inxglb+1) - nproccol*nsxy1
    do i = 0, nproccol-1
       if(i.le.(nsxy2-1)) then
          xpystable(i) = nsxy1+1
       else
          xpystable(i) = nsxy1
       endif
    enddo

    nxpylc2 = xpystable(myidy)

    ! Open boundary conditions!
    call openBC2D1G(innx,inny,source,hx,hy,&
         nxpylc2,myidy,nproccol,commcol,pytable,xpystable,inxglb,inyglb,&
         shift)

    !call MPI_BARRIER(mpicommwd,ierr)
    t_field = t_field + elapsedtime_Timer(t0)

  end subroutine fieldsolver2d1G



  !----------------------------------------------------------------------
  ! Solving Poisson's equation with open BCs.
  subroutine openBC2D1G(innx,inny,rho,hx,hy,&
       nxpylc2,myidy,npy,commcol,pytable,xpystable,inxglb,inyglb,&
       shift)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innx,inny,inxglb,inyglb
    integer, intent(in) :: nxpylc2
    integer, intent(in) :: myidy,npy,commcol
    integer, dimension(0:npy-1), intent(in) :: pytable,xpystable
    double precision, dimension(innx,inny), intent(inout) :: rho
    double precision, intent(in) :: hx, hy
    double precision, dimension(2),intent(in) :: shift
    double precision :: scalex,scaley
    integer :: i,j,n1,n2,nylc22
    double complex, dimension(2*inyglb,nxpylc2) :: rho2out,grn
    integer :: ginny

    n1 = 2*inxglb
    n2 = 2*inyglb

    nylc22 = nxpylc2

    scalex = 1.0d0
    scaley = 1.0d0

    call fft2d1G_FFT(n1,n2,inny,nylc22,1,scalex,scaley,rho,&
         xpystable,pytable,npy,commcol,myidy,rho2out)

    !c compute FFT of the Green function on the grid:
    ! here the +1 is from the unsymmetry of green function
    ! on double-sized grid.
    if(myidy.eq.(npy-1)) then
       ginny = inny + 1
    else
       ginny = inny
    endif
    !call greenf2d1G(inxglb,inyglb,ginny,nylc22, &
    !       hx,hy,myidy,npy,commcol,xpystable,pytable,grn)
    !call greenf2d1G(inxglb,inyglb,nylc22,myidy,npy,xpystable,hx,hy,grn,&
    !              shift)
    !using parallel integrated green function.
    call greenf2d1G(inxglb,inyglb,nylc22,hx,hy,myidy,npy,commcol,&
         xpystable,grn,shift)

    ! multiply transformed charge density and transformed Green 
    ! function:
    do j = 1, nylc22
       do i = 1, n2
          rho2out(i,j) = rho2out(i,j)*grn(i,j)
       enddo
    enddo

    ! inverse FFT:
    scalex = 1.0d0/dble(n1)
    scaley = 1.0d0/dble(n2)
    call invfft2d1G_FFT(n2,n1,nylc22,inny,-1,scalex,scaley,rho2out,&
         pytable,xpystable,npy,commcol,myidy,rho)

  end subroutine openBC2D1G



  ! parallel integrated green function for extended array.
  subroutine greenf2d1G(nx,ny,nsizexy,&
       hx,hy,myidy,npy,commcol,ystable,grnout,shift)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: nx,ny,nsizexy,commcol
    integer, intent(in) :: myidy,npy
    integer, dimension(0:npy-1),intent(in) :: ystable
    double precision, intent(in) :: hx, hy
    double precision, dimension(2), intent(in) :: shift
    double complex, intent (out), &
         dimension (2*ny,nsizexy) :: grnout
    integer :: i,j,iii,jjj,n1,n2
    integer :: nblocky,ksign,nx1,ny1
    double precision :: scalex,scaley
    double precision :: t0
    double precision, dimension(2*nx,2*ny/npy) :: tmp2
    double complex, dimension(nx+1,2*ny/npy) :: tmp10
    double precision :: xcent,ycent,rr1
    integer :: flag0,inny,npbc,jgb
    integer, dimension(0:npy-1) :: yrtable
    double precision :: xtmp1,ytmp1,fgrn1,fgrn2,fgrn3,fgrn4,xoy,&
         hxhalf,hyhalf,hxi,hyi

    call starttime_Timer(t0)

    npbc = 0
    inny = 2*ny/npy
    nblocky = (myidy-npbc)*inny 
    n1 = 2*nx
    n2 = 2*ny
    nx1 = n1/2 + 1
    ny1 = n2/2 + 1

    flag0 = 0
    xcent = shift(1)
    ycent = shift(2)
    hxhalf = hx/2
    hyhalf = hy/2
    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    do j = 1, inny
       jgb = j + nblocky
       if(jgb.lt.ny1) then
          do i = 1, nx
             jjj = jgb - 1
             iii = i - 1
             xtmp1 = xcent+hx*iii+hxhalf
             ytmp1 = ycent+hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn1 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent+hx*iii-hxhalf
             ytmp1 = ycent+hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn2 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent+hx*iii+hxhalf
             ytmp1 = ycent+hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn3 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent+hx*iii-hxhalf
             ytmp1 = ycent+hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn4 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             tmp2(i,j) =-(fgrn1-fgrn2-fgrn3+fgrn4)*hxi*hyi/2
          enddo
          do i = nx1,n1
             jjj = jgb - 1
             iii = n1-i+1
             xtmp1 = xcent-hx*iii+hxhalf
             ytmp1 = ycent+hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn1 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent-hx*iii-hxhalf
             ytmp1 = ycent+hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn2 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent-hx*iii+hxhalf
             ytmp1 = ycent+hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn3 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent-hx*iii-hxhalf
             ytmp1 = ycent+hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn4 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             tmp2(i,j) =-(fgrn1-fgrn2-fgrn3+fgrn4)*hxi*hyi/2
          enddo
       else
          do i = 1, nx
             jjj = n2-jgb+1
             iii = i - 1
             xtmp1 = xcent+hx*iii+hxhalf
             ytmp1 = ycent-hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn1 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent+hx*iii-hxhalf
             ytmp1 = ycent-hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn2 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent+hx*iii+hxhalf
             ytmp1 = ycent-hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn3 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent+hx*iii-hxhalf
             ytmp1 = ycent-hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn4 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             tmp2(i,j) =-(fgrn1-fgrn2-fgrn3+fgrn4)*hxi*hyi/2
          enddo
          do i = nx1,n1
             jjj = n2-jgb+1
             iii = n1-i+1
             xtmp1 = xcent-hx*iii+hxhalf
             ytmp1 = ycent-hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn1 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent-hx*iii-hxhalf
             ytmp1 = ycent-hy*jjj+hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn2 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent-hx*iii+hxhalf
             ytmp1 = ycent-hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn3 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             xtmp1 = xcent-hx*iii-hxhalf
             ytmp1 = ycent-hy*jjj-hyhalf
             rr1 = xtmp1**2 + ytmp1**2
             xoy = xtmp1/ytmp1
             fgrn4 = -3*xtmp1*ytmp1+xtmp1**2*atan(1.0d0/xoy)+ &
                  ytmp1**2*atan(xoy)+xtmp1*ytmp1*log(rr1)
             tmp2(i,j) =-(fgrn1-fgrn2-fgrn3+fgrn4)*hxi*hyi/2
          enddo
       endif
    enddo

    scalex = 1.0d0
    scaley = 1.0d0
    ksign = 1

    ! FFTs along x dimensions:
    call fftrclocal_FFT(ksign,scalex,tmp2,n1,inny,tmp10)

    ! FFTs along y dimensions:
    yrtable = n2/npy
    call trans2d21G_TRANSP(nx1,n2,nsizexy,inny,tmp10,grnout,npy,&
         ystable,yrtable,commcol,myidy)

    call fftlocal_FFT(ksign,scaley,grnout,n2,nsizexy)

    t_greenf = t_greenf + elapsedtime_Timer(t0)

  end subroutine greenf2d1G

  ! FFT for 2D open boundary conditions. 
  ! The original computational domain is doubled in each dimension
  ! to apply the FFT for the new domain.
  ! 2_D FFT.
  subroutine fft2d1G_FFT(nx,ny,nsizey,nsizexy,ksign,scalex,scaley,x,&
       ystable,yrtable,nproccol,commcol,myidy,xout)
    implicit none
    !include 'mpif.h'
    integer,intent(in) :: nx,ny,nsizey,nsizexy
    integer,intent(in) :: nproccol,commcol
    integer,intent(in) :: ksign,myidy
    double precision, intent(in) :: scalex,scaley
    integer,dimension(0:nproccol-1),intent(in) :: ystable,yrtable
    double precision, dimension(nx/2,nsizey), intent(in) &
         :: x
    double complex, dimension(ny,nsizexy), intent(out) &
         :: xout
    double precision, dimension(nx,nsizey) :: tmp1
    double complex, dimension(nx/2+1,nsizey) :: tmp10
    double complex, dimension(ny,nsizexy) :: tmp2
    double complex, dimension(ny/2,nsizexy) :: x1
    double complex, dimension(nx/2+1,nsizey) :: x0
    integer :: i,j
    double precision :: t0
    integer :: nxx !,ierr

    call starttime_Timer(t0)

    nxx = nx/2 + 1
    !FFTs along x dimensions: could be a lot of cache miss.
    do j = 1, nsizey 
       do i = 1, nx/2
          tmp1(i,j) = x(i,j)
       enddo
       do i = nx/2+1, nx
          tmp1(i,j) = 0.0
       enddo
    enddo

    ! FFTs along x dimensions:
    call fftrclocal_FFT(ksign,scalex,tmp1,nx,nsizey,tmp10)

    do j = 1, nsizey 
       do i = 1, nx/2+1
          x0(i,j) = tmp10(i,j)
       enddo
    enddo

    ! FFTs along y dimensions:
    !        call MPI_BARRIER(commcol,ierr)
    ! yrtable needs to be changed.
    ! transpose between the x and y dimension.
    !        call trans2dnew_TRANSP(nxx,ny/2,nsizexy,nsizey,x0,x1,nproccol,&
    !                     ystable,yrtable,commcol)
    call trans2d21G_TRANSP(nxx,ny/2,nsizexy,nsizey,x0,x1,nproccol,&
         ystable,yrtable,commcol,myidy)

    do j = 1, nsizexy 
       do i = 1, ny/2
          tmp2(i,j) = x1(i,j)
       enddo
       do i = ny/2+1,ny
          tmp2(i,j) = (0.0,0.0)
       enddo
    enddo

    call fftlocal_FFT(ksign,scaley,tmp2,ny,nsizexy)

    do j = 1, nsizexy 
       do i = 1, ny
          xout(i,j) = tmp2(i,j) 
       enddo
    enddo

    t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

    return
  end subroutine fft2d1G_FFT



  ! 2_D inverse FFT for open BCs.
  subroutine invfft2d1G_FFT(ny,nx,nsizexy,nsizey,ksign,scalex,scaley,x,&
       ystable,yrtable,nproccol,commcol,&
       myidy,xout)
    implicit none
    !include 'mpif.h'
    integer,intent(in) :: nx,ny,nsizey,nsizexy
    integer,intent(in) :: nproccol,commcol
    integer,intent(in) :: ksign,myidy
    double precision, intent(in) :: scalex,scaley
    integer,dimension(0:nproccol-1),intent(in) :: ystable,yrtable
    double complex, dimension(ny,nsizexy), intent(inout) &
         :: x
    double precision, dimension(nx/2,nsizey), intent(out) &
         :: xout
    double complex, dimension(ny,nsizexy) :: tmp2
    double complex, dimension(nx/2+1,nsizey) :: tmp3
    double precision, dimension(nx,nsizey) :: tmp30
    double complex, dimension(ny/2,nsizexy) :: x0
    double complex, dimension(nx/2+1,nsizey) :: x1
    integer :: i,j !,k
    double precision :: t0
    integer :: nxx !,ierr
    !double complex, allocatable, dimension(:,:) :: x0
    !double complex, allocatable, dimension(:,:) :: x1

    call starttime_Timer(t0)

    nxx = nx/2 + 1

    do j = 1, nsizexy 
       do i = 1, ny
          tmp2(i,j) = x(i,j)
       enddo
    enddo

    call fftlocal_FFT(ksign,scaley,tmp2,ny,nsizexy)

    do j = 1, nsizexy 
       do i = 1, ny/2
          x0(i,j) = tmp2(i,j) 
          !            x0(i,j,k) = tmp2(i,j)*scaley 
       enddo
    enddo

    !        call trans2dnew_TRANSP(ny/2,nxx,nsizey,nsizexy,x0,x1,nproccol,&
    !                     ystable,yrtable,commcol)
    call trans2d21G_TRANSP(ny/2,nxx,nsizey,nsizexy,x0,x1,nproccol,&
         ystable,yrtable,commcol,myidy)
    !        call MPI_BARRIER(commcol,ierr)

    do j = 1, nsizey
       do i = 1, nxx
          tmp3(i,j) = x1(i,j) 
       enddo
    enddo

    call fftcrlocal_FFT(ksign,scalex,tmp3,nx,nsizey,tmp30)

    do j = 1, nsizey
       do i = 1, nx/2
          xout(i,j) = tmp30(i,j)
          !            xout(i,j) = tmp30(i,j)*scalex*2
       enddo
    enddo

    t_fft2dhpf = t_fft2dhpf + elapsedtime_Timer(t0)

    return
  end subroutine invfft2d1G_FFT

end module FieldSolverclass
