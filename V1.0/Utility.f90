module Utilityclass
  use Pgrid2dclass

contains

  ! calculate the centroid shift from beam2 to beam 1, 2 group procs.
  subroutine findshift_Utility(Pts1,innp,nptot,myidy,npyhalf,shift)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot,myidy,npyhalf
    double precision, dimension(2), intent(out) :: shift
    double precision, dimension(2) :: center1,center2
    double precision:: x0lc,y0lc
    integer :: i,ierr
    double precision, dimension(4) :: tmplc,tmpgl

    x0lc = 0.0
    y0lc = 0.0

    do i = 1, innp
       x0lc = x0lc + Pts1(1,i)
       y0lc = y0lc + Pts1(3,i)
    enddo
    !        x0lc = x0lc/nptot
    !        y0lc = y0lc/nptot

    !try to obtain the centorid of bunch from each group of processors
    if(myidy.lt.npyhalf) then
       tmplc(1) = x0lc
       tmplc(2) = y0lc
       tmplc(3) = 0.0
       tmplc(4) = 0.0
    else
       tmplc(1) = 0.0
       tmplc(2) = 0.0
       tmplc(3) = x0lc
       tmplc(4) = y0lc
    endif

    call MPI_ALLREDUCE(tmplc,tmpgl,4,MPI_DOUBLE_PRECISION,&
         MPI_SUM,mpicommwd,ierr)

    center1(1) = tmpgl(1)/nptot
    center1(2) = tmpgl(2)/nptot
    center2(1) = tmpgl(3)/nptot
    center2(2) = tmpgl(4)/nptot
    shift = center2 - center1

  end subroutine findshift_Utility



  !//calculate the centroid of one beam
  subroutine findcenter_Utility(Pts1,innp,nptot,center)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot
    double precision, dimension(3), intent(out) :: center
    double precision:: x0lc,y0lc
    integer :: i,ierr
    double precision, dimension(2) :: tmplc

    x0lc = 0.0
    y0lc = 0.0
    center = 0.0

    do i = 1, innp
       x0lc = x0lc + Pts1(1,i)
       y0lc = y0lc + Pts1(3,i)
    enddo
    x0lc = x0lc/nptot
    y0lc = y0lc/nptot

    !try to obtain the centorid of bunch from each group of processors
    tmplc(1) = x0lc
    tmplc(2) = y0lc

    call MPI_ALLREDUCE(tmplc,center,2,MPI_DOUBLE_PRECISION,&
         MPI_SUM,mpicommwd,ierr)

  end subroutine findcenter_Utility



  ! calculate the 6D centroids of both beams
  subroutine findcentroids_Utility(Pts, innp, nptot, myidy, npyhalf, center1, center2)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:), intent(in) :: Pts
    integer, intent(in) :: innp, nptot, myidy, npyhalf
    double precision, dimension(6), intent(out) :: center1, center2
    integer :: i, j, ierr
    double precision, dimension(12) :: tmplc, tmpgl

    tmplc = 0.0 
    tmpgl = 0.0 
    !// assign tmplc(1:6) to coordinates of first beam and tmplc(7:12) to second beam
    if(myidy.lt.npyhalf) then
       do i = 1, innp
          do j = 1, 6
             tmplc(j) = tmplc(j) + Pts(j,i)
          enddo
       enddo
    else
       do i = 1, innp
          do j = 1, 6
             tmplc(j+6) = tmplc(j+6) + Pts(j,i)
          enddo
       enddo
    endif

    call MPI_ALLREDUCE(tmplc, tmpgl, 12, MPI_DOUBLE_PRECISION, MPI_SUM, mpicommwd, ierr)

    center1(1:6) = tmpgl(1:6)/nptot
    center2(1:6) = tmpgl(7:12)/nptot

  end subroutine findcentroids_Utility



  !// Convenience function to print centroids easily
  subroutine printcentroids_Utility(Pts, innp, nptot, myidy, npyhalf)
    double precision, pointer, dimension(:,:), intent(in) :: Pts
    integer, intent(in) :: innp, nptot, myidy, npyhalf

    double precision, dimension(6) :: c1, c2
    integer :: myid, ierr

    call MPI_COMM_RANK(mpicommwd, myid, ierr)

    call findcentroids_Utility(Pts, innp, nptot, myidy, npyhalf, c1, c2)
    if(myid==0) then
       print*, "centroid1:", c1
       print*, "centroid2:", c2
    endif

  end subroutine printcentroids_Utility



  ! calculate the centroid shift from beam2 to beam 1,centroid,sigma
  subroutine shiftsigma_Utility(Pts1,innp,nptot,myidy,npyhalf,&
       shift,center,sigma)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot,myidy,npyhalf
    double precision, dimension(2), intent(out) :: shift
    double precision, dimension(3), intent(out) :: center
    double precision, dimension(6), intent(out) :: sigma
    double precision, dimension(2) :: center1,center2
    double precision:: x0lc,px0lc,y0lc,py0lc,z0lc,pz0lc,xrlc,pxrlc,yrlc,pyrlc, &
         zrlc,pzrlc
    integer :: i,ierr
    double precision, dimension(24) :: tmplc,tmpgl

    x0lc = 0.0
    px0lc = 0.0
    y0lc = 0.0
    py0lc = 0.0
    z0lc = 0.0
    pz0lc = 0.0
    xrlc = 0.0
    pxrlc = 0.0
    yrlc = 0.0
    pyrlc = 0.0
    zrlc = 0.0
    pzrlc = 0.0
    do i = 1, innp
       x0lc = x0lc + Pts1(1,i)
       px0lc = px0lc + Pts1(2,i)
       y0lc = y0lc + Pts1(3,i)
       py0lc = py0lc + Pts1(4,i)
       z0lc = z0lc + Pts1(5,i)
       pz0lc = pz0lc + Pts1(6,i)
       xrlc = xrlc + Pts1(1,i)*Pts1(1,i)
       pxrlc = pxrlc + Pts1(2,i)*Pts1(2,i)
       yrlc = yrlc + Pts1(3,i)*Pts1(3,i)
       pyrlc = pyrlc + Pts1(4,i)*Pts1(4,i)
       zrlc = zrlc + Pts1(5,i)*Pts1(5,i)
       pzrlc = pzrlc + Pts1(6,i)*Pts1(6,i)
    enddo

    !try to obtain the centorid and rms of bunch from each group of processors
    tmplc = 0.0 
    tmpgl = 0.0 
    if(myidy.lt.npyhalf) then
       tmplc(1) = x0lc
       tmplc(2) = px0lc
       tmplc(3) = y0lc
       tmplc(4) = py0lc
       tmplc(5) = z0lc
       tmplc(6) = pz0lc
       tmplc(7) = xrlc
       tmplc(8) = pxrlc
       tmplc(9) = yrlc
       tmplc(10) = pyrlc
       tmplc(11) = zrlc
       tmplc(12) = pzrlc
    else
       tmplc(13) = x0lc
       tmplc(14) = px0lc
       tmplc(15) = y0lc
       tmplc(16) = py0lc
       tmplc(17) = z0lc
       tmplc(18) = pz0lc
       tmplc(19) = xrlc
       tmplc(20) = pxrlc
       tmplc(21) = yrlc
       tmplc(22) = pyrlc
       tmplc(23) = zrlc
       tmplc(24) = pzrlc
    endif

    call MPI_ALLREDUCE(tmplc,tmpgl,24,MPI_DOUBLE_PRECISION,&
         MPI_SUM,mpicommwd,ierr)

    center1(1) = tmpgl(1)/nptot
    center1(2) = tmpgl(3)/nptot
    center2(1) = tmpgl(13)/nptot
    center2(2) = tmpgl(15)/nptot
    shift = center2 - center1

    if(myidy.lt.npyhalf) then
       center(1) = tmpgl(1)/nptot
       center(2) = tmpgl(3)/nptot
       center(3) = tmpgl(5)/nptot
       sigma(1) = sqrt(tmpgl(7)/nptot - (tmpgl(1)/nptot)**2)
       sigma(2) = sqrt(tmpgl(8)/nptot - (tmpgl(2)/nptot)**2)
       sigma(3) = sqrt(tmpgl(9)/nptot - (tmpgl(3)/nptot)**2)
       sigma(4) = sqrt(tmpgl(10)/nptot - (tmpgl(4)/nptot)**2)
       sigma(5) = sqrt(tmpgl(11)/nptot - (tmpgl(5)/nptot)**2)
       sigma(6) = sqrt(tmpgl(12)/nptot - (tmpgl(6)/nptot)**2)
    else
       center(1) = tmpgl(13)/nptot
       center(2) = tmpgl(15)/nptot
       center(3) = tmpgl(17)/nptot
       sigma(1) = sqrt(tmpgl(19)/nptot - (tmpgl(13)/nptot)**2)
       sigma(2) = sqrt(tmpgl(20)/nptot - (tmpgl(14)/nptot)**2)
       sigma(3) = sqrt(tmpgl(21)/nptot - (tmpgl(15)/nptot)**2)
       sigma(4) = sqrt(tmpgl(22)/nptot - (tmpgl(16)/nptot)**2)
       sigma(5) = sqrt(tmpgl(23)/nptot - (tmpgl(17)/nptot)**2)
       sigma(6) = sqrt(tmpgl(24)/nptot - (tmpgl(18)/nptot)**2)
    endif

  end subroutine shiftsigma_Utility


  ! calculate the centroid shift from beam2 to beam 1,centroid,sigma
  subroutine shiftsigma2_Utility(Pts1,innp,nptot,myidy,npyhalf,&
       shift,center,sigma)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot,myidy,npyhalf
    double precision, dimension(2), intent(out) :: shift
    double precision, dimension(6), intent(out) :: center
    double precision, dimension(6), intent(out) :: sigma
    double precision, dimension(2) :: center1,center2
    double precision:: x0lc,px0lc,y0lc,py0lc,z0lc,pz0lc,xrlc,pxrlc,yrlc,pyrlc, &
         zrlc,pzrlc
    integer :: i,ierr
    double precision, dimension(24) :: tmplc,tmpgl

    x0lc = 0.0
    px0lc = 0.0
    y0lc = 0.0
    py0lc = 0.0
    z0lc = 0.0
    pz0lc = 0.0
    xrlc = 0.0
    pxrlc = 0.0
    yrlc = 0.0
    pyrlc = 0.0
    zrlc = 0.0
    pzrlc = 0.0
    do i = 1, innp
       x0lc = x0lc + Pts1(1,i)
       px0lc = px0lc + Pts1(2,i)
       y0lc = y0lc + Pts1(3,i)
       py0lc = py0lc + Pts1(4,i)
       z0lc = z0lc + Pts1(5,i)
       pz0lc = pz0lc + Pts1(6,i)
       xrlc = xrlc + Pts1(1,i)*Pts1(1,i)
       pxrlc = pxrlc + Pts1(2,i)*Pts1(2,i)
       yrlc = yrlc + Pts1(3,i)*Pts1(3,i)
       pyrlc = pyrlc + Pts1(4,i)*Pts1(4,i)
       zrlc = zrlc + Pts1(5,i)*Pts1(5,i)
       pzrlc = pzrlc + Pts1(6,i)*Pts1(6,i)
    enddo

    !try to obtain the centorid and rms of bunch from each group of processors
    tmplc = 0.0 
    tmpgl = 0.0 
    if(myidy.lt.npyhalf) then
       tmplc(1) = x0lc
       tmplc(2) = px0lc
       tmplc(3) = y0lc
       tmplc(4) = py0lc
       tmplc(5) = z0lc
       tmplc(6) = pz0lc
       tmplc(7) = xrlc
       tmplc(8) = pxrlc
       tmplc(9) = yrlc
       tmplc(10) = pyrlc
       tmplc(11) = zrlc
       tmplc(12) = pzrlc
    else
       tmplc(13) = x0lc
       tmplc(14) = px0lc
       tmplc(15) = y0lc
       tmplc(16) = py0lc
       tmplc(17) = z0lc
       tmplc(18) = pz0lc
       tmplc(19) = xrlc
       tmplc(20) = pxrlc
       tmplc(21) = yrlc
       tmplc(22) = pyrlc
       tmplc(23) = zrlc
       tmplc(24) = pzrlc
    endif

    call MPI_ALLREDUCE(tmplc,tmpgl,24,MPI_DOUBLE_PRECISION,&
         MPI_SUM,mpicommwd,ierr)

    center1(1) = tmpgl(1)/nptot
    center1(2) = tmpgl(3)/nptot
    center2(1) = tmpgl(13)/nptot
    center2(2) = tmpgl(15)/nptot
    shift = center2 - center1

    if(myidy.lt.npyhalf) then
       center(1) = tmpgl(1)/nptot
       center(2) = tmpgl(2)/nptot
       center(3) = tmpgl(3)/nptot
       center(4) = tmpgl(4)/nptot
       center(5) = tmpgl(5)/nptot
       center(6) = tmpgl(6)/nptot
       sigma(1) = sqrt(tmpgl(7)/nptot - (tmpgl(1)/nptot)**2)
       sigma(2) = sqrt(tmpgl(8)/nptot - (tmpgl(2)/nptot)**2)
       sigma(3) = sqrt(tmpgl(9)/nptot - (tmpgl(3)/nptot)**2)
       sigma(4) = sqrt(tmpgl(10)/nptot - (tmpgl(4)/nptot)**2)
       sigma(5) = sqrt(tmpgl(11)/nptot - (tmpgl(5)/nptot)**2)
       sigma(6) = sqrt(tmpgl(12)/nptot - (tmpgl(6)/nptot)**2)
    else
       center(1) = tmpgl(13)/nptot
       center(2) = tmpgl(14)/nptot
       center(3) = tmpgl(15)/nptot
       center(4) = tmpgl(16)/nptot
       center(5) = tmpgl(17)/nptot
       center(6) = tmpgl(18)/nptot
       sigma(1) = sqrt(tmpgl(19)/nptot - (tmpgl(13)/nptot)**2)
       sigma(2) = sqrt(tmpgl(20)/nptot - (tmpgl(14)/nptot)**2)
       sigma(3) = sqrt(tmpgl(21)/nptot - (tmpgl(15)/nptot)**2)
       sigma(4) = sqrt(tmpgl(22)/nptot - (tmpgl(16)/nptot)**2)
       sigma(5) = sqrt(tmpgl(23)/nptot - (tmpgl(17)/nptot)**2)
       sigma(6) = sqrt(tmpgl(24)/nptot - (tmpgl(18)/nptot)**2)
    endif

  end subroutine shiftsigma2_Utility


  ! Calculate the average centroid and sigma of both beams.
  subroutine findmoments12_Utility(Pts1,innp,nptot,center,sigma)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Pts1
    integer, intent(in) :: innp,nptot
    double precision, dimension(3), intent(out) :: center
    double precision, dimension(6), intent(out) :: sigma
    integer :: i,ierr
    double precision, dimension(12) :: tmplc,tmp

    tmplc = 0.0
    tmp = 0.0 

    do i = 1, innp
       tmplc(1) = tmplc(1) + Pts1(1,i)
       tmplc(2) = tmplc(2) + Pts1(2,i)
       tmplc(3) = tmplc(3) + Pts1(3,i)
       tmplc(4) = tmplc(4) + Pts1(4,i)
       tmplc(5) = tmplc(5) + Pts1(5,i)
       tmplc(6) = tmplc(6) + Pts1(6,i)
       tmplc(7) = tmplc(7) + Pts1(1,i)*Pts1(1,i)
       tmplc(8) = tmplc(8) + Pts1(2,i)*Pts1(2,i)
       tmplc(9) = tmplc(9) + Pts1(3,i)*Pts1(3,i)
       tmplc(10) = tmplc(10) + Pts1(4,i)*Pts1(4,i)
       tmplc(11) = tmplc(11) + Pts1(5,i)*Pts1(5,i)
       tmplc(12) = tmplc(12) + Pts1(6,i)*Pts1(6,i)
    enddo

    tmplc = tmplc/nptot

    !try to obtain the centorid of bunch from each group of processors
    call MPI_ALLREDUCE(tmplc,tmp,12,MPI_DOUBLE_PRECISION,&
         MPI_SUM,mpicommwd,ierr)

    center(1) = tmp(1)
    center(2) = tmp(3)
    center(3) = tmp(5)
    sigma(1) = sqrt(tmp(7)-tmp(1)*tmp(1))
    sigma(2) = sqrt(tmp(8)-tmp(2)*tmp(2))
    sigma(3) = sqrt(tmp(9)-tmp(3)*tmp(3))
    sigma(4) = sqrt(tmp(10)-tmp(4)*tmp(4))
    sigma(5) = sqrt(tmp(11)-tmp(5)*tmp(5))
    sigma(6) = sqrt(tmp(12)-tmp(6)*tmp(6))

  end subroutine findmoments12_Utility



  ! calculate the transverse centroid particles
  subroutine findmomT_Utility(Pts1,maxpt,innp,nptot,center,sigma)
    implicit none
    include 'mpif.h'
    double precision, dimension(6,maxpt) :: Pts1
    integer, intent(in) :: maxpt,innp,nptot
    double precision, dimension(2), intent(out) :: center
    double precision, dimension(2), intent(out) :: sigma
    integer :: i,ierr
    double precision, dimension(4) :: tmplc,tmp

    tmplc = 0.0
    tmp = 0.0 

    do i = 1, innp
       tmplc(1) = tmplc(1) + Pts1(1,i)
       tmplc(2) = tmplc(2) + Pts1(3,i)
       tmplc(3) = tmplc(3) + Pts1(1,i)*Pts1(1,i)
       tmplc(4) = tmplc(4) + Pts1(3,i)*Pts1(3,i)
    enddo

    tmplc = tmplc/nptot

    !try to obtain the centorid of bunch from each group of processors
    call MPI_ALLREDUCE(tmplc,tmp,4,MPI_DOUBLE_PRECISION,&
         MPI_SUM,mpicommwd,ierr)

    center(1) = tmp(1)
    center(2) = tmp(2)
    sigma(1) = sqrt(tmp(3)-tmp(1)*tmp(1))
    sigma(2) = sqrt(tmp(4)-tmp(2)*tmp(2))

  end subroutine findmomT_Utility



  ! calculate the transverse centroids and rms sizes of "opposite beam"
  subroutine findmomO_Utility(Pts1,maxpt,innp,nptot,npyhalf,myidy,center,sigma)
    implicit none
    include 'mpif.h'
    double precision, dimension(6,innp) :: Pts1
    integer, intent(in) :: maxpt,innp,nptot,npyhalf,myidy
    double precision, dimension(2), intent(out) :: center
    double precision, dimension(2), intent(out) :: sigma
    integer :: i,ierr
    double precision, dimension(8) :: tmplc2,tmp2
    double precision, dimension(4) :: tmplc

    tmplc = 0.0d0
    tmp2 = 0.0d0

    do i = 1, innp
       tmplc(1) = tmplc(1) + Pts1(1,i)
       tmplc(2) = tmplc(2) + Pts1(3,i)
       tmplc(3) = tmplc(3) + Pts1(1,i)*Pts1(1,i)
       tmplc(4) = tmplc(4) + Pts1(3,i)*Pts1(3,i)
    enddo

    tmplc = tmplc/nptot

    tmplc2 = 0.0d0
    if(myidy.lt.npyhalf) then
       tmplc2(1:4) = tmplc
    else
       tmplc2(5:8) = tmplc
    endif

    !try to obtain the centorid of bunch from each group of processors
    call MPI_ALLREDUCE(tmplc2,tmp2,8,MPI_DOUBLE_PRECISION,&
         MPI_SUM,mpicommwd,ierr)

    if(myidy.lt.npyhalf) then
       center(1) = tmp2(5)
       center(2) = tmp2(6)
       sigma(1) = sqrt(tmp2(7)-tmp2(5)*tmp2(5))
       sigma(2) = sqrt(tmp2(8)-tmp2(6)*tmp2(6))
    else
       center(1) = tmp2(1)
       center(2) = tmp2(2)
       sigma(1) = sqrt(tmp2(3)-tmp2(1)*tmp2(1))
       sigma(2) = sqrt(tmp2(4)-tmp2(2)*tmp2(2))
    endif

  end subroutine findmomO_Utility



  subroutine Ptlost1_Utility(Ptcls,xmin,xmax,ymin,ymax,nptlc,nptot,npyhalf,&
       myidy,shift21)
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, intent(in) :: xmin,xmax,ymin,ymax
    double precision, dimension(2) :: shift21
    integer, intent(inout) :: nptlc,nptot
    integer, intent(in) :: npyhalf,myidy
    integer :: i,i0,ilost
    integer, dimension(2) :: tmplc, tmp
    double precision :: x,y

    ilost = 0
    do i0 = 1, nptlc
       i = i0 - ilost
       !go to the coordinate used by Hirata.
       if(myidy.lt.npyhalf) then
          x = Ptcls(1,i0) 
          y = Ptcls(3,i0)
       else
          x = Ptcls(1,i0) - shift21(1) 
          y = Ptcls(3,i0) - shift21(2)
       endif
       if(x.le.xmin .or. x.ge.xmax .or. y.le.ymin .or. y.ge.ymax) then
          ilost = ilost + 1
       endif
    enddo
    nptlc = nptlc - ilost
    if(ilost.gt.0) then
       if(myidy.lt.npyhalf) then
          tmplc(1) = nptlc
          tmplc(2) = 0
       else
          tmplc(1) = 0
          tmplc(2) = nptlc
       endif
       call MPI_ALLREDUCE(tmplc,tmp,2,MPI_INTEGER,MPI_SUM,&
            mpicommwd,ierr)
       if(myidy.lt.npyhalf) then
          nptot = tmp(1)
       else
          nptot = tmp(2)
       endif
    endif

  end subroutine Ptlost1_Utility



  subroutine Ptlost2_Utility(Ptcls,xmin,xmax,ymin,ymax,nptlc,nptot,shift21,id)
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, dimension(2) :: shift21
    double precision, intent(in) :: xmin,xmax,ymin,ymax
    integer, intent(inout) :: nptlc,nptot,id
    integer :: i,i0,ilost,sign
    double precision :: x,y

    if(id.eq.1) then
       sign = 0
    else if(id.eq.2) then
       sign = 1
    endif

    ilost = 0
    do i0 = 1, nptlc
       i = i0 - ilost
       !go to the coordinate used by Hirata.
       x = Ptcls(1,i0) - sign*shift21(1)
       y = Ptcls(3,i0) - sign*shift21(2)
       if(x.le.xmin .or. x.ge.xmax .or. y.le.ymin .or. y.ge.ymax) then
          ilost = ilost + 1
       endif
    enddo
    nptlc = nptlc - ilost
    call MPI_ALLREDUCE(nptlc,nptot,1,MPI_INTEGER,MPI_SUM,&
         mpicommwd,ierr)

  end subroutine Ptlost2_Utility



  !//Check the particle loss from transverse aperture size. 
  !//Here, the xmin, xmax are defined as the physical aperture size
  !//with respect to the axis (may not be the close orbit) of the machine.
  !//When the particle gets lost, it is no longer tracked.
  subroutine ptlost2G_Utility(Ptcls,xmin,xmax,ymin,ymax,nptlc,nptot,npyhalf,myidy)
    implicit none
    include 'mpif.h'
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, intent(in) :: xmin,xmax,ymin,ymax
    integer, intent(inout) :: nptlc,nptot
    integer, intent(in) :: npyhalf,myidy
    integer :: i,i0,ilost,ierr
    integer, dimension(2) :: tmplc,tmp
    double precision :: x,y

    ilost = 0
    do i0 = 1, nptlc
       i = i0 - ilost
       x = Ptcls(1,i0) 
       y = Ptcls(3,i0)
       if(x.le.xmin .or. x.ge.xmax .or. y.le.ymin .or. y.ge.ymax) then
          ilost = ilost + 1
       endif
       Ptcls(1,i) = Ptcls(1,i0)
       Ptcls(2,i) = Ptcls(2,i0)
       Ptcls(3,i) = Ptcls(3,i0)
       Ptcls(4,i) = Ptcls(4,i0)
       Ptcls(5,i) = Ptcls(5,i0)
       Ptcls(6,i) = Ptcls(6,i0)
    enddo

    nptlc = nptlc - ilost
    if(myidy.lt.npyhalf) then
       tmplc(1) = nptlc
       tmplc(2) = 0
    else
       tmplc(1) = 0
       tmplc(2) = nptlc
    endif
    call MPI_ALLREDUCE(tmplc,tmp,2,MPI_INTEGER,MPI_SUM,mpicommwd,ierr)
    if(myidy.lt.npyhalf) then
       nptot = tmp(1)
    else
       nptot = tmp(2)
    endif

  end subroutine ptlost2G_Utility



  !//Check the particle loss from transverse aperture size. 
  !//Here, the xmin, xmax are defined as the physical aperture size
  !//with repect to the axis (may not be the close orbit) of the machine.
  !//When the particle gets lost, it is no longer been tracked.
  subroutine ptlost1G_Utility(Ptcls,xmin,xmax,ymin,ymax,nptlc,nptot)
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, intent(in) :: xmin,xmax,ymin,ymax
    integer, intent(inout) :: nptlc,nptot
    integer :: i,i0,ilost
    double precision :: x,y

    ilost = 0
    do i0 = 1, nptlc
       i = i0 - ilost
       x = Ptcls(1,i0) 
       y = Ptcls(3,i0)
       if(x.le.xmin .or. x.ge.xmax .or. y.le.ymin .or. y.ge.ymax) then
          ilost = ilost + 1
       endif
       Ptcls(1,i) = Ptcls(1,i0)
       Ptcls(2,i) = Ptcls(2,i0)
       Ptcls(3,i) = Ptcls(3,i0)
       Ptcls(4,i) = Ptcls(4,i0)
       Ptcls(5,i) = Ptcls(5,i0)
       Ptcls(6,i) = Ptcls(6,i0)
    enddo

    nptlc = nptlc - ilost
    if(ilost.gt.0) then
       call MPI_ALLREDUCE(nptlc,nptot,1,MPI_INTEGER,MPI_SUM,&
            mpicommwd,ierr)
    endif

  end subroutine ptlost1G_Utility



  ! generate "nsize" uniformly distributed random number betwee 0 and 1.
  ! from Numerical Recipes.
  ! initial random seed has to be a negative integer number.
  subroutine ran2(idum,randarray,nsize)
    integer, intent(in) :: nsize
    integer, intent(inout) :: idum
    double precision, dimension(nsize), intent(out) :: randarray
    INTEGER IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    REAL AM,EPS,RNMX
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
         IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,   &
         NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    INTEGER idum2,j,k,iv(NTAB),iy, ii
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/

    do ii = 1, nsize
       if (idum.le.0) then
          idum=max(-idum,1)
          idum2=idum
          do 11 j=NTAB+8,1,-1
             k=idum/IQ1
             idum=IA1*(idum-k*IQ1)-k*IR1
             if (idum.lt.0) idum=idum+IM1
             if (j.le.NTAB) iv(j)=idum
11           continue
          iy=iv(1)
       endif
       k=idum/IQ1
       idum=IA1*(idum-k*IQ1)-k*IR1
       if (idum.lt.0) idum=idum+IM1
       k=idum2/IQ2
       idum2=IA2*(idum2-k*IQ2)-k*IR2
       if (idum2.lt.0) idum2=idum2+IM2
       j=1+iy/NDIV
       iy=iv(j)-idum2
       iv(j)=idum
       if(iy.lt.1)iy=iy+IMM1
       randarray(ii)=min(AM*iy,RNMX)
    enddo

  END subroutine ran2



  !matrix inversion.
  subroutine invmt(amt,np,bmt)
    integer :: np
    double precision, dimension(np,np) :: amt, bmt
    integer :: i, j
    integer, dimension(np) :: indx
    double precision :: d
    double precision, dimension(np) :: tmp

    do i = 1, np
       do j = 1, np
          bmt(i,j) = 0.0
       enddo
       bmt(i,i) = 1.0
    enddo

    call ludcmp(amt,np,np,indx,d)

    do j = 1, np
       do i = 1, np
          tmp(i) = bmt(i,j)
       enddo
       call lubksb(amt,np,np,indx,tmp)
       do i = 1, np
          bmt(i,j) = tmp(i)
       enddo
    enddo

  end subroutine invmt



  !LU decomposition.
  SUBROUTINE ludcmp(a,n,np,indx,d)
    INTEGER n,np,indx(n),NMAX
    double precision d,a(np,np),TINY
    PARAMETER (NMAX=500,TINY=1.0d-20)
    INTEGER i,imax,j,k
    double precision aamax,dum,sum,vv(NMAX)

    d=1.0d0
    do 12 i=1,n
       aamax=0.0d0
       do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11        continue
       if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
          vv(i)=1./aamax
12        continue
    do 19 j=1,n
       do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
             sum=sum-a(i,k)*a(k,j)
13           continue
          a(i,j)=sum
14           continue
       aamax=0.d0
       do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
             sum=sum-a(i,k)*a(k,j)
15           continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
             imax=i
             aamax=dum
          endif
16        continue
       if (j.ne.imax)then
          do 17 k=1,n
             dum=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=dum
17           continue
          d=-d
          vv(imax)=vv(j)
       endif
       indx(j)=imax
       if(a(j,j).eq.0.0d0)a(j,j)=TINY
       if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
             a(i,j)=a(i,j)*dum
18           continue
       endif
19     continue
   !  return
  END subroutine ludcmp



  SUBROUTINE lubksb(a,n,np,indx,b)
    INTEGER n,np,indx(n)
    double precision a(np,np),b(n)
    INTEGER i,ii,j,ll
    double precision sum
    ii=0
    do 12 i=1,n
       ll=indx(i)
       sum=b(ll)
       b(ll)=b(i)
       if (ii.ne.0)then
          do 11 j=ii,i-1
             sum=sum-a(i,j)*b(j)
11           continue
       else if (sum.ne.0.) then
          ii=i
       endif
       b(i)=sum
12     continue
    do 14 i=n,1,-1
       sum=b(i)
       do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13        continue
       b(i)=sum/a(i,i)
14     continue
   !  return
  END subroutine lubksb


  
  subroutine matrixmul(mt1,mt2,mt3,l,m,n)
    integer :: l,m,n
    double precision, dimension(l,m) :: mt1
    double precision, dimension(m,n) :: mt2
    double precision, dimension(l,n) :: mt3
    integer :: i,j,k

    do j = 1, n
       do i = 1, l
          mt3(i,j) = 0.0
          do k = 1, m
             mt3(i,j) = mt3(i,j) + mt1(i,k)*mt2(k,j)
          enddo
       enddo
    enddo

  end subroutine matrixmul



  !// Generate random list of indices of a random subset of particles
  subroutine randomIndexList(Nplocal,nSamples,sampleInd)
    integer, intent(in) :: Nplocal, nSamples
    integer, intent(out), dimension(nSamples) :: sampleInd
    integer :: ii
    integer, dimension(Nplocal) :: indList
    double precision, allocatable, dimension(:) :: randNum

    if(nSamples>Nplocal) print*,"Warning: Sample larger than number of available particles in Utility::randomIndexList."
    do ii=1,Nplocal
       indList(ii) = ii
    end do
    if(nSamples==Nplocal)then
       sampleInd = indList
    else
       allocate(randNum(nSamples))
       call random_number(randNum)
       do ipt=1, nSamples
          ii = 1+int(randNum(ipt)*(Nplocal+1-ipt))  !// pick random index in remaining set; the range is decremented every step
          sampleInd(ipt) = indList(ii)  !// copy picked index to sample collection
          indList(ii) = indList(Nplocal+1-ipt)  !// copy last index to position of picked one
       enddo
       deallocate(randNum)
    endif
  end subroutine randomIndexList

end module Utilityclass
