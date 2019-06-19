!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Beambeam class:
! Version: 1.0
! Author: Ji Qiang, LBNL
! Description: 
! Comments:
!----------------------------------------------------------------
module Beambeamclass
  use CompDomclass
  use Pgrid2dclass
  use Timerclass
  use Outputclass
  use Communclass
  use DepoScatclass
  use FieldSolverclass
  use Utilityclass
  !include "mpif.h"
  !type beambeam
  !  double precision :: tunex,tuney,tunez
  !  double precision :: alphax,betax,emitx,alphay,betay,emity 
  !  double precision, dimension(2,2) :: matx,maty,matz
  !end type beambeam


contains

  !//find the strong-strong beam-beam kick in a 3D multi-slice model.
  !//Here, we will use the integrated Green phi function to handle
  !//large aspect ratio. The slice of field is slightly different the
  !//particle slice for interpolation.
  subroutine bbmSlcnew_Beambeam(nprocrow,npyhalf,innx,inny,&
       nslice1,nslice2,nslice,nsliceO,jendmax,nptlc,myidx,myidy,Pts,&
       grid2d,Ageom,zmin,zmax,nptot,coef,nytot,commrow,commcol,nz,&
       maxPtclSlice,shift21,curin,flglum,lum3d)
    implicit none
    integer, intent(in) :: nprocrow,npyhalf,innx,nslice1,&
         nslice2,myidx,myidy,nptot,nytot,commrow,commcol,nz,&
         maxPtclSlice,nslice,nsliceO,jendmax,flglum
    integer, intent(inout) :: inny,nptlc
    double precision, pointer, dimension(:,:) :: Pts
    double precision, intent(out) :: lum3d
    type (Pgrid2d), intent(in) :: grid2d
    type (CompDom), intent(inout) :: Ageom
    double precision, intent(in) :: zmin,coef,zmax,curin
    double precision, dimension(2), intent(in) :: shift21
    double precision, dimension(2) :: shift12,shift
    double precision, dimension(nsliceO) :: zslice2,hzslice2
    double precision, dimension(nslice) :: weight,zslice,hzslice
    double precision, dimension(6) :: range
    double precision, dimension(4) :: range1
    double precision, dimension(3) :: msize
    double precision, dimension(innx,nytot) :: rho,phiout,phiout2
    double precision, dimension(innx,inny) :: rholc
    double precision, dimension(innx,inny,jendmax) :: rhotot,rhotot2
    double precision, dimension(6,maxPtclSlice) :: ray
    double precision, dimension(7,maxPtclSlice,nslice) :: beamslice
    double precision :: xmin,ymin,hx,hy
    double complex, allocatable, dimension (:,:) :: grn
    integer, dimension(nslice) :: countlc,count
    integer, dimension(2,0:nprocrow-1,0:npyhalf-1) :: table
    integer :: i,j,nptlcslice,sign,ii,jj,ipt,k,ierr,&
         i1,j1,jend,istart,jstart,myjend,itmp1,itmp2,rowend
    integer, dimension(0:npyhalf-1) :: xpystable,pytable
    integer :: npbc,nsxy1,nsxy2,nxpylc2,nxtot,itb,nxlum,nylum
    double precision :: xtmptmp,lumtmp
    double precision :: tmpx1,tmpx2,tmpy1,tmpy2,distance1,distance2
    integer :: ibf,sign2
    real*8 :: zbk,hzi1,tmpx3,tmpy3,distance3

    !//shift12 - beam 1 with respect to beam 2
    !//shift21 - beam 2 with respect to beam 1
    shift12 = -shift21
    nxtot = innx
    xtmptmp = 1.0
    lum3d = 0.0
    nxlum = 128
    nylum = 128

    !//set up the parameters for the lower group and upper group processors.
    if(myidy.lt.npyhalf) then
       shift = shift21
    else
       shift = shift12
    endif
    if(myidy.lt.npyhalf) then
       npbc = 0
    else
       npbc = npyhalf
    endif

    ! +1 is from the real to complex fft.
    nsxy1 = (nxtot+1)/npyhalf
    nsxy2 = (nxtot+1) - npyhalf*nsxy1
    do i = 0, npyhalf-1
       if(i.le.(nsxy2-1)) then
          xpystable(i) = nsxy1+1
       else
          xpystable(i) = nsxy1
       endif
    enddo

    nxpylc2 = xpystable(myidy-npbc)

    allocate(grn(2*nytot,nxpylc2))

    !//here, we have used equal area weight
    !//copy particles from a bunch into slices. Each slice has
    !//roughly same number of particles.
    call slicernew(Pts,nptlc,nslice,nsliceO,maxPtclSlice,nz,zmin,&
         zmax,beamslice,weight,zslice,zslice2,hzslice,hzslice2,nptot,count,countlc,myidy,&
         npyhalf,commcol,commrow)
    call MPI_BARRIER(mpicommwd,ierr)

    !//loop through the total overlaping steps.
    do i = 1, nslice1+nslice2-1
       !//find the number of colliding slices 
       if(i.le.nslice1) then
          jend = min(i,nslice2) !//jend is the number of the overlaped slices.
          istart = i+1  !//istart is the starting slice + 1 in beam1
          jstart = 0    !//jstart is is the starting slice - 1 in beam2
       else
          jend = min((nslice2-i+nslice1),nslice1)
          istart = nslice1+1
          jstart = i-nslice1
       endif

       !//find the range of the colliding slices
       !find the range xmin,xmax,ymin,ymax for all the colliding slices,jend.
       if(myidy.lt.npyhalf) then
          !distance2 = (beamslice(5,ipt,i1)-zslice2(j1))/2
          !distance1 = (zslice(istart-1)-zslice2(jstart+1))/2
          distance1 = (beamslice(5,1,istart-1)-zslice2(jstart+1))/2
          tmpx1 = beamslice(1,1,istart-1) + distance1*beamslice(2,1,istart-1)
          tmpy1 = beamslice(3,1,istart-1) + distance1*beamslice(4,1,istart-1)
          range1(1) = tmpx1
          range1(2) = tmpx1
          range1(3) = tmpy1
          range1(4) = tmpy1
       else
          distance1 = (beamslice(5,1,jstart+1)-zslice2(istart-1))/2
          !distance1 = (zslice(jstart+1)-zslice2(istart-1))/2
          tmpx1 = beamslice(1,1,jstart+1) + distance1*beamslice(2,1,jstart+1)
          tmpy1 = beamslice(3,1,jstart+1) + distance1*beamslice(4,1,jstart+1)
          range1(1) = tmpx1
          range1(2) = tmpx1
          range1(3) = tmpy1
          range1(4) = tmpy1
       endif
       do j = 1, jend !//loop through the number of overlaped slices.
          if(myidy.lt.npyhalf) then
             i1 = istart - j
             j1 = jstart + j
          else
             i1 = jstart + j
             j1 = istart - j
          endif
          !distance1 = (zslice(i1)-zslice2(j1))/2
          distance1 = (zslice(i1)+0.5*hzslice2(j1)-zslice2(j1))/2
          distance3 = (zslice(i1)-0.5*hzslice2(j1)-zslice2(j1))/2
          do ipt = 1, countlc(i1)
             distance2 = (beamslice(5,ipt,i1)-zslice2(j1))/2
             tmpx1 = beamslice(1,ipt,i1) + distance1*beamslice(2,ipt,i1)
             tmpx2 = beamslice(1,ipt,i1) + distance2*beamslice(2,ipt,i1)
             tmpx3 = beamslice(1,ipt,i1) + distance3*beamslice(2,ipt,i1)
             tmpy1 = beamslice(3,ipt,i1) + distance1*beamslice(4,ipt,i1)
             tmpy2 = beamslice(3,ipt,i1) + distance2*beamslice(4,ipt,i1)
             tmpy3 = beamslice(3,ipt,i1) + distance3*beamslice(4,ipt,i1)
             if(range1(1).gt.tmpx1) range1(1) = tmpx1
             if(range1(1).gt.tmpx2) range1(1) = tmpx2
             if(range1(1).gt.tmpx3) range1(1) = tmpx3
             if(range1(2).lt.tmpx1) range1(2) = tmpx1
             if(range1(2).lt.tmpx2) range1(2) = tmpx2
             if(range1(2).lt.tmpx3) range1(2) = tmpx3
             if(range1(3).gt.tmpy1) range1(3) = tmpy1
             if(range1(3).gt.tmpy2) range1(3) = tmpy2
             if(range1(3).gt.tmpy3) range1(3) = tmpy3
             if(range1(4).lt.tmpy1) range1(4) = tmpy1
             if(range1(4).lt.tmpy2) range1(4) = tmpy2
             if(range1(4).lt.tmpy3) range1(4) = tmpy3
          enddo
       enddo

       !//assume the Green function is same for all colliding slices at each step.
       !//assume all colliding slices at each step are contained in the same box.
       !//shift the range to local range, ie. the range with respect
       !//to the center of the beam.
       if(myidy.ge.npyhalf) then
          range1(1) = range1(1) - shift21(1)
          range1(2) = range1(2) - shift21(1)
          range1(3) = range1(3) - shift21(2)
          range1(4) = range1(4) - shift21(2)
       endif
       !
       !//update the computational domain after drift
       call update2d_CompDom(Ageom,range1,grid2d)
       call getmsize_CompDom(Ageom,msize)
       hx = msize(1)
       hy = msize(2)
       call getlctabnm_CompDom(Ageom,table)
       do itb = 0, npyhalf-1
          pytable(itb) = table(2,0,itb)
       enddo
       !//find the Green function for the domain containing "jend"
       !//slice.
       call greenf2d(nxtot,nytot,nxpylc2,hx,hy,myidy,npyhalf,commcol,&
            xpystable,grn,shift)
       call getrange_CompDom(Ageom,range)
       if(myidy.lt.npyhalf) then
          xmin = range(1)
          ymin = range(3)
       else
          xmin = range(1) + shift21(1)
          ymin = range(3) + shift21(2)
       endif
       !endif

       rhotot = 0.0d0
       rhotot2 = 0.0d0
       do ibf = 1, 2 

          !--------------------------------------------------------------------------------
          !//drift the collision slice particles to the collision point.
          sign = 1
          if(ibf.eq.1) then
             sign2 = -1
          else
             sign2 = 1
          endif
          if(nslice.gt.1) then
             call driftCPfixnew(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
                  istart,jstart,jend,zslice,zslice2,hzslice2,countlc,maxPtclSlice,&
                  sign2)
          endif

          !//find the transverse charge density of overlaped slices
          do j = 1, jend
             if(myidy.lt.npyhalf) then
                i1 = istart - j 
             else
                i1 = jstart + j
             endif

             nptlcslice = countlc(i1)
             do k = 1, nptlcslice
                ray(:,k) = beamslice(1:6,k,i1)
             enddo

             !            !//calculate the 3d luminosity
             !            if(flglum.eq.1) then
             !              call luminosity2G3d_Output(ray,nptlcslice,nxlum,nylum,&
             !                    curin,weight(i1),count(i1),myidy,npyhalf,lumtmp)
             !              lum3d = lum3d + lumtmp
             !            endif

             !//deposition for each slice.
             call deposit2d(maxPtclSlice,nptlcslice,innx,nytot,hx,hy,xmin,&
                  ymin,ray,rho)

             !//collect the contribution from the other processors in the same column.
             !//This function and the function guardsum2drow can be replaced by
             !//a Allreduce in the whole domain with double sized container.
             call guardsum2d(rho,rholc,innx,inny,nytot,npyhalf,myidy,commcol)

             if(ibf.eq.1) then
                do jj = 1, inny
                   do ii = 1, innx
                      rhotot(ii,jj,j) = rholc(ii,jj)*weight(i1)/count(i1)
                   enddo
                enddo
             else
                do jj = 1, inny
                   do ii = 1, innx
                      rhotot2(ii,jj,j) = rholc(ii,jj)*weight(i1)/count(i1)
                   enddo
                enddo
             endif
          enddo

          !//collect the density onto the colliding slices
          myjend = (jend-1)/nprocrow + 1
          if(ibf.eq.1) then
             if(nslice.gt.1) then
                call guardsum2drow(rhotot,innx,inny,myjend,jendmax,nprocrow,&
                     myidx,commrow)
             endif
          else
             if(nslice.gt.1) then
                call guardsum2drow(rhotot2,innx,inny,myjend,jendmax,nprocrow,&
                     myidx,commrow)
             endif
          endif

          !//distribute "jend" slice along "nprocrow" row processors.
          itmp1 = jend/nprocrow
          itmp2 = jend - nprocrow*itmp1
          if(myidx.lt.itmp2) then
             rowend = itmp1 + 1
          else
             rowend = itmp1
          endif

          !//solve the Poisson's equation in parallel for processors along row.
          !//transpose will be used along column to solve the Poisson equation. 
          do j = 1, rowend

             if(ibf.eq.1) then
                do jj = 1, inny
                   do ii = 1, innx
                      rholc(ii,jj) =  rhotot(ii,jj,j)
                   enddo
                enddo
             else
                do jj = 1, inny
                   do ii = 1, innx
                      rholc(ii,jj) =  rhotot2(ii,jj,j)
                   enddo
                enddo
             endif

             !//solve the potential for each slice of distributed "jend" slice.
             !call fieldsolver2d(rholc,innx,inny,nprocrow,npyhalf,&
             !nytot,nxpylc2,xpystable,pytable,hx,hy,myidy,commcol,shift,grn)
             call fieldsolver2d(rholc,innx,inny,nprocrow,npyhalf,&
                  nytot,nxpylc2,xpystable,pytable,hx,hy,myidy,commcol,grn)


             if(ibf.eq.1) then !store the potential at back of the slice
                do jj = 1, inny
                   do ii = 1, innx
                      rhotot(ii,jj,j) =  rholc(ii,jj)
                   enddo
                enddo
             else !store the potential at front of the slice
                do jj = 1, inny
                   do ii = 1, innx
                      rhotot2(ii,jj,j) =  rholc(ii,jj)
                   enddo
                enddo
             endif

          enddo

          !drift back the fixed distance
          sign = -1
          if(ibf.eq.1) then
             sign2 = -1
          else
             sign2 = 1
          endif
          if(nslice.gt.1) then
             call driftCPfixnew(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
                  istart,jstart,jend,zslice,zslice2,hzslice2,countlc,maxPtclSlice,sign2)
          endif
          !--------------------------------------------------------------------------------

       enddo

       !drift to the field interpolation point
       sign = 1
       if(nslice.gt.1) then
          call driftCP(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
               istart,jstart,jend,zslice,zslice2,countlc,maxPtclSlice)
       endif

       !//scatter potentials along row so that each processor contain 
       !//"jend" slices potential. This function can be replaced by
       !//Allreduce command along row.
       if(nslice.gt.1) then
          call guardexch2drow(rhotot,innx,inny,jendmax,myjend,nprocrow,&
               myidx,commrow)
          call guardexch2drow(rhotot2,innx,inny,jendmax,myjend,nprocrow,&
               myidx,commrow)
       endif

       !interpolate the beam-beam field onto particles on each slice.
       do j = 1, jend 
          if(myidy.lt.npyhalf) then
             i1 = istart - j 
             j1 = jstart + j
          else
             i1 = jstart + j
             j1 = istart -j
          endif

          nptlcslice = countlc(i1)
          do k = 1, nptlcslice
             ray(:,k) = beamslice(1:6,k,i1)
          enddo

          do jj = 1, inny
             do ii = 1, innx
                rholc(ii,jj) =  rhotot(ii,jj,j)
             enddo
          enddo

          !//collect all the potential along column to a single proc. This
          !//function can be replaced by Allgather along column with doubled size.
          call guardexch2d(rholc,phiout,innx,inny,nytot,npyhalf,myidy,commcol)

          do jj = 1, inny
             do ii = 1, innx
                rholc(ii,jj) =  rhotot2(ii,jj,j)
             enddo
          enddo
          call guardexch2d(rholc,phiout2,innx,inny,nytot,npyhalf,myidy,commcol)

          !//interpolate field of each slice to slice particles
          zbk = zslice(i1)-0.5*hzslice(i1)
          hzi1 = hzslice(i1)
          call scatter2dnew(maxPtclSlice,nptlcslice,innx,&
               hx,hy,xmin,ymin,ray,phiout,phiout2,xtmptmp,coef,&
               nytot,zbk,hzi1)

          do k = 1, nptlcslice
             beamslice(1:6,k,i1) = ray(:,k)
          enddo

          !//calculate the 3d luminosity
          if(flglum.eq.1) then
             call luminosity2G3d_Output(ray,nptlcslice,nxlum,nylum,&
                  curin,weight(i1),count(i1),myidy,npyhalf,lumtmp)
             lum3d = lum3d + lumtmp
          endif
       enddo

       !//drift all particles back to their original longitudinal location
       if(nslice.gt.1) then
          sign = -1
          call driftCP(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
               istart,jstart,jend,zslice,zslice2,countlc,maxPtclSlice)
       endif
    enddo

    deallocate(grn)

    call MPI_BARRIER(mpicommwd,ierr)
    !//copy each slice particle to the beam bunch 
    ipt = 0
    do i = 1, nslice
       do j = 1, countlc(i)
          !ipt = ipt + 1
          ipt = int(beamslice(7,j,i)+0.01)
          Pts(:,ipt) = beamslice(1:6,j,i)
       enddo
    enddo

  end subroutine bbmSlcnew_Beambeam

  subroutine bbmSlcfixnew_Beambeam(nprocrow,npyhalf,innx,inny,&
       nslice1,nslice2,nslice,nsliceO,jendmax,nptlc,myidx,myidy,Pts,&
       grid2d,Ageom,zmin,zmax,nptot,coef,nytot,commrow,commcol,nz,&
       maxPtclSlice,nxpylc2,grn,xpystable,pytable,hx,hy,xmin,ymin,&
       curin,flglum,lum3d)
    implicit none
    integer, intent(in) :: nprocrow,npyhalf,innx,nslice1,&
         nslice2,myidx,myidy,nptot,nytot,commrow,commcol,nz,&
         maxPtclSlice,nslice,nsliceO,jendmax,nxpylc2,flglum
    integer, intent(inout) :: inny,nptlc
    double precision, pointer, dimension(:,:) :: Pts
    double precision, intent(out) :: lum3d
    type (Pgrid2d), intent(in) :: grid2d
    type (CompDom), intent(inout) :: Ageom
    double precision, intent(in) :: zmin,coef,zmax,curin
    !        double precision, dimension(2), intent(in) :: shift21
    !        double precision, dimension(2) :: shift12,shift
    double precision, dimension(nsliceO) :: zslice2,hzslice2
    double precision, dimension(nslice) :: weight,zslice,hzslice
    !        double precision, dimension(6) :: range
    !        double precision, dimension(4) :: range1
    !        double precision, dimension(3) :: msize
    double precision, dimension(innx,nytot) :: rho,phiout,phiout2
    double precision, dimension(innx,inny) :: rholc
    double precision, dimension(innx,inny,jendmax) :: rhotot,rhotot2
    double precision, dimension(6,maxPtclSlice) :: ray
    double precision, dimension(6,maxPtclSlice,nslice) :: beamslice
    double precision :: xmin,ymin,hx,hy
    double complex, dimension (2*nytot,nxpylc2) :: grn
    integer, dimension(nslice) :: countlc,count
!    integer, dimension(2,0:nprocrow-1,0:npyhalf-1) :: table
    integer :: i,j,nptlcslice,sign,ii,jj,ipt,k,ierr,&
         i1,j1,jend,istart,jstart,myjend,itmp1,itmp2,rowend
    integer, dimension(0:npyhalf-1) :: xpystable,pytable
    integer :: nxtot,nxlum,nylum !,npbc,nsxy1,nsxy2
    double precision :: xtmptmp,lumtmp
!    double precision :: tmpx1,tmpx2,tmpx3,tmpy1,tmpy2,tmpy3
    integer :: ibf,sign2
    real*8 :: zbk,hzi1 !,distance1,distance2distance3

    !//shift12 - beam 1 with respect to beam 2
    !//shift21 - beam 2 with respect to beam 1
    !        shift12 = -shift21
    nxtot = innx
    xtmptmp = 1.0
    lum3d = 0.0d0
    nxlum = 128
    nylum = 128

    !        !//set up the parameters for the lower group and upper group processors.
    !        if(myidy.lt.npyhalf) then
    !           shift = shift21
    !        else
    !           shift = shift12
    !        endif
    !        if(myidy.lt.npyhalf) then
    !          npbc = 0
    !        else
    !          npbc = npyhalf
    !        endif
    !
    !        ! +1 is from the real to complex fft.
    !        nsxy1 = (nxtot+1)/npyhalf
    !        nsxy2 = (nxtot+1) - npyhalf*nsxy1
    !        do i = 0, npyhalf-1
    !          if(i.le.(nsxy2-1)) then
    !            xpystable(i) = nsxy1+1
    !          else
    !            xpystable(i) = nsxy1
    !          endif
    !        enddo
    !
    !        nxpylc2 = xpystable(myidy-npbc)
    !
    !        allocate(grn(2*nytot,nxpylc2))

    !//here, we have used equal area weight
    !//copy particles from a bunch into slices. Each slice has
    !//roughly same number of particles.
    call slicernew(Pts,nptlc,nslice,nsliceO,maxPtclSlice,nz,zmin,&
         zmax,beamslice,weight,zslice,zslice2,hzslice,hzslice2,nptot,count,countlc,myidy,&
         npyhalf,commcol,commrow)

    call MPI_BARRIER(mpicommwd,ierr)

    !//loop through the total overlaping steps.
    do i = 1, nslice1+nslice2-1
       !//find the number of colliding slices 
       if(i.le.nslice1) then
          jend = min(i,nslice2) !//jend is the number of the overlaped slices.
          istart = i+1  !//istart is the starting slice + 1 in beam1
          jstart = 0    !//jstart is is the starting slice - 1 in beam2
       else
          jend = min((nslice2-i+nslice1),nslice1)
          istart = nslice1+1
          jstart = i-nslice1
       endif

       !//find the range of the colliding slices
       !          !find the range xmin,xmax,ymin,ymax for all the colliding slices,jend.
       !          if(myidy.lt.npyhalf) then
       !              distance2 = (beamslice(5,ipt,i1)-zslice2(j1))/2
       !              !distance1 = (zslice(istart-1)-zslice2(jstart+1))/2
       !              tmpx1 = beamslice(1,1,istart-1) + distance1*beamslice(2,1,istart-1)
       !              tmpy1 = beamslice(3,1,istart-1) + distance1*beamslice(4,1,istart-1)
       !              range1(1) = tmpx1
       !              range1(2) = tmpx1
       !              range1(3) = tmpy1
       !              range1(4) = tmpy1
       !          else
       !              distance1 = (beamslice(5,1,jstart+1)-zslice2(istart-1))/2
       !              !distance1 = (zslice(jstart+1)-zslice2(istart-1))/2
       !              tmpx1 = beamslice(1,1,jstart+1) + distance1*beamslice(2,1,jstart+1)
       !              tmpy1 = beamslice(3,1,jstart+1) + distance1*beamslice(4,1,jstart+1)
       !              range1(1) = tmpx1
       !              range1(2) = tmpx1
       !              range1(3) = tmpy1
       !              range1(4) = tmpy1
       !          endif
       !          do j = 1, jend !//loop through the number of overlaped slices.
       !            if(myidy.lt.npyhalf) then
       !              i1 = istart - j
       !              j1 = jstart + j
       !            else
       !              i1 = jstart + j
       !              j1 = istart - j
       !            endif
       !            !distance1 = (zslice(i1)-zslice2(j1))/2
       !            distance1 = (zslice(i1)+0.5*hzslice2(j1)-zslice2(j1))/2
       !            distance3 = (zslice(i1)-0.5*hzslice2(j1)-zslice2(j1))/2
       !            do ipt = 1, countlc(i1)
       !              distance2 = (beamslice(5,ipt,i1)-zslice2(j1))/2
       !              tmpx1 = beamslice(1,ipt,i1) + distance1*beamslice(2,ipt,i1)
       !              tmpx2 = beamslice(1,ipt,i1) + distance2*beamslice(2,ipt,i1)
       !              tmpx3 = beamslice(1,ipt,i1) + distance3*beamslice(2,ipt,i1)
       !              tmpy1 = beamslice(3,ipt,i1) + distance1*beamslice(4,ipt,i1)
       !              tmpy2 = beamslice(3,ipt,i1) + distance2*beamslice(4,ipt,i1)
       !              tmpy3 = beamslice(3,ipt,i1) + distance3*beamslice(4,ipt,i1)
       !              if(range1(1).gt.tmpx1) range1(1) = tmpx1
       !              if(range1(1).gt.tmpx2) range1(1) = tmpx2
       !              if(range1(1).gt.tmpx3) range1(1) = tmpx3
       !              if(range1(2).lt.tmpx1) range1(2) = tmpx1
       !              if(range1(2).lt.tmpx2) range1(2) = tmpx2
       !              if(range1(2).lt.tmpx3) range1(2) = tmpx3
       !              if(range1(3).gt.tmpy1) range1(3) = tmpy1
       !              if(range1(3).gt.tmpy2) range1(3) = tmpy2
       !              if(range1(3).gt.tmpy3) range1(3) = tmpy3
       !              if(range1(4).lt.tmpy1) range1(4) = tmpy1
       !              if(range1(4).lt.tmpy2) range1(4) = tmpy2
       !              if(range1(4).lt.tmpy3) range1(4) = tmpy3
       !            enddo
       !          enddo
       !
       !          !//assume the Green function is same for all colliding slices at each step.
       !          !//assume all colliding slices at each step are contained in the same box.
       !            !//shift the range to local range, ie. the range with respect
       !            !//to the center of the beam.
       !            if(myidy.ge.npyhalf) then
       !              range1(1) = range1(1) - shift21(1)
       !              range1(2) = range1(2) - shift21(1)
       !              range1(3) = range1(3) - shift21(2)
       !              range1(4) = range1(4) - shift21(2)
       !            endif
       !
       !            !//update the computational domain after drift
       !            call update2d_CompDom(Ageom,range1,grid2d)
       !            call getmsize_CompDom(Ageom,msize)
       !            hx = msize(1)
       !            hy = msize(2)
       !            call getlctabnm_CompDom(Ageom,table)
       !            do itb = 0, npyhalf-1
       !              pytable(itb) = table(2,0,itb)
       !            enddo
       !            !//find the Green function for the domain containing "jend"
       !            !//slice.
       !            call greenf2d(nxtot,nytot,nxpylc2,hx,hy,myidy,npyhalf,commcol,&
       !                          xpystable,grn,shift)
       !            call getrange_CompDom(Ageom,range)
       !            if(myidy.lt.npyhalf) then
       !              xmin = range(1)
       !              ymin = range(3)
       !            else
       !              xmin = range(1) + shift21(1)
       !              ymin = range(3) + shift21(2)
       !            endif
       !          !endif

       rhotot = 0.0
       rhotot2 = 0.0
       do ibf = 1, 2 

          !----------------------------------------------------------------------------
          !//drift the collision slice particles to the collision point.
          sign = 1
          if(ibf.eq.1) then
             sign2 = -1
          else
             sign2 = 1
          endif
          if(nslice.gt.1) then
             call driftCPfixnew(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
                  istart,jstart,jend,zslice,zslice2,hzslice2,countlc,maxPtclSlice,&
                  sign2)
          endif

          !//find the transverse charge density of overlaped slices
          do j = 1, jend
             if(myidy.lt.npyhalf) then
                i1 = istart - j 
             else
                i1 = jstart + j
             endif

             nptlcslice = countlc(i1)
             do k = 1, nptlcslice
                ray(:,k) = beamslice(:,k,i1)
             enddo

             !            !//calculate the 3d luminosity
             !            if(flglum.eq.1) then
             !              call luminosity2G3d_Output(ray,nptlcslice,nxlum,nylum,&
             !                    curin,weight(i1),count(i1),myidy,npyhalf,lumtmp)
             !              lum3d = lum3d + lumtmp
             !            endif

             !//deposition for each slice.
             call deposit2d(maxPtclSlice,nptlcslice,innx,nytot,hx,hy,xmin,&
                  ymin,ray,rho)

             !//collect the contribution from the other processors in the same column.
             !//This function and the function guardsum2drow can be replaced by
             !//a Allreduce in the whole domain with double sized container.
             call guardsum2d(rho,rholc,innx,inny,nytot,npyhalf,myidy,commcol)

             if(ibf.eq.1) then
                do jj = 1, inny
                   do ii = 1, innx
                      rhotot(ii,jj,j) = rholc(ii,jj)*weight(i1)/count(i1)
                   enddo
                enddo
             else
                do jj = 1, inny
                   do ii = 1, innx
                      rhotot2(ii,jj,j) = rholc(ii,jj)*weight(i1)/count(i1)
                   enddo
                enddo
             endif
          enddo

          !//collect the density onto the colliding slices
          myjend = (jend-1)/nprocrow + 1
          if(ibf.eq.1) then
             if(nslice.gt.1) then
                call guardsum2drow(rhotot,innx,inny,myjend,jendmax,nprocrow,&
                     myidx,commrow)
             endif
          else
             if(nslice.gt.1) then
                call guardsum2drow(rhotot2,innx,inny,myjend,jendmax,nprocrow,&
                     myidx,commrow)
             endif
          endif

          !//distribute "jend" slice along "nprocrow" row processors.
          itmp1 = jend/nprocrow
          itmp2 = jend - nprocrow*itmp1
          if(myidx.lt.itmp2) then
             rowend = itmp1 + 1
          else
             rowend = itmp1
          endif

          !//solve the Poisson's equation in parallel for processors along row.
          !//transpose will be used along column to solve the Poisson equation. 
          do j = 1, rowend

             if(ibf.eq.1) then
                do jj = 1, inny
                   do ii = 1, innx
                      rholc(ii,jj) =  rhotot(ii,jj,j)
                   enddo
                enddo
             else
                do jj = 1, inny
                   do ii = 1, innx
                      rholc(ii,jj) =  rhotot2(ii,jj,j)
                   enddo
                enddo
             endif

             !//solve the potential for each slice of distributed "jend" slice.
             !call fieldsolver2d(rholc,innx,inny,nprocrow,npyhalf,&
             !nytot,nxpylc2,xpystable,pytable,hx,hy,myidy,commcol,shift,grn)
             call fieldsolver2d(rholc,innx,inny,nprocrow,npyhalf,&
                  nytot,nxpylc2,xpystable,pytable,hx,hy,myidy,commcol,grn)


             if(ibf.eq.1) then !store the potential at back of the slice
                do jj = 1, inny
                   do ii = 1, innx
                      rhotot(ii,jj,j) =  rholc(ii,jj)
                   enddo
                enddo
             else !store the potential at front of the slice
                do jj = 1, inny
                   do ii = 1, innx
                      rhotot2(ii,jj,j) =  rholc(ii,jj)
                   enddo
                enddo
             endif

          enddo

          !drift back the fixed distance
          sign = -1
          if(ibf.eq.1) then
             sign2 = -1
          else
             sign2 = 1
          endif
          if(nslice.gt.1) then
             call driftCPfixnew(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
                  istart,jstart,jend,zslice,zslice2,hzslice2,countlc,maxPtclSlice,sign2)
          endif
          !----------------------------------------------------------------------------

       enddo

       !drift to the field interpolation point
       sign = 1
       if(nslice.gt.1) then
          call driftCP(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
               istart,jstart,jend,zslice,zslice2,countlc,maxPtclSlice)
       endif

       !//scatter potentials along row so that each processor contain 
       !//"jend" slices potential. This function can be replaced by
       !//Allreduce command along row.
       if(nslice.gt.1) then
          call guardexch2drow(rhotot,innx,inny,jendmax,myjend,nprocrow,&
               myidx,commrow)
          call guardexch2drow(rhotot2,innx,inny,jendmax,myjend,nprocrow,&
               myidx,commrow)
       endif

       !interpolate the beam-beam field onto particles on each slice.
       do j = 1, jend 
          if(myidy.lt.npyhalf) then
             i1 = istart - j 
             j1 = jstart + j
          else
             i1 = jstart + j
             j1 = istart -j
          endif

          nptlcslice = countlc(i1)
          do k = 1, nptlcslice
             ray(:,k) = beamslice(:,k,i1)
          enddo

          do jj = 1, inny
             do ii = 1, innx
                rholc(ii,jj) =  rhotot(ii,jj,j)
             enddo
          enddo

          !//collect all the potential along column to a single proc. This
          !//function can be replaced by Allgather along column with doubled size.
          call guardexch2d(rholc,phiout,innx,inny,nytot,npyhalf,myidy,commcol)

          do jj = 1, inny
             do ii = 1, innx
                rholc(ii,jj) =  rhotot2(ii,jj,j)
             enddo
          enddo
          call guardexch2d(rholc,phiout2,innx,inny,nytot,npyhalf,myidy,commcol)

          !//interpolate field of each slice to slice particles
          zbk = zslice(i1)-0.5*hzslice(i1)
          hzi1 = hzslice(i1)
          call scatter2dnew(maxPtclSlice,nptlcslice,innx,&
               hx,hy,xmin,ymin,ray,phiout,phiout2,xtmptmp,coef,&
               nytot,zbk,hzi1)

          do k = 1, nptlcslice
             beamslice(:,k,i1) = ray(:,k)
          enddo

          !//calculate the 3d luminosity
          if(flglum.eq.1) then
             call luminosity2G3d_Output(ray,nptlcslice,nxlum,nylum,&
                  curin,weight(i1),count(i1),myidy,npyhalf,lumtmp)
             lum3d = lum3d + lumtmp
          endif

       enddo

       !//drift all particles back to their original longitudinal location
       if(nslice.gt.1) then
          sign = -1
          call driftCP(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
               istart,jstart,jend,zslice,zslice2,countlc,maxPtclSlice)
       endif
    enddo

    !        deallocate(grn)

    call MPI_BARRIER(mpicommwd,ierr)
    !//copy each slice particle to the beam bunch 
    ipt = 0
    do i = 1, nslice
       do j = 1, countlc(i)
          ipt = ipt + 1
          Pts(:,ipt) = beamslice(:,j,i)
       enddo
    enddo

  end subroutine bbmSlcfixnew_Beambeam


  subroutine slicernew(Pts1in,nptlc,nslice,nslice2,maxsptcl,nz,zmin,zmax,beamslice,&
       weight,zslice,zslice2,hzslice,hzslice2,nptot,count,countlc,myidy,npyhalf,commcol,commrow)
    implicit none
    include 'mpif.h'
    integer, parameter :: Ndiv = 10
    integer, intent(in) :: nptlc,nslice,nslice2,maxsptcl,nptot,nz,myidy,&
         commcol,commrow,npyhalf
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, intent(in) :: zmin,zmax
    double precision, dimension(7,maxsptcl,nslice),&
         intent(inout) :: beamslice
    double precision, dimension(nslice),intent(out) :: weight,zslice,hzslice
    double precision, dimension(nslice2),intent(out) :: zslice2,hzslice2
    integer, dimension(nslice), intent(out) :: count,countlc
    double precision, dimension(nslice+nslice2)  :: fcount,fcountlc
    double precision, dimension(0:nslice) :: rangez
    double precision, dimension(2*nz) :: rhozlc,rhoz
    double precision, dimension(nz*Ndiv) :: rhoznew
    integer :: i,ierr,ii,iz,iz1,ip,op,mydes,mysou
    double precision :: ab,sumint,hnew1,hnew2,interval0,hz,hzi
    integer status(MPI_STATUS_SIZE)
    real*8 :: sum0,sum1,sumip,z0,z1,zip
    integer :: nsendata,ip0

    op = 1
    if(nslice.gt.1) then
       !//find the charge density along z
       hz = (zmax-zmin)/(nz-1)
       hzi = 1.0/hz
       rhozlc = 0.0
       if(myidy.lt.npyhalf) then
          do i = 1, nptlc
             iz=int((Pts1in(5,i)-zmin)*hzi + 1)
             ab=(zmin-Pts1in(5,i)+iz*hz)*hzi
             iz1 = iz + 1
             rhozlc(iz) = rhozlc(iz) + ab
             rhozlc(iz1) = rhozlc(iz1) + (1.0-ab)
          enddo
       else
          do i = 1, nptlc
             iz=int((Pts1in(5,i)-zmin)*hzi + 1) !//nz is for group 2
             ab=(zmin-Pts1in(5,i)+iz*hz)*hzi
             iz1 = iz + 1 
             rhozlc(iz+nz) = rhozlc(iz+nz) + ab
             rhozlc(iz1+nz) = rhozlc(iz1+nz) + (1.0-ab)
          enddo
       endif
       rhozlc = rhozlc/nptot

       call MPI_ALLREDUCE(rhozlc,rhoz,2*nz,MPI_DOUBLE_PRECISION,MPI_SUM,mpicommwd,ierr)

       !//find total density, should be 1
       sumint = 0.0
       if(myidy.lt.npyhalf) then
          do i = 1, nz
             sumint = sumint + rhoz(i)
          enddo
       else
          do i = 1, nz
             sumint = sumint + rhoz(i+nz)
          enddo
       endif

       !//subdivided the charge density for finer resolution.
       if(myidy.lt.npyhalf) then
          do i = 1, nz*Ndiv
             ii = (i-1)/Ndiv + 1
             rhoznew(i) = rhoz(ii)/Ndiv
          enddo
       else
          do i = 1, nz*Ndiv
             ii = (i-1)/Ndiv + 1 + nz
             rhoznew(i) = rhoz(ii)/Ndiv
          enddo
       endif
       hnew1 = hz/(Ndiv-1)/2
       hnew2 = hz/Ndiv

       interval0 = sumint/nslice
       sumint = 0.0

       !//find the range of each slice
       !          do i = 1, Ndiv*nz
       !            sumint = sumint + rhoznew(i)
       !            ip = sumint/interval0
       !            if(ip.ne.nslice) then
       !              if(i.le.Ndiv) then
       !                rangez(ip+1) = zmin + (i-1)*hnew1 + 0.5*hnew1
       !              else if((i.gt.Ndiv).and.(i.le.Ndiv*nz-Ndiv)) then
       !                rangez(ip+1) = zmin +0.5*hz+(i-Ndiv)*hnew2 + 0.5*hnew2
       !              else
       !                rangez(ip+1) = zmin + (nz-1)*hz - 0.5*hz + &
       !                                 (i+Ndiv-Ndiv*nz)*hz/Ndiv
       !              endif
       !            endif
       !          enddo

       ip0 = 0
       if(myidy.lt.npyhalf) then
          do i = 1, nz
             sumint = sumint + rhoz(i)
             ip = int(sumint/interval0)
             if((ip.ne.ip0).and.(ip.lt.nslice)) then
                z1 = zmin + (i-1)*hz
                z0 = zmin + (i-2)*hz
                sum1 = sumint
                sum0 = sumint-rhoz(i)
                sumip = ip*interval0
                zip = z0 + (sumip-sum0)*(z1-z0)/(sum1-sum0)
                rangez(ip) = zip
                ip0 = ip
             endif
          enddo
       else
          do i = 1, nz
             sumint = sumint + rhoz(i+nz)
             ip = int(sumint/interval0)
             if((ip.ne.ip0).and.(ip.lt.nslice)) then
                z1 = zmin + (i-1)*hz
                z0 = zmin + (i-2)*hz
                sum1 = sumint
                sum0 = sumint-rhoz(i+nz)
                sumip = ip*interval0
                zip = z0 + (sumip-sum0)*(z1-z0)/(sum1-sum0)
                rangez(ip) = zip
                ip0 = ip
             endif
          enddo
       endif

       rangez(nslice) = zmax
       rangez(0) = zmin


       !//find the # of particle per slice and put the particles into
       !//single slice. Here, zmax corresponds to slice # 1.
       do i = 1, nslice
          countlc(i) = 0
       enddo
       do i = 1,nptlc
          if(Pts1in(5,i).gt.rangez(nslice-1)) then
             !countlc(nslice) = countlc(nslice) + 1
             countlc(1) = countlc(1) + 1
             if(countlc(1).gt.maxsptcl) then
                print*,"overflow local slice limit in nslice: ",countlc(1)
                stop
             endif
             beamslice(1:6,countlc(1),1) = Pts1in(:,i)
             beamslice(7,countlc(1),1) = i
          else if(Pts1in(5,i).le.rangez(1)) then
             !countlc(1) = countlc(1) + 1
             countlc(nslice) = countlc(nslice) + 1
             if(countlc(nslice).gt.maxsptcl) then
                print*,"overflow local slice limit in nslice: ",countlc(nslice)
                stop
             endif
             beamslice(1:6,countlc(nslice),nslice) = Pts1in(:,i)
             beamslice(7,countlc(nslice),nslice) = i
          else
             do ii = 2,nslice-1
                if((Pts1in(5,i).gt.rangez(ii-1)).and.&
                     (Pts1in(5,i).le.rangez(ii)) ) then
                   countlc(nslice-ii+1) = countlc(nslice-ii+1) + 1
                   if(countlc(nslice-ii+1).gt.maxsptcl) then
                      print*,"overflow local slice limit in ii: ",ii,&
                           countlc(nslice-ii+1)
                      stop
                   endif
                   beamslice(1:6,countlc(nslice-ii+1),nslice-ii+1)=Pts1in(:,i)
                   beamslice(7,countlc(nslice-ii+1),nslice-ii+1)=i
                endif
             enddo
          endif
       enddo

       fcountlc = 0
       if(myidy.lt.npyhalf) then
          do i = 1, nslice
             fcountlc(i) = countlc(i)
          enddo
       else
          do i = 1, nslice
             fcountlc(i+nslice2) = countlc(i)
          enddo
       endif
       nsendata = nslice+nslice2
       call MPI_ALLREDUCE(fcountlc,fcount,nsendata,MPI_DOUBLE_PRECISION,&
            MPI_SUM,mpicommwd,ierr)
       if(myidy.lt.npyhalf) then
          do i = 1, nslice
             count(i) = int(fcount(i) + 0.0001) !//0.0001 for roundoff
             weight(i) = fcount(i)/nptot
          enddo
       else
          do i = 1, nslice
             count(i) = int(fcount(i+nslice2) + 0.0001) !//0.0001 for roundoff
             weight(i) = fcount(i+nslice2)/nptot
          enddo
       endif

       do i = 1, nslice
          !zslice(i) = (rangez(i-1)+rangez(i))/2 !//this is a terrible bug.
          zslice(nslice-i+1) = (rangez(i-1)+rangez(i))/2
          hzslice(nslice-i+1) = rangez(i)-rangez(i-1)
       enddo

       !//send local zslice to the another group zslice2
       if(myidy.lt.npyhalf) then
          mydes = myidy + npyhalf
          mysou = myidy + npyhalf
       else
          mydes = myidy - npyhalf
          mysou = myidy - npyhalf
       endif
       call MPI_SEND(zslice,nslice,MPI_DOUBLE_PRECISION,mydes,1,commcol,&
            ierr)
       call MPI_RECV(zslice2,nslice2,MPI_DOUBLE_PRECISION,mysou,1,commcol,&
            status,ierr)
       call MPI_SEND(hzslice,nslice,MPI_DOUBLE_PRECISION,mydes,1,commcol,&
            ierr)
       call MPI_RECV(hzslice2,nslice2,MPI_DOUBLE_PRECISION,mysou,1,commcol,&
            status,ierr)

    else !//single slice case !not sure whether it is correct or not if one beam has 1 slice the other beam has N slice
       do i = 1, nptlc
          beamslice(1:6,i,nslice) = Pts1in(:,i)
          beamslice(7,i,nslice) = i
       enddo
       weight(nslice) = 1.0
       rangez(nslice) = zmax
       rangez(0) = zmin
       zslice(1) = 0.0
       zslice2(1) = 0.0
       countlc = nptlc
       !call MPI_ALLREDUCE(countlc,count,nslice,MPI_INTEGER,&
       !                   MPI_SUM,mpicommwd,ierr)
       count(1) = nptot
       hzslice(1) = zmax - zmin
       hzslice2(1) = zmax - zmin !not correct but not used
    endif

  end subroutine slicernew



  !//zslice is the z location of itself and zslice2 is
  !the z location of the opposite beam
  subroutine slicer(Pts1in,nptlc,nslice,nslice2,maxsptcl,&
       nz,zmin,&
       zmax,beamslice,weight,zslice,zslice2,nptot,count,countlc,myidy,&
       npyhalf,commcol,commrow)
    implicit none
    include 'mpif.h'
    integer, parameter :: Ndiv = 10
    integer, intent(in) :: nptlc,nslice,nslice2,maxsptcl,nptot,nz,myidy,&
         commcol,commrow,npyhalf
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, intent(in) :: zmin,zmax
    double precision, dimension(6,maxsptcl,nslice),&
         intent(inout) :: beamslice
    double precision, dimension(nslice),intent(out) :: weight,zslice
    double precision, dimension(nslice2),intent(out) :: zslice2
    integer, dimension(nslice), intent(out) :: count,countlc
    double precision, dimension(2*nslice)  :: fcount,fcountlc
    double precision, dimension(0:nslice) :: rangez
    double precision, dimension(2*nz) :: rhozlc,rhoz
    double precision, dimension(nz*Ndiv) :: rhoznew
    integer :: i,ierr,ii,iz,iz1,ip,op,mydes,mysou
    double precision :: ab,sumint,hnew1,hnew2,interval0,hz,hzi
    integer status(MPI_STATUS_SIZE)

    op = 1
    if(nslice.gt.1) then
       !//find the charge density along z
       hz = (zmax-zmin)/(nz-1)
       hzi = 1.0/hz
       rhozlc = 0.0
       if(myidy.lt.npyhalf) then
          do i = 1, nptlc
             iz=int((Pts1in(5,i)-zmin)*hzi + 1)
             ab=(zmin-Pts1in(5,i)+iz*hz)*hzi
             iz1 = iz + 1
             rhozlc(iz) = rhozlc(iz) + ab
             rhozlc(iz1) = rhozlc(iz1) + (1.0-ab)
          enddo
       else
          do i = 1, nptlc
             iz=int((Pts1in(5,i)-zmin)*hzi + 1) !//nz is for group 2
             ab=(zmin-Pts1in(5,i)+iz*hz)*hzi
             iz1 = iz + 1 
             rhozlc(iz+nz) = rhozlc(iz+nz) + ab
             rhozlc(iz1+nz) = rhozlc(iz1+nz) + (1.0-ab)
          enddo
       endif
       rhozlc = rhozlc/nptot

       call MPI_ALLREDUCE(rhozlc,rhoz,2*nz,MPI_DOUBLE_PRECISION,&
            MPI_SUM,mpicommwd,ierr)

       !//find totol density, should be 1
       sumint = 0.0
       if(myidy.lt.npyhalf) then
          do i = 1, nz
             sumint = sumint + rhoz(i)
          enddo
       else
          do i = 1, nz
             sumint = sumint + rhoz(i+nz)
          enddo
       endif

       !//subdivided the charge density for finer resolution.
       if(myidy.lt.npyhalf) then
          do i = 1, nz*Ndiv
             ii = (i-1)/Ndiv + 1
             rhoznew(i) = rhoz(ii)/Ndiv
          enddo
       else
          do i = 1, nz*Ndiv
             ii = (i-1)/Ndiv + 1 + nz
             rhoznew(i) = rhoz(ii)/Ndiv
          enddo
       endif
       hnew1 = hz/(Ndiv-1)/2
       hnew2 = hz/Ndiv

       interval0 = sumint/nslice
       sumint = 0.0
       !//find the range of each slice
       do i = 1, Ndiv*nz
          sumint = sumint + rhoznew(i)
          ip = int(sumint/interval0)
          if(ip.ne.nslice) then
             if(i.le.Ndiv) then
                rangez(ip+1) = zmin + (i-1)*hnew1 + 0.5*hnew1
             else if((i.gt.Ndiv).and.(i.le.Ndiv*nz-Ndiv)) then
                rangez(ip+1) = zmin +0.5*hz+(i-Ndiv)*hnew2 + 0.5*hnew2
             else
                rangez(ip+1) = zmin + (nz-1)*hz - 0.5*hz + &
                     (i+Ndiv-Ndiv*nz)*hz/Ndiv
             endif
          endif
       enddo
       rangez(nslice) = zmax
       rangez(0) = zmin


       !//find the # of particle per slice and put the particles into
       !//single slice. Here, zmax corresponds to slice # 1.
       do i = 1, nslice
          countlc(i) = 0
       enddo
       do i = 1,nptlc
          if(Pts1in(5,i).gt.rangez(nslice-1)) then
             !countlc(nslice) = countlc(nslice) + 1
             countlc(1) = countlc(1) + 1
             if(countlc(1).gt.maxsptcl) then
                print*,"overflow local slice limit in nslice: ",countlc(1)
                stop
             endif
             beamslice(:,countlc(1),1) = Pts1in(:,i)
          else if(Pts1in(5,i).le.rangez(1)) then
             !countlc(1) = countlc(1) + 1
             countlc(nslice) = countlc(nslice) + 1
             if(countlc(nslice).gt.maxsptcl) then
                print*,"overflow local slice limit in nslice: ",countlc(nslice)
                stop
             endif
             beamslice(:,countlc(nslice),nslice) = Pts1in(:,i)
          else
             do ii = 2,nslice-1
                if((Pts1in(5,i).gt.rangez(ii-1)).and.&
                     (Pts1in(5,i).le.rangez(ii)) ) then
                   countlc(nslice-ii+1) = countlc(nslice-ii+1) + 1
                   if(countlc(nslice-ii+1).gt.maxsptcl) then
                      print*,"overflow local slice limit in ii: ",ii,&
                           countlc(nslice-ii+1)
                      stop
                   endif
                   beamslice(:,countlc(nslice-ii+1),nslice-ii+1)=Pts1in(:,i)
                endif
             enddo
          endif
       enddo

       fcountlc = 0
       if(myidy.lt.npyhalf) then
          do i = 1, nslice
             fcountlc(i) = countlc(i)
          enddo
       else
          do i = 1, nslice
             fcountlc(i+nslice) = countlc(i)
          enddo
       endif
       call MPI_ALLREDUCE(fcountlc,fcount,2*nslice,MPI_DOUBLE_PRECISION,&
            MPI_SUM,mpicommwd,ierr)

       if(myidy.lt.npyhalf) then
          do i = 1, nslice
             count(i) = int(fcount(i) + 0.0001) !//0.0001 for roundoff
             weight(i) = fcount(i)/nptot
          enddo
       else
          do i = 1, nslice
             count(i) = int(fcount(i+nslice) + 0.0001) !//0.0001 for roundoff
             weight(i) = fcount(i+nslice)/nptot
          enddo
       endif

       do i = 1, nslice
          !zslice(i) = (rangez(i-1)+rangez(i))/2 !//this is a terrible bug.
          zslice(nslice-i+1) = (rangez(i-1)+rangez(i))/2
       enddo

       !//send local zslice to the another group zslice2
       if(myidy.lt.npyhalf) then
          mydes = myidy + npyhalf
          mysou = myidy + npyhalf
       else
          mydes = myidy - npyhalf
          mysou = myidy - npyhalf
       endif
       call MPI_SEND(zslice,nslice,MPI_DOUBLE_PRECISION,mydes,1,commcol,&
            ierr)
       call MPI_RECV(zslice2,nslice2,MPI_DOUBLE_PRECISION,mysou,1,commcol,&
            status,ierr)

    else !//single slice case
       do i = 1, nptlc
          beamslice(:,i,nslice) = Pts1in(:,i)
       enddo
       weight(nslice) = 1.0
       rangez(nslice) = zmax
       rangez(0) = zmin
       zslice(1) = 0.0
       zslice2(1) = 0.0
       countlc = nptlc
       !call MPI_ALLREDUCE(countlc,count,nslice,MPI_INTEGER,&
       !                   MPI_SUM,mpicommwd,ierr)
       count(1) = nptot
    endif

  end subroutine slicer



  subroutine slicerold(Pts1in,nptlc,nslice,maxsptcl,countlcin,zmin,&
       hzslice,beamslice,weight,zslice,nptot,count)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc,nslice,maxsptcl,nptot
    integer, dimension(nslice), intent(in) :: countlcin
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, intent(in) :: zmin,hzslice
    double precision, dimension(6,maxsptcl,nslice),&
         intent(inout) :: beamslice
    double precision, dimension(nslice),intent(out) :: weight,zslice
    integer, dimension(nslice), intent(out) :: count
    integer :: i,ierr,islice
    integer, dimension(nslice) :: countlc

    call MPI_ALLREDUCE(countlcin,count,nslice,MPI_INTEGER,&
         MPI_SUM,mpicommwd,ierr)

    do i = 1, nslice
       weight(i) = dble(count(i))/nptot
    enddo

    if(nslice.gt.1) then
       do i = 1, nslice
          zslice(i) = zmin + (i-1)*hzslice + hzslice/2
       enddo
    else
       zslice(1) = 0.0
    endif

    do i = 1, nslice
       countlc(i) = 0
    enddo
    do i = 1, nptlc
       islice = int((Pts1in(5,i) - zmin)/hzslice + 1)
       countlc(islice) = countlc(islice) + 1
       beamslice(:,countlc(islice),islice) = Pts1in(:,i)
    enddo

  end subroutine slicerold



  !zslice is the z location of itself and zslice2 is
  !the z location of the opposite beam
  subroutine slicerold2(Pts1in,nptlc,nslice,nslice2,maxsptcl,nz,zmin,&
       zmax,beamslice,weight,zslice,zslice2,nptot,count,countlc,myidy,&
       npyhalf,commcol,commrow)
    implicit none
    include 'mpif.h'
    integer, parameter :: Ndiv = 10
    integer, intent(in) :: nptlc,nslice,nslice2,maxsptcl,nptot,nz,myidy,&
         commcol,commrow,npyhalf
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, intent(in) :: zmin,zmax
    double precision, dimension(6,maxsptcl,nslice),&
         intent(inout) :: beamslice
    double precision, dimension(nslice),intent(out) :: weight,zslice
    double precision, dimension(nslice2),intent(out) :: zslice2
    integer, dimension(nslice), intent(out) :: count,countlc
    double precision, dimension(nslice)  :: fcount,fcountlc
    double precision, dimension(0:nslice) :: rangez
    double precision, dimension(nz) :: rhozlc,rhoz
    double precision, dimension(nz*Ndiv) :: rhoznew
    integer :: i,ierr,ii,iz,iz1,ip,op,mydes,mysou
    double precision :: ab,sumint,hnew1,hnew2,interval0,hz,hzi
    integer status(MPI_STATUS_SIZE)

    op = 1
    if(nslice.gt.1) then
       hz = (zmax-zmin)/(nz-1)
       hzi = 1.0/hz
       rhozlc = 0.0
       do i = 1, nptlc
          iz=int((Pts1in(5,i)-zmin)*hzi + 1)
          ab=(zmin-Pts1in(5,i)+iz*hz)*hzi
          iz1 = iz + 1
          rhozlc(iz) = rhozlc(iz) + ab
          rhozlc(iz1) = rhozlc(iz1) + (1.0-ab)
       enddo
       rhozlc = rhozlc/nptot

       call MPI_ALLREDUCE(rhozlc,rhoz,nz,MPI_DOUBLE_PRECISION,&
            MPI_SUM,commrow,ierr)
       rhozlc = rhoz
       call allreduce(rhozlc,rhoz,nz,npyhalf,myidy,commcol,op)

       sumint = 0.0
       do i = 1, nz
          sumint = sumint + rhoz(i)
       enddo

       do i = 1, nz*Ndiv
          ii = (i-1)/Ndiv + 1
          rhoznew(i) = rhoz(ii)/Ndiv
       enddo
       hnew1 = hz/(Ndiv-1)/2
       hnew2 = hz/Ndiv

       interval0 = sumint/nslice
       sumint = 0.0

       do i = 1, Ndiv*nz
          sumint = sumint + rhoznew(i)
          ip = int(sumint/interval0)
          if(ip.ne.nslice) then
             if(i.le.Ndiv) then
                rangez(ip+1) = zmin + (i-1)*hnew1 + 0.5*hnew1
             else if((i.gt.Ndiv).and.(i.le.Ndiv*nz-Ndiv)) then
                rangez(ip+1) = zmin +0.5*hz+(i-Ndiv)*hnew2 + 0.5*hnew2
             else
                rangez(ip+1) = zmin + (nz-1)*hz - 0.5*hz + &
                     (i+Ndiv-Ndiv*nz)*hz/Ndiv
             endif
          endif
       enddo
       rangez(nslice) = zmax
       rangez(0) = zmin


       do i = 1, nslice
          countlc(i) = 0
       enddo
       do i = 1,nptlc
          if(Pts1in(5,i).gt.rangez(nslice-1)) then
             !countlc(nslice) = countlc(nslice) + 1
             countlc(1) = countlc(1) + 1
             if(countlc(1).gt.maxsptcl) then
                print*,"overflow local slice limit in nslice: ",countlc(1)
                stop
             endif
             beamslice(:,countlc(1),1) = Pts1in(:,i)
          else if(Pts1in(5,i).le.rangez(1)) then
             !countlc(1) = countlc(1) + 1
             countlc(nslice) = countlc(nslice) + 1
             if(countlc(nslice).gt.maxsptcl) then
                print*,"overflow local slice limit in nslice: ",countlc(nslice)
                stop
             endif
             beamslice(:,countlc(nslice),nslice) = Pts1in(:,i)
          else
             do ii = 2,nslice-1
                if((Pts1in(5,i).gt.rangez(ii-1)).and.&
                     (Pts1in(5,i).le.rangez(ii)) ) then
                   countlc(nslice-ii+1) = countlc(nslice-ii+1) + 1
                   if(countlc(nslice-ii+1).gt.maxsptcl) then
                      print*,"overflow local slice limit in ii: ",ii,&
                           countlc(nslice-ii+1)
                      stop
                   endif
                   beamslice(:,countlc(nslice-ii+1),nslice-ii+1)=Pts1in(:,i)
                endif
             enddo
          endif
       enddo

       call MPI_ALLREDUCE(countlc,count,nslice,MPI_INTEGER,&
            MPI_SUM,commrow,ierr)
       fcountlc = count
       call allreduce(fcountlc,fcount,nslice,npyhalf,myidy,commcol,op)

       !The 0.1 is try to avoid the error in double precision
       count = int(fcount + 0.1)

       do i = 1, nslice
          weight(i) = fcount(i)/nptot
       enddo

       do i = 1, nslice
          zslice(i) = (rangez(i-1)+rangez(i))/2
       enddo

       if(myidy.lt.npyhalf) then
          mydes = myidy + npyhalf
          mysou = myidy + npyhalf
       else
          mydes = myidy - npyhalf
          mysou = myidy - npyhalf
       endif
       call MPI_SEND(zslice,nslice,MPI_DOUBLE_PRECISION,mydes,1,commcol,&
            ierr)
       call MPI_RECV(zslice2,nslice2,MPI_DOUBLE_PRECISION,mysou,1,commcol,&
            status,ierr)

    else
       do i = 1, nptlc
          beamslice(:,i,nslice) = Pts1in(:,i)
       enddo
       weight(nslice) = 1.0
       rangez(nslice) = zmax
       rangez(0) = zmin
       zslice(1) = 0.0
       zslice2(1) = 0.0
       countlc = nptlc
       !call MPI_ALLREDUCE(countlc,count,nslice,MPI_INTEGER,&
       !                   MPI_SUM,mpicommwd,ierr)
       count(1) = nptot
    endif

  end subroutine slicerold2



  subroutine driftCPfixnew(ray,sign,nslice,nslice2,myidy,npyhalf,&
       istart,jstart,jcross,zslice,zslice2,hzslice2,countlc,maxnpt,sign2)
    integer, intent(in) :: sign,nslice,myidy,npyhalf,istart,jstart,&
         jcross,nslice2,maxnpt,sign2
    integer, dimension(nslice) :: countlc
    double precision, dimension(nslice) :: zslice
    double precision, dimension(nslice2) :: zslice2,hzslice2
    double precision :: distance
    double precision, dimension(7,maxnpt,nslice), intent(inout) :: ray
    integer :: i,j,i1,j1

    do j = 1, jcross !//loop through the number of overlaped slices.
       if(myidy.lt.npyhalf) then
          i1 = istart - j
          j1 = jstart + j
       else
          i1 = jstart + j
          j1 = istart - j
       endif
       !//the drift distance from each beam has the same value but
       !//opposite sign. This means that if one beam slice particle moves
       !//forward, the other beam slice particle is going to move
       !//backward.
       !distance = (zslice(i1)+sign2*0.5*hzslice(i1)-zslice2(j1))/2
       distance = ( zslice(i1)- (sign2*0.5*hzslice2(j1)+zslice2(j1)) )/2
       do i = 1, countlc(i1)
          ray(1,i,i1) = ray(1,i,i1) + distance*sign*ray(2,i,i1)
          ray(3,i,i1) = ray(3,i,i1) + distance*sign*ray(4,i,i1)
       enddo
    enddo

  end subroutine driftCPfixnew



  !//drift slice particles to the collision points for charge density 
  !//(here, all particles drift zj - zi) and find the local range of the beam.
  !//for field interpolation, the drift is delta zi - zj. It is is different
  !//for different particles.
  subroutine driftCPfix(ray,sign,nslice,nslice2,myidy,npyhalf,&
       istart,jstart,jcross,zslice,zslice2,countlc,maxnpt)
    integer, intent(in) :: sign,nslice,myidy,npyhalf,istart,jstart,&
         jcross,nslice2,maxnpt
    integer, dimension(nslice) :: countlc
    double precision, dimension(nslice) :: zslice
    double precision, dimension(nslice2) :: zslice2
    double precision :: distance
    double precision, dimension(6,maxnpt,nslice), intent(inout) :: ray
    integer :: i,j,i1,j1

    do j = 1, jcross !//loop through the number of overlaped slices.
       if(myidy.lt.npyhalf) then
          i1 = istart - j
          j1 = jstart + j
       else
          i1 = jstart + j
          j1 = istart - j
       endif
       !//the drift distance from each beam has the same value but
       !//opposite sign. This means that if one beam slice particle moves
       !//forward, the other beam slice particle is going to move
       !//backward.
       distance = (zslice(i1)-zslice2(j1))/2
       do i = 1, countlc(i1)
          ray(1,i,i1) = ray(1,i,i1) + distance*sign*ray(2,i,i1)
          ray(3,i,i1) = ray(3,i,i1) + distance*sign*ray(4,i,i1)
       enddo
    enddo

  end subroutine driftCPfix



  !//drift slice particles to the collision field interpolation point.
  subroutine driftCP(ray,sign,nslice,nslice2,myidy,npyhalf,&
       istart,jstart,jcross,zslice,zslice2,countlc,maxnpt)
    integer, intent(in) :: sign,nslice,myidy,npyhalf,istart,jstart,&
         jcross,nslice2,maxnpt
    integer, dimension(nslice) :: countlc
    double precision, dimension(nslice) :: zslice
    double precision, dimension(nslice2) :: zslice2
    double precision :: distance
    double precision, dimension(7,maxnpt,nslice), intent(inout) :: ray
    integer :: i,j,i1,j1

    do j = 1, jcross !//loop through the number of overlaped slices.
       if(myidy.lt.npyhalf) then
          i1 = istart - j
          j1 = jstart + j
       else
          i1 = jstart + j
          j1 = istart - j
       endif
       !//the drift distance from each beam has the same value but
       !//opposite sign. This means that if one beam slice particle moves
       !//forward, the other beam slice particle is going to move
       !//backward.
       do i = 1, countlc(i1)
          distance = (ray(5,i,i1)-zslice2(j1))/2
          ray(1,i,i1) = ray(1,i,i1) + distance*sign*ray(2,i,i1)
          ray(3,i,i1) = ray(3,i,i1) + distance*sign*ray(4,i,i1)
       enddo
    enddo

  end subroutine driftCP



  ! drift each slice to the IP point.
  subroutine driftIP(ray,distance,sign,slnptlc)
    integer, intent(in) :: sign,slnptlc
    double precision, intent(in) :: distance
    double precision, dimension(6,slnptlc), intent(inout) :: ray
    integer :: i

    do i = 1, slnptlc
       ray(1,i) = ray(1,i) + distance*sign*ray(2,i)
       ray(3,i) = ray(3,i) + distance*sign*ray(4,i)
    enddo

  end subroutine driftIP



  ! drift each slice to the IP point.
  subroutine driftIPold(ray,distance,range,sign,slnptlc)
    integer, intent(in) :: sign,slnptlc
    double precision, intent(in) :: distance
    double precision, dimension(6,slnptlc), intent(inout) :: ray
    double precision, dimension(4) :: range
    integer :: i

    range(1) = 1.0e10
    range(2) = -1.0e10
    range(3) = 1.0e10
    range(4) = -1.0e10
    do i = 1, slnptlc
       ray(1,i) = ray(1,i) + distance*sign*ray(2,i)
       if(range(1).gt.ray(1,i)) then
          range(1) = ray(1,i)
       else if(range(2).le.ray(1,i)) then
          range(2) = ray(1,i)
       else
       endif
       ray(3,i) = ray(3,i) + distance*sign*ray(4,i)
       if(range(3).gt.ray(3,i)) then
          range(3) = ray(3,i)
       else if(range(4).le.ray(3,i)) then
          range(4) = ray(3,i)
       else
       endif
    enddo

  end subroutine driftIPold



  !//--------------------------------------------------------------------------
  !//find the strong-strong beam-beam kick in a 2d one slice model.
  subroutine bb1Slc_Beambeam(nprocrow,npyhalf,innx,inny,&
       nptlc,myidy,Pts1,grid2d,Ageom,&
       nptot,coef,nytot,commcol,shift21)
    implicit none
    integer, intent(in) :: nprocrow,npyhalf,innx,&
         myidy,nptot,nytot,commcol
    integer, intent(inout) :: inny
    integer, intent(inout) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1
    type (Pgrid2d), intent(in) :: grid2d
    type (CompDom), intent(inout) :: Ageom
    double precision, intent(in) :: coef
    double precision, dimension(2), intent(in) :: shift21
    double precision, dimension(2) :: shift12,shift
    double precision, dimension(4) :: range
    double precision, dimension(6) :: glrange
    double precision, dimension(3) :: msize
    double precision, dimension(innx,nytot) :: rho,phiout
    double precision, dimension(innx,inny) :: rholc
    double precision :: xmin,ymin,hx,hy,weight
    integer, dimension(2,0:nprocrow-1,0:npyhalf-1) :: table
    integer :: innyr,k,ierr

    shift12 = -shift21
    innyr = inny
    range = 0.0
    range(1) = Pts1(1,1)
    range(2) = Pts1(1,1)
    range(3) = Pts1(3,1)
    range(4) = Pts1(3,1)
    if(myidy.lt.npyhalf) then
       shift = shift21 
    else
       shift = shift12 
    endif

    !//here, we have used equal area weight
    weight = 1.0d0

    call MPI_BARRIER(mpicommwd,ierr)

    !//find the local range of slice particle after drift.
    do k = 1, nptlc
       if(range(1).gt.Pts1(1,k)) then
          range(1) = Pts1(1,k)
       else if(range(2).le.Pts1(1,k)) then
          range(2) = Pts1(1,k)
       else
       endif
       if(range(3).gt.Pts1(3,k)) then
          range(3) = Pts1(3,k)
       else if(range(4).le.Pts1(3,k)) then
          range(4) = Pts1(3,k)
       else
       endif
    enddo

    if(myidy.ge.npyhalf) then
       range(1) = range(1) - shift21(1)
       range(2) = range(2) - shift21(1)
       range(3) = range(3) - shift21(2)
       range(4) = range(4) - shift21(2)
    endif

    !//update the computational domain after drift
    !call update2d_CompDom(Ageom,range1,range2,grid2d)
    call update2d_CompDom(Ageom,range,grid2d)
    call getmsize_CompDom(Ageom,msize)
    hx = msize(1)
    hy = msize(2)
    call getrange_CompDom(Ageom,glrange)
    if(myidy.lt.npyhalf) then
       xmin = glrange(1)
       ymin = glrange(3)
    else
       xmin = glrange(1) + shift21(1)
       ymin = glrange(3) + shift21(2)
    endif
    call getlctabnm_CompDom(Ageom,table)

    !//get the charge density on the grid.
    call deposit2dtsc(nptlc,nptlc,innx,nytot,hx,hy,xmin,&
                   ymin,Pts1,rho)
    !call deposit2d(nptlc,nptlc,innx,nytot,hx,hy,xmin,&
    !     ymin,Pts1,rho)
    !//sum together the contributions from the particles outside the subdomain
    call guardsum2d(rho,rholc,innx,inny,nytot,npyhalf,myidy,commcol)
    rholc = rholc/nptot

    !//solve the 2d Poisson equation in parallel.
    call fieldsolver2d1slc(rholc,innx,inny,nprocrow,npyhalf,&
         innyr,nytot,table,hx,hy,myidy,commcol,shift)

    !//collect the potential from local domain to the global for all PEs.
    call guardexch2d(rholc,phiout,innx,inny,nytot,npyhalf,myidy,commcol)

    !//calculate the field and interpolate back to opposite particles.
    !call scatter2d(nptlc,nptlc,innx,&
    !     hx,hy,xmin,ymin,Pts1,phiout,weight,coef,nytot)
    call scatter2dtsc(nptlc,nptlc,innx,&
         hx,hy,xmin,ymin,Pts1,phiout,weight,coef,nytot,myidy,npyhalf)
    !call scatter2dtsc(nptlc,nptlc,innx,&
    !hx,hy,xmin,ymin,Pts1,phiout,weight,coef,nytot,myidy,npyhalf)

    call MPI_BARRIER(mpicommwd,ierr)

  end subroutine bb1Slc_Beambeam



  !//find the strong-strong beam-beam kick in a 2d one slice model
  !//using a fixed computational domain, i.e.
  !//reuse Green function of the domain, no shift.
  subroutine bb1Slcfix_Beambeam(nprocrow,npyhalf,innx,inny,nptlc,&
       myidy,Pts1,grid2d,nptot,coef,nytot,commcol,nxpylc2,grn,&
       xpystable,pytable,hx,hy,xmin,ymin)
    implicit none
    integer, intent(in) :: nprocrow,npyhalf,innx,&
         myidy,nptot,nytot,commcol,nxpylc2
    integer, intent(inout) :: inny
    integer, intent(inout) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1
    type (Pgrid2d), intent(in) :: grid2d
    double precision, intent(in) :: coef
    double precision, dimension(innx,nytot) :: rho,phiout
    double precision, dimension(innx,inny) :: rholc
    double precision :: xmin,ymin,hx,hy,weight
    integer :: ierr
    double complex, dimension (2*nytot,nxpylc2) :: grn
    integer, dimension(0:npyhalf-1) :: xpystable,pytable

    !//here, we have used equal area weight
    weight = 1.0

    call MPI_BARRIER(mpicommwd,ierr)

    !//get the charge density on the grid.
    call deposit2d(nptlc,nptlc,innx,nytot,hx,hy,xmin,&
         ymin,Pts1,rho)
    !//sum together the contributions from the particles outside the subdomain
    call guardsum2d(rho,rholc,innx,inny,nytot,npyhalf,myidy,commcol)
    rholc = rholc/nptot

    !//solve the 2d Poisson equation in parallel using the input Green function.
    call fieldsolver2d(rholc,innx,inny,nprocrow,npyhalf,&
         nytot,nxpylc2,xpystable,pytable,hx,hy,myidy,commcol,grn)

    !//collect the potential from local domain to the global for all PEs.
    call guardexch2d(rholc,phiout,innx,inny,nytot,npyhalf,myidy,commcol)

    !//calculate the field and interpolate back to opposite particles.
    call scatter2d(nptlc,nptlc,innx,&
         hx,hy,xmin,ymin,Pts1,phiout,weight,coef,&
         nytot)

    call MPI_BARRIER(mpicommwd,ierr)

  end subroutine bb1Slcfix_Beambeam



  !//find the strong-strong beam-beam kick in a 2d one slice model.
  subroutine bb1Slcold_Beambeam(nprocrow,nproccol,innx,inny,&
       nptlc,myidy,Pts1,grid2d,Ageom,&
       nptot,coef,nytot,commcol,shift21)
    implicit none
    integer, intent(in) :: nprocrow,nproccol,innx,&
         myidy,nptot,nytot,commcol
    integer, intent(inout) :: inny
    integer, intent(inout) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1
    type (Pgrid2d), intent(in) :: grid2d
    type (CompDom), intent(inout) :: Ageom
    double precision, intent(in) :: coef
    double precision, dimension(2), intent(in) :: shift21
    double precision, dimension(2) :: shift12,shift
    double precision, dimension(4) :: range
    double precision, dimension(6) :: glrange
    double precision, dimension(3) :: msize
    double precision, dimension(innx,nytot) :: rho,phiout
    double precision, dimension(innx,inny) :: rholc
    double precision :: xmin,ymin,hx,hy,weight
    integer, dimension(2,0:nprocrow-1,0:nproccol/2-1) :: table
    integer :: innyr,k,ierr,npyhalf

    shift12 = -shift21
    innyr = inny
    range = 0.0
    if(mod(nproccol,2).ne.0) then
       print*,"number of nproccol has to be even"
       stop
    endif
    npyhalf = nproccol/2
    range(1) = Pts1(1,1)
    range(2) = Pts1(1,1)
    range(3) = Pts1(3,1)
    range(4) = Pts1(3,1)
    if(myidy.lt.npyhalf) then
       shift = shift21 
    else
       shift = shift12 
    endif

    !here, we have used equal area weight
    weight = 1.0

    call MPI_BARRIER(mpicommwd,ierr)

    !find the range of slice particle after drift.
    do k = 1, nptlc
       if(range(1).gt.Pts1(1,k)) then
          range(1) = Pts1(1,k)
       else if(range(2).le.Pts1(1,k)) then
          range(2) = Pts1(1,k)
       else
       endif
       if(range(3).gt.Pts1(3,k)) then
          range(3) = Pts1(3,k)
       else if(range(4).le.Pts1(3,k)) then
          range(4) = Pts1(3,k)
       else
       endif
    enddo

    if(myidy.ge.npyhalf) then
       range(1) = range(1) - shift21(1)
       range(2) = range(2) - shift21(1)
       range(3) = range(3) - shift21(2)
       range(4) = range(4) - shift21(2)
    endif

    !update the computational domain after drift
    !call update2d_CompDom(Ageom,range1,range2,grid2d)
    call update2d_CompDom(Ageom,range,grid2d)
    call getmsize_CompDom(Ageom,msize)
    hx = msize(1)
    hy = msize(2)
    call getrange_CompDom(Ageom,glrange)
    if(myidy.lt.npyhalf) then
       xmin = glrange(1)
       ymin = glrange(3)
    else
       xmin = glrange(1) + shift21(1)
       ymin = glrange(3) + shift21(2)
    endif
    call getlctabnm_CompDom(Ageom,table)

    call deposit2d(nptlc,nptlc,innx,nytot,hx,hy,xmin,&
         ymin,Pts1,rho)
    call guardsum2d(rho,rholc,innx,inny,nytot,npyhalf,myidy,commcol)
    rholc = rholc/nptot

    call fieldsolver2d1slc(rholc,innx,inny,nprocrow,npyhalf,&
         innyr,nytot,table,hx,hy,myidy,commcol,shift)

    call guardexch2d(rholc,phiout,innx,inny,nytot,npyhalf,myidy,commcol)


    call scatter2d(nptlc,nptlc,innx,&
         hx,hy,xmin,ymin,Pts1,phiout,weight,coef,&
         nytot)

    call MPI_BARRIER(mpicommwd,ierr)

  end subroutine bb1Slcold_Beambeam



  !//--------------------------------------------------------------------------
  !find the strong-weak beam-beam kick in a 3D slice model.
  !here, beam 1 is strong beam; beam 2 is weak beam.
  !There are three type of weak-strong calculations:
  !1)The kick from strong beam is computed only one turn and
  !  stored for use. This works correctly only with single slice
  !  model. With multiple slice model, it is just an approximation,
  !  which assume that the SC does not change for strong beam.
  !2)The kick from strong beam is computed every turn. 
  !3)The kick from strong beam is computed using a soft-Gaussian
  !  approximation.  
  subroutine bbwkst(nprocrow,nproccol,innx,inny,&
       nslice1,nslice2,nslicephi,nptlc1,nptlc2,myidy,Pts1,Pts2,grid2d,Ageom,&
       zmin1,zmax1,zmin2,zmax2,nptot1,nptot2,&
       coef1,coef2,nytot,commcol,nz1,nz2,xmin,xmax,ymin,ymax,&
       maxPtclSlice1,maxPtclSlice2,shift21,weight1,zslice1,phistrg,&
       flagwk,stepid,sigx,sigy,centx,centy)
    implicit none
    integer, intent(in) :: nprocrow,nproccol,innx,nslice1,&
         nslice2,nslicephi,myidy,nytot,commcol,nz1,nz2,&
         maxPtclSlice1,maxPtclSlice2,inny
    integer, intent(inout) :: nptlc1,nptlc2,nptot1,nptot2
    integer, intent(in) :: flagwk,stepid
    double precision, pointer, dimension(:,:) :: Pts1,Pts2
    type (Pgrid2d), intent(in) :: grid2d
    type (CompDom), intent(inout) :: Ageom
    double precision, intent(in) :: zmin1,zmin2,&
         coef1,coef2,zmax1,zmax2
    double precision, intent(inout) :: xmin,xmax,ymin,ymax
    double precision, dimension(2), intent(in) :: shift21
    double precision, dimension(nslice1), intent(inout) :: weight1,zslice1
    double precision, dimension(nslice1), intent(inout) :: sigx,sigy
    double precision, dimension(nslice1), intent(inout) :: centx,centy
    double precision, dimension(innx,nytot,nslicephi),&
         intent(inout)::phistrg
    double precision, dimension(2) :: shift12
    double precision, dimension(nslice2) :: weight2,zslice2
    double precision, dimension(6) :: range
    double precision, dimension(4) :: range1,range2
    double precision, dimension(3) :: msize
    double precision, dimension(innx,nytot) :: tmp2d
    double precision, dimension(innx,inny) :: rho1lc
    double precision, dimension(6,maxPtclSlice1) :: ray1 
    double precision, dimension(6,maxPtclSlice2) :: ray2     
    !double precision, pointer, dimension(:,:) :: ray1,ray2     
    double precision, dimension(6,maxPtclSlice1,nslice1) :: beamslice1
    double precision, dimension(6,maxPtclSlice2,nslice2) :: beamslice2
    double precision :: sc,hx,hy,xmin2,ymin2,tmpx,tmpy,dist1
    integer, dimension(nslice1) :: countlc1,count1
    integer, dimension(nslice2) :: countlc2,count2
    integer, dimension(2,0:nprocrow-1,0:nproccol-1) :: table
    integer :: i,j,innyr,nptlcslice1,i1,j1,nptlcslice2,sign,ipt,k,ierr,ilost,k0,nptslice1
    double precision, dimension(2) :: sigma,center
    double precision :: sc1, sc2

    shift12 = -shift21
    innyr = inny

    !here, we have used equal area weight
    !copy particles from a bunch into slices. Each slice has
    !roughly same number of particles.
    !For strong beam flagwk=1 or 3, this needs to be done only once.
    if((flagwk.eq.2).or.(flagwk.eq.4).or.(stepid.eq.1)) then
       call slicer1G(Pts1,nptlc1,nslice1,maxPtclSlice1,nz1,zmin1,&
            zmax1,beamslice1,weight1,zslice1,nptot1,count1,countlc1)
    endif

    !get beam slice 2
    !here, we have used equal area weight
    call slicer1G(Pts2,nptlc2,nslice2,maxPtclSlice2,nz2,zmin2,&
         zmax2,beamslice2,weight2,zslice2,nptot2,count2,countlc2)

    call MPI_BARRIER(mpicommwd,ierr)

    do i = 1, nslice1 !loop througth strong beam.
       do j = 1, nslice2 !loop througth weak beam.

          sign = 1
          !find the collision point between two slices
          !in this version, we assume that the collision points
          !does not change a lot during weak strong interaction.
          sc1 = (zslice1(i) - zslice2(j))/2
          sc2 = -sc1

          !copy each slice into a single slice particle array
          !and move to the collision point
          nptlcslice2 = countlc2(j)
          do k = 1, nptlcslice2
             ray2(:,k) = beamslice2(:,k,j)
          enddo

          if(nslice2.gt.1) then
             call driftIP(ray2,sc2,sign,nptlcslice2)
          endif

          !for weak-strong type 1, check weak beam is inside the 
          !computational domain after the 1st step. 
          if( (flagwk.eq.1).and.(stepid.ne.1) ) then
             ilost = 0
             do k0 = 1, nptlcslice2
                k = k0 - ilost
                tmpx = ray2(1,k0) - shift21(1)
                tmpy = ray2(3,k0) - shift21(2)
                if( (tmpx.ge.xmax).or.(tmpx.le.xmin) ) then
                   ilost = ilost + 1
                else if( (tmpy.ge.ymax).or.(tmpy.le.ymin) ) then
                   ilost = ilost + 1
                endif
                ray2(1,k) = ray2(1,k0)
                ray2(2,k) = ray2(2,k0)
                ray2(3,k) = ray2(3,k0)
                ray2(4,k) = ray2(4,k0)
                ray2(5,k) = ray2(5,k0)
                ray2(6,k) = ray2(6,k0)
             enddo
             nptlcslice2 = nptlcslice2 - ilost
             countlc2(j) = nptlcslice2
          endif

          !for weak-strong type 1, enlarge the domain by 120%.
          if( (flagwk.eq.1).and.(stepid.eq.1)) then
             range2(1) = ray2(1,1)
             range2(2) = ray2(1,1)
             range2(3) = ray2(3,1)
             range2(4) = ray2(3,1)
             do k = 1, nptlcslice2
                if(range2(1).gt.ray2(1,k)) then
                   range2(1) = ray2(1,k)
                else if(range2(2).le.ray2(1,k)) then
                   range2(2) = ray2(1,k)
                else
                endif
                if(range2(3).gt.ray2(3,k)) then
                   range2(3) = ray2(3,k)
                else if(range2(4).le.ray2(3,k)) then
                   range2(4) = ray2(3,k)
                else
                endif
             enddo
             range2(1) = range2(1) - shift21(1)
             range2(2) = range2(2) - shift21(1)
             range2(3) = range2(3) - shift21(2)
             range2(4) = range2(4) - shift21(2)
             !enlarge the computational domain by 120%.
             dist1 = range2(2)-range2(1)
             range2(1) = range2(1) - 0.1*dist1
             range2(2) = range2(2) + 0.1*dist1
             dist1 = range2(4)-range2(3)
             range2(3) = range2(3) - 0.1*dist1
             range2(4) = range2(4) + 0.1*dist1
          endif

          !for weak-strong type 2, compute the new domain range.
          if( (flagwk.ne.1) ) then
             range2(1) = ray2(1,1)
             range2(2) = ray2(1,1)
             range2(3) = ray2(3,1)
             range2(4) = ray2(3,1)
             do k = 1, nptlcslice2
                if(range2(1).gt.ray2(1,k)) then
                   range2(1) = ray2(1,k)
                else if(range2(2).le.ray2(1,k)) then
                   range2(2) = ray2(1,k)
                else
                endif
                if(range2(3).gt.ray2(3,k)) then
                   range2(3) = ray2(3,k)
                else if(range2(4).le.ray2(3,k)) then
                   range2(4) = ray2(3,k)
                else
                endif
             enddo
             range2(1) = range2(1) - shift21(1)
             range2(2) = range2(2) - shift21(1)
             range2(3) = range2(3) - shift21(2)
             range2(4) = range2(4) - shift21(2)
          endif

          !find the electric potential from strong beam 1.
          !for weak-strong type 1, this is done only in the 1st step.
          !for weak-strong type 2, this has to be done every step.
          if( (flagwk.ne.2) .or. (stepid.eq.1) ) then
             nptlcslice1 = countlc1(i)
             do k = 1, nptlcslice1
                ray1(:,k) = beamslice1(:,k,i)
             enddo

             !drift strong beam 1 to the location of collision
             if(nslice1.gt.1) then
                call driftIP(ray1,sc1,sign,nptlcslice1)
             endif

             !find the range of slice particle after drift.
             range1(1) = ray1(1,1)
             range1(2) = ray1(1,1)
             range1(3) = ray1(3,1)
             range1(4) = ray1(3,1)
             do k = 1, nptlcslice1
                if(range1(1).gt.ray1(1,k)) then
                   range1(1) = ray1(1,k)
                else if(range1(2).le.ray1(1,k)) then
                   range1(2) = ray1(1,k)
                else
                endif
                if(range1(3).gt.ray1(3,k)) then
                   range1(3) = ray1(3,k)
                else if(range1(4).le.ray1(3,k)) then
                   range1(4) = ray1(3,k)
                else
                endif
             enddo
             !use an enlarged domain for type 1 weak-strong interaction
             if(flagwk.eq.1) then
                !enlarge the computational domain by 120%.
                dist1 = range1(2)-range1(1)
                range1(1) = range1(1) - 0.1*dist1
                range1(2) = range1(2) + 0.1*dist1
                dist1 = range1(4)-range1(3)
                range1(3) = range1(3) - 0.1*dist1
                range1(4) = range1(4) + 0.1*dist1
             endif

             call MPI_BARRIER(mpicommwd,ierr)

             !update the computational domain after drift
             call update2d1G_CompDom(Ageom,range1,range2,grid2d)
             call getmsize_CompDom(Ageom,msize)
             hx = msize(1)
             hy = msize(2)
             call getrange_CompDom(Ageom,range)
             xmin = range(1)
             xmax = range(2)
             ymin = range(3)
             ymax = range(4)
             call getlctabnm_CompDom(Ageom,table)

             call deposit2d(nptlcslice1,nptlcslice1,innx,nytot,hx,hy,xmin,&
                  ymin,ray1,tmp2d)
             !tmp2d = tmp2d/count1(i)
             call guardsum2d1G(tmp2d,rho1lc,innx,inny,nytot,nproccol,commcol)
             rho1lc = rho1lc/count1(i)

             !                   nptlcslice2,nptlc1,nptlc2,innx,inny,myidy
             call fieldsolver2d1G(rho1lc,innx,inny,nprocrow,nproccol,&
                  innyr,nytot,table,hx,hy,myidy,commcol,shift21)
             call guardexch2d1G(rho1lc,tmp2d,innx,inny,nytot,commcol)

             !store the potential from strong beam.
             if(flagwk.eq.1) then
                do j1 = 1, nytot
                   do i1=1, innx
                      phistrg(i1,j1,(i-1)*nslice1+j) = tmp2d(i1,j1)
                   enddo
                enddo
             endif
             !you do not need to copy back to the strong beam since
             !strong beam will not be modified by beam-beam force.
             !if(nslice1.gt.1) then
             !  !drift back the strong beam 1
             !  sign = -1
             !  call driftIP(ray1,sc,sign,nptlcslice1)
             !endif
             !do k = 1, nptlcslice1
             !  beamslice1(:,k,i) = ray1(:,k)
             !enddo
          else
             !extract the potential from strong beam 1.
             do j1 = 1, nytot 
                do i1=1, innx
                   tmp2d(i1,j1) = phistrg(i1,j1,(i-1)*nslice1+j)
                enddo
             enddo
          endif
          !using Gaussian model
          if( (flagwk.eq.4) .or. (stepid.eq.1) ) then
             nptlcslice1 = countlc1(i)
             do k = 1, nptlcslice1
                ray1(:,k) = beamslice1(:,k,i)
             enddo

             !drift strong beam 1 to the location of collision
             if(nslice1.gt.1) then
                call driftIP(ray1,sc1,sign,nptlcslice1)
             endif

             nptslice1 = count1(i) !round to nearest integer
             call findmomT_Utility(ray1,maxptclSlice1,nptlcslice1,nptslice1,&
                  center,sigma)

             if(flagwk.eq.3) then !store transverse rms size and centeroid
                sigx(i) = sigma(1)
                sigy(i) = sigma(2)
                centx(i) = center(1)
                centy(i) = center(2)
             endif
          endif

          !apply the kick from strong beam to weak beam.
          xmin2 = xmin + shift21(1)
          ymin2 = ymin + shift21(2)

          if(flagwk.eq.1 .or. flagwk.eq.2 ) then
             call scatter2d(nptlcslice2,nptlcslice2,innx,&
                  hx,hy,xmin2,ymin2,ray2,tmp2d,weight1(i),coef2,&
                  nytot)
          else if(flagwk.eq.3) then
             sigma(1) = sigx(i)
             sigma(2) = sigy(i) 
             center(1) = centx(i) 
             center(2) = centy(i) 
             call scatter2dgauss(nptlcslice2,nptlcslice2,&
                  ray2,weight1(i),coef2,sigma,center)
          else if(flagwk.eq.4) then
             call scatter2dgauss(nptlcslice2,nptlcslice2,&
                  ray2,weight1(i),coef2,sigma,center)
          else
          endif

          if(nslice2.gt.1) then
             !move weak beam back to its normal position.
             sign = -1
             call driftIP(ray2,sc,sign,nptlcslice2)
          endif

          do k = 1, nptlcslice2
             beamslice2(:,k,j) = ray2(:,k)
          enddo

       enddo
    enddo

    ipt = 0
    do i = 1, nslice2
       do j = 1, countlc2(i)
          ipt = ipt + 1
          Pts2(:,ipt) = beamslice2(:,j,i)
       enddo
    enddo

    !we do not need to update the total number of weak beam since
    !it does not affect the beam-beam force on weak beam.
    !call MPI_BARRIER(mpicommwd,ierr)
    if(nptlc2.ne.ipt) then
       nptlc2 = ipt
       !  call MPI_ALLREDUCE(nptlc2,nptot2,1,MPI_INTEGER,MPI_SUM,&
       !                     mpicommwd,ierr)
    endif

  end subroutine bbwkst



  !//find the strong-weak beam-beam kick in a 3D multi-slice model.
  subroutine bbmslsw_Beambeam(nprocrow,nproccol,innx,inny,&
       nslice1,nslice2,jendmax,nptlc1,nptlc2,myidx,myidy,Pts1,Pts2,&
       grid2d,Ageom,zmin1,zmax1,zmin2,zmax2,nptot1,nptot2,coef2,nytot,commrow,commcol,&
       nz1,nz2,maxPtclSlice1,maxPtclSlice2,shift21,curin1,curin2,flglum,lum3d)
    implicit none
    integer, intent(in) :: nprocrow,nproccol,innx,nslice1,&
         nslice2,myidx,myidy,nptot1,nptot2,nytot,commrow,commcol,nz1,nz2,&
         maxPtclSlice1,maxPtclSlice2,jendmax,flglum
    integer, intent(inout) :: inny,nptlc1,nptlc2
    double precision, pointer, dimension(:,:) :: Pts1,Pts2
    double precision, intent(out) :: lum3d
    type (Pgrid2d), intent(in) :: grid2d
    type (CompDom), intent(inout) :: Ageom
    double precision, intent(in) :: zmin1,zmin2,coef2,zmax1,zmax2,curin1,curin2
    double precision, dimension(2), intent(in) :: shift21
    double precision, dimension(2) :: shift12
    double precision, dimension(nslice2) :: weight2,zslice2
    double precision, dimension(nslice1) :: weight1,zslice1
    double precision, dimension(6) :: range
    double precision, dimension(4) :: range1,range2
    double precision, dimension(3) :: msize
    double precision, dimension(innx,nytot) :: rho,phiout
    double precision, dimension(innx,inny) :: rholc
    double precision, dimension(innx,inny,jendmax) :: rhotot
    double precision, dimension(6,maxPtclSlice1) :: ray1
    double precision, dimension(6,maxPtclSlice1,nslice1) :: beamslice1
    double precision, dimension(6,maxPtclSlice2) :: ray2
    double precision, dimension(6,maxPtclSlice2,nslice2) :: beamslice2
    double precision :: xmin,ymin,hx,hy
    double complex, allocatable, dimension (:,:) :: grn
    integer, dimension(nslice1) :: countlc1,count1
    integer, dimension(nslice2) :: countlc2,count2
    integer, dimension(2,0:nprocrow-1,0:nproccol-1) :: table
    integer :: i,j,nptlcslice1,nptlcslice2,sign,ii,jj,ipt,k,ierr,&
         i1,j1,jend,istart,jstart,myjend,itmp1,itmp2,rowend,flagbeam
    integer, dimension(0:nproccol-1) :: xpystable,pytable
    integer :: nsxy1,nsxy2,nxpylc2,nxtot,itb,nxlum,nylum,c1,c2,i2
    double precision :: xtmptmp,lumtmp,w1,w2

    !//shift12 - beam 1 with respect to beam 2
    !//shift21 - beam 2 with respect to beam 1
    shift12 = -shift21
    nxtot = innx
    xtmptmp = 1.0
    lum3d = 0.0
    nxlum = 128
    nylum = 128

    ! set up the storage for the Green function
    ! +1 is from the real to complex fft.
    nsxy1 = (nxtot+1)/nproccol
    nsxy2 = (nxtot+1) - nproccol*nsxy1
    do i = 0, nproccol-1
       if(i.le.(nsxy2-1)) then
          xpystable(i) = nsxy1+1
       else
          xpystable(i) = nsxy1
       endif
    enddo

    nxpylc2 = xpystable(myidy)

    allocate(grn(2*nytot,nxpylc2))

    !//here, we have used equal area weight
    !//copy particles from a bunch into slices. Each slice has
    !//roughly same number of particles.
    call slicer1G(Pts1,nptlc1,nslice1,maxPtclSlice1,nz1,zmin1,&
         zmax1,beamslice1,weight1,zslice1,nptot1,count1,countlc1)
    !get beam slice 2
    !here, we have used equal area weight
    call slicer1G(Pts2,nptlc2,nslice2,maxPtclSlice2,nz2,zmin2,&
         zmax2,beamslice2,weight2,zslice2,nptot2,count2,countlc2)

    call MPI_BARRIER(mpicommwd,ierr)
    !//loop through the total overlaping steps.
    do i = 1, nslice1+nslice2-1
       !//find the number of colliding slices 
       if(i.le.nslice1) then
          jend = min(i,nslice2) !//jend is the number of the overlaped slices.
          istart = i+1  !//istart is the starting slice + 1 in beam1
          jstart = 0    !//jstart is is the starting slice - 1 in beam2
       else
          jend = min((nslice2-i+nslice1),nslice1)
          istart = nslice1+1
          jstart = i-nslice1
       endif

       !//drift the collision slice particles to the collision point.
       sign = 1
       if(nslice1.gt.1) then
          flagbeam = 1
          call driftCP1Gfix(beamslice1,sign,nslice1,nslice2,flagbeam,&
               istart,jstart,jend,zslice1,zslice2,countlc1,maxPtclSlice1)
       endif
       if(nslice2.gt.1) then
          flagbeam = 2
          call driftCP1G(beamslice2,sign,nslice2,nslice1,flagbeam,&
               istart,jstart,jend,zslice2,zslice1,countlc2,maxPtclSlice2)
       endif

       !//assume the Green function is same for all colliding slices at each step.
       !//assume all colliding slices at each step are contained in the same box.
       !find the range xmin,xmax,ymin,ymax for all the colliding jend slices.
       range1(1) = beamslice1(1,1,istart-1)
       range1(2) = beamslice1(1,1,istart-1)
       range1(3) = beamslice1(3,1,istart-1)
       range1(4) = beamslice1(3,1,istart-1)

       do j = 1, jend
          !//find the slice id number in each beam.
          i1 = istart - j

          !//find the range of slice particle after drift.
          nptlcslice1 = countlc1(i1)
          do k = 1, nptlcslice1
             if(range1(1).gt.beamslice1(1,k,i1)) then
                range1(1) = beamslice1(1,k,i1)
             else if(range1(2).le.beamslice1(1,k,i1)) then
                range1(2) = beamslice1(1,k,i1)
             else
             endif
             if(range1(3).gt.beamslice1(3,k,i1)) then
                range1(3) = beamslice1(3,k,i1)
             else if(range1(4).le.beamslice1(3,k,i1)) then
                range1(4) = beamslice1(3,k,i1)
             else
             endif
          enddo
       enddo

       range2(1) = beamslice2(1,1,jstart+1)
       range2(2) = beamslice2(1,1,jstart+1)
       range2(3) = beamslice2(3,1,jstart+1)
       range2(4) = beamslice2(3,1,jstart+1)
       do j = 1, jend
          !//find the slice id number in each beam.
          i1 = jstart + j

          !//find the range of slice particle after drift.
          nptlcslice2 = countlc2(i1)
          do k = 1, nptlcslice2
             if(range2(1).gt.beamslice2(1,k,i1)) then
                range2(1) = beamslice2(1,k,i1)
             else if(range2(2).le.beamslice2(1,k,i1)) then
                range2(2) = beamslice2(1,k,i1)
             else
             endif
             if(range2(3).gt.beamslice2(3,k,i1)) then
                range2(3) = beamslice2(3,k,i1)
             else if(range2(4).le.beamslice2(3,k,i1)) then
                range2(4) = beamslice2(3,k,i1)
             else
             endif
          enddo
       enddo
       !//shift the range to local range, ie. the range with respect
       !//to the center of the beam.
       range2(1) = range2(1) - shift21(1)
       range2(2) = range2(2) - shift21(1)
       range2(3) = range2(3) - shift21(2)
       range2(4) = range2(4) - shift21(2)

       !//update the computational domain after drift
       call update2d1G_CompDom(Ageom,range1,range2,grid2d)

       call getmsize_CompDom(Ageom,msize)
       hx = msize(1)
       hy = msize(2)
       call getlctabnm_CompDom(Ageom,table)
       do itb = 0, nproccol-1
          pytable(itb) = table(2,0,itb)
       enddo
       !//find the Green function for the domain containing "jend"
       !//slice.
       !            call greenf2d1G(nxtot,nytot,nxpylc2,myidy,nproccol,xpystable,hx,hy,&
       !                            grn,shift21)
       !using parallel Green function
       call greenf2d1G(nxtot,nytot,nxpylc2,hx,hy,myidy,nproccol,commcol,&
            xpystable,grn,shift21)

       call getrange_CompDom(Ageom,range)
       xmin = range(1)
       ymin = range(3)

       !//find the transverse charge density of overlaped slices
       !//and store in rhotot.
       rhotot = 0.0
       do j = 1, jend
          i1 = istart - j 
          i2 = jstart + j

          nptlcslice1 = countlc1(i1)
          do k = 1, nptlcslice1
             ray1(:,k) = beamslice1(:,k,i1)
          enddo
          nptlcslice2 = countlc2(i2)
          do k = 1, nptlcslice2
             ray2(:,k) = beamslice2(:,k,i2)
          enddo
          w1 = weight1(i1)
          w2 = weight2(i2)
          c1 = count1(i1)
          c2 = count2(i2)

          !//calculate the 3d luminosity
          if(flglum.eq.1) then
             call luminosity1G3d_Output(ray1,nptlcslice1,ray2,nptlcslice2,&
                  nxlum,nylum,curin1,curin2,w1,w2,c1,c2,lumtmp)
             lum3d = lum3d + lumtmp
          endif

          !//deposition for each slice.
          call deposit2d(maxPtclSlice1,nptlcslice1,innx,nytot,hx,hy,xmin,&
               ymin,ray1,rho)

          !//collect the contribution from the other processors in the same column.
          !//This function and the function guardsum2drow can be replaced by
          !//a Allreduce in the whole domain with double sized container.
          call guardsum2d1G(rho,rholc,innx,inny,nytot,nproccol,commcol)

          do jj = 1, inny
             do ii = 1, innx
                rhotot(ii,jj,j) = rholc(ii,jj)*weight1(i1)/count1(i1)
             enddo
          enddo
       enddo

       !//collect the density onto the colliding slices
       myjend = (jend-1)/nprocrow + 1
       if(nslice1.gt.1) then
          call guardsum2drow(rhotot,innx,inny,myjend,jendmax,nprocrow,&
               myidx,commrow)
       endif

       !//distribute "jend" slice along "nprocrow" row processors.
       itmp1 = jend/nprocrow
       itmp2 = jend - nprocrow*itmp1
       if(myidx.lt.itmp2) then
          rowend = itmp1 + 1
       else
          rowend = itmp1
       endif

       !//solve the Poisson's equation in parallel for processors along row.
       !//transpose will be used along column to solve the Poisson equation. 
       do j = 1, rowend

          !//store the density on one slice
          do jj = 1, inny
             do ii = 1, innx
                rholc(ii,jj) =  rhotot(ii,jj,j)
             enddo
          enddo
          !//solve the potential for each slice of distributed "jend" slice.
          call fieldsolver2d1Gwk(rholc,innx,inny,nprocrow,nproccol,&
               nytot,nxpylc2,xpystable,pytable,hx,hy,myidy,commcol,grn)

          do jj = 1, inny
             do ii = 1, innx
                rhotot(ii,jj,j) =  rholc(ii,jj)
             enddo
          enddo

       enddo

       !//scatter potentials along row so that each processor contain 
       !//"jend" slices potential. This function can be replaced by
       !//Allreduce command along row.
       if(nslice1.gt.1) then
          call guardexch2drow(rhotot,innx,inny,jendmax,myjend,nprocrow,&
               myidx,commrow)
       endif

       !interpolate the beam-beam field from strong beam onto 
       !particles of each slice of weak beam.
       xmin = range(1) + shift21(1)
       ymin = range(3) + shift21(2)
       do j = 1, jend 
          i1 = jstart + j
          j1 = istart -j

          nptlcslice2 = countlc2(i1)
          do k = 1, nptlcslice2
             ray2(:,k) = beamslice2(:,k,i1)
          enddo

          do jj = 1, inny
             do ii = 1, innx
                rholc(ii,jj) =  rhotot(ii,jj,j)
             enddo
          enddo

          !//collect all the potential along column to a single proc. This
          !//function can be replaced by Allgather along column with doubled size.
          call guardexch2d1G(rholc,phiout,innx,inny,nytot,commcol)

          !//interpolate field of each slice to slice particles
          call scatter2d(maxPtclSlice2,nptlcslice2,innx,&
               hx,hy,xmin,ymin,ray2,phiout,xtmptmp,coef2,&
               nytot)

          do k = 1, nptlcslice2
             beamslice2(:,k,i1) = ray2(:,k)
          enddo
       enddo

       !//drift all particles back to their original longitudinal location
       sign = -1
       if(nslice1.gt.1) then
          flagbeam = 1
          call driftCP1Gfix(beamslice1,sign,nslice1,nslice2,flagbeam,&
               istart,jstart,jend,zslice1,zslice2,countlc1,maxPtclSlice1)
       endif
       if(nslice2.gt.1) then
          flagbeam = 2
          call driftCP1G(beamslice2,sign,nslice2,nslice1,flagbeam,&
               istart,jstart,jend,zslice2,zslice1,countlc2,maxPtclSlice2)
       endif
    enddo

    deallocate(grn)

    call MPI_BARRIER(mpicommwd,ierr)
    !//copy each slice particle to the beam bunch 
    ipt = 0
    do i = 1, nslice2
       do j = 1, countlc2(i)
          ipt = ipt + 1
          Pts2(:,ipt) = beamslice2(:,j,i)
       enddo
    enddo

  end subroutine bbmslsw_Beambeam



  !//find the strong-weak beam-beam kick in a 3D multi-slice model
  !//using a fixed domain Green function.
  subroutine bbmslswfix_Beambeam(nprocrow,nproccol,innx,inny,&
       nslice1,nslice2,jendmax,nptlc1,nptlc2,myidx,myidy,Pts1,&
       Pts2,grid2d,Ageom,zmin1,zmax1,zmin2,zmax2,nptot1,nptot2,&
       coef2,nytot,commrow,commcol,nz1,nz2,maxPtclSlice1,maxPtclSlice2,&
       nxpylc2,grn,xpystable,pytable,hx,hy,xmin,ymin,curin1,&
       curin2,flglum,lum3d)
    implicit none
    integer, intent(in) :: nprocrow,nproccol,innx,nslice1,&
         nslice2,myidx,myidy,nptot1,nptot2,nytot,commrow,commcol,&
         nz1,nz2,maxPtclSlice1,maxPtclSlice2,jendmax,nxpylc2,flglum
    integer, intent(inout) :: inny,nptlc1,nptlc2
    double precision, pointer, dimension(:,:) :: Pts1,Pts2
    double precision, intent(out) :: lum3d
    type (Pgrid2d), intent(in) :: grid2d
    type (CompDom), intent(inout) :: Ageom
    double precision, intent(in) :: zmin1,zmin2,coef2,zmax1,zmax2,&
         curin1,curin2
    double precision, dimension(nslice1) :: weight1,zslice1
    double precision, dimension(nslice2) :: weight2,zslice2
    double precision, dimension(innx,nytot) :: rho,phiout
    double precision, dimension(innx,inny) :: rholc
    double precision, dimension(innx,inny,jendmax) :: rhotot
    double precision, dimension(6,maxPtclSlice1) :: ray1
    double precision, dimension(6,maxPtclSlice2) :: ray2
    double precision, dimension(6,maxPtclSlice1,nslice1) :: beamslice1
    double precision, dimension(6,maxPtclSlice2,nslice2) :: beamslice2
    double precision :: xmin,ymin,hx,hy
    double complex, dimension (2*nytot,nxpylc2) :: grn
    integer, dimension(nslice1) :: countlc1,count1
    integer, dimension(nslice2) :: countlc2,count2
    integer :: i,j,nptlcslice1,nptlcslice2,sign,ii,jj,ipt,k,ierr,&
         i1,j1,jend,istart,jstart,myjend,itmp1,itmp2,rowend,flagbeam
    integer, dimension(0:nproccol-1) :: xpystable,pytable
    integer :: nxtot,nxlum,nylum,c1,c2,i2
    double precision :: xtmptmp,lumtmp,w1,w2

    nxtot = innx
    xtmptmp = 1.0
    lum3d = 0.0
    nxlum = 128
    nylum = 128

    !//here, we have used equal area weight
    !//copy particles from a bunch into slices. Each slice has
    !//roughly same number of particles.
    call slicer1G(Pts1,nptlc1,nslice1,maxPtclSlice1,nz1,zmin1,&
         zmax1,beamslice1,weight1,zslice1,nptot1,count1,countlc1)
    !get beam slice 2
    !here, we have used equal area weight
    call slicer1G(Pts2,nptlc2,nslice2,maxPtclSlice2,nz2,zmin2,&
         zmax2,beamslice2,weight2,zslice2,nptot2,count2,countlc2)

    call MPI_BARRIER(mpicommwd,ierr)

    !//loop through the total overlaping steps.
    do i = 1, nslice1+nslice2-1
       !//find the number of colliding slices 
       if(i.le.nslice1) then
          jend = min(i,nslice2)
          istart = i+1
          jstart = 0
       else
          jend = min((nslice2-i+nslice1),nslice1)
          istart = nslice1+1
          jstart = i-nslice1
       endif

       !//drift the collision slice particles to the collision point.
       sign = 1
       if(nslice1.gt.1) then
          flagbeam = 1
          call driftCP1Gfix(beamslice1,sign,nslice1,nslice2,flagbeam,&
               istart,jstart,jend,zslice1,zslice2,countlc1,maxPtclSlice1)
       endif
       if(nslice2.gt.1) then
          flagbeam = 2
          call driftCP1G(beamslice2,sign,nslice2,nslice1,flagbeam,&
               istart,jstart,jend,zslice2,zslice1,countlc2,maxPtclSlice2)
       endif

       !find the transverse charge density of all slices
       rhotot = 0.0
       do j = 1, jend
          i1 = istart - j
          i2 = jstart + j

          nptlcslice1 = countlc1(i1)
          do k = 1, nptlcslice1
             ray1(:,k) = beamslice1(:,k,i1)
          enddo
          nptlcslice2 = countlc2(i2)
          do k = 1, nptlcslice2
             ray2(:,k) = beamslice2(:,k,i2)
          enddo
          w1 = weight1(i1)
          w2 = weight2(i2)
          c1 = count1(i1)
          c2 = count2(i2)

          !//calculate the 3d luminosity
          if(flglum.eq.1) then
             call luminosity1G3d_Output(ray1,nptlcslice1,ray2,nptlcslice2,&
                  nxlum,nylum,curin1,curin2,w1,w2,c1,c2,lumtmp)
             lum3d = lum3d + lumtmp
          endif

          !//deposition for each slice.
          call deposit2d(maxPtclSlice1,nptlcslice1,innx,nytot,hx,hy,xmin,&
               ymin,ray1,rho)

          !//collect the contribution from the other processors in the same column.
          !//This function and the function guardsum2drow can be replaced by
          !//a Allreduce in the whole domain with double sized container.
          call guardsum2d1G(rho,rholc,innx,inny,nytot,nproccol,commcol)

          do jj = 1, inny
             do ii = 1, innx
                rhotot(ii,jj,j) = rholc(ii,jj)*weight1(i1)/count1(i1)
             enddo
          enddo
       enddo

       !collect the density onto the colliding slices
       myjend = (jend-1)/nprocrow + 1
       if(nslice1.gt.1) then
          call guardsum2drow(rhotot,innx,inny,myjend,jendmax,nprocrow,&
               myidx,commrow)
       endif

       itmp1 = jend/nprocrow
       itmp2 = jend - nprocrow*itmp1
       if(myidx.lt.itmp2) then
          rowend = itmp1 + 1
       else
          rowend = itmp1
       endif

       !solve the Poisson's equation in parallel for processors along row.
       !transpose will be used along column to solve the Poisson equation. 
       do j = 1, rowend

          do jj = 1, inny
             do ii = 1, innx
                rholc(ii,jj) =  rhotot(ii,jj,j)
             enddo
          enddo

          call fieldsolver2d1Gwk(rholc,innx,inny,nprocrow,nproccol,&
               nytot,nxpylc2,xpystable,pytable,hx,hy,myidy,commcol,grn)

          do jj = 1, inny
             do ii = 1, innx
                rhotot(ii,jj,j) =  rholc(ii,jj)
             enddo
          enddo

       enddo

       !scatter potentials along row
       if(nslice1.gt.1) then
          call guardexch2drow(rhotot,innx,inny,jendmax,myjend,nprocrow,&
               myidx,commrow)
       endif

       !interpolate the beam-beam field onto particles on each slice.
       do j = 1, jend 
          i1 = jstart + j
          j1 = istart -j

          nptlcslice2 = countlc2(i1)
          do k = 1, nptlcslice2
             ray2(:,k) = beamslice2(:,k,i1)
          enddo

          do jj = 1, inny
             do ii = 1, innx
                rholc(ii,jj) =  rhotot(ii,jj,j)
             enddo
          enddo

          !//collect all the potential along column to a single proc. This
          !//function can be replaced by Allgather along column with doubled size.
          call guardexch2d1G(rholc,phiout,innx,inny,nytot,commcol)

          !//interpolate field of each slice to slice particles
          call scatter2d(maxPtclSlice2,nptlcslice2,innx,&
               hx,hy,xmin,ymin,ray2,phiout,xtmptmp,coef2,&
               nytot)

          do k = 1, nptlcslice2
             beamslice2(:,k,i1) = ray2(:,k)
          enddo

       enddo

       !drift all particles back to their original longitudinal location
       sign = -1
       if(nslice1.gt.1) then
          flagbeam = 1
          call driftCP1Gfix(beamslice1,sign,nslice1,nslice2,flagbeam,&
               istart,jstart,jend,zslice1,zslice2,countlc1,maxPtclSlice1)
       endif
       if(nslice2.gt.1) then
          flagbeam = 2
          call driftCP1G(beamslice2,sign,nslice2,nslice1,flagbeam,&
               istart,jstart,jend,zslice2,zslice1,countlc2,maxPtclSlice2)
       endif
    enddo

    call MPI_BARRIER(mpicommwd,ierr)
    !determine the new local particle number on each PE.
    !copy each slice particle to the beam bunch 
    ipt = 0
    do i = 1, nslice2
       do j = 1, countlc2(i)
          ipt = ipt + 1
          Pts2(:,ipt) = beamslice2(:,j,i)
       enddo
    enddo

  end subroutine bbmslswfix_Beambeam



  !//find the strong-weak beam-beam kick in a 3D multi-slice model
  !//using Gaussian approxmiation of the strong beam.
  subroutine bbmslswgauss_Beambeam(nprocrow,nproccol,innx,inny,&
       nslice1,nslice2,jendmax,nptlc1,nptlc2,myidx,myidy,Pts1,Pts2,&
       grid2d,Ageom,zmin1,zmax1,zmin2,zmax2,nptot1,nptot2,coef2,nytot,commrow,commcol,&
       nz1,nz2,maxPtclSlice1,maxPtclSlice2,curin1,curin2,flglum,lum3d)
    implicit none
    integer, intent(in) :: nprocrow,nproccol,innx,nslice1,&
         nslice2,myidx,myidy,nptot1,nptot2,nytot,commrow,commcol,nz1,nz2,&
         maxPtclSlice1,maxPtclSlice2,jendmax,flglum
    integer, intent(inout) :: inny,nptlc1,nptlc2
    double precision, pointer, dimension(:,:) :: Pts1,Pts2
    double precision, intent(out) :: lum3d
    type (Pgrid2d), intent(in) :: grid2d
    type (CompDom), intent(inout) :: Ageom
    double precision, intent(in) :: zmin1,zmin2,coef2,zmax1,zmax2,curin1,curin2
    double precision, dimension(nslice2) :: weight2,zslice2
    double precision, dimension(nslice1) :: weight1,zslice1
    double precision, dimension(6,maxPtclSlice1) :: ray1
    double precision, dimension(6,maxPtclSlice1,nslice1) :: beamslice1
    double precision, dimension(6,maxPtclSlice2) :: ray2
    double precision, dimension(6,maxPtclSlice2,nslice2) :: beamslice2
    double precision, dimension(2) :: center,sigma
    integer, dimension(nslice1) :: countlc1,count1
    integer, dimension(nslice2) :: countlc2,count2
    integer :: i,j,nptlcslice1,nptlcslice2,sign,ipt,k,ierr,&
         i1,jend,istart,jstart,flagbeam
    integer :: nxtot,nxlum,nylum,c1,c2,i2
    double precision :: xtmptmp,lumtmp,w1,w2

    !//shift12 - beam 1 with respect to beam 2
    !//shift21 - beam 2 with respect to beam 1
    nxtot = innx
    xtmptmp = 1.0
    lum3d = 0.0
    nxlum = 128
    nylum = 128

    !//here, we have used equal area weight
    !//copy particles from a bunch into slices. Each slice has
    !//roughly same number of particles.
    call slicer1G(Pts1,nptlc1,nslice1,maxPtclSlice1,nz1,zmin1,&
         zmax1,beamslice1,weight1,zslice1,nptot1,count1,countlc1)
    !get beam slice 2
    !here, we have used equal area weight
    call slicer1G(Pts2,nptlc2,nslice2,maxPtclSlice2,nz2,zmin2,&
         zmax2,beamslice2,weight2,zslice2,nptot2,count2,countlc2)

    call MPI_BARRIER(mpicommwd,ierr)
    !//loop through the total overlaping steps.
    do i = 1, nslice1+nslice2-1
       !//find the number of colliding slices 
       if(i.le.nslice1) then
          jend = min(i,nslice2) !//jend is the number of the overlaped slices.
          istart = i+1  !//istart is the starting slice + 1 in beam1
          jstart = 0    !//jstart is is the starting slice - 1 in beam2
       else
          jend = min((nslice2-i+nslice1),nslice1)
          istart = nslice1+1
          jstart = i-nslice1
       endif

       !//drift the collision slice particles to the collision point.
       sign = 1
       if(nslice1.gt.1) then
          flagbeam = 1
          call driftCP1Gfix(beamslice1,sign,nslice1,nslice2,flagbeam,&
               istart,jstart,jend,zslice1,zslice2,countlc1,maxPtclSlice1)
       endif
       if(nslice2.gt.1) then
          flagbeam = 2
          call driftCP1G(beamslice2,sign,nslice2,nslice1,flagbeam,&
               istart,jstart,jend,zslice2,zslice1,countlc2,maxPtclSlice2)
       endif

       !//loop through the overlaping slices
       do j = 1, jend 
          i1 = istart - j
          i2 = jstart + j

          nptlcslice1 = countlc1(i1)
          do k = 1, nptlcslice1
             ray1(:,k) = beamslice1(:,k,i1)
          enddo
          nptlcslice2 = countlc2(i2)
          do k = 1, nptlcslice2
             ray2(:,k) = beamslice2(:,k,i2)
          enddo
          w1 = weight1(i1)
          w2 = weight2(i2)
          c1 = count1(i1)
          c2 = count2(i2)

          !//calculate the 3d luminosity
          if(flglum.eq.1) then
             call luminosity1G3d_Output(ray1,nptlcslice1,ray2,nptlcslice2,&
                  nxlum,nylum,curin1,curin2,w1,w2,c1,c2,lumtmp)
             lum3d = lum3d + lumtmp
          endif

          !find the centoid and rms size of the strong beam slice
          call findmomT_Utility(ray1,maxPtclSlice1,nptlcslice1,c1,&
               center,sigma)

          !calculate the field on each particle assuming the Gaussian
          !distribution of the strong beam.
          call scatter2dgauss(maxPtclSlice2,nptlcslice2,&
               ray2,w1,coef2,sigma,center)

          do k = 1, nptlcslice2
             beamslice2(:,k,i2) = ray2(:,k)
          enddo
       enddo

       !//drift all particles back to their original longitudinal location
       sign = -1
       if(nslice1.gt.1) then
          flagbeam = 1
          call driftCP1Gfix(beamslice1,sign,nslice1,nslice2,flagbeam,&
               istart,jstart,jend,zslice1,zslice2,countlc1,maxPtclSlice1)
       endif
       if(nslice2.gt.1) then
          flagbeam = 2
          call driftCP1G(beamslice2,sign,nslice2,nslice1,flagbeam,&
               istart,jstart,jend,zslice2,zslice1,countlc2,maxPtclSlice2)
       endif
    enddo

    call MPI_BARRIER(mpicommwd,ierr)
    !//copy each slice particle to the beam bunch 
    ipt = 0
    do i = 1, nslice2
       do j = 1, countlc2(i)
          ipt = ipt + 1
          Pts2(:,ipt) = beamslice2(:,j,i)
       enddo
    enddo

  end subroutine bbmslswgauss_Beambeam


  !cut a bunch into slices using 1 group PEs.
  !//zslice is the z location of itself.
  subroutine slicer1G(Pts1in,nptlc,nslice,maxsptcl,nz,zmin,&
       zmax,beamslice,weight,zslice,nptot,count,countlc)
    implicit none
    include 'mpif.h'
    integer, parameter :: Ndiv = 10
    integer, intent(in) :: nptlc,nslice,maxsptcl,nptot,nz
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, intent(in) :: zmin,zmax
    double precision, dimension(6,maxsptcl,nslice),&
         intent(inout) :: beamslice
    double precision, dimension(nslice),intent(out) :: weight,zslice
    integer, dimension(nslice), intent(out) :: count,countlc
    double precision, dimension(nslice)  :: fcount,fcountlc
    double precision, dimension(0:nslice) :: rangez
    double precision, dimension(nz) :: rhozlc,rhoz
    double precision, dimension(nz*Ndiv) :: rhoznew
    integer :: i,ierr,ii,iz,iz1,ip
    double precision :: ab,sumint,hnew1,hnew2,interval0,hz,hzi
    !integer status(MPI_STATUS_SIZE)

    if(nslice.gt.1) then
       !//find the charge density along z
       hz = (zmax-zmin)/(nz-1)
       hzi = 1.0/hz
       rhozlc = 0.0
       do i = 1, nptlc
          iz=int((Pts1in(5,i)-zmin)*hzi + 1)
          ab=(zmin-Pts1in(5,i)+iz*hz)*hzi
          iz1 = iz + 1
          rhozlc(iz) = rhozlc(iz) + ab
          rhozlc(iz1) = rhozlc(iz1) + (1.0-ab)
       enddo
       rhozlc = rhozlc/nptot

       call MPI_ALLREDUCE(rhozlc,rhoz,nz,MPI_DOUBLE_PRECISION,&
            MPI_SUM,mpicommwd,ierr)

       !//find totol density, should be 1
       sumint = 0.0
       do i = 1, nz
          sumint = sumint + rhoz(i)
       enddo

       !//subdivided the charge density for finer resolution.
       do i = 1, nz*Ndiv
          ii = (i-1)/Ndiv + 1
          rhoznew(i) = rhoz(ii)/Ndiv
       enddo
       hnew1 = hz/(Ndiv-1)/2
       hnew2 = hz/Ndiv

       interval0 = sumint/nslice
       sumint = 0.0
       !//find the range of each slice
       do i = 1, Ndiv*nz
          sumint = sumint + rhoznew(i)
          ip = int(sumint/interval0)
          if(ip.ne.nslice) then
             if(i.le.Ndiv) then
                rangez(ip+1) = zmin + (i-1)*hnew1 + 0.5*hnew1
             else if((i.gt.Ndiv).and.(i.le.Ndiv*nz-Ndiv)) then
                rangez(ip+1) = zmin +0.5*hz+(i-Ndiv)*hnew2 + 0.5*hnew2
             else
                rangez(ip+1) = zmin + (nz-1)*hz - 0.5*hz + &
                     (i+Ndiv-Ndiv*nz)*hz/Ndiv
             endif
          endif
       enddo
       rangez(nslice) = zmax
       rangez(0) = zmin


       !//find the # of particle per slice and put the particles into
       !//single slice. Here, zmax corresponds to slice # 1.
       do i = 1, nslice
          countlc(i) = 0
       enddo
       do i = 1,nptlc
          if(Pts1in(5,i).gt.rangez(nslice-1)) then
             !countlc(nslice) = countlc(nslice) + 1
             countlc(1) = countlc(1) + 1
             if(countlc(1).gt.maxsptcl) then
                print*,"overflow local slice limit in nslice: ",countlc(1)
                stop
             endif
             beamslice(:,countlc(1),1) = Pts1in(:,i)
          else if(Pts1in(5,i).le.rangez(1)) then
             !countlc(1) = countlc(1) + 1
             countlc(nslice) = countlc(nslice) + 1
             if(countlc(nslice).gt.maxsptcl) then
                print*,"overflow local slice limit in nslice: ",countlc(nslice)
                stop
             endif
             beamslice(:,countlc(nslice),nslice) = Pts1in(:,i)
          else
             do ii = 2,nslice-1
                if((Pts1in(5,i).gt.rangez(ii-1)).and.&
                     (Pts1in(5,i).le.rangez(ii)) ) then
                   countlc(nslice-ii+1) = countlc(nslice-ii+1) + 1
                   if(countlc(nslice-ii+1).gt.maxsptcl) then
                      print*,"overflow local slice limit in ii: ",ii,&
                           countlc(nslice-ii+1)
                      stop
                   endif
                   beamslice(:,countlc(nslice-ii+1),nslice-ii+1)=Pts1in(:,i)
                endif
             enddo
          endif
       enddo

       fcountlc = 0
       do i = 1, nslice
          fcountlc(i) = countlc(i)
       enddo
       call MPI_ALLREDUCE(fcountlc,fcount,nslice,MPI_DOUBLE_PRECISION,&
            MPI_SUM,mpicommwd,ierr)

       do i = 1, nslice
          count(i) = int(fcount(i) + 0.0001) !//0.0001 for roundoff
          weight(i) = fcount(i)/nptot
       enddo

       do i = 1, nslice
          !zslice(i) = (rangez(i-1)+rangez(i))/2 !//this is a terrible bug.
          zslice(nslice-i+1) = (rangez(i-1)+rangez(i))/2
       enddo

    else !//single slice case
       do i = 1, nptlc
          beamslice(:,i,nslice) = Pts1in(:,i)
       enddo
       weight(nslice) = 1.0
       rangez(nslice) = zmax
       rangez(0) = zmin
       zslice(1) = 0.0
       countlc = nptlc
       count(1) = nptot
    endif

  end subroutine slicer1G



  !//drift the whole slice particles same distance to the collision point.
  subroutine driftCP1Gfix(ray,sign,nslice,nslice2,flagbeam,&
       istart,jstart,jcross,zslice,zslice2,countlc,maxnpt)
    integer, intent(in) :: sign,nslice,flagbeam,istart,jstart,&
         jcross,nslice2,maxnpt
    integer, dimension(nslice) :: countlc
    double precision, dimension(nslice) :: zslice
    double precision, dimension(nslice2) :: zslice2
    double precision :: distance
    double precision, dimension(6,maxnpt,nslice), intent(inout) :: ray
    integer :: i,j,i1,j1

    do j = 1, jcross !//loop through the number of overlaped slices.
       if(flagbeam.eq.1) then
          i1 = istart - j
          j1 = jstart + j
       else if(flagbeam.eq.2) then
          i1 = jstart + j
          j1 = istart - j
       endif
       !//the drift distance from each beam has the same value but
       !//opposite sign. This means that if one beam slice particle moves
       !//forward, the other beam slice particle is going to move
       !//backward.
       distance = (zslice(i1)-zslice2(j1))/2
       do i = 1, countlc(i1)
          ray(1,i,i1) = ray(1,i,i1) + distance*sign*ray(2,i,i1)
          ray(3,i,i1) = ray(3,i,i1) + distance*sign*ray(4,i,i1)
       enddo
    enddo

  end subroutine driftCP1Gfix



  !//drift slice particles to the collision field interpolation point.
  subroutine driftCP1G(ray,sign,nslice,nslice2,flagbeam,&
       istart,jstart,jcross,zslice,zslice2,countlc,maxnpt)
    integer, intent(in) :: sign,nslice,flagbeam,istart,jstart,&
         jcross,nslice2,maxnpt
    integer, dimension(nslice) :: countlc
    double precision, dimension(nslice) :: zslice
    double precision, dimension(nslice2) :: zslice2
    double precision :: distance
    double precision, dimension(6,maxnpt,nslice), intent(inout) :: ray
    integer :: i,j,i1,j1

    do j = 1, jcross !//loop through the number of overlaped slices.
       if(flagbeam.eq.1) then
          i1 = istart - j
          j1 = jstart + j
       else if(flagbeam.eq.2) then
          i1 = jstart + j
          j1 = istart - j
       endif
       !//the drift distance from each beam has the same value but
       !//opposite sign. This means that if one beam slice particle moves
       !//forward, the other beam slice particle is going to move
       !//backward.
       !distance = (zslice(i1)-zslice2(j1))/2
       do i = 1, countlc(i1)
          distance = (ray(5,i,i1)-zslice2(j1))/2
          ray(1,i,i1) = ray(1,i,i1) + distance*sign*ray(2,i,i1)
          ray(3,i,i1) = ray(3,i,i1) + distance*sign*ray(4,i,i1)
       enddo
    enddo

  end subroutine driftCP1G


  !--------------------------------------------------------------------------
  !//find the strong-strong beam-beam kick in a 3D multi-slice model.
  !//Here, the beam-beam kick is based on the linear interpolation along
  !//longitudinal direction. This will slow down the code significantly.
  !//some test runs using KEKB parameters does not show that using
  !//the linear interpolation can reduce the requirement of the # of slice.
  subroutine bbmSlcLinear_Beambeam(nprocrow,npyhalf,innx,inny,&
       nslice1,nslice2,nslice,nsliceO,jendmax,nptlc,myidx,myidy,Pts,&
       grid2d,Ageom,zmin,zmax,nptot,coef,nytot,commrow,commcol,nz,&
       maxPtclSlice,shift21,curin,flglum,lum3d)
    implicit none
    integer, intent(in) :: nprocrow,npyhalf,innx,nslice1,&
         nslice2,myidx,myidy,nptot,nytot,commrow,commcol,nz,&
         maxPtclSlice,nslice,nsliceO,jendmax,flglum
    integer, intent(inout) :: inny,nptlc
    double precision, pointer, dimension(:,:) :: Pts
    double precision, intent(out) :: lum3d
    type (Pgrid2d), intent(in) :: grid2d
    type (CompDom), intent(inout) :: Ageom
    double precision, intent(in) :: zmin,coef,zmax,curin
    double precision, dimension(2), intent(in) :: shift21
    double precision, dimension(2) :: shift12,shift
    double precision, dimension(nsliceO) :: zslice2,zhalf2,ztmp2
    double precision, dimension(nslice) :: weight,zslice,zhalf
    double precision, dimension(6) :: range
    double precision, dimension(4) :: range1
    double precision, dimension(3) :: msize
    double precision, dimension(innx,nytot) :: rho,phiout,phiout2
    double precision, dimension(innx,inny) :: rholc
    double precision, dimension(innx,inny,jendmax) :: rhotot,rhotot2
    double precision, dimension(6,maxPtclSlice) :: ray
    double precision, dimension(6,maxPtclSlice,nslice) :: beamslice
    double precision :: xmin,ymin,hx,hy
    double complex, allocatable, dimension (:,:) :: grn
    integer, dimension(nslice) :: countlc,count
    integer, dimension(2,0:nprocrow-1,0:npyhalf-1) :: table
    integer :: i,j,nptlcslice,sign,ii,jj,ipt,k,&
         i1,j1,jend,istart,jstart,myjend,itmp1,itmp2,rowend !,ierr
    integer, dimension(0:npyhalf-1) :: xpystable,pytable
    integer :: npbc,nsxy1,nsxy2,nxpylc2,nxtot,itb,nxlum,nylum
    double precision :: xtmptmp,lumtmp,z2,zthk,sbk
    double precision :: distance1,distance2,distance3,tmpx1,&
         tmpx2,tmpx3,tmpy1,tmpy2,tmpy3

    !//shift12 - beam 1 with respect to beam 2
    !//shift21 - beam 2 with respect to beam 1
    shift12 = -shift21
    nxtot = innx
    xtmptmp = 1.0
    lum3d = 0.0
    nxlum = 128
    nylum = 128

    !//set up the parameters for the lower group and upper group processors.
    if(myidy.lt.npyhalf) then
       shift = shift21
    else
       shift = shift12
    endif
    if(myidy.lt.npyhalf) then
       npbc = 0
    else
       npbc = npyhalf
    endif

    ! +1 is from the real to complex fft.
    nsxy1 = (nxtot+1)/npyhalf
    nsxy2 = (nxtot+1) - npyhalf*nsxy1
    do i = 0, npyhalf-1
       if(i.le.(nsxy2-1)) then
          xpystable(i) = nsxy1+1
       else
          xpystable(i) = nsxy1
       endif
    enddo

    nxpylc2 = xpystable(myidy-npbc)

    allocate(grn(2*nytot,nxpylc2))

    !//here, we have used equal area weight
    !//copy particles from a bunch into slices. Each slice has
    !//roughly same number of particles.
    call slicer2z(Pts,nptlc,nslice,nsliceO,maxPtclSlice,nz,zmin,&
         zmax,beamslice,weight,zslice,zslice2,zhalf,zhalf2,nptot,count,&
         countlc,myidy,npyhalf,commcol,commrow)

    !        call MPI_BARRIER(mpicommwd,ierr)

    !//loop through the total overlaping steps.
    do i = 1, nslice1+nslice2-1
       !//find the number of colliding slices 
       if(i.le.nslice1) then
          jend = min(i,nslice2) !//jend is the number of the overlaped slices.
          istart = i+1  !//istart is the starting slice + 1 in beam1
          jstart = 0    !//jstart is is the starting slice - 1 in beam2
       else
          jend = min((nslice2-i+nslice1),nslice1)
          istart = nslice1+1
          jstart = i-nslice1
       endif

       !//find the range of the colliding slices
       !find the range xmin,xmax,ymin,ymax for all the colliding slices,jend.
       if(myidy.lt.npyhalf) then
          distance3 = (beamslice(5,1,istart-1)-zslice2(jstart+1))/2
          !distance3 = (zslice(istart-1)-zslice2(jstart+1))/2
          tmpx1 = beamslice(1,1,istart-1) + distance3*beamslice(2,1,istart-1)
          tmpy1 = beamslice(3,1,istart-1) + distance3*beamslice(4,1,istart-1)
          range1(1) = tmpx1
          range1(2) = tmpx1
          range1(3) = tmpy1
          range1(4) = tmpy1
       else
          distance3 = (beamslice(5,1,jstart+1)-zslice2(istart-1))/2
          !distance3 = (zslice(jstart+1)-zslice2(istart-1))/2
          tmpx1 = beamslice(1,1,jstart+1) + distance3*beamslice(2,1,jstart+1)
          tmpy1 = beamslice(3,1,jstart+1) + distance3*beamslice(4,1,jstart+1)
          range1(1) = tmpx1
          range1(2) = tmpx1
          range1(3) = tmpy1
          range1(4) = tmpy1
       endif
       do j = 1, jend !//loop through the number of overlaped slices.
          if(myidy.lt.npyhalf) then
             i1 = istart - j
             j1 = jstart + j
          else
             i1 = jstart + j
             j1 = istart - j
          endif
          distance1 = (zslice(i1)-(zslice2(j1)-zhalf2(j1)))/2
          distance2 = (zslice(i1)-(zslice2(j1)+zhalf2(j1)))/2
          do ipt = 1, countlc(i1)
             distance3 = (beamslice(5,ipt,i1)-zslice2(j1))/2
             tmpx1 = beamslice(1,ipt,i1) + distance1*beamslice(2,ipt,i1)
             tmpx2 = beamslice(1,ipt,i1) + distance2*beamslice(2,ipt,i1)
             tmpx3 = beamslice(1,ipt,i1) + distance3*beamslice(2,ipt,i1)
             tmpy1 = beamslice(3,ipt,i1) + distance1*beamslice(4,ipt,i1)
             tmpy2 = beamslice(3,ipt,i1) + distance2*beamslice(4,ipt,i1)
             tmpy3 = beamslice(3,ipt,i1) + distance3*beamslice(4,ipt,i1)
             if(range1(1).gt.tmpx1) range1(1) = tmpx1
             if(range1(1).gt.tmpx2) range1(1) = tmpx2
             if(range1(1).gt.tmpx3) range1(1) = tmpx3
             if(range1(2).lt.tmpx1) range1(2) = tmpx1
             if(range1(2).lt.tmpx2) range1(2) = tmpx2
             if(range1(2).lt.tmpx3) range1(2) = tmpx3
             if(range1(3).gt.tmpy1) range1(3) = tmpy1
             if(range1(3).gt.tmpy2) range1(3) = tmpy2
             if(range1(3).gt.tmpy3) range1(3) = tmpy3
             if(range1(4).lt.tmpy1) range1(4) = tmpy1
             if(range1(4).lt.tmpy2) range1(4) = tmpy2
             if(range1(4).lt.tmpy3) range1(4) = tmpy3
          enddo
       enddo
       !//shift the range to local range, ie. the range with respect
       !//to the center of the beam.
       if(myidy.ge.npyhalf) then
          range1(1) = range1(1) - shift21(1)
          range1(2) = range1(2) - shift21(1)
          range1(3) = range1(3) - shift21(2)
          range1(4) = range1(4) - shift21(2)
       endif
       !//update the computational domain after drift
       call update2d_CompDom(Ageom,range1,grid2d)
       call getmsize_CompDom(Ageom,msize)
       hx = msize(1)
       hy = msize(2)
       call getlctabnm_CompDom(Ageom,table)
       do itb = 0, npyhalf-1
          pytable(itb) = table(2,0,itb)
       enddo
       !//find the Green function for the domain containing "jend"
       !//slice.
       !//assume the Green function is same for all colliding slices at each step.
       !//assume all colliding slices at each step are contained in the same box.
       call greenf2d(nxtot,nytot,nxpylc2,hx,hy,myidy,npyhalf,commcol,&
            xpystable,grn,shift)
       call getrange_CompDom(Ageom,range)
       if(myidy.lt.npyhalf) then
          xmin = range(1)
          ymin = range(3)
       else
          xmin = range(1) + shift21(1)
          ymin = range(3) + shift21(2)
       endif

       !//drift the collision slice particles to the back of the collision point.
       ztmp2 = zslice2 - zhalf2
       sign = 1
       if(nslice.gt.1) then
          call driftCPfix(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
               istart,jstart,jend,zslice,ztmp2,countlc,maxPtclSlice)
       endif

       !//find the transverse charge density of overlaped slices
       rhotot = 0.0
       do j = 1, jend
          if(myidy.lt.npyhalf) then
             i1 = istart - j 
          else
             i1 = jstart + j
          endif

          nptlcslice = countlc(i1)
          do k = 1, nptlcslice
             ray(:,k) = beamslice(:,k,i1)
          enddo

          !//deposition for each slice.
          call deposit2d(maxPtclSlice,nptlcslice,innx,nytot,hx,hy,xmin,&
               ymin,ray,rho)

          !//collect the contribution from the other processors in the same column.
          !//This function and the function guardsum2drow can be replaced by
          !//a Allreduce in the whole domain with double sized container.
          call guardsum2d(rho,rholc,innx,inny,nytot,npyhalf,myidy,commcol)

          do jj = 1, inny
             do ii = 1, innx
                rhotot(ii,jj,j) = rholc(ii,jj)*weight(i1)/count(i1)
             enddo
          enddo
       enddo

       !//collect the density onto the colliding slices
       myjend = (jend-1)/nprocrow + 1
       if(nslice.gt.1) then
          call guardsum2drow(rhotot,innx,inny,myjend,jendmax,nprocrow,&
               myidx,commrow)
       endif

       !//distribute "jend" slice along "nprocrow" row processors.
       itmp1 = jend/nprocrow
       itmp2 = jend - nprocrow*itmp1
       if(myidx.lt.itmp2) then
          rowend = itmp1 + 1
       else
          rowend = itmp1
       endif

       !//solve the Poisson's equation in parallel for processors along row.
       !//transpose will be used along column to solve the Poisson equation. 
       do j = 1, rowend
          do jj = 1, inny
             do ii = 1, innx
                rholc(ii,jj) =  rhotot(ii,jj,j)
             enddo
          enddo

          !//solve the potential for each slice of distributed "jend" slice.
          !call fieldsolver2d(rholc,innx,inny,nprocrow,npyhalf,&
          !nytot,nxpylc2,xpystable,pytable,hx,hy,myidy,commcol,shift,grn)
          call fieldsolver2d(rholc,innx,inny,nprocrow,npyhalf,&
               nytot,nxpylc2,xpystable,pytable,hx,hy,myidy,commcol,grn)

          do jj = 1, inny
             do ii = 1, innx
                rhotot(ii,jj,j) =  rholc(ii,jj)
             enddo
          enddo
       enddo

       !drift back the fixed distance
       sign = -1
       if(nslice.gt.1) then
          call driftCPfix(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
               istart,jstart,jend,zslice,ztmp2,countlc,maxPtclSlice)
       endif

       !//calculate the 3d luminosity, luminosity is calculated at i-j.
       if(flglum.eq.1) then
          sign = 1
          if(nslice.gt.1) then
             call driftCPfix(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
                  istart,jstart,jend,zslice,zslice2,countlc,maxPtclSlice)
          endif
          do j = 1, jend
             if(myidy.lt.npyhalf) then
                i1 = istart - j
             else
                i1 = jstart + j
             endif

             nptlcslice = countlc(i1)
             do k = 1, nptlcslice
                ray(:,k) = beamslice(:,k,i1)
             enddo

             call luminosity2G3d_Output(ray,nptlcslice,nxlum,nylum,&
                  curin,weight(i1),count(i1),myidy,npyhalf,lumtmp)
             lum3d = lum3d + lumtmp
          enddo
          sign = -1
          if(nslice.gt.1) then
             call driftCPfix(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
                  istart,jstart,jend,zslice,zslice2,countlc,maxPtclSlice)
          endif
       endif

       !drift to the front of the collision point
       ztmp2 = zslice2 + zhalf2
       sign = 1
       if(nslice.gt.1) then
          call driftCPfix(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
               istart,jstart,jend,zslice,ztmp2,countlc,maxPtclSlice)
       endif
       !//find the transverse charge density of overlaped slices
       rhotot2 = 0.0
       do j = 1, jend
          if(myidy.lt.npyhalf) then
             i1 = istart - j
          else
             i1 = jstart + j
          endif

          nptlcslice = countlc(i1)
          do k = 1, nptlcslice
             ray(:,k) = beamslice(:,k,i1)
          enddo

          !//deposition for each slice.
          call deposit2d(maxPtclSlice,nptlcslice,innx,nytot,hx,hy,xmin,&
               ymin,ray,rho)

          call guardsum2d(rho,rholc,innx,inny,nytot,npyhalf,myidy,commcol)

          do jj = 1, inny
             do ii = 1, innx
                rhotot2(ii,jj,j) = rholc(ii,jj)*weight(i1)/count(i1)
             enddo
          enddo
       enddo

       !//collect the density onto the colliding slices
       myjend = (jend-1)/nprocrow + 1
       if(nslice.gt.1) then
          call guardsum2drow(rhotot2,innx,inny,myjend,jendmax,nprocrow,&
               myidx,commrow)
       endif

       !//distribute "jend" slice along "nprocrow" row processors.
       itmp1 = jend/nprocrow
       itmp2 = jend - nprocrow*itmp1
       if(myidx.lt.itmp2) then
          rowend = itmp1 + 1
       else
          rowend = itmp1
       endif

       !//solve the Poisson's equation in parallel for processors along row.
       !//transpose will be used along column to solve the Poisson equation.
       do j = 1, rowend
          do jj = 1, inny
             do ii = 1, innx
                rholc(ii,jj) =  rhotot2(ii,jj,j)
             enddo
          enddo
          !//solve the potential for each slice of distributed "jend" slice.
          call fieldsolver2d(rholc,innx,inny,nprocrow,npyhalf,&
               nytot,nxpylc2,xpystable,pytable,hx,hy,myidy,commcol,grn)

          do jj = 1, inny
             do ii = 1, innx
                rhotot2(ii,jj,j) =  rholc(ii,jj)
             enddo
          enddo
       enddo

       !drift back the fixed distance
       sign = -1
       if(nslice.gt.1) then
          call driftCPfix(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
               istart,jstart,jend,zslice,ztmp2,countlc,maxPtclSlice)
       endif
       !drift to the field interpolation point
       sign = 1
       if(nslice.gt.1) then
          call driftCP(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
               istart,jstart,jend,zslice,zslice2,countlc,maxPtclSlice)
       endif


       !//scatter potentials along row so that each processor contain 
       !//"jend" slices potential. This function can be replaced by
       !//Allreduce command along row.
       if(nslice.gt.1) then
          call guardexch2drow(rhotot,innx,inny,jendmax,myjend,nprocrow,&
               myidx,commrow)
       endif
       if(nslice.gt.1) then
          call guardexch2drow(rhotot2,innx,inny,jendmax,myjend,nprocrow,&
               myidx,commrow)
       endif


       !interpolate the beam-beam field onto particles on each slice.
       do j = 1, jend 
          if(myidy.lt.npyhalf) then
             i1 = istart - j 
             j1 = jstart + j
          else
             i1 = jstart + j
             j1 = istart -j
          endif

          nptlcslice = countlc(i1)
          do k = 1, nptlcslice
             ray(:,k) = beamslice(:,k,i1)
          enddo

          !//for back collision 
          do jj = 1, inny
             do ii = 1, innx
                rholc(ii,jj) =  rhotot(ii,jj,j)
             enddo
          enddo
          !//collect all the potential along column to a single proc. This
          !//function can be replaced by Allgather along column with doubled size.
          call guardexch2d(rholc,phiout,innx,inny,nytot,npyhalf,myidy,commcol)
          !//for front collision
          do jj = 1, inny
             do ii = 1, innx
                rholc(ii,jj) =  rhotot2(ii,jj,j)
             enddo
          enddo
          !//collect all the potential along column to a single proc. This
          !//function can be replaced by Allgather along column with doubled size.
          call guardexch2d(rholc,phiout2,innx,inny,nytot,npyhalf,myidy,commcol)

          !//interpolate field of each slice to slice particles
          !//here, we have used a linear interpolation between the 
          !//slice front and the slice back collision field.
          z2 = zslice2(j1)
          zthk = zhalf(i1)
          sbk = (zslice(i1)-zthk-z2)/2
          call scatter2d2zTsc(maxPtclSlice,nptlcslice,innx,&
               hx,hy,xmin,ymin,ray,phiout,phiout2,xtmptmp,coef,&
               nytot,z2,zthk,sbk)

          do k = 1, nptlcslice
             beamslice(:,k,i1) = ray(:,k)
          enddo
       enddo
       !drift to the field interpolation point
       sign = -1
       if(nslice.gt.1) then
          call driftCP(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
               istart,jstart,jend,zslice,zslice2,countlc,maxPtclSlice)
       endif

    enddo

    deallocate(grn)

    !        call MPI_BARRIER(mpicommwd,ierr)
    !//copy each slice particle to the beam bunch 
    ipt = 0
    do i = 1, nslice
       do j = 1, countlc(i)
          ipt = ipt + 1
          Pts(:,ipt) = beamslice(:,j,i)
       enddo
    enddo

  end subroutine bbmSlcLinear_Beambeam



  !//zslice is the z location of itself and zslice2 is
  !the z location of the opposite beam
  !//zhlf is the half thickness of the slice.
  subroutine slicer2z(Pts1in,nptlc,nslice,nslice2,maxsptcl,&
       nz,zmin,zmax,beamslice,weight,zslice,zslice2,zhlf,zhlf2,&
       nptot,count,countlc,myidy,npyhalf,commcol,commrow)
    implicit none
    include 'mpif.h'
    integer, parameter :: Ndiv = 10
    integer, intent(in) :: nptlc,nslice,nslice2,maxsptcl,nptot,nz,myidy,&
         commcol,commrow,npyhalf
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, intent(in) :: zmin,zmax
    double precision, dimension(6,maxsptcl,nslice),&
         intent(inout) :: beamslice
    double precision, dimension(nslice),intent(out) :: weight,zslice,zhlf
    double precision, dimension(nslice2),intent(out) :: zslice2,zhlf2
    integer, dimension(nslice), intent(out) :: count,countlc
    double precision, dimension(2*nslice)  :: fcount,fcountlc
    double precision, dimension(0:nslice) :: rangez
    double precision, dimension(2*nz) :: rhozlc,rhoz
    double precision, dimension(nz*Ndiv) :: rhoznew
    integer :: i,ierr,ii,iz,iz1,ip,op,mydes,mysou
    double precision :: ab,sumint,hnew1,hnew2,interval0,hz,hzi
    integer status(MPI_STATUS_SIZE)

    op = 1
    if(nslice.gt.1) then
       !//find the charge density along z
       hz = (zmax-zmin)/(nz-1)
       hzi = 1.0/hz
       rhozlc = 0.0
       if(myidy.lt.npyhalf) then
          do i = 1, nptlc
             iz=int((Pts1in(5,i)-zmin)*hzi + 1)
             ab=(zmin-Pts1in(5,i)+iz*hz)*hzi
             iz1 = iz + 1
             rhozlc(iz) = rhozlc(iz) + ab
             rhozlc(iz1) = rhozlc(iz1) + (1.0-ab)
          enddo
       else
          do i = 1, nptlc
             iz=int((Pts1in(5,i)-zmin)*hzi + 1) !//nz is for group 2
             ab=(zmin-Pts1in(5,i)+iz*hz)*hzi
             iz1 = iz + 1 
             rhozlc(iz+nz) = rhozlc(iz+nz) + ab
             rhozlc(iz1+nz) = rhozlc(iz1+nz) + (1.0-ab)
          enddo
       endif
       rhozlc = rhozlc/nptot

       call MPI_ALLREDUCE(rhozlc,rhoz,2*nz,MPI_DOUBLE_PRECISION,&
            MPI_SUM,mpicommwd,ierr)

       !//find totol density, should be 1
       sumint = 0.0
       if(myidy.lt.npyhalf) then
          do i = 1, nz
             sumint = sumint + rhoz(i)
          enddo
       else
          do i = 1, nz
             sumint = sumint + rhoz(i+nz)
          enddo
       endif

       !//subdivided the charge density for finer resolution.
       if(myidy.lt.npyhalf) then
          do i = 1, nz*Ndiv
             ii = (i-1)/Ndiv + 1
             rhoznew(i) = rhoz(ii)/Ndiv
          enddo
       else
          do i = 1, nz*Ndiv
             ii = (i-1)/Ndiv + 1 + nz
             rhoznew(i) = rhoz(ii)/Ndiv
          enddo
       endif
       hnew1 = hz/(Ndiv-1)/2
       hnew2 = hz/Ndiv

       interval0 = sumint/nslice
       sumint = 0.0
       !//find the range of each slice
       do i = 1, Ndiv*nz
          sumint = sumint + rhoznew(i)
          ip = int(sumint/interval0)
          if(ip.ne.nslice) then
             if(i.le.Ndiv) then
                rangez(ip+1) = zmin + (i-1)*hnew1 + 0.5*hnew1
             else if((i.gt.Ndiv).and.(i.le.Ndiv*nz-Ndiv)) then
                rangez(ip+1) = zmin +0.5*hz+(i-Ndiv)*hnew2 + 0.5*hnew2
             else
                rangez(ip+1) = zmin + (nz-1)*hz - 0.5*hz + &
                     (i+Ndiv-Ndiv*nz)*hz/Ndiv
             endif
          endif
       enddo
       rangez(nslice) = zmax
       rangez(0) = zmin


       !//find the # of particle per slice and put the particles into
       !//single slice. Here, zmax corresponds to slice # 1.
       do i = 1, nslice
          countlc(i) = 0
       enddo
       do i = 1,nptlc
          if(Pts1in(5,i).gt.rangez(nslice-1)) then
             !countlc(nslice) = countlc(nslice) + 1
             countlc(1) = countlc(1) + 1
             if(countlc(1).gt.maxsptcl) then
                print*,"overflow local slice limit in nslice: ",countlc(1)
                stop
             endif
             beamslice(:,countlc(1),1) = Pts1in(:,i)
          else if(Pts1in(5,i).le.rangez(1)) then
             !countlc(1) = countlc(1) + 1
             countlc(nslice) = countlc(nslice) + 1
             if(countlc(nslice).gt.maxsptcl) then
                print*,"overflow local slice limit in nslice: ",countlc(nslice)
                stop
             endif
             beamslice(:,countlc(nslice),nslice) = Pts1in(:,i)
          else
             do ii = 2,nslice-1
                if((Pts1in(5,i).gt.rangez(ii-1)).and.&
                     (Pts1in(5,i).le.rangez(ii)) ) then
                   countlc(nslice-ii+1) = countlc(nslice-ii+1) + 1
                   if(countlc(nslice-ii+1).gt.maxsptcl) then
                      print*,"overflow local slice limit in ii: ",ii,&
                           countlc(nslice-ii+1)
                      stop
                   endif
                   beamslice(:,countlc(nslice-ii+1),nslice-ii+1)=Pts1in(:,i)
                endif
             enddo
          endif
       enddo

       fcountlc = 0
       if(myidy.lt.npyhalf) then
          do i = 1, nslice
             fcountlc(i) = countlc(i)
          enddo
       else
          do i = 1, nslice
             fcountlc(i+nslice) = countlc(i)
          enddo
       endif
       call MPI_ALLREDUCE(fcountlc,fcount,2*nslice,MPI_DOUBLE_PRECISION,&
            MPI_SUM,mpicommwd,ierr)

       if(myidy.lt.npyhalf) then
          do i = 1, nslice
             count(i) = int(fcount(i) + 0.0001) !//0.0001 for roundoff
             weight(i) = fcount(i)/nptot
          enddo
       else
          do i = 1, nslice
             count(i) = int(fcount(i+nslice) + 0.0001) !//0.0001 for roundoff
             weight(i) = fcount(i+nslice)/nptot
          enddo
       endif

       do i = 1, nslice
          !zslice(i) = (rangez(i-1)+rangez(i))/2 !//this is a terrible bug.
          zslice(nslice-i+1) = (rangez(i-1)+rangez(i))/2 !//centroid of slice
          zhlf(nslice-i+1) = (rangez(i)-rangez(i-1))/2 !//half thickness of slice
       enddo

       !//send local zslice to the another group zslice2
       if(myidy.lt.npyhalf) then
          mydes = myidy + npyhalf
          mysou = myidy + npyhalf
       else
          mydes = myidy - npyhalf
          mysou = myidy - npyhalf
       endif
       call MPI_SEND(zslice,nslice,MPI_DOUBLE_PRECISION,mydes,1,commcol,&
            ierr)
       call MPI_RECV(zslice2,nslice2,MPI_DOUBLE_PRECISION,mysou,1,commcol,&
            status,ierr)
       call MPI_SEND(zhlf,nslice,MPI_DOUBLE_PRECISION,mydes,1,commcol,&
            ierr)
       call MPI_RECV(zhlf2,nslice2,MPI_DOUBLE_PRECISION,mysou,1,commcol,&
            status,ierr)

    else !//single slice case
       do i = 1, nptlc
          beamslice(:,i,nslice) = Pts1in(:,i)
       enddo
       weight(nslice) = 1.0
       rangez(nslice) = zmax
       rangez(0) = zmin
       zslice(1) = 0.0
       zslice2(1) = 0.0
       zhlf(1) = 0.0
       zhlf2(1) = 0.0
       countlc = nptlc
       !call MPI_ALLREDUCE(countlc,count,nslice,MPI_INTEGER,&
       !                   MPI_SUM,mpicommwd,ierr)
       count(1) = nptot
    endif

  end subroutine slicer2z


  !calculate strong-strong beam-beam interaction using a soft-Gaussian model
  subroutine bbmSlcnewGauss_Beambeam(nprocrow,npyhalf,&
       nslice1,nslice2,nslice,nsliceO,jendmax,nptlc,myidx,myidy,Pts,&
       zmin,zmax,nptot,coef,commrow,commcol,nz,&
       maxPtclSlice,curin,flglum,lum3d,nptot2)
    implicit none
    integer, intent(in) :: nprocrow,npyhalf,nslice1,&
         nslice2,myidx,myidy,nptot,commrow,commcol,nz,&
         maxPtclSlice,nslice,nsliceO,jendmax,flglum,nptot2
    integer, intent(inout) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts
    double precision, intent(out) :: lum3d
    double precision, intent(in) :: zmin,coef,zmax,curin
    double precision, dimension(nsliceO) :: zslice2,hzslice2,count2
    double precision, dimension(nsliceO) :: xslcent,yslcent,sigxsl,sigysl
    double precision, dimension(nslice) :: weight,zslice,hzslice
    double precision, dimension(6,maxPtclSlice) :: ray
    double precision, dimension(7,maxPtclSlice,nslice) :: beamslice
    integer, dimension(nslice) :: countlc,count
    integer :: i,j,nptlcslice,sign,ii,jj,ipt,k,ierr,&
         i1,j1,jend,istart,jstart,myjend,itmp1,itmp2,rowend
    integer :: npbc,nsxy1,nsxy2,nxpylc2,nxtot,itb,nxlum,nylum
    double precision :: xtmptmp,lumtmp,ww
    real*8, dimension(2) :: center,sigma

    lum3d = 0.0
    nxlum = 128
    nylum = 128

    !//set up the parameters for the lower group and upper group processors.

    !//here, we have used equal area weight
    !//copy particles from a bunch into slices. Each slice has
    !//roughly same number of particles.
    call slicernew4gaus(Pts,nptlc,nslice,nsliceO,maxPtclSlice,nz,zmin,&
         zmax,beamslice,weight,zslice,zslice2,hzslice,hzslice2,nptot,count,countlc,myidy,&
         npyhalf,commcol,commrow,xslcent,yslcent,sigxsl,sigysl,count2)

    call MPI_BARRIER(mpicommwd,ierr)

    !//loop through the total overlaping steps.
    do i = 1, nslice1+nslice2-1
       !//find the number of colliding slices 
       if(i.le.nslice1) then
          jend = min(i,nslice2) !//jend is the number of the overlaped slices.
          istart = i+1  !//istart is the starting slice + 1 in beam1
          jstart = 0    !//jstart is is the starting slice - 1 in beam2
       else
          jend = min((nslice2-i+nslice1),nslice1)
          istart = nslice1+1
          jstart = i-nslice1
       endif

       !drift to the field interpolation point
       if(nslice.gt.1) then
          sign = 1
          call driftCP(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
               istart,jstart,jend,zslice,zslice2,countlc,maxPtclSlice)
       endif

       !interpolate the beam-beam field onto particles on each slice.
       do j = 1, jend
          if(myidy.lt.npyhalf) then
             i1 = istart - j
             j1 = jstart + j
          else
             i1 = jstart + j
             j1 = istart - j
          endif

          nptlcslice = countlc(i1)
          do k = 1, nptlcslice
             ray(:,k) = beamslice(1:6,k,i1)
          enddo

          sigma(1) = sigxsl(j1) 
          sigma(2) = sigysl(j1)  
          center(1) = xslcent(j1) 
          center(2) = yslcent(j1)  
          ww = count2(j1)/nptot2
          call scatter2dgauss(nptlcslice,nptlcslice,& 
               ray,ww,coef,sigma,center) 

          do k = 1, nptlcslice
             beamslice(1:6,k,i1) = ray(:,k)
          enddo

          !//calculate the 3d luminosity
          if(flglum.eq.1) then
             call luminosity2G3d_Output(ray,nptlcslice,nxlum,nylum,&
                  curin,weight(i1),count(i1),myidy,npyhalf,lumtmp)
             lum3d = lum3d + lumtmp
          endif

       enddo

       !//drift all particles back to their original longitudinal location
       if(nslice.gt.1) then 
          sign = -1 
          call driftCP(beamslice,sign,nslice,nsliceO,myidy,npyhalf,&
               istart,jstart,jend,zslice,zslice2,countlc,maxPtclSlice)
       endif
    enddo

    call MPI_BARRIER(mpicommwd,ierr)
    !//copy each slice particle to the beam bunch
    ipt = 0
    do i = 1, nslice
       do j = 1, countlc(i)
          !ipt = ipt + 1
          ipt = beamslice(7,j,i)+0.01
          Pts(:,ipt) = beamslice(1:6,j,i)
       enddo
    enddo

  end subroutine bbmSlcnewGauss_Beambeam



  !slicing the beams longitudinal with non-uniform grid size trying to
  !maintain the same weight for each slice.
  !The centers and the sigmas of each slice from the opposite beam are
  !caclulated for the soft-Gaussian model usage.
  subroutine slicernew4gaus(Pts1in,nptlc,nslice,nslice2,maxsptcl,&
       nz,zmin,&
       zmax,beamslice,weight,zslice,zslice2,hzslice,hzslice2,nptot,count,countlc,myidy,&
       npyhalf,commcol,commrow,xslcent,yslcent,sigxsl,sigysl,count2)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc,nslice,nslice2,maxsptcl,nptot,nz,myidy,&
         commcol,commrow,npyhalf
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, intent(in) :: zmin,zmax
    double precision, dimension(7,maxsptcl,nslice),&
         intent(inout) :: beamslice
    double precision, dimension(nslice),intent(out) :: weight,zslice,hzslice
    double precision, dimension(nslice2),intent(out) :: xslcent,yslcent,sigxsl,sigysl
    double precision, dimension(nslice2),intent(out) :: zslice2,hzslice2,count2
    integer, dimension(nslice), intent(out) :: count,countlc
    double precision, dimension(5,(nslice+nslice2))  :: fcount,fcountlc
    double precision, dimension(0:nslice) :: rangez
    double precision, dimension(2*nz) :: rhozlc,rhoz
    real*8, dimension(nslice2) :: sigslx,sigslx2,sigslxglb,sigslx2glb
    real*8, dimension(nslice2) :: sigsly,sigsly2,sigslyglb,sigsly2glb
    integer :: i,ierr,islice,ii,iz,iz1,ip,op,mydes,mysou,ip0
    double precision :: ab,sumint,hnew1,hnew2,interval0,hz,hzi
    real*8 :: sum0,sum1,sumip,z0,z1,zip
    integer status(MPI_STATUS_SIZE)
    integer :: nsendata

    !print*,"zmin: ",zmin,zmax
    op = 1
    if(nslice.gt.1) then
       !//find the charge density along z
       hz = (zmax-zmin)/(nz-1)
       hzi = 1.0/hz
       rhozlc = 0.0
       if(myidy.lt.npyhalf) then
          do i = 1, nptlc
             iz=(Pts1in(5,i)-zmin)*hzi + 1
             ab=(zmin-Pts1in(5,i)+iz*hz)*hzi
             iz1 = iz + 1
             rhozlc(iz) = rhozlc(iz) + ab
             rhozlc(iz1) = rhozlc(iz1) + (1.0-ab)
          enddo
       else
          do i = 1, nptlc
             iz=(Pts1in(5,i)-zmin)*hzi + 1 
             ab=(zmin-Pts1in(5,i)+iz*hz)*hzi
             iz1 = iz + 1 
             !//nz is for group 2, and nz must be the same for 2 beams
             rhozlc(iz+nz) = rhozlc(iz+nz) + ab
             rhozlc(iz1+nz) = rhozlc(iz1+nz) + (1.0-ab)
          enddo
       endif
       rhozlc = rhozlc/nptot

       call MPI_ALLREDUCE(rhozlc,rhoz,2*nz,MPI_DOUBLE_PRECISION,&
            MPI_SUM,mpicommwd,ierr)

       !//find totol density, should be 1
       sumint = 0.0
       if(myidy.lt.npyhalf) then
          do i = 1, nz
             sumint = sumint + rhoz(i)
          enddo
       else
          do i = 1, nz
             sumint = sumint + rhoz(i+nz)
          enddo
       endif
       !print*,"sumint: ",sumint,myidy,npyhalf,nptlc
       !print*,"rhoz2: ",rhoz

       !using the following way, the nz needs to be sufficiently large
       interval0 = sumint/nslice
       sumint = 0.0
       ip0 = 0
       do i = 1, nz
          sumint = sumint + rhoz(i)
          ip = sumint/interval0
          if((ip.ne.ip0).and.(ip.lt.nslice)) then
             z1 = zmin + (i-1)*hz
             z0 = zmin + (i-2)*hz
             sum1 = sumint
             sum0 = sumint-rhoz(i)
             sumip = ip*interval0
             zip = z0 + (sumip-sum0)*(z1-z0)/(sum1-sum0)
             rangez(ip) = zip 
             ip0 = ip
          endif
       enddo

       rangez(nslice) = zmax
       rangez(0) = zmin

       !print*,"rangez: ",rangez

       !//find the # of particle per slice and put the particles into
       !//single slice. Here, zmax corresponds to slice # 1.
       do i = 1, nslice
          countlc(i) = 0
       enddo
       sigslx = 0.0d0
       sigslx2 = 0.0d0
       sigsly = 0.0d0
       sigsly2 = 0.0d0
       do i = 1,nptlc
          if(Pts1in(5,i).gt.rangez(nslice-1)) then
             !countlc(nslice) = countlc(nslice) + 1
             countlc(1) = countlc(1) + 1
             if(countlc(1).gt.maxsptcl) then
                print*,"overflow local slice limit in nslice: ",countlc(1)
                stop
             endif
             beamslice(1:6,countlc(1),1) = Pts1in(:,i)
             beamslice(7,countlc(1),1) = i
             sigslx(1) = sigslx(1) + Pts1in(1,i)
             sigslx2(1) = sigslx2(1) + Pts1in(1,i)**2
             sigsly(1) = sigsly(1) + Pts1in(3,i)
             sigsly2(1) = sigsly2(1) + Pts1in(3,i)**2
          else if(Pts1in(5,i).le.rangez(1)) then
             !countlc(1) = countlc(1) + 1
             countlc(nslice) = countlc(nslice) + 1
             if(countlc(nslice).gt.maxsptcl) then
                print*,"overflow local slice limit in nslice: ",countlc(nslice)
                stop
             endif
             beamslice(1:6,countlc(nslice),nslice) = Pts1in(:,i)
             beamslice(7,countlc(nslice),nslice) = i
             sigslx(nslice) = sigslx(nslice) + Pts1in(1,i)
             sigslx2(nslice) = sigslx2(nslice) + Pts1in(1,i)**2
             sigsly(nslice) = sigsly(nslice) + Pts1in(3,i)
             sigsly2(nslice) = sigsly2(nslice) + Pts1in(3,i)**2
          else
             do ii = 2,nslice-1
                if((Pts1in(5,i).gt.rangez(ii-1)).and.&
                     (Pts1in(5,i).le.rangez(ii)) ) then
                   countlc(nslice-ii+1) = countlc(nslice-ii+1) + 1
                   if(countlc(nslice-ii+1).gt.maxsptcl) then
                      print*,"overflow local slice limit in ii: ",ii,&
                           countlc(nslice-ii+1)
                      stop
                   endif
                   beamslice(1:6,countlc(nslice-ii+1),nslice-ii+1)=Pts1in(:,i)
                   beamslice(7,countlc(nslice-ii+1),nslice-ii+1)=i
                   sigslx(nslice-ii+1) = sigslx(nslice-ii+1) + Pts1in(1,i)
                   sigslx2(nslice-ii+1) = sigslx2(nslice-ii+1) + Pts1in(1,i)**2
                   sigsly(nslice-ii+1) = sigsly(nslice-ii+1) + Pts1in(3,i)
                   sigsly2(nslice-ii+1) = sigsly2(nslice-ii+1) + Pts1in(3,i)**2
                endif
             enddo
          endif
       enddo
       !print*,"countlc: ",myidy,countlc,sigslx2

       fcountlc = 0.0d0
       if(myidy.lt.npyhalf) then
          do i = 1, nslice
             fcountlc(1,i) = countlc(i)
             fcountlc(2,i) = sigslx(i)
             fcountlc(3,i) = sigslx2(i)
             fcountlc(4,i) = sigsly(i)
             fcountlc(5,i) = sigsly2(i)
          enddo
       else
          do i = 1, nslice
             fcountlc(1,i+nslice2) = countlc(i)
             fcountlc(2,i+nslice2) = sigslx(i)
             fcountlc(3,i+nslice2) = sigslx2(i)
             fcountlc(4,i+nslice2) = sigsly(i)
             fcountlc(5,i+nslice2) = sigsly2(i)
          enddo
       endif
       nsendata = 5*(nslice+nslice2)
       call MPI_ALLREDUCE(fcountlc,fcount,nsendata,MPI_DOUBLE_PRECISION,&
            MPI_SUM,mpicommwd,ierr)

       if(myidy.lt.npyhalf) then
          do i = 1, nslice
             count(i) = fcount(1,i) + 0.0001 !//0.0001 for roundoff
             weight(i) = fcount(1,i)/nptot
             !              print*,"count1: ",count(i),fcountlc(i),myidy,nptot
          enddo
          do i = 1, nslice2
             count2(i) = fcount(1,i+nslice) + 0.0001 !//0.0001 for roundoff
             if(count2(i).eq.0) print*,"something wrong in slicernew4gaus"
             !here the center and sigma are the opposite beam not itself
             sigslxglb(i) = fcount(2,i+nslice)/count2(i)
             sigslx2glb(i) = fcount(3,i+nslice)/count2(i)
             sigslyglb(i) = fcount(4,i+nslice)/count2(i)
             sigsly2glb(i) = fcount(5,i+nslice)/count2(i)
          enddo
       else
          do i = 1, nslice
             count(i) = fcount(1,i+nslice2) + 0.0001 !//0.0001 for roundoff
             weight(i) = fcount(1,i+nslice2)/nptot
             if(count(i).eq.0) print*,"something wrong in slicernew4gaus"
             !              print*,"count2: ",count(i),fcountlc(i+nslice),myidy,nptot
          enddo
          do i = 1, nslice2
             count2(i) = fcount(1,i) + 0.0001 !//0.0001 for roundoff
             !here the center and sigma are the opposite beam not itself
             if(count2(i).eq.0) print*,"something wrong in slicernew4gaus"
             sigslxglb(i) = fcount(2,i)/count2(i)
             sigslx2glb(i) = fcount(3,i)/count2(i)
             sigslyglb(i) = fcount(4,i)/count2(i)
             sigsly2glb(i) = fcount(5,i)/count2(i)
          enddo
       endif
       !print*,"sig: ",sqrt(sigslx2glb)

       do i = 1, nslice
          zslice(nslice-i+1) = (rangez(i-1)+rangez(i))/2
          hzslice(nslice-i+1) = rangez(i)-rangez(i-1)
       enddo
       do i = 1, nslice2
          sigxsl(i)  = sqrt(sigslx2glb(i)-sigslxglb(i)**2)
          sigysl(i)  = sqrt(sigsly2glb(i)-sigslyglb(i)**2)
          xslcent(i) = sigslxglb(i)
          yslcent(i) = sigslyglb(i)
       enddo

       !//send local zslice to the another group zslice2
       if(myidy.lt.npyhalf) then
          mydes = myidy + npyhalf
          mysou = myidy + npyhalf
       else
          mydes = myidy - npyhalf
          mysou = myidy - npyhalf
       endif
       !print*,"myidy: ",myidy,mydes,mysou,npyhalf
       call MPI_SEND(zslice,nslice,MPI_DOUBLE_PRECISION,mydes,1,commcol,&
            ierr)
       call MPI_RECV(zslice2,nslice2,MPI_DOUBLE_PRECISION,mysou,1,commcol,&
            status,ierr)
       call MPI_SEND(hzslice,nslice,MPI_DOUBLE_PRECISION,mydes,1,commcol,&
            ierr)
       call MPI_RECV(hzslice2,nslice2,MPI_DOUBLE_PRECISION,mysou,1,commcol,&
            status,ierr)

    else !//single slice case !not sure whether it is correct or not if one beam has 1 slice the other beam has N slice
       print*,"this option in slicernew4gaus is not valid yet!"
       stop
       do i = 1, nptlc
          beamslice(1:6,i,nslice) = Pts1in(:,i)
          beamslice(7,i,nslice) = i
       enddo
       weight(nslice) = 1.0
       rangez(nslice) = zmax
       rangez(0) = zmin
       zslice(1) = 0.0
       zslice2(1) = 0.0
       countlc = nptlc
       !call MPI_ALLREDUCE(countlc,count,nslice,MPI_INTEGER,&
       !                   MPI_SUM,mpicommwd,ierr)
       count(1) = nptot
       hzslice(1) = zmax - zmin
       hzslice2(1) = zmax - zmin !not correct but not used
    endif
    !print*,"rangez: ",rangez

  end subroutine slicernew4gaus



  !//find the strong-strong soft Gaussian beam-beam kick in a 2d one slice model.
  subroutine bb1SlcGaus_Beambeam(npyhalf,nptlc,myidy,Pts1,nptot,coef)
    implicit none
    integer, intent(in) :: npyhalf,myidy,nptot
    integer, intent(inout) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1
    double precision, intent(in) :: coef
    double precision :: weight
    real*8, dimension(2) :: center,sigma

    call findmomO_Utility(Pts1,nptlc,nptlc,nptot,npyhalf,myidy,center,sigma)
    weight = 1.0d0
    call scatter2dgauss(nptlc,nptlc,Pts1,weight,coef,sigma,center)

  end subroutine bb1SlcGaus_Beambeam


end module Beambeamclass
