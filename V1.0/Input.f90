!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! Inputclass: Input class in I/O module of CONTROL layer. 
! Version: 1.0
! Author: Ji Qiang, LBNL 
! Description: This class defines functions to input the global
!              beam and computational parameters and the lattice input
!              parameters in the accelerator.
! Comments:
!----------------------------------------------------------------

module Inputclass
  use Pgrid2dclass

contains
  ! Start MPI
  subroutine init_Input(time)
    implicit none
    include 'mpif.h'
    double precision, intent(out) :: time
    integer :: ierr

    ! start MPI library.
    call MPI_INIT(ierr)
    time = MPI_WTIME()
    ! for measurement of memory usage.
    !        call system_stats()

  end subroutine init_Input



  ! Input all parameters except beam line element parameters.
  subroutine in_Input(filename,odim,onp,onx,ony,onz,onslice,&
       onturn,oflagdist,distparam,nparam,obcurr,obkenergy,obmass,&
       obcharge,otunex,otuney,otunez,oax,oay,obx,oby,oemitx,oemity,&
       otauy,odampart,osigz,ospop,onpcol,onprow,oclose,oclosess,&
       oifdbk,oisweep,oswrad,oswtune,oiclosq,onmeet,osigmaz,osigmapz,&
       oalpha,ophi,oflgrange,oflgslice,oflgmp2nd,oflgmp3rd,oflgmp4th,oncl,&
       oextrange,oflagint,oflagwk,oqx,oqy,ohcv,onfrqlum,otauz,odampflag,&
       onbunch,onsteps,oistart,otmax,saveAll,oNresamp,ptRate,ptFrac,momOutRate,gain,betaCrab,crabFreq,&
       crabvnorm,cfreq2,cvnorm2,np_xy,len_xy,repRate_xy, &
       ogainx,ogainy,ofrqhighx,ofrqhighy)

    implicit none
    include 'mpif.h'
    character*8 :: filename
    integer, intent(out) :: odim,onp,onx,ony,onz,onslice,onturn,oflagdist
    integer, intent(out) :: onpcol,onprow,oifdbk,oisweep,oiclosq,onmeet,oNresamp 
    integer, intent(in) :: nparam
    double precision, dimension(nparam), intent(out) :: distparam
    double precision, dimension(2), intent(out) :: oclose,oclosess
    double precision, intent(out) :: obcurr,obkenergy,obmass
    double precision, intent(out) :: obcharge,otunex,otuney,otunez,&
         oax,oay,obx,oby,oemitx,oemity,otauy,odampart,otauz,osigz,ospop,&
         oswrad,oswtune,osigmaz,osigmapz,oalpha,ophi,oqx,oqy,ohcv,otmax
    integer, intent(out) :: oflgrange,oflgslice,oflgmp2nd,oflgmp3rd,&
         oflgmp4th,oncl,oflagint,oflagwk,onfrqlum,odampflag,onbunch,onsteps,oistart,saveAll
    double precision, dimension(6), intent(out) :: oextrange
    integer, intent(out) :: ptRate,ptFrac,momOutRate
    double precision, intent(out), dimension(4) :: gain
    double precision, intent(out) :: betaCrab,crabFreq,crabVnorm,cfreq2,cvnorm2
    integer, intent(out) :: np_xy,len_xy,repRate_xy
    real*8, intent(out) :: ogainx,ogainy,ofrqhighx,ofrqhighy

    integer :: my_rank,ierr,np

    call MPI_COMM_RANK(mpicommwd,my_rank,ierr)
    call MPI_COMM_SIZE(mpicommwd,np,ierr)

    if(my_rank.eq.0) then
       open(unit=13,file=filename,status='old')

       read(13,*)onpcol,onprow
       !onprow has always to be 1 in the beam-beam simulation.
       !if(onprow.ne.1) then
       !  onprow = 1
       !  onpcol = np
       !endif
       read(13,*)odim, onp 
       read(13,*)onx, ony, onz 
       read(13,*)onslice
       read(13,*)onturn
       read(13,*)oflagdist
       distparam = 0.0
       read(13,*)distparam(1),distparam(2),distparam(3),distparam(4),&
            distparam(5),distparam(6),distparam(7)
       read(13,*)distparam(8),distparam(9),distparam(10),&
            distparam(11),distparam(12),distparam(13),distparam(14)
       read(13,*)distparam(15),distparam(16),&
            distparam(17),distparam(18),distparam(19),distparam(20),distparam(21)
       read(13,*)obcurr,obkenergy,obmass,obcharge
       read(13,*)otunex,otuney,otunez
       read(13,*)oax,oay
       read(13,*)obx,oby
       read(13,*)oemitx,oemity
       read(13,*)otauy,odampart,otauz,odampflag
       read(13,*)osigz,ospop
       read(13,*)oisweep,oswrad,oswtune
       read(13,*)oifdbk,ogainx,ogainy,ofrqhighx,ofrqhighy
       read(13,*)oiclosq,onmeet
       read(13,*)oclosess(1),oclosess(2)
       read(13,*)osigmaz,osigmapz
       read(13,*)oalpha,ophi
       read(13,*)oflgrange,oflgslice,oflgmp2nd,oflgmp3rd,oflgmp4th
       read(13,*)oncl
       read(13,*)oextrange(1),oextrange(2),oextrange(3),oextrange(4),&
            oextrange(5),oextrange(6)
       read(13,*)oflagint,oflagwk
       read(13,*)oqx,oqy !chromaticity
       read(13,*)ohcv !curvature
       read(13,*)onfrqlum !luminosity sample frequency
       read(13,*)onbunch,onsteps
       read(13,*)oistart,otmax,saveAll
       read(13,*)oNresamp
       read(13,*)ptRate,ptFrac,momOutRate
       read(13,*)gain(1:4)
       read(13,*)betaCrab,crabFreq,crabVnorm,cfreq2,cvnorm2
       read(13,*)np_xy,len_xy,repRate_xy
       close(13)
    endif

    call MPI_BCAST(onpcol,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(onprow,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(odim,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(onp,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(onx,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(ony,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(onz,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(onslice,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(onturn,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(oflagdist,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(oisweep,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(oifdbk,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(oiclosq,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(onmeet,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(distparam(1),nparam,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    oclose(1) = distparam(6)
    oclose(2) = distparam(13)
    !        oclose(1) = 0.0
    !        oclose(2) = 0.0
    call MPI_BCAST(oclosess(1),2,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(obcurr,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(obkenergy,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(obmass,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(obcharge,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(otunex,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(otuney,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(otunez,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(oax,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(oay,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(obx,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(oby,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(oemitx,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(oemity,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(otauy,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(odampart,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(otauz,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(odampflag,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(osigz,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(ospop,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(oswrad,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(oswtune,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(osigmaz,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(osigmapz,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(oalpha,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(ophi,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(oflgrange,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(oflgslice,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(oflgmp2nd,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(oflgmp3rd,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(oflgmp4th,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(oncl,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(oextrange,6,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(oflagint,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(oflagwk,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(onfrqlum,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(oqx,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(oqy,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(ohcv,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(onbunch,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(onsteps,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(oistart,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(saveAll,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(otmax,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(oNresamp,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(ptRate,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(ptFrac,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(momOutRate,1,MPI_INTEGER,0,mpicommwd,ierr)
    call MPI_BCAST(gain,4,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(betaCrab,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(crabFreq,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(crabVnorm,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(cfreq2,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(cvnorm2,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(np_xy,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(len_xy,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(repRate_xy,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(ogainx,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(ogainy,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(ofrqhighx,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)
    call MPI_BCAST(ofrqhighy,1,MPI_DOUBLE_PRECISION,0,mpicommwd,ierr)

  end subroutine in_Input


end module Inputclass
