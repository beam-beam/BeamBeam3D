!----------------------------------------------------------------
! (c) Copyright, 2016 by the Regents of the University of California.
! AccSimulatorclass: BeamBeam3D.
! Version: 1.0
! Developer: Ji Qiang 
! Contributors: 
!...............................................
! Description: This class defines functions to set up the initial beam 
!              particle distribution, field information, computational
!              domain, beam line element lattice and run the dynamics
!              simulation through the system.
! Comments: The coordinate of beam 1 is x,y,z right hand system.
!           The coordinate of beam 2 is x,y,z left hand system.
!           Beam 1 and beam 2 have the same transverse x, y but
!           different longitudinal direction. This follow the 
!           notation of MAD. The coordinate system also follows
!           MAD definition. The # of particles has to be multiple of
!           the half processor. The # of grid has to be multiple
!           of the half processor. For strong/weak beam-beam, the
!           first beam has to be the strong beam, the second beam has
!           to be the weak beam.
!----------------------------------------------------------------

module AccSimulatorclass
  use Pgrid2dclass
  use CompDomclass
  use Timerclass
  use Inputclass
  use Outputclass
  use Distributionclass
  use Linearmapclass
  use Extmap2ndclass
  use Extmap3rdclass
  use Extmap4thclass
  use Dampflucclass
  use Beambeamclass
  use Orbitclass
  use Transferclass
  use DepoScatclass
  use Utilityclass
  use Feedbackclass

  !include 'mpif.h'
  implicit none

  double precision, parameter :: epsilon0 = 8.854187817d-12
  double precision, parameter :: qe = 1.60217733d-19
  double precision, parameter :: clight = 2.99792458d8
  double precision, parameter :: pmass = 1.67262158d-27
  double precision, parameter :: emass = 9.10938188d-31
  double precision, parameter :: z0 = 1/(epsilon0*clight)
  double precision, parameter :: pi = 2d0*asin(1.d0)

  !maximum # of collisions per step
  integer, parameter :: MaxCol = 100
  !maximum # of steps per turn
  integer, parameter :: MaxStep = 100
  !# of phase dim., num. total and local particles, int. dist. 
  !and restart switch, error study switch, substep for space-charge
  !switch
  !6d phase space, # of turns, # of PEs in y, # of PEs in z
  integer, private :: Dim,Nturn,npcol,nprow
  !# of total particles, # of local particles, switch of
  !initial distribution, global # of mesh along x, y, z,
  !local # of mesh along y, # of slice 
  integer, private :: Np1,Nplocal1,Flagdist1,Nx1,Ny1,Nz1,&
       Nylocal1,Nslice1
  integer, private :: Np2,Nplocal2,Flagdist2,Nx2,Ny2,Nz2,&
       Nylocal2,Nslice2
  !global # of particles, loal # of particles, switch of
  !initial distribution, # of sampling frequency for luminosity output.
  integer, private :: Np,NpO,Nplocal,Flagdist,nfrqlum

  !beam current, kin. energy, part. mass, and charge for one bunch.
  double precision, private :: Bcurr1,Bkenergy1,Bmass1,Bcharge1
  !beam current, kin. energy, part. mass, and charge for a list of bunches
  double precision, allocatable,dimension(:) :: Bcurrlist1,Bmasslist1,Bchargelist1
  !horizontal tune, vertical tune, synchrotron tune, alpha x, alpha y, beta x,
  !beta y, emittance x, emittance y, damping time y, ...............
  double precision, private :: tunex1,tuney1,tunez1,ax1,ay1,&
       bx1,by1,emitx1,emity1,tauy1,dampart1,tauz1,sigz1,spop1,Bcurr
  double precision, private :: Bcurr2,Bkenergy2,Bmass2,Bcharge2
  double precision, allocatable,dimension(:) :: Bcurrlist2,Bmasslist2,Bchargelist2
  double precision, private :: tunex2,tuney2,tunez2,ax2,ay2,&
       bx2,by2,emitx2,emity2,tauy2,dampart2,tauz2,sigz2,spop2

  double precision, private :: tunex,tuney,tunez,ax,ay,bx,by,tauy,dampart,tauz
  double precision, private :: bet1,bet2,bet

  !x chromaticity, y chromaticity, curvature of ring
  double precision, private :: qx1,qx2,qy1,qy2,hcuv1,hcuv2,qx,qy,hcuv

  !conts. in init. dist.
  double precision, private, dimension(21) :: distparam1,distparam2,&
       distparam
  !average rms sigma size of the opposite beam.
  double precision :: sigopp

  !6 rms information for each beam 
  double precision, private, dimension(6) :: sigma1,sigma2,sigma

  !static closed orbit
  double precision, private, dimension(2) :: close1ss,close2ss,closess
  !closed orbit
  double precision, private, dimension(2) :: close1,close2,close2g
  !beam sweeping amplitude, tune
  double precision :: swrad1,swrad2,swtune1,swtune2
  double precision :: swrad,swtune
  !rms z and pz
  double precision :: sigmaz1,sigmapz1,sigmaz2,sigmapz2

  !switches of feedback, beam sweeping, orbit closing, # of turns for closing
  integer :: ifdbk,isweep1,isweep2,iclosq1,iclosq2,nmeet,isweep,iclosq

  !coefficients before the beam-beam kick (from 2D poisson solver).
  double precision :: coef1,coef2,coef

  !collision angles ....................
  double precision :: alpha1,alpha2,alpha,phi1,phi2,phi

  !1d logical processor array.
  type (Pgrid2d), private :: grid2d

  !beam particle object for one bunch
  double precision, pointer, dimension(:,:) :: Bpts,Bpts1,Bpts2
  !temporary beam object for 2nd order and 3rd map
  double precision, pointer, dimension(:,:) :: Bptstmp
  !beam particle object for an array of bunches
  double precision, pointer, dimension(:,:,:) :: BptsArry

  !common geometry object for both beams.
  type (CompDom), private :: Ageom

  !linear map object.
  type (Linearmap), private :: Alinearmap1,Alinearmap2,Alinearmap
  !type (Extmap2nd), dimension(MaxCol) :: Amap2nd,Amap2nd1,Amap2nd2
  type (Extmap2nd),allocatable,dimension(:) :: Amap2nd,Amap2nd1,Amap2nd2
  type (Extmap2nd),allocatable,dimension(:,:) :: Amap2ndArry
  !3rd tensor map object
  !type (Extmap3rd), dimension(MaxCol) :: Amap3rd,Amap3rd1,Amap3rd2
  type (Extmap3rd),allocatable,dimension(:) :: Amap3rd,Amap3rd1,Amap3rd2
  type (Extmap3rd),allocatable,dimension(:,:) :: Amap3rdArry
  !4th order tensor map object
  !type (Extmap4th), dimension(MaxCol) :: Amap4th,Amap4th1,Amap4th2
  type (Extmap4th),allocatable,dimension(:) :: Amap4th,Amap4th1,Amap4th2
  type (Extmap4th),allocatable,dimension(:,:) :: Amap4thArry

  !a list of slices at each collision points.
  integer, dimension(MaxCol) :: Nslice1list,Nslice2list
  !a list of the bunch id at each collision point
  integer, dimension(MaxCol,MaxStep) :: idblist,idblistO
  !# of global particle per bunch, # of local particle per bunch.
  integer, allocatable,dimension(:) :: Nplist,Nplclist

  !radiation damping and quantumn fluctuation object.
  type (Dampfluc), private :: Adampfluc,Adampfluc1,Adampfluc2

  !switch flag for using "fixed range", "variable slice",
  !"linear map", "3rd tensor map", group of PEs, strong-weak interaction
  integer :: flagrange,flagslice1,flagslice2,flagmap2nd,flagmap3rd,&
       flagmap4th,flagGroup,flagwk
  !# of collisions per step, # of bunches, # of steps per turn
  integer :: Ncl,nbunch1,nbunch2,nbunch,nsteps

  !external supplied range of domain. 
  !when "flagrange" is true, the computation will use this range
  !find find the electro-magnetic force from opposite beam.
  double precision, dimension(6) :: Extrange1,Extrange2,Extrange

  !switch flag for the choice of the radiation damping
  !flagdamp = 0, no damping; 1, damping using tauy and damping part;
  !           2, damping with taux, tauy, tauz
  integer :: flagdamp1,flagdamp2,flagdamp

  !parameters for restart function
  integer :: iendleft,idrstart,iticleft,iseedend
  double precision :: tmax
  integer*8 :: icount

  !Interval (in turns) for resampling the particle distribution of beam 1
  integer :: Nresamp
  integer :: SaveAll

  !// turn sample frequency and particle sample frequency for phase space output
  !// sampling rate for ouptut of moments
  integer :: ptRate,ptFrac,momOutRate

  integer :: myid,myidx,myidy,comm2d,commcol,commrow

  !flag for initial offset: 0 - no offset, otherwise - offset
  integer :: flagoffset

  double precision, dimension(4) :: gain1, gain2

  integer :: np_xy, repRate_xy, len_xy  !// variables to control data output for spectral analysis

  double precision, private :: crabFreq, crabVnorm, betaCrab, cfreq2,cvnorm2

  real*8 :: gainx1,gainy1,frqhighx1,frqhighy1
  real*8 :: gainx2,gainy2,frqhighx2,frqhighy2
  real*8 :: gainx,gainy,frqhighx,frqhighy

contains

  !set up objects and parameters.
  subroutine init_AccSimulator(time)
    implicit none
    integer :: ierr,Flagbc,npyhalf
    double precision :: time
    double precision :: t0,ga1,ga2
    double precision :: xrad,yrad,Perdlen,s1,s2,s3,gx,gy
    double precision, dimension(6,6) :: tmpmap2nd
    double precision, dimension(6,6,6) :: tmpmap3rd
    double precision, dimension(6,6,6,6) :: tmpmap4th
    integer :: i,j,ii,jj,kk,ngroup,ipt,i6,ii0 !,ib
    double precision :: tmpvalue
    double precision :: tmpvalue1,tmpvalue2,tmpvalue3,tmpvalue4,&
         tmpvalue5,tmpvalue6,error
    !        double precision, dimension(2) :: tmpsum, tmpsumg
    double precision, dimension(2) :: closetmp
    double precision, dimension(6) :: tmp1
    double precision :: betaCrab1,betaCrab2,crabFreq1,crabFreq2,crabVnorm1,crabVnorm2
    double precision :: cfreq12,cfreq22,cvnorm12,cvnorm22
    double precision :: bbparm1x,bbparm2x,rree1,rree2,bbparm1y,bbparm2y
    integer :: comm_size
    character :: message

    !start up MPI.
    call init_Input(time)
    !set up the MPI communicator
    mpicommwd = MPI_COMM_WORLD
    call MPI_COMM_RANK(mpicommwd,myid,ierr)
    call MPI_COMM_SIZE(mpicommwd,comm_size,ierr)

    ! initialize Timer.
    call construct_Timer(0.0d0)
    call starttime_Timer(t0)

    !flag for offset
    flagoffset = 0

    !-------------------------------------------------------------------
    ! get all global input parameters.
    call in_Input("beam1.in",Dim,Np1,Nx1,Ny1,Nz1,Nslice1,Nturn,&
         Flagdist1,distparam1,21,Bcurr1,Bkenergy1,Bmass1,Bcharge1,&
         tunex1,tuney1,tunez1,ax1,ay1,bx1,by1,emitx1,emity1,tauy1,&
         dampart1,sigz1,spop1,npcol,nprow,close1,close1ss,&
         ifdbk,isweep1,swrad1,swtune1,iclosq1,nmeet,sigmaz1,sigmapz1,&
         alpha1,phi1,flagrange,flagslice1,flagmap2nd,flagmap3rd,flagmap4th,Ncl,&
         Extrange1,flagGroup,flagwk,qx1,qy1,hcuv1,nfrqlum,tauz1,flagdamp1,&
         nbunch1,nsteps,idrstart,tmax,saveAll,Nresamp,ptRate,ptFrac,momOutRate,gain1,betaCrab1,crabFreq1,&
         crabVnorm1,cfreq12,cvnorm12,np_xy,len_xy,repRate_xy, &
         gainx1,gainy1,frqhighx1,frqhighy1)

    !// ensure number of macro particles fits number of cores (per beam)
    tmpvalue = modulo(Np1,npcol*nprow/flagGroup)
    if(myid==0 .and. tmpvalue/=0) then
       Np1 = int(Np1 - tmpvalue)
       print*, "Warning: The number of macro particles in beam 1 was truncated to", &
            Np1, " for an equal particle number on each core."
    endif

    gx = (1.0+ax1*ax1)/bx1
    sigma1(1) = sqrt(emitx1*bx1)
    sigma1(2) = sqrt(emitx1*gx)
    gy = (1.0+ay1*ay1)/by1
    sigma1(3) = sqrt(emity1*by1)
    sigma1(4) = sqrt(emity1*gy)
    sigma1(5) = distparam1(15)
    sigma1(6) = distparam1(16)
    if(close1(1).gt.1.0d-4*sigma1(1)) then  !// check initial offset
       flagoffset = 1
    endif
    if(close1(2).gt.1.0d-4*sigma1(3)) then
       flagoffset = 1
    endif
    ga1 = Bkenergy1/Bmass1 
    bet1 = sqrt(1.0-1.0/ga1/ga1)

    call in_Input("beam2.in",Dim,Np2,Nx2,Ny2,Nz2,Nslice2,Nturn,&
         Flagdist2,distparam2,21,Bcurr2,Bkenergy2,Bmass2,Bcharge2,&
         tunex2,tuney2,tunez2,ax2,ay2,bx2,by2,emitx2,emity2,tauy2,&
         dampart2,sigz2,spop2,npcol,nprow,close2,close2ss,&
         ifdbk,isweep2,swrad2,swtune2,iclosq2,nmeet,sigmaz2,sigmapz2,&
         alpha2,phi2,flagrange,flagslice2,flagmap2nd,flagmap3rd,flagmap4th,Ncl,&
         Extrange2,flagGroup,flagwk,qx2,qy2,hcuv2,nfrqlum,tauz2,flagdamp2,&
         nbunch2,nsteps,idrstart,tmax,saveAll,Nresamp,ptRate,ptFrac,momOutRate,gain2,betaCrab2,crabFreq2,&
         crabVnorm2,cfreq22,cvnorm22,np_xy,len_xy,repRate_xy,&
         gainx2,gainy2,frqhighx2,frqhighy2)

    !// Ensure number of macro particles fits number of cores (per beam)
    tmpvalue = modulo(Np2,npcol*nprow/flagGroup)
    if(tmpvalue/=0) then
       Np2 = int(Np2 - tmpvalue)
       if(myid==0) print*, "Warning: The number of macro particles in beam 2 was truncated to", &
            Np2, " for an equal particle number on each core."
    endif
    if(Nturn<len_xy) then
       if(myid==0) print*, "Warning: Requested length of xy series is larger than simulation time. &
            &Shortened length of xy series."
       len_xy=Nturn
    endif
    if(len_xy>repRate_xy) then
       if(myid==0) print*, "Warning: Requested length of xy series is too large for given repetition rate.&
            & Increased repetition rate."
       repRate_xy=len_xy
    endif
    if(np_xy>Np1 .or. np_xy>Np2) then
       if(myid==0) print*, "Warning: Requested number of particles for xy series exceeds number of &
            macro particles. Reduced number of particles to track."
       np_xy = min(Np1,Np2)
    endif
    tmpvalue = np_xy*2/comm_size
    if(tmpvalue<np_xy) then
       np_xy = tmpvalue
       if(myid==0) print*,"Warning: Reduced number particles tracked for xy series to make it divisible &
            without remainder by the number of processors."
    endif

    gx = (1.0d0+ax2*ax2)/bx1
    sigma2(1) = sqrt(emitx2*bx2)
    sigma2(2) = sqrt(emitx2*gx)
    gy = (1.0d0+ay2*ay2)/by2
    sigma2(3) = sqrt(emity2*by2)
    sigma2(4) = sqrt(emity2*gy)
    sigma2(5) = distparam2(15)
    sigma2(6) = distparam2(16)
    if(close2(1).gt.1.0d-4*sigma2(1)) then  !//  check initial offset
       flagoffset = 1
    endif
    if(close2(2).gt.1.0d-4*sigma2(3)) then
       flagoffset = 1
    endif
    ga2 = Bkenergy2/Bmass2 
    bet2 = sqrt(1.0d0-1.0d0/ga2/ga2)

    coef1 = Bcharge1*Bcharge2*qe*Bcurr2*(1.0+bet1*bet2) &
         /(ga1*bet1*(bet1+bet2)*2*pi*epsilon0*Bmass1)
    coef2 = Bcharge1*Bcharge2*qe*Bcurr1*(1.0+bet1*bet2) &
         /(ga2*bet2*(bet1+bet2)*2*pi*epsilon0*Bmass2)
    rree1 = 1.0d0/(4*pi*8.854187817620d-12)*1.6d-19/Bmass1
    rree2 = 1.0d0/(4*pi*8.854187817620d-12)*1.6d-19/Bmass2
    bbparm1x = (Bcurr2*rree1/ga1)*&
         (bx1/(2*pi*sigma2(1)*(sigma2(1)+sigma2(3))))
    bbparm1y = (Bcurr2*rree1/ga1)*&
         (by1/(2*pi*sigma2(3)*(sigma2(1)+sigma2(3))))
    bbparm2x = (Bcurr1*rree2/ga2)*&
         (bx2/(2*pi*sigma1(1)*(sigma1(1)+sigma1(3))))
    bbparm2y = (Bcurr1*rree2/ga2)*&
         (by2/(2*pi*sigma1(3)*(sigma1(1)+sigma1(3))))

    !-------------------------------------------------------------------
    ! construct 2D logical processor Cartesian coordinate
    call construct_Pgrid2d(grid2d,mpicommwd,nprow,npcol)
    call getpost_Pgrid2d(grid2d,myid,myidy,myidx)
    call getcomm_Pgrid2d(grid2d,comm2d,commcol,commrow)

    if(myid.eq.0) then
       print*,"Start simulation:"
       !generate the log file
       open(1,file="BeamBeam3D.log",status="unknown")
       write(1,*)"Simulation Using the BeamBeam3D Code:"
       write(1,*)"Copyright by the Regents of the University of California."
       write(1,*)"The number of processors used in this simulation is ",&
            nprow*npcol
       write(1,*)"The number of turns to be tracked: ",Nturn
       write(1,*)"The frequency of particle output: ",ptFrac
       write(1,*)"The fraction of particle to sampled for output: ",ptRate
       write(1,*)"Moment output rate: ",momOutRate
       write(1,*)"The particle distribution is resampled after (turns): ",Nresamp
       write(1,*)"The flag for strong-strong or weak strong simulation: ",flagwk
       write(1,*)"The group of processors ",flagGroup
       write(1,*)"The maximum simulation time for particle output: ",tmax
       write(1,*)"The flag for restarting simulation from previous stopping point: ",&
            idrstart
       write(1,*)"The number of frequency for luminosity calculation (turn): ",nfrqlum
       write(1,*)"The flag of using feedback control of bunch center: ",ifdbk
       write(1,*)"The number of turn for close-orbit squeeze: ",nmeet
       write(1,*)"The horizontal half crossing angle: ",alpha1, alpha2
       write(1,*)"The vertical half crossing angle: ",phi1, phi2
       write(1,*)"The number of collision points: ",Ncl
       write(1,*)"The number of steps per bunch in one turn: ",nsteps
       write(1,*)"The flag for external linear transfer map: ",flagmap2nd
       write(1,*)"The flag for external 2nd order transfer map: ",flagmap3rd
       write(1,*)"The flag for external 3rd order transfer map: ",flagmap4th
       write(1,*)"-----------------------------------------------------------"
       write(1,*)"Parameters for beam 1:"
       write(1,*)"-----------------------------------------------------------"
       write(1,*)"The number of transverse grids is ",Nx1,Ny1
       write(1,*)"The number of longitudinal slices ",Nslice1
       write(1,*)"The number of macroparticles ",Np1
       write(1,*)"The initial distribution type is ",Flagdist1
       write(1,*)"The parameters of initial distribution: "
       write(1,101)distparam1(1:7)
       write(1,101)distparam1(8:14)
       write(1,101)distparam1(15:21)
       write(1,*)"The real number of particles is: ",Bcurr1
       write(1,*)"The total energy (eV) of the beam: ",Bkenergy1
       write(1,*)"The mass of the particle (eV): ",Bmass1
       write(1,*)"The charge the particle: ",Bcharge1
       write(1,*)"The horizontal tune: ",tunex1
       write(1,*)"The vertical tune: ",tuney1
       write(1,*)"The synchrotron tune: ",tunez1
       write(1,*)"The horizontal beam alpha at IP: ",ax1
       write(1,*)"The vertical beam alpha at IP: ",ay1
       write(1,*)"The horizontal beam beta at IP: ",bx1
       write(1,*)"The vertical beam beta at IP: ",by1
       write(1,*)"The horizontal unnormalized beam emittance at IP: ",emitx1
       write(1,*)"The vertical unnormalized beam emittance at IP: ",emity1
       write(1,*)"The longitudinal rms bunch length  at IP: ",sigmaz1
       write(1,*)"The longitudinal rms relative momentum deviation  at IP: ",sigmapz1
       write(1,*)"The flag for radiation damping, quantum fluctuation calculation: ",&
            flagdamp1
       write(1,*)"The horizontal damping time: ",tauy1
       write(1,*)"The vertical damping time: ",dampart1
       write(1,*)"The longitudinal damping time: ",tauz1
       write(1,*)"The horizontal chromaticity: ",qx1
       write(1,*)"The vertical chromaticity: ",qy1
       write(1,*)"The flag for fixed computational domain or adaptive computational &
            &domain: ",flagrange
       write(1,*)"The range of external computational domain: "
       write(1,102)Extrange1
       write(1,*)"The number of bunches per beam: ",nbunch1
       write(1,*)"The curvature of the ring: ",hcuv1
       write(1,*)"The flag of using readin number of slice: ",flagslice1
       write(1,*)"The flag of doing close-orbit squeeze: ",iclosq1
       write(1,*)"The final close-orbit coordinates: ",close1ss
       write(1,*)"The flag of doing close-orbit sweeping: ",isweep1
       write(1,*)"The tune of close-orbit sweeping: ",swtune1
       write(1,*)"The amplitude of close-orbit sweeping: ",swrad1
       write(1,*)"Gains of the transverse feedback system (damper): ",gain1
       write(1,*)"beta function at crab cavity and crab cavity frequency: ",betaCrab1,crabFreq1
       write(1,*)"particle number, length and repetition rate of output in xy data files: ",np_xy,len_xy,repRate_xy
       write(1,*)"-----------------------------------------------------------"
       write(1,*)"Parameters for beam 2:"
       write(1,*)"-----------------------------------------------------------"
       write(1,*)"The number of transverse grids is ",Nx2,Ny2
       write(1,*)"The number of longitudinal slices ",Nslice2
       write(1,*)"The number of macroparticles ",Np2
       write(1,*)"The initial distribution type is ",Flagdist2
       write(1,*)"The parameters of initial distribution: "
       write(1,101)distparam2(1:7)
       write(1,101)distparam2(8:14)
       write(1,101)distparam2(15:21)
       write(1,*)"The real number of particles is: ",Bcurr2
       write(1,*)"The total energy (eV) of the beam: ",Bkenergy2
       write(1,*)"The mass of the particle (eV): ",Bmass2
       write(1,*)"The charge the particle: ",Bcharge2
       write(1,*)"The horizontal tune: ",tunex2
       write(1,*)"The vertical tune: ",tuney2
       write(1,*)"The synchrotron tune: ",tunez2
       write(1,*)"The horizontal beam alpha at IP: ",ax2
       write(1,*)"The vertical beam alpha at IP: ",ay2
       write(1,*)"The horizontal beam beta at IP: ",bx2
       write(1,*)"The vertical beam beta at IP: ",by2
       write(1,*)"The horizontal unnormalized beam emittance at IP: ",emitx2
       write(1,*)"The vertical unnormalized beam emittance at IP: ",emity2
       write(1,*)"The longitudinal rms bunch length  at IP: ",sigmaz2
       write(1,*)"The longitudinal rms relative momentum deviation  at IP: ",sigmapz2
       write(1,*)"The flag for radiation damping, quantum fluctuation calculation: ", &
            flagdamp2
       write(1,*)"The horizontal damping time: ",tauy2
       write(1,*)"The vertical damping time: ",dampart2
       write(1,*)"The longitudinal damping time: ",tauz2
       write(1,*)"The horizontal chromaticity: ",qx2
       write(1,*)"The vertical chromaticity: ",qy2
       write(1,*)"The flag for fixed computational domain or adaptive computational &
            &domain: ",flagrange
       write(1,*)"The range of external computational domain: "
       write(1,102)Extrange2
       write(1,*)"The number of bunches per beam: ",nbunch2
       write(1,*)"The curvature of the ring: ",hcuv2
       write(1,*)"The flag of using readin number of slice: ",flagslice2
       write(1,*)"The flag of doing close-orbit squeeze: ",iclosq2
       write(1,*)"The final close-orbit coordinates: ",close2ss
       write(1,*)"The flag of doing close-orbit sweeping: ",isweep2
       write(1,*)"The tune of close-orbit sweeping: ",swtune2
       write(1,*)"The amplitude of close-orbit sweeping: ",swrad2
       write(1,*)"Gains of the transverse feedback system (damper): ", gain2
       write(1,*)"beta function at crab cavity and crab cavity frequency: ",betaCrab2,crabFreq2
       write(1,*)"particle number, length and repetition rate of output in xy data files: ",np_xy,len_xy,repRate_xy
       write(1,*)"----------------------------------------------------------"
       write(1,*)"The beam-beam parameters per IP are (beam 1 x,y, beam 2 x,y):"
       write(1,103)bbparm1x,bbparm1y,bbparm2x,bbparm2y
       close(1)
    endif
101 format(7(1x,es15.7))
102 format(6(1x,es15.7))
103 format(4(1x,es15.7))

    !-------------------------------------------------------------------
    if(flagGroup.eq.2) then !strong-strong simulation using 2 group PEs.
       ngroup = 2
       if(mod(npcol,2).ne.0) then
          print*,"npcol has to be even number"
          stop
       endif
       npyhalf = npcol/2
       if(myidy.lt.npyhalf) then
          Nplocal = Np1/(npyhalf*nprow)
          distparam = distparam1
          Flagdist = Flagdist1
          Np = Np1
          NpO = Np2
          tunex = tunex1
          tuney = tuney1
          tunez = tunez1
          ax = ax1
          bx = bx1
          ay = ay1
          by = by1
          phi = phi1
          alpha = alpha1
          tauy = tauy1
          dampart = dampart1
          tauz = tauz1
          flagdamp = flagdamp1
          coef = coef1
          swrad = swrad1
          swtune = swtune1
          isweep = isweep1
          close2g = close1
          iclosq = iclosq1
          closess = close1ss
          sigma = sigma1
          sigopp = (sigma2(1)+sigma2(3))/2
          Extrange = Extrange1
          qx = qx1
          qy = qy1
          hcuv = hcuv1
          Bcurr = Bcurr1
          nbunch = nbunch1
          betaCrab = betaCrab1
          crabFreq = crabFreq1 
          crabVnorm = crabVnorm1
          cfreq2 = cfreq12
          cvnorm2 = cvnorm12
          frqhighx = frqhighx1
          frqhighy = frqhighy1
          gainx = gainx1
          gainy = gainy1
       else
          Nplocal = Np2/(npyhalf*nprow)
          distparam = distparam2
          Flagdist = Flagdist2
          Np = Np2
          NpO = Np1
          tunex = tunex2
          tuney = tuney2
          tunez = tunez2
          ax = ax2
          bx = bx2
          ay = ay2
          by = by2
          phi = phi2
          alpha = alpha2
          tauy = tauy2
          dampart = dampart2
          tauz = tauz2
          flagdamp = flagdamp2
          coef = coef2
          swrad = swrad2
          swtune = swtune2
          isweep = isweep2
          close2g = close2
          iclosq = iclosq2
          closess = close2ss
          sigma = sigma2
          sigopp = (sigma1(1)+sigma1(3))/2
          Extrange = Extrange2
          qx = qx2
          qy = qy2
          hcuv = hcuv2
          Bcurr = Bcurr2
          nbunch = nbunch2
          betaCrab = betaCrab2
          crabFreq = crabFreq2 
          crabVnorm = crabVnorm2
          cfreq2 = cfreq22
          cvnorm2 = cvnorm22
          frqhighx = frqhighx2
          frqhighy = frqhighy2
          gainx = gainx2
          gainy = gainy2
       endif

       ! crab cavity parameters
       if(crabfreq < 1.0d-8) then  ! caclulate crab voltage
          crabVnorm = 0.d0  ! cavity off
       endif

       !assign the charge, current to each bunch
       !The following needs to be modified for more robust cases.
       allocate(Bcurrlist1(nbunch1))
       allocate(Bmasslist1(nbunch1))
       allocate(Bchargelist1(nbunch1))
       allocate(Bcurrlist2(nbunch2))
       allocate(Bmasslist2(nbunch2))
       allocate(Bchargelist2(nbunch2))
       !// Here, we have assumed that all bunches have the same mass, charge and number of macro particles
       Bmasslist1 = Bmass1
       Bchargelist1 = Bcharge1
       Bmasslist2 = Bmass2
       Bchargelist2 = Bcharge2
       do i = 1, nbunch1
          Bcurrlist1(i) = Bcurr1
       enddo
       do i = 1, nbunch2
          Bcurrlist2(i) = Bcurr2
       enddo

       !Here, we have assumed that all bunches have the same number of particles.
       allocate(Nplist(nbunch))
       allocate(Nplclist(nbunch))
       Nplist = Np
       Nplclist = Nplocal
       allocate(BptsArry(6,Nplocal,nbunch))
       allocate(Bpts(6,Nplocal))
       allocate(Bptstmp(6,Nplocal))
       !generate initial distribution for each bunch from restart
       iendleft = 0
       iticleft = 0
       iseedend = -10-27*myid
       icount = 0
       if(idrstart.eq.1) then
          !read in the distribution from the last simulation 
          do i = 1, nbunch
             call phaseinBB_Output(Bpts,myid,iendleft,i,Nplocal,closetmp,iticleft,icount)
             close2g = closetmp
             BptsArry(:,:,i) = Bpts(:,:)
          enddo

          !initialization of random number generator
          call ran2(iseedend,closetmp,2)
          do i = 1, int(icount)
             call ran2(iseedend,tmp1,6)
          enddo
       else
          !initialization of random number generator
          call ran2(iseedend,closetmp,2)
          if(myidy.lt.npyhalf) then  !// initialize Distribution:sobseq
             call sobseq(-2,Bpts(1:6,1),303)  !// Bpts is just a dummy here
          else
             call sobseq(-2,Bpts(1:6,1),3030)
          endif
          !generate the initial distribution from given parameters
          do i = 1, nbunch
             call sample_Dist(Bpts,distparam,21,Flagdist,grid2d,Np,Nplocal,ngroup,i)
             BptsArry(:,:,i) = Bpts(:,:)
          enddo
       endif
    else if(flagGroup.eq.1) then !weak-strong simulation using 1 group PEs.
       ngroup = 1
       Nplocal1 = Np1/(npcol*nprow)
       allocate(Bpts1(6,Nplocal1))
       !this Bptstmp is not used
       allocate(Bptstmp(6,Nplocal1))
       if(Flagdist1.ne.11) then
          call sample_Dist(Bpts1,distparam1,21,Flagdist1,grid2d,Np1,Nplocal1,&
               ngroup,ngroup)
       else
          call read1g_Dist(Bpts1,21,distparam1,grid2d,&
               Np1,Nplocal1,ngroup,"partcl1.data")

       endif
       Nplocal2 = Np2/(npcol*nprow)
       allocate(Bpts2(6,Nplocal2))
       if(Flagdist1.ne.11) then
          call sample_Dist(Bpts2,distparam2,21,Flagdist2,grid2d,Np2,Nplocal2,&
               ngroup,ngroup)
       else
          call read1g_Dist(Bpts2,21,distparam2,grid2d,&
               Np2,Nplocal2,ngroup,"partcl2.data")
       endif
    endif

    if(myid.eq.0) print*,"finish sampling..."
    !construct the computation domain using the maximum
    !rms size of all bunches.
    s1 = max(distparam1(1),distparam2(1))
    s2 = max(distparam1(8),distparam2(8))
    s3 = max(distparam1(15),distparam2(15))

    distparam1(1) = s1
    distparam1(8) = s2
    distparam1(15) = s3
    xrad = 0.0
    yrad = 0.0
    Perdlen = 0.0

    !-------------------------------------------------------------------
    ! construct computational domain CompDom class and get local geometry 
    ! information on each processor.
    if( (Nx1.ne.Nx2).or.(Ny1.ne.Ny2) ) then
       if(myid.eq.0) then
          print*,"In this version beam 1 and beam 2 should have the"
          print*,"same computational domain!!!"
       endif
       stop
    endif
    Flagbc = 1
    if(myid.eq.0) print*,"before construct the domain...."
    !here we have assumed that Nx1=Nx2,Ny1=Ny2,Nz1=Nz2. Hence,
    !we have only one computational domain instead of two
    !computational domain for each beam.
    if(flagGroup.eq.2) then
       call construct_CompDom(Ageom,distparam1,21,Flagdist1,&
            Nx1,Ny1,Nz1,grid2d,nprow,npcol,Flagbc,xrad,yrad,Perdlen)
    else if(flagGroup.eq.1) then
       call construct1G_CompDom(Ageom,distparam1,21,Flagdist1,&
            Nx1,Ny1,Nz1,grid2d,nprow,npcol,Flagbc,xrad,yrad,Perdlen)
    endif
    if(myid.eq.0) print*,"finish construct the domain...."

    !------------------------------------------------------------------
    if(flagGroup.eq.2) then !using two group PEs
       !//construct linear lattice map using alpha, beta, tune or 
       !//read in from MAD output
       if(flagmap2nd.ne.1) then
          call construct_Linearmap(Alinearmap,tunex/dble(Ncl),tuney/dble(Ncl),&
               tunez/dble(Ncl),ax,bx,ay,by)
       else
          allocate(Amap2nd(Ncl))
          allocate(Amap2ndArry(Ncl,nbunch))
          if(myidy.lt.npyhalf) then
             open(1,file="map2nd1.in",status="old")
             !loop through the collision points
             do i = 1, Ncl
                call construct_Extmap2nd(Amap2nd(i),tunex/dble(Ncl),tuney/dble(Ncl),&
                     tunez/dble(Ncl),ax,bx,ay,by)
                do j = 1, 6
                   read(1,*)tmpvalue1,tmpvalue2,tmpvalue3,tmpvalue4,tmpvalue5,tmpvalue6
                   tmpmap2nd(j,1) = tmpvalue1
                   tmpmap2nd(j,2) = tmpvalue2
                   tmpmap2nd(j,3) = tmpvalue3
                   tmpmap2nd(j,4) = tmpvalue4
                   tmpmap2nd(j,5) = tmpvalue5
                   tmpmap2nd(j,6) = tmpvalue6
                enddo

                !                  call symplect(tmpmap2nd,6,error)
                call setmap_Extmap2nd(Amap2nd(i),tmpmap2nd)
             enddo
             close(1)
          else
             open(2,file="map2nd2.in",status="old")
             do i = 1, Ncl
                call construct_Extmap2nd(Amap2nd(i),tunex/dble(Ncl),tuney/dble(Ncl),&
                     tunez/dble(Ncl),ax,bx,ay,by)
                do j = 1, 6
                   read(2,*)tmpvalue1,tmpvalue2,tmpvalue3,tmpvalue4,tmpvalue5,tmpvalue6
                   tmpmap2nd(j,1) = tmpvalue1
                   tmpmap2nd(j,2) = tmpvalue2
                   tmpmap2nd(j,3) = tmpvalue3
                   tmpmap2nd(j,4) = tmpvalue4
                   tmpmap2nd(j,5) = tmpvalue5
                   tmpmap2nd(j,6) = tmpvalue6
                enddo

                !                  call symplect(tmpmap2nd,6,error)
                call setmap_Extmap2nd(Amap2nd(i),tmpmap2nd)
             enddo
             close(2)
          endif
       endif
       if(flagmap3rd.eq.1) then
          !//construct nonlinear lattice map.
          allocate(Amap3rd(Ncl))
          allocate(Amap3rdArry(Ncl,nbunch))
          if(myidy.lt.npyhalf) then
             open(3,file="map3rd1.in",status="old")
             do i = 1, Ncl
                call construct_Extmap3rd(Amap3rd(i))
                do j = 1, 216
                   ii = (j-1)/36 + 1
                   jj = (j - (ii-1)*36 - 1)/6 + 1
                   kk = j - (ii-1)*36 - (jj-1)*6
                   read(3,*)tmpvalue
                   tmpmap3rd(ii,jj,kk) = tmpvalue
                enddo
                call setmap_Extmap3rd(Amap3rd(i),tmpmap3rd)
             enddo
             close(3)
          else
             open(4,file="map3rd2.in",status="old")
             do i = 1, Ncl
                call construct_Extmap3rd(Amap3rd(i))
                do j = 1, 216
                   ii = (j-1)/36 + 1
                   jj = (j - (ii-1)*36 - 1)/6 + 1
                   kk = j - (ii-1)*36 - (jj-1)*6
                   read(4,*)tmpvalue
                   tmpmap3rd(ii,jj,kk) = tmpvalue
                enddo
                call setmap_Extmap3rd(Amap3rd(i),tmpmap3rd)
             enddo
             close(4)
          endif
       endif
       if(flagmap4th.eq.1) then
          !//construct nonlinear lattice map.
          allocate(Amap4th(Ncl))
          allocate(Amap4thArry(Ncl,nbunch))
          if(myidy.lt.npyhalf) then
             open(3,file="map4th1.in",status="old")
             do i = 1, Ncl
                call construct_Extmap4th(Amap4th(i))
                do j = 1, 1296
                   ii0 = (j-1)/216 + 1
                   ii = (j-(ii0-1)*216-1)/36 + 1
                   jj = (j - (ii0-1)*216 - (ii-1)*36 - 1)/6 + 1
                   kk = j - (ii0-1)*216 - (ii-1)*36 - (jj-1)*6
                   read(3,*)tmpvalue
                   tmpmap4th(ii0,ii,jj,kk) = tmpvalue
                enddo
                call setmap_Extmap4th(Amap4th(i),tmpmap4th)
             enddo
             close(3)
          else
             open(4,file="map4th2.in",status="old")
             do i = 1, Ncl
                call construct_Extmap4th(Amap4th(i))
                do j = 1, 1296
                   ii0 = (j-1)/216 + 1
                   ii = (j-(ii0-1)*216-1)/36 + 1
                   jj = (j - (ii0-1)*216 - (ii-1)*36 - 1)/6 + 1
                   kk = j - (ii0-1)*216 - (ii-1)*36 - (jj-1)*6
                   read(4,*)tmpvalue
                   tmpmap4th(ii0,ii,jj,kk) = tmpvalue
                enddo
                call setmap_Extmap4th(Amap4th(i),tmpmap4th)
             enddo
             close(4)
          endif
       endif

       !//read in the bunch id at each collision point.
       if(Ncl>1 .or. nsteps>1) then
          if(myidy.lt.npyhalf) then
             open(3,file="idblist1.in",status="old")
             do i = 1, Ncl
                read(3,*)idblist(i,1:nsteps)
             enddo
             close(3)
             open(4,file="idblist2.in",status="old")
             do i = 1, Ncl
                read(4,*)idblistO(i,1:nsteps)
             enddo
             close(4)
          else
             open(4,file="idblist2.in",status="old")
             do i = 1, Ncl
                read(4,*)idblist(i,1:nsteps)
             enddo
             close(4)
             open(3,file="idblist1.in",status="old")
             do i = 1, Ncl
                read(3,*)idblistO(i,1:nsteps)
             enddo
             close(3)
          endif
       else
          idblist = 1
          idblistO = 1
       endif

    else if(flagGroup.eq.1) then !use 1 group PE
       ! construct linear lattice map
       if(myid.eq.0) print*,"strong-weak beam-beam....",flagGroup,flagwk
       if(flagmap2nd.ne.1) then
          call construct_Linearmap(Alinearmap1,tunex1/dble(Ncl),tuney1/dble(Ncl),&
               tunez1/dble(Ncl),ax1,bx1,ay1,by1)
          call construct_Linearmap(Alinearmap2,tunex2/dble(Ncl),tuney2/dble(Ncl),&
               tunez2/dble(Ncl),ax2,bx2,ay2,by2)
       else
          allocate(Amap2nd1(Ncl))
          if(myid.eq.0) print*,"start construct map2dn1...",Ncl
          open(1,file="map2nd1.in",status="old")
          do i = 1, Ncl
             call construct_Extmap2nd(Amap2nd1(i),tunex1/dble(Ncl),tuney1/dble(Ncl),&
                  tunez1/dble(Ncl),ax1,bx1,ay1,by1)
             do j = 1, 6
                read(1,*)tmpvalue1,tmpvalue2,tmpvalue3,tmpvalue4,tmpvalue5,tmpvalue6
                tmpmap2nd(j,1) = tmpvalue1
                tmpmap2nd(j,2) = tmpvalue2
                tmpmap2nd(j,3) = tmpvalue3
                tmpmap2nd(j,4) = tmpvalue4
                tmpmap2nd(j,5) = tmpvalue5
                tmpmap2nd(j,6) = tmpvalue6
             enddo

             call symplect(tmpmap2nd,6,error)
             call setmap_Extmap2nd(Amap2nd1(i),tmpmap2nd)
          enddo
          close(1)
          if(myid.eq.0) print*,"pass construct map2dn1..."
          allocate(Amap2nd2(Ncl))
          open(2,file="map2nd2.in",status="old")
          do i = 1, Ncl
             call construct_Extmap2nd(Amap2nd2(i),tunex2/dble(Ncl),tuney2/dble(Ncl),&
                  tunez2/dble(Ncl),ax2,bx2,ay2,by2)
             do j = 1, 6
                read(2,*)tmpvalue1,tmpvalue2,tmpvalue3,tmpvalue4,tmpvalue5,tmpvalue6
                tmpmap2nd(j,1) = tmpvalue1
                tmpmap2nd(j,2) = tmpvalue2
                tmpmap2nd(j,3) = tmpvalue3
                tmpmap2nd(j,4) = tmpvalue4
                tmpmap2nd(j,5) = tmpvalue5
                tmpmap2nd(j,6) = tmpvalue6
             enddo

             call symplect(tmpmap2nd,6,error)
             call setmap_Extmap2nd(Amap2nd2(i),tmpmap2nd)
          enddo
          close(2)
       endif
       if(flagmap3rd.eq.1) then
          !//construct nonlinear lattice map.
          !//allocate(Amap3rd1(Ncl+1))
          allocate(Amap3rd1(Ncl))
          open(3,file="map3rd1.in",status="old")
          do i = 1, Ncl
             call construct_Extmap3rd(Amap3rd1(i))
             do j = 1, 216
                ii = (j-1)/36 + 1
                jj = (j - (ii-1)*36 - 1)/6 + 1 
                kk = j - (ii-1)*36 - (jj-1)*6
                read(3,*)tmpvalue
                tmpmap3rd(ii,jj,kk) = tmpvalue
             enddo
             call setmap_Extmap3rd(Amap3rd1(i),tmpmap3rd)
          enddo
          close(3)
          !//allocate(Amap3rd2(Ncl+1))
          allocate(Amap3rd2(Ncl))
          open(4,file="map3rd2.in",status="old")
          do i = 1, Ncl
             call construct_Extmap3rd(Amap3rd2(i))
             do j = 1, 216
                ii = (j-1)/36 + 1
                jj = (j - (ii-1)*36 - 1)/6 + 1 
                kk = j - (ii-1)*36 - (jj-1)*6
                read(4,*)tmpvalue
                tmpmap3rd(ii,jj,kk) = tmpvalue
             enddo
             call setmap_Extmap3rd(Amap3rd2(i),tmpmap3rd)
          enddo
          close(4)
       endif
       if(flagmap4th.eq.1) then
          !//construct nonlinear lattice map.
          !//allocate(Amap4th1(Ncl+1))
          allocate(Amap4th1(Ncl))
          open(3,file="map4th1.in",status="old")
          do i = 1, Ncl
             call construct_Extmap4th(Amap4th1(i))
             do j = 1, 1296
                ii0 = (j-1)/216 + 1
                ii = (j-(ii0-1)*216-1)/36 + 1
                jj = (j - (ii0-1)*216 - (ii-1)*36 - 1)/6 + 1
                kk = j - (ii0-1)*216 - (ii-1)*36 - (jj-1)*6
                read(3,*)tmpvalue
                tmpmap4th(ii0,ii,jj,kk) = tmpvalue
             enddo
             call setmap_Extmap4th(Amap4th1(i),tmpmap4th)
          enddo
          close(3)
          !//allocate(Amap4th2(Ncl+1))
          allocate(Amap4th2(Ncl))
          open(4,file="map4th2.in",status="old")
          do i = 1, Ncl
             call construct_Extmap4th(Amap4th2(i))
             do j = 1, 1296
                ii0 = (j-1)/216 + 1
                ii = (j-(ii0-1)*216-1)/36 + 1
                jj = (j - (ii0-1)*216 - (ii-1)*36 - 1)/6 + 1
                kk = j - (ii0-1)*216 - (ii-1)*36 - (jj-1)*6
                read(3,*)tmpvalue
                tmpmap4th(ii0,ii,jj,kk) = tmpvalue
             enddo
             call setmap_Extmap4th(Amap4th2(i),tmpmap4th)
          enddo
          close(4)
       endif
    endif

    !-------------------------------------------------------------------
    ! read in number of slices at each collision point
    if(flagslice1.eq.1) then
       open(4,file="Nslice1.in",status="old")
       do i = 1, Ncl
          read(4,*)Nslice1list(i)
       enddo
       close(4)
    else !//all collision points have the same # of slices
       Nslice1list = Nslice1
    endif
    if(flagslice2.eq.1) then
       open(8,file="Nslice2.in",status="old")
       do i = 1, Ncl
          read(8,*)Nslice2list(i)
       enddo
       close(8)
    else
       Nslice2list = Nslice2
    endif

    !-------------------------------------------------------------------
    ! construct radiation damping and quantumn fluctuation map
    if(flagGroup.eq.2) then
       if(flagdamp.eq.1) then
          call construct_Dampfluc(Adampfluc,tauy,dampart) 
       else if(flagdamp.eq.2) then
          !in this case, tauy is actually taux,dampart is tauy
          call construct_Dampfluc(Adampfluc,tauy,dampart,tauz) 
       endif
    else if(flagGroup.eq.1) then
       if(flagdamp1.eq.1) then
          call construct_Dampfluc(Adampfluc1,tauy1,dampart1) 
       else if(flagdamp1.eq.2) then
          !in this case, tauy is actually taux,dampart is tauy
          call construct_Dampfluc(Adampfluc1,tauy1,dampart1,tauz1) 
       endif
       if(flagdamp2.eq.1) then
          call construct_Dampfluc(Adampfluc2,tauy2,dampart2) 
       else if(flagdamp2.eq.2) then
          !in this case, tauy is actually taux,dampart is tauy
          call construct_Dampfluc(Adampfluc2,tauy2,dampart2,tauz2) 
       endif
    endif

    !-------------------------------------------------------------------
    call MPI_BARRIER(mpicommwd,ierr)
    if(myid.eq.0) print*,"finish initalization....."

    t_init = t_init + elapsedtime_Timer(t0)

  end subroutine init_AccSimulator

  !//Run a beam-beam simulation through accelerator.
  subroutine run_AccSimulator()
    implicit none
    !include 'mpif.h'
    integer :: i,n,innx,inny
    double precision :: zmin1,zmin2,zmax1,zmax2,z
    double precision :: t0,eps
    double precision, dimension(2) :: tmp1,tmp2,tmplc1,tmplc2,close1init,close2init,&
         delta1,delta2,delta,closeold1,closeold2,shift,closeold,closeinit
    double precision :: factor,szspz
    double precision, dimension(3) :: center,centrl1,centrl2
    !only for CERN benchmark
    double precision, dimension(6) :: centerCern
    integer :: ierr,maxpts1,maxpts2,ipt
    integer :: nxlum,nylum,ic,npyhalf,jendmax,nz,NsliceO,Nslice
    double precision :: zmin,zmax
    !//for sweeping and squeeze only
    double precision :: ang0
    double precision, dimension(6) :: sigtmp1,sigtmp2
    double precision :: betas1,betas2,ga2,ga1,betas
    double precision :: xmin,xmax,ymin,ymax
    !//for fixed external domain Green function.
    double precision, dimension(3) :: msize
    double precision, dimension(4) :: range1,range2
    double precision, dimension(6) :: tmprange
    double precision  :: hx,hy
    double complex, allocatable, dimension (:,:) :: grn
    integer, allocatable, dimension(:,:,:) :: table
    integer, allocatable, dimension(:) :: xpystable,pytable
    integer :: nxtot,nytot,nsxy1,nsxy2,nxpylc2,npbc,flaglum
    !//luminosity 3d
    double precision :: lum3d,Bcurr0,Bcurr01,Bcurr02,pi
    integer :: itic,Np0,Np01,Np02,idbunch,i6,istp,idbunchO
    double precision :: time0,time1,timepass
    integer :: iend,iseedinit !,ntmp1,ntmp2
    double precision :: gambet
    double precision :: timepasslc,phi_s
    integer :: errSeed,j,ptOut
    type(Extmap2nd) :: tmpMap2nd2
    integer, allocatable, dimension(:) :: ptSampleInd
    integer, allocatable, target, dimension(:) :: QsampleInd
    integer, pointer, dimension(:) :: ptr => NULL()
    integer :: comm_size,lumFile
    double precision :: T_0
    double precision, dimension(3) :: offset
    integer :: iseedglb
    real*8 :: xerr,crabvnorm0,xerr2
    real*8, dimension(2) :: rry
    real*8 :: a1,a2,a3,a4,a5,a6,a7
    real*8 :: b2,b3,b4,b5,b6,b7,eng0

    call starttime_Timer(t0)
    call MPI_BARRIER(mpicommwd,ierr)
    time0 = MPI_WTIME()
    iseedinit = iseedend

    call MPI_BARRIER(mpicommwd,ierr)

    if(ptFrac/=0) then
       ptOut = nplocal/ptFrac
    else
       ptOut = 0
    endif

    lumFile = 11

    !// Initialize crab cavity noise variables
    T_0 =  26.659d3/2.99792458d8  ! revolution period
    offset = [0.d0, 0.d0, 0.d0]  ! initial offset

    pi = 2*asin(1.d0)
    !ga1 = Bkenergy1/Bmass1 + 1.d0
    ga1 = Bkenergy1/Bmass1 
    betas1 = sqrt(1.d0-1.d0/ga1/ga1)
    !ga2 = Bkenergy2/Bmass2 + 1.d0
    ga2 = Bkenergy2/Bmass2 
    betas2 = sqrt(1.d0-1.d0/ga2/ga2)

    xmin = Extrange(1)
    xmax = Extrange(2)
    ymin = Extrange(3)
    ymax = Extrange(4)

    !//Used for the oscillation of beam 2.
    swrad2 = sqrt(close2(1)**2+close2(2)**2)
    if(close2(1).gt.0) then
       ang0 = asin(close2(2)/swrad2)
    else
       ang0 = 2*asin(1.0d0) - asin(close2(2)/swrad2)
    endif

    nxlum = 256 !//Grid number of x for luminosity calculation
    nylum = 256 !//Grid number of y for luminosity calculation
    lum3d = 0.0

    eps = 1.0d-4 !//extra boundary
    factor = 8.0d0 !//scaling factor for the maximum # of particles per slice (tolerance for unequal distribution)
    if(idrstart.eq.1) then
       iend = iendleft
    else
       iend = 0
    endif
    z = iend
    closeold = close2g

    crabvnorm0 = crabVnorm

    if(flagGroup.eq.2) then !//Strong-strong using 2 group!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
       !//Strong-strong beam-beam interaction below.
       !//To do a strong-strong simulation using only 2 bunches but multiple
       !//collision points works correctly only for 2 collision points.
       if(myid.eq.0) print*,"strong-strong beam-beam...."
       Np0 = Np !//initial total # of ptcls.
       Bcurr0 = Bcurr !//initial current
       closeinit = close2g
       !//Get local grid number in each computaional domain
       npyhalf = npcol/2
       innx = Nx1
       inny = Ny1/npyhalf

       !initialize the Gaussian random generator
       if(myidy.lt.npyhalf) then 
          errSeed = -1000 
          iseedglb = -47
       else 
          iseedglb = -477
          errSeed = -10000 
       endif

       if(myidy.lt.npyhalf) then
          nz = Nz1
          betas = betas1
          gambet = ga1*betas1
          eng0 = Bkenergy1 
       else
          nz = Nz2
          betas = betas2
          gambet = ga2*betas2
          eng0 = Bkenergy2 
       endif
       sigtmp1 = sigma

       !// Configuration of samples for output
       allocate(ptSampleInd(ptOut))
       call randomIndexList(Nplocal,ptOut,ptSampleInd)
       allocate(QsampleInd(np_xy))
       call randomIndexList(Nplocal,np_xy,QsampleInd)
       ptr => QsampleInd

       phi_s = 0  !// phase shift of crab cavities

       do idbunch = 1, nbunch
          Bpts(:,:) = BptsArry(:,:,idbunch)

          !// Initial output
          if(ptOut>0) call ptclout_Output(Bpts,myid,0,idbunch,Nplocal,ptOut,nturn,ptRate,ptSampleInd)
          if(np_xy>0) call xySeries_Output(Bpts,myid,0,idbunch,Nplocal,np_xy,repRate_xy,&
               len_xy,Nturn,ptr)
          if(momOutRate /= 0) call diagnostic2GMBeic_Output(z,Bpts,Nplocal,Np,myidx,myidy,nprow,npcol,&
               commrow,commcol,comm2d,idbunch,ax,bx,ay,by)

          !// Crab the bunch
          call crabTransform(Bpts,Nplocal,myidy,npyhalf,0,betaCrab,crabFreq,crabVnorm,&
               phi_s,bx,by,close2g,alpha,cfreq2,cvnorm2) 
!---

          BptsArry(:,:,idbunch)=Bpts(:,:)
       enddo
       call MPI_BARRIER(mpicommwd,ierr)

       !//Output initial luminosity
       if(Nslice1.eq.1 .and. Nslice2.eq.1) then
          call luminosity2G_Output(Bpts,Nplocal,nxlum,nylum,itic,Bcurr,Np,myidy,npyhalf)
       endif

       !//find the Green function for the fixed computation domain.
       !//The shift between two beam has to be 0 or fixed.
       if(flagrange.eq.1) then
          shift = 0.0d0

          if(myidy.lt.npyhalf) then
             npbc = 0
          else
             npbc = npyhalf
          endif

          nxtot = Nx1
          nytot = Ny1
          !// +1 is from the real to complex fft.
          nsxy1 = (nxtot+1)/npyhalf
          nsxy2 = (nxtot+1) - npyhalf*nsxy1
          allocate(xpystable(0:npyhalf-1))
          do i = 0, npyhalf-1
             if(i.le.(nsxy2-1)) then
                xpystable(i) = nsxy1+1
             else
                xpystable(i) = nsxy1
             endif
          enddo
          nxpylc2 = xpystable(myidy-npbc)

          allocate(grn(2*nytot,nxpylc2))

          range1(1) = Extrange(1)
          range1(2) = Extrange(2)
          range1(3) = Extrange(3)
          range1(4) = Extrange(4)
          call update2d_CompDom(Ageom,range1,grid2d)
          call getmsize_CompDom(Ageom,msize)
          hx = msize(1)
          hy = msize(2)
          allocate(table(2,0:nprow-1,0:npyhalf-1))
          call getlctabnm_CompDom(Ageom,table)
          allocate(pytable(0:npyhalf-1))
          do i = 0, npyhalf-1
             pytable(i) = table(2,0,i)
          enddo

          !//get the Green function for the fixed domain.
          !//here, you can not have the domain shift
          call greenf2d(nxtot,nytot,nxpylc2,hx,hy,myidy,npyhalf,commcol,&
               xpystable,grn,shift)

          call getrange_CompDom(Ageom,tmprange)
          if(myidy.lt.npyhalf) then
             xmin = tmprange(1)
             ymin = tmprange(3)
          else
             xmin = tmprange(1) + shift(1)
             ymin = tmprange(3) + shift(2)
          endif
       endif !//finish the Green function for the fixed domain.

       itic = iticleft
       call MPI_BARRIER(mpicommwd,ierr)
       time0 = MPI_WTIME() 

       !// Open file for luminosity data
       if(Nslice1list(1)>1)then
          if(idrstart==0) then
             open(lumFile,file="luminosity.data",status="replace",iostat=ierr)
          else
             open(lumFile,file="luminosity.data",status="unknown",position="append",iostat=ierr)
          end if
          if(ierr/=0) print*,"Warning: 'luminosity.data' could not be opened."
       end if

!-----------------------------------------------------------------------------------
       !// Loop through "Nturn" simulation.
!-----------------------------------------------------------------------------------
       do i = iend + 1, Nturn
          !resample the particle distribution every Nresamp turns
          if(mod(i,Nresamp).eq.0) then
             do idbunch = 1, nbunch
                !generate the initial distribution from given parameters
                call sample_Dist(Bpts,distparam,21,Flagdist,grid2d,Np,Nplocal,&
                     flagGroup,idbunch)
                !replace the beam 1 by the new sampled particles
                if(myidy.lt.npyhalf) then
                   BptsArry(:,:,idbunch) = Bpts(:,:)
                endif
             enddo
          endif

          !// Loop through nsteps and Ncl collision points for each turn. After nsteps,
          !/  a bunch returns to its original location.
          do istp = 1, nsteps
             !// Loop through "Ncl" collisions for each step.
             do ic = 1, Ncl
                !get the bunch 1ds that are involved in the collision.
                idbunch = idblist(ic,istp)
                !get the bunch ID of the other beam
                idbunchO = idblistO(ic,istp)

                if(idbunch.le.0) goto 500
                itic = itic + 1
                !//We need to find the shift between the slice at each
                !//collision point: center2 - center1
                Nslice1 = Nslice1list(ic)
                Nslice2 = Nslice2list(ic)
                jendmax = ((min(Nslice1,Nslice2)-1)/nprow + 1)*nprow
                if(myidy.lt.npyhalf) then
                   Nslice = Nslice1list(ic) !//local slice number
                   NsliceO = Nslice2list(ic)!//slice number of the opposite beam
                else
                   Nslice = Nslice2list(ic)
                   NsliceO = Nslice1list(ic)
                endif
                Np = Nplist(idbunch)
                maxpts1 = int(Np/Nslice*factor) !//maximum # of ptl. per slice

                Nplocal = Nplclist(idbunch)
                Bpts(:,:) = BptsArry(:,:,idbunch)

                if(flagoffset.eq.0) then
                   shift = 0.0d0
                endif

!-----------------------------------------------
                !//from Lab coordinate to boost coordinate
                !//when dealing with crossing angle collision.
                !//the particle transverse size is checked to
                !//decide the lost particle.
                !//phi is the half crossing angle in x-z, alpha is the cross angle in x-y
                call transfer1old2(Bpts,alpha,phi,Nplocal)

                !//find zmin and zmax of each beam
                !//The following can be simplified using 2 communicators
                !//The following can be replaced by fixed maximum range
                tmplc1 = 0.0
                tmplc2 = 0.0
                if(Nslice.gt.1) then
                   !if(flagrange.eq.1) then !// fixed external range
                   !  zmin = Extrange(5)
                   !  zmax = Extrange(6)
                   !else !// this could be improved by using MPI_ALLREC
                   if(myidy.lt.npyhalf) then
                      do n = 2, Nplocal
                         if(tmplc1(1).gt.Bpts(5,n)) then
                            tmplc1(1) = Bpts(5,n)
                         else if(tmplc2(1).lt.Bpts(5,n)) then
                            tmplc2(1) = Bpts(5,n)
                         endif
                      enddo
                   else
                      do n = 2, Nplocal
                         if(tmplc1(2).gt.Bpts(5,n)) then
                            tmplc1(2) = Bpts(5,n)
                         else if(tmplc2(2).lt.Bpts(5,n)) then
                            tmplc2(2) = Bpts(5,n)
                         endif
                      enddo
                   endif
                   call MPI_ALLREDUCE(tmplc1,tmp1,2,MPI_DOUBLE_PRECISION,&
                        MPI_MIN,mpicommwd,ierr)
                   call MPI_ALLREDUCE(tmplc2,tmp2,2,MPI_DOUBLE_PRECISION,&
                        MPI_MAX,mpicommwd,ierr)
                   if(myidy.lt.npyhalf) then
                      zmin = tmp1(1) - eps*(tmp2(1)-tmp1(1))
                      zmax = tmp2(1) + eps*(tmp2(1)-tmp1(1))
                   else
                      zmin = tmp1(2) - eps*(tmp2(2)-tmp1(2))
                      zmax = tmp2(2) + eps*(tmp2(2)-tmp1(2))
                   endif
                else
                   zmin = 0.0
                   zmax = 0.0
                endif

                !// Set flag for 3d luminosity output
                if(mod(i-1,nfrqlum).eq.0) then
                   flaglum = 1
                else
                   flaglum = 0
                endif

                if(myidy.lt.npyhalf) then
                   coef = Bchargelist1(idbunch)*Bchargelist2(idbunchO)*qe*&
                        Bcurrlist2(idbunchO)*(1.0+betas1*betas2)/ &
                        (ga1*betas1*(betas1+betas2)*2*pi*epsilon0*Bmasslist1(idbunch))
                else
                   coef = Bchargelist2(idbunch)*Bchargelist1(idbunchO)*qe*&
                        Bcurrlist1(idbunchO)*(1.0+betas1*betas2)/ &
                        (ga2*betas2*(betas1+betas2)*2*pi*epsilon0*Bmasslist2(idbunch))
                endif

                !//do strong-strong beam-beam kick
                if((Nslice1.eq.1).and.(Nslice2.eq.1)) then
                   !//single slice model
                   if(flagrange.eq.1) then !//fixed comp. domain, reuse Green function.
                      call bb1Slcfix_Beambeam(nprow,npyhalf,innx,inny,&
                           Nplocal,myidy,Bpts,grid2d,Np,coef,Ny1,commcol,nxpylc2,grn,&
                           xpystable,pytable,hx,hy,xmin,ymin)
                   else !//dynamic computation range.
                      if(flagwk.eq.10) then !soft-Gaussian model
                         call bb1SlcGaus_Beambeam(npyhalf,Nplocal,myidy,Bpts,Np,coef)
                      else
                         call bb1Slc_Beambeam(nprow,npyhalf,innx,inny,&
                              Nplocal,myidy,Bpts,grid2d,Ageom,Np,coef,Ny1,commcol,shift)
                      endif
                   endif
                else
                   !//multiple slice model
                   !//To include the effects of Ez, just replace bbmSlcfix by
                   !//bbmSlcfixEz or bbmSlc by bbmSlcEz. However, the effects
                   !//of Ez is very small. The computational time with Ez will
                   !//two third more than the case without Ez.
                   if(flagrange.eq.1) then !//fixed comp. domain, reuse Green function.
                      call bbmSlcfixnew_Beambeam(nprow,npyhalf,innx,inny,Nslice1,Nslice2,&
                           Nslice,NsliceO,jendmax,Nplocal,myidx,myidy,Bpts,grid2d,Ageom,zmin,&
                           zmax,Np,coef,Ny1,commrow,commcol,nz,maxpts1,nxpylc2,&
                           grn,xpystable,pytable,hx,hy,xmin,ymin,Bcurr,flaglum,lum3d)
                   else
                      if(flagwk.eq.10) then !soft-Gaussian model
                         call bbmSlcnewGauss_Beambeam(nprow,npyhalf,Nslice1,Nslice2,Nslice,&
                              NsliceO,jendmax,Nplocal,myidx,myidy,Bpts,zmin,zmax,&
                              Np,coef,commrow,commcol,nz,maxpts1,Bcurr,flaglum,lum3d,NpO)
                      else
                         call bbmSlcnew_Beambeam(nprow,npyhalf,innx,inny,Nslice1,Nslice2,&
                              Nslice,NsliceO,jendmax,Nplocal,myidx,myidy,Bpts,grid2d,&
                              Ageom,zmin,zmax,Np,coef,Ny1,commrow,commcol,nz,maxpts1,shift,&
                              Bcurr,flaglum,lum3d)
                      endif
                   endif
                endif

                !//transfer back to the lab frame.
                call transfer2(Bpts,alpha,phi,Nplocal)

                !// Apply the restore crab kick (for local scheme)
                call crabTransform(Bpts,Nplocal,myidy,npyhalf,1,betaCrab,crabFreq,crabVnorm,&
                     phi_s,bx,by,close2g,alpha,cfreq2,cvnorm2)

                Bptstmp = Bpts
                BptsArry(:,:,idbunch) = Bpts(:,:)
                !// Transfer to next IP
                !// Apply the "ic_th" linear map.
                if(flagmap2nd.eq.1) then
                   call kick_Extmap2nd(Bpts,Nplocal,Amap2nd(ic),close2g)
                else !//this works only for one collision pt per turn
                   !no kick for a test
                   szspz = sigtmp1(5)/sigtmp1(6)
                   call kick_Linearmap(Bpts,Nplocal,Alinearmap,close2g,szspz)
                endif

                !// Apply the "ic_th" nonlinear map.
                if(flagmap3rd.eq.1) then
                   call kick_Extmap3rd(Bptstmp,Nplocal,Amap3rd(ic),close2g)
                   Bpts = Bpts + Bptstmp
                   Bptstmp(:,:) = BptsArry(:,:,idbunch)
                endif
                if(flagmap4th.eq.1) then
                   call kick_Extmap4th(Bptstmp,Nplocal,Amap4th(ic),close2g)
                   Bpts = Bpts + Bptstmp
                endif

                Bcurr = dble(Np)*Bcurr0/Np0
                Nplist(idbunch) = Np
                Nplclist(idbunch) = Nplocal

                !// Save 3D luminosity
                if(flaglum==1) then
                   if(myid.eq.0) then 
                      write(lumFile,'(1x,i9.1,1x,i9.1,1x,i9.1,1x,i9.1,es16.6)')i,ic,idbunch,idbuncho,lum3d
                      call flush(lumFile)
                   endif
                   !// Save 2d luminosity for single slice model
                   if(Nslice1.eq.1 .and. Nslice2.eq.1) then
                      call luminosity2G_Output(Bpts,Nplocal,nxlum,nylum,i,Bcurr,Np,myidy,npyhalf)
                   endif
                endif

                !// Apply the crab kick for next step
                call crabTransform(Bpts,Nplocal,myidy,npyhalf,0,betaCrab,crabFreq,crabVnorm,&
                     phi_s,bx,by,close2g,alpha,cfreq2,cvnorm2)

                BptsArry(:,:,idbunch) = Bpts(:,:)
500             continue
             enddo !//end loop through collisions points per step
          enddo !//end loop through collision steps per turn


          !// Loop through the bunches for once-per-turn actions
          do idbunch = 1, nbunch
             !// Copy coordinates of particles in current bunch into temporary array
             Bpts(:,:) = BptsArry(:,:,idbunch)

             !// Apply the restore crab kick to perform following actions with untilted bunches
             !// This is a "virtual" kick, i. e. not taking place in the machine and therefore without noise
             call crabTransform(Bpts,Nplocal,myidy,npyhalf,1,betaCrab,crabFreq,crabVnorm,&
                  phi_s,bx,by,close2g,alpha,cfreq2,cvnorm2)

             !//apply chromaticity kick
             !//ax - alpha x, bx - beta x, qx - linear chromaticity in x.
             if(abs(qx).gt.1.0d-6 .or. abs(qy).gt.1.0d-6) then
                call chromkick(Bpts,Nplocal,close2g,ax,bx,ay,by,qx,qy)
             endif

             !//apply radiation damping and quantumn excitation.
             !//sigtmp1 is the equlibrium sigma or initial sigma.
             if(flagdamp.ne.0) then
                call kickold_Dampfluc(Bpts,Nplocal,Adampfluc,sigtmp1,z,myid,&
                     close2g,iseedinit,icount)
             endif

             !//do close-orbit sweeping
             if(isweep.eq.1) then
                if(i.eq.1) then
                   closeold = close2g
                endif
                call sweep(Bpts,close2g,swrad,swtune,closeold,Nplocal,i)
             endif

             !//do close-orbit squeeze
             if(iclosq.eq.1) then
                delta = (closess - closeinit)/nmeet
                call clo_orb_squeeze(Bpts,close2g,delta,Nplocal,nmeet,i)
             endif

             !//do feedback control
             if(ifdbk.eq.1) then
                call shiftsigma_Utility(Bpts,Nplocal,Np,myidy,npyhalf,shift,center,sigma)
                call feedback(Bpts,close2g,center,Nplocal)
             endif

             z = i
             !//calculate moments of distribution.
             if(momOutRate /= 0 .and. mod(i,momOutRate).eq.0) then
                call diagnostic2GMBeic_Output(z,Bpts,Nplocal,Np,myidx,myidy,nprow,npcol,&
             commrow,commcol,comm2d,idbunch,ax,bx,ay,by)
             endif

             !// Save particle distribution in phase space
             if(ptOut>0 .and. mod(i,ptRate).eq.0) then
                call ptclout_Output(Bpts,myid,i,idbunch,Nplocal,ptOut,nturn,ptRate,ptSampleInd)
             endif
             if(np_xy>0) call xySeries_Output(Bpts,myid,i,idbunch,Nplocal,np_xy,&
                  repRate_xy,len_xy,Nturn,ptr)

             !//calculate 2D luminosity
             !            call luminosity2G_Output(Bpts,Nplocal,nxlum,nylum,i,&
             !                                Bcurr1,Bcurr2,Np1,Np2,myidy,npyhalf)

             !// Apply the crab kick for next turn
             !// This is a "virtual" kick, i. e. not taking place in the machine and therefore without noise
             call crabTransform(Bpts,Nplocal,myidy,npyhalf,0,betaCrab,crabFreq,crabVnorm,&
                  phi_s,bx,by,close2g,alpha,cfreq2,cvnorm2)

             !// Update list with coordinates
             BptsArry(:,:,idbunch) = Bpts(:,:)
          enddo

          call MPI_BARRIER(mpicommwd,ierr)
          time1 = MPI_WTIME()
          timepasslc = time1 - time0
          call MPI_ALLREDUCE(timepasslc,timepass,1,MPI_DOUBLE_PRECISION,MPI_SUM,mpicommwd,ierr)
          timepass = timepass/npcol/nprow

          !// Eventually save particle distribution to resume simulation
          if((timepass.gt.tmax).and.(saveAll.eq.1)) then
             do idbunch = 1, nbunch
                Bpts(:,:) = BptsArry(:,:,idbunch)
                call phaseoutBB_Output(Bpts,myid,i,idbunch,Nplocal,close2g,itic,icount)
             enddo
             saveAll = 0
          endif
       enddo !loop through Nturn

       !// Save final particle distribution to resume simulation (if not saved before)
       if(saveAll.eq.1) then
          do idbunch = 1, nbunch
             Bpts(:,:) = BptsArry(:,:,idbunch)
             call phaseoutBB_Output(Bpts,myid,i,idbunch,Nplocal,close2g,itic,icount)
          enddo
       endif

       deallocate(Nplist)
       deallocate(Nplclist)
       deallocate(Bpts)
       deallocate(Bptstmp)
       deallocate(BptsArry)
       deallocate(Bcurrlist1)
       deallocate(Bcurrlist2)
       deallocate(Bchargelist1)
       deallocate(Bchargelist2)
       deallocate(Bmasslist1)
       deallocate(Bmasslist2)
       nullify(ptr)
       deallocate(ptSampleInd)
       deallocate(QsampleInd)
       if(flagmap2nd.eq.1) then
          deallocate(Amap2nd)
          deallocate(Amap2ndArry)
       endif
       if(flagmap3rd.eq.1) then
          deallocate(Amap3rd)
          deallocate(Amap3rdArry)
       endif
       if(flagmap4th.eq.1) then
          deallocate(Amap4th)
          deallocate(Amap4thArry)
       endif
       if(flagrange.eq.1) then
          deallocate(grn)
          deallocate(table)
          deallocate(pytable)
          deallocate(xpystable)
       endif
       close(1)

    else if(flagGroup.eq.1) then !weak-strong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !//Strong-weak beam-beam interaction below.
       !//To do a strong-strong simulation using only 2 bunches with multiple
       !//collision points.
       if(myid.eq.0) print*,"strong-weak beam-beam...."
       Np01 = Np1 !//initial total # of ptcls.
       Bcurr01 = Bcurr1 !//initial current
       Np02 = Np2 !//initial total # of ptcls.
       Bcurr02 = Bcurr2 !//initial current
       close1init = close1
       close2init = close2
       !//Get local grid number in each computaional domain
       innx = Nx1
       inny = Ny1/npcol

       sigtmp1 = sigma1
       !//Output initial moments
       call findmoments12_Utility(Bpts1,Nplocal1,Np1,centrl1,sigma1)
       !//rescale the particles coordinates to the prescribed value
       !          do i = 1, Nplocal1
       !            Bpts1(1,i) = Bpts1(1,i)*sigtmp1(1)/sigma1(1)
       !            Bpts1(2,i) = Bpts1(2,i)*sigtmp1(2)/sigma1(2)
       !            Bpts1(3,i) = Bpts1(3,i)*sigtmp1(3)/sigma1(3)
       !            Bpts1(4,i) = Bpts1(4,i)*sigtmp1(4)/sigma1(4)
       !            Bpts1(5,i) = Bpts1(5,i)*sigtmp1(5)/sigma1(5)
       !            Bpts1(6,i) = Bpts1(6,i)*sigtmp1(6)/sigma1(6)
       !          enddo
       call diagnostic1G1_Output(z,Bpts1,Nplocal1,Np1,centrl1)
       sigtmp2 = sigma2
       !//Output initial moments
       call findmoments12_Utility(Bpts2,Nplocal2,Np2,centrl2,sigma2)
       !//rescale the particles coordinates to the prescribed value
       !          do i = 1, Nplocal2
       !            Bpts2(1,i) = Bpts2(1,i)*sigtmp2(1)/sigma2(1)
       !            Bpts2(2,i) = Bpts2(2,i)*sigtmp2(2)/sigma2(2)
       !            Bpts2(3,i) = Bpts2(3,i)*sigtmp2(3)/sigma2(3)
       !            Bpts2(4,i) = Bpts2(4,i)*sigtmp2(4)/sigma2(4)
       !            Bpts2(5,i) = Bpts2(5,i)*sigtmp2(5)/sigma2(5)
       !            Bpts2(6,i) = Bpts2(6,i)*sigtmp2(6)/sigma2(6)
       !          enddo
       call diagnostic2G1_Output(z,Bpts2,Nplocal2,Np2,centrl2)

       shift(1) = centrl2(1)-centrl1(1)
       shift(2) = centrl2(2)-centrl1(2)

       !initialization required for the Gaussian approximation, flag 1
       call tablew

       call MPI_BARRIER(mpicommwd,ierr)
       !//Output initial 2D luminosity
       i = 0
       flaglum = 0  !//flag for 3d luminosity calculation.
       open(1,file="luminosity.data",status="unknown")

       !//find the Green function for the fixed computation domain.
       !//The shift between two beam has to be 0 or fixed.
       if(flagrange.eq.1) then
          shift = 0.0d0

          nxtot = Nx1
          nytot = Ny1
          !// +1 is from the real to complex fft.
          nsxy1 = (nxtot+1)/npcol
          nsxy2 = (nxtot+1) - npcol*nsxy1
          allocate(xpystable(0:npcol-1))
          do i = 0, npcol-1
             if(i.le.(nsxy2-1)) then
                xpystable(i) = nsxy1+1
             else
                xpystable(i) = nsxy1
             endif
          enddo
          nxpylc2 = xpystable(myidy)

          allocate(grn(2*nytot,nxpylc2))

          range1(1) = Extrange1(1)
          range1(2) = Extrange1(2)
          range1(3) = Extrange1(3)
          range1(4) = Extrange1(4)
          range2(1) = Extrange2(1)
          range2(2) = Extrange2(2)
          range2(3) = Extrange2(3)
          range2(4) = Extrange2(4)
          call update2d1G_CompDom(Ageom,range1,range2,grid2d)
          call getmsize_CompDom(Ageom,msize)
          hx = msize(1)
          hy = msize(2)
          allocate(table(2,0:nprow-1,0:npcol-1))
          call getlctabnm_CompDom(Ageom,table)
          allocate(pytable(0:npcol-1))
          do i = 0, npcol-1
             pytable(i) = table(2,0,i)
          enddo

          !//get the Green function for the fixed domain.
          !//here, you can not have the domain shift
          !             call greenf2d1G(nxtot,nytot,nxpylc2,myidy,npcol,xpystable,hx,hy,&
          !                            grn,shift)
          !using parallel green function
          call greenf2d1G(nxtot,nytot,nxpylc2,hx,hy,myidy,npcol,commcol,&
               xpystable,grn,shift)


          call getrange_CompDom(Ageom,tmprange)
          xmin = tmprange(1)
          ymin = tmprange(3)
       endif !//finish the Green function for the fixed domain.

       itic = 0
       !//Loop through "Nturn" simulation.
       do i = 1, Nturn
          !//Loop through "Ncl" collision points
          !//The maximum Ncl is 2 for single bunch S-S beam-beam
          do ic = 1, Ncl
             itic = itic + 1
             !//We need to find the shift between the slice at each
             !//collision point: center2 - center1
             Nslice1 = Nslice1list(ic)
             Nslice2 = Nslice2list(ic)
             jendmax = ((min(Nslice1,Nslice2)-1)/nprow + 1)*nprow
             maxpts1 = int(Np1/Nslice1*factor) !//maximum # of ptl. per slice
             maxpts2 = int(Np2/Nslice2*factor) !//maximum # of ptl. per slice

             !//both shift and sigma are needed
             call findmoments12_Utility(Bpts1,Nplocal1,Np1,centrl1,sigma1)
             call findmoments12_Utility(Bpts2,Nplocal2,Np2,centrl2,sigma2)
             shift(1) = centrl2(1)-centrl1(1)
             shift(2) = centrl2(2)-centrl1(2)

             !//from Lab coordinate to boost coordinate
             !//when deal with crossing angle collision.
             !//the particle transverse size is checked to
             !//decide the lost particle.
             !//phi is the half crossing angle in x-z, alpha is the cross angle in x-y
             call transfer1old2(Bpts1,alpha,phi,Nplocal1)
             call transfer1old2(Bpts2,alpha,phi,Nplocal2)

             !//find zmin and zmax of each beam
             !//The following can be simplified using 2 communicators
             !//The following can be replaced by fixed maximum range
             tmplc1 = 1.0e10
             tmplc2 = -1.0e10
             if(Nslice1.gt.1 .or. Nslice2.gt.1) then
                if(flagrange.eq.1) then !// fixed external range
                   zmin1 = Extrange1(5)
                   zmax1 = Extrange1(6) 
                   zmin2 = Extrange2(5)
                   zmax2 = Extrange2(6) 
                else !// this could be improved by using MPI_ALLREC
                   do n = 1, Nplocal1
                      if(tmplc1(1).gt.Bpts1(5,n)) then
                         tmplc1(1) = Bpts1(5,n)
                      else if(tmplc2(1).le.Bpts1(5,n)) then
                         tmplc2(1) = Bpts1(5,n)
                      else
                      endif
                   enddo
                   do n = 1, Nplocal2
                      if(tmplc1(2).gt.Bpts2(5,n)) then
                         tmplc1(2) = Bpts2(5,n)
                      else if(tmplc2(2).le.Bpts2(5,n)) then
                         tmplc2(2) = Bpts2(5,n)
                      else
                      endif
                   enddo
                   call MPI_ALLREDUCE(tmplc1,tmp1,2,MPI_DOUBLE_PRECISION,&
                        MPI_MIN,mpicommwd,ierr)
                   call MPI_ALLREDUCE(tmplc2,tmp2,2,MPI_DOUBLE_PRECISION,&
                        MPI_MAX,mpicommwd,ierr)
                   zmin1 = tmp1(1) - eps*(tmp2(1)-tmp1(1))
                   zmax1 = tmp2(1) + eps*(tmp2(1)-tmp1(1))
                   zmin2 = tmp1(2) - eps*(tmp2(2)-tmp1(2))
                   zmax2 = tmp2(2) + eps*(tmp2(2)-tmp1(2))
                endif
             else
                zmin1 = 0.0
                zmax1 = 0.0
                zmin2 = 0.0
                zmax2 = 0.0
             endif

             !//choose 3d luminosity
             if(mod(itic,nfrqlum).eq.0) then
                flaglum = 1
             else
                flaglum = 0
             endif
             !//do weak-strong beam-beam kick
             if(flagwk.eq.1) then !//arbitray distribution of strong beam
                !//multiple slice model
                if(flagrange.eq.1) then !//fixed comp. domain, reuse Green function.
                   call bbmslswfix_Beambeam(nprow,npcol,innx,inny,Nslice1,Nslice2,&
                        jendmax,Nplocal1,Nplocal2,myidx,myidy,Bpts1,Bpts2,grid2d,Ageom,&
                        zmin1,zmax1,zmin2,zmax2,Np1,Np2,coef,Ny1,commrow,commcol,&
                        nz1,nz2,maxpts1,maxpts2,nxpylc2,grn,xpystable,pytable,hx,hy,&
                        xmin,ymin,Bcurr1,Bcurr2,flaglum,lum3d)
                else
                   call bbmslsw_Beambeam(nprow,npcol,innx,inny,Nslice1,Nslice2,&
                        jendmax,Nplocal1,Nplocal2,myidx,myidy,Bpts1,Bpts2,grid2d,Ageom,&
                        zmin1,zmax1,zmin2,zmax2,Np1,Np2,coef,Ny1,commrow,commcol,&
                        nz1,nz2,maxpts1,maxpts2,shift,Bcurr1,Bcurr2,flaglum,lum3d)
                endif
             else if(flagwk.eq.2) then !//Gaussian approximation of strong beam
                call bbmslswgauss_Beambeam(nprow,npcol,innx,inny,Nslice1,Nslice2,&
                     jendmax,Nplocal1,Nplocal2,myidx,myidy,Bpts1,Bpts2,grid2d,Ageom,&
                     zmin1,zmax1,zmin2,zmax2,Np1,Np2,coef,Ny1,commrow,commcol,&
                     nz1,nz2,maxpts1,maxpts2,Bcurr1,Bcurr2,flaglum,lum3d)
             endif

             !//transfer back to the lab frame.
             call transfer2(Bpts1,alpha,phi,Nplocal1)
             call transfer2(Bpts2,alpha,phi,Nplocal2)
             
             !//apply the "ic_th" linear map.
             if(flagmap2nd.eq.1) then
                call kick_Extmap2nd(Bpts1,Nplocal1,Amap2nd1(ic),close1)
                call kick_Extmap2nd(Bpts2,Nplocal2,Amap2nd2(ic),close2)
             else !//this works only for one collision pt per turn
                szspz = sigtmp1(5)/sigtmp1(6)
                call kick_Linearmap(Bpts1,Nplocal1,Alinearmap1,close1,szspz)
                szspz = sigtmp2(5)/sigtmp2(6)
                call kick_Linearmap(Bpts2,Nplocal2,Alinearmap2,close2,szspz)
             endif
             !//apply the "ic_th" nonlinear map.
             if(flagmap3rd.eq.1) then
                call kick_Extmap3rd(Bpts1,Nplocal1,Amap3rd1(ic),close1)
                call kick_Extmap3rd(Bpts2,Nplocal2,Amap3rd2(ic),close2)
             endif
             if(flagmap4th.eq.1) then
                call kick_Extmap4th(Bpts1,Nplocal1,Amap4th1(ic),close1)
                call kick_Extmap4th(Bpts2,Nplocal2,Amap4th2(ic),close2)
             endif
             !//update the close orbit to current location.
             !//whenever you update the particles using external map, you
             !//should update the close orbit.
             !call update_clorbit(close2g)
             !//Check particle lost ?????????
             !              call ptlost1_Utility(Bpts1,xmin,xmax,ymin,ymax,Nplocal1,Np1)
             !              Bcurr1 = dble(Np1)*Bcurr01/Np01
             !              call ptlost2_Utility(Bpts2,xmin,xmax,ymin,ymax,Nplocal2,Np2)
             !              Bcurr2 = dble(Np2)*Bcurr02/Np02

             !//choose the 3d luminosity printout.
             if(mod(itic,nfrqlum).eq.0) then
                !test purpose
                !                flaglum = 1
                if(myid.eq.0) then 
                   write(1,1000)itic*1.0,lum3d
                   call flush(1)
                endif
             endif

          enddo !//end loop through collisions per turn

          !//apply chromaticity kick
          !//ax - alpha x, bx - beta x, qx - linear chromaticity in x.
          !            call chromkick(Bpts1,Nplocal1,close1,ax1,bx1,ay1,by1,qx1,qy1)
          !            call chromkick(Bpts2,Nplocal2,close2,ax2,bx2,ay2,by2,qx2,qy2)

          !//apply radiation damping and quantumn excitation.
          !//sigtmp1 is the equlibrium sigma or initial sigma.
          if(flagdamp1.ne.0) then
             call kickold_Dampfluc(Bpts1,Nplocal1,Adampfluc1,sigtmp1,z,myid,&
                  close1,iseedinit,icount)
             !              call kick_Dampfluc(Bpts1,Nplocal1,Adampfluc1,sigtmp1,z,myid,close1,&
             !                               ax1,bx1,ay1,by1)
          endif
          if(flagdamp2.ne.0) then
             call kickold_Dampfluc(Bpts2,Nplocal2,Adampfluc2,sigtmp2,z,myid,&
                  close2,iseedinit,icount)
             !              call kick_Dampfluc(Bpts2,Nplocal2,Adampfluc2,sigtmp2,z,myid,close2,&
             !                               ax2,bx2,ay2,by2)
          endif

          !//do close-orbit sweeping
          if(isweep1.eq.1) then
             if(i.eq.1) then
                closeold1(1) = swrad1
                closeold1(2) = 0.0
             endif
             call sweep(Bpts1,close1,swrad1,swtune1,closeold1,Nplocal1,i)
          endif
          if(isweep2.eq.1) then
             if(i.eq.1) then
                closeold2(1) = swrad2
                closeold2(2) = 0.0
             endif
             call sweep(Bpts2,close2,swrad2,swtune2,closeold2,Nplocal2,i)
          endif
          !//do close-orbit squeeze
          if(iclosq1.eq.1) then
             delta1 = (close1ss - close1init)/nmeet
             call clo_orb_squeeze(Bpts1,close1,delta1,Nplocal1,nmeet,i)
          endif
          if(iclosq2.eq.1) then
             delta2 = (close2ss - close2init)/nmeet
             call clo_orb_squeeze(Bpts2,close2,delta2,Nplocal2,nmeet,i)
          endif

          !//do feedback control
          if(ifdbk.eq.1) then
             call findmoments12_Utility(Bpts1,Nplocal1,Np1,centrl1,sigma1)
             call feedback(Bpts1,close1,centrl1,Nplocal1)
             call findmoments12_Utility(Bpts2,Nplocal2,Np2,centrl2,sigma2)
             call feedback(Bpts2,close2,centrl2,Nplocal2)
          endif

          z = z + 1
          !//calculate moments of distribution.
          call diagnostic1G1_Output(z,Bpts1,Nplocal1,Np1,centrl1)
          call diagnostic2G1_Output(z,Bpts2,Nplocal2,Np2,centrl2)

       enddo !loop through Nturn

       deallocate(Bpts1)
       deallocate(Bpts2)
       if(flagmap2nd.eq.1) then
          deallocate(Amap2nd1)
          deallocate(Amap2nd2)
       endif
       if(flagmap3rd.eq.1) then
          deallocate(Amap3rd1)
          deallocate(Amap3rd2)
       endif
       if(flagmap4th.eq.1) then
          deallocate(Amap4th1)
          deallocate(Amap4th2)
       endif
       if(flagrange.eq.1) then
          deallocate(grn)
          deallocate(table)
          deallocate(pytable)
          deallocate(xpystable)
       endif
       close(1)
    endif

    ! final output.
1000 format(2(1x,es14.7))
1001 format(2(1x,i13.1),5(1x,es13.6))
    call MPI_BARRIER(mpicommwd,ierr)
    t_integ = t_integ + elapsedtime_Timer(t0)
    call showtime_Timer()

  end subroutine run_AccSimulator

  subroutine destruct_AccSimulator(time)
    implicit none
    include 'mpif.h'
    double precision :: time

    call destruct_CompDom(Ageom)
    call end_Output(time)

  end subroutine destruct_AccSimulator

  !given the linear transfer matrix "amt", output symplectic matrix "bmt"
  !testsum: a check of symplectic condition. It is 0 for exact symplectic matrix.
  subroutine symplect(amt,np,testsum)
    integer, intent(in) :: np
    double precision, intent(out) :: testsum
    double precision, dimension(np,np), intent(inout) :: amt
    double precision, dimension(np,np) :: tmp1,tmp2,tmp3,ss,bmt,tmppp
    integer :: i,j,k
    double precision :: testinv

    !form factor I + F
    do j = 1, np
       do i = 1, np
          tmp1(i,j) = amt(i,j)
       enddo
       tmp1(j,j) = 1.0 + amt(j,j)
    enddo

    !find (I + F)^-1
    call invmt(tmp1,np,tmp2)

    do j = 1, np
       do i = 1, np
          tmp1(i,j) = amt(i,j)
       enddo
       tmp1(j,j) = 1.0 + amt(j,j)
    enddo
    do j = 1, np
       do i = 1, np
          tmppp(i,j) = 0.0
          do k = 1, np
             tmppp(i,j) = tmppp(i,j) + tmp1(i,k)*tmp2(k,j)
          enddo
       enddo
    enddo
    testinv = 0.0
    do j = 1, np
       do i = 1,np
          testinv = testinv + abs(tmppp(i,j))
       enddo
    enddo

    !form factor (I - F)
    do j = 1, np
       do i = 1, np
          tmp3(i,j) = -amt(i,j)
       enddo
       tmp3(j,j) = 1.0 - amt(j,j)
    enddo

    !find (I-F)(I+F)^(-1)
    do j = 1, np
       do i = 1, np
          tmp1(i,j) = 0.0
          do k = 1, np
             tmp1(i,j) = tmp1(i,j) + tmp3(i,k)*tmp2(k,j)
          enddo
       enddo
    enddo

    if(mod(np,2).ne.0) then
       print*,"np has to be even: ",np
       stop
    endif

    !form S matrix
    do j = 1, np
       do i = 1, np
          ss(i,j) = 0.0
       enddo
    enddo
    do i = 1, np/2
       ss(2*i-1,2*i) = 1.0
       ss(2*i,2*i-1) = -1.0
    enddo

    !find V = S(I-F)(I+F)^(-1)
    do j = 1, np
       do i = 1, np
          tmp2(i,j) = 0.0
          do k = 1, np
             tmp2(i,j) = tmp2(i,j) + ss(i,k)*tmp1(k,j)
          enddo
       enddo
    enddo

    do i = 1, np
       do j = 1, np
          tmp3(i,j) = tmp2(j,i)
       enddo
    enddo
    !form W = (V + V^(t))/2
    do j = 1, np
       do i = 1, np
          tmp1(i,j) = (tmp2(i,j)+tmp3(i,j))/2
       enddo
    enddo

    !form S W
    do j = 1, np
       do i = 1, np
          tmp2(i,j) = 0.0
          do k = 1, np
             tmp2(i,j) = tmp2(i,j) + ss(i,k)*tmp1(k,j)
          enddo
       enddo
    enddo

    tmp1 = tmp2

    !form (I-SW)
    do j = 1, np
       do i = 1, np
          tmp2(i,j) = -tmp1(i,j)
       enddo
       tmp2(j,j) = 1.0 - tmp1(j,j)
    enddo

    !form (I+SW) 
    do j = 1, np
       do i = 1, np
          tmp3(i,j) = tmp1(i,j)
       enddo
       tmp3(j,j) = 1.0 + tmp1(j,j)
    enddo

    !find (I-W)^(-1)
    call invmt(tmp2,np,tmp1)

    !find symplectic matrix, (I+SW)(I-SW)^(-1)
    do j = 1, np
       do i = 1, np
          bmt(i,j) = 0.0
          do k = 1, np
             bmt(i,j) = bmt(i,j) + tmp3(i,k)*tmp1(k,j)
          enddo
       enddo
    enddo

    !test the symplectic condition: R^t S R = S
    do j = 1, np
       do i = 1, np
          tmp1(i,j) = 0.0
          do k = 1, np
             tmp1(i,j) = tmp1(i,j) + ss(i,k)*bmt(k,j)
          enddo
       enddo
    enddo

    do i = 1, np
       do j = 1, np
          tmp3(i,j) = bmt(j,i)
       enddo
    enddo
    do j = 1, np
       do i = 1, np
          tmp2(i,j) = 0.0
          do k = 1, np
             !tmp2(i,j) = tmp2(i,j) + bmt(k,i)*tmp1(k,j)
             tmp2(i,j) = tmp2(i,j) + tmp3(i,k)*tmp1(k,j)
          enddo
       enddo
    enddo

    testsum = 0.0
    do j = 1, np
       do i = 1, np
          testsum = testsum + abs(ss(i,j)-tmp2(i,j))
       enddo
    enddo
    amt = bmt

  end subroutine symplect

  !kick from linear machine chromaticity.
  !ax - alpha x, bx - beta x, qx1 - x chromaticity.
  !close: close orbit
  subroutine chromkick(Pts1in,nptlc,close2g,ax,bx,ay,by,qx1,qy1)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, dimension(2), intent(in) :: close2g
    double precision :: ax,bx,ay,by,qx1,qy1
    integer :: i
    double precision :: tmp1,tmp2,pi,gx,gy,dphi
    double precision, dimension(2,2) :: matx1,maty1

    pi = 2*asin(1.0d0)
    gx = (1.0d0+ax*ax)/bx
    gy = (1.0d0+ay*ay)/by
    do i = 1, nptlc
       tmp1 = Pts1in(1,i) - close2g(1)
       tmp2 = Pts1in(2,i)
       dphi = 2*pi*qx1*Pts1in(6,i)
       matx1(1,1) = cos(dphi)+ax*sin(dphi)
       matx1(2,1) = -gx*sin(dphi)
       matx1(1,2) = bx*sin(dphi)
       matx1(2,2) = cos(dphi)-ax*sin(dphi)
       Pts1in(1,i) = close2g(1) + matx1(1,1)*tmp1+matx1(1,2)*tmp2
       Pts1in(2,i) = matx1(2,1)*tmp1+matx1(2,2)*tmp2
       tmp1 = Pts1in(3,i) - close2g(2)
       tmp2 = Pts1in(4,i)
       dphi = 2*pi*qy1*Pts1in(6,i)
       maty1(1,1) = cos(dphi)+ay*sin(dphi)
       maty1(2,1) = -gy*sin(dphi)
       maty1(1,2) = by*sin(dphi)
       maty1(2,2) = cos(dphi)-ay*sin(dphi)
       Pts1in(3,i) = close2g(2) + maty1(1,1)*tmp1+maty1(1,2)*tmp2
       Pts1in(4,i) = maty1(2,1)*tmp1+maty1(2,2)*tmp2
    enddo

  end subroutine chromkick

  !horizontal transfer map for 90 degree phase advance
  subroutine map90_x(Pts1in,nptlc,beta1,beta2,close2g)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, intent(in) :: beta1,beta2
    double precision, dimension(2),intent(in) :: close2g
    double precision ::tmp1,tmp2
    integer :: i

    do i = 1, nptlc
       tmp1 = pts1in(1,i) - close2g(1)
       tmp2 = pts1in(2,i)
       pts1in(1,i) = sqrt(beta1*beta2)*tmp2 + close2g(1)
       pts1in(2,i) = -tmp1/sqrt(beta1*beta2)
    enddo

  end subroutine map90_x

  !horizontal inverse transfer map for 90 degree phase advance
  subroutine invmap90_x(Pts1in,nptlc,beta1,beta2,close2g) 
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, intent(in) :: beta1,beta2
    double precision, dimension(2),intent(in) :: close2g
    double precision :: tmp1,tmp2
    integer :: i

    do i = 1, nptlc
       tmp1 = pts1in(1,i) - close2g(1)
       tmp2 = pts1in(2,i)
       pts1in(1,i) = close2g(1) - sqrt(beta1*beta2)*tmp2
       pts1in(2,i) = tmp1/sqrt(beta1*beta2) 
    enddo

  end subroutine invmap90_x

  !// vertical transfer map for 90 degree phase advance
  subroutine map90_y(Pts1in,nptlc,beta1,beta2,close2g)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, intent(in) :: beta1,beta2
    double precision, dimension(2),intent(in) :: close2g
    double precision ::tmp1,tmp2
    integer :: i

    do i = 1, nptlc
       tmp1 = pts1in(3,i) - close2g(2)
       tmp2 = pts1in(4,i)
       pts1in(3,i) = sqrt(beta1*beta2)*tmp2 + close2g(2)
       pts1in(4,i) = -tmp1/sqrt(beta1*beta2)
    enddo

  end subroutine map90_y

  !// vertical inverse transfer map for 90 degree phase advance
  subroutine invmap90_y(Pts1in,nptlc,beta1,beta2,close2g) 
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision :: beta1,beta2,tmp1,tmp2
    double precision, dimension(2) :: close2g
    integer :: i

    do i = 1, nptlc
       tmp1 = pts1in(3,i) - close2g(2)
       tmp2 = pts1in(4,i)
       pts1in(3,i) = close2g(2) - sqrt(beta1*beta2)*tmp2
       pts1in(4,i) = tmp1/sqrt(beta1*beta2) 
    enddo

  end subroutine invmap90_y

  !thin lens crab cavity model in X-Z plane
  subroutine crabcavitythinXZorg(Pts1in,nptlc,freq,vnorm,phi_s,freq2,vnorm2)
    implicit none
    include 'mpif.h'
    double precision, intent(in) :: vnorm,freq,phi_s,vnorm2,freq2
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:), intent(inout) :: Pts1in
    double precision :: woc,woc2,v22,v11

    woc = 4*asin(1.0d0)*freq/2.99792458d8
    woc2 = 4*asin(1.0d0)*freq2/2.99792458d8

    if (woc2>0) then
       v22 = vnorm*vnorm2*woc/woc2
       v11 = vnorm*(1.0-vnorm2)
    else
       v22 = 0.0
       v11 = vnorm
    endif
   
    Pts1in(2,1:nptlc) = Pts1in(2,1:nptlc) + v11*sin(woc*Pts1in(5,1:nptlc)+phi_s) + &
                                            v22*sin(woc2*Pts1in(5,1:nptlc)+phi_s)
    Pts1in(6,1:nptlc) = Pts1in(6,1:nptlc) + (v11*woc*cos(woc*Pts1in(5,1:nptlc)+phi_s)+&
                        v22*woc2*cos(woc2*Pts1in(5,1:nptlc)+phi_s))*Pts1in(1,1:nptlc)

  end subroutine crabcavitythinXZorg

  !// thin lens crab cavity model in Y-Z plane
  subroutine crabcavitythinYZorg(Pts1in,nptlc,freq,vnorm,phi_s,freq2,vnorm2)
    implicit none
    include 'mpif.h'
    double precision, intent(in) :: vnorm,freq,phi_s,vnorm2,freq2
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:), intent(inout) :: Pts1in
    double precision :: woc,woc2,v22,v11

    woc = 4*asin(1.0d0)*freq/2.99792458d8
    woc2 = 4*asin(1.0d0)*freq2/2.99792458d8

    if (woc2>0) then
       v22 = vnorm*vnorm2*woc/woc2
       v11 = vnorm*(1.0-vnorm2)
    else
       v22 = 0.0
       v11 = vnorm
    endif

    Pts1in(4,1:nptlc) = Pts1in(4,1:nptlc) + v11*sin(woc*Pts1in(5,1:nptlc)+phi_s) + &
                                            v22*sin(woc2*Pts1in(5,1:nptlc)+phi_s)
    Pts1in(6,1:nptlc) = Pts1in(6,1:nptlc) + (v11*woc*cos(woc*Pts1in(5,1:nptlc)+phi_s)+&
                        v22*woc2*cos(woc2*Pts1in(5,1:nptlc)+phi_s))*Pts1in(3,1:nptlc)

  end subroutine crabcavitythinYZorg

  subroutine crabTransform(Bpts,Nplocal,myidy,npyhalf,reverseFlag,betaCrab,crabFreq,crabVnorm,&
    phi_s,bx,by,close2g,alpha,cfreq2,cvnorm2)
    integer, intent(in) :: Nplocal,myidy,npyhalf,reverseFlag
    double precision, intent(in) :: betaCrab,crabFreq,crabVnorm,phi_s,bx,by,alpha,&
                                    cfreq2,cvnorm2
    double precision, dimension(2), intent(in) :: close2g
    double precision, pointer, dimension(:,:), intent(inout) :: Bpts
    double precision :: pi

    pi = 2*asin(1.d0)

    if(reverseFlag==0) then  !// Crab bunch
       if(abs(alpha)<=1.0d-9) then
          call invmap90(Bpts,Nplocal,betaCrab,bx,close2g)
          call crabcavitythinXZorg(Bpts,Nplocal,crabfreq,crabVnorm,phi_s,cfreq2,cvnorm2)
          call map90(Bpts,Nplocal,betaCrab,bx,close2g)
       else
          call invmap90(Bpts,Nplocal,betaCrab,by,close2g)
          call crabcavitythinYZorg(Bpts,Nplocal,crabfreq,crabVnorm,phi_s,cfreq2,cvnorm2)
          call map90(Bpts,Nplocal,betaCrab,by,close2g)
       endif
    endif

    !// Uncrab bunch
    if(reverseFlag==1) then
       if(abs(alpha)<=1.0d-9) then
          call map90(Bpts,Nplocal,betaCrab,bx,close2g)
          call crabcavitythinXZorg(Bpts,Nplocal,crabfreq,crabVnorm,phi_s,cfreq2,cvnorm2)
          call invmap90(Bpts,Nplocal,betaCrab,bx,close2g)
       else
          call map90(Bpts,Nplocal,betaCrab,by,close2g)
          call crabcavitythinYZorg(Bpts,Nplocal,crabfreq,crabVnorm,phi_s,cfreq2,cvnorm2)
          call invmap90(Bpts,Nplocal,betaCrab,by,close2g)
       endif
    endif

  end subroutine crabTransform

  subroutine normdv2(iseed,y)
    implicit none
    include 'mpif.h'
    integer, intent(inout) :: iseed
    double precision, dimension(2), intent(out) :: y
    double precision :: twopi,x1,x2,epsilon

    epsilon = 1.0d-18

    twopi = 4.0d0*asin(1.0d0)
    x2 = ran222(iseed)
10    x1 = ran222(iseed)
!    call random_number(x2)
!10  call random_number(x1)
    if(x1.eq.0.0d0) goto 10
    !        if(x1.eq.0.0) x1 = epsilon
    y(1) = sqrt(-2.0d0*log(x1))*cos(twopi*x2)
    y(2) = sqrt(-2.0d0*log(x1))*sin(twopi*x2)

  end subroutine normdv2

      FUNCTION ran222(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran222,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2d-7,RNMX=1.d0-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
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
      ran222=min(AM*iy,RNMX)
      return
      END Function ran222


  !horizontal transfer map for 90 degree phase advance
  !assume betax = betay
  subroutine map90(Pts1in,nptlc,beta1,beta2,close2g)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, intent(in) :: beta1,beta2
    double precision, dimension(2),intent(in) :: close2g
    double precision ::tmp1,tmp2
    integer :: i
 
    do i = 1, nptlc
       tmp1 = pts1in(1,i) - close2g(1)
       tmp2 = pts1in(2,i)
       pts1in(1,i) = sqrt(beta1*beta2)*tmp2 + close2g(1)
       pts1in(2,i) = -tmp1/sqrt(beta1*beta2)
       tmp1 = pts1in(3,i) - close2g(2)
       tmp2 = pts1in(4,i)
       pts1in(3,i) = sqrt(beta1*beta2)*tmp2 + close2g(2)
       pts1in(4,i) = -tmp1/sqrt(beta1*beta2)
    enddo
 
  end subroutine map90

  !// vertical inverse transfer map for 90 degree phase advance
  subroutine invmap90(Pts1in,nptlc,beta1,beta2,close2g)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision :: beta1,beta2,tmp1,tmp2
    double precision, dimension(2) :: close2g
    integer :: i
 
    do i = 1, nptlc
       tmp1 = pts1in(1,i) - close2g(1)
       tmp2 = pts1in(2,i)
       pts1in(1,i) = close2g(1) - sqrt(beta1*beta2)*tmp2
       pts1in(2,i) = tmp1/sqrt(beta1*beta2)
       tmp1 = pts1in(3,i) - close2g(2)
       tmp2 = pts1in(4,i)
       pts1in(3,i) = close2g(2) - sqrt(beta1*beta2)*tmp2
       pts1in(4,i) = tmp1/sqrt(beta1*beta2)
    enddo
 
  end subroutine invmap90

end module AccSimulatorclass
