module Dampflucclass
  use Utilityclass

  type Dampfluc
     double precision :: lambdax,lambday,lambdaz
  end type Dampfluc
  interface construct_Dampfluc
     module procedure init1_Dampfluc, init2_Dampfluc
  end interface
contains
  ! generate damping-related quantities (reference: "Comparisons of 
  ! beam-beam simulations"). Definitions:
  ! rlx=exp(-1/taux), rly=exp(-1/tauy), rlz=exp(-1/tauz)
  ! deloe=DeltaE/E=synchr. radiation relative energy loss per turn
  ! 1/taux=del*Jx/2, 1/tauy=del*Jy/2, 1/tauz=del*Jz/2
  ! where taux,tauy and tauz are the damping times in units of turns and
  ! Jx=dpjx, Jy=dpjy, Jz=dpjz, and Jx+Jy+Jz=4. 
  ! For isomagnetic lattice, Jx=1-dd, Jy=1, Jz=2+dd
  ! where dd=damping partition no. (normally, dd=0), so that 
  ! 1/taux+1/tauz=3/tauy
  ! NOTE: tauy and dampart are assumed to be input quantities
  subroutine init1_Dampfluc(this,tauy,dampart)
    implicit none
    include "mpif.h"
    type (Dampfluc), intent(out) :: this
    double precision, intent(in) :: tauy,dampart
    double precision :: deleoe,dpjx,dpjy,dpjz,taux,tauz

    dpjx = 1-dampart
    dpjy = 1
    dpjz = 2+dampart
    deleoe = 2/(tauy*dpjy)
    taux = 2/(deleoe*dpjx)
    tauz = 2/(deleoe*dpjz)
    this%lambdax=exp(-1/taux)
    this%lambday=exp(-1/tauy)
    this%lambdaz=exp(-1/tauz)

  end subroutine init1_Dampfluc



  ! generate damping-related quantities with input damping time
  subroutine init2_Dampfluc(this,taux,tauy,tauz)
    implicit none
    include "mpif.h"
    type (Dampfluc), intent(out) :: this
    double precision, intent(in) :: taux,tauy,tauz

    this%lambdax=exp(-1/taux)
    this%lambday=exp(-1/tauy)
    this%lambdaz=exp(-1/tauz)

  end subroutine init2_Dampfluc



  !in this subroutine, we use z, pz/pz0 instead of z/sigmaz, pz/sigmapz
  subroutine kickold_Dampfluc(Pts1in,nptlc,dampfluc1,sigma,t,myid,close2g,&
       iseedinp,icount)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc,myid
    integer, intent(inout) :: iseedinp
    integer*8, intent(inout) :: icount
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, dimension(6), intent(in) :: sigma
    double precision, dimension(2), intent(in) :: close2g
    double precision, intent(in) :: t
    type (Dampfluc), intent(in) :: dampfluc1
    double precision, dimension(6) :: tmp
    !double precision, dimension(6*nptlc) :: rannum6
    double precision, dimension(6) :: rannum6
    double precision :: sqrt3,tmpx,tmpy
    integer :: i,seedsize,ir
    !integer, allocatable, dimension(:) :: seedarray
    !double precision :: xcheck

    !        call random_seed(SIZE=seedsize)
    !        allocate(seedarray(seedsize))
    !        do i = 1, seedsize
    !          seedarray(i) = (1000+5*myid)*(myid+7)+i-1+t
    !        enddo
    !        call random_seed(PUT=seedarray)

    !        call random_number(xcheck)
    !        if(t.lt.10.0) then
    !        endif


    !        seedsize = 6*nptlc
    !        call ran2(iseedinp,rannum6,seedsize)
    !        icount = icount + seedsize
    !        call random_number(rannum6)

    seedsize = 6

    sqrt3 = sqrt(3.0d0)
    !        rannum6 = (2*rannum6-1.0d0)*sqrt3

    tmp(1) = dampfluc1%lambdax
    tmp(2) = sqrt(1.0 - (dampfluc1%lambdax)**2)
    tmp(3) = dampfluc1%lambday
    tmp(4) = sqrt(1.0 - (dampfluc1%lambday)**2)
    tmp(5) = dampfluc1%lambdaz
    tmp(6) = sqrt(1.0 - (dampfluc1%lambdaz)**2)

    !        if(t.lt.10.0) then
    !        endif


    do i = 1, nptlc
       ir = 0
       call ran2(iseedinp,rannum6,seedsize)
       rannum6 = (2*rannum6-1.0d0)*sqrt3
       ir = ir + 1
       tmpx = Pts1in(1,i) - close2g(1)
       Pts1in(1,i) = close2g(1)+tmp(1)*tmpx + tmp(2)*rannum6(ir)*sigma(1)
       !Pts1in(1,i) = close2g(1)+ tmp(2)*rannum6(ir)*sigma(1)
       !Pts1in(1,i) = close2g(1)+tmp(1)*tmpx 
       ir = ir + 1
       Pts1in(2,i) = tmp(1)*Pts1in(2,i) + tmp(2)*rannum6(ir)*sigma(2)
       !Pts1in(2,i) =  tmp(2)*rannum6(ir)*sigma(2)
       !Pts1in(2,i) = tmp(1)*Pts1in(2,i) 
       ir = ir + 1
       tmpy = Pts1in(3,i) - close2g(2)
       Pts1in(3,i) = close2g(2)+tmp(3)*tmpy + tmp(4)*rannum6(ir)*sigma(3)
       !Pts1in(3,i) = close2g(2)+ tmp(4)*rannum6(ir)*sigma(3)
       !Pts1in(3,i) = close2g(2)+tmp(3)*tmpy 
       ir = ir + 1
       Pts1in(4,i) = tmp(3)*Pts1in(4,i) + tmp(4)*rannum6(ir)*sigma(4)
       !Pts1in(4,i) = tmp(4)*rannum6(ir)*sigma(4)
       !Pts1in(4,i) = tmp(3)*Pts1in(4,i) 
       !the following has to be double checked for the definition of
       !z and pz
       ir = ir + 1
       Pts1in(5,i) = tmp(5)*Pts1in(5,i) + tmp(6)*rannum6(ir)*sigma(5)
       !Pts1in(5,i) = tmp(5)*Pts1in(5,i) + tmp(6)*rannum6(ir)
       ir = ir + 1
       Pts1in(6,i) = tmp(5)*Pts1in(6,i) + tmp(6)*rannum6(ir)*sigma(6)
       !Pts1in(6,i) = tmp(5)*Pts1in(6,i) + tmp(6)*rannum6(ir)
    enddo
    icount = icount + nptlc

    !        deallocate(seedarray)

  end subroutine kickold_Dampfluc



  !//in this subroutine, we use z, pz/pz0 instead of z/sigmaz, pz/sigmapz
  !//we need to go into the system with alpha = 0.
  !//(X,Px)' = (1 0 // alpha beta)(x,px)'
  subroutine kick_Dampfluc(Pts1in,nptlc,dampfluc1,sigma,t,myid,close2g,&
       ax,bx,ay,by)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc,myid
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, dimension(6), intent(in) :: sigma
    double precision, dimension(2), intent(in) :: close2g
    double precision, intent(in) :: t,ax,bx,ay,by
    type (Dampfluc), intent(in) :: dampfluc1
    double precision, dimension(6) :: tmp
    double precision, dimension(6*nptlc) :: rannum6
    double precision :: sqrt3,tmpx,tmpy,tmppx,tmppy,epsx,epsy,&
         a11x,a21x,a22x,a11y,a21y,a22y
    integer :: i,ir
    !integer, allocatable, dimension(:) :: seedarray
    !integer :: seedsize
    double precision :: xcheck

    !        call random_seed(SIZE=seedsize)
    !        allocate(seedarray(seedsize))
    !        do i = 1, seedsize
    !          seedarray(i) = (1000+5*myid)*(myid+7)+i-1+t
    !        enddo
    !        call random_seed(PUT=seedarray)

    call random_number(xcheck)
    call random_number(rannum6)
    sqrt3 = sqrt(3.0)
    rannum6 = (2*rannum6-1)*sqrt3

    epsx = sigma(1)*sigma(2)/sqrt(1.0+ax*ax)
    epsy = sigma(3)*sigma(4)/sqrt(1.0+ay*ay)
    tmp(1) = dampfluc1%lambdax
    tmp(2) = sqrt(epsx)*sqrt(1.0 - (dampfluc1%lambdax)**2)
    tmp(3) = dampfluc1%lambday
    tmp(4) = sqrt(epsy)*sqrt(1.0 - (dampfluc1%lambday)**2)
    tmp(5) = dampfluc1%lambdaz
    tmp(6) = sqrt(1.0 - (dampfluc1%lambdaz)**2)

    !//prepare for transformation matrix A
    a11x = sqrt(bx)
    a21x = -ax/sqrt(bx)
    a22x = 1.0/sqrt(bx)
    a11y = sqrt(by)
    a21y = -ay/sqrt(by)
    a22y = 1.0/sqrt(by)

    ir = 0
    do i = 1, nptlc
       !//radation in x-px
       tmpx = (Pts1in(1,i)-close2g(1))*a22x
       tmppx = a11x*Pts1in(2,i) - (Pts1in(1,i)-close2g(1))*a21x
       ir = ir + 1
       tmpx = tmp(1)*tmpx + tmp(2)*rannum6(ir)
       ir = ir + 1
       tmppx = tmp(1)*tmppx + tmp(2)*rannum6(ir)
       Pts1in(1,i) = close2g(1)+tmpx*a11x
       Pts1in(2,i) = tmpx*a21x + tmppx*a22x 

       !//radation in y-py
       tmpy = (Pts1in(3,i)-close2g(2))*a22y
       tmppy = a11y*Pts1in(4,i) - (Pts1in(3,i)-close2g(2))*a21y
       ir = ir + 1
       tmpy = tmp(3)*tmpy + tmp(4)*rannum6(ir)
       ir = ir + 1
       tmppy = tmp(3)*tmppy + tmp(4)*rannum6(ir)
       Pts1in(3,i) = close2g(2)+tmpy*a11y
       Pts1in(4,i) = tmpy*a21y + tmppy*a22y 

       !the following has to be double checked for the definition of
       !z and pz
       ir = ir + 1
       Pts1in(5,i) = tmp(5)*Pts1in(5,i) + tmp(6)*rannum6(ir)*sigma(5)
       ir = ir + 1
       Pts1in(6,i) = tmp(5)*Pts1in(6,i) + tmp(6)*rannum6(ir)*sigma(6)
    enddo

    !        deallocate(seedarray)

  end subroutine kick_Dampfluc



end module Dampflucclass
