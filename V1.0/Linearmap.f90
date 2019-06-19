module Linearmapclass
  type Linearmap
     double precision, dimension(2,2) :: matx,maty,matz
  end type Linearmap

contains

  subroutine construct_Linearmap(this,tunex,tuney,tunez,alphax,betax,alphay,betay)
    implicit none
    include "mpif.h"
    type (Linearmap), intent(out) :: this
    double precision, intent(in) :: tunex,tuney,tunez,alphax,betax,alphay,betay
    double precision :: twopi,cx,sx,gx,cy,sy,gy,cz,sz

    twopi = 4*dasin(1.0d0)
    gx = (1+alphax*alphax)/betax
    cx=dcos(twopi*tunex)
    sx=dsin(twopi*tunex)
    this%matx(1,1) = cx + alphax*sx
    this%matx(2,1) = -gx*sx
    this%matx(1,2) = betax*sx
    this%matx(2,2) = cx - alphax*sx
    gy = (1+alphay*alphay)/betay
    cy=dcos(twopi*tuney)
    sy=dsin(twopi*tuney)
    this%maty(1,1) = cy + alphay*sy
    this%maty(2,1) = -gy*sy
    this%maty(1,2) = betay*sy
    this%maty(2,2) = cy - alphay*sy
    cz=dcos(twopi*tunez)
    sz=dsin(twopi*tunez)
    this%matz(1,1) = cz
    this%matz(2,1) = -sz
    this%matz(1,2) = sz
    this%matz(2,2) = cz

  end subroutine construct_Linearmap



  subroutine kickold_Linearmap(Pts1in,nptlc,linearmap1,close2g)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, dimension(2) :: close2g
    type (Linearmap), intent(in) :: linearmap1
    integer :: i
    double precision :: tmp1,tmp2
    double precision, dimension(2,2) :: matx1,maty1,matz1

    matx1 = linearmap1%matx
    maty1 = linearmap1%maty
    matz1 = linearmap1%matz

    do i = 1, nptlc
       tmp1 = Pts1in(1,i) - close2g(1)
       tmp2 = Pts1in(2,i)
       Pts1in(1,i) = close2g(1) + matx1(1,1)*tmp1+matx1(1,2)*tmp2
       Pts1in(2,i) = matx1(2,1)*tmp1+matx1(2,2)*tmp2
       tmp1 = Pts1in(3,i) - close2g(2)
       tmp2 = Pts1in(4,i)
       Pts1in(3,i) = close2g(2) + maty1(1,1)*tmp1+maty1(1,2)*tmp2
       Pts1in(4,i) = maty1(2,1)*tmp1+maty1(2,2)*tmp2
       tmp1 = Pts1in(5,i)
       tmp2 = Pts1in(6,i)
       Pts1in(5,i) = matz1(1,1)*tmp1+matz1(1,2)*tmp2
       Pts1in(6,i) = matz1(2,1)*tmp1+matz1(2,2)*tmp2
    enddo

  end subroutine kickold_Linearmap



  !in this subroutine, we use z, pz/pz0 instead of z/sigmaz, pz/sigmapz
  subroutine kick_Linearmap(Pts1in,nptlc,linearmap1,close2g,szspz)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, dimension(2) :: close2g
    double precision :: szspz 
    type (Linearmap), intent(in) :: linearmap1
    integer :: i
    double precision :: tmp1,tmp2
    double precision, dimension(2,2) :: matx1,maty1,matz1

    matx1 = linearmap1%matx
    maty1 = linearmap1%maty
    matz1 = linearmap1%matz

    do i = 1, nptlc
       tmp1 = Pts1in(1,i) - close2g(1)
       tmp2 = Pts1in(2,i)
       Pts1in(1,i) = close2g(1) + matx1(1,1)*tmp1+matx1(1,2)*tmp2
       Pts1in(2,i) = matx1(2,1)*tmp1+matx1(2,2)*tmp2
       tmp1 = Pts1in(3,i) - close2g(2)
       tmp2 = Pts1in(4,i)
       Pts1in(3,i) = close2g(2) + maty1(1,1)*tmp1+maty1(1,2)*tmp2
       Pts1in(4,i) = maty1(2,1)*tmp1+maty1(2,2)*tmp2
       tmp1 = Pts1in(5,i)
       tmp2 = Pts1in(6,i)
       !use z, pz/pz0
       Pts1in(5,i) = matz1(1,1)*tmp1+matz1(1,2)*tmp2*szspz
       Pts1in(6,i) = matz1(2,1)*tmp1/szspz+matz1(2,2)*tmp2
    enddo

  end subroutine kick_Linearmap



  !in this subroutine, we use z, pz/pz0 instead of z/sigmaz, pz/sigmapz
  subroutine kicktranv_Linearmap(Pts1in,nptlc,linearmap1,close2g,szspz)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, dimension(2) :: close2g
    double precision :: szspz 
    type (Linearmap), intent(in) :: linearmap1
    integer :: i
    double precision :: tmp1,tmp2
    double precision, dimension(2,2) :: matx1,maty1,matz1

    matx1 = linearmap1%matx
    maty1 = linearmap1%maty
    matz1 = linearmap1%matz

    do i = 1, nptlc
       tmp1 = Pts1in(1,i) - close2g(1)
       tmp2 = Pts1in(2,i)
       Pts1in(1,i) = close2g(1) + matx1(1,1)*tmp1+matx1(1,2)*tmp2
       Pts1in(2,i) = matx1(2,1)*tmp1+matx1(2,2)*tmp2
       tmp1 = Pts1in(3,i) - close2g(2)
       tmp2 = Pts1in(4,i)
       Pts1in(3,i) = close2g(2) + maty1(1,1)*tmp1+maty1(1,2)*tmp2
       Pts1in(4,i) = maty1(2,1)*tmp1+maty1(2,2)*tmp2
    enddo

  end subroutine kicktranv_Linearmap



  !in this subroutine, we use z, pz/pz0 instead of z/sigmaz, pz/sigmapz
  !subroutine kickLong_Linearmap(Pts1in,nptlc,linearmap1,close2g,szspz)
  subroutine kickLong_Linearmap(Pts1in,nptlc,tunez,close2g,szspz)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, dimension(2) :: close2g
    double precision :: szspz,tunez
    !type (Linearmap), intent(in) :: linearmap1
    integer :: i
    double precision :: tmp1,tmp2
    double precision, dimension(2,2) :: matz1 !matx1,maty1
    real*8 :: twopi,cz,sz

    !matx1 = linearmap1%matx
    !maty1 = linearmap1%maty
    !matz1 = linearmap1%matz
    twopi = 4*asin(1.0d0)
    cz=dcos(twopi*tunez)
    sz=dsin(twopi*tunez)
    matz1(1,1) = cz
    matz1(2,1) = -sz
    matz1(1,2) = sz
    matz1(2,2) = cz

    do i = 1, nptlc
       !tmp1 = Pts1in(1,i) - close2g(1)
       !tmp2 = Pts1in(2,i)
       !Pts1in(1,i) = close2g(1) + matx1(1,1)*tmp1+matx1(1,2)*tmp2
       !Pts1in(2,i) = matx1(2,1)*tmp1+matx1(2,2)*tmp2
       !tmp1 = Pts1in(3,i) - close2g(2)
       !tmp2 = Pts1in(4,i)
       !Pts1in(3,i) = close2g(2) + maty1(1,1)*tmp1+maty1(1,2)*tmp2
       !Pts1in(4,i) = maty1(2,1)*tmp1+maty1(2,2)*tmp2
       tmp1 = Pts1in(5,i)
       tmp2 = Pts1in(6,i)
       !use z, pz/pz0
       Pts1in(5,i) = matz1(1,1)*tmp1+matz1(1,2)*tmp2*szspz
       Pts1in(6,i) = matz1(2,1)*tmp1/szspz+matz1(2,2)*tmp2
    enddo

  end subroutine kickLong_Linearmap


end module Linearmapclass
