module Extmap2ndclass

  type Extmap2nd
     double precision, dimension(6,6) :: linearmap
  end type Extmap2nd
  interface kick_Extmap2nd
     module procedure kick1_Extmap2nd,kick2_Extmap2nd
  end interface


contains

  subroutine construct_Extmap2nd(this,tunex,tuney,tunez,alphax,betax,alphay,betay)
    implicit none
    include "mpif.h"
    type (Extmap2nd), intent(out) :: this
    double precision, intent(in) :: tunex,tuney,tunez,alphax,betax,alphay,betay
    double precision :: twopi,cx,sx,gx,cy,sy,gy,cz,sz

    this%linearmap = 0.0d0
    twopi = 4*dasin(1.0d0)
    gx = (1+alphax*alphax)/betax
    cx=dcos(twopi*tunex)
    sx=dsin(twopi*tunex)
    this%linearmap(1,1) = cx + alphax*sx
    this%linearmap(2,1) = -gx*sx
    this%linearmap(1,2) = betax*sx
    this%linearmap(2,2) = cx - alphax*sx
    gy = (1+alphay*alphay)/betay
    cy=dcos(twopi*tuney)
    sy=dsin(twopi*tuney)
    this%linearmap(3,3) = cy + alphay*sy
    this%linearmap(4,3) = -gy*sy
    this%linearmap(3,4) = betay*sy
    this%linearmap(4,4) = cy - alphay*sy
    cz=dcos(twopi*tunez)
    sz=dsin(twopi*tunez)
    this%linearmap(5,5) = cz
    this%linearmap(6,5) = -sz
    this%linearmap(5,6) = sz
    this%linearmap(6,6) = cz

  end subroutine construct_Extmap2nd



  !apply linear map kick
  !in this subroutine, we use z, pz/pz0 instead of z/sigmaz, pz/sigmapz
  subroutine kick1_Extmap2nd(Pts1in,nptlc,this,close2g,szspz)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, dimension(2) :: close2g
    double precision :: szspz 
    type (Extmap2nd), intent(in) :: this
    integer :: i,j,k
    double precision, dimension(6) :: tmp

    do i = 1, nptlc
       tmp(1) = Pts1in(1,i) - close2g(1)
       tmp(2) = Pts1in(2,i)
       tmp(3) = Pts1in(3,i) - close2g(2)
       tmp(4) = Pts1in(4,i)
       tmp(5) = Pts1in(5,i)
       tmp(6) = Pts1in(6,i)
       do j = 1, 6
          Pts1in(j,i) = 0.0
       enddo
       do k = 1, 6
          do j = 1, 6
             Pts1in(j,i) = Pts1in(j,i) + tmp(k)*this%linearmap(j,k)
          enddo
       enddo

       Pts1in(1,i) = Pts1in(1,i) + close2g(1) 
       Pts1in(3,i) = Pts1in(3,i) + close2g(2) 
       !use z, pz/pz0
       Pts1in(5,i) = this%linearmap(5,5)*tmp(5) + &
            this%linearmap(5,6)*tmp(6)*szspz
       Pts1in(6,i) = this%linearmap(6,5)*tmp(5)/szspz + &
            this%linearmap(6,6)*tmp(6)
    enddo

  end subroutine kick1_Extmap2nd



  !apply linear map kick
  !in this subroutine, we use z, pz/pz0 instead of z/sigmaz, pz/sigmapz
  subroutine kick2_Extmap2nd(Pts1in,nptlc,this,close2g)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, dimension(2) :: close2g
    type (Extmap2nd), intent(in) :: this
    integer :: i,j,k
    double precision, dimension(6) :: tmp

    do i = 1, nptlc
       tmp(1) = Pts1in(1,i) - close2g(1)
       tmp(2) = Pts1in(2,i)
       tmp(3) = Pts1in(3,i) - close2g(2)
       tmp(4) = Pts1in(4,i)
       tmp(5) = Pts1in(5,i)
       tmp(6) = Pts1in(6,i)
       do j = 1, 6
          Pts1in(j,i) = 0.0d0
       enddo
       do k = 1, 6
          do j = 1, 6
             Pts1in(j,i) = Pts1in(j,i) + tmp(k)*this%linearmap(j,k)
          enddo
       enddo

       Pts1in(1,i) = Pts1in(1,i) + close2g(1) 
       Pts1in(3,i) = Pts1in(3,i) + close2g(2) 
    enddo

  end subroutine kick2_Extmap2nd



  subroutine setij_Extmap2nd(this,i,j,value)
    implicit none
    include "mpif.h"
    type (Extmap2nd), intent(out) :: this
    double precision, intent(in) :: value
    integer, intent(in) :: i,j

    this%linearmap(i,j) = value

  end subroutine setij_Extmap2nd



  subroutine setmap_Extmap2nd(this,value)
    implicit none
    include "mpif.h"
    type (Extmap2nd), intent(out) :: this
    double precision, dimension(6,6), intent(in) :: value
    integer :: i,j

    do j = 1, 6
       do i = 1, 6
          this%linearmap(i,j) = value(i,j)
       enddo
    enddo

  end subroutine setmap_Extmap2nd



  !apply linear map kick
  !in this subroutine, we use z, pz/pz0 instead of z/sigmaz, pz/sigmapz
  subroutine kick1test_Extmap2nd(Pts1in,nptlc,this,close2g,szspz)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, dimension(2) :: close2g
    double precision :: szspz 
    type (Extmap2nd), intent(in) :: this
    integer :: i!,j,k
    double precision, dimension(6) :: tmp

    do i = 1, nptlc
       tmp(1) = Pts1in(1,i) - close2g(1)
       tmp(2) = Pts1in(2,i)
       tmp(3) = Pts1in(3,i) - close2g(2)
       tmp(4) = Pts1in(4,i)
       tmp(5) = Pts1in(5,i)
       tmp(6) = Pts1in(6,i)

       Pts1in(1,i) = this%linearmap(1,1)*tmp(1) + &
            this%linearmap(1,2)*tmp(2) + close2g(1)
       Pts1in(2,i) = this%linearmap(2,1)*tmp(1) + &
            this%linearmap(2,2)*tmp(2)
       Pts1in(3,i) = this%linearmap(3,3)*tmp(3) + &
            this%linearmap(3,4)*tmp(4) + close2g(2)
       Pts1in(4,i) = this%linearmap(4,3)*tmp(3) + &
            this%linearmap(4,4)*tmp(4)
       !use z, pz/pz0
       !Pts1in(5,i) = this%linearmap(5,5)*tmp(5) + &
       !              this%linearmap(5,6)*tmp(6)*szspz
       Pts1in(5,i) = this%linearmap(5,5)*tmp(5) + &
            this%linearmap(5,6)*tmp(6)
       !Pts1in(6,i) = this%linearmap(6,5)*tmp(5)/szspz + &
       !              this%linearmap(6,6)*tmp(6)
       Pts1in(6,i) = this%linearmap(6,5)*tmp(5) + &
            this%linearmap(6,6)*tmp(6)
    enddo

  end subroutine kick1test_Extmap2nd



end module Extmap2ndclass
