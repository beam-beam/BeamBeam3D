!3rd order transfer map (4th order Hamiltonia)
module Extmap4thclass
  type Extmap4th
     double precision, dimension(6,6,6,6) :: map4th
  end type Extmap4th
contains
  subroutine construct_Extmap4th(this)
    implicit none
    include "mpif.h"
    type (Extmap4th), intent(out) :: this
    integer :: i,j,k,l

    do l = 1, 6
       do k = 1, 6
          do j = 1, 6
             do i = 1, 6
                this%map4th(i,j,k,l) = 0.0d0
             enddo
          enddo
       enddo
    enddo

  end subroutine construct_Extmap4th



  !apply linear map kick
  !in this subroutine, we use z, pz/pz0 instead of z/sigmaz, pz/sigmapz
  subroutine kick_Extmap4th(Pts1in,nptlc,this,close2g)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, dimension(2), intent(in) :: close2g
    type (Extmap4th), intent(in) :: this
    integer :: i,j,k,l,m
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
       do m = 1, 6
          do l = 1, 6
             do k = 1, 6
                do j = 1, 6
                   Pts1in(j,i) = Pts1in(j,i) + tmp(k)*tmp(l)*tmp(m)*this%map4th(j,k,l,m)
                enddo
             enddo
          enddo
       enddo

       Pts1in(1,i) = Pts1in(1,i) + close2g(1) 
       Pts1in(3,i) = Pts1in(3,i) + close2g(2) 
    enddo

  end subroutine kick_Extmap4th



  subroutine setijk_Extmap4th(this,i,j,k,l,value)
    implicit none
    include "mpif.h"
    type (Extmap4th), intent(out) :: this
    double precision, intent(in) :: value
    integer, intent(in) :: i,j,k,l

    this%map4th(i,j,k,l) = value

  end subroutine setijk_Extmap4th



  subroutine setmap_Extmap4th(this,value)
    implicit none
    include "mpif.h"
    type (Extmap4th), intent(out) :: this
    double precision, dimension(6,6,6,6), intent(in) :: value
    integer :: i,j,k,l

    do l = 1, 6
       do k = 1, 6
          do j = 1, 6
             do i = 1, 6
                this%map4th(i,j,k,l) = value(i,j,k,l)
             enddo
          enddo
       enddo
    enddo

  end subroutine setmap_Extmap4th



end module Extmap4thclass
