module Extmap3rdclass
  type Extmap3rd
     double precision, dimension(6,6,6) :: map3rd
  end type Extmap3rd
contains
  subroutine construct_Extmap3rd(this)
    implicit none
    include "mpif.h"
    type (Extmap3rd), intent(out) :: this
    integer :: i,j,k

    do k = 1, 6
       do j = 1, 6
          do i = 1, 6
             this%map3rd(i,j,k) = 0.0
          enddo
       enddo
    enddo

  end subroutine construct_Extmap3rd



  !apply linear map kick
  !in this subroutine, we use z, pz/pz0 instead of z/sigmaz, pz/sigmapz
  subroutine kick_Extmap3rd(Pts1in,nptlc,this,close2g)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: nptlc
    double precision, pointer, dimension(:,:) :: Pts1in
    double precision, dimension(2), intent(in) :: close2g
    type (Extmap3rd), intent(in) :: this
    integer :: i,j,k,l
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
       do l = 1, 6
          do k = 1, 6
             do j = 1, 6
                Pts1in(j,i) = Pts1in(j,i) + tmp(k)*tmp(l)*this%map3rd(j,k,l)
             enddo
          enddo
       enddo

       Pts1in(1,i) = Pts1in(1,i) + close2g(1) 
       Pts1in(3,i) = Pts1in(3,i) + close2g(2) 
    enddo

  end subroutine kick_Extmap3rd



  subroutine setijk_Extmap3rd(this,i,j,k,value)
    implicit none
    include "mpif.h"
    type (Extmap3rd), intent(out) :: this
    double precision, intent(in) :: value
    integer, intent(in) :: i,j,k

    this%map3rd(i,j,k) = value

  end subroutine setijk_Extmap3rd



  subroutine setmap_Extmap3rd(this,value)
    implicit none
    include "mpif.h"
    type (Extmap3rd), intent(out) :: this
    double precision, dimension(6,6,6), intent(in) :: value
    integer :: i,j,k

    do k = 1, 6
       do j = 1, 6
          do i = 1, 6
             this%map3rd(i,j,k) = value(i,j,k)
          enddo
       enddo
    enddo

  end subroutine setmap_Extmap3rd



end module Extmap3rdclass
