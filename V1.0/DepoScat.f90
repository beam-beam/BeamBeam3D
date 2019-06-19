module DepoScatclass
  use Timerclass

contains

  ! deposit particles onto grid.
  subroutine deposit2d(maxnpt,innp,innx,inny,hx,hy, &
       xmin,ymin,rays,rho)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,inny,maxnpt
    double precision, intent (in), dimension (6, maxnpt) :: rays
    double precision, intent (out), dimension (innx,inny) :: rho
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1
    double precision :: ab,cd
    integer :: i !, j
    double precision :: t0
    double precision :: hxi,hyi

    call starttime_Timer( t0 )

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    rho=0.d0
    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=((xmin-rays(1,i))+ix*hx)*hxi
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=((ymin-rays(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.inny)) &
       !then
       !   innp
       !   stop
       !endif

       ! (i,j):
       rho(ix,jx) = rho(ix,jx) + ab*cd
       ! (i,j+1):
       rho(ix,jx1) = rho(ix,jx1) + ab*(1.0d0-cd)
       ! (i+1,j):
       rho(ix1,jx) = rho(ix1,jx)+(1.0d0-ab)*cd
       ! (i+1,j+1):
       rho(ix1,jx1) = rho(ix1,jx1)+(1.0d0-ab)*(1.0d0-cd)
    enddo

    t_rhofas = t_rhofas + elapsedtime_Timer( t0 )

  end subroutine deposit2d



  !interpolate the field from grid onto particles.
  !linear interpolation in x and y too. 
  !take out the dipole force of the offset beams
  subroutine scatter2dAorg(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,weight,&
       coef,nytot,myidy,npyhalf)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,nytot,maxnpt,myidy,npyhalf
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(in), dimension(innx,nytot) :: phi
    double precision, intent(in) :: weight,coef
    double precision, dimension(innx,nytot) :: egx,egy
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1
    integer :: i, j
    double precision :: ab,cd
    double precision :: hxi,hyi,eyn,exn,ex0,tmpx,tmpy

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    !here, the beam 2 is located at +0.15 sigma away from x-axis
    !we need to substract the dipole kick which already been taken into account when the
    !the close orbit is calculated.
    if(myidy.lt.npyhalf) then
       ex0 = -(1.0d0-exp(-0.15*0.15/2.0))/(0.15*1.5899686e-5)
       !!          !beam 2 has 10% rms size different.
       !          ex0 = -(1.0d0-exp(-0.15*0.15/2.0))/(0.15*1.43097174e-5)
    else
       ex0 = (1.0d0-exp(-0.15*0.15/2.0))/(0.15*1.5899686e-5)
    endif
    !       nominal case
    ex0 = 0.0d0

    !Ex
    egx = 0.0
    do j = 1, nytot
       !egx(1,j) = hxi*(phi(1,j)-phi(2,j))
       egx(1,j) = hxi*(1.5d0*phi(1,j)-2*phi(2,j)+0.5d0*phi(3,j)) - ex0
       do i = 2, innx-1
          egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j)) - ex0
       enddo
       !egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
       egx(innx,j) = hxi*(-0.5d0*phi(innx-2,j)+2*phi(innx-1,j)-1.5d0*phi(innx,j)) - ex0
    enddo

    !Ey
    egy = 0.0
    do i = 1, innx
       !egy(i,1) = hyi*(phi(i,1)-phi(i,2))
       egy(i,1) = hyi*(1.5d0*phi(i,1)-2*phi(i,2)+0.5d0*phi(i,3))
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1))
       enddo
    enddo
    do i = 1, innx
       !egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
       egy(i,nytot) = hyi*(-0.5d0*phi(i,nytot-2)+2*phi(i,nytot-1)-1.5d0*phi(i,nytot))
    enddo

    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=((xmin-rays(1,i))+ix*hx)*hxi
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=((ymin-rays(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.nytot)) &
       !then
       !                          innp
       !   stop
       !endif
       exn = egx(ix,jx)*ab*cd+egx(ix,jx1)*ab*(1.0d0-cd) &
            +egx(ix1,jx)*(1.0d0-ab)*cd+egx(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       eyn = egy(ix,jx)*ab*cd+egy(ix,jx1)*ab*(1.0d0-cd) &
            +egy(ix1,jx)*(1.0d0-ab)*cd+egy(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4
    enddo

  end subroutine scatter2dAorg



  !interpolate the field from grid onto particles.
  !linear interpolation in x and y too.
  subroutine scatter2d(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,weight,&
       coef,nytot)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,nytot,maxnpt
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(in), dimension(innx,nytot) :: phi
    double precision, intent(in) :: weight,coef
    double precision, dimension(innx,nytot) :: egx,egy
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1
    integer :: i, j
    double precision :: ab,cd
    double precision :: hxi,hyi,eyn,exn,tmpx,tmpy

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    !Ex
    egx = 0.0
    do j = 1, nytot
       !egx(1,j) = hxi*(phi(1,j)-phi(2,j))
       egx(1,j) = hxi*(1.5d0*phi(1,j)-2*phi(2,j)+0.5d0*phi(3,j))
       do i = 2, innx-1
          egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j))
       enddo
       !egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
       egx(innx,j) = hxi*(-0.5d0*phi(innx-2,j)+2*phi(innx-1,j)-1.5d0*phi(innx,j))
    enddo

    !Ey
    egy = 0.0
    do i = 1, innx
       !egy(i,1) = hyi*(phi(i,1)-phi(i,2))
       egy(i,1) = hyi*(1.5d0*phi(i,1)-2*phi(i,2)+0.5d0*phi(i,3))
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1))
       enddo
    enddo
    do i = 1, innx
       !egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
       egy(i,nytot) = hyi*(-0.5d0*phi(i,nytot-2)+2*phi(i,nytot-1)-1.5d0*phi(i,nytot))
    enddo

    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=((xmin-rays(1,i))+ix*hx)*hxi
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=((ymin-rays(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.nytot)) &
       !then
       !                          innp
       !   stop
       !endif
       exn = egx(ix,jx)*ab*cd+egx(ix,jx1)*ab*(1.0d0-cd) &
            +egx(ix1,jx)*(1.0d0-ab)*cd+egx(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       eyn = egy(ix,jx)*ab*cd+egy(ix,jx1)*ab*(1.0d0-cd) &
            +egy(ix1,jx)*(1.0d0-ab)*cd+egy(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4
    enddo

  end subroutine scatter2d



  !interpolate the Ex and Ey field from grid onto particles.
  !linear interpolation in x and y too.
  subroutine scatter2dExEy(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,weight,&
       coef,nytot,egx,egy)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,nytot,maxnpt
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(in) :: weight,coef
    double precision, dimension(innx,nytot) :: egx,egy
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1
    integer :: i
    double precision :: ab,cd
    double precision :: hxi,hyi,eyn,exn,tmpx,tmpy

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=((xmin-rays(1,i))+ix*hx)*hxi
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=((ymin-rays(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.nytot)) &
       !then
       !                          innp
       !   stop
       !endif
       exn = egx(ix,jx)*ab*cd+egx(ix,jx1)*ab*(1.0d0-cd) &
            +egx(ix1,jx)*(1.0d0-ab)*cd+egx(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       eyn = egy(ix,jx)*ab*cd+egy(ix,jx1)*ab*(1.0d0-cd) &
            +egy(ix1,jx)*(1.0d0-ab)*cd+egy(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       !rays(2,i) = rays(2,i) + exn*weight*coef
       !rays(4,i) = rays(4,i) + eyn*weight*coef
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4
    enddo

  end subroutine scatter2dExEy



  subroutine scatter2dnewold(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,phi2,weight,&
       coef,nytot,zbk,hzi1)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,nytot,maxnpt
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(inout), dimension(innx,nytot) :: phi,phi2
    double precision, intent(in) :: weight,coef,zbk,hzi1
    double precision, dimension(innx,nytot) :: egx,egy,egx2,egy2
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1
    integer :: i, j
    double precision :: ab,cd
    double precision :: hxi,hyi,eyn,exn,ss,tmpx,tmpy !,distance

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    !Ex
    egx = 0.0
    do j = 1, nytot
       !egx(1,j) = hxi*(phi(1,j)-phi(2,j))
       egx(1,j) = hxi*(1.5d0*phi(1,j)-2*phi(2,j)+0.5d0*phi(3,j))
       do i = 2, innx-1
          egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j))
       enddo
       !egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
       egx(innx,j) = hxi*(-0.5d0*phi(innx-2,j)+2*phi(innx-1,j)-1.5d0*phi(innx,j))
    enddo

    !Ey
    egy = 0.0
    do i = 1, innx
       !egy(i,1) = hyi*(phi(i,1)-phi(i,2))
       egy(i,1) = hyi*(1.5d0*phi(i,1)-2*phi(i,2)+0.5d0*phi(i,3))
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1))
       enddo
    enddo
    do i = 1, innx
       !egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
       egy(i,nytot) = hyi*(-0.5d0*phi(i,nytot-2)+2*phi(i,nytot-1)-1.5d0*phi(i,nytot))
    enddo

    !Ex2
    egx2 = 0.0
    do j = 1, nytot
       !egx2(1,j) = hxi*(phi2(1,j)-phi2(2,j))
       egx2(1,j) = hxi*(1.5d0*phi2(1,j)-2*phi2(2,j)+0.5d0*phi2(3,j))
       do i = 2, innx-1
          egx2(i,j) = 0.5d0*hxi*(phi2(i-1,j)-phi2(i+1,j))
       enddo
       !egx2(innx,j) = hxi*(phi2(innx-1,j)-phi2(innx,j))
       egx2(innx,j) = hxi*(-0.5d0*phi2(innx-2,j)+2*phi2(innx-1,j)-1.5d0*phi2(innx,j))
    enddo

    !Ey2
    egy2 = 0.0
    do i = 1, innx
       !egy2(i,1) = hyi*(phi2(i,1)-phi2(i,2))
       egy2(i,1) = hyi*(1.5d0*phi2(i,1)-2*phi2(i,2)+0.5d0*phi2(i,3))
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy2(i,j) = 0.5d0*hyi*(phi2(i,j-1)-phi2(i,j+1))
       enddo
    enddo
    do i = 1, innx
       !egy2(i,nytot) = hyi*(phi2(i,nytot-1)-phi2(i,nytot))
       egy2(i,nytot) = hyi*(-0.5d0*phi2(i,nytot-2)+2*phi2(i,nytot-1)-1.5d0*phi2(i,nytot))
    enddo

    !here, phi and phi2 have been used to store the E field
    do j = 1, nytot
       do i = 1, innx
          phi(i,j) = (egx2(i,j)-egx(i,j))/hzi1
          phi2(i,j) = (egy2(i,j)-egy(i,j))/hzi1
       enddo
    enddo

    do i = 1, innp
       !          distance = 0.0
       !          if(i.lt.5) then
       !          endif
       !          tmpx = rays(1,i)+distance*rays(2,i)
       !          ix=(tmpx-xmin)*hxi + 1
       !          ab=((xmin-tmpx)+ix*hx)*hxi
       !          tmpy = rays(3,i)+distance*rays(4,i)
       !          jx=(tmpy-ymin)*hyi + 1
       !          cd=((ymin-tmpy)+jx*hy)*hyi
       !          ix1=ix+1
       !          jx1=jx+1

       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=((xmin-rays(1,i))+ix*hx)*hxi
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=((ymin-rays(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.nytot)) &
       !then
       !                          innp
       !   stop
       !endif
       ss = rays(5,i) - zbk
       !          ss = 0.0
       !          ss = zthk
       exn = (egx(ix,jx)+ss*phi(ix,jx))*ab*cd+(egx(ix,jx1)+ss*phi(ix,jx1))*&
            ab*(1.0d0-cd)+(egx(ix1,jx)+ss*phi(ix1,jx))*(1.0d0-ab)*cd+&
            (egx(ix1,jx1)+ss*phi(ix1,jx1))*(1.0d0-ab)*(1.0d0-cd)
       eyn = (egy(ix,jx)+ss*phi2(ix,jx))*ab*cd+(egy(ix,jx1)+ss*phi2(ix,jx1))*&
            ab*(1.0d0-cd)+(egy(ix1,jx)+ss*phi2(ix1,jx))*(1.0d0-ab)*cd+&
            (egy(ix1,jx1)+ss*phi2(ix1,jx1))*(1.0d0-ab)*(1.0d0-cd)
       !rays(2,i) = rays(2,i) + exn*weight*coef
       !rays(4,i) = rays(4,i) + eyn*weight*coef
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4
    enddo

  end subroutine scatter2dnewold



  !including the longitudinal Ez
  subroutine scatter2dnew(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,phi2,weight,&
       coef,nytot,zbk,hzi1)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,nytot,maxnpt
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(inout), dimension(innx,nytot) :: phi,phi2
    double precision, intent(in) :: weight,coef,zbk,hzi1
    double precision, dimension(innx,nytot) :: egx,egy,egx2,egy2,egz
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1
    integer :: i, j
    double precision :: ab,cd
    double precision :: hxi,hyi,eyn,exn,ss,tmpx,tmpy,ezn !,distance

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    !Ez
    do j = 1, nytot
       do i = 1, innx
          egz(i,j) = (phi(i,j) - phi2(i,j))/hzi1 
       enddo
    enddo
    !Ex
    egx = 0.0
    do j = 1, nytot
       !egx(1,j) = hxi*(phi(1,j)-phi(2,j))
       egx(1,j) = hxi*(1.5d0*phi(1,j)-2*phi(2,j)+0.5d0*phi(3,j))
       do i = 2, innx-1
          egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j))
       enddo
       !egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
       egx(innx,j) = hxi*(-0.5d0*phi(innx-2,j)+2*phi(innx-1,j)-1.5d0*phi(innx,j))
    enddo

    !Ey
    egy = 0.0
    do i = 1, innx
       !egy(i,1) = hyi*(phi(i,1)-phi(i,2))
       egy(i,1) = hyi*(1.5d0*phi(i,1)-2*phi(i,2)+0.5d0*phi(i,3))
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1))
       enddo
    enddo
    do i = 1, innx
       !egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
       egy(i,nytot) = hyi*(-0.5d0*phi(i,nytot-2)+2*phi(i,nytot-1)-1.5d0*phi(i,nytot))
    enddo

    !Ex2
    egx2 = 0.0
    do j = 1, nytot
       !egx2(1,j) = hxi*(phi2(1,j)-phi2(2,j))
       egx2(1,j) = hxi*(1.5d0*phi2(1,j)-2*phi2(2,j)+0.5d0*phi2(3,j))
       do i = 2, innx-1
          egx2(i,j) = 0.5d0*hxi*(phi2(i-1,j)-phi2(i+1,j))
       enddo
       !egx2(innx,j) = hxi*(phi2(innx-1,j)-phi2(innx,j))
       egx2(innx,j) = hxi*(-0.5d0*phi2(innx-2,j)+2*phi2(innx-1,j)-1.5d0*phi2(innx,j))
    enddo

    !Ey2
    egy2 = 0.0
    do i = 1, innx
       !egy2(i,1) = hyi*(phi2(i,1)-phi2(i,2))
       egy2(i,1) = hyi*(1.5d0*phi2(i,1)-2*phi2(i,2)+0.5d0*phi2(i,3))
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy2(i,j) = 0.5d0*hyi*(phi2(i,j-1)-phi2(i,j+1))
       enddo
    enddo
    do i = 1, innx
       !egy2(i,nytot) = hyi*(phi2(i,nytot-1)-phi2(i,nytot))
       egy2(i,nytot) = hyi*(-0.5d0*phi2(i,nytot-2)+2*phi2(i,nytot-1)-1.5d0*phi2(i,nytot))
    enddo

    !here, phi and phi2 have been used to store the E field
    do j = 1, nytot
       do i = 1, innx
          phi(i,j) = (egx2(i,j)-egx(i,j))/hzi1
          phi2(i,j) = (egy2(i,j)-egy(i,j))/hzi1
       enddo
    enddo

    do i = 1, innp
       !          distance = 0.0
       !          if(i.lt.5) then
       !          endif
       !          tmpx = rays(1,i)+distance*rays(2,i)
       !          ix=(tmpx-xmin)*hxi + 1
       !          ab=((xmin-tmpx)+ix*hx)*hxi
       !          tmpy = rays(3,i)+distance*rays(4,i)
       !          jx=(tmpy-ymin)*hyi + 1
       !          cd=((ymin-tmpy)+jx*hy)*hyi
       !          ix1=ix+1
       !          jx1=jx+1

       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=((xmin-rays(1,i))+ix*hx)*hxi
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=((ymin-rays(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.nytot)) &
       !then
       !                          innp
       !   stop
       !endif
       ss = rays(5,i) - zbk
       !          ss = 0.0
       !          ss = zthk
       exn = (egx(ix,jx)+ss*phi(ix,jx))*ab*cd+(egx(ix,jx1)+ss*phi(ix,jx1))*&
            ab*(1.0d0-cd)+(egx(ix1,jx)+ss*phi(ix1,jx))*(1.0d0-ab)*cd+&
            (egx(ix1,jx1)+ss*phi(ix1,jx1))*(1.0d0-ab)*(1.0d0-cd)
       eyn = (egy(ix,jx)+ss*phi2(ix,jx))*ab*cd+(egy(ix,jx1)+ss*phi2(ix,jx1))*&
            ab*(1.0d0-cd)+(egy(ix1,jx)+ss*phi2(ix1,jx))*(1.0d0-ab)*cd+&
            (egy(ix1,jx1)+ss*phi2(ix1,jx1))*(1.0d0-ab)*(1.0d0-cd)
       ezn = egz(ix,jx)*ab*cd+egz(ix,jx1)*ab*(1.0d0-cd) &
            +egz(ix1,jx)*(1.0d0-ab)*cd+egz(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       tmpx = rays(2,i) + exn*weight*coef/2
       tmpy = rays(4,i) + eyn*weight*coef/2
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + 0.5d0*weight*coef*(exn*tmpx+eyn*tmpy+2*ezn)
    enddo

  end subroutine scatter2dnew



  !including the longitudinal Ez and substract offset
  subroutine scatter2dnewoffset(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,&
       phi,phi2,weight,coef,nytot,zbk,hzi1,shift,sigopp,myidy,npyhalf)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,nytot,maxnpt,myidy,npyhalf
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(inout), dimension(innx,nytot) :: phi,phi2
    double precision, intent(in) :: weight,coef,zbk,hzi1
    double precision, intent(in), dimension(2) :: shift
    double precision, dimension(innx,nytot) :: egx,egy,egx2,egy2,egz
    double precision :: xmin,ymin,hx,hy,sigopp
    integer :: ix,jx,ix1,jx1
    integer :: i, j
    double precision :: ab,cd
    double precision :: hxi,hyi,eyn,exn,ss,tmpx,tmpy,ezn !,distance
    real*8 :: ex0,ey0,rr2,sigavg2

    !the close orbit is calculated.
    rr2 = shift(1)**2+shift(2)**2
    sigavg2 = sigopp**2
    ex0 = 0.0d0
    ey0 = 0.0d0
    !here shift is defined as beam2-beam1
    if(myidy.lt.npyhalf) then
       ex0 = -shift(1)/rr2*(1.0d0-exp(-rr2/sigavg2/2.0d0))
       ey0 = -shift(2)/rr2*(1.0d0-exp(-rr2/sigavg2/2.0d0))
    else
       ex0 = shift(1)/rr2*(1.0d0-exp(-rr2/sigavg2/2.0d0))
       ey0 = shift(2)/rr2*(1.0d0-exp(-rr2/sigavg2/2.0d0))
    endif

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    !Ez
    do j = 1, nytot
       do i = 1, innx
          egz(i,j) = (phi(i,j) - phi2(i,j))/hzi1 
       enddo
    enddo
    !Ex
    egx = 0.0
    do j = 1, nytot
       !egx(1,j) = hxi*(phi(1,j)-phi(2,j))
       egx(1,j) = hxi*(1.5d0*phi(1,j)-2*phi(2,j)+0.5d0*phi(3,j)) -ex0
       do i = 2, innx-1
          egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j)) -ex0
       enddo
       !egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
       egx(innx,j) = hxi*(-0.5d0*phi(innx-2,j)+2*phi(innx-1,j)-1.5d0*phi(innx,j)) -ex0
    enddo

    !Ey
    egy = 0.0
    do i = 1, innx
       !egy(i,1) = hyi*(phi(i,1)-phi(i,2))
       egy(i,1) = hyi*(1.5d0*phi(i,1)-2*phi(i,2)+0.5d0*phi(i,3)) -ey0
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1)) -ey0
       enddo
    enddo
    do i = 1, innx
       !egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
       egy(i,nytot) = hyi*(-0.5d0*phi(i,nytot-2)+2*phi(i,nytot-1)-1.5d0*phi(i,nytot)) -ey0
    enddo

    !Ex2
    egx2 = 0.0
    do j = 1, nytot
       !egx2(1,j) = hxi*(phi2(1,j)-phi2(2,j))
       egx2(1,j) = hxi*(1.5d0*phi2(1,j)-2*phi2(2,j)+0.5d0*phi2(3,j)) -ex0
       do i = 2, innx-1
          egx2(i,j) = 0.5d0*hxi*(phi2(i-1,j)-phi2(i+1,j)) -ex0
       enddo
       !egx2(innx,j) = hxi*(phi2(innx-1,j)-phi2(innx,j))
       egx2(innx,j) = hxi*(-0.5d0*phi2(innx-2,j)+2*phi2(innx-1,j)-1.5d0*phi2(innx,j)) -ex0
    enddo

    !Ey2
    egy2 = 0.0
    do i = 1, innx
       !egy2(i,1) = hyi*(phi2(i,1)-phi2(i,2))
       egy2(i,1) = hyi*(1.5d0*phi2(i,1)-2*phi2(i,2)+0.5d0*phi2(i,3)) -ey0
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy2(i,j) = 0.5d0*hyi*(phi2(i,j-1)-phi2(i,j+1)) -ey0
       enddo
    enddo
    do i = 1, innx
       !egy2(i,nytot) = hyi*(phi2(i,nytot-1)-phi2(i,nytot))
       egy2(i,nytot) = hyi*(-0.5d0*phi2(i,nytot-2)+2*phi2(i,nytot-1)-1.5d0*phi2(i,nytot)) -ey0
    enddo

    !here, phi and phi2 have been used to store the E field
    do j = 1, nytot
       do i = 1, innx
          phi(i,j) = (egx2(i,j)-egx(i,j))/hzi1
          phi2(i,j) = (egy2(i,j)-egy(i,j))/hzi1
       enddo
    enddo

    do i = 1, innp
       !          distance = 0.0
       !          if(i.lt.5) then
       !          endif
       !          tmpx = rays(1,i)+distance*rays(2,i)
       !          ix=(tmpx-xmin)*hxi + 1
       !          ab=((xmin-tmpx)+ix*hx)*hxi
       !          tmpy = rays(3,i)+distance*rays(4,i)
       !          jx=(tmpy-ymin)*hyi + 1
       !          cd=((ymin-tmpy)+jx*hy)*hyi
       !          ix1=ix+1
       !          jx1=jx+1

       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=((xmin-rays(1,i))+ix*hx)*hxi
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=((ymin-rays(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.nytot)) &
       !then
       !                          innp
       !   stop
       !endif
       ss = rays(5,i) - zbk
       !          ss = 0.0
       !          ss = zthk
       exn = (egx(ix,jx)+ss*phi(ix,jx))*ab*cd+(egx(ix,jx1)+ss*phi(ix,jx1))*&
            ab*(1.0d0-cd)+(egx(ix1,jx)+ss*phi(ix1,jx))*(1.0d0-ab)*cd+&
            (egx(ix1,jx1)+ss*phi(ix1,jx1))*(1.0d0-ab)*(1.0d0-cd)
       eyn = (egy(ix,jx)+ss*phi2(ix,jx))*ab*cd+(egy(ix,jx1)+ss*phi2(ix,jx1))*&
            ab*(1.0d0-cd)+(egy(ix1,jx)+ss*phi2(ix1,jx))*(1.0d0-ab)*cd+&
            (egy(ix1,jx1)+ss*phi2(ix1,jx1))*(1.0d0-ab)*(1.0d0-cd)
       ezn = egz(ix,jx)*ab*cd+egz(ix,jx1)*ab*(1.0d0-cd) &
            +egz(ix1,jx)*(1.0d0-ab)*cd+egz(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       tmpx = rays(2,i) + exn*weight*coef/2
       tmpy = rays(4,i) + eyn*weight*coef/2
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + 0.5d0*weight*coef*(exn*tmpx+eyn*tmpy+2*ezn)
    enddo

  end subroutine scatter2dnewoffset



  !interpolate the field onto particles. Here, we have used a linear
  !interpolation along z. We have found the potential at front and back
  !of the slice and do a linear interpolation when we think particle
  !drift to the center of the opposite slice.
  !linear interpolation in x and y too.
  subroutine scatter2d2z(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,phi2,weight,&
       coef,nytot,z2,zthk,sbk)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,nytot,maxnpt
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(inout), dimension(innx,nytot) :: phi,phi2
    double precision, intent(in) :: weight,coef,z2,zthk,sbk
    double precision, dimension(innx,nytot) :: egx,egy,egx2,egy2
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1
    integer :: i, j
    double precision :: ab,cd
    double precision :: hxi,hyi,eyn,exn,ss,distance,tmpx,tmpy

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    !Ex
    egx = 0.0
    do j = 1, nytot
       !egx(1,j) = hxi*(phi(1,j)-phi(2,j))
       egx(1,j) = hxi*(1.5d0*phi(1,j)-2*phi(2,j)+0.5d0*phi(3,j))
       do i = 2, innx-1
          egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j))
       enddo
       !egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
       egx(innx,j) = hxi*(-0.5d0*phi(innx-2,j)+2*phi(innx-1,j)-1.5d0*phi(innx,j))
    enddo

    !Ey
    egy = 0.0
    do i = 1, innx
       !egy(i,1) = hyi*(phi(i,1)-phi(i,2))
       egy(i,1) = hyi*(1.5d0*phi(i,1)-2*phi(i,2)+0.5d0*phi(i,3))
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1))
       enddo
    enddo
    do i = 1, innx
       !egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
       egy(i,nytot) = hyi*(-0.5d0*phi(i,nytot-2)+2*phi(i,nytot-1)-1.5d0*phi(i,nytot))
    enddo

    !Ex2
    egx2 = 0.0
    do j = 1, nytot
       !egx2(1,j) = hxi*(phi2(1,j)-phi2(2,j))
       egx2(1,j) = hxi*(1.5d0*phi2(1,j)-2*phi2(2,j)+0.5d0*phi2(3,j))
       do i = 2, innx-1
          egx2(i,j) = 0.5d0*hxi*(phi2(i-1,j)-phi2(i+1,j))
       enddo
       !egx2(innx,j) = hxi*(phi2(innx-1,j)-phi2(innx,j))
       egx2(innx,j) = hxi*(-0.5d0*phi2(innx-2,j)+2*phi2(innx-1,j)-1.5d0*phi2(innx,j))
    enddo

    !Ey2
    egy2 = 0.0
    do i = 1, innx
       !egy2(i,1) = hyi*(phi2(i,1)-phi2(i,2))
       egy2(i,1) = hyi*(1.5d0*phi2(i,1)-2*phi2(i,2)+0.5d0*phi2(i,3))
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy2(i,j) = 0.5d0*hyi*(phi2(i,j-1)-phi2(i,j+1))
       enddo
    enddo
    do i = 1, innx
       !egy2(i,nytot) = hyi*(phi2(i,nytot-1)-phi2(i,nytot))
       egy2(i,nytot) = hyi*(-0.5d0*phi2(i,nytot-2)+2*phi2(i,nytot-1)-1.5d0*phi2(i,nytot))
    enddo

    !here, phi and phi2 have been used to store the E field
    do j = 1, nytot
       do i = 1, innx
          phi(i,j) = (egx2(i,j)-egx(i,j))/zthk
          phi2(i,j) = (egy2(i,j)-egy(i,j))/zthk
       enddo
    enddo

    do i = 1, innp
       !          distance = 0.0
       !          if(i.lt.5) then
       !          endif
       !          tmpx = rays(1,i)+distance*rays(2,i)
       !          ix=(tmpx-xmin)*hxi + 1
       !          ab=((xmin-tmpx)+ix*hx)*hxi
       !          tmpy = rays(3,i)+distance*rays(4,i)
       !          jx=(tmpy-ymin)*hyi + 1
       !          cd=((ymin-tmpy)+jx*hy)*hyi
       !          ix1=ix+1
       !          jx1=jx+1

       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=((xmin-rays(1,i))+ix*hx)*hxi
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=((ymin-rays(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.nytot)) &
       !then
       !                          innp
       !   stop
       !endif
       distance = (rays(5,i)-z2)/2
       ss = distance - sbk
       !          ss = 0.0
       !          ss = zthk
       exn = (egx(ix,jx)+ss*phi(ix,jx))*ab*cd+(egx(ix,jx1)+ss*phi(ix,jx1))*&
            ab*(1.0d0-cd)+(egx(ix1,jx)+ss*phi(ix1,jx))*(1.0d0-ab)*cd+&
            (egx(ix1,jx1)+ss*phi(ix1,jx1))*(1.0d0-ab)*(1.0d0-cd)
       eyn = (egy(ix,jx)+ss*phi2(ix,jx))*ab*cd+(egy(ix,jx1)+ss*phi2(ix,jx1))*&
            ab*(1.0d0-cd)+(egy(ix1,jx)+ss*phi2(ix1,jx))*(1.0d0-ab)*cd+&
            (egy(ix1,jx1)+ss*phi2(ix1,jx1))*(1.0d0-ab)*(1.0d0-cd)
       !rays(2,i) = rays(2,i) + exn*weight*coef
       !rays(4,i) = rays(4,i) + eyn*weight*coef
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4
    enddo

  end subroutine scatter2d2z



  !interpolate the field onto particles. Here, we have used a linear
  !interpolation along z. We have found the potential at front and back
  !of the slice and do a linear interpolation when we think particle
  !drift to the center of the opposite slice.
  !TSC interpolation in x and y too.
  subroutine scatter2d2zTsc(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,phi2,weight,&
       coef,nytot,z2,zthk,sbk)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,nytot,maxnpt
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(inout), dimension(innx,nytot) :: phi,phi2
    double precision, intent(in) :: weight,coef,z2,zthk,sbk
    double precision, dimension(innx,nytot) :: egx,egy,egx2,egy2
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1,ix2,jx2
    integer :: i, j
    double precision :: wix,wix1,wix2,wjx,wjx1,wjx2,ab,cd
    double precision :: hxi,hyi,eyn,exn,ss,distance,tmpx,tmpy

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    !Ex
    egx = 0.0
    do j = 1, nytot
       !egx(1,j) = hxi*(phi(1,j)-phi(2,j))
       egx(1,j) = hxi*(1.5d0*phi(1,j)-2*phi(2,j)+0.5d0*phi(3,j))
       do i = 2, innx-1
          egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j))
       enddo
       !egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
       egx(innx,j) = hxi*(-0.5d0*phi(innx-2,j)+2*phi(innx-1,j)-1.5d0*phi(innx,j))
    enddo

    !Ey
    egy = 0.0
    do i = 1, innx
       !egy(i,1) = hyi*(phi(i,1)-phi(i,2))
       egy(i,1) = hyi*(1.5d0*phi(i,1)-2*phi(i,2)+0.5d0*phi(i,3))
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1))
       enddo
    enddo
    do i = 1, innx
       !egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
       egy(i,nytot) = hyi*(-0.5d0*phi(i,nytot-2)+2*phi(i,nytot-1)-1.5d0*phi(i,nytot))
    enddo

    !Ex2
    egx2 = 0.0
    do j = 1, nytot
       !egx2(1,j) = hxi*(phi2(1,j)-phi2(2,j))
       egx2(1,j) = hxi*(1.5d0*phi2(1,j)-2*phi2(2,j)+0.5d0*phi2(3,j))
       do i = 2, innx-1
          egx2(i,j) = 0.5d0*hxi*(phi2(i-1,j)-phi2(i+1,j))
       enddo
       !egx2(innx,j) = hxi*(phi2(innx-1,j)-phi2(innx,j))
       egx2(innx,j) = hxi*(-0.5d0*phi2(innx-2,j)+2*phi2(innx-1,j)-1.5d0*phi2(innx,j))
    enddo

    !Ey2
    egy2 = 0.0
    do i = 1, innx
       !egy2(i,1) = hyi*(phi2(i,1)-phi2(i,2))
       egy2(i,1) = hyi*(1.5d0*phi2(i,1)-2*phi2(i,2)+0.5d0*phi2(i,3))
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy2(i,j) = 0.5d0*hyi*(phi2(i,j-1)-phi2(i,j+1))
       enddo
    enddo
    do i = 1, innx
       !egy2(i,nytot) = hyi*(phi2(i,nytot-1)-phi2(i,nytot))
       egy2(i,nytot) = hyi*(-0.5d0*phi2(i,nytot-2)+2*phi2(i,nytot-1)-1.5d0*phi2(i,nytot))
    enddo

    !here, phi and phi2 have been used to store the E field
    do j = 1, nytot
       do i = 1, innx
          phi(i,j) = (egx2(i,j)-egx(i,j))/zthk
          phi2(i,j) = (egy2(i,j)-egy(i,j))/zthk
       enddo
    enddo

    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=(rays(1,i)-xmin-(ix-1)*hx)*hxi
       !I am not so sure that we should use "ge" or "gt"
       if(ab.le.0.5) then
          ix2 = ix - 1
          wix = 0.75-ab*ab
          wix1 = (0.5+ab)**2/2
          wix2 = (0.5-ab)**2/2
       else
          ix2 = ix + 2
          wix = (1.5d0-ab)**2/2
          wix1 = 0.75 - (1-ab)*(1-ab)
          wix2 = (ab-0.5)**2/2
       endif
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=(rays(3,i)-ymin-(jx-1)*hy)*hyi
       if(cd.le.0.5) then
          jx2 = jx - 1
          wjx = 0.75-cd*cd
          wjx1 = (0.5+cd)**2/2
          wjx2 = (0.5-cd)**2/2
       else
          jx2 = jx + 2
          wjx = (1.5d0-cd)**2/2
          wjx1 = 0.75 - (1-cd)*(1-cd)
          wjx2 = (cd-0.5)**2/2
       endif
       ix1=ix+1
       jx1=jx+1

       distance = (rays(5,i)-z2)/2
       ss = distance - sbk
       !          ss = 0.0
       !          ss = zthk

       exn = (egx(ix,jx)+ss*phi(ix,jx))*wix*wjx+&
            (egx(ix,jx1)+ss*phi(ix,jx1))*wix*wjx1+&
            (egx(ix1,jx)+ss*phi(ix1,jx))*wix1*wjx+&
            (egx(ix1,jx1)+ss*phi(ix1,jx1))*wix1*wjx1+&
            (egx(ix,jx2)+ss*phi(ix,jx2))*wix*wjx2+&
            (egx(ix1,jx2)+ss*phi(ix1,jx2))*wix1*wjx2+&
            (egx(ix2,jx2)+ss*phi(ix2,jx2))*wix2*wjx2+&
            (egx(ix2,jx)+ss*phi(ix2,jx))*wix2*wjx+&
            (egx(ix2,jx1)+ss*phi(ix2,jx1))*wix2*wjx1
       eyn = (egy(ix,jx)+ss*phi2(ix,jx))*wix*wjx+&
            (egy(ix,jx1)+ss*phi2(ix,jx1))*wix*wjx1+&
            (egy(ix1,jx)+ss*phi2(ix1,jx))*wix1*wjx+&
            (egy(ix1,jx1)+ss*phi2(ix1,jx1))*wix1*wjx1+&
            (egy(ix,jx2)+ss*phi2(ix,jx2))*wix*wjx2+&
            (egy(ix1,jx2)+ss*phi2(ix1,jx2))*wix1*wjx2+&
            (egy(ix2,jx2)+ss*phi2(ix2,jx2))*wix2*wjx2+&
            (egy(ix2,jx)+ss*phi2(ix2,jx))*wix2*wjx+&
            (egy(ix2,jx1)+ss*phi2(ix2,jx1))*wix2*wjx1

       !rays(2,i) = rays(2,i) + exn*weight*coef
       !rays(4,i) = rays(4,i) + eyn*weight*coef
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4
    enddo

  end subroutine scatter2d2zTsc



  ! deposit particles onto grid using quadratic form. 
  ! this is not TSC since we do not want to do correction
  ! However, this could cause discontinunity across boundary
  ! there should be 1 extra-grid on each side of domain.
  !subroutine deposit2d(maxnpt,innp,innx,inny,hx,hy, &
  subroutine deposit2dqd(maxnpt,innp,innx,inny,hx,hy, &
       xmin,ymin,rays,rho)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,inny,maxnpt
    double precision, intent (in), dimension (6, maxnpt) :: rays
    double precision, intent (out), dimension (innx,inny) :: rho
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1,ix2,jx2
    double precision :: wix,wix1,wix2,wjx,wjx1,wjx2,ab,cd
    integer :: i !, j
    double precision :: t0
    double precision :: hxi,hyi

    call starttime_Timer( t0 )

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    rho=0.
    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=(rays(1,i)-xmin-(ix-1)*hx)*hxi
       !I am not so sure that we should use "ge" or "gt"
       if(ab.le.0.5) then
          ix2 = ix - 1
          wix = 1.0d0-ab*ab
          wix1 = 1.0d0 + 1.5d0*(ab-1.0d0) + 0.5d0*(ab - 1.0d0)**2
          wix2 = 1.0d0 - 1.5d0*(ab+1.0d0) + 0.5d0*(ab + 1.0d0)**2
       else
          ix2 = ix + 2 
          wix = 1.0d0 - 1.5d0*ab + 0.5d0*ab**2
          wix1 = 1.0-(1.0d0-ab)*(1.0d0-ab)
          wix2 = 1.0d0 + 1.5d0*(ab-2) + 0.5d0*(ab-2)**2
       endif
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=(rays(3,i)-ymin-(jx-1)*hy)*hyi
       if(cd.le.0.5) then
          jx2 = jx - 1
          wjx = 1.0d0-cd*cd
          wjx1 = 1.0d0 + 1.5d0*(cd-1.0d0) + 0.5d0*(cd - 1.0d0)**2
          wjx2 = 1.0d0 - 1.5d0*(cd+1.0d0) + 0.5d0*(cd + 1.0d0)**2
       else
          jx2 = jx + 2 
          wjx = 1.0d0 - 1.5d0*cd + 0.5d0*cd**2
          wjx1 = 1.0d0-(1.0d0-cd)*(1.0d0-cd)
          wjx2 = 1.0d0 + 1.5d0*(cd-2) + 0.5d0*(cd-2)**2
       endif
       ix1=ix+1
       jx1=jx+1

       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.inny)) &
       !then
       !   innp
       !   stop
       !endif

       ! (i,j):
       rho(ix,jx) = rho(ix,jx) + wix*wjx
       ! (i,j+1):
       rho(ix,jx1) = rho(ix,jx1) + wix*wjx1
       ! (i+1,j):
       rho(ix1,jx) = rho(ix1,jx) + wix1*wjx
       ! (i+1,j+1):
       rho(ix1,jx1) = rho(ix1,jx1) + wix1*wjx1
       ! (i,j2)
       rho(ix,jx2) = rho(ix,jx2) + wix*wjx2
       ! (i+1,j2)
       rho(ix1,jx2) = rho(ix1,jx2) + wix1*wjx2
       ! (i2,j2)
       rho(ix2,jx2) = rho(ix2,jx2) + wix2*wjx2
       ! (i2,j)
       rho(ix2,jx) = rho(ix2,jx) + wix2*wjx
       ! (i2,j+1)
       rho(ix2,jx1) = rho(ix2,jx1) + wix2*wjx1
    enddo

    t_rhofas = t_rhofas + elapsedtime_Timer( t0 )

    !end subroutine deposit2d


  end subroutine deposit2dqd



  ! deposit particles onto grid using quadratic TSC form. 
  ! correction is needed to use this subroutine.
  ! there should be 1 extra-grid on each side of domain.
  subroutine deposit2dtsc(maxnpt,innp,innx,inny,hx,hy, &
       xmin,ymin,rays,rho)
    !subroutine deposit2d(maxnpt,innp,innx,inny,hx,hy, &
    !           xmin,ymin,rays,rho)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,inny,maxnpt
    double precision, intent (in), dimension (6, maxnpt) :: rays
    double precision, intent (out), dimension (innx,inny) :: rho
    !double precision, dimension (innx) :: aa,bb
    !double precision, dimension (inny) :: ay,by
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1,ix2,jx2
    double precision :: wix,wix1,wix2,wjx,wjx1,wjx2,ab,cd
    integer :: i !, j
    double precision :: t0
    double precision :: hxi,hyi

    call starttime_Timer( t0 )

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    rho=0.d0
    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=(rays(1,i)-xmin-(ix-1)*hx)*hxi
       !I am not so sure that we should use "ge" or "gt"
       if(ab.le.0.5d0) then
          ix2 = ix - 1
          wix = 0.75d0-ab*ab
          wix1 = (0.5d0+ab)**2/2 
          wix2 = (0.5d0-ab)**2/2 
       else
          ix2 = ix + 2 
          wix = (1.5d0-ab)**2/2 
          wix1 = 0.75d0 - (1-ab)*(1-ab) 
          wix2 = (ab-0.5d0)**2/2 
       endif
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=(rays(3,i)-ymin-(jx-1)*hy)*hyi
       if(cd.le.0.5d0) then
          jx2 = jx - 1
          wjx = 0.75d0-cd*cd
          wjx1 = (0.5d0+cd)**2/2 
          wjx2 = (0.5d0-cd)**2/2 
       else
          jx2 = jx + 2 
          wjx = (1.5d0-cd)**2/2 
          wjx1 = 0.75d0 - (1-cd)*(1-cd) 
          wjx2 = (cd-0.5d0)**2/2 
       endif
       ix1=ix+1
       jx1=jx+1

       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.inny)) &
       !then
       !   innp
       !   stop
       !endif

       ! (i,j):
       rho(ix,jx) = rho(ix,jx) + wix*wjx
       ! (i,j+1):
       rho(ix,jx1) = rho(ix,jx1) + wix*wjx1
       ! (i+1,j):
       rho(ix1,jx) = rho(ix1,jx) + wix1*wjx
       ! (i+1,j+1):
       rho(ix1,jx1) = rho(ix1,jx1) + wix1*wjx1
       ! (i,j2)
       rho(ix,jx2) = rho(ix,jx2) + wix*wjx2
       ! (i+1,j2)
       rho(ix1,jx2) = rho(ix1,jx2) + wix1*wjx2
       ! (i2,j2)
       rho(ix2,jx2) = rho(ix2,jx2) + wix2*wjx2
       ! (i2,j)
       rho(ix2,jx) = rho(ix2,jx) + wix2*wjx
       ! (i2,j+1)
       rho(ix2,jx1) = rho(ix2,jx1) + wix2*wjx1
    enddo

    ! correction along x - i
    !        do j = 1, inny
    !          aa(1) = 0.75
    !          bb(1) = rho(1,j)
    !          do i = 2, innx
    !            aa(i) = aa(i) - 0.125*0.125/aa(i-1)
    !            bb(i) = rho(i,j) - bb(i-1)*0.125/aa(i-1)      
    !          enddo
    !          rho(innx,j) = bb(innx)/aa(innx)
    !          do i = innx-1,1,-1
    !            rho(i,j) = (bb(i) - 0.125*rho(i+1,j))/aa(i)
    !          enddo
    !        enddo
    !        ! correction along y - j
    !        do i = 1, innx
    !          ay(1) = 0.75
    !          by(1) = rho(i,1)
    !          do j = 2, inny
    !            ay(j) = ay(j) - 0.125*0.125/ay(j-1)
    !            by(j) = rho(i,j) - by(j-1)*0.125/ay(j-1)      
    !          enddo
    !          rho(i,inny) = by(inny)/ay(inny)
    !          do j = inny-1,1,-1
    !            rho(i,j) = (by(j) - 0.125*rho(i,j+1))/ay(j)
    !          enddo
    !        enddo

    t_rhofas = t_rhofas + elapsedtime_Timer( t0 )

  end subroutine deposit2dtsc

  !end subroutine deposit2d



  ! interpolate the field from grid to particles using a linear function.
  ! here the fields are calculated from the potential using 4th order alg.
  subroutine scatter2dln(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,weight,&
       coef,nytot)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,nytot,maxnpt
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(in), dimension(innx,nytot) :: phi
    double precision, intent(in) :: weight,coef
    double precision, dimension(innx,nytot) :: egx,egy
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1
    integer :: i, j
    double precision :: ab,cd
    double precision :: hxi,hyi,eyn,exn,tmpx,tmpy

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    !Ex
    egx = 0.0
    !2nd order approximation
    !do j = 1, nytot
    !  egx(1,j) = hxi*(phi(1,j)-phi(2,j))
    !  do i = 2, innx-1
    !    egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j))
    !  enddo
    !  egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
    !enddo

    !4th order approximation of E from potential
    egx(1,1) = hxi*(phi(1,1)-phi(2,1))
    egx(1,nytot) = hxi*(phi(1,nytot)-phi(2,nytot))
    do i = 2, innx - 1
       egx(i,1) = 0.5d0*hxi*(phi(i-1,1)-phi(i+1,1))
       egx(i,nytot) = 0.5d0*hxi*(phi(i-1,nytot)-phi(i+1,nytot))
    enddo
    egx(innx,1) = hxi*(phi(innx-1,1)-phi(innx,1))
    egx(innx,nytot) = hxi*(phi(innx-1,nytot)-phi(innx,nytot))
    do j = 2, nytot-1
       egx(1,j) = hxi*(phi(1,j)-phi(2,j))
       do i = 2, innx-1
          egx(i,j) = hxi*((phi(i-1,j+1)-phi(i+1,j+1)) + &
               4*(phi(i-1,j)-phi(i+1,j))+(phi(i-1,j-1)-phi(i+1,j-1)) )/12
       enddo
       egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
    enddo

    !Ey
    egy = 0.0
    !2nd order approximation
    !do i = 1, innx
    !  egy(i,1) = hyi*(phi(i,1)-phi(i,2))
    !enddo
    !do j = 2, nytot-1
    !  do i = 1, innx
    !    egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1))
    !  enddo
    !enddo
    !do i = 1, innx
    !  egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
    !enddo

    !4th order approximation Ey from potential
    do i = 1, innx
       egy(i,1) = hyi*(phi(i,1)-phi(i,2))
    enddo
    do j = 2, nytot-1
       do i = 2, innx-1
          egy(i,j) = hyi*((phi(i+1,j-1)-phi(i+1,j+1))+ &
               4*(phi(i,j-1)-phi(i,j+1)) + (phi(i-1,j-1)-phi(i-1,j+1)) )/12
       enddo
    enddo
    do i = 1, innx
       egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
    enddo
    do j = 2, nytot-1
       egy(1,j) = hyi*(phi(1,j-1)-phi(1,j+1))/2
       egy(innx,j) = hyi*(phi(innx,j-1)-phi(innx,j+1))/2
    enddo

    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=((xmin-rays(1,i))+ix*hx)*hxi
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=((ymin-rays(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.inny)) &
       !then
       !                          innp
       !   stop
       !endif
       exn = egx(ix,jx)*ab*cd+egx(ix,jx1)*ab*(1.0d0-cd) &
            +egx(ix1,jx)*(1.0d0-ab)*cd+egx(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       eyn = egy(ix,jx)*ab*cd+egy(ix,jx1)*ab*(1.0d0-cd) &
            +egy(ix1,jx)*(1.0d0-ab)*cd+egy(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       !rays(2,i) = rays(2,i) + exn*weight*coef
       !rays(4,i) = rays(4,i) + eyn*weight*coef
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4
    enddo

  end subroutine scatter2dln



  ! interpolate field onto particles using quadratic function.
  ! you need to have extra-grid on each side of domain to use it.
  subroutine scatter2dqd(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,weight,&
       coef,nytot)
    !subroutine scatter2dA(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,weight,&
    !                     coef,nytot,myidy,npyhalf)
    implicit none
    !include 'mpif.h'
    !integer, intent(in) :: innp,innx,nytot,maxnpt,myidy,npyhalf
    integer, intent(in) :: innp,innx,nytot,maxnpt
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(in), dimension(innx,nytot) :: phi
    double precision, intent(in) :: weight,coef
    double precision, dimension(innx,nytot) :: egx,egy
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1,ix2,jx2
    integer :: i, j
    double precision :: ab,cd,wix,wix1,wix2,wjx,wjx1,wjx2
    double precision :: hxi,hyi,eyn,exn,tmpx,tmpy

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    !Ex
    egx = 0.0
    !2nd order approximation
    !do j = 1, nytot
    !  egx(1,j) = hxi*(phi(1,j)-phi(2,j))
    !  do i = 2, innx-1
    !    egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j))
    !  enddo
    !  egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
    !enddo

    !4th order approximation of E from potential
    egx(1,1) = hxi*(phi(1,1)-phi(2,1))
    egx(1,nytot) = hxi*(phi(1,nytot)-phi(2,nytot))
    do i = 2, innx - 1
       egx(i,1) = 0.5d0*hxi*(phi(i-1,1)-phi(i+1,1))
       egx(i,nytot) = 0.5d0*hxi*(phi(i-1,nytot)-phi(i+1,nytot))
    enddo
    egx(innx,1) = hxi*(phi(innx-1,1)-phi(innx,1))
    egx(innx,nytot) = hxi*(phi(innx-1,nytot)-phi(innx,nytot))
    do j = 2, nytot-1
       egx(1,j) = hxi*(phi(1,j)-phi(2,j))
       do i = 2, innx-1
          egx(i,j) = hxi*((phi(i-1,j+1)-phi(i+1,j+1)) + &
               4*(phi(i-1,j)-phi(i+1,j))+(phi(i-1,j-1)-phi(i+1,j-1)) )/12
       enddo
       egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
    enddo

    !Ey
    egy = 0.0
    !2nd order approximation
    !do i = 1, innx
    !  egy(i,1) = hyi*(phi(i,1)-phi(i,2))
    !enddo
    !do j = 2, nytot-1
    !  do i = 1, innx
    !    egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1))
    !  enddo
    !enddo
    !do i = 1, innx
    !  egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
    !enddo

    !4th order approximation Ey from potential
    do i = 1, innx
       egy(i,1) = hyi*(phi(i,1)-phi(i,2))
    enddo
    do j = 2, nytot-1
       do i = 2, innx-1
          egy(i,j) = hyi*((phi(i+1,j-1)-phi(i+1,j+1))+ &
               4*(phi(i,j-1)-phi(i,j+1)) + (phi(i-1,j-1)-phi(i-1,j+1)) )/12
       enddo
    enddo
    do i = 1, innx
       egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
    enddo
    do j = 2, nytot-1
       egy(1,j) = hyi*(phi(1,j-1)-phi(1,j+1))/2
       egy(innx,j) = hyi*(phi(innx,j-1)-phi(innx,j+1))/2
    enddo

    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=(rays(1,i)-xmin-(ix-1)*hx)*hxi
       !I am not so sure that we should use "ge" or "gt"
       if(ab.le.0.5) then
          ix2 = ix - 1
          wix = 1.0d0-ab*ab
          wix1 = 1.0d0 + 1.5d0*(ab-1.0d0) + 0.5d0*(ab - 1.0d0)**2
          wix2 = 1.0d0 - 1.5d0*(ab+1.0d0) + 0.5d0*(ab + 1.0d0)**2
       else
          ix2 = ix + 2
          wix = 1.0d0 - 1.5d0*ab + 0.5d0*ab**2
          wix1 = 1.0d0-(1.0d0-ab)*(1.0d0-ab)
          wix2 = 1.0d0 + 1.5d0*(ab-2) + 0.5d0*(ab-2)**2
       endif
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=(rays(3,i)-ymin-(jx-1)*hy)*hyi
       if(cd.le.0.5) then
          jx2 = jx - 1
          wjx = 1.0d0-cd*cd
          wjx1 = 1.0d0 + 1.5d0*(cd-1.0d0) + 0.5d0*(cd - 1.0d0)**2
          wjx2 = 1.0d0 - 1.5d0*(cd+1.0d0) + 0.5d0*(cd + 1.0d0)**2
       else
          jx2 = jx + 2
          wjx = 1.0d0 - 1.5d0*cd + 0.5d0*cd**2
          wjx1 = 1.0d0-(1.0d0-cd)*(1.0d0-cd)
          wjx2 = 1.0d0 + 1.5d0*(cd-2) + 0.5d0*(cd-2)**2
       endif
       ix1=ix+1
       jx1=jx+1

       exn = egx(ix,jx)*wix*wjx+egx(ix,jx1)*wix*wjx1+egx(ix1,jx)* &
            wix1*wjx+egx(ix1,jx1)*wix1*wjx1+egx(ix,jx2)*wix*wjx2+&
            egx(ix1,jx2)*wix1*wjx2+egx(ix2,jx2)*wix2*wjx2+&
            egx(ix2,jx)*wix2*wjx+egx(ix2,jx1)*wix2*wjx1
       eyn = egy(ix,jx)*wix*wjx+egy(ix,jx1)*wix*wjx1+egy(ix1,jx)* &
            wix1*wjx+egy(ix1,jx1)*wix1*wjx1+egy(ix,jx2)*wix*wjx2+&
            egy(ix1,jx2)*wix1*wjx2+egy(ix2,jx2)*wix2*wjx2+&
            egy(ix2,jx)*wix2*wjx+egy(ix2,jx1)*wix2*wjx1
       !rays(2,i) = rays(2,i) + exn*weight*coef
       !rays(4,i) = rays(4,i) + eyn*weight*coef
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4
    enddo

    !end subroutine scatter2dA

  end subroutine scatter2dqd



  ! interpolate field onto particles using TSC function.
  ! you need to have extra-grid on each side of domain to use it.
  !subroutine scatter2dtsc(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,&
  !                        weight,coef,nytot)
  !subroutine scatter2dA(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,&
  !     weight,coef,nytot,myidy,npyhalf)
  subroutine scatter2dtsc(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,&
                          weight,coef,nytot,myidy,npyhalf)
    implicit none
    !include 'mpif.h'
    !integer, intent(in) :: innp,innx,nytot,maxnpt
    integer, intent(in) :: innp,innx,nytot,maxnpt,myidy,npyhalf
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(in), dimension(innx,nytot) :: phi
    double precision, intent(in) :: weight,coef
    double precision, dimension(innx,nytot) :: egx,egy
    !double precision, dimension(innx) :: aa,bb
    !double precision, dimension(nytot) :: ay,by
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1,ix2,jx2
    integer :: i, j
    double precision :: ab,cd,wix,wix1,wix2,wjx,wjx1,wjx2
    double precision :: hxi,hyi,eyn,exn,tmpx,tmpy,ey0

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    !here, the beam 2 is located at 0.1 sigma away from y-axis
    !we need to substract the dipole kick which already been taken into account when the         !the close orbit is calculated.
    ey0 = 0.0d0
    if(myidy.lt.npyhalf) then 
       !          ey0 = -(1.0d0-exp(-0.1*0.1/2.0))/(0.1*1.6e-5)
       !          ey0 = -(1.0d0-exp(-0.2*0.2/2.0))/(0.2*1.6d-5)
       !          ey0 = -(1.0d0-exp(-0.3*0.3/2.0))/(0.3*1.6d-5)
       !          ey0 = -(1.0d0-exp(-0.4*0.4/2.0))/(0.4*1.6d-5)
       !          ey0 = -(1.0d0-exp(-0.6*0.6/2.0))/(0.6*1.6d-5)
       !          ey0 = -(1.0d0-exp(-0.8*0.8/2.0))/(0.8*1.6d-5)
       !          ey0 = -(1.0d0-exp(-1.0*1.0d0/2.0))/(1.0*3.17d-5)
       !          ey0 = -(1.0d0-exp(-1.2*1.2/2.0))/(1.2*1.6d-5)
       !          ey0 = -(1.0d0-exp(-1.4*1.4/2.0))/(1.4*1.6d-5)
       !          ey0 = -(1.0d0-exp(-1.6*1.6/2.0))/(1.6*1.6d-5)
       !          ey0 = -(1.0d0-exp(-1.8*1.8/2.0))/(1.8*1.6d-5)
       !          ey0 = -(1.0d0-exp(-2.0*2.0d0/2.0))/(2.0*1.6d-5)
    else 
       !          ey0 = (1.0d0-exp(-0.1*0.1/2.0))/(0.1*1.6e-5)
       !          ey0 = (1.0d0-exp(-0.2*0.2/2.0))/(0.2*1.6d-5)
       !          ey0 = (1.0d0-exp(-0.3*0.3/2.0))/(0.3*1.6d-5)
       !          ey0 = (1.0d0-exp(-0.4*0.4/2.0))/(0.4*1.6d-5)
       !          ey0 = (1.0d0-exp(-0.6*0.6/2.0))/(0.6*1.6d-5)
       !          ey0 = (1.0d0-exp(-0.8*0.8/2.0))/(0.8*1.6d-5)
       !          ey0 = (1.0d0-exp(-1.0*1.0d0/2.0))/(1.0*3.17d-5)
       !          ey0 = (1.0d0-exp(-1.2*1.2/2.0))/(1.2*1.6d-5)
       !          ey0 = (1.0d0-exp(-1.4*1.4/2.0))/(1.4*1.6d-5)
       !          ey0 = (1.0d0-exp(-1.6*1.6/2.0))/(1.6*1.6d-5)
       !          ey0 = (1.0d0-exp(-1.8*1.8/2.0))/(1.8*1.6d-5)
       !          ey0 = (1.0d0-exp(-2.0*2.0d0/2.0))/(2.0*1.6d-5)
    endif
    !nominal case
    ey0 = 0.0d0

    !Ex
    egx = 0.0d0
    !2nd order approximation
    do j = 1, nytot
       !egx(1,j) = hxi*(phi(1,j)-phi(2,j))
       egx(1,j) = hxi*(1.5d0*phi(1,j)-2*phi(2,j)+0.5d0*phi(3,j))
       do i = 2, innx-1
          egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j))
       enddo
       !egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
       egx(innx,j) = hxi*(-0.5d0*phi(innx-2,j)+2*phi(innx-1,j)-1.5d0*phi(innx,j))
    enddo

    !4th order approximation of E from potential
    !egx(1,1) = hxi*(phi(1,1)-phi(2,1))
    !egx(1,nytot) = hxi*(phi(1,nytot)-phi(2,nytot))
    !do i = 2, innx - 1
    !  egx(i,1) = 0.5d0*hxi*(phi(i-1,1)-phi(i+1,1))
    !  egx(i,nytot) = 0.5d0*hxi*(phi(i-1,nytot)-phi(i+1,nytot))
    !enddo
    !egx(innx,1) = hxi*(phi(innx-1,1)-phi(innx,1))
    !egx(innx,nytot) = hxi*(phi(innx-1,nytot)-phi(innx,nytot))
    !do j = 2, nytot-1
    !  egx(1,j) = hxi*(phi(1,j)-phi(2,j))
    !  do i = 2, innx-1
    !    egx(i,j) = hxi*((phi(i-1,j+1)-phi(i+1,j+1)) + &
    !    4*(phi(i-1,j)-phi(i+1,j))+(phi(i-1,j-1)-phi(i+1,j-1)) )/12
    !  enddo
    !  egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
    !enddo

    ! correction along x - i for Ex
    !        do j = 1, nytot
    !          aa(1) = 0.75
    !          bb(1) = egx(1,j)
    !          do i = 2, innx
    !            aa(i) = aa(i) - 0.125*0.125/aa(i-1)
    !            bb(i) = egx(i,j) - bb(i-1)*0.125/aa(i-1)      
    !          enddo
    !          egx(innx,j) = bb(innx)/aa(innx)
    !          do i = innx-1,1,-1
    !            egx(i,j) = (bb(i) - 0.125*egx(i+1,j))/aa(i)
    !          enddo
    !        enddo
    !        ! correction along y - j for Ex
    !        do i = 1, innx
    !          ay(1) = 0.75
    !          by(1) = egx(i,1)
    !          do j = 2, nytot
    !            ay(j) = ay(j) - 0.125*0.125/ay(j-1)
    !            by(j) = egx(i,j) - by(j-1)*0.125/ay(j-1)      
    !          enddo
    !          egx(i,nytot) = by(nytot)/ay(nytot)
    !          do j = nytot-1,1,-1
    !            egx(i,j) = (by(j) - 0.125*egx(i,j+1))/ay(j)
    !          enddo
    !        enddo

    !Ey
    egy = 0.0
    !2nd order approximation
    do i = 1, innx
       !egy(i,1) = hyi*(phi(i,1)-phi(i,2))
       egy(i,1) = hyi*(1.5d0*phi(i,1)-2*phi(i,2)+0.5d0*phi(i,3)) - ey0
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1)) - ey0
       enddo
    enddo
    do i = 1, innx
       !egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
       egy(i,nytot) = hyi*(-0.5d0*phi(i,nytot-2)+2*phi(i,nytot-1)-1.5d0*phi(i,nytot)) - ey0
    enddo
    !4th order approximation Ey from potential
    !do i = 1, innx
    !  egy(i,1) = hyi*(phi(i,1)-phi(i,2))
    !enddo
    !do j = 2, nytot-1
    !  do i = 2, innx-1
    !    egy(i,j) = hyi*((phi(i+1,j-1)-phi(i+1,j+1))+ &
    !    4*(phi(i,j-1)-phi(i,j+1)) + (phi(i-1,j-1)-phi(i-1,j+1)) )/12
    !  enddo
    !enddo
    !do i = 1, innx
    !  egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
    !enddo
    !do j = 2, nytot-1
    !  egy(1,j) = hyi*(phi(1,j-1)-phi(1,j+1))/2
    !  egy(innx,j) = hyi*(phi(innx,j-1)-phi(innx,j+1))/2
    !enddo

    ! correction along x - i for Ey
    !        do j = 1, nytot
    !          aa(1) = 0.75
    !          bb(1) = egy(1,j)
    !          do i = 2, innx
    !            aa(i) = aa(i) - 0.125*0.125/aa(i-1)
    !            bb(i) = egy(i,j) - bb(i-1)*0.125/aa(i-1)      
    !          enddo
    !          egy(innx,j) = bb(innx)/aa(innx)
    !          do i = innx-1,1,-1
    !            egy(i,j) = (bb(i) - 0.125*egy(i+1,j))/aa(i)
    !          enddo
    !        enddo
    !        ! correction along y - j for Ey
    !        do i = 1, innx
    !          ay(1) = 0.75
    !          by(1) = egy(i,1)
    !          do j = 2, nytot
    !            ay(j) = ay(j) - 0.125*0.125/ay(j-1)
    !            by(j) = egy(i,j) - by(j-1)*0.125/ay(j-1)      
    !          enddo
    !          egy(i,nytot) = by(nytot)/ay(nytot)
    !          do j = nytot-1,1,-1
    !            egy(i,j) = (by(j) - 0.125*egy(i,j+1))/ay(j)
    !          enddo
    !        enddo

    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=(rays(1,i)-xmin-(ix-1)*hx)*hxi
       !I am not so sure that we should use "ge" or "gt"
       if(ab.le.0.5d0) then
          ix2 = ix - 1
          wix = 0.75d0-ab*ab
          wix1 = (0.5d0+ab)**2/2 
          wix2 = (0.5d0-ab)**2/2 
       else
          ix2 = ix + 2
          wix = (1.5d0-ab)**2/2
          wix1 = 0.75d0 - (1-ab)*(1-ab)
          wix2 = (ab-0.5d0)**2/2
       endif
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=(rays(3,i)-ymin-(jx-1)*hy)*hyi
       if(cd.le.0.5d0) then
          jx2 = jx - 1
          wjx = 0.75d0-cd*cd
          wjx1 = (0.5d0+cd)**2/2 
          wjx2 = (0.5d0-cd)**2/2 
       else
          jx2 = jx + 2
          wjx = (1.5d0-cd)**2/2 
          wjx1 = 0.75d0 - (1-cd)*(1-cd) 
          wjx2 = (cd-0.5d0)**2/2 
       endif
       ix1=ix+1
       jx1=jx+1

       exn = egx(ix,jx)*wix*wjx+egx(ix,jx1)*wix*wjx1+egx(ix1,jx)* &
            wix1*wjx+egx(ix1,jx1)*wix1*wjx1+egx(ix,jx2)*wix*wjx2+&
            egx(ix1,jx2)*wix1*wjx2+egx(ix2,jx2)*wix2*wjx2+&
            egx(ix2,jx)*wix2*wjx+egx(ix2,jx1)*wix2*wjx1
       eyn = egy(ix,jx)*wix*wjx+egy(ix,jx1)*wix*wjx1+egy(ix1,jx)* &
            wix1*wjx+egy(ix1,jx1)*wix1*wjx1+egy(ix,jx2)*wix*wjx2+&
            egy(ix1,jx2)*wix1*wjx2+egy(ix2,jx2)*wix2*wjx2+&
            egy(ix2,jx)*wix2*wjx+egy(ix2,jx1)*wix2*wjx1
       !rays(2,i) = rays(2,i) + exn*weight*coef
       !rays(4,i) = rays(4,i) + eyn*weight*coef
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4
    enddo

    end subroutine scatter2dtsc

!  end subroutine scatter2dA



  ! interpolate field onto particles using TSC function.
  ! you need to have extra-grid on each side of domain to use it.
  ! this subroutine also substract the dipole kick from the oppossite beam
  ! with arbitrary separation.
  subroutine scatter2dtscoffset(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,&
       weight,coef,nytot,myidy,npyhalf,close1,close2,sigavgopp)
    implicit none
    !include 'mpif.h'
    !integer, intent(in) :: innp,innx,nytot,maxnpt
    integer, intent(in) :: innp,innx,nytot,maxnpt,myidy,npyhalf
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(in), dimension(2) :: close1,close2
    double precision, intent(in), dimension(innx,nytot) :: phi
    double precision, intent(in) :: weight,coef
    double precision, dimension(innx,nytot) :: egx,egy
    !double precision, dimension(innx) :: aa,bb
    !double precision, dimension(nytot) :: ay,by
    double precision :: xmin,ymin,hx,hy,sigavgopp
    integer :: ix,jx,ix1,jx1,ix2,jx2
    integer :: i, j
    double precision :: ab,cd,wix,wix1,wix2,wjx,wjx1,wjx2
    double precision :: hxi,hyi,eyn,exn,tmpx,tmpy,ey0,ex0,rr2,sigavg2

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    !here, the beam 2 is located at 0.1 sigma away from y-axis
    !we need to substract the dipole kick which already been taken into account when the         !the close orbit is calculated.
    rr2 = (close1(1)-close2(1))**2+(close1(2)-close2(2))**2
    sigavg2 = sigavgopp**2
    ex0 = 0.0d0
    ey0 = 0.0d0
    if(myidy.lt.npyhalf) then 
       !          ey0 = -(1.0d0-exp(-0.1*0.1/2.0))/(0.1*1.6e-5)
       !          ey0 = -(1.0d0-exp(-0.2*0.2/2.0))/(0.2*1.6d-5)
       !          ey0 = -(1.0d0-exp(-0.3*0.3/2.0))/(0.3*1.6d-5)
       !          ey0 = -(1.0d0-exp(-0.4*0.4/2.0))/(0.4*1.6d-5)
       !          ey0 = -(1.0d0-exp(-0.6*0.6/2.0))/(0.6*1.6d-5)
       !          ey0 = -(1.0d0-exp(-0.8*0.8/2.0))/(0.8*1.6d-5)
       !          ey0 = -(1.0d0-exp(-1.0*1.0d0/2.0))/(1.0*1.6d-5)
       !          ey0 = -(1.0d0-exp(-1.2*1.2/2.0))/(1.2*1.6d-5)
       !          ey0 = -(1.0d0-exp(-1.4*1.4/2.0))/(1.4*1.6d-5)
       !          ey0 = -(1.0d0-exp(-1.6*1.6/2.0))/(1.6*1.6d-5)
       !          ey0 = -(1.0d0-exp(-1.8*1.8/2.0))/(1.8*1.6d-5)
       !          ey0 = -(1.0d0-exp(-2.0*2.0d0/2.0))/(2.0*1.6d-5)
       ex0 = (close1(1)-close2(1))/rr2*(1.0d0-exp(-rr2/sigavg2/2.0d0))
       ey0 = (close1(2)-close2(2))/rr2*(1.0d0-exp(-rr2/sigavg2/2.0d0))
    else 
       !          ey0 = (1.0d0-exp(-0.1*0.1/2.0))/(0.1*1.6e-5)
       !          ey0 = (1.0d0-exp(-0.2*0.2/2.0))/(0.2*1.6d-5)
       !          ey0 = (1.0d0-exp(-0.3*0.3/2.0))/(0.3*1.6d-5)
       !          ey0 = (1.0d0-exp(-0.4*0.4/2.0))/(0.4*1.6d-5)
       !          ey0 = (1.0d0-exp(-0.6*0.6/2.0))/(0.6*1.6d-5)
       !          ey0 = (1.0d0-exp(-0.8*0.8/2.0))/(0.8*1.6d-5)
       !          ey0 = (1.0d0-exp(-1.0*1.0d0/2.0))/(1.0*1.6d-5)
       !          ey0 = (1.0d0-exp(-1.2*1.2/2.0))/(1.2*1.6d-5)
       !          ey0 = (1.0d0-exp(-1.4*1.4/2.0))/(1.4*1.6d-5)
       !          ey0 = (1.0d0-exp(-1.6*1.6/2.0))/(1.6*1.6d-5)
       !          ey0 = (1.0d0-exp(-1.8*1.8/2.0))/(1.8*1.6d-5)
       !          ey0 = (1.0d0-exp(-2.0*2.0d0/2.0))/(2.0*1.6d-5)
       ex0 = (close2(1)-close1(1))/rr2*(1.0d0-exp(-rr2/sigavg2/2.0d0))
       ey0 = (close2(2)-close1(2))/rr2*(1.0d0-exp(-rr2/sigavg2/2.0d0))
    endif
    !nominal case
    !        ey0 = 0.0d0

    !Ex
    egx = 0.0
    !2nd order approximation
    do j = 1, nytot
       !egx(1,j) = hxi*(phi(1,j)-phi(2,j))
       egx(1,j) = hxi*(1.5d0*phi(1,j)-2*phi(2,j)+0.5d0*phi(3,j)) - ex0
       do i = 2, innx-1
          egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j)) - ex0
       enddo
       !egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
       egx(innx,j) = hxi*(-0.5d0*phi(innx-2,j)+2*phi(innx-1,j)-1.5d0*phi(innx,j)) - ex0
    enddo

    !4th order approximation of E from potential
    !egx(1,1) = hxi*(phi(1,1)-phi(2,1))
    !egx(1,nytot) = hxi*(phi(1,nytot)-phi(2,nytot))
    !do i = 2, innx - 1
    !  egx(i,1) = 0.5d0*hxi*(phi(i-1,1)-phi(i+1,1))
    !  egx(i,nytot) = 0.5d0*hxi*(phi(i-1,nytot)-phi(i+1,nytot))
    !enddo
    !egx(innx,1) = hxi*(phi(innx-1,1)-phi(innx,1))
    !egx(innx,nytot) = hxi*(phi(innx-1,nytot)-phi(innx,nytot))
    !do j = 2, nytot-1
    !  egx(1,j) = hxi*(phi(1,j)-phi(2,j))
    !  do i = 2, innx-1
    !    egx(i,j) = hxi*((phi(i-1,j+1)-phi(i+1,j+1)) + &
    !    4*(phi(i-1,j)-phi(i+1,j))+(phi(i-1,j-1)-phi(i+1,j-1)) )/12
    !  enddo
    !  egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
    !enddo

    ! correction along x - i for Ex
    !        do j = 1, nytot
    !          aa(1) = 0.75
    !          bb(1) = egx(1,j)
    !          do i = 2, innx
    !            aa(i) = aa(i) - 0.125*0.125/aa(i-1)
    !            bb(i) = egx(i,j) - bb(i-1)*0.125/aa(i-1)      
    !          enddo
    !          egx(innx,j) = bb(innx)/aa(innx)
    !          do i = innx-1,1,-1
    !            egx(i,j) = (bb(i) - 0.125*egx(i+1,j))/aa(i)
    !          enddo
    !        enddo
    !        ! correction along y - j for Ex
    !        do i = 1, innx
    !          ay(1) = 0.75
    !          by(1) = egx(i,1)
    !          do j = 2, nytot
    !            ay(j) = ay(j) - 0.125*0.125/ay(j-1)
    !            by(j) = egx(i,j) - by(j-1)*0.125/ay(j-1)      
    !          enddo
    !          egx(i,nytot) = by(nytot)/ay(nytot)
    !          do j = nytot-1,1,-1
    !            egx(i,j) = (by(j) - 0.125*egx(i,j+1))/ay(j)
    !          enddo
    !        enddo

    !Ey
    egy = 0.0
    !2nd order approximation
    do i = 1, innx
       !egy(i,1) = hyi*(phi(i,1)-phi(i,2))
       egy(i,1) = hyi*(1.5d0*phi(i,1)-2*phi(i,2)+0.5d0*phi(i,3)) - ey0
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1)) - ey0
       enddo
    enddo
    do i = 1, innx
       !egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
       egy(i,nytot) = hyi*(-0.5d0*phi(i,nytot-2)+2*phi(i,nytot-1)-1.5d0*phi(i,nytot)) - ey0
    enddo
    !4th order approximation Ey from potential
    !do i = 1, innx
    !  egy(i,1) = hyi*(phi(i,1)-phi(i,2))
    !enddo
    !do j = 2, nytot-1
    !  do i = 2, innx-1
    !    egy(i,j) = hyi*((phi(i+1,j-1)-phi(i+1,j+1))+ &
    !    4*(phi(i,j-1)-phi(i,j+1)) + (phi(i-1,j-1)-phi(i-1,j+1)) )/12
    !  enddo
    !enddo
    !do i = 1, innx
    !  egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
    !enddo
    !do j = 2, nytot-1
    !  egy(1,j) = hyi*(phi(1,j-1)-phi(1,j+1))/2
    !  egy(innx,j) = hyi*(phi(innx,j-1)-phi(innx,j+1))/2
    !enddo

    ! correction along x - i for Ey
    !        do j = 1, nytot
    !          aa(1) = 0.75
    !          bb(1) = egy(1,j)
    !          do i = 2, innx
    !            aa(i) = aa(i) - 0.125*0.125/aa(i-1)
    !            bb(i) = egy(i,j) - bb(i-1)*0.125/aa(i-1)      
    !          enddo
    !          egy(innx,j) = bb(innx)/aa(innx)
    !          do i = innx-1,1,-1
    !            egy(i,j) = (bb(i) - 0.125*egy(i+1,j))/aa(i)
    !          enddo
    !        enddo
    !        ! correction along y - j for Ey
    !        do i = 1, innx
    !          ay(1) = 0.75
    !          by(1) = egy(i,1)
    !          do j = 2, nytot
    !            ay(j) = ay(j) - 0.125*0.125/ay(j-1)
    !            by(j) = egy(i,j) - by(j-1)*0.125/ay(j-1)      
    !          enddo
    !          egy(i,nytot) = by(nytot)/ay(nytot)
    !          do j = nytot-1,1,-1
    !            egy(i,j) = (by(j) - 0.125*egy(i,j+1))/ay(j)
    !          enddo
    !        enddo

    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=(rays(1,i)-xmin-(ix-1)*hx)*hxi
       !I am not so sure that we should use "ge" or "gt"
       if(ab.le.0.5) then
          ix2 = ix - 1
          wix = 0.75-ab*ab
          wix1 = (0.5+ab)**2/2 
          wix2 = (0.5-ab)**2/2 
       else
          ix2 = ix + 2
          wix = (1.5d0-ab)**2/2
          wix1 = 0.75 - (1-ab)*(1-ab)
          wix2 = (ab-0.5)**2/2
       endif
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=(rays(3,i)-ymin-(jx-1)*hy)*hyi
       if(cd.le.0.5) then
          jx2 = jx - 1
          wjx = 0.75-cd*cd
          wjx1 = (0.5+cd)**2/2 
          wjx2 = (0.5-cd)**2/2 
       else
          jx2 = jx + 2
          wjx = (1.5d0-cd)**2/2 
          wjx1 = 0.75 - (1-cd)*(1-cd) 
          wjx2 = (cd-0.5)**2/2 
       endif
       ix1=ix+1
       jx1=jx+1

       exn = egx(ix,jx)*wix*wjx+egx(ix,jx1)*wix*wjx1+egx(ix1,jx)* &
            wix1*wjx+egx(ix1,jx1)*wix1*wjx1+egx(ix,jx2)*wix*wjx2+&
            egx(ix1,jx2)*wix1*wjx2+egx(ix2,jx2)*wix2*wjx2+&
            egx(ix2,jx)*wix2*wjx+egx(ix2,jx1)*wix2*wjx1
       eyn = egy(ix,jx)*wix*wjx+egy(ix,jx1)*wix*wjx1+egy(ix1,jx)* &
            wix1*wjx+egy(ix1,jx1)*wix1*wjx1+egy(ix,jx2)*wix*wjx2+&
            egy(ix1,jx2)*wix1*wjx2+egy(ix2,jx2)*wix2*wjx2+&
            egy(ix2,jx)*wix2*wjx+egy(ix2,jx1)*wix2*wjx1
       !rays(2,i) = rays(2,i) + exn*weight*coef
       !rays(4,i) = rays(4,i) + eyn*weight*coef
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4
    enddo

  end subroutine scatter2dtscoffset

  !        end subroutine scatter2dA



  ! The following soft Gaussian model is from Miguel Furman's TSC code
  ! scatter field onto particles by soft-Gaussian approximation..
  subroutine scatter2dgauss(maxnpt,innp,rays,weight,coef,sigma,cent)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,maxnpt
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(in), dimension(2) :: cent
    double precision, intent(in), dimension(2) :: sigma
    double precision, intent(in) :: weight,coef
    !real :: sigmax,sigmay
    double precision :: sigmax,sigmay,tmpx,tmpy

    integer :: i, icerrf, ipgamma
    !real :: x,y
    double precision :: x,y
    !complex :: cf
    double complex :: cf

    !icerrf = 1 !require initialization for this option
    icerrf = 2 !Pade approximation
    ipgamma = 1

    !here, the factor 2 is from Miguel exp(-x^2/sigma^2).
    sigmax = sigma(1)*sqrt(2.0)
    sigmay = sigma(2)*sqrt(2.0)
    do i = 1, innp
       x = rays(1,i) - cent(1)
       y = rays(3,i) - cent(2)
       cf = cfMAF(x,y,sigmax,sigmay,icerrf,ipgamma)
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = rays(2,i) + weight*coef*real(cf)/2
       rays(4,i) = rays(4,i) + weight*coef*aimag(cf)/2
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4
    enddo

  end subroutine scatter2dgauss



  ! scatter field onto particles by soft-Gaussian approximation..
  subroutine scatter2dgaussfnal(maxnpt,innp,rays,weight,coef,sigma,cent)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,maxnpt
    double precision, pointer, dimension(:,:) :: rays
    double precision, intent(in), dimension(2) :: cent
    double precision, intent(in), dimension(2) :: sigma
    double precision, intent(in) :: weight,coef
    !real :: sigmax,sigmay
    double precision :: sigmax,sigmay,tmpx,tmpy

    integer :: i, icerrf, ipgamma
    !real :: x,y
    double precision :: x,y
    !complex :: cf
    double complex :: cf
    double precision :: hcur,tmp1,tmp2,tmp0 

    icerrf = 2
    !icerrf = 1  !require initialization
    ipgamma = 1
    hcur = 1.0d-3 !for Tevatron

    sigmax = sigma(1)
    sigmay = sigma(2)
    !the coef has to be double checked. 
    !there is a factor of 2 difference between Miguel's original code.
    do i = 1, innp
       x = rays(1,i) - cent(1)
       y = rays(3,i) - cent(2)
       cf = cfMAF(x,y,sigmax,sigmay,icerrf,ipgamma)
       !the following correction from px to x' is from MAD physical manual.
       tmp0 = 1.0d0+hcur*rays(1,i)-rays(6,i)
       tmp1 = rays(2,i)*tmp0
       tmp2 = rays(4,i)*tmp0
       tmp1 = tmp1 + weight*coef*real(cf)/2
       tmp2 = tmp2 + weight*coef*aimag(cf)/2
       !rays(2,i) = tmp1/tmp0
       !rays(4,i) = tmp2/tmp0
       !rays(2,i) = rays(2,i) + weight*coef*real(cf)/2
       !rays(4,i) = rays(4,i) + weight*coef*aimag(cf)/2
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = tmp1/tmp0
       rays(4,i) = tmp2/tmp0
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4

    enddo

  end subroutine scatter2dgaussfnal



  double complex function cfMAF(x,y,a,b,icerrf,ipgamma)
    !complex function cfMAF(x,y,a,b,icerrf,ipgamma)
    parameter(rootpi=1.7724538509055160272981674833d0)
    double precision :: x, y, a, b
    !	complex z,ef,z1,z2,z3,z4,cpad10,cerf,cerfe,cfroundMAF,w1,w2
    double complex z,ef,z1,z2,z3,w1,w2
    double precision :: xx,yy,aa,xa,ya,ba,root,rcsi,fx,fy,rootsq
    ! ipgamma=1 is the gaussian case
    ! ipgamma=2 is the gamma distribution with p=2
    ! switch ICERRF=
    ! 1: table interpolation to 4th order 
    ! 2: Pade approximant (Talman via Ziemann)
    ! 3: IMSL10.0 function CERFE(Z)
    if((ipgamma.lt.1).or.(ipgamma.gt.2)) then
       stop 'ERROR in cfMAF: ipgamma<1 or >2 not done; exiting...'
    end if

    cfMAF=(0.,0.)
    if(a.eq.b) then             !use round beam (exact)
       cfMAF=cfroundMAF(x,y,a,ipgamma)
       return
    end if

    if(x.ne.0.) signx=x/abs(x)  !calculate in the first quadrant only
    if(y.ne.0.) signy=y/abs(y)
    xx=abs(x)
    yy=abs(y)
    if(a.gt.b) then
       xa=xx
       ya=yy
       aa=a
       ba=b
    else if(b.gt.a) then       !reverse x <--> y
       xa=yy
       ya=xx
       aa=b
       ba=a
    end if

    z=cmplx(xa,ya)
    root=sqrt(aa**2-ba**2)
    z1=z/root
    z2=cmplx(xa*ba/aa,ya*aa/ba)/root
    rcsi=(xa/aa)**2+(ya/ba)**2

    if(icerrf.eq.1) then           !lookup table 4th order
       w1=cerf(z1)                 !Lookup table
       w2=cerf(z2)                 !Lookup table
    else if(icerrf.eq.2) then      !Pade approximant
       w1=cpad10(z1)               !Pade appr. (Talman via Ziemann)
       w2=cpad10(z2)               !Pade appr. (Talman via Ziemann)
    else if(icerrf.eq.3) then      !IMSL 10.0 function CERFE(z)
       !	 w1=CERFE(z1)		!IMSL 10.0 routine single prec.
       !	 w2=CERFE(z2)		!IMSL 10.0 routine single prec.
       !	 w1=ZERFE(z1)		!IMSL 10.0 routine double prec.
       !	 w2=ZERFE(z2)		!IMSL 10.0 routine double prec.
    end if

    ef=w1-exp(-rcsi)*w2
    ef=conjg(ef)
    ef=ef*2*(0.,1.0d0)*rootpi/root
    if(a.gt.b) then
       fx=signx*real(ef)           !restore signs to fx and fy
       fy=signy*aimag(ef)
    else if(b.gt.a) then           !reverse x <--> y
       fy=signy*real(ef)
       fx=signx*aimag(ef)
    end if
    cfMAF=cmplx(fx,fy)
    if(ipgamma.eq.1) return
    rootsq=root**2
    z=cmplx(x,y)
    z2=cmplx(x*b/a,y*a/b)
    rcsi=(x/a)**2+(y/b)**2
    z1=conjg(1+z**2/rootsq)
    z3=conjg(z2*exp(-rcsi)-z)
    cfMAF=z1*cfMAF+2*z3/rootsq
    return
  end function cfMAF



  !complex function cfroundMAF(x,y,a,ipgamma)         !round gaussian or gamma
  double complex function cfroundMAF(x,y,a,ipgamma)   !round gaussian or gamma
    double precision :: x,y,a
    !complex z,zp
    double complex z,zp
    parameter(small=1e-6)
    double precision :: rw
    ! ipgamma=1 is the gaussian case
    ! ipgamma=2 is the gamma distribution with p=2
    if((ipgamma.lt.1).or.(ipgamma.gt.2)) then
       stop 'ERROR in cfroundMAF: ipgamma<1 or >2 not done; exiting...'
    end if
    rw=(x**2+y**2)/a**2
    z=cmplx(x,y)
    if(rw.le.small) then
       cfroundMAF=(2*z/a**2)*(1-rw/2)   !Taylor exp. about the origin
    else
       zp=conjg(z)
       cfroundMAF=(2/zp)*(1-exp(-rw))
    endif
    if(ipgamma.eq.1) return
    if(rw.le.small) then
       cfroundMAF=(z/a**2)*(1-rw**2/6)  !Taylor exp. about the origin
    else
       cfroundMAF=(2/zp)*(1-(1+rw/2)*exp(-rw))
    endif
    return
  end function cfroundMAF



  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !-----------------------------------------------------------------------
  ! ***** function CERF **************************************************
  ! lifted intact from Jeff Tennyson's program TRS.
  ! This function is the implementation of Abramowits and Stegun for the
  ! complex error function for z in the first quadrant
  !-----------------------------------------------------------------------
  !complex function cerf(z)
  double complex function cerf(z)
    !complex z,dz,w,ii,z2
    double complex z,dz,w,ii,z2
    double precision :: x,y
    common/CerfTabl/w(0:4,0:39,0:30)     ! Table of error functs
    data &
         c1/0.4613135d0/,c2/0.1901635d0/,c3/0.09999216d0/,c4/1.7844927d0/, &
         c5/.002883894d0/,c6/5.5253437d0/,c7/0.5124242d0/,c8/0.2752551d0/, &
         c9/.05176536d0/,c10/2.724745d0/

    ii=(0.0,1.0d0)
    x=real(z)
    y=aimag(z)
    if(x.lt.0. .or. y.lt.0.)  &
         stop ' CERF: z not in the first quadrant! exiting...'
    ! ***** if real(x) greater than 4.5 or imag(z) greater than 3.5, use a *
    ! ***** polynomial approximation for the complex error function. *******

    if(x.ge.3.8999.or.y.ge.2.9999) go to 20

    ! ***** if not, get w(z) using third order taylor expansion to *********
    ! ***** interpolate in lookup table ************************************

    x=10.0*x
    y=10.0*y
    ix=int(x+.5)
    iy=int(y+.5)
    dz=0.1*cmplx(x-ix,y-iy)
    ! 4th-order table interpolation:
    cerf=w(0,ix,iy)+dz*(w(1,ix,iy) &
         +dz*(w(2,ix,iy)/2. &
         +dz*(w(3,ix,iy)/6. &
         +dz*(w(4,ix,iy)/24.))))
    return
20  z2=z*z
    if(x.ge.6.0.or.y.ge.6.0) go to 30
    ! Abramowitz and Stegun p. 328:
    cerf=ii*z*(c1/(z2-c2)+c3/(z2-c4)+c5/(z2-c6))
    return
    ! A*S p.328:
30  cerf=ii*z*(c7/(z2-c8)+c9/(z2-c10))
    return
  end function cerf



  ! This routine calculates the complex error function w(z) by its
  ! Pade approximation as described by Talman (AIP, 1987). The co-
  ! efficients aa and ab were calculated separately in program pade2.
  ! The polynomial in the numerator is of order 10, in the denominator
  ! of order 11.
  ! In the range, where no overflow occurs the values coincide with
  ! those of cerfw to about 5 or 6 significant figures. This is
  ! sufficient for single precision.  However, the accuracy becomes less
  ! for larger arguments. There are problems for negative imaginary
  ! values, therefore we will calculate the values for positive imaginary
  ! values and transform by eq. 7.1.12 from Abramowitz * Stegun.
  ! internal calculations are performed in double precision.
  ! single precision

  !complex function cpad10(cz)
  double complex function cpad10(cz)
    implicit double complex (c)
    real aa(10),ab(11)
    data aa/-0.183720683d1,0.173008270d1,-0.104473469d1,0.440660164d0, &
         -0.134239810d0,0.297440683d-1,-0.471679983d-2, &
         0.511415009d-3,-0.342080840d-4,0.107209367d-5/
    data ab/-0.296558599d1,0.407638815d1,-0.343111294d1,0.196673867d1, &
         -0.806936246d0,0.242087060d0,-0.531718509d-1, &
         0.839066844d-2,-0.907406326d-3,0.606324265d-4, &
         -0.190023170d-5/
    data cone/(1.0d0,0.0d0)/

    if (aimag(cz).ge.0.) then
       cu=cmplx(-aimag(cz),real(cz))
       cnum=cone+cu*(aa(1)+cu*(aa(2)+cu*(aa(3)+cu*(aa(4)+cu*(aa(5) &
            +cu*(aa(6)+cu*(aa(7)+cu*(aa(8)+cu*(aa(9)+cu*aa(10))))))))))
       cden=cone+cu*(ab(1)+cu*(ab(2)+cu*(ab(3)+cu*(ab(4)+cu*(ab(5) &
            +cu*(ab(6)+cu*(ab(7)+cu*(ab(8)+cu*(ab(9)+cu*(ab(10)+ &
            cu*ab(11)))))))))))
       cw=cnum/cden
    else
       cu=cmplx(aimag(cz),real(cz))
       cnum=cone+cu*(aa(1)+cu*(aa(2)+cu*(aa(3)+cu*(aa(4)+cu*(aa(5) &
            +cu*(aa(6)+cu*(aa(7)+cu*(aa(8)+cu*(aa(9)+cu*aa(10))))))))))
       cden=cone+cu*(ab(1)+cu*(ab(2)+cu*(ab(3)+cu*(ab(4)+cu*(ab(5) &
            +cu*(ab(6)+cu*(ab(7)+cu*(ab(8)+cu*(ab(9) &
            +cu*(ab(10)+cu*ab(11)))))))))))
       cw=cnum/cden
       !cw=2*cexp(-cz*cz)-conjg(cw)
       cw=2*cdexp(-cz*cz)-conjg(cw)
    endif
    cpad10=cw
    return
  end function cpad10



  subroutine tablew
    complex ii,z !,CERFE
    double complex w
    common/CerfTabl/w(0:4,0:39,0:30)        !table of error functs
    parameter(sqrtpi=1.7724538509055160272981674833d0)
    ii=(0.0,1.0d0)
    ! Abramowitz and Stegun pp. 325:
    do 20 i=0,39
       do 20 j=0,30
          x=0.1*i
          y=0.1*j
          !seaborg has something wrong with imsl
          !        w(0,i,j)=CERFE(cmplx(x,y))                !IMSLSFUN 10.0 single prec.
          w(0,i,j)= 1.0d0                     !IMSLSFUN 10.0 single prec.
20        continue
    ! first, second, third and fourth derivatives of w(z):
    do 30 i=0,39
       do 30 j=0,30
          z=0.1*cmplx(i,j)
          w(1,i,j)=-2*(z*w(0,i,j)-ii/sqrtpi)      !1st derivative
          w(2,i,j)=-2*(z*w(1,i,j)+1*w(0,i,j))     !2nd derivative
          w(3,i,j)=-2*(z*w(2,i,j)+2*w(1,i,j))     !3rd derivative
          w(4,i,j)=-2*(z*w(3,i,j)+3*w(2,i,j))     !4th derivative
30        continue
    return
  end subroutine tablew

  
  !---------------------------------------------------------------------------
  ! deposit particles*x' onto grid.
  subroutine deposit2dpx(maxnpt,innp,innx,inny,hx,hy, &
       xmin,ymin,rays,rho)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,inny,maxnpt
    double precision, intent (in), dimension (6, maxnpt) :: rays
    double precision, intent (out), dimension (innx,inny) :: rho
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1
    double precision :: ab,cd
    integer :: i !,j
    double precision :: t0
    double precision :: hxi,hyi

    call starttime_Timer( t0 )

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    rho=0.
    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=((xmin-rays(1,i))+ix*hx)*hxi
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=((ymin-rays(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.inny)) &
       !then
       !   innp
       !   stop
       !endif

       ! (i,j):
       rho(ix,jx) = rho(ix,jx) + ab*cd*rays(2,i)
       ! (i,j+1):
       rho(ix,jx1) = rho(ix,jx1) + ab*(1.0d0-cd)*rays(2,i)
       ! (i+1,j):
       rho(ix1,jx) = rho(ix1,jx)+(1.0d0-ab)*cd*rays(2,i)
       ! (i+1,j+1):
       rho(ix1,jx1) = rho(ix1,jx1)+(1.0d0-ab)*(1.0d0-cd)*rays(2,i)
    enddo

    t_rhofas = t_rhofas + elapsedtime_Timer( t0 )

  end subroutine deposit2dpx



  ! deposit particles*y' onto grid.
  subroutine deposit2dpy(maxnpt,innp,innx,inny,hx,hy, &
       xmin,ymin,rays,rho)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,inny,maxnpt
    double precision, intent (in), dimension (6, maxnpt) :: rays
    double precision, intent (out), dimension (innx,inny) :: rho
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1
    double precision :: ab,cd
    integer :: i !, j
    double precision :: t0
    double precision :: hxi,hyi

    call starttime_Timer( t0 )

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    rho=0.
    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=((xmin-rays(1,i))+ix*hx)*hxi
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=((ymin-rays(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.inny)) &
       !then
       !   innp
       !   stop
       !endif

       ! (i,j):
       rho(ix,jx) = rho(ix,jx) + ab*cd*rays(4,i)
       ! (i,j+1):
       rho(ix,jx1) = rho(ix,jx1) + ab*(1.0d0-cd)*rays(4,i)
       ! (i+1,j):
       rho(ix1,jx) = rho(ix1,jx)+(1.0d0-ab)*cd*rays(4,i)
       ! (i+1,j+1):
       rho(ix1,jx1) = rho(ix1,jx1)+(1.0d0-ab)*(1.0d0-cd)*rays(4,i)
    enddo

    t_rhofas = t_rhofas + elapsedtime_Timer( t0 )

  end subroutine deposit2dpy



  !wrong one
  !//interpolate field onto particles including Ez.
  subroutine scatter2dEz(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,&
       phi,phipx,phipy,weight,coef,nytot)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,nytot,maxnpt
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(in), dimension(innx,nytot) :: phi,&
         phipx,phipy
    double precision, intent(in) :: weight,coef
    double precision, dimension(innx,nytot) :: egx,egy,&
         egxpx,egypy,egz
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1
    integer :: i, j
    double precision :: ab,cd
    double precision :: hxi,hyi,eyn,exn,ezn,wtcf,tmppx,tmppy

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy
    wtcf = weight*coef

    !Ex
    egx = 0.0
    do j = 1, nytot
       egx(1,j) = hxi*(phi(1,j)-phi(2,j))
       do i = 2, innx-1
          egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j))
       enddo
       egx(innx,j) = hxi*(phi(innx-1,j)-phi(innx,j))
    enddo

    !ExPx
    egxpx = 0.0
    do j = 1, nytot
       egxpx(1,j) = hxi*(phipx(1,j)-phipx(2,j))
       do i = 2, innx-1
          egxpx(i,j) = 0.5d0*hxi*(phipx(i-1,j)-phipx(i+1,j))
       enddo
       egxpx(innx,j) = hxi*(phipx(innx-1,j)-phipx(innx,j))
    enddo

    !Ey
    egy = 0.0
    do i = 1, innx
       egy(i,1) = hyi*(phi(i,1)-phi(i,2))
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1))
       enddo
    enddo
    do i = 1, innx
       egy(i,nytot) = hyi*(phi(i,nytot-1)-phi(i,nytot))
    enddo

    !EyPy
    egypy = 0.0
    do i = 1, innx
       egypy(i,1) = hyi*(phipy(i,1)-phipy(i,2))
    enddo
    do j = 2, nytot-1
       do i = 1, innx
          egypy(i,j) = 0.5d0*hyi*(phipy(i,j-1)-phipy(i,j+1))
       enddo
    enddo
    do i = 1, innx
       egypy(i,nytot) = hyi*(phipy(i,nytot-1)-phipy(i,nytot))
    enddo

    !Ez: the calculation of Ez follows Hirata et. al.
    !//Part. Accel. 40, p.205 1993.
    do j = 1, nytot
       do i = 1, innx
          egz(i,j) = -egxpx(i,j) - egypy(i,j)
       enddo
    enddo

    do i = 1, innp
       ix=int((rays(1,i)-xmin)*hxi + 1)
       ab=((xmin-rays(1,i))+ix*hx)*hxi
       jx=int((rays(3,i)-ymin)*hyi + 1)
       cd=((ymin-rays(3,i))+jx*hy)*hyi
       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.inny)) &
       !then
       !                          innp
       !   stop
       !endif
       exn = egx(ix,jx)*ab*cd+egx(ix,jx1)*ab*(1.0d0-cd) &
            +egx(ix1,jx)*(1.0d0-ab)*cd+egx(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       eyn = egy(ix,jx)*ab*cd+egy(ix,jx1)*ab*(1.0d0-cd) &
            +egy(ix1,jx)*(1.0d0-ab)*cd+egy(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       ezn = egz(ix,jx)*ab*cd+egz(ix,jx1)*ab*(1.0d0-cd) &
            +egz(ix1,jx)*(1.0d0-ab)*cd+egz(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       rays(2,i) = rays(2,i) + exn*wtcf
       rays(4,i) = rays(4,i) + eyn*wtcf
       tmppx = rays(2,i) + 0.5d0*exn*wtcf
       tmppy = rays(4,i) + 0.5d0*eyn*wtcf
       rays(6,i) = rays(6,i) + 0.5d0*exn*wtcf*tmppx+&
            0.5d0*eyn*wtcf*tmppy + 0.5d0*ezn*wtcf
    enddo

  end subroutine scatter2dEz



  ! deposit particles onto fixed domain grid in r-theta coordinate.
  ! here, "inny" is actually the same as "nytot"
  subroutine deposit2dR(innp,innx,inny,rays,&
       rho,hx,hy,xmin,ymin)
    implicit none
    include 'mpif.h'
    integer, intent(in) :: innp,innx,inny
    double precision, intent(in) :: hx,hy,xmin,ymin
    double precision, dimension(6,innp), intent(in) :: rays
    double precision, dimension(innx,inny),intent(inout)::rho
    integer :: ix,jx,ix1,jx1
    double precision :: ab,cd
    integer :: i !,j
    double precision :: t0,ri,thi,pi
    double precision :: hxi,hyi,twopi
    double precision :: htmp,hx2

    call starttime_Timer( t0 )

    pi = 2*asin(1.0d0)
    twopi = 2*pi

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy
    rho = 0.0d0
    htmp = twopi/(innp-1)
    hx2 = hx*hx
    do i = 1, innp
       ri = sqrt(rays(1,i)*rays(1,i)+rays(3,i)*rays(3,i))
       !if(ri.gt.0.0d0) then
       !  if(rays(1,i).gt.0.0) then
       !    if(rays(3,i).gt.0.0) then
       !      thi = asin(rays(3,i)/ri)
       !    else
       !      thi = 2*pi+asin(rays(3,i)/ri)
       !    endif
       !  else
       !    thi = pi - asin(rays(3,i)/ri)
       !  endif
       !else
       !  thi = 0.0d0
       !endif
       if(rays(3,i).gt.0.0) then
          thi = acos(rays(1,i)/ri)
       else
          thi = twopi - acos(rays(1,i)/ri)
       endif
       !thi = (i-1)*htmp !test of speed
       !ix=(ri-xmin)*hxi + 1
       !ab=(ix*ix*hx*hx-(ri-xmin)*(ri-xmin))/ &
       !   (ix*ix*hx*hx-(ix-1)*(ix-1)*hx*hx)
       !jx=(thi-ymin)*hyi + 1 
       !cd=((ymin-thi)+jx*hy)*hyi
       ix=int(ri*hxi + 1)
       !ab=(ix*ix*hx*hx-ri*ri)/ &
       !   (ix*ix*hx*hx-(ix-1)*(ix-1)*hx*hx)
       ab=(ix*ix*hx2-ri*ri)/ &
            (ix*ix*hx2-(ix-1)*(ix-1)*hx2)
       jx=int(thi*hyi + 1)
       cd=(jx*hy-thi)*hyi

       ix1=ix+1
       jx1=jx+1
       !          if((ix1.le.0).or.(ix1.gt.innx)) then
       !          endif
       !          if((jx1.le.0).or.(jx1.gt.inny)) then
       !          endif
       ! (i,j):
       rho(ix,jx) = rho(ix,jx) + ab*cd
       ! (i,j+1):
       rho(ix,jx1) = rho(ix,jx1) + ab*(1.0d0-cd)
       ! (i+1,j):
       rho(ix1,jx) = rho(ix1,jx)+(1.0d0-ab)*cd
       ! (i+1,j+1):
       rho(ix1,jx1) = rho(ix1,jx1)+(1.0d0-ab)*(1.0d0-cd)
    enddo

    !this is due to the periodic BC along theta. 0 = 2pi
    do i = 1, innx
       rho(i,1) = rho(i,1)+rho(i,inny)
    enddo
    do i = 1, innx
       rho(i,inny) = rho(i,1)
    enddo

    t_rhofas = t_rhofas + elapsedtime_Timer( t0 )

  end subroutine deposit2dR



  !interpolation in r-theta coordinate
  subroutine scatter2dRorg(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,weight,&
       coef,nytot)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,nytot,maxnpt
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(in), dimension(innx,nytot) :: phi
    double precision, intent(in) :: weight,coef
    double precision, dimension(innx,nytot) :: egx,egy
    double precision :: xmin,ymin,hx,hy
    integer :: ix,jx,ix1,jx1
    integer :: i, j
    double precision :: ab,cd
    double precision :: hxi,hyi,eyn,exn,ri,thi,ern,etn,xx,pi,twopi
    double precision :: htmp,hx2,cs,ss,tmpx,tmpy

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    !Ex
    egx = 0.0
    do j = 1, nytot
       egx(1,j) = -hxi*(2*phi(2,j)-phi(1,j)*3/2-phi(3,j)/2)
       do i = 2, innx-1
          egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j))
       enddo
       egx(innx,j) = -hxi*(-2*phi(innx-1,j)+phi(innx-2,j)/2+phi(innx,j)*3/2)
    enddo

    !Ey
    egy = 0.0
    do i = 2, innx
       xx = xmin + (i-1)*hx
       egy(i,1) = 0.5d0*hyi*(phi(i,nytot-1)-phi(i,2))/xx
    enddo
    do j = 2, nytot-1
       do i = 2, innx
          xx = xmin + (i-1)*hx
          egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1))/xx
       enddo
    enddo
    do i = 2, innx
       xx = xmin + (i-1)*hx
       egy(i,nytot) = 0.5d0*hyi*(phi(i,nytot-1)-phi(i,2))/xx
    enddo

    pi = 2*asin(1.0d0)
    twopi = 2*pi
    htmp = twopi/(innp-1)
    hx2 = hx*hx
    do i = 1, innp
       ri = sqrt(rays(1,i)*rays(1,i)+rays(3,i)*rays(3,i))
       !if(ri.gt.0.0d0) then
       !  if(rays(1,i).gt.0.0) then
       !    if(rays(3,i).gt.0.0) then
       !      thi = asin(rays(3,i)/ri)
       !    else
       !      thi = 2*pi+asin(rays(3,i)/ri)
       !    endif
       !  else
       !    thi = pi - asin(rays(3,i)/ri)
       !  endif
       !else
       !  thi = 0.0d0
       !endif
       if(rays(3,i).gt.0.0) then
          thi = acos(rays(1,i)/ri)
       else
          thi = twopi - acos(rays(1,i)/ri)
       endif
       !thi = (i-1)*htmp !test of speed
       !ix=(ri-xmin)*hxi + 1
       !ab=(ix*ix*hx*hx-(ri-xmin)*(ri-xmin))/ &
       !   (ix*ix*hx*hx-(ix-1)*(ix-1)*hx*hx)
       !jx=(thi-ymin)*hyi + 1
       !cd=((ymin-thi)+jx*hy)*hyi
       ix=int(ri*hxi + 1)
       !ab=(ix*ix*hx*hx-ri*ri)/ &
       !   (ix*ix*hx*hx-(ix-1)*(ix-1)*hx*hx)
       ab=(ix*ix*hx2-ri*ri)/ &
            (ix*ix*hx2-(ix-1)*(ix-1)*hx2)
       jx=int(thi*hyi + 1)
       cd=(jx*hy-thi)*hyi

       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.nytot)) &
       !then
       !                          innp
       !   stop
       !endif
       ern = egx(ix,jx)*ab*cd+egx(ix,jx1)*ab*(1.0d0-cd) &
            +egx(ix1,jx)*(1.0d0-ab)*cd+egx(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       etn = egy(ix,jx)*ab*cd+egy(ix,jx1)*ab*(1.0d0-cd) &
            +egy(ix1,jx)*(1.0d0-ab)*cd+egy(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       cs = rays(1,i)/ri
       ss = rays(3,i)/ri
       exn = ern*cs - etn*ss
       eyn = ern*ss + etn*cs
       !exn = ern*cos(thi) - etn*sin(thi)
       !eyn = ern*sin(thi) + etn*cos(thi)
       !rays(2,i) = rays(2,i) + exn*weight*coef
       !rays(4,i) = rays(4,i) + eyn*weight*coef
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4
    enddo

  end subroutine scatter2dRorg



  !interpolation in r-theta coordinate
  subroutine scatter2dR(maxnpt,innp,innx,hx,hy,xmin,ymin,rays,phi,weight,&
       coef,nytot,myidy,npyhalf)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: innp,innx,nytot,maxnpt,myidy,npyhalf
    double precision, intent(inout), dimension(6, maxnpt) :: rays
    double precision, intent(in), dimension(innx,nytot) :: phi
    double precision, intent(in) :: weight,coef
    double precision, dimension(innx,nytot) :: egx,egy
    double precision :: xmin,ymin,hx,hy,tmpx,tmpy
    integer :: ix,jx,ix1,jx1
    integer :: i, j
    double precision :: ab,cd
    double precision :: hxi,hyi,eyn,exn,ri,thi,ern,etn,xx,pi,twopi
    double precision :: htmp,hx2,cs,ss,ex0

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy

    !Ex
    egx = 0.0
    do j = 1, nytot
       egx(1,j) = -hxi*(2*phi(2,j)-phi(1,j)*3/2-phi(3,j)/2)
       do i = 2, innx-1
          egx(i,j) = 0.5d0*hxi*(phi(i-1,j)-phi(i+1,j))
       enddo
       egx(innx,j) = -hxi*(-2*phi(innx-1,j)+phi(innx-2,j)/2+phi(innx,j)*3/2)
    enddo

    !Ey
    egy = 0.0
    do i = 2, innx
       xx = xmin + (i-1)*hx
       egy(i,1) = 0.5d0*hyi*(phi(i,nytot-1)-phi(i,2))/xx
    enddo
    do j = 2, nytot-1
       do i = 2, innx
          xx = xmin + (i-1)*hx
          egy(i,j) = 0.5d0*hyi*(phi(i,j-1)-phi(i,j+1))/xx
       enddo
    enddo
    do i = 2, innx
       xx = xmin + (i-1)*hx
       egy(i,nytot) = 0.5d0*hyi*(phi(i,nytot-1)-phi(i,2))/xx
    enddo

    !here, the beam 2 is located at +0.15 sigma away from x-axis 
    !we need to substract the dipole kick which already been taken into account when the
    !the close orbit is calculated.
    !        if(myidy.lt.npyhalf) then
    !!          ex0 = -(1.0d0-exp(-0.15*0.15/2.0))/(0.15*1.5899686e-5)
    !!          !beam 2 has 10% rms size different.
    !          ex0 = -(1.0d0-exp(-0.15*0.15/2.0))/(0.15*1.43097174e-5) 
    !        else
    !          ex0 = (1.0d0-exp(-0.15*0.15/2.0))/(0.15*1.5899686e-5)
    !        endif
    !       nominal case
    ex0 = 0.0d0
    pi = 2*asin(1.0d0)
    twopi = 2*pi
    htmp = twopi/(innp-1)
    hx2 = hx*hx
    do i = 1, innp
       ri = sqrt(rays(1,i)*rays(1,i)+rays(3,i)*rays(3,i))
       !if(ri.gt.0.0d0) then
       !  if(rays(1,i).gt.0.0) then
       !    if(rays(3,i).gt.0.0) then
       !      thi = asin(rays(3,i)/ri)
       !    else
       !      thi = 2*pi+asin(rays(3,i)/ri)
       !    endif
       !  else
       !    thi = pi - asin(rays(3,i)/ri)
       !  endif
       !else
       !  thi = 0.0d0
       !endif
       if(rays(3,i).gt.0.0) then
          thi = acos(rays(1,i)/ri)
       else
          thi = twopi - acos(rays(1,i)/ri)
       endif
       !thi = (i-1)*htmp !test of speed
       !ix=(ri-xmin)*hxi + 1
       !ab=(ix*ix*hx*hx-(ri-xmin)*(ri-xmin))/ &
       !   (ix*ix*hx*hx-(ix-1)*(ix-1)*hx*hx)
       !jx=(thi-ymin)*hyi + 1
       !cd=((ymin-thi)+jx*hy)*hyi
       ix=int(ri*hxi + 1)
       ab=(ix*ix*hx2-ri*ri)/(ix*ix*hx2-(ix-1)*(ix-1)*hx2)
       jx=int(thi*hyi + 1)
       cd=(jx*hy-thi)*hyi

       ix1=ix+1
       jx1=jx+1
       !if( (ix.lt.1).or.(ix.ge.innx).or.(jx.lt.1).or.(jx.ge.nytot)) &
       !then
       !                          innp
       !   stop
       !endif
       ern = egx(ix,jx)*ab*cd+egx(ix,jx1)*ab*(1.0d0-cd) &
            +egx(ix1,jx)*(1.0d0-ab)*cd+egx(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       etn = egy(ix,jx)*ab*cd+egy(ix,jx1)*ab*(1.0d0-cd) &
            +egy(ix1,jx)*(1.0d0-ab)*cd+egy(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd)
       cs = rays(1,i)/ri
       ss = rays(3,i)/ri
       exn = ern*cs - etn*ss - ex0
       eyn = ern*ss + etn*cs
       !exn = ern*cos(thi) - etn*sin(thi)
       !eyn = ern*sin(thi) + etn*cos(thi)
       !rays(2,i) = rays(2,i) + exn*weight*coef
       !rays(4,i) = rays(4,i) + eyn*weight*coef
       tmpx = rays(2,i)
       tmpy = rays(4,i)
       rays(2,i) = rays(2,i) + exn*weight*coef
       rays(4,i) = rays(4,i) + eyn*weight*coef
       rays(6,i) = rays(6,i) + (-tmpx**2-tmpy**2+rays(2,i)**2+rays(4,i)**2)/4
    enddo

  end subroutine scatter2dR



  subroutine interpol1(nxtot,nytot,rho,rhotran,hx,hy,xmin,ymin,rmax)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: nxtot,nytot
    double precision, intent(inout), dimension(nxtot,nytot) :: rho
    double precision, intent(inout), dimension(nytot,nxtot) ::rhotran
    double precision :: xmin,ymin,hx,hy,rmax
    integer :: ix,jx,ix1,jx1
    integer :: i, j
    double precision :: ab,cd
    double precision :: hxi,hyi,xx,yy,rr,thi,hr,hth
    double precision :: hri,hthi,twopi,hr2

    hxi = 1.0d0/hx
    hyi = 1.0d0/hy
    hr = rmax/(nxtot-1)
    hr2 = hr*hr
    hth = 4*asin(1.0d0)/(nytot-1)
    hri = 1.0d0/hr
    hthi = 1.0d0/hth
    twopi = 4*asin(1.0d0)

    rhotran = 0.0d0
    do i = 1, nxtot
       do j = 1, nytot
          !            rr=(i-1)*hr
          !            thi = (j-1)*hth
          !            xx = rr*cos(thi)
          !            yy = rr*sin(thi)
          !            ix=(xx-xmin)*hxi + 1
          !            ab=((xmin-xx)+ix*hx)*hxi
          !            jx=(yy-ymin)*hyi + 1
          !            cd=((ymin-yy)+jx*hy)*hyi
          xx = xmin + (i-1)*hx
          yy = ymin + (j-1)*hy
          rr = sqrt(xx*xx+yy*yy)
          if(yy.ge.0.0) then
             thi = acos(xx/rr)
          else
             thi = twopi - acos(xx/rr)
          endif
          ix=int(rr*hri + 1)
          ab=(ix*ix*hr2-rr*rr)/ &
               (ix*ix*hr2-(ix-1)*(ix-1)*hr2)
          jx=int(thi*hthi + 1)
          cd=(jx*hth-thi)*hthi
          ix1=ix+1
          jx1=jx+1
          !has been transposed
          !            rhotran(j,i) = rho(ix,jx)*ab*cd+rho(ix,jx1)*ab*(1.0d0-cd)+&
          !                rho(ix1,jx)*(1.0d0-ab)*cd+rho(ix1,jx1)*(1.0d0-ab)*(1.0d0-cd) 
          rhotran(jx,ix) = rhotran(jx,ix) + ab*cd*rho(i,j)
          rhotran(jx1,ix) = rhotran(jx1,ix) + ab*(1.0d0-cd)*rho(i,j)
          rhotran(jx,ix1) = rhotran(jx,ix1) + (1.0d0-ab)*cd*rho(i,j)
          rhotran(jx1,ix1) = rhotran(jx1,ix1) + (1.0d0-ab)*(1.0d0-cd)*rho(i,j)
       enddo
    enddo

  end subroutine interpol1



  subroutine interpol2(nxtot,nytot,rho,rhotran,hx,hy,xmin,ymin,rmax)
    implicit none
    !include 'mpif.h'
    integer, intent(in) :: nxtot,nytot
    double precision, intent(inout), dimension(nxtot,nytot) :: rho
    double precision, intent(inout), dimension(nytot,nxtot) ::rhotran
    double precision :: xmin,ymin,hx,hy,rmax
    integer :: ix,jx,ix1,jx1
    integer :: i, j
    double precision :: ab,cd
    double precision :: xx,yy,rr,thi,hr,hth,hr2,hri,hthi,twopi

    hr = rmax/(nxtot-1)
    hr2 = hr*hr
    hri = 1.d0/hr
    hth = 4*asin(1.0d0)/(nytot-1)
    hthi = 1.d0/hth
    twopi = 4*asin(1.0d0)

    do j = 1, nytot
       do i = 1, nxtot
          xx = xmin + (i-1)*hx
          yy = ymin + (j-1)*hy
          rr = sqrt(xx*xx+yy*yy)
          if(yy.ge.0.0) then
             thi = acos(xx/rr)
          else
             thi = twopi - acos(xx/rr)
          endif
          ix=int(rr*hri + 1)
          ab=(ix*ix*hr2-rr*rr)/ &
               (ix*ix*hr2-(ix-1)*(ix-1)*hr2)
          jx=int(thi*hthi + 1)
          cd=(jx*hth-thi)*hthi
          ix1 = ix+1
          jx1 = jx+1
          !rhtran ohas been transposed
          rho(i,j) = rhotran(jx,ix)*ab*cd+rhotran(jx1,ix)*ab*(1.0d0-cd)+&
               rhotran(jx,ix1)*(1.0d0-ab)*cd+rhotran(jx1,ix1)*(1.0d0-ab)*(1.0d0-cd)
       enddo
    enddo

  end subroutine interpol2



end module DepoScatclass
