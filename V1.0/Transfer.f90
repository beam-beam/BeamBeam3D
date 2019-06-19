module Transferclass
  use Communclass

contains

  !use new coordinate x,x',y,y',z,p/p0, i.e. the coodinate used
  !by Hirata et. al.
  !from lab frame to the boost frame:
  !alpha is the angle with respect to x on x-y plane, 
  !phi is the angle with respect to s on s-x plane
  subroutine transfer1(Ptcls,alpha,phi,rad,nptlc,nptot,npy,myidy,commcol)
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, intent(in) :: alpha,phi,rad
    integer, intent(inout) :: nptlc,nptot
    integer, intent(in) :: npy,myidy,commcol
    integer :: i,i0,ilost
    double precision :: x,px,y,py,z,pz,h0,pxstart,pystart,pzstart,&
         h0start,tmp,xstart,ystart,zstart,h0spx,h0spy,h0spz
    double precision, dimension(1) :: rholc,rhotot
    integer :: ndata, op

    ilost = 0
    do i0 = 1, nptlc
       i = i0 - ilost
       !go to the coordinate used by Hirata.
       x = Ptcls(1,i0) 
       px = Ptcls(2,i0)   
       y = Ptcls(3,i0)
       py = Ptcls(4,i0)
       z = Ptcls(5,i0)
       pz = Ptcls(6,i0)
       h0 = pz+1-sqrt(abs((pz+1)**2-px**2-py**2)) 
       pxstart = px/cos(phi) - h0*cos(alpha)*tan(phi)/cos(phi)
       pystart = py/cos(phi) - h0*sin(alpha)*tan(phi)/cos(phi)
       pzstart = pz-px*cos(alpha)*tan(phi)-py*sin(alpha)*tan(phi)+&
            h0*tan(phi)**2
       h0start = pzstart+1-sqrt(abs((pzstart+1)**2-pxstart**2-pystart**2))
       tmp = sqrt(abs((pzstart+1)**2-pxstart**2-pystart**2))
       h0spx = pxstart/tmp
       h0spy = pystart/tmp
       h0spz = 1 - (pzstart+1)/tmp
       xstart=x*(1+h0spx*cos(alpha)*sin(phi))+y*h0spx*sin(alpha)*sin(phi)+&
            z*cos(alpha)*tan(phi)
       ystart=y*(1+h0spy*sin(alpha)*sin(phi))+x*h0spy*cos(alpha)*sin(phi)+&
            z*sin(alpha)*tan(phi)
       zstart = z/cos(phi)+h0spz*(x*cos(alpha)*sin(phi)+y*sin(alpha)*&
            sin(phi))
       Ptcls(1,i) = xstart
       Ptcls(2,i) = pxstart
       Ptcls(3,i) = ystart
       Ptcls(4,i) = pystart
       Ptcls(5,i) = zstart
       Ptcls(6,i) = pzstart
       radtest = sqrt(Ptcls(1,i)*Ptcls(1,i)+Ptcls(3,i)*Ptcls(3,i))
       if(radtest.ge.rad) then
          ilost = ilost + 1
       endif
    enddo
    nptlc = nptlc - ilost
    if(ilost.gt.0) then
       op = 1
       ndata = 1
       rholc(1) = nptlc
       call allreduce(rholc,rhotot,ndata,npy,myidy,commcol,op)
       nptot = int(rhotot(1) + 0.1)
    endif

  end subroutine transfer1



  !use new coordinate x,x',y,y',z,p/p0, i.e. the coodinate used
  !by Hirata et. al.
  !from lab frame to the boost frame:
  !alpha is the angle with respect to x on x-y plane, 
  !phi is the angle with respect to s on s-x plane
  subroutine transfer1old2(Ptcls,alpha,phi,nptlc)
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, intent(in) :: alpha,phi
    integer, intent(in) :: nptlc
    integer :: i

    double precision :: x,px,y,py,z,pz,h0,pxstart,pystart,pzstart,&
         h0start,tmp,xstart,ystart,zstart,h0spx,h0spy,h0spz
    do i = 1, nptlc
       !go to the coordinate used by Hirata.
       x = Ptcls(1,i) 
       px = Ptcls(2,i)   
       y = Ptcls(3,i)
       py = Ptcls(4,i)
       z = Ptcls(5,i)
       pz = Ptcls(6,i)
       h0 = pz+1-sqrt(abs((pz+1)**2-px**2-py**2)) 
       pxstart = px/cos(phi) - h0*cos(alpha)*tan(phi)/cos(phi)
       pystart = py/cos(phi) - h0*sin(alpha)*tan(phi)/cos(phi)
       pzstart = pz-px*cos(alpha)*tan(phi)-py*sin(alpha)*tan(phi)+&
            h0*tan(phi)**2
       h0start = pzstart+1-sqrt(abs((pzstart+1)**2-pxstart**2-pystart**2))
       tmp = sqrt(abs((pzstart+1)**2-pxstart**2-pystart**2))
       h0spx = pxstart/tmp
       h0spy = pystart/tmp
       h0spz = 1 - (pzstart+1)/tmp
       xstart=x*(1+h0spx*cos(alpha)*sin(phi))+y*h0spx*sin(alpha)*sin(phi)+&
            z*cos(alpha)*tan(phi)
       ystart=y*(1+h0spy*sin(alpha)*sin(phi))+x*h0spy*cos(alpha)*sin(phi)+&
            z*sin(alpha)*tan(phi)
       zstart = z/cos(phi)+h0spz*(x*cos(alpha)*sin(phi)+y*sin(alpha)*&
            sin(phi))
       Ptcls(1,i) = xstart
       Ptcls(2,i) = pxstart
       Ptcls(3,i) = ystart
       Ptcls(4,i) = pystart
       Ptcls(5,i) = zstart
       Ptcls(6,i) = pzstart
    enddo

  end subroutine transfer1old2



  !from boost frame to the lab frame:
  !alpha is the angle with respect to x on x-y plane, 
  !phi is the angle with respect to s on s-x plane
  subroutine transfer2(Ptcls,alpha,phi,nptlc)
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, intent(in) :: alpha,phi
    integer, intent(in) :: nptlc
    integer :: i
    double precision :: px,py,pz,pxstart,pystart,pzstart,&
         h0start
    double precision :: a11,a12,a13,a21,a22,a23,a31,a32,a33,&
         b11,b12,b13,b21,b22,b23,b31,b32,b33,tmp,h0spx,h0spy,h0spz,det,&
         tmpx,tmpy,tmpz

    do i = 1, nptlc
       !go to the coordinate used by Hirata.
       pxstart = Ptcls(2,i)  
       pystart = Ptcls(4,i)
       pzstart = Ptcls(6,i)
       h0start = pzstart+1-sqrt(abs((pzstart+1)**2-pxstart**2-pystart**2)) 
       px = (pxstart + sin(phi)*cos(alpha)*h0start)*cos(phi)
       py = (pystart + sin(alpha)*sin(phi)*h0start)*cos(phi)
       pz = (pxstart*cos(alpha)*tan(phi)+pystart*sin(alpha)*tan(phi)+&
            pzstart/cos(phi))*cos(phi)
       tmp = sqrt(abs((pzstart+1)**2-pxstart**2-pystart**2))
       h0spx = pxstart/tmp
       h0spy = pystart/tmp
       h0spz = 1 - (pzstart+1)/tmp
       a11 = 1+h0spx*cos(alpha)*sin(phi)
       a12 = h0spx*sin(alpha)*sin(phi)
       a13 = cos(alpha)*tan(phi)
       a21 = h0spy*cos(alpha)*sin(phi)
       a22 = 1+h0spy*sin(alpha)*sin(phi)
       a23 = sin(alpha)*tan(phi)
       a31 = h0spz*cos(alpha)*sin(phi)
       a32 = h0spz*sin(alpha)*sin(phi)
       a33 = 1/cos(phi)
       det = a11*(a22*a33-a32*a23)-a12*(a21*a33-a31*a23)+&
            a13*(a21*a32-a22*a31)
       b11 = a22*a33-a23*a32
       b12 = -a12*a33+a13*a32
       b13 = a12*a23-a22*a13
       b21 = -a21*a33+a31*a23
       b22 = a11*a33-a13*a31
       b23 = -a11*a23+a21*a13
       b31 = a21*a32-a22*a31
       b32 = -a11*a32+a12*a31
       b33 = a11*a22-a12*a21
       tmpx = (b11*Ptcls(1,i)+b12*Ptcls(3,i)+b13*Ptcls(5,i))/det
       tmpy = (b21*Ptcls(1,i)+b22*Ptcls(3,i)+b23*Ptcls(5,i))/det
       tmpz = (b31*Ptcls(1,i)+b32*Ptcls(3,i)+b33*Ptcls(5,i))/det
       Ptcls(1,i) = tmpx
       Ptcls(2,i) = px
       Ptcls(3,i) = tmpy
       Ptcls(4,i) = py
       Ptcls(5,i) = tmpz
       Ptcls(6,i) = pz
    enddo

  end subroutine transfer2



  !use old coordinate, ie. x,x',y,y',z/sigmaz,pz/sigmapz
  subroutine transfer1old(Ptcls,sigmaz,sigmapz,pz0,alpha,phi,nptlc)
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, intent(in) :: sigmaz,sigmapz,pz0,alpha,phi
    integer, intent(in) :: nptlc
    integer :: i
    double precision :: x,px,y,py,z,pz,h0,pxstart,pystart,pzstart,&
         h0start,tmp,xstart,ystart,zstart,h0spx,h0spy,h0spz

    do i = 1, nptlc
       !go to the coordinate used by Hirata.
       x = Ptcls(1,i) 
       px = Ptcls(2,i)   
       y = Ptcls(3,i)
       py = Ptcls(4,i)
       z = Ptcls(5,i)*sigmaz
       pz = Ptcls(6,i)*sigmapz/pz0
       h0 = pz+1-sqrt(abs((pz+1)**2-px**2-py**2)) 
       pxstart = px/cos(phi) - h0*cos(alpha)*tan(phi)/cos(phi)
       pystart = py/cos(phi) - h0*sin(alpha)*tan(phi)/cos(phi)
       pzstart = pz-px*cos(alpha)*tan(phi)-py*sin(alpha)*tan(phi)+&
            h0*tan(phi)**2
       h0start = pzstart+1-sqrt(abs((pzstart+1)**2-pxstart**2-pystart**2))
       tmp = sqrt(abs((pzstart+1)**2-pxstart**2-pystart**2))
       h0spx = pxstart/tmp
       h0spy = pystart/tmp
       h0spz = 1 - (pzstart+1)/tmp
       xstart=x*(1+h0spx*cos(alpha)*sin(phi))+y*h0spx*sin(alpha)*sin(phi)+&
            z*cos(alpha)*tan(phi)
       ystart=y*(1+h0spy*sin(alpha)*sin(phi))+x*h0spy*cos(alpha)*sin(phi)+&
            z*sin(alpha)*tan(phi)
       zstart = z/cos(phi)+h0spz*(x*cos(alpha)*sin(phi)+y*sin(alpha)*&
            sin(phi))
       Ptcls(1,i) = xstart
       Ptcls(2,i) = pxstart
       Ptcls(3,i) = ystart
       Ptcls(4,i) = pystart
       Ptcls(5,i) = zstart/sigmaz
       Ptcls(6,i) = pzstart*pz0/sigmapz
    enddo

  end subroutine transfer1old



  subroutine transfer2old(Ptcls,sigmaz,sigmapz,pz0,alpha,phi,nptlc)
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, intent(in) :: sigmaz,sigmapz,pz0,alpha,phi
    integer, intent(in) :: nptlc
    integer :: i
    double precision :: px,py,pz,pxstart,pystart,pzstart,&
         h0start
    double precision :: a11,a12,a13,a21,a22,a23,a31,a32,a33,&
         b11,b12,b13,b21,b22,b23,b31,b32,b33,tmp,h0spx,h0spy,h0spz,det,&
         tmpx,tmpy,tmpz

    do i = 1, nptlc
       !go to the coordinate used by Hirata.
       pxstart = Ptcls(2,i)  
       pystart = Ptcls(4,i)
       pzstart = Ptcls(6,i)*sigmapz/pz0
       h0start = pzstart+1-sqrt(abs((pzstart+1)**2-pxstart**2-pystart**2)) 
       px = (pxstart + sin(phi)*cos(alpha)*h0start)*cos(phi)
       py = (pystart + sin(alpha)*sin(phi)*h0start)*cos(phi)
       pz = (pxstart*cos(alpha)*tan(phi)+pystart*sin(alpha)*tan(phi)+&
            pzstart/cos(phi))*cos(phi)
       tmp = sqrt(abs((pzstart+1)**2-pxstart**2-pystart**2))
       h0spx = pxstart/tmp
       h0spy = pystart/tmp
       h0spz = 1 - (pzstart+1)/tmp
       a11 = 1+h0spx*cos(alpha)*sin(phi)
       a12 = h0spx*sin(alpha)*sin(phi)
       a13 = cos(alpha)*tan(phi)
       a21 = h0spy*cos(alpha)*sin(phi)
       a22 = 1+h0spy*sin(alpha)*sin(phi)
       a23 = sin(alpha)*tan(phi)
       a31 = h0spz*cos(alpha)*sin(phi)
       a32 = h0spz*sin(alpha)*sin(phi)
       a33 = 1/cos(phi)
       det = a11*(a22*a33-a32*a23)-a12*(a21*a33-a31*a23)+&
            a13*(a21*a32-a22*a31)
       b11 = a22*a33-a23*a32
       b12 = -a12*a33+a13*a32
       b13 = a12*a23-a22*a13
       b21 = -a21*a33+a31*a23
       b22 = a11*a33-a13*a31
       b23 = -a11*a23+a21*a13
       b31 = a21*a32-a22*a31
       b32 = -a11*a32+a12*a31
       b33 = a11*a22-a12*a21
       tmpx = (b11*Ptcls(1,i)+b12*Ptcls(3,i)+b13*Ptcls(5,i)*sigmaz)/det
       tmpy = (b21*Ptcls(1,i)+b22*Ptcls(3,i)+b23*Ptcls(5,i)*sigmaz)/det
       tmpz = (b31*Ptcls(1,i)+b32*Ptcls(3,i)+b33*Ptcls(5,i)*sigmaz)/det
       Ptcls(1,i) = tmpx
       Ptcls(2,i) = px
       Ptcls(3,i) = tmpy
       Ptcls(4,i) = py
       Ptcls(5,i) = tmpz/sigmaz
       Ptcls(6,i) = pz*pz0/sigmapz
    enddo

  end subroutine transfer2old


end module Transferclass
