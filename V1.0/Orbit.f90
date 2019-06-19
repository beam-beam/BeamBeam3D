module Orbitclass
  use Utilityclass

contains
  !---------------------------------------------------------------------
  subroutine sweep(Ptcls,close2g,swrad,swtune,close2gold,nptlc,iturn)
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, dimension(2), intent(inout) :: close2g,close2gold
    double precision, intent(in) :: swrad,swtune
    integer, intent(in) :: nptlc,iturn
    double precision :: ang,twopi
    double precision, dimension(2) :: tmp,delta

    if(swrad.eq.0.0) return
    if(swtune.eq.0.0) return
    twopi = 4*asin(1.0d0)
    ang = twopi*swtune*iturn
    tmp(1) = swrad*cos(ang)
    tmp(2) = swrad*sin(ang)
    delta = tmp - close2gold
    close2g = close2g + delta
    do i = 1, nptlc
       Ptcls(1,i) = Ptcls(1,i) + delta(1)
       Ptcls(3,i) = Ptcls(3,i) + delta(2)
    enddo
    close2gold = tmp

  end subroutine sweep



  subroutine oscil(Ptcls,close2g,swrad,swtune,close2gold,nptlc,iturn,ang0)
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, dimension(2), intent(inout) :: close2g,close2gold
    double precision, intent(in) :: swrad,swtune,ang0
    integer, intent(in) :: nptlc,iturn
    double precision :: ang,twopi,cosang0,sinang0
    double precision, dimension(2) :: tmp,delta

    if(swrad.eq.0.0) return
    if(swtune.eq.0.0) return
    !cosang0 = close2gold(1)/sqrt(close2gold(1)**2+close2gold(2)**2) 
    !sinang0 = close2gold(2)/sqrt(close2gold(1)**2+close2gold(2)**2) 
    cosang0 = cos(ang0)
    sinang0 = sin(ang0)
    twopi = 4*asin(1.0d0)
    ang = twopi*swtune*iturn
    tmp(1) = swrad*cosang0*cos(ang)
    !for oscillation purpose
    tmp(2) = swrad*sinang0*cos(ang)
    delta = tmp - close2gold
    close2g = close2g + delta
    do i = 1, nptlc
       Ptcls(1,i) = Ptcls(1,i) + delta(1)
       Ptcls(3,i) = Ptcls(3,i) + delta(2)
    enddo
    close2gold = tmp

  end subroutine oscil



  subroutine oscil2(Ptcls,close2g,swrad,swtune,close2gold,nptlc,iturn,ang0,&
       close2g0)
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, dimension(2), intent(inout) :: close2g,close2gold
    double precision, dimension(2), intent(in) :: close2g0
    double precision, intent(in) :: swrad,swtune,ang0
    integer, intent(in) :: nptlc,iturn
    double precision :: ang,twopi,cosang0,sinang0
    double precision, dimension(2) :: tmp,delta

    if(swrad.eq.0.0) return
    if(swtune.eq.0.0) return
    !cosang0 = close2gold(1)/sqrt(close2gold(1)**2+close2gold(2)**2) 
    !sinang0 = close2gold(2)/sqrt(close2gold(1)**2+close2gold(2)**2) 
    cosang0 = cos(ang0)
    sinang0 = sin(ang0)
    twopi = 4*asin(1.0d0)
    ang = twopi*swtune*iturn
    tmp(1) = close2g0(1)+swrad*cosang0*cos(ang)
    !for oscillation purpose
    tmp(2) = close2g0(2)+swrad*sinang0*cos(ang)
    delta = tmp - close2gold
    close2g = close2g + delta
    do i = 1, nptlc
       Ptcls(1,i) = Ptcls(1,i) + delta(1)
       Ptcls(3,i) = Ptcls(3,i) + delta(2)
    enddo
    close2gold = tmp

  end subroutine oscil2



  subroutine clo_orb_squeeze(Ptcls,close2g,delta,nptlc,nmeet,iturn)
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, dimension(2), intent(inout) :: close2g
    double precision, dimension(2), intent(in) :: delta
    integer, intent(in) :: nptlc,nmeet,iturn

    if(iturn.gt.nmeet) return
    close2g = close2g + delta
    do i = 1, nptlc
       Ptcls(1,i) = Ptcls(1,i) + delta(1)
       Ptcls(3,i) = Ptcls(3,i) + delta(2)
    enddo

  end subroutine clo_orb_squeeze



  subroutine feedback(Ptcls,close2g,centrl,nptlc)
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, dimension(2), intent(in) :: centrl,close2g
    integer, intent(in) :: nptlc

    do i = 1, nptlc
       Ptcls(1,i) = Ptcls(1,i) - centrl(1) + close2g(1)
       Ptcls(3,i) = Ptcls(3,i) - centrl(2) + close2g(2)
    enddo

  end subroutine feedback



  !emulate the random centroid orbit jitter due to bunch-bunch variation,
  !parasitic collision, etc.
  subroutine orbitjitter(Ptcls,swrad,close2g,close2ginit,nptlc,iturn)
    implicit none
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, intent(in) :: swrad
    double precision, dimension(2), intent(inout) :: close2g,close2ginit
    integer, intent(in) :: nptlc,iturn
    double precision :: ang2,twopi,rd
    double precision, dimension(2) :: tmp,x,delta,tmp0
    integer :: i,iseed

    if(swrad.eq.0.0) return
    twopi = 4*asin(1.0d0)

    !      call random_number(x)
    !//random jitter of closed orbit
    iseed = -300-iturn
    call ran2(iseed,x,2)
    rd = 0.3*swrad*x(1)
    !      rd = 0.0*swrad*x(1)
    ang2 = twopi*x(2)
    tmp0(1) = rd*cos(ang2)
    tmp0(2) = rd*sin(ang2)

    !//new positions of closed orbit
    tmp(1) = close2ginit(1) + tmp0(1)
    tmp(2) = close2ginit(2) + tmp0(2)

    !shift of closed orbit during each turn
    delta = tmp - close2g
    close2g = tmp

    do i = 1, nptlc
       Ptcls(1,i) = Ptcls(1,i) + delta(1)
       Ptcls(3,i) = Ptcls(3,i) + delta(2)
    enddo

  end subroutine orbitjitter



  subroutine sweepjitter(Ptcls,close2g,swrad,swtune,close2gold,nptlc,iturn)
    implicit none
    double precision, pointer, dimension(:,:) :: Ptcls
    double precision, dimension(2), intent(inout) :: close2g,close2gold
    double precision, intent(in) :: swrad,swtune
    integer, intent(in) :: nptlc,iturn
    double precision :: ang,twopi,rd,ang2
    double precision, dimension(2) :: tmp,delta,x,tmp0
    integer :: i,iseed

    if(swrad.eq.0.0) return
    if(swtune.eq.0.0) return
    twopi = 4*asin(1.0d0)
    ang = twopi*swtune*iturn

    !//random jitter of closed orbit
    iseed = -50-iturn
    call ran2(iseed,x,2)
    !      rd = 0.3*swrad*x(1)
    !      rd = swrad*x(1)
    rd = 0.0*swrad*x(1)
    ang2 = twopi*x(2)
    tmp0(1) = rd*cos(ang2)
    tmp0(2) = rd*sin(ang2)

    !//new positions of closed orbit
    tmp(1) = swrad*cos(ang) + tmp0(1)
    tmp(2) = swrad*sin(ang) + tmp0(2)

    !shift of closed orbit during each turn
    delta = tmp - close2gold
    close2gold = tmp

    close2g = close2gold

    !      close2g = close2g + delta

    do i = 1, nptlc
       Ptcls(1,i) = Ptcls(1,i) + delta(1)
       Ptcls(3,i) = Ptcls(3,i) + delta(2)
    enddo

  end subroutine sweepjitter


end module Orbitclass
