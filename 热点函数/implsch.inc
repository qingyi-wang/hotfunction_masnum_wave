!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
!*DeckYinxq: implsch

  subroutine implsch

  implicit none

!-------------------------------------------------------------------------------

!  real :: alog10pwk, awkss, vx, vy, ww, wkpm, wkpmt, wakt, beta
  real :: awkss, vx, vy, ww, wkpm, wkpmt, wakt, beta
!  real :: deltt, deltt5, xx, wp11, wp12, wp21, wp22, wm11, wm12, wm21, wm22
  real :: xx, wp11, wp12, wp21, wp22, wm11, wm12, wm21, wm22

  real :: wp112, wp122, wp212, wp222, wm112, wm122, wm212, wm222, ffacp, ffacp1
  real :: cwks17, rij, eij, ea1, ea2, ea3, ea4, ea5, ea6, ea7, ea8, up, up1
  real :: um, um1, sap, sam, eij2, zua, ead1, ead2, fcen, ad, adp, adm, delad
  real :: deladm, cd, theta0, costh, sinth, wl, wlstar, wk0, wf0, ws0, bett
  real :: awfss, arkss, aess, eks, ekspm, sds, asiss, ssds, d0, dk, ssbo, duxdx0
  real :: duydx0, duydy0, th0, cost2, sint2, cg, cp, cgdc, cu1, cu2, cu3, sscu
  real :: wstar, deltee, eef, gadiag, eefab, amax1, sig, deladp, duxdy0

  integer :: kpmt, kakt, ks1, ks, ksp1, k, j, kh, ip, ip1, im, im1
  integer :: kp, kp1, kp2, kp3, mr, js, j11, j12, j21, j22, i


  real :: sbo

!-------------------------------------------------------------------------------

!yinxq:  alog10pwk=alog10(pwk)

  call mean2

!-------------------------------------------------------------------------------

  do 100 ic=iys,iyl
  do 100 ia=ixs,ixl

    if(nsp(ia,ic).ne.1) cycle

    awkss=awk(ia,ic)
    vx=wx(ia,ic)
    vy=wy(ia,ic)
    ww=vx**2+vy**2
    if (ww.eq.0.) ww=0.05
    wkpm=gc2/ww
    wkpmt=cksp*wkpm
    kpmt=alog10(wkpmt/wk(1))/alog10pwk+1
    wakt=cksa*awkss
    kakt=alog10(wakt/wk(1))/alog10pwk+1
    ks1=max0(kpmt,kakt)
    ks=min0(ks1,kl)
    ksp1=ks+1
    ks0(ia,ic)=ks
    kpmt0(ia,ic)=kpmt
    kakt0(ia,ic)=kakt
    do 103 k=1,ks
      fconst0(k,ia,ic)=1.0
    103      continue
    do 104 k=ksp1,kl
      fconst0(k,ia,ic)=0.0
    104      continue
  100      continue

!yinxq:  deltt=delttm*60.
!yinxq:  deltt5=deltt*0.5

  do 20000 ia=ixs,ixl
  do 20000 ic=iys,iyl

    if(nsp(ia,ic).ne.1) goto 20000

    do j=1,jl
    do k=1,klp1
      se(k,j)=0.
      dse(k,j)=0.
    enddo
    enddo

!*************************************
!      call snonlin(e)
!*************************************

!         x=0.75*d(ia,ic)*akmean(ia,ic)
  xx=0.75*d(ia,ic)*awk(ia,ic)
  if (xx.lt.0.5) xx=0.5
  enh(ia,ic)=1.+(5.5/xx)*(1.-0.833*xx)*exp(-1.25*xx)

  kh=0
  do 1001 k=1,kl
    ks=k
    wp11=wp(k,1,1)
    wp12=wp(k,1,2)
    wp21=wp(k,2,1)
    wp22=wp(k,2,2)
    wm11=wm(k,1,1)
    wm12=wm(k,1,2)
    wm21=wm(k,2,1)
    wm22=wm(k,2,2)

    wp112=wp11**2
    wp122=wp12**2
    wp212=wp21**2
    wp222=wp22**2
    wm112=wm11**2
    wm122=wm12**2
    wm212=wm21**2
    wm222=wm22**2

    ffacp=1.
    ffacp1=1.
    ip=ikp(k)
    ip1=ikp1(k)
    im=ikm(k)
    im1=ikm1(k)
    cwks17=cong*wks17(k)
    kp=ip
    kp1=ip1
    kp2=kp
    kp3=kp1
    !if(kp.lt.kl)go to 1004
    if(kp >= kl)then
      kh=kh+1
      kp2=kl+1
      if(kp.eq.kl)kp2=kl
      kp=kl
      kp1=kl
      kp3=kl+1
      ffacp=wkh(kh)
      ffacp1=wkh(kh+1)
    endif
    !1004      continue

    do 1101 mr=1,2

!*      1.2 angular loop

      1200      continue
      do 1201 j=1,jl
        js=j
        j11=jp1(mr,j)
        j12=jp2(mr,j)
        j21=jm1(mr,j)
        j22=jm2(mr,j)

!****************************************************************

!      do 1202 ia=ix1,ix2
!      do 1202 ic=iy1,iy2

        eij=e(ks,js,ia,ic)
        !if (eij.lt.1.e-20) goto 1201
        if (eij.lt.1.e-20) cycle
        ea1=e(kp ,j11,ia,ic)
        ea2=e(kp ,j12,ia,ic)
        ea3=e(kp1,j11,ia,ic)
        ea4=e(kp1,j12,ia,ic)
        ea5=e(im ,j21,ia,ic)
        ea6=e(im ,j22,ia,ic)
        ea7=e(im1,j21,ia,ic)
        ea8=e(im1,j22,ia,ic)
        up =(wp11*ea1+wp12*ea2)*ffacp
        up1=(wp21*ea3+wp22*ea4)*ffacp1
        um =wm11*ea5+wm12*ea6
        um1=wm21*ea7+wm22*ea8
        sap=up+up1
        sam=um+um1
        eij2=eij**2
        zua=2.*eij/al31
        ead1=sap/al11+sam/al21
        ead2=-2.*sap*sam/al31
        !      fcen=fcnss(k,ia,ic)*enh(ia,ic)
        fcen=fconst0(k,ia,ic)*enh(ia,ic)
        ad=cwks17*(eij2*ead1+ead2*eij)*fcen
        adp=ad/al13
        adm=ad/al23
        delad =cwks17*(eij*2.*ead1+ead2) *fcen
        deladp=cwks17*(eij2/al11-zua*sam)*fcen/al13
        deladm=cwks17*(eij2/al21-zua*sap)*fcen/al23

        !*      nonlinear transfer

        se(ks ,js )= se(ks ,js )-2.0*ad
        se(kp2,j11)= se(kp2,j11)+adp*wp11
        se(kp2,j12)= se(kp2,j12)+adp*wp12
        se(kp3,j11)= se(kp3,j11)+adp*wp21
        se(kp3,j12)= se(kp3,j12)+adp*wp22
        se(im ,j21)= se(im ,j21)+adm*wm11
        se(im ,j22)= se(im ,j22)+adm*wm12
        se(im1,j21)= se(im1,j21)+adm*wm21
        se(im1,j22)= se(im1,j22)+adm*wm22

        dse(ks ,js )=dse(ks ,js )-2.0*delad
        dse(kp2,j11)=dse(kp2,j11)+deladp*wp112
        dse(kp2,j12)=dse(kp2,j12)+deladp*wp122
        dse(kp3,j11)=dse(kp3,j11)+deladp*wp212
        dse(kp3,j12)=dse(kp3,j12)+deladp*wp222
        dse(im ,j21)=dse(im ,j21)+deladm*wm112
        dse(im ,j22)=dse(im ,j22)+deladm*wm122
        dse(im1,j21)=dse(im1,j21)+deladm*wm212
        dse(im1,j22)=dse(im1,j22)+deladm*wm222

!        1202      continue
!        1212      continue
      1201      continue
    1101      continue
  1001      continue

!      write(*,*)'endcall snon'

!*************************************
!      call sinput(e)
!*************************************

  ks=ks0(ia,ic)
  vx=wx(ia,ic)
  vy=wy(ia,ic)
  ww=vx**2+vy**2
  w(ia,ic)=sqrt(ww)
  cd=(0.80+0.065*w(ia,ic))*0.001

  do 2200 j=1,jl
    theta0=thet(j)
    costh=cos(theta0)
    sinth=sin(theta0)
    wl=vx*costh+vy*sinth
    wlstar=wl*sqrt(cd)
    do 2500 k=1,ks
      wk0=wk(k)
      wf0=wf(k,ia,ic)
      ws0=zpi*wf0
!yinxq:      bett=beta10*(1.+beta1*winc(ia,ic))
      bett=beta10
      beta=amax1(0.,bett*(wk0*28.*wlstar-ws0))

      sein(k, j) = beta * e(k, j, ia, ic) ! yinxq: 2011.05.06

      se(k,j)= se(k,j)+beta*e(k,j,ia,ic)
      dse(k,j)=dse(k,j)+beta
    2500      continue
  2200      continue

!      write(*,*)'endcall sinp'
!***************************************************
!      goto 3009
!      if (logsdiss.eq.1) goto 3009
!      call sdissip(e)--old
!***************************************************
  3009      continue
!      call sdissip(e)--new
!***************************************************

  ks=ks0(ia,ic)
  awfss=awf(ia,ic)
  asiss=asi(ia,ic)
  arkss=ark(ia,ic)
  aess = ae(ia,ic)
  eks=aess*arkss*arkss
  ekspm=eks/0.0030162
!      sds=2.36e-5*asiss*arkss**3*aess**2/alpm2
!      2.36e-5/alpm2=2.587605807
!      sds=2.60*(zpi*awfss)*arkss**3*aess**2
3010      sds=d1*asiss/arkss*sqrt(ekspm)*exp(-d2*0.64/ekspm)

  do 3505 k=1,ks
    wk0=wk(k)
    do 3202 j=1,jl

      ssds=-ads*sds*wk0

      seds(k, j) = ssds * e(k,j,ia,ic) ! yinxq: 2011.05.06

      se(k,j)= se(k,j)+ssds*e(k,j,ia,ic)
      dse(k,j)=dse(k,j)+ssds

    3202      continue
  3505      continue

3999      continue
!*************************************
!      call sbottom(e)
!*************************************

  sbo=0.038*2./g
  ks=ks0(ia,ic)
  d0=d(ia,ic)

  do 4500 k=1,ks
    wk0=wk(k)
    dk=d0*wk0
    do 4200 j=1,jl

      if (dk.ge.30.)then
      ssbo=0.
      else
      ssbo=-abo*sbo*wk0/sinh(2.*dk)
      endif

  	  sebo(k, j) = ssbo * e(k,j,ia,ic) ! yinxq: 2011.05.06

      se(k,j)= se(k,j)+ssbo*e(k,j,ia,ic)
      dse(k,j)=dse(k,j)+ssbo

    4200      continue
  4500      continue

!*************************************
!      call scurrent(e)
!*************************************

  duxdx0=uxx(ia,ic)
  duxdy0=uxy(ia,ic)
  duydx0=uyx(ia,ic)
  duydy0=uyy(ia,ic)

  ks=ks0(ia,ic)
  do 5200 j=1,jl
    th0=thet(j)
    sinth=sin(th0)
    costh=cos(th0)
    cost2=costh*costh
    sint2=sinth*sinth
    do 5500 k=1,ks
      wk0=wk(k)
      ws0=zpi*wf(k,ia,ic)
      cg=ccg(k,ia,ic)
      cp=ws0/wk0
      cgdc=cg/cp

      cu1=(cgdc*(1.+cost2)-0.5)*duxdx0
      cu2= cgdc*sinth*costh*(duxdy0+duydx0)
      cu3=(cgdc*(1.+sint2)-0.5)*duydy0
      sscu=-acu*(cu1+cu2+cu3)

      se(k,j)= se(k,j)+sscu*e(k,j,ia,ic)
      dse(k,j)=dse(k,j)+sscu

    5500      continue
  5200      continue

!*************************************
!      call source terms end
!*************************************
!      write(*,*)'endcall source terms'
  ks=ks0(ia,ic)
  vx=wx(ia,ic)
  vy=wy(ia,ic)
  ww=vx**2+vy**2
  w(ia,ic)=sqrt(ww)
  cd=(0.80+0.065*w(ia,ic))*0.001
  wstar=w(ia,ic)*sqrt(cd)

  do 6001 j=1,jl
    do 6002 k=1,ks

      deltee=wstar*grolim(k)
      gadiag=1.-deltt5*dse(k,j)
      gadiag=amax1(gadiag,1.)
      eef=deltt*se(k,j)/gadiag
      eefab=abs(eef)
      eefab=amax1(eefab,0.1e-19)
      sig=eef/eefab
      eefab=amin1(eefab,deltee)
      eef=e(k,j,ia,ic)+sig*eefab
!yinxq:      ee(k,j,ia,ic)=amax1(eef,0.)
      ee(k,j,ia,ic)=max(small, amax1(eef,0.))
    6002      continue
    do 6012 k=ks+1,kl
      i=k-ks+1
!yinxq:      ee(k,j,ia,ic)=ee(ks,j,ia,ic)*wkh(i)
      ee(k,j,ia,ic)=max(small, ee(ks,j,ia,ic)*wkh(i))
    6012      continue
  6001      continue

! yinxq: 2011.05.06
	  pein(ia,ic) = 0.0
	  peds(ia,ic) = 0.0
	  pebo(ia,ic) = 0.0

	  do k=1,kl
	  do j=1,jl
		  pein(ia,ic)=pein(ia,ic)+sein(k,j)*2.0*dwk(k)
		  pebo(ia,ic)=pebo(ia,ic)+sebo(k,j)*2.0*dwk(k)
		  peds(ia,ic)=peds(ia,ic)+seds(k,j)*2.0*dwk(k)
	  enddo !	  
	  enddo !	  
! yinxq: 2011.05.06

  20000      continue

!      write(*,*)'call water boundary'
!********************************
!      set water boundary
!********************************

  call setspec(2)

!      write(*,*)'endcall water boundary'

  return

  end subroutine implsch

!-------------------------------------------------------------------------------
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-------------------------------------------------------------------------------
