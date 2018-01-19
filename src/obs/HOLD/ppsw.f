!http://sam.ucsd.edu/sio210/propseawater/ppsw_fortran/ppsw.f

c separate set of subroutines, adapted from WHOI CTD group
      real function sal78(cnd,t,p,m)
c function to convert conductivity ratio to salinity (m = 0)
c internal functions
c  dsal(xr,xt)  function for derivative of sal(xr,xt) with xr.
c  function rt35 :  c(35,t,0)/c(35,15,0) variation with temperature
      real function svan(s,t,p0,sigma)
      real function depth(p,lat)
      real function tf(s,p)
c   function to compute the freezing point of seawater
      real function cpsw(s,t,p0)
      real function atg(s,t,p)
      real function theta(s,t0,p0,pr)
      real function svel(s,t,p0)
      function bvfrq(s,t,p,nobs,pav,e)
c the function bvfrq
      function grady(y,p,nobs,pav,ybar)
c  function compute least squares slope 'grady' of y versus p
       function p80(dpth,xlat)
      function gravity(xlat)
      function coriol(xlat)
c
c
c n fofonoff & r millard
c
c sal78 fcn ********** mar 28 1983 *******
      real function sal78(cnd,t,p,m)
c ********************************************************************
c     the conductivity ratio (cnd) = 1.0000000 for salinity = 35 pss-78
c     temperature = 15.0 deg. celsius , and atmospheric pressure.
c********************************************************************
c
c function to convert conductivity ratio to salinity (m = 0)
c salinity to conductivity ratio (m = 1,cnd becomes input salinity)
c*****************************************************************
c   references:   also located in unesco report # 37 1981
c  practical salinity scale 1978: e.l. lewis ieee ocean eng. jan. 1980
c *******************************************************************
c units:
c       pressure        p        decibars
c       temperature     t        deg celsius (ipts-68)
c       conductivity    cnd      ratio     (m=0)
c       salinity        sal78    (pss-78)  (m=0)
c  checkvalues:
c      sal78=1.888091 :cnd= 40.0000,t=40 deg c,p= 10000 decibars:   m= 1
c      sal78=40.00000 :cnd=1.888091,t=40 deg c,p=10000 decibars:    m=0
c*******************************************************************
c sal78 ratio: returns zero for conductivity ratio:  < 0.0005
c sal78: returns zero for salinity:  < 0.02
c ****************************************************************
c internal functions
c ****************************************************************
c  practical salinity scale 1978 definition with temperature correction
c  xt=t-15.0 : xr=sqrt(rt)
      sal(xr,xt) =((((2.7081*xr-7.0261)*xr+14.0941)*xr+25.3851)*xr
     x-0.1692)*  xr+0.0080
     x  +(xt/(1.0+0.0162*xt))*(((((-0.0144*xr+
     x   0.0636)*xr-0.0375)*xr-0.0066)*xr-0.0056)*xr+0.0005)
c  dsal(xr,xt)  function for derivative of sal(xr,xt) with xr.
      dsal(xr,xt) =((((13.5405*xr-28.1044)*xr+42.2823)*xr+50.7702)*xr
     x   -0.1692)+(xt/(1.0+0.0162*xt))*((((-0.0720*xr+0.2544)*xr
     x   -0.1125)*xr-0.0132)*xr-0.0056)
c  function rt35 :  c(35,t,0)/c(35,15,0) variation with temperature
c  with temperature.
      rt35(xt) = (((1.0031e-9*xt-6.9698e-7)*xt+1.104259e-4)*xt
     x   + 2.00564e-2)*xt + 0.6766097
c polynomials of rp: c(s,t,p)/c(s,t,0) variation with pressure
c  c(xp) polynomial corresponds to a1-a3 constants: lewis 1980
      c(xp) = ((3.989e-15*xp-6.370e-10)*xp+2.070e-5)*xp
      b(xt) = (4.464e-4*xt+3.426e-2)*xt + 1.0
c  a(xt) polynomial corresponds to b3 and b4 constants: lewis 1980
      a(xt) = -3.107e-3*xt + 0.4215
c*******************************************************************
c zero salinity/conductivity trap
      sal78=0.0
      if((m.eq.0).and.(cnd.le.5e-4)) return
      if((m.eq.1).and.(cnd.le.0.02)) return
c ***********************************************
      dt = t - 15.0
c select branch for salinity (m=0) or conductivity (m=1)
      if(m.eq.1)  go to 10
c ************************************************
c convert conductivity to salinity
      r = cnd
      rt = r/(rt35(t)*(1.0 + c(p)/(b(t) + a(t)*r)))
      rt = sqrt(abs(rt))
      sal78 = sal(rt,dt)
      return
c *********  end of conductivity to salinity section ***********
c *******************************************************
c invert salinity to conductivity by the
c  newton-raphson iterative method.
c ******************************************************
c first approximation
   10 rt = sqrt(cnd/35.0)
      si = sal(rt,dt)
      n = 0
c
c  iteration loop begins here with a maximum of 10 cycles
c
   15 rt = rt + (cnd - si)/dsal(rt,dt)
      si = sal(rt,dt)
      n = n + 1
      dels = abs(si - cnd)
      if((dels.gt.1.0e-4).and.(n.lt.10))go to 15
c
c ****************************end of iteration loop ********
c
c compute conductivity ratio
      rtt = rt35(t)*rt*rt
      at = a(t)
      bt = b(t)
      cp = c(p)
      cp = rtt*(cp + bt)
      bt = bt - rtt*at
c
c solve quadratic equation for r: r=rt35*rt*(1+c/ar+b)
c
      r = sqrt(abs(bt*bt + 4.0*at*cp)) - bt
c conductivity return
      sal78 = 0.5*r/at
      return
      end
      real function svan(s,t,p0,sigma)
c  modified rcm
c ******************************************************
c specific volume anomaly (steric anomaly) based on 1980 equation
c of state for seawater and 1978 practical salinity scale.
c references
c millero, et al (1980) deep-sea res.,27a,255-264
c millero and poisson 1981,deep-sea res.,28a pp 625-629.
c both above references are also found in unesco report 38 (1981)
c units:
c       pressure        p0       decibars
c       temperature     t        deg celsius (ipts-68)
c       salinity        s        (ipss-78)
c       spec. vol. ana. svan     m**3/kg *1.0e-8
c       density ana.    sigma    kg/m**3
c ******************************************************************
c check value: svan=981.3021 e-8 m**3/kg.  for s = 40 (ipss-78) ,
c t = 40 deg c, p0= 10000 decibars.
c check value: sigma = 59.82037  kg/m**3 for s = 40 (ipss-78) ,
c t = 40 deg c, p0= 10000 decibars.
c *******************************************************
      real p,t,s,sig,sr,r1,r2,r3,r4
      real a,b,c,d,e,a1,b1,aw,bw,k,k0,kw,k35
c equiv
      equivalence (e,d,b1),(bw,b,r3),(c,a1,r2)
      equivalence (aw,a,r1),(kw,k0,k)
c ********************
c data
      data r3500,r4/1028.1063,4.8314e-4/
      data dr350/28.106331/
c   r4 is refered to as  c  in millero and poisson 1981
c convert pressure to bars and take square root salinity.
      p=p0/10.
      sr = sqrt(abs(s))
c *********************************************************
c pure water density at atmospheric pressure
c   bigg p.h.,(1967) br. j. applied physics 8 pp 521-537.
c
      r1 = ((((6.536332e-9*t-1.120083e-6)*t+1.001685e-4)*t
     x-9.095290e-3)*t+6.793952e-2)*t-28.263737
c seawater density atm press.
c  coefficients involving salinity
c  r2 = a   in notation of millero and poisson 1981
      r2 = (((5.3875e-9*t-8.2467e-7)*t+7.6438e-5)*t-4.0899e-3)*t
     x+8.24493e-1
c  r3 = b  in notation of millero and poisson 1981
      r3 = (-1.6546e-6*t+1.0227e-4)*t-5.72466e-3
c  international one-atmosphere equation of state of seawater
      sig = (r4*s + r3*sr + r2)*s + r1
c specific volume at atmospheric pressure
      v350p = 1.0/r3500
      sva = -sig*v350p/(r3500+sig)
      sigma=sig+dr350
c  scale specific vol. anamoly to normally reported units
      svan=sva*1.0e+8
      if(p.eq.0.0) return
c ******************************************************************
c ******  new high pressure equation of state for seawater ********
c ******************************************************************
c        millero, et al , 1980 dsr 27a, pp 255-264
c               constant notation follows article
c********************************************************
c compute compression terms
      e = (9.1697e-10*t+2.0816e-8)*t-9.9348e-7
      bw = (5.2787e-8*t-6.12293e-6)*t+3.47718e-5
      b = bw + e*s
c
      d = 1.91075e-4
      c = (-1.6078e-6*t-1.0981e-5)*t+2.2838e-3
      aw = ((-5.77905e-7*t+1.16092e-4)*t+1.43713e-3)*t
     x-0.1194975
      a = (d*sr + c)*s + aw
c
      b1 = (-5.3009e-4*t+1.6483e-2)*t+7.944e-2
      a1 = ((-6.1670e-5*t+1.09987e-2)*t-0.603459)*t+54.6746
      kw = (((-5.155288e-5*t+1.360477e-2)*t-2.327105)*t
     x+148.4206)*t-1930.06
      k0 = (b1*sr + a1)*s + kw
c evaluate pressure polynomial
c ***********************************************
c   k equals the secant bulk modulus of seawater
c   dk=k(s,t,p)-k(35,0,p)
c  k35=k(35,0,p)
c ***********************************************
      dk = (b*p + a)*p + k0
      k35  = (5.03217e-5*p+3.359406)*p+21582.27
      gam=p/k35
      pk = 1.0 - gam
      sva = sva*pk + (v350p+sva)*p*dk/(k35*(k35+dk))
c  scale specific vol. anamoly to normally reported units
      svan=sva*1.0e+8
      v350p = v350p*pk
c  ****************************************************
c compute density anamoly with respect to 1000.0 kg/m**3
c  1) dr350: density anamoly at 35 (ipss-78), 0 deg. c and 0 decibars
c  2) dr35p: density anamoly 35 (ipss-78), 0 deg. c ,  pres. variation
c  3) dvan : density anamoly variations involving specfic vol. anamoly
c ********************************************************************
c check value: sigma = 59.82037  kg/m**3 for s = 40 (ipss-78),
c t = 40 deg c, p0= 10000 decibars.
c *******************************************************
      dr35p=gam/v350p
      dvan=sva/(v350p*(v350p+sva))
      sigma=dr350+dr35p-dvan
      return
      end
      real function depth(p,lat)
c ********************************
c depth in meters from pressure in decibars using
c saunders and fofonoff's method.
c deep-sea res., 1976,23,109-111.
c formula refitted for 1980 equation of state
c units:
c       pressure        p        decibars
c       latitude        lat      degrees
c       depth           depth    meters
c checkvalue: depth = 9712.653 m for p=10000 decibars, latitude=30 deg
c     above for standard ocean: t=0 deg. celsuis ; s=35 (ipss-78)
c
      real lat
c
      x = sin(lat/57.29578)
c**************************
      x = x*x
c gr= gravity variation with latitude: anon (1970) bulletin geodesique
      gr = 9.780318*(1.0+(5.2788e-3+2.36e-5*x)*x) + 1.092e-6*p
      depth = (((-1.82e-15*p+2.279e-10)*p-2.2512e-5)*p+9.72659)*p
      depth=depth/gr
      return
      end
      real function tf(s,p)
c   function to compute the freezing point of seawater
c
c   reference: unesco tech. papers in the marine science no. 28. 1978
c   eighth report jpots
c   annex 6 freezing point of seawater f.j. millero pp.29-35.
c
c  units:
c         pressure      p          decibars
c         salinity      s          pss-78
c         temperature   tf         degrees celsius
c         freezing pt.
c************************************************************
c  checkvalue: tf= -2.588567 deg. c for s=40.0, p=500. decibars
      tf=(-.0575+1.710523e-3*sqrt(abs(s))-2.154996e-4*s)*s-7.53e-4*p
      return
      end
c************************************
      real function cpsw(s,t,p0)
c ****************************
c units:
c       pressure        p0       decibars
c       temperature     t        deg celsius (ipts-68)
c       salinity        s        (ipss-78)
c       specific heat   cpsw     j/(kg deg c)
c********************************************************
c ref: millero et al,1973,jgr,78,4499-4507
c       millero et al, unesco report no. 38 1981 pp. 99-188.
c pressure variation from least squares polynomial
c developed by fofonoff 1980.
c check value: cpsw = 3849.500 j/(kg deg. c) for s = 40 (ipss-78),
c t = 40 deg c, p0= 10000 decibars
c   scale pressure to bars
      p=p0/10.
c**************************
c sqrt salinity for fractional terms
      sr = sqrt(abs(s))
c specific heat cp0 for p=0 (millero et al ,unesco 1981)
      a = (-1.38385e-3*t+0.1072763)*t-7.643575
      b = (5.148e-5*t-4.07718e-3)*t+0.1770383
      c = (((2.093236e-5*t-2.654387e-3)*t+0.1412855)*t
     x    -3.720283)*t+4217.4
      cp0 = (b*sr + a)*s + c
c cp1 pressure and temperature terms for s = 0
      a = (((1.7168e-8*t+2.0357e-6)*t-3.13885e-4)*t+1.45747e-2)*t
     x   -0.49592
      b = (((2.2956e-11*t-4.0027e-9)*t+2.87533e-7)*t-1.08645e-5)*t
     x   +2.4931e-4
      c = ((6.136e-13*t-6.5637e-11)*t+2.6380e-9)*t-5.422e-8
      cp1 = ((c*p+b)*p+a)*p
c cp2 pressure and temperature terms for s > 0
      a = (((-2.9179e-10*t+2.5941e-8)*t+9.802e-7)*t-1.28315e-4)*t
     x   +4.9247e-3
      b = (3.122e-8*t-1.517e-6)*t-1.2331e-4
      a = (a+b*sr)*s
      b = ((1.8448e-11*t-2.3905e-9)*t+1.17054e-7)*t-2.9558e-6
      b = (b+9.971e-8*sr)*s
      c = (3.513e-13*t-1.7682e-11)*t+5.540e-10
      c = (c-1.4300e-12*t*sr)*s
      cp2 = ((c*p+b)*p+a)*p
c specific heat return
      cpsw = cp0 + cp1 + cp2
      return
      end
      real function atg(s,t,p)
c ****************************
c adiabatic temperature gradient deg c per decibar
c ref: bryden,h.,1973,deep-sea res.,20,401-408
c units:
c       pressure        p        decibars
c       temperature     t        deg celsius (ipts-68)
c       salinity        s        (ipss-78)
c       adiabatic       atg      deg. c/decibar
c checkvalue: atg=3.255976e-4 c/dbar for s=40 (ipss-78),
c t=40 deg c,p0=10000 decibars
      ds = s - 35.0
      atg = (((-2.1687e-16*t+1.8676e-14)*t-4.6206e-13)*p
     x+((2.7759e-12*t-1.1351e-10)*ds+((-5.4481e-14*t
     x+8.733e-12)*t-6.7795e-10)*t+1.8741e-8))*p
     x+(-4.2393e-8*t+1.8932e-6)*ds
     x+((6.6228e-10*t-6.836e-8)*t+8.5258e-6)*t+3.5803e-5
      return
      end
      real function theta(s,t0,p0,pr)
c ***********************************
c to compute local potential temperature at pr
c using bryden 1973 polynomial for adiabatic lapse rate
c and runge-kutta 4-th order integration algorithm.
c ref: bryden,h.,1973,deep-sea res.,20,401-408
c fofonoff,n.,1977,deep-sea res.,24,489-491
c units:
c       pressure        p0       decibars
c       temperature     t0       deg celsius (ipts-68)
c       salinity        s        (ipss-78)
c       reference prs   pr       decibars
c       potential tmp.  theta    deg celsius
c checkvalue: theta= 36.89073 c,s=40 (ipss-78),t0=40 deg c,
c p0=10000 decibars,pr=0 decibars
c
c      set-up intermediate temperature and pressure variables
      p=p0
      t=t0
c**************
      h = pr - p
      xk = h*atg(s,t,p)
      t = t + 0.5*xk
      q = xk
      p = p + 0.5*h
      xk = h*atg(s,t,p)
      t = t + 0.29289322*(xk-q)
      q = 0.58578644*xk + 0.121320344*q
      xk = h*atg(s,t,p)
      t = t + 1.707106781*(xk-q)
      q = 3.414213562*xk - 4.121320344*q
      p = p + 0.5*h
      xk = h*atg(s,t,p)
      theta = t + (xk-2.0*q)/6.0
      return
      end
      real function svel(s,t,p0)
c *******************************
c sound speed seawater chen and millero 1977,jasa,62,1129-1135
c units:
c       pressure        p0       decibars
c       temperature     t        deg celsius (ipts-68)
c       salinity        s        (ipss-78)
c       sound speed     svel     meters/second
c checkvalue: svel=1731.995 m/s, s=40 (ipss-78),t=40 deg c,p=10000 dbar
c
      equivalence (a0,b0,c0),(a1,b1,c1),(a2,c2),(a3,c3)
c
c   scale pressure to bars
      p=p0/10.
c**************************
      sr = sqrt(abs(s))
c s**2 term
      d = 1.727e-3 - 7.9836e-6*p
c s**3/2 term
      b1 = 7.3637e-5 +1.7945e-7*t
      b0 = -1.922e-2 -4.42e-5*t
      b = b0 + b1*p
c s**1 term
      a3 = (-3.389e-13*t+6.649e-12)*t+1.100e-10
      a2 = ((7.988e-12*t-1.6002e-10)*t+9.1041e-9)*t-3.9064e-7
      a1 = (((-2.0122e-10*t+1.0507e-8)*t-6.4885e-8)*t-1.2580e-5)*t
     x     +9.4742e-5
      a0 = (((-3.21e-8*t+2.006e-6)*t+7.164e-5)*t-1.262e-2)*t
     x     +1.389
      a = ((a3*p+a2)*p+a1)*p+a0
c s**0 term
      c3 = (-2.3643e-12*t+3.8504e-10)*t-9.7729e-9
      c2 = (((1.0405e-12*t-2.5335e-10)*t+2.5974e-8)*t-1.7107e-6)*t
     x     +3.1260e-5
      c1 = (((-6.1185e-10*t+1.3621e-7)*t-8.1788e-6)*t+6.8982e-4)*t
     x     +0.153563
      c0 = ((((3.1464e-9*t-1.47800e-6)*t+3.3420e-4)*t-5.80852e-2)*t
     x     +5.03711)*t+1402.388
      c = ((c3*p+c2)*p+c1)*p+c0
c sound speed return
      svel = c + (a+b*sr+d*s)*s
      return
      end
c
c brval ***** brunt-vaisala freq *****
c  uses 1980 equation of state
      function bvfrq(s,t,p,nobs,pav,e)
c ************************************
c units:
c       pressure        p0       decibars
c       temperature     t        deg celsius (ipts-68)
c       salinity        s        (ipss-78)
c       buoyancy freq   bvfrq    cph
c       n**2            e        radians/second
c checkvalue: bvfrq=14.57836 cph e=6.4739928e-4 rad/sec.
c            s(1)=35.0, t(1)=5.0, p(1)=1000.0
c            s(2)=35.0, t(2)=4.0, p(2)=1002.0
c  ********note result centered at pav=1001.0 dbars **********
c inputs:  s     salinity
c          t     temperature
c          p     pressure
c          nobs  number of observations to be weighted
c  note that s,t,p should be pointed, in the call to
c       bvfrq, at the beginning of the nobs to be
c       used (e.g. s(150) with nobs = 11 will do fit
c       to points 150 through 160, centered at 155.
c*****************************************************
c
c  r millard
c july 12 1982
c computes n in cycles per hour, and e=n**2 in rad/sec**2
c after formulation of breck owen's & n.p. fofonoff
c
c modified jan 1987 LDT to incorporate gaussian weighting
c maximum number of weights is 200
c minimum number of points in filter (nobs) is 3
c
c*******************************************************
      parameter ( maxwgt = 200 )
      real*4 p(1),t(1),s(1)
      dimension weight(maxwgt)
      e=0.0
      bvfrq=0.0
      if(nobs.lt.3) return
      cxx=0.0
      cx=0.0
      cxy=0.0
      cy=0.0
c compute least squares estimate of specific volume anomaly gradient
c compute weights for least squares fit
c weights are already normalized
      call bvwgt(nobs,weight)
c compute mean and slope
      do 20 k=1,nobs
  20  cx =cx+weight(k)*p(k)
      pav = cx
      do 35 k=1,nobs
      data= svan(s(k),theta(s(k),t(k),p(k),pav),pav,sig)*1.0e-8
      cxy=cxy+weight(k)*data*(p(k)-pav)
      cy =cy+weight(k)*data
      cxx=cxx+weight(k)*(p(k)-pav)**2
   35 continue
      if(cxx.eq.0.0) return
      a0=cxy/cxx
      v350p=(1./(sig+1000.))-data
      vbar=v350p+cy
      dvdp=a0
c
      if(vbar.eq.0.0) return
      e   = -.96168423e-2*dvdp/(vbar)**2
      bvfrq = 572.9578*sign(sqrt(abs(e)),e)
      return
      end
c----------------------------------------------------
c subroutine to compute weights for least squares fit in
c the function bvfrq
c uses the gaussian weights from the whoi ctd routine gau1
c
      subroutine bvwgt(nobs,w)
      dimension w(1)
      data pi/3.1415926535898/
      n = nobs/2
      alp = ((pi/n)**2)/4.5
      wsum = 0.
      do 10 i=1,n
       w(i) = (n+1) - i
       w(i) = exp( -alp * ( w(i)**2 ) )
c      w(i) = 1.0
       wsum = wsum + w(i) 
 10   continue
      w(n+1) = 1.0
      wsum = 2*wsum + w(n+1)
       do 20 i = 1,n+1
       w(i) = w(i) / wsum
 20   continue
      nn = 2 * n + 1
      do 30 i=n+2,nn
       w(i) = w(nn+1-i)
 30   continue
      return
      end
c
c   compute gradient of y versus p
c   july 15 1982
      function grady(y,p,nobs,pav,ybar)
c  function compute least squares slope 'grady' of y versus p
c  the gradient is representive of the interval centered at pav
      real*4 p(1),y(1)
      grady=0.0
      a0=0.0
      cxx=0.0
      cx=0.0
      cxy=0.0
      cy=0.0
      if(nobs.le.1) go to 30
      do 20 k=1,nobs
  20  cx =cx+p(k)
      pav=cx/nobs
      do 35 k=1,nobs
      cxy=cxy+y(k)*(p(k)-pav)
      cy =cy+y(k)
      cxx=cxx+(p(k)-pav)**2
   35 continue
      if(cxx.eq.0.0) return
      a0=cxy/cxx
      ybar=cy/nobs
   30 continue
      grady=a0
      return
      end
c
c pressure from depth from saunder's formula with eos80.
c reference: saunders,peter m., practical conversion of pressure
c            to depth., j.p.o. , april 1981.
c r millard
c march 9, 1983
c check value: p80=7500.004 dbars;for lat=30 deg., depth=7321.45 meters
       function p80(dpth,xlat)
      parameter pi=3.141592654
      plat=abs(xlat*pi/180.)
      d=sin(plat)
      c1=5.92e-3+d**2*5.25e-3
      p80=((1-c1)-sqrt(((1-c1)**2)-(8.84e-6*dpth)))/4.42e-6
      return
      end

c compute gravity & coriolius parameters
      function gravity(xlat)
c compute gravity
c rcm jan. 7 1981
c correct sin(plat) coef. from 5.032 to 5.3024
      parameter pi=3.141592654
      plat=abs(xlat*pi/180.)
      gravity=978.0318*(1.0+5.3024e-3*sin(plat)**2-5.9e-6*sin(2.*plat)
     &**2)
      return
      end
      function coriol(xlat)
      parameter pi=3.141592654
c  compute coriolius parameter from beginning latitude
      plat=abs(xlat*pi/180.)
      coriol=14.5842e-5*sin(plat)
      return
      end

