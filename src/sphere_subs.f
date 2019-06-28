! SPHERE_SUBS is set of FORTRAN subroutines to compute distances
! and angles on a spherical Earth.  All angles are in degrees.
! Latitude used on standard maps is geographic latitude; this may be
! converted to the geocentric latitude used in these routines
! by using the SPH_GEOCENTRIC subroutine.  Longitude input to
! these routines may be from either -180 to 180 or from 0 to 360.
! Longitude returned from these routines will be from 0 to 360.
!
!                                                Peter Shearer
!                                                pshearer@ucsd.edu

! SPH_MID finds midpoint between two surface points on sphere
! and azimuth at midpoint to second point
!
! Requires:  SPH_AZI, SPH_LOC
!
! Inputs:  flat1  =  latitude of first point (degrees) 
!          flon1  =  longitude of first point (degrees)
!          flat2  =  latitude of second point (degrees)
!          flon2  =  longitude of second point (degrees)
! Returns: del    =  angular separation between points (degrees)
!          flat3  =  latitude of midpoint (degrees)
!          flon3  =  longitude of midpoint (degrees)
!          azi    =  azimuth at midpoint to first point
!
      subroutine SPH_MID(flat1,flon1,flat2,flon2,del,flat3,flon3,azi)
      call SPH_AZI(flat1,flon1,flat2,flon2,del,azi0)
      del0=del/2.
      call SPH_LOC(flat1,flon1,del0,azi0,flat3,flon3)
      call SPH_AZI(flat3,flon3,flat2,flon2,del2,azi)
      return
      end


! SPH_LOC finds location of second point on sphere, given range 
! and azimuth at first point.
!
! Inputs:  flat1  =  latitude of first point (degrees) 
!          flon2  =  longitude of first point (degrees)
!          del    =  angular separation between points (degrees)
!          azi    =  azimuth at 1st point to 2nd point, from N (deg.)
! Returns: flat2  =  latitude of second point (degrees)
!          flon2  =  longitude of second point (degrees)
!
      subroutine SPH_LOC(flat1,flon1,del,azi,flat2,flon2)
      if (del.eq.0.) then
         flat2=flat1
         flon2=flon1
         return
      end if
      pi=3.141592654
      raddeg=pi/180.
      delr=del*raddeg
      azr=azi*raddeg
      theta1=(90.-flat1)*raddeg
      phi1=flon1*raddeg      
      ctheta2=sin(delr)*sin(theta1)*cos(azr)+cos(theta1)*cos(delr)
      theta2=acos(ctheta2)
      if (theta1.eq.0.) then
         phi2=azr
      else if (theta2.eq.0.) then
         phi2=0.0
      else
         sphi2=sin(delr)*sin(azr)/sin(theta2)
         cphi2=(cos(delr)-cos(theta1)*ctheta2)/(sin(theta1)*sin(theta2))
         phi2=phi1+atan2(sphi2,cphi2)
      end if
      flat2=90.-theta2/raddeg
      flon2=phi2/raddeg
      if (flon2.gt.360.) flon2=flon2-360.
      if (flon2.lt.0.) flon2=flon2+360.
      return
      end


! SPH_DIST computes angular separation of two points on sphere 
!
! Inputs:  flat1  =  latitude of first point (degrees)
!          flon1  =  longitude of first point (degrees)
!          flat2  =  latitude of second point (degrees)
!          flon2  =  longitude of second point (degrees)
! Returns: del    =  angular separation between points (degrees)
!
! Note:  This routine is inaccurate for del less than about 0.5 degrees. 
!        For greater accuracy, use SPH_DISTDP or perform a separate
!        calculation for close ranges using Cartesian geometry.
!
      subroutine SPH_DIST(flat1,flon1,flat2,flon2,del)
      pi=3.141592654
      raddeg=pi/180.
      theta1=(90.-flat1)*raddeg
      theta2=(90.-flat2)*raddeg
      phi1=flon1*raddeg
      phi2=flon2*raddeg
      ang=acos( sin(theta1)*sin(theta2)*cos(phi2-phi1)+
     &          cos(theta1)*cos(theta2) )
      del=ang/raddeg
      return
      end


! SPH_AZI computes distance and azimuth between two points on sphere
!
! Inputs:  flat1  =  latitude of first point (degrees) 
!          flon2  =  longitude of first point (degrees)
!          flat2  =  latitude of second point (degrees)
!          flon2  =  longitude of second point (degrees)
! Returns: del    =  angular separation between points (degrees)
!          azi    =  azimuth at 1st point to 2nd point, from N (deg.)
!
! Note:  This routine is inaccurate for del less than about 0.5 degrees. 
!        For greater accuracy, use SPH_AZIDP or perform a separate
!        calculation for close ranges using Cartesian geometry.
!
      subroutine SPH_AZI(flat1,flon1,flat2,flon2,del,azi)
      if ((flat1.eq.flat2.and.flon1.eq.flon2).or.
     &    (flat1.eq.90..and.flat2.eq.90.).or.
     &    (flat1.eq.-90..and.flat2.eq.-90.))  then
         del=0.
         azi=0.
         return
      end if
      pi=3.141592654
      raddeg=pi/180.
      theta1=(90.-flat1)*raddeg
      theta2=(90.-flat2)*raddeg
      phi1=flon1*raddeg
      phi2=flon2*raddeg
      stheta1=sin(theta1)
      stheta2=sin(theta2)
      ctheta1=cos(theta1)
      ctheta2=cos(theta2)
      cang=stheta1*stheta2*cos(phi2-phi1)+ctheta1*ctheta2
      ang=acos(cang)
      del=ang/raddeg
      sang=sqrt(1.-cang*cang)
      caz=(ctheta2-ctheta1*cang)/(sang*stheta1)
      saz=-stheta2*sin(phi1-phi2)/sang
      az=atan2(saz,caz)
      azi=az/raddeg
      if (azi.lt.0.) azi=azi+360.
      return
      end


! SPH_DISTDP computes angular separation of two points on sphere
!          (double precision version with real*4 i/o)
!
! Inputs:  flat1  =  latitude of first point (degrees) 
!          flon2  =  longitude of first point (degrees)
!          flat2  =  latitude of second point (degrees)
!          flon2  =  longitude of second point (degrees)
! Returns: del    =  angular separation between points (degrees)
!
! Note:  This routine is inaccurate for del less than about 0.01 degrees.
!
      subroutine SPH_DISTDP(flat1,flon1,flat2,flon2,del)
      implicit double precision (a-h,o-z)
      real*4 flat1,flon1,flat2,flon2,del
      pi=3.14159265358979267
      raddeg=pi/dble(180.)
      theta1=(dble(90.)-dble(flat1))*raddeg
      theta2=(dble(90.)-dble(flat2))*raddeg
      phi1=dble(flon1)*raddeg
      phi2=dble(flon2)*raddeg
      ang=dacos( dsin(theta1)*dsin(theta2)*dcos(phi2-phi1)+
     &          dcos(theta1)*dcos(theta2) )
      del=real(ang/raddeg)
      return
      end


! SPH_AZIDP computes distance and azimuth between two points on sphere
!          (double precision version with real*4 i/o)
!
! Inputs:  flat1  =  latitude of first point (degrees) 
!          flon2  =  longitude of first point (degrees)
!          flat2  =  latitude of second point (degrees)
!          flon2  =  longitude of second point (degrees)
! Returns: del    =  angular separation between points (degrees)
!          azi    =  azimuth at 1st point to 2nd point, from N (deg.)
!
! Note:  This routine is inaccurate for del less than about 0.01 degrees.
!
      subroutine SPH_AZIDP(flat1,flon1,flat2,flon2,del,azi)
      implicit double precision (a-h,o-z)
      real*4 flat1,flon1,flat2,flon2,del,azi
      if ((flat1.eq.flat2.and.flon1.eq.flon2).or.
     &    (flat1.eq.90..and.flat2.eq.90.).or.
     &    (flat1.eq.-90..and.flat2.eq.-90.))  then
         del=0.
         azi=0.
         return
      end if
      pi=3.14159265358979267
      raddeg=pi/dble(180.)
      theta1=(dble(90.)-dble(flat1))*raddeg
      theta2=(dble(90.)-dble(flat2))*raddeg
      phi1=dble(flon1)*raddeg
      phi2=dble(flon2)*raddeg
      stheta1=dsin(theta1)
      stheta2=dsin(theta2)
      ctheta1=dcos(theta1)
      ctheta2=dcos(theta2)
      cang=stheta1*stheta2*cos(phi2-phi1)+ctheta1*ctheta2
      ang=dacos(cang)
      del=real(ang/raddeg)
      sang=dsqrt(1.d0-cang*cang)
      caz=(ctheta2-ctheta1*cang)/(sang*stheta1)
      saz=-stheta2*dsin(phi1-phi2)/sang
      az=datan2(saz,caz)
      azi=real(az/raddeg)
      if (azi.lt.0.) azi=azi+360.
      return
      end


! SPH_GEOCENTRIC converts Earth's geographic latitude to geocentric latitude
!
!  Inputs:  flat1  =  geographic latitude (degrees) 
!                     (relative to local vertical)
!  Returns: flat2  =  geocentric latitude (degrees)
!                     (relative to Earth's center)
!
      subroutine SPH_GEOCENTRIC(flat1,flat2)
      if (flat1.eq.90.) then
         flat2=90.
         return
      end if
      pi=3.141592654
      degrad=180./pi
      factor=0.9933056      !=(1-1/298.256)**2
      phi=flat1/degrad
      theta=atan(factor*tan(phi))
      flat2=theta*degrad
      return
      end


! SPH_GEOGRAPHIC converts Earth's geocentric latitude to geographic latitude
!
!  Inputs:  flat1  =  geocentric latitude (degrees)
!                     (relative to Earth's center)
!  Returns: flat2  =  geographic latitude (degrees) 
!                     (relative to local vertical)
!
      subroutine SPH_GEOGRAPHIC(flat1,flat2)
      if (flat1.eq.90.) then
         flat2=90.
         return
      end if
      pi=3.141592654
      degrad=180./pi
      factor=0.9933056      !=(1-1/298.256)**2
      theta=flat1/degrad
      phi=atan(tan(theta)/factor)
      flat2=phi*degrad
      return
      end


! SPH_CART converts from spherical (lat,lon,r) to cartesian (x,y,z)
!
! Inputs:   flat  =  latitude (degrees)
!           flon  =  longitude (degrees)
!           r     =  radius
! Returns:  x,y,z =  x,y,z coordinates (z=up)
!
      subroutine SPH_CART(flat,flon,r,x,y,z)
      degrad=180./3.1415927
      theta=(90.-flat)/degrad
      phi=flon/degrad
      stheta=sin(theta)
      x=stheta*cos(phi)*r
      y=stheta*sin(phi)*r
      z=cos(theta)*r
      return
      end


! SPH_SPH converts from cartesian (x,y,z) to spherical (lat,lon,r)
!
! Inputs:   x,y,z =  x,y,z coordinates (z=up)
! Returns:  flat  =  latitude (degrees)
!           flon  =  longitude (degrees)
!           r     =  radius
!
      subroutine SPH_SPH(x,y,z,flat,flon,r)
      degrad=180./3.1415927
      r=sqrt(x*x+y*y+z*z)
      if (r.eq.0.) then
         flat=0.
         flon=0.
         return
      end if
      d=sqrt(x*x+y*y)
      theta=atan2(d,z)
      phi=atan2(y,x)
      flat=90.-theta*degrad
      flon=phi*degrad
      if (flon.lt.0.) flon=flon+360.
      return
      end


! SPH_POLE finds pole of great circle joining two points on sphere
!
! Requires:  SPH_CART, SPH_SPH
!      
! Inputs:  flat1  =  latitude of first point (degrees) 
!          flon2  =  longitude of first point (degrees)
!          flat2  =  latitude of second point (degrees)
!          flon2  =  longitude of second point (degrees)
! Returns: flat3  =  latitude of great circle pole (degrees)
!          flon3  =  longitude of great circle pole (degrees)
!
      subroutine SPH_POLE(flat1,flon1,flat2,flon2,flat3,flon3)
      call SPH_CART(flat1,flon1,1.,x1,y1,z1)
      call SPH_CART(flat2,flon2,1.,x2,y2,z2)
      x3=y1*z2-y2*z1
      y3=x2*z1-x1*z2
      z3=x1*y2-x2*y1
      call SPH_SPH(x3,y3,z3,flat3,flon3,r)
      return
      end


! SPH_TO_GCIRC rotates (lat,lon) to new coordinates using
! specified great circle as new equator 
!
! Requires:  SPH_CART, SPH_SPH
!
! Inputs:   flat  =  latitude of point (degrees)
!           flon  =  longitude of point (degrees)
!           plat  =  latitude of great circle pole (degrees)
!           plon  =  longitude of great circle pole (degrees)
! Returns:  glat  =  latitude of point relative to great circle (degrees)
!           glon  =  longitude of point relative to great circle (degrees)
!
      subroutine SPH_TO_GCIRC(flat,flon,plat,plon,glat,glon)
      degrad=180./3.1415927
      flon1=flon-plon
      if (flon1.lt.0.) flon1=flon1+360.
      call SPH_CART(flat,flon1,1.,x,y,z)
      pcolat=90.-plat      
      the=-pcolat/degrad
      costhe=cos(the)
      sinthe=sin(the)
      yr=y
      xr=costhe*x+sinthe*z
      zr=costhe*z-sinthe*x
      call SPH_SPH(xr,yr,zr,glat,glon,r)
      return
      end


! SPH_FROM_GCIRC rotates (lat,lon) from great circle coordinates
! back to original coordinates
!
! Requires:  SPH_CART, SPH_SPH
!
! Inputs:   glat  =  latitude of point relative to great circle (degrees)
!           glon  =  longitude of point relative to great circle (degrees)
!           plat  =  latitude of great circle pole (degrees)
!           plon  =  longitude of great circle pole (degrees)
! Returns:  flat  =  latitude of point (degrees)
!           flon  =  longitude of point (degrees)
!
      subroutine SPH_FROM_GCIRC(glat,glon,plat,plon,flat,flon)
      degrad=180./3.1415927
      call SPH_CART(glat,glon,1.,xr,yr,zr)
      pcolat=90.-plat
      the=-pcolat/degrad
      costhe=cos(the)
      sinthe=sin(the)
      y=yr
      x=costhe*xr-sinthe*zr
      z=costhe*zr+sinthe*xr
      call SPH_SPH(x,y,z,flat,flon1,r)
      flon=flon1+plon
      if (flon.gt.360.) flon=flon-360.
      return
      end
