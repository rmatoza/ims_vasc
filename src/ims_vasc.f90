! IMS_VASC: combined infrasound signal association and source location
! using a brute-force, grid-search, cross-bearings approach
!
! written by Robin S. Matoza, University of California, Santa Barbara
! email: rmatoza@ucsb.edu
!
! If you use this software, please cite:
!
! Matoza, R.S., D.N. Green, A. Le Pichon, P.M. Shearer, D. Fee,
!     P. Mialle, and L. Ceranna (2017), Automated detection and
!     cataloging of global explosive volcanism using the International
!     Monitoring System infrasound network, J. Geophys. Res.
!     Solid Earth, 122, 2946–2971, doi:10.1002/2016JB013356
!
!
! Modified BSD-2 License - for Non-Commercial Research and Educational Use Only
!
! Copyright (c) 2017, The Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted for non-commercial research and educational
! use only provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! For permission to use for commercial purposes, please contact
! UCSB’s Office of Technology & Industry Alliances at 805-893-5180
! or info@tia.ucsb.edu.
!
      
      module bulltype ! define types for storing infrasound array detection files (PMCC bulletin files)
      
      ! array dimensions
      integer, parameter :: nfam_max = 300000  ! maximum number of families (lines in bulletin input file)
      integer, parameter :: nsta_max = 60      ! maximum number of stations
      type bulldata
         character (len=4) :: fstname
         integer :: nfam
         real, dimension(nfam_max) :: azi
         real, dimension(nfam_max) :: fmin, fmax, fmean
         real*8, dimension(nfam_max) :: tday  ! days since 00:00 1 Jan 1600
         integer, dimension(nfam_max) :: npix ! for bulletin files, this is npix, the number of pixels in each family.
      end type bulldata
      
      end module bulltype
      
      program ims_vasc
	        
      use bulltype
      implicit none
      
      ! array dimensions
      integer nsta, nbullfile      ! actual sizes
	  
      ! general
      character (len=100) :: stafile, bulllistfile
      character (len=100) :: outfile_all, outfile_all_pref
      integer i, j, k, p, isel, isnap
      integer ios
	  
      ! stations
      real*4, dimension(nsta_max) :: stlat, stlon
      character (len=4), dimension(nsta_max) :: stname
	  
      ! family files
      character (len=100), dimension(nsta_max) :: bullfile
      character (len=4), dimension(nsta_max) :: stbullfile
	  
      ! grid search variables
      real*4 lat_min, lat_max, dlat
      real*4 lon_min, lon_max, dlon
      real az_dev
      real celer
      real distmax
      integer npixcross
      integer nlat, nlon
      real freq_selmin, freq_selmax, freq_selmean
      integer, dimension(:,:), allocatable :: grid_stack, grid_sum
      integer isdet
      integer, dimension(:,:), allocatable :: grid_stacount
      real*4, dimension(:), allocatable :: latvec, lonvec
      real*4 lat_trial, lon_trial, azi_trial, del_trial, del_trialkm, azi_hi, azi_lo
      logical, dimension(:), allocatable :: ismatch
      
      ! azimuthal gap 
      real*4, dimension(nsta_max) :: azivec=0.
      real*4 azigap
      integer iazi, mazi
      real*4 azi_back, del_back
      real, dimension(:,:), allocatable :: grid_azigap
      
      !distance sorting
      integer nnearsta, kdum
      real*4, dimension(nsta_max) :: delstavec=0.
      integer, dimension(nsta_max) :: indxdel
      
      ! time-evolving
      integer yr_on, mon_on, dy_on, hr_on, mn_on
      real sec_on
      integer yr_off, mon_off, dy_off, hr_off, mn_off
      real sec_off
      real*8 tday_son, tday_soff, t1, t2
      real*8 tback
      character (len=10) charsnap 

      ! Bulletin file arrays
      integer yr, mon, dy, hr, mn
      real hr_r, mn_r
      real sec, fmin, fmax, azi, vapp, arms, fmean
      integer npix
      integer itday
      real*8 dfrac
      
      type(bulldata), dimension(nsta_max) :: STADET, STARAW ! station detection lists (bulletin files)
      
      ! distances
      real degrad, degkm
      real*8 tshift
      real*4 latgeocen1, latgeocen2
      
      ! WAITBAR
      real pct
      integer pct_i
      character (len=13), dimension(11) :: waitbar
      waitbar(1)  = "|          |"
      waitbar(2)  = "|#         |"
      waitbar(3)  = "|##        |"
      waitbar(4)  = "|###       |"
      waitbar(5)  = "|####      |"
      waitbar(6)  = "|#####     |"
      waitbar(7)  = "|######    |"
      waitbar(8)  = "|#######   |"
      waitbar(9)  = "|########  |"
      waitbar(10) = "|######### |"
      waitbar(11) = "|##########|"
      
      degrad = 180./3.1415927
      degkm = 111.19
      
      ! FORMATS
100   format (a4,1x,a100) ! bulletin file list
200   format (a12,5x,f5.1,a3,i10,a4,i10) ! waitbar
300   format (f8.3,1x,f10.3,1x,f5.0,1x,i6,a50) ! volcano database file
400   format(1000i10) ! matrix output 1000 is a large number > nlon; better for ifort
500   format(f8.3,1x,f10.3,1x,i10)  ! xyz file, easier for GMT xyz2grd
600   format(f8.3,1x,f10.3,1x,f5.1) ! xyz file, easier for GMT xyz2grd
700   format(f8.3,1x,f10.3,1x,i10,1x,i10,1x,f10.3) ! combined xyz file (ascii)


      print *, ''
      print *, ' ___ __  __ ____ __     ___    ____   ____  '
      print *, '|_ _|  \/  / ___|\ \   / / \  / ___| / ___| '
      print *, ' | || |\/| \___ \ \ \ / / _ \ \___ \| |     '
      print *, ' | || |  | |___) | \ V / ___ \ ___) | |___  '
      print *, '|___|_|  |_|____/___\_/_/   \_\____/ \____| '
      print *, '               |_____|                      '
      print *, ''

      print *, ''
      print *, 'Combined infrasound signal association and source location &
                 using a brute-force, grid-search, cross-bearings approach.'
      print *, ''
      print *, 'Robin S. Matoza, University of California, Santa Barbara'
      print *, 'email: rmatoza@ucsb.edu'
      print *, ''
      print *, ''
      print *, 'Modified BSD-2 License - for Non-Commercial Research and Educational Use Only'
      print *, 'Copyright (c) 2017, The Regents of the University of California'
      print *, 'All rights reserved.'
      print *, ''
      print *, ''

      !*****************************************************************
      ! READ PARAMETERS/INPUT FILES                                    *
      !*****************************************************************

      print *,'Enter station list file'
      read (*,'(a)') stafile
      print *, 'Reading: ', trim(stafile)
      open (11,file=stafile,status='old')
	  
      ! Read the station file
      i = 0
      ios = 0
      do while (ios==0)	
          i=i+1	
          read(11,*,iostat=ios) stlat(i),stlon(i),stname(i)
          if (ios < 0) then
             exit
          else if (ios > 0) then
             print *, '***ERROR, read error on line: ', i, trim(stafile)
             stop
          endif
      enddo
      close(11)
      nsta = i - 1
      print *, 'Number of stations read: ', nsta	
      if (nsta .gt. nsta_max) then
          print *, '**ERROR: nsta .gt. nsta_max, incease nsta_max', nsta, nsta_max
          stop
      endif
      print *,'Enter min lat, max lat, dlat, for search area: '
      read *, lat_min, lat_max, dlat
      print *, 'Read: ', lat_min, lat_max, dlat
	  
      print *,'Enter min lon, max lon, dlon, for search area: '
      read *, lon_min, lon_max, dlon
      print *, 'Read: ', lon_min, lon_max, dlon
	  
      print *,'Enter allowed azimuth deviation [deg]: '
      read *, az_dev
      print *, 'Read: ', az_dev
	  
      print *,'Enter association celerity [km/s]: '
      read *, celer
      print *, 'Read: ', celer
      
      print *,'Enter association max distance [km]: '
      read *, distmax
      print *, 'Read: ', distmax
      
      print *,'Enter min, max, mean frequency [Hz], [enter 0 0 0 for all frequencies]: '
      read *, freq_selmin, freq_selmax, freq_selmean
      
      print *,'Enter min number of pixels to count a station in intersection: '
      read *, npixcross
	  print *, 'Read: ', npixcross
      
      print *,'Enter number of nearest stations that must detect (e.g., 2) [enter 0 to not use this option]'
      read *, nnearsta
      print *, nnearsta
      if (nnearsta .eq. 0) print *, 'Nearest stations constraint not used'
	  
      print *,'Enter bulletin list file'
      read (*,'(a)') bulllistfile
      print *, 'Reading: ', trim(bulllistfile)
      open (11,file=bulllistfile,status='old')
	  
      ! Read the bulletin list file
      i = 0
      ios = 0
      do while (ios==0)	
          i=i+1	
          read(11,100,iostat=ios) stbullfile(i),bullfile(i)
          if (ios < 0) then
             exit
          else if (ios > 0) then
             print *, '***ERROR, read error on line: ', i, trim(bulllistfile)
             stop
          endif
          bullfile(i) = trim(bullfile(i))
          if (stname(i) .ne. stbullfile(i)) then
             print *, '**ERROR: station coordinates and bulletin file orders do not match: ', stname(i), '  ', stbullfile(i)
             stop
          endif
      enddo
      close(11)
      nbullfile = i - 1
      print *, 'Number of bulletin files in list: ', nbullfile
      if (nbullfile .ne. nsta) print *, '***WARNING: nbullfile .ne. nsta', nbullfile, nsta	
      if (nbullfile .gt. nsta) then
          print *, '**ERROR: nbullfile .gt. nsta_max, incease nsta_max', nbullfile, nsta_max
          stop
      endif	
      
      print *,'Enter stacking start time [yr, mon, dy, hr, mn, sec]: '
      read *, yr_on, mon_on, dy_on, hr_on, mn_on, sec_on
      print *, 'Read: ', yr_on, '[year]'
      print *, 'Read: ', mon_on, '[month]'
      print *, 'Read: ', dy_on, '[day]'
      print *, 'Read: ', hr_on, '[hour]'
      print *, 'Read: ', mn_on, '[min]'
      print *, 'Read: ', sec_on, '[sec]'
      CALL DT_GET_TDAY(yr_on, mon_on, dy_on, itday)
      dfrac = float(hr_on)/24 + float(mn_on)/(24*60) + sec_on/(24*60*60)
      tday_son = float(itday) + dfrac
      
      print *,'Enter stacking stop time [yr, mon, dy, hr, mn, sec]: '
      read *, yr_off, mon_off, dy_off, hr_off, mn_off, sec_off
      print *, 'Read: ', yr_off, '[year]'
      print *, 'Read: ', mon_off, '[month]'
      print *, 'Read: ', dy_off, '[day]'
      print *, 'Read: ', hr_off, '[hour]'
      print *, 'Read: ', mn_off, '[min]'
      print *, 'Read: ', sec_off, '[sec]'
      CALL DT_GET_TDAY(yr_off, mon_off, dy_off, itday)
      dfrac = float(hr_off)/24 + float(mn_off)/(24*60) + sec_off/(24*60*60)
      tday_soff = float(itday) + dfrac
      
      ! Read the bulletin files
      do j = 1,nbullfile
         print *, 'READING: ', trim(bullfile(j))
         STARAW(j)%fstname = stbullfile(j)       
         open (12,file=trim(bullfile(j)),status='old')    
         i = 0
         ios = 0
         do while (ios==0)	
            i=i+1	
            if (i .gt. nfam_max) then
               print *, '**ERROR: too many families, increase nfam_max', i, nfam_max
               stop
            endif    
            read(12,*,iostat=ios) yr, mon, dy, hr, mn, sec, fmin, fmax, azi, vapp, arms, npix, fmean
            if (ios < 0) then
               exit
            else if (ios > 0) then
               print *, '***ERROR, read error on line: ', i, trim(bullfile(j))
               stop
            endif
            CALL DT_GET_TDAY(yr, mon, dy, itday)
            dfrac = float(hr)/24 + float(mn)/(24*60) + sec/(24*60*60)
            STARAW(j)%tday(i) = float(itday) + dfrac
            STARAW(j)%azi(i) = azi
            STARAW(j)%fmin(i) = fmin
            STARAW(j)%fmax(i) = fmax
            STARAW(j)%fmean(i) = fmean
            STARAW(j)%npix(i) = npix
         enddo
         STARAW(j)%nfam = i-1
         if (STARAW(j)%nfam .gt. nfam_max) then
            print *, '**ERROR: too many input data, increase nfam_max', STARAW(j)%fstname, STARAW(j)%nfam, nfam_max
            stop
         endif
         print *, 'READ # families: ', STARAW(j)%fstname, STARAW(j)%nfam
         close(12)
      enddo
      
      ! Parse the detections based on frequency content
      if ((freq_selmin.ne.0) .and. (freq_selmax.ne.0)) then
         print *, 'Applying frequency parsing [Hz]: ', freq_selmin, freq_selmax, freq_selmean
         CALL FCLNDET(STARAW,freq_selmin,freq_selmax,freq_selmean,nbullfile,STADET)
      else if ((freq_selmin.eq.0) .and. (freq_selmax.eq.0)) then
         print *, 'No frequency parsing, all detections will be used'
         do j = 1,nbullfile
            STADET(j) = STARAW(j)
         enddo
      endif

      ! example of how to access
      !do i =1,100
      !   print *, 'Enter index (-99 to quit): '
      !   read *, j
      !   if (j .eq. -99) exit
      !  CALL DT_GET_YMD(floor(STADET(1)%tday(j)),yr, mon, dy)
      !  print *, STADET(1)%fstname, STADET(1)%azi(j), STADET(1)%tday(j), STADET(1)%nfam, yr, mon, dy
      !enddo
      
      print *,'Enter binary output filename prefix (lat, lon, G, stadet, azigap)'
      read (*,'(a)') outfile_all_pref
      outfile_all = trim(outfile_all_pref)  // ".dat"
      
      !*****************************************************************
      ! GRID SEARCH                                                    *
      !*****************************************************************
      
      ! Allocate grid_stack, latvec, lonvec
      nlat = (lat_max-lat_min)/dlat + 1
      nlon = (lon_max-lon_min)/dlon + 1
      print *, 'Allocating grid_stack nlat x nlon: ', nlat, nlon
      allocate (grid_stack(nlat, nlon))
      allocate (grid_sum(nlat, nlon))
      allocate (grid_stacount(nlat, nlon)) ! counts the number of stations detecting 
      allocate (grid_azigap(nlat, nlon))   ! grid of the azimuthal gap values 
      allocate (latvec(nlat))
      allocate (lonvec(nlon))

      ! initialize
      grid_stack = 0
      grid_sum = 0 
      grid_stacount = 0
      grid_azigap = -999.
      do i = 1,nlat
        latvec(i) = lat_min+dlat*(i-1)
      enddo
      do i = 1,nlon
        lonvec(i) = lon_min+dlon*(i-1)
      enddo

      t1 = tday_son
      t2 = tday_soff
      grid_stack = 0
      grid_stacount = 0

      ! Start the main calculation
      print *, 'Calculating...'
      call timer

      do i = 1,nlat ! loop over latitude
         ! waitbar
         pct = (float(i)/float(nlat))*100
         pct_i = int(pct/10) + 1
         if(mod(i,10).eq.0) write(6,200) waitbar(pct_i), pct, ' % ', i, ' of ', nlat
         do j = 1,nlon ! loop over longitude
            lat_trial = latvec(i)
            lon_trial = lonvec(j)
            iazi = 0 ! initialize variables for azimuthal gap calculation at this (lat_trial, lon_trial)
            mazi = 0
            azivec = 0.

            ! sort stations in order of distance to *this* (lat,lon) source node
            do kdum = 1,nsta ! loop over stations
               CALL SPH_GEOCENTRIC(stlat(kdum),latgeocen1)
               CALL SPH_GEOCENTRIC(lat_trial,latgeocen2)
               CALL SPH_AZIDP(latgeocen1,stlon(kdum),latgeocen2,lon_trial,del_trial,azi_trial) ! azimuth from station to trial source location
               delstavec(kdum) = del_trial
            enddo
            CALL INDEXX(nsta, delstavec, indxdel)

            do kdum = 1,nsta ! loop over stations, in order of increasing distance from (lat_trial, lon_trial)
               k = indxdel(kdum)
               CALL SPH_GEOCENTRIC(stlat(k),latgeocen1)
               CALL SPH_GEOCENTRIC(lat_trial,latgeocen2)
               CALL SPH_AZIDP(latgeocen1,stlon(k),latgeocen2,lon_trial,del_trial,azi_trial) ! azimuth from station to trial source location
               del_trialkm = degkm*del_trial
               tshift = del_trialkm/celer ! time-shift in seconds
               azi_hi = azi_trial + az_dev
               azi_lo = azi_trial - az_dev
               isdet = 0
               if (del_trialkm .le. distmax) then ! check point within max dist
                  do p = 1,STADET(k)%nfam
                     tback = STADET(k)%tday(p) - tshift/(24*60*60)
                     if ((tback .ge. t1) .and. (tback .le. t2)) then ! check if time-shifted detection time fits in stack time window
                        if ((STADET(k)%azi(p) .ge. azi_lo) .and. (STADET(k)%azi(p) .le. azi_hi)) then ! check if azimuth is within az_dev
                           grid_stack(i,j) = grid_stack(i,j) + STADET(k)%npix(p)  ! count detection in time snapshot: note not just 1, but npix for bulletin files
                           grid_sum(i,j) = grid_sum(i,j) + STADET(k)%npix(p)      ! count cumulative detections over full data duration: note not just 1, but npix for bulletin files
                           isdet = isdet + STADET(k)%npix(p)                      ! count on *just* this station k for just this loc(i, j); note isdet=0 within loop above
                        endif
                     endif
                  enddo
               endif

               if (nnearsta .ne. 0) then            ! if nnearsta is zero, we don't do any of this
                  if (kdum .le. nnearsta) then      ! only do this check on the nnearsta closest stations (e.g., 2 closest stations)
                    if (isdet .lt. npixcross) exit ! if any of the nnearsta stations did not detect, stop working on this trial source (end loop over station)
                  endif
               endif

               if (isdet .ge. npixcross) then                 ! Check if this station associated at least npixcross pixels to this location
                  grid_stacount(i,j) = grid_stacount(i,j) + 1 ! increase grid_stacount by +1 if *this* station detected enough npixcross for this trial source i, j
                  iazi = iazi + 1
                  CALL SPH_GEOCENTRIC(lat_trial,latgeocen1)
                  CALL SPH_GEOCENTRIC(stlat(k),latgeocen2)
                  CALL SPH_AZIDP(latgeocen1,lon_trial,latgeocen2,stlon(k),del_back,azi_back) ! note this azi is from source to station
                  azivec(iazi) = azi_back                     ! store source-station azimuth values in vector for azigap calculation
               endif
            enddo  ! end loop over station
            mazi = iazi           ! azimuthal gap calculation for this (lat, lon)
            if (mazi .ge. 3) then ! only calculate azimuthal gap if .ge. 3 stations
               call getmaxazigap(azivec,mazi,azigap)
               grid_azigap(i,j) = azigap
               !print *, lat_trial, lon_trial, "gap", azigap, mazi, "vec", azivec(1:mazi)
            endif
         enddo     ! end loop over longitude
      enddo        ! end loop over latitude
      write(6,200) waitbar(11), 100., ' % ', i-1, ' of ', nlat


      open (13,file=trim(outfile_all),status='new', form='unformatted')
      print *, 'Writing binary file: ', trim(outfile_all)
      write(13) nlat, nlon, npixcross
      write(13) latvec(1:nlat), lonvec(1:nlon)
      write(13) grid_stack(1:nlat,1:nlon)
      write(13) grid_stacount(1:nlat,1:nlon)
      write(13) grid_azigap(1:nlat,1:nlon)
      write(13) yr_on, mon_on, dy_on, hr_on, mn_on, sec_on, tday_son
      write(13) yr_off, mon_off, dy_off, hr_off, mn_off, sec_off, tday_soff


      close(13)
          

      call timer
      print *, 'IMS_VASC complete.'
end program ims_vasc


   subroutine FCLNDET(STAIN,freq_selmin,freq_selmax,freq_selmean,nbullfile,STAOUT)
   ! Parse a family type variable based on freq_selmin,freq_selmax
   ! updated 09/10/2015 to include mean freq upper limit

   use bulltype
   implicit none
   
   type(bulldata), dimension(nsta_max) :: STAIN, STAOUT ! station detection lists (family files)
   real freq_selmin,freq_selmax,freq_selmean
   integer nbullfile
   integer i, j, k, icount
   
   do j=1,nbullfile
      print *, 'FCLNDET working on: ', trim(STAIN(j)%fstname)
      STAOUT(j)%fstname = STAIN(j)%fstname
      icount = 0
      do i=1,STAIN(j)%nfam
         if((STAIN(j)%fmin(i) .ge. freq_selmin) .and. (STAIN(j)%fmax(i) .le. freq_selmax) &
            .and. (STAIN(j)%fmean(i) .le. freq_selmean)) then
            icount = icount + 1
            STAOUT(j)%tday(icount) = STAIN(j)%tday(i)
            STAOUT(j)%azi(icount) = STAIN(j)%azi(i)
            STAOUT(j)%fmin(icount) = STAIN(j)%fmin(i)
            STAOUT(j)%fmax(icount) = STAIN(j)%fmax(i)
            STAOUT(j)%npix(icount) = STAIN(j)%npix(i)
         endif
      enddo
      STAOUT(j)%nfam = icount
      print *, 'Number of detections in chosen band: ', STAOUT(j)%nfam
   enddo
   
   end


   subroutine getmaxazigap(azivec,mazi,azigap)
   ! Azimuthal gap
   ! R.S. Matoza 06/09/2015
   ! Obtain the largest open azimuth (azigap) between recording stations, e.g., Sweeney (1996) 
   ! azivec is vector of azimuths of stations to the source
   ! mazi2 is 2*mazi
   
   implicit none
   integer, parameter :: nsta_max = 200
   integer mazi
   real*4, dimension(mazi) :: azivec
   real*4, dimension(nsta_max) :: azivecrep=0., azidiff=0.
   real*4 azigap
   integer i
   
   ! repeat on the circle
   azivecrep(1:mazi) = azivec(1:mazi)
   azivecrep((mazi+1):(mazi*2)) = azivec(1:mazi) + 360.
   
   ! sort values, twice around circle
   call SORT4(mazi*2,azivecrep)
    
   ! Vector of differences between sorted azimuth value and next one round circle 
   do i = 1,((mazi*2)-1)
      azidiff(i) = azivecrep(i+1)-azivecrep(i)
   enddo
   
   azigap = maxval(azidiff(1:(mazi*2))) ! get max gap between adjacent stations
   
   return 
   end
   


! Numerical recipes routines
!-----------------------------------------------------------------------
!
      SUBROUTINE MEDIAN(X,N,XMED)
      IMPLICIT NONE

      INTEGER N
      INTEGER N2

      REAL*8 XMED
      REAL*8 X(N)

      if (n.eq.0) then
!         print *,'**WARNING in MEDIAN, n = 0'
         xmed=0.
         return
      else if (n.eq.1) then
         xmed=x(1)
         return
      end if
      CALL SORT(N,X)
      N2=N/2
      IF(2*N2.EQ.N)THEN
        XMED=0.5*(X(N2)+X(N2+1))
      ELSE
        XMED=X(N2+1)
      ENDIF
      RETURN
      END
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      SUBROUTINE SORT(N,RA)
      IMPLICIT NONE 
    
      INTEGER I
      INTEGER IR
      INTEGER J
      INTEGER L
      INTEGER N

      REAL*8 RRA
      REAL*8 RA(N)
!
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END
!       
!-----------------------------------------------------------------------
!

   
   !-----------------------------------------------------------------------
      SUBROUTINE SORT4(N,RA)
      IMPLICIT NONE 
    
      INTEGER I
      INTEGER IR
      INTEGER J
      INTEGER L
      INTEGER N

      REAL*4 RRA
      REAL*4 RA(N)
!
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
        ELSE
          RRA=RA(IR)
          RA(IR)=RA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
      GO TO 10
      END

!
!-----------------------------------------------------------------------
!
subroutine mean(x,N,xmean)
implicit none
integer i, N
real*8 x(N)
real*8 xsum, xmean

xmean = 0
xsum = 0

do i = 1,N
	xsum = xsum + x(i)
enddo

xmean = xsum/N

end
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine timer
implicit none

integer if1

real diftim
real etime
real oldtim
real tnow
real tarray(2)

save oldtim
data if1/1/
!

if(if1.eq.0) goto 5
if1=0
oldtim=etime(tarray)
return
5 continue
tnow=etime(tarray)
diftim=tnow-oldtim
oldtim=tnow
print *,'elapsed seconds =',diftim
return
end
!
!-----------------------------------------------------------------------
      SUBROUTINE INDEXX(N,ARRIN,INDX)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END





