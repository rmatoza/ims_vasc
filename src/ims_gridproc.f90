      program ims_gridproc
      ! This file is part of IMS_VASC: combined infrasound signal
      ! association and source location using a brute-force,
      ! grid-search, cross-bearings approach
      !
      ! IMS_GRIDPROC performs a basic clutter correction by grid differencing
      ! a current time snapshot with a prior grid [Matoza et al., 2017]
      !
      ! written by Robin S. Matoza, University of California, Santa Barbara
      ! email: rmatoza@ucsb.edu

      implicit none
      
      integer i, j
      real  tau_p_des, tau_d_des
      real  tau_p, tau_d
      real alpha
      character (len=100) :: outfile_during, outfile_prior
      character (len=100) :: outfile_clean, outfile_mask
      character (len=100) :: outfile_catalog
      integer ncross_usr
      real azgap_usr
      integer, dimension(:,:), allocatable :: grid_clean, grid_mask, grid_thresh
      integer iwrite 
      integer pixthreshdet
      integer, dimension(2) :: iloc
      integer igmax
      real*4 slat, slon
      integer isexceed
      
      ! during grid
      integer nlat_d, nlon_d
      real*4, dimension(:), allocatable :: latvec_d, lonvec_d
      integer, dimension(:,:), allocatable :: grid_stack_d
      integer, dimension(:,:), allocatable :: grid_stacount_d
      real, dimension(:,:), allocatable :: grid_azigap_d
      integer npixcross_d
      integer yr_on_d, mon_on_d, dy_on_d, hr_on_d, mn_on_d
      real sec_on_d
      integer yr_off_d, mon_off_d, dy_off_d, hr_off_d, mn_off_d
      real sec_off_d
      real*8 tday_son_d, tday_soff_d
      
      ! prior grid
      integer nlat_p, nlon_p
      real*4, dimension(:), allocatable :: latvec_p, lonvec_p
      integer, dimension(:,:), allocatable :: grid_stack_p
      integer, dimension(:,:), allocatable :: grid_stacount_p
      real, dimension(:,:), allocatable :: grid_azigap_p
      integer npixcross_p
      integer yr_on_p, mon_on_p, dy_on_p, hr_on_p, mn_on_p
      real sec_on_p
      integer yr_off_p, mon_off_p, dy_off_p, hr_off_p, mn_off_p
      real sec_off_p
      real*8 tday_son_p, tday_soff_p
      
      ! volcano catalaog file
      character (len=100) :: volcfile
      integer, parameter :: nvolc_max = 2000   ! maximum number of volcanoes
      integer nvolc
      real*4, dimension(nvolc_max) :: vlat, vlon, velev
      integer, dimension(nvolc_max) ::  vnum
      character (len=50), dimension(nvolc_max) :: vname
      integer ios, isel
      integer, dimension(1) :: ivmin
      real*4 del_volc, azi_volc, vmindist
      real*4, dimension(nvolc_max) :: del_volckm=0.
      real degrad, degkm
      

      degrad = 180./3.1415927
      degkm = 111.19
      
      ! FORMATS
300   format (f8.3,1x,f10.3,1x,f5.0,1x,i6,a50) ! volcano database file      
500   format(f8.3,1x,f10.3,1x,i10) ! combined xyz file (ascii)
800   format(i4,1x,i0.2,1x,i0.2,1x,i0.2,1x,i0.2,1x,f7.3,1x,f10.3,1x,&
      f10.3,1x,f10.3,1x,i10,1x,i3,1x,f10.3,4x,'|*g*:|',1x,i10,1x,&
      i10,1x,i10,1x,i10,1x,f5.1,1x,f5.1,4x,'|*v*:|',1x,f10.3,1x,f10.3,&
      1x,f8.1,1x,i6,1x,i4,a50) ! output catalog file

      print *,'Enter binary input filename (during grid)'
      read (*,'(a)') outfile_during
      print *, trim(outfile_during)
      
      print *,'Enter binary input filename (prior grid)'
      read (*,'(a)') outfile_prior
      print *, trim(outfile_prior)
      
      print *,'Enter duration of during grid (days)'
      read *, tau_d_des
      print *, tau_d_des
      
      print *,'Enter duration of prior grid (days)'
      read *, tau_p_des
      print *, tau_p_des

      print *,'Enter alpha, background rate exceed ratio, e.g., 1.5 means detection must exceed background by 1.5 times'
      read *, alpha
      print *, alpha

      print *,'Enter minimum number of interesecting stations (e.g., 3)'
      read *, ncross_usr
      print *, ncross_usr
      
      print *,'Enter maximum azimuthal gap [deg] (e.g., 200.)'
      read *, azgap_usr
      print *, azgap_usr
      
      print *,'Enter binary output filename (cleaned and masked grid) or [none]'
      read (*,'(a)') outfile_clean
      print *, trim(outfile_clean)
      
      iwrite = 1
      if (trim(outfile_clean) .eq. 'none') iwrite = 0
      
      print *,'Enter # pixels threshold for detecion (e.g., 30000)'
      read *, pixthreshdet
      print *, pixthreshdet
      
      ! Read the volcano list file (Smithsonian Institution database)
      print *,'Enter volcano list file'
      read (*,'(a)') volcfile
      print *, 'Reading: ', trim(volcfile)
      open (11,file=volcfile,status='old')
      i = 0
      ios = 0
      do while (ios==0)	
		i=i+1	
        read(11,300,iostat=ios) vlat(i), vlon(i), velev(i), vnum(i), vname(i)
      enddo
      close(11)
      nvolc = i - 1
      print *, 'Number of volcanoes in database: ', nvolc
      if (nvolc .gt. nvolc_max) then
          print *, '**ERROR: nvolc .gt. nvolc_max, incease nvolc_max', nvolc, nvolc_max
          stop
      endif	
	  
      ! Example access to volcano database
      !isel = 766 ! should be Sarychev Peak
      !print *, vlat(isel), vlon(isel), velev(isel), vnum(isel), trim(vname(isel))
      
      print *,'Event detection file'
      read (*,'(a)') outfile_catalog
      print *, trim(outfile_catalog)
      
      ! During grid
      print *, 'DURING GRID'
      open (13,file=trim(outfile_during),form='unformatted')
      print *, 'Reading binary file: ', trim(outfile_during)
      read(13) nlat_d, nlon_d, npixcross_d
      print *, "nlat, nlon, npixcross: ", nlat_d, nlon_d, npixcross_d
      
      print *, "Allocating arrays"
      allocate (latvec_d(nlat_d))
      allocate (lonvec_d(nlon_d))
      allocate (grid_stack_d(nlat_d, nlon_d))
      allocate (grid_stacount_d(nlat_d, nlon_d)) 
      allocate (grid_azigap_d(nlat_d, nlon_d))   
      
      print *, "Reading arrays"
      read(13) latvec_d(1:nlat_d), lonvec_d(1:nlon_d)
      read(13) grid_stack_d(1:nlat_d,1:nlon_d)
      read(13) grid_stacount_d(1:nlat_d,1:nlon_d)
      read(13) grid_azigap_d(1:nlat_d,1:nlon_d)
      read(13) yr_on_d, mon_on_d, dy_on_d, hr_on_d, mn_on_d, sec_on_d, tday_son_d
      read(13) yr_off_d, mon_off_d, dy_off_d, hr_off_d, mn_off_d, sec_off_d, tday_soff_d
          
      close(13)
      
      tau_d = tday_soff_d - tday_son_d !last day doesn't count; buffer on celerity, no +1
      print *, yr_on_d, mon_on_d, dy_on_d, hr_on_d, mn_on_d, sec_on_d
      print *, yr_off_d, mon_off_d, dy_off_d, hr_off_d, mn_off_d, sec_off_d
      print *, tday_son_d, tday_soff_d, tau_d ! note that we can use this for grid processing
      
      if (tau_d .ne. tau_d_des) then
        print *, '**Error: tau_d does not match tau_d_des: ', tau_d, tau_d_des
        stop
      endif
      
      ! Prior grid
      print *, 'PRIOR GRID'
      open (14,file=trim(outfile_prior),form='unformatted')
      print *, 'Reading binary file: ', trim(outfile_prior)
      read(14) nlat_p, nlon_p, npixcross_p
      print *, "nlat, nlon, npixcross: ", nlat_p, nlon_p, npixcross_p
      
      if (nlat_p .ne. nlat_d) then
        print *, '**Error: nlat_p does not match nlat_d: ', nlat_p, nlat_d
        stop
      endif
      
      if (nlon_p .ne. nlon_d) then
        print *, '**Error: nlon_p does not match nlon_d: ', nlon_p, nlon_d
        stop
      endif
      
      if (npixcross_p .ne. npixcross_d) then
        print *, '**Error: npixcross_p does not match npixcross_d: ', npixcross_p, npixcross_d
        stop
      endif
      
      print *, "Allocating arrays"
      allocate (latvec_p(nlat_p))
      allocate (lonvec_p(nlon_p))
      allocate (grid_stack_p(nlat_p, nlon_p))
      allocate (grid_stacount_p(nlat_p, nlon_p)) 
      allocate (grid_azigap_p(nlat_p, nlon_p))   
      
      print *, "Reading arrays"
      read(14) latvec_p(1:nlat_p), lonvec_p(1:nlon_p)
      read(14) grid_stack_p(1:nlat_p,1:nlon_p)
      read(14) grid_stacount_p(1:nlat_p,1:nlon_p)
      read(14) grid_azigap_p(1:nlat_p,1:nlon_p)
      read(14) yr_on_p, mon_on_p, dy_on_p, hr_on_p, mn_on_p, sec_on_p, tday_son_p
      read(14) yr_off_p, mon_off_p, dy_off_p, hr_off_p, mn_off_p, sec_off_p, tday_soff_p
          
      close(14)
      
      tau_p = tday_soff_p - tday_son_p!last day doesn't count; buffer on celerity, no +1
      print *, yr_on_p, mon_on_p, dy_on_p, hr_on_p, mn_on_p, sec_on_p
      print *, yr_off_p, mon_off_p, dy_off_p, hr_off_p, mn_off_p, sec_off_p
      print *, tday_son_p, tday_soff_p, tau_p ! note that we can use this for grid processing
      
      if (tau_p .ne. tau_p_des) then
        print *, '**Error: tau_p does not match tau_p_des: ', tau_p, tau_p_des
        stop
      endif
      
      ! CLEAN and MASK GRID
      allocate (grid_clean(nlat_p, nlon_p))
      allocate (grid_mask(nlat_p, nlon_p))
      allocate (grid_thresh(nlat_p, nlon_p))
      grid_clean = 0
      grid_mask = 0
      grid_thresh = 0
      
      grid_clean = grid_stack_d - alpha*(tau_d/tau_p)*grid_stack_p
      
      isexceed = 0
      do i = 1,nlat_d
         do j = 1,nlon_d
            if ((grid_stacount_d(i,j) .ge. ncross_usr) .and. &
                   (grid_azigap_d(i,j) .le. azgap_usr) .and. (grid_azigap_d(i,j) .gt. 0.)) then
               
               grid_mask(i,j) = grid_clean(i,j)
               if (grid_mask(i,j) .ge. pixthreshdet) then
                  grid_thresh(i,j) = grid_mask(i,j)
                  isexceed = 1
               endif
            endif
         enddo
      enddo
      
      ! OUTPUT FILE
      !ASCII FORMAT OUTPUT FILES
      !open (15,file=trim(outfile_clean),status='new')
      !print *, 'Writing: ', trim(outfile_clean)
      !do i = 1,nlat_d
      !  do j = 1,nlon_d
      !     write(15,500) latvec_d(i), lonvec_d(j), grid_clean(i,j)
      !  enddo
      !enddo 
      !close(15)
      
      !open (16,file=trim(outfile_mask),status='new')
      !print *, 'Writing: ', trim(outfile_mask)
      !do i = 1,nlat_d
      !  do j = 1,nlon_d
      !     write(16,500) latvec_d(i), lonvec_d(j), grid_mask(i,j)
      !  enddo
      !enddo 
      !close(16)
      
      ! binary output file
      if (iwrite .eq. 1) then
          open (15,file=trim(outfile_clean),status='new', form='unformatted')
          print *, 'Writing binary file: ', trim(outfile_clean)
          write(15) nlat_d, nlon_d, npixcross_d
          write(15) latvec_d(1:nlat_d), lonvec_d(1:nlon_d)
          write(15) grid_thresh(1:nlat_d,1:nlon_d)
          write(15) grid_mask(1:nlat_d,1:nlon_d)
          write(15) grid_clean(1:nlat_d,1:nlon_d)
          write(15) grid_stacount_d(1:nlat_d,1:nlon_d)
          write(15) grid_azigap_d(1:nlat_d,1:nlon_d)
          write(15) yr_on_d, mon_on_d, dy_on_d, hr_on_d, mn_on_d, sec_on_d
          write(15) yr_off_d, mon_off_d, dy_off_d, hr_off_d, mn_off_d, sec_off_d
          write(15) yr_on_p, mon_on_p, dy_on_p, hr_on_p, mn_on_p, sec_on_p
          write(15) yr_off_p, mon_off_p, dy_off_p, hr_off_p, mn_off_p, sec_off_p
          close(15)
      endif

      ! event detection
      if (isexceed .eq. 1) then
          iloc = maxloc(grid_thresh)
          igmax = maxval(grid_thresh)
          
          print *, 'Max loc: ', iloc
          print *, 'Max val: ', igmax
          slat = latvec_d(iloc(1))
          slon = lonvec_d(iloc(2))
          print *, 'Location: ', slat, slon
          
          ! look for nearby volcanoes:
          del_volckm=99999.
          do isel=1,nvolc
             CALL SPH_AZIDP(slat,slon,vlat(isel),vlon(isel),del_volc,azi_volc) ! note this azi is from source to station
             del_volckm(isel) = degkm*del_volc
          enddo
          vmindist = minval(del_volckm(1:nvolc))
          ivmin = minloc(del_volckm(1:nvolc))
          
          print *, 'Min distance to volcano [km]: ', vmindist 
          print *, 'Min loc: ', ivmin
          print *, vlat(ivmin), vlon(ivmin), velev(ivmin), vnum(ivmin), trim(vname(ivmin(1)))
          
          open (17,file=trim(outfile_catalog),status='new')
          print *, 'Writing: ', trim(outfile_catalog)
          write(17,800) yr_on_d, mon_on_d, dy_on_d, hr_on_d, mn_on_d, sec_on_d, &
                      slat, slon, vmindist, igmax, grid_stacount_d(iloc(1),iloc(2)), &
                      grid_azigap_p(iloc(1),iloc(2)), grid_thresh(iloc(1),iloc(2)), &
                      grid_clean(iloc(1),iloc(2)), grid_stack_d(iloc(1),iloc(2)), &
                      grid_stack_p(iloc(1),iloc(2)), &
                      tau_d, tau_p, &
                      vlat(ivmin), vlon(ivmin), velev(ivmin), &
                      vnum(ivmin), ivmin(1), trim(vname(ivmin(1)))
          close(17)
      endif
      print *, isexceed
      if (isexceed .eq. 0) print *, 'No event detected above threshold: ', pixthreshdet      

end program ims_gridproc

