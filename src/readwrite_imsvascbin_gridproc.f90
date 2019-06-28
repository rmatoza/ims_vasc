      program readwrite_imsvasscbin_gridproc

      ! Utility program to read binary output files produced by
      ! IMS_GRIDPROC and produce an ASCII file.
      !
      ! This file is part of IMS_VASC
      !
      ! written by Robin S. Matoza, University of California, Santa Barbara
      ! email: rmatoza@ucsb.edu

      implicit none
      

      character (len=100) :: outfile_all, outfile_ascii
      integer nlat, nlon
      real*4, dimension(:), allocatable :: latvec, lonvec
      integer, dimension(:,:), allocatable :: grid_stack
      integer, dimension(:,:), allocatable :: grid_stacount
      real, dimension(:,:), allocatable :: grid_azigap
      integer, dimension(:,:), allocatable :: grid_clean, grid_mask, grid_thresh
      integer npixcross
      integer i, j
      
      integer yr_on_d, mon_on_d, dy_on_d, hr_on_d, mn_on_d
      real sec_on_d
      integer yr_off_d, mon_off_d, dy_off_d, hr_off_d, mn_off_d
      real sec_off_d
      
      integer yr_on_p, mon_on_p, dy_on_p, hr_on_p, mn_on_p
      real sec_on_p
      integer yr_off_p, mon_off_p, dy_off_p, hr_off_p, mn_off_p
      real sec_off_p
      
900   format(f8.3,1x,f10.3,1x,i10,1x,i10,1x,i10,1x,f10.3) ! combined xyz file (ascii)


      print *,'Enter binary input filename (lat, lon, G, stadet, azigap)'
      read (*,'(a)') outfile_all
      
      print *,'Enter ascii output filename'
      read (*,'(a)') outfile_ascii
      

      open (13,file=trim(outfile_all),form='unformatted')
      print *, 'Reading binary file: ', trim(outfile_all)
      read(13) nlat, nlon, npixcross
      print *, "nlat, nlon, npixcross: ", nlat, nlon, npixcross 
      
      print *, "Allocating arrays"
      allocate (latvec(nlat))
      allocate (lonvec(nlon))
      allocate (grid_stack(nlat, nlon))
      allocate (grid_stacount(nlat, nlon)) 
      allocate (grid_azigap(nlat, nlon))   
      allocate (grid_thresh(nlat, nlon))
      allocate (grid_mask(nlat, nlon))
      allocate (grid_clean(nlat, nlon))
      
      print *, "Reading arrays"
      read(13) latvec(1:nlat), lonvec(1:nlon)
      read(13) grid_thresh(1:nlat,1:nlon)
      read(13) grid_mask(1:nlat,1:nlon)
      read(13) grid_clean(1:nlat,1:nlon)
      read(13) grid_stacount(1:nlat,1:nlon)
      read(13) grid_azigap(1:nlat,1:nlon)
      read(13) yr_on_d, mon_on_d, dy_on_d, hr_on_d, mn_on_d, sec_on_d
      read(13) yr_off_d, mon_off_d, dy_off_d, hr_off_d, mn_off_d, sec_off_d
      read(13) yr_on_p, mon_on_p, dy_on_p, hr_on_p, mn_on_p, sec_on_p
      read(13) yr_off_p, mon_off_p, dy_off_p, hr_off_p, mn_off_p, sec_off_p
      
      close(13)
      
      print *, yr_on_d, mon_on_d, dy_on_d, hr_on_d, mn_on_d, sec_on_d
      print *, yr_off_d, mon_off_d, dy_off_d, hr_off_d, mn_off_d, sec_off_d
     
      print *, yr_on_p, mon_on_p, dy_on_p, hr_on_p, mn_on_p, sec_on_p
      print *, yr_off_p, mon_off_p, dy_off_p, hr_off_p, mn_off_p, sec_off_p
      
      print *, 'Writing ascii file: ', trim(outfile_ascii)
      open (14,file=trim(outfile_ascii),status='new')
      do i = 1,nlat
        do j = 1,nlon
            write(14,900) latvec(i), lonvec(j), grid_mask(i,j), &
            grid_clean(i,j), grid_stacount(i,j), grid_azigap(i,j)
        enddo
      enddo 
      close(14)



end program readwrite_imsvasscbin_gridproc

