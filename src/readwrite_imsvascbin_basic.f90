      program readwrite_imsvasscbin_basic

      ! Utility program to read binary output files produced by
      ! IMS_VASC and produce an ASCII file.
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
      integer npixcross
      integer yr_on, mon_on, dy_on, hr_on, mn_on
      real sec_on
      integer yr_off, mon_off, dy_off, hr_off, mn_off
      real sec_off
      real*8 tday_son, tday_soff
      integer i, j
      
700   format(f8.3,1x,f10.3,1x,i10,1x,i10,1x,f10.3) ! combined xyz file (ascii)


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
      
      print *, "Reading arrays"
      read(13) latvec(1:nlat), lonvec(1:nlon)
      read(13) grid_stack(1:nlat,1:nlon)
      read(13) grid_stacount(1:nlat,1:nlon)
      read(13) grid_azigap(1:nlat,1:nlon)
      read(13) yr_on, mon_on, dy_on, hr_on, mn_on, sec_on, tday_son
      read(13) yr_off, mon_off, dy_off, hr_off, mn_off, sec_off, tday_soff
          
      close(13)
      
      print *, yr_on, mon_on, dy_on, hr_on, mn_on, sec_on
      print *, yr_off, mon_off, dy_off, hr_off, mn_off, sec_off
      print *, tday_son, tday_soff, tday_soff-tday_son ! note that we can use this for grid processing
     
      print *, 'Writing ascii file: ', trim(outfile_ascii)
      open (14,file=trim(outfile_ascii),status='new')
      do i = 1,nlat
        do j = 1,nlon
            write(14,700) latvec(i), lonvec(j), grid_stack(i,j), grid_stacount(i,j), grid_azigap(i,j)
        enddo
      enddo 
      close(14)


end program readwrite_imsvasscbin_basic

