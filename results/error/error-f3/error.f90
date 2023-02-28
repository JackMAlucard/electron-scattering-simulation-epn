program error
implicit none

integer, parameter :: dp = selected_real_kind(15, 307)
!PARAMETERS
!Directories' adresses
character(*), parameter :: dir_m = '../../main/main-f3/'
character(*), parameter :: dir_mo = '../../main-opt/main-opt-f3/'
integer, parameter :: N = 5000
!Neq is for when embedded and scattered electrons are equal
integer :: Ne_m, Ns_m, Neq_m
integer :: Ne_mo, Ns_mo, Neq_mo
integer :: i, aux_int!, N_comp! N_comp number up until which comparison is valid
real(dp), dimension(2) :: aux_dp
real(dp), dimension(3) :: r_fep, r_emb, r_sct
real(dp), dimension(2) :: r_sea, r_sea_mean_m, r_sea_mean_mo
real(dp) :: r_sea_sd_m, r_sea_sd_mo
real(dp) :: d_sea_mean_m, d_sea_mean_mo
real(dp) :: d_sea_sd_m, d_sea_sd_mo
real(dp), dimension(3) :: r_emb_mean_m, r_emb_mean_mo

!TIME DISTRIBUTION**************************************************************
open(unit=11,file=dir_m//'td.dat',status='unknown')
open(unit=12,file=dir_mo//'td.dat',status='unknown')
open(unit=13,file='td-info.dat',status='unknown')

	do i=1, N
		read(11,*) Ne_m, Ns_m, aux_int, aux_dp
		if (Ne_m .eq. Ns_m) Neq_m = aux_int
		read(12,*) Ne_mo, Ns_mo, aux_int, aux_dp
		if (Ne_mo .eq. Ns_mo) Neq_mo = aux_int
!		if ((Ne_m .eq. Ne_mo) .and. (Ns_m .eq. Ns_mo)) N_comp = aux_int
	end do
	!TO CONSOLE
	print*, 'Regular'
	print*, 'Ne', Ne_m
	print*, 'Ns', Ns_m
	print*, 'Neq', Neq_m
	print*
	print*, 'Optimized'
	print*, 'Ne', Ne_mo
	print*, 'Ns', Ns_mo
	print*, 'Neq', Neq_mo
	print*
	print*, 'Error'
	print*, 'Ne', (abs(Ne_mo-Ne_m)/dfloat(Ne_m))*100._dp, '%'
	print*, 'Ns', (abs(Ns_mo-Ns_m)/dfloat(Ns_m))*100._dp, '%'
	print*, 'Neq', (abs(Neq_mo-Neq_m)/dfloat(Neq_m))*100._dp, '%'
	print*
!	print*, 'Comparison is valid until', N_comp
	!TO FILE
	write(13,*) '!Number of electrons: embedded, scattered and when both are equal'
	write(13,*) '!Regular Simulation'
	write(13,*) Ne_m
	write(13,*) Ns_m
	write(13,*) Neq_m
	write(13,*) '!Optimized Simulation'
	write(13,*) Ne_mo
	write(13,*) Ns_mo
	write(13,*) Neq_mo
	write(13,*) '!Percent error of each number [%]'
	write(13,*) (abs(Ne_mo-Ne_m)/dfloat(Ne_m))*100._dp
	write(13,*) (abs(Ns_mo-Ns_m)/dfloat(Ns_m))*100._dp
	write(13,*) (abs(Neq_mo-Neq_m)/dfloat(Neq_m))*100._dp

close(11)
close(12)
close(13)
!SCATTERED ANGLES***************************************************************
open(unit=11,file=dir_m//'sea.dat',status='unknown')
open(unit=12,file=dir_mo//'sea.dat',status='unknown')

	!Mean (2D)
	!Regular simulation data
	r_sea_mean_m = 0
	do i=1, Ns_m
		read(11,*) r_sea, aux_int
		r_sea_mean_m = r_sea_mean_m + r_sea
	end do
	r_sea_mean_m = r_sea_mean_m/Ns_m
	
	!Optimized simulation data
	r_sea_mean_mo = 0
	do i=1, Ns_mo
		read(12,*) r_sea, aux_int
		r_sea_mean_mo = r_sea_mean_mo + r_sea
	end do
	r_sea_mean_mo = r_sea_mean_mo/Ns_mo

close(11)
close(12)
open(unit=11,file=dir_m//'sea.dat',status='unknown')
open(unit=12,file=dir_mo//'sea.dat',status='unknown')

	!Standard deviation (2D dispersion)
	!Regular simulation data
	r_sea_sd_m = 0
	do i=1, Ns_m
		read(11,*) r_sea, aux_int
		r_sea_sd_m = r_sea_sd_m + norm2(r_sea - r_sea_mean_m)
	end do
	r_sea_sd_m = r_sea_sd_m/Ns_m
	
	!Optimized simulation data
	r_sea_sd_mo = 0
	do i=1, Ns_mo
		read(12,*) r_sea, aux_int
		r_sea_sd_mo = r_sea_sd_mo + norm2(r_sea - r_sea_mean_mo)
	end do
	r_sea_sd_mo = r_sea_sd_mo/Ns_mo

close(11)
close(12)
open(unit=11,file=dir_m//'sea.dat',status='unknown')
open(unit=12,file=dir_mo//'sea.dat',status='unknown')

	!Mean (distance)
	!Regular simulation data
	d_sea_mean_m = 0
	do i=1, Ns_m
		read(11,*) r_sea, aux_int
		d_sea_mean_m = d_sea_mean_m + norm2(r_sea)
	end do
	d_sea_mean_m = d_sea_mean_m/Ns_m
	
	!Optimized simulation data
	d_sea_mean_mo = 0
	do i=1, Ns_mo
		read(12,*) r_sea, aux_int
		d_sea_mean_mo = d_sea_mean_mo + norm2(r_sea)
	end do
	d_sea_mean_mo = d_sea_mean_mo/Ns_mo

close(11)
close(12)
open(unit=11,file=dir_m//'sea.dat',status='unknown')
open(unit=12,file=dir_mo//'sea.dat',status='unknown')

	!Standard deviation (distance)
	!Regular simulation data
	d_sea_sd_m = 0
	do i=1, Ns_m
		read(11,*) r_sea, aux_int
		d_sea_sd_m = d_sea_sd_m + (norm2(r_sea) - d_sea_mean_m)**2
	end do
	d_sea_sd_m = d_sea_sd_m/Ns_m
	d_sea_sd_m = dsqrt(d_sea_sd_m)
	
	!Optimized simulation data
	d_sea_sd_mo = 0
	do i=1, Ns_mo
		read(12,*) r_sea, aux_int
		d_sea_sd_mo = d_sea_sd_mo + (norm2(r_sea) - d_sea_mean_mo)**2
	end do
	d_sea_sd_mo = d_sea_sd_mo/Ns_mo
	d_sea_sd_mo = dsqrt(d_sea_sd_mo)

close(11)
close(12)

	!PRINT RESULTS
	!TO CONSOLE
	print*, 'Scattered angles mean (2D)'
	print*, 'Regular', r_sea_mean_m, norm2(r_sea_mean_m)
	print*, 'Optimized', r_sea_mean_mo, norm2(r_sea_mean_mo)
	print*
	print*, 'Error [%]'
!	print*, 'r_sea_mean_m-r_sea_mean_mo', r_sea_mean_m-r_sea_mean_mo
!	print*, 'r_sea_mean_m', r_sea_mean_m
	print*, dabs((r_sea_mean_m-r_sea_mean_mo)/(r_sea_mean_m))*100._dp
	print*, dabs((norm2(r_sea_mean_m)-norm2(r_sea_mean_mo))/(norm2(r_sea_mean_m)))*100._dp
!	(abs(Ne_mo-Ne_m)/dfloat(Ne_m))*100._dp
	print*, 'Error (distance)'
	print*, norm2(r_sea_mean_m-r_sea_mean_mo)
	print*
	print*, 'Standard deviation (2D dispersion)'
	print*, 'Regular', r_sea_sd_m
	print*, 'Optimized', r_sea_sd_mo
	print*
	print*, 'Error [%]'
	print*, (dabs(r_sea_sd_m-r_sea_sd_mo)/(r_sea_sd_m))*100._dp
	print*
	
	print*, 'Mean (distance r)'
	print*, 'Regular', d_sea_mean_m
	print*, 'Optimized', d_sea_mean_mo
	print*
	print*, 'Error [%]'
	print*, (dabs(d_sea_mean_m-d_sea_mean_mo)/(d_sea_mean_m))*100._dp
	print*
	print*, 'Standard deviation (distance r)'
	print*, 'Regular', d_sea_sd_m
	print*, 'Optimized', d_sea_sd_mo
	print*
	print*, 'Error [%]'
	print*, (dabs(d_sea_sd_m-d_sea_sd_mo)/(d_sea_sd_m))*100._dp
	print*
	

open(unit=13,file='sea-error.dat',status='unknown')

close(13)

!CHARGE PATCH*******************************************************************
!open(unit=11,file=dir_m//'cp.dat',status='unknown')
!open(unit=12,file=dir_mo//'cp.dat',status='unknown')
!open(unit=13,file='cp-error.dat',status='unknown')

!	do i=1, Ne
	
!	end do

!close(11)
!close(12)
!close(13)

end program