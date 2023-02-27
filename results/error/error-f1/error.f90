program error
implicit none

integer, parameter :: dp = selected_real_kind(15, 307)
!PARAMETERS
!Directories' adresses
character(*), parameter :: dir_m = '../../main/main-f1/'
character(*), parameter :: dir_mo = '../../main-opt/main-opt-f1/'
integer, parameter :: N = 5000
!Neq is for when embedded and scattered electrons are equal
integer :: Ne_m, Ns_m, Neq_m
integer :: Ne_mo, Ns_mo, Neq_mo
integer :: i, aux_int, N_comp! N_comp number up until which comparison is valid
real(dp), dimension(2) :: aux_dp
real(dp), dimension(3) :: r_fep, r_emb, r_sct

!TIME DISTRIBUTION**************************************************************
open(unit=11,file=dir_m//'td.dat',status='unknown')
open(unit=12,file=dir_mo//'td.dat',status='unknown')
open(unit=13,file='td-info.dat',status='unknown')

	do i=1, N
		read(11,*) Ne_m, Ns_m, aux_int, aux_dp
		if (Ne_m .eq. Ns_m) Neq_m = aux_int
		read(12,*) Ne_mo, Ns_mo, aux_int, aux_dp
		if (Ne_mo .eq. Ns_mo) Neq_mo = aux_int
		if ((Ne_m .eq. Ne_mo) .and. (Ns_m .eq. Ns_mo)) N_comp = aux_int
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
	print*, 'Comparison is valid until', N_comp
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

!CHARGE PATCH*******************************************************************
!open(unit=11,file=dir_m//'cp.dat',status='unknown')
!open(unit=12,file=dir_mo//'cp.dat',status='unknown')
!open(unit=13,file='cp-error.dat',status='unknown')

!	do i=1, Ne
	
!	end do

!close(11)
!close(12)
!close(13)

!SCATTERED ANGLES***************************************************************

!open(unit=13,file='sea.dat',status='unknown')

!close(13)

!open(unit=14,file='td.dat',status='unknown')

!close(14)

!	i = 0
!	Ne = 0
!	Ns = 0
!	read(12,*) r_emb, aux
!	read(13,*) r_sct, aux
!	do k=1, 5000
!		read(11,*) r_fep, i		
!		if (all(r_fep .eq. r_emb)) then
!			Ne = Ne + 1
!			if (k .lt. 5000) read(12,*,IOSTAT=IOstatus) r_emb, aux
!		else
!			Ns = Ns + 1
!			if (k .lt. 5000) read(13,*,IOSTAT=IOstatus) r_sct, aux
!		end if
!		print*, Ne, Ns, i, r_emb, r_sct
!		write(14,*) Ne, Ns, i, real(Ne,dp)/i, real(Ns,dp)/i
!	end do
!close(11)
!close(12)
!close(13)
!close(14)

end program