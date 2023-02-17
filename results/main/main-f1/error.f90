program error
implicit none

integer, parameter :: dp = selected_real_kind(15, 307)

!integer :: Ns, Ne, i
!integer :: aux, k, iostatus
real(dp), dimension(3) :: r_fep, r_emb, r_sct

open(unit=11,file='error.dat',status='unknown')

open(unit=12,file='emb.dat',status='unknown')

close(12)

open(unit=13,file='sea.dat',status='unknown')

close(13)

open(unit=14,file='td.dat',status='unknown')

close(14)

	i = 0
	Ne = 0
	Ns = 0
	read(12,*) r_emb, aux
	read(13,*) r_sct, aux
	do k=1, 5000
		read(11,*) r_fep, i		
		if (all(r_fep .eq. r_emb)) then
			Ne = Ne + 1
			if (k .lt. 5000) read(12,*,IOSTAT=IOstatus) r_emb, aux
		else
			Ns = Ns + 1
			if (k .lt. 5000) read(13,*,IOSTAT=IOstatus) r_sct, aux
		end if
		print*, Ne, Ns, i, r_emb, r_sct
		write(14,*) Ne, Ns, i, real(Ne,dp)/i, real(Ns,dp)/i
	end do
close(11)
close(12)
close(13)
close(14)

end program