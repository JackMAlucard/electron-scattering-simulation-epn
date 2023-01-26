program main

use m0!kind_parameters, only: sp, dp, i8
!use units, only: EinAU, dinAUfromSI
use m1!MS1, only: electron_beam, trnsf_to_surf_ref_sys
use m2!MS2! only???
use m3!MS3

implicit none
!*******************************************************************************
!MS1 Test variables
!electron_beam subroutine variables
real(dp), allocatable :: r(:,:), v(:,:), a(:,:) !Electron arrays
integer(i8) :: N_eb	!Number of electrons in the beam
real(dp) :: E, R_eb	!Energy (keV) and radius (Å) of the beam
integer :: pd	!Probability distribution
!trnsf_to_surf_ref_sys
real(dp) :: theta, d0	!incident angle (º) and initial electron beam distance (Å)
!*******************************************************************************
!MS2 Test variables
!simple_silica_model variables
real(dp) :: d!, Lx, Ly, Lz?	!Interatomic distance and dimensions of the cuboid
integer(i8) :: Nx, Ny, Nz!, N_scc
real(dp), allocatable :: ssm(:,:,:,:)	!Indexes and positions array
integer, allocatable :: Z(:,:,:) !Indexes and charges array
!*******************************************************************************
!MS3 Test variables
!sort_electron and add_electron_to_bin variables
integer(i8) :: Nxb, Nyb, Nzb, N_bins	!Bins boundary indices and number of bins
integer(i8) :: max_eq		!Maximum possible equivalent charge, 180 from tests (max 182)
!Maximum possible number of equivalent electrons per bin
integer(i8) :: max_neq	!At most int(N_beam/max_eq)+1 equivalent electrons
!Number of equivalent electrons present in each bin, between 0 and max_neq
integer(i8), allocatable :: b_neq(:,:,:)	!(-Nxb:Nxb,-Nyb:-1,-Nzb:Nzb)
!Equivalent charge of each equivalent electron of each bin, between 0 and max_eq
integer(i8), allocatable :: b_ec(:,:,:,:)	!(-Nxb:Nxb,-Nyb:-1,-Nzb:Nzb,max_neq)
!Equivalent (mean) position of each equivalent electron of each bin
real(dp), allocatable :: b_er(:,:,:,:,:)	!(-Nxb:Nxb,-Nyb:-1,-Nzb:Nzb,max_neq)
real(dp) :: FFd	!Far Field (effective) distance
integer :: bin_size
!trajectory subroutine variables
integer(i8) :: P, T
real(dp) :: dt, a0, v0, tf, l
integer(i8) :: Ne, Ns	!Number of embedded and scattered electrons
real(dp), allocatable :: emb(:,:)	!Simulation embedded electrons array
real(dp), allocatable :: sct(:,:) !Simulation scattered electrons array
!scattered electron angles
real(dp) :: xs, ys, zs
real(dp) :: s, alfa, beta
!*******************************************************************************
!General variables
integer(i8) :: i, j, k
integer(i8), parameter :: ou = 13, oiu = 24	!Output file for data and info
real(sp) :: startT, endT, execTime			!Program timer variables
integer :: seed_size
integer, allocatable :: seed(:)
real(dp) :: u
integer(i8) :: Ne_p, Ns_p	!Count of embedded or scattered electrons to plot
	
	!FIX RANDOM NUMBERS FOR REPRODUCIBILITY***************************************
	call random_seed(size=seed_size)
	allocate(seed(seed_size))
	seed = 23466992! putting arbitrary seed to all elements
	call random_seed(put=seed)
	deallocate(seed)

	!MS1 TEST: ELECTRON BEAM GENERATION*******************************************	
	!Input parameters
	N_eb = 150
	E = 10	!keV
	R_eb = 3*1.77 !Å, considering the width of the standard cuboid in scc is Nx*1.77
	pd = 1 !Gaussian !4 Uniform distribution
	call electron_beam(r, v, a, N_eb, E, R_eb, pd)
	!e_pos, e_vel, e_acc, N_in, E_in, R_in, prob_dist_in
	open(unit=11,file='haz.dat',status='unknown')
		do i=1,size(r)/3
			write(11,*) r(i,:), v(i,:), a(i,:), 0
		end do

		write(11,*)
		write(11,*)
		
		!Input parameters
		theta = 84!º
		d0 = 80!Å
		call trnsf_to_surf_ref_sys(r, v, a, theta_in=theta, d_in=d0)
		!e_pos, e_vel, e_acc, theta_in, d_in
		
		do i=1,size(r)/3
			write(11,*) r(i,:), v(i,:), a(i,:), 0
		end do

	close(11)

	!UNCOMMENT/COMMENT TO SEE/HIDE PLOTS
!	call system('gnuplot "haz.gp"')
!	call system('gnuplot "vel.gp"')


	!MS2 TEST: ATOMIC POSITIONS GENERATION****************************************

	!FOR SIMPLE SILICA MODEL
	!Input parameters
	Nx = 30
	Ny = 30
	Nz = 75
	d = 1.77!Å, sum of Silicon and Oxygen covalent radii
	call simple_silica_model(ssm, Z, Nx, Ny, Nz, d)
	
	open(unit=11,file='ssm_Si.dat',status='unknown')
	open(unit=12,file='ssm_O.dat',status='unknown')
	
	do i = -Nx, Nx
		do j = 0, -Ny, -1
			do k = -Nz, Nz
				if (Z(i,j,k) .eq. 14) write(11, *) ssm(i,j,k,:)
				if (Z(i,j,k) .eq. 8) write(12, *) ssm(i,j,k,:)
			end do
		end do
	end do

	close(11)
	close(12)

	!Print to console charge positions
!	do j = 0, -Ny, -1
!		!print*, "LAYER:", j
!		do i = -Nx, Nx
!			!print*, Z(i,j,:)
!			!print*
!		end do
!	end do
	!read*

	!PLOT THAT ONLY SHOWS SSM STRUCTURE 
!	call system('gnuplot "ssm.gp"')
!	print*, "Stopping point, Ctrl+C to end program"
!	read*
	
	
	!EMBEDDED AND SCATTERED ELECTRONS' POSITIONS ARRAYS
	!Allocating and initializing embedded and scattered electron positions arrays
	allocate(emb(N_eb,3))
	allocate(sct(N_eb,3))
	emb = 0
	sct = 0


	!MS3 TEST: TRAJECTORIES*******************************************************
	print*, "Interatomic distance (a0):", d
!	print*, "Bin distance (a0)", d_bin, d_bin_inv
	print*, "Boundaries"
	print*, "Lx: +/-", Nx*d
	print*, "Ly:    ", -Ny*d
	print*, "Lz: +/-", Nz*d
	
	!Trajectory subroutine input parameters
	P = 5000	!Number of points to be plotted
	dt = 1.d-5	!Time step size
	l = 1/33.7465!SiO2 mean free using weighted scattering cross sections
	print*, "l", l
	print*, "Mean free path:", 1/l
	Ne = 0	!Initially, zero electrons are embedded
	Ns = 0	!Initially, zero electrons are scattered
	Ne_p = 0
	Ns_p = 0	
	!Open file to write electron trajectories
	open(unit=ou, file="beam_trajectories.dat", status='unknown')
	!Open file to write final position of embedded electrons
	open(unit=ou+1,file='emb.dat',status='unknown')
	!Open file to write final position of scattered electrons
	open(unit=ou+2,file='sct.dat',status='unknown')
	!Open file to write final position of electrons in general
	open(unit=ou+3,file='fep.dat',status='unknown')
	!Open file to write scattered electrons angles
	open(unit=ou+4,file='sea.dat',status='unknown')
	!Open file to write time distribution of # of embedded and scattered electrons
	open(unit=ou+5,file='td.dat',status='unknown')
	!Simulation time start
	call cpu_time(startT)
		!Computing the trajectories of each electron in the beam
		do i=1, N_eb
			!Final time estimation
			d0 = norm2(r(i,:))
			a0 = 1._dp/(d0**2)	!Initial acceleration due to N_e electrons, in au
			v0 = dsqrt(2*E)	!Initial speed taken from E, main energy of the beam, in au
			tf = (-v0 + dsqrt(v0*v0 - 2*a0*d0))/a0
			tf = -2*tf
			!Estimated number of iterations
			T = dint(tf/dt)

			!Printing simulation info on console (1/2)
			!print*, "Estimated number of iterations T:", T
			!print*, "Number of points to be plotted P:", P
			if (T .lt. P) P = T !Can't plot less points than the number of simulated ones
			!print*, "Final estimated number of iterations T:", T
			print*, "Simulating", i, "out of", N_eb, "electron trajectories"

			!Computing i-th trajectory
			call trajectory(i,N_eb,P,ou,Ne,Ns,emb,sct,Nx,Ny,Nz,ssm,Z,d0,T,d,dt,l,r,v,a)
			write(ou,*)
			write(ou,*)
			!Write final electron positions to file 
			if (Ne .gt. Ne_p) then
				Ne_p = Ne
				write(ou+1, *) emb(Ne,:), i
			else if (Ns .gt. Ns_p) then
				Ns_p = Ns
				write(ou+2, *) sct(Ns,:), i
				xs = sct(Ns,1)
				ys = sct(Ns,2)
				zs = sct(Ns,3)
				s = dsqrt(xs**2 + zs**2)
				alfa = datan2(ys,s)*180/PI
				beta = -datan2(xs,zs)*180/PI
				write(ou+4, *) alfa, beta, i
			end if
			write(ou+3, *) r(i,:), i
			!Write time distribution of number of embedded and scattered electrons
			write(ou+5,*) Ne, Ns, i, real(Ne,dp)/i, real(Ns,dp)/i

			!Simulation time from beginning to this electron
			call cpu_time(endT)
			execTime = endT - startT
			!Printing simulation info on console (2/2)
			print*, "  Number of electrons embedded: ", Ne
			print*, "  Number of electrons scattered:", Ns
			print*, "Simulation time thus far:", execTime
			!Writing simulation to file
			open(unit=oiu,file='sim-info.dat',status='unknown')
				write(oiu,*) i, "out of", N_eb, "electron trajectories simulated"
				write(oiu,*) "  Number of electrons embedded: ", Ne
				write(oiu,*) "  Number of electrons scattered:", Ns
				write(oiu,*) "Simulation time thus far:", execTime
			close(oiu)
		end do
	close(ou)
	close(ou+1)
	close(ou+2)
	close(ou+3)
	close(ou+4)
	close(ou+5)
	
	call cpu_time(endT)
	execTime = endT - startT
	!SIMULATION END MESSAGES
	!print*
	print*, "SIMULATION COMPLETE!"
	print*, "Total simulation time:", execTime
	open(unit=oiu,file='sim-info.dat',status='unknown')
		write(oiu,*) "SIMULATION COMPLETE!"
		write(oiu,*) "  Number of electrons embedded:", Ne
		write(oiu,*) "  Number of electrons scattered:", Ns
		write(oiu,*) "Total simulation time:", execTime
	close(oiu)
	
	!WRITE POSITION OF SIMULATED EMBEDDED ELECTRONS TO FILE
!	open(unit=11,file='emb.dat',status='unknown')
!	do i = 1, Ne
!		write(11, *) emb(i,:)
!	end do
!	close(11)

	!WRITE POSITION OF SIMULATED SCATTERED ELECTRONS TO FILE
!	open(unit=11,file='sct.dat',status='unknown')
!	do i = 1, Ns
!		write(11, *) sct(i,:)
!	end do
!	close(11)

	!PLOT THAT SHOWS SIMULATED EMBEDDED, SCATTERED ELECTRONS AND SCC STRUCTURE 	
	!call system('gnuplot "particles.gp"')
	
	!PLOT THAT SHOWS ELECTRON TRAJECTORIES 	
	!call system('gnuplot "beam_trajectories_zx.gp"')
	!call system('gnuplot "beam_trajectories_zy.gp"')


	if (allocated(r)) deallocate(r)
	if (allocated(v)) deallocate(v)
	if (allocated(a)) deallocate(a)
	if (allocated(ssm)) deallocate(ssm)
	if (allocated(Z)) deallocate(Z)
	if (allocated(emb)) deallocate(emb)
	if (allocated(sct)) deallocate(sct)
	if (allocated(b_neq)) deallocate(b_neq)
	if (allocated(b_ec)) deallocate(b_ec)
	if (allocated(b_er)) deallocate(b_er)

end program