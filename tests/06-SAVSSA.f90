!Scattering Angle vs Step Size Automatic (SAVSSA)
!This program aims to obtain data relating step size and change (or percent
!change) in scattering angle.
program SAVSS
use kind_parameters, only: sp, dp, i8
use units!, only: !EinAU, dinAU
use, intrinsic:: iso_fortran_env, only: stdin=>input_unit
implicit none
!General
real(dp), parameter :: PI = dacos(-1.d0)
real(dp), parameter :: t0 = 0._dp
integer(i8) :: i, j, option	!Program structure parameters
integer(i8) :: N 						!Number of points to be plotted/simulated
real(dp) :: K0, K0f, hd, hdf	!Fixed simulation parameters
real(dp) :: b, dt							!Variable simulation parameters
real(dp) :: rt(3)						!Target electron position
real(dp) :: r0(3), v0(3)		!Simulation position and velocity
integer(i8) :: T						!Number of iterations
real(dp) :: a, c, e!, b			!Theoretical trajectory geometric parameters
real(dp) :: xf, yf		!Last theoretical trajectory point generated
real(dp) :: TSA, PSA, SSA	!Theoretical, Plot, and Simulation Scattering Angles
real(dp) :: ti, ri(3)	!Time and space simulation variables
real(dp) :: vi(3)			!Velocity simulation variables
real(dp) :: ai(3) 		!Acceleration simulation variables
real(dp) :: U0, L0, Ui, Li	!Conserved quantities variables
logical :: estimated_time
character(len=24) :: current_time
real(sp) :: startT, endT, execTime			!Program timer variables
logical :: approaching_center		!Approaching Center
real(dp) :: cdsc, dsc, crsc(3) 	!Closest Distance to Scattering Center
integer(i8) :: NS, k
character(len=*), parameter :: input_file = "06-SAVSSA.in", aux_file = "aux.in"
character(len=10) :: K0c, hdc, bc, dtc!AT MOST 10 CHARACTERS FOR INPUT VALUES!!!
character(len=:), allocatable :: K0ct, hdct, bct, dtct
character(len=:), allocatable :: output_file, output_file_info
integer(i8), parameter :: inu = 11, icu = 12	!Input file as numbers and chars
integer(i8), parameter :: ou = 13, oiu = 14	!Output file for data and info
character(len=80) :: FMTS		!Format string

	print*, "Scattering Angle vs Step Size Automatic (SAVSSA)"
	print*, "This program originally aimed to obtain data relating step size and"
	print*, "change (or percent change) in scattering angle."
	print*, ""
	print*, "Now it is used to test when the VV algorithm breaks by varying"
	print*, "mainly the impact parameter and the time step size used."
	print*, ""
!*******************************************************************************
	!Reading fixed parameters from input file
	!Open input file twice: unit 11 to use as numbers, unit 12 to use as chars
	call system ("cp "//input_file//" "//aux_file)
	open(unit=inu, file=input_file, status='unknown')
	!open(unit=icu, file=input_file, status='unknown')
	open(unit=icu, file=aux_file, status='unknown')

	!Skip the first 5 lines of the input file
	do i=1,5
		read(inu, *)
		read(icu, *)
	end do

	!# of simulations/parameter sets, NS
	read(inu, *) NS
	read(icu, *)
	!# of points to be plotted, N
	read(inu, *) N
	read(icu, *)
	!Initial kinetic energy, K0 [keV]
	read(inu, *) K0f
	read(icu, *) K0c
	K0ct = trim(K0c)
	!Initial horizontal distance, hd [Å]
	read(inu, *) hdf
	read(icu, *) hdc
	hdct = trim(hdc)

	!Skip the next 4 lines of the input file
	do i=1,4
		read(inu, *)
		read(icu, *)
	end do

!*******************************************************************************
	do k=1, NS
		!Fixed parameter values
		K0 = K0f
		hd = hdf
		!Reading variable parameters from input file
		!Impact parameter, b [Å]
		read(inu, *) b
		read(icu, *) bc
		bct = trim(bc)
		!Time step size, dt [aut]
		read(inu, *) dt
		read(icu, *) dtc
		dtct = trim(dtc)

		!Naming and opening output files
		output_file = K0ct//'_'//hdct//'_'//bct//'_'//dtct//'.dat'
		output_file_info = K0ct//'_'//hdct//'_'//bct//'_'//dtct//'_info.dat'
		open(unit=ou, file=output_file, status='unknown')
		open(unit=oiu, file=output_file_info, status='unknown')

		!Simulation info to print in console
		write (oiu, "('*** ELECTRON-ELECTRON SCATTERING ***')")
		print "('SIMULATION ', i3, ' OUT OF ', i3)", k, NS
		print*, 'K0: '//K0ct//'[keV], hd: '//hdct//'[Å], b: '//bct//'[Å], dt: '//dtct//'[aut]'
		!Unit conversion of simulation parameters to au and printing to info file
		call parameter_init(hd, b, K0, dt, r0, v0, T, N, oiu)
		!Add blank space on Output Info File
		write (oiu, *)

		!***************************************************************************
		!Theoretical trajectory
		!Plotting N points
		call theoretical_trajectory(N, r0, v0, a, b, c, e, xf, yf, ou, oiu)
		!Add blank space on Output Info File
		write (oiu, *)

		!***************************************************************************
		!Simulation initialization
		rt = 0

		ti = t0
		ri = r0
		vi = v0
		call akP(rt, ri, ai)

		!Initial values of conserved quantities
		call conserved_quantities(r0, v0, U0, L0)

		Ui = U0
		Li = L0

		write(ou,*) ti, ri, Ui, 100*(dabs(Ui-U0)/U0), Li, 100*(dabs(Li-L0)/L0)
		j = 1

		!Always checking for the Closest Distance to the Scattering Center
		approaching_center = .true.
		cdsc = norm2(ri)
		crsc = ri

		estimated_time = .true.
		!Start timer RIGHT before the first iteration
		call cpu_time(startT)

		do i=1,T
			!Plotting only N points
			if ( (mod(i,T/N) .eq. 0) .and. j .lt. N ) then
				!Computing conserved quantities
				call conserved_quantities(ri, vi, Ui, Li)

				!Write values to file
				write(ou,*) ti, ri, Ui, 100*(dabs(Ui-U0)/U0), Li, 100*(dabs(Li-L0)/L0)
				j = j + 1

				!For estimation of simulation time (it runs only once)
				if (estimated_time) then
					!Compute iteration time after 1/N-enth of the simulation
					call cpu_time(endT)
					execTime = endT - startT
					call cpu_time(startT)
					!Estimate simulation time, print to console and output file
					call fdate(current_time)
					FMTS = "('Estimated simulation time: ', f8.2, '[s]', f8.2, '[min]')"
					!Print to console
					print "('Start time:                ', a)", current_time
					print FMTS, execTime*N, execTime*N/60
					!Write to Output Info File
					write(oiu, "('*** SIMULATION TIME ***')")
					write(oiu, "('Start time:                ', a)") current_time
					write(oiu, FMTS) execTime*N, execTime*N/60
					!To only run this statement once
					estimated_time = .false.
				end if

			end if

			!Velocity Verlet step calculation
			call vv_step(i, rt, t0, dt, ti, ri, vi, ai)

			!Checking if still approching center, searching closest distance to center
			if (approaching_center) then
				dsc = norm2(ri)
				if (dsc .gt. cdsc) then
					approaching_center = .false.
				else
					cdsc = dsc
					crsc = ri
				end if
			end if

		end do

		call cpu_time(endT)
		execTime = execTime + (endT - startT)

		call fdate(current_time)
		FMTS = "('Total simulation time:     ', f8.2, '[s]', f8.2, '[min]')"
		!Print to console
		print "('End time:                  ', a)", current_time
		print FMTS, execTime, execTime/60
		!Write to Output Info File
		write(oiu, "('End time:                  ', a)") current_time
		write(oiu, FMTS) execTime, execTime/60

		!Add blank space on Output Info File
		write (oiu, *)

		!***************************************************************************
		!Scattering Angles computation and comparison
		write(oiu, "('*** SCATTERING ANGLE COMPARISONS ***')")

		!Theoretical Scattering Angle
		TSA = 2*datan(1/(2*K0*b))
		FMTS = "('Theoretical Scattering Angle (TSA): ', f20.16, 'º')"
		write(oiu, FMTS) TSA*180/PI
		!Plot Scattering Angle
		PSA = datan2(yf,xf)
		if (yf .lt. 0._dp) PSA = 2*PI + PSA
		FMTS = "('Plot Scattering Angle (PSA):        ', f20.16, 'º')"
		write(oiu, FMTS) PSA*180/PI
		!Simulation Scattering Angle
		SSA = datan2(ri(2),ri(1))
		if (ri(2) .lt. 0._dp) SSA = 2*PI + SSA
		FMTS = "('Simulation Scattering Angle (SSA):  ', f20.16, 'º')"
		write(oiu, FMTS) SSA*180/PI

		!Comparison between SSA and PSA
		FMTS = "('Percent error between SSA and PSA: ', d12.4, '%')"
		write(oiu, FMTS) 100*dabs(SSA-PSA)/PSA
		!Comparison between SSA and TSA
		FMTS = "('Percent error between SSA and TSA: ', d12.4, '%')"
		write(oiu, FMTS) 100*dabs(SSA-TSA)/TSA
		!Comparison between PSA and TSA
		FMTS = "('Percent error between PSA and TSA: ', d12.4, '%')"
		write(oiu, FMTS) 100*dabs(PSA-TSA)/TSA

		!Add blank space on Output Info File
		write (oiu, *)

		!***************************************************************************
		!Closest distance and position to scattering center

		write(oiu, "('*** CLOSEST DISTANCE AND POSITION TO SCATTERING CENTER ***')")
		FMTS = "('Closest distance to scattering center:   ', d12.4, '[au]')"
		write(oiu, FMTS) cdsc
		FMTS = "('Closest postion to scattering center, x: ', d12.4, '[au]')"
		write(oiu, FMTS) crsc(1)
		FMTS = "('Closest postion to scattering center, y: ', d12.4, '[au]')"
		write(oiu, FMTS) crsc(2)
		FMTS = "('Closest postion to scattering center, z: ', d12.4, '[au]')"
		write(oiu, FMTS) crsc(3)

		close(ou)
		close(oiu)

		print*

	end do

	close(inu)
	close(icu)
	call system ("rm "//aux_file)

contains

!This subroutine reads values for initial horizontal distance, hd (Å), impact
!parameter, b (Å), initial kinetic energy, K0 (KeV), and time step size,
!dt (aut). It uses these values to initialize the simulation.
!It also uses the number of points to be plotted, N, and the output info file
!unit, oifu.
!Here, the final simulation time is estimated assuming that acceleration is
!constant.
subroutine parameter_init(hd, b, K0, dt, r0, v0, T, N, oifu)
implicit none
integer(i8), intent(in) :: oifu
integer(i8), intent(inout) :: N
real(dp), intent(inout) :: K0, hd, b, dt			!Variable simulation parameters
real(dp), intent(out) :: r0(3), v0(3)	!Simulation variables
integer(i8), intent (out) :: T						!Number of iterations
real(dp) :: conv_aux				!Conversion auxiliar to print values
real(dp) :: v, a, tf				!Acceleration for time estimation
character(len=80) :: FMTS		!Format string

	write(oifu, "('*** SIMULATION PARAMETERS ***')")
	write(oifu, "('Number of points to be plotted, N: ', i6)") N

	!Initial kinetic energy, K0
	conv_aux = K0
	K0 = EinAU(K0)
	FMTS = "('Initial kinetic energy, K0:     ', d12.4, '[keV] =', d12.4, '[Eh]')"
	write(oifu, FMTS) conv_aux, K0

	!Initial horizontal distance, hd
	conv_aux = hd
	hd = dinAU(hd)
	FMTS = "('Initial horizontal distance, hd:', d12.4, '[Å]   =', d12.4, '[a0]')"
	write(oifu, FMTS) conv_aux, hd

	!Impact parameter, b
	conv_aux = b
	b = dinAU(b)
	FMTS = "('Impact parameter, b:            ', d12.4, '[Å]   =', d12.4, '[a0]')"
	write(oifu, FMTS) conv_aux, b

	!Time step size, dt
	conv_aux = dt*2.418d-17
	FMTS = "('Time step size, dt:             ', d12.4, '[aut] =', d12.4, '[s]')"
	write(oifu, FMTS) dt, conv_aux

	!Parameter initialization
	!Initial position
	r0 = (/-hd, b, 0._dp/)

	!Initial velocity
	v0 = (/dsqrt(2*K0), 0._dp, 0._dp/)

	!Final time estimation
	a = norm2(r0)
	a = 1/a
	v = norm2(v0)
	tf = (-v + dsqrt(v*v - 2*a*hd))/a
	tf = -2*tf
	!Number of iterations
	T = dint(tf/dt)
	if (T .lt. N) N = T !Can't plot less points than the number of simulated ones
end subroutine parameter_init

!Number of points to be plotted, initial position vector, initial velocity
!vector, hyperbola geometrical parameters.
!The resulting coordinates are written into the Output file Unit, ou, and
!the trajectory parameters to Output Info file Unit, oiu which must be
!previously opened in the program.
subroutine theoretical_trajectory(N, r0, v0, a, b, c, e, xf, yf, ou, oiu)
implicit none
real(dp), parameter :: PI = dacos(-1.d0)
integer(i8), intent(in) :: N						!Number of points to be plotted
integer(i8), intent(in) :: ou, oiu			!Output files Unit's
real(dp), intent(in) :: r0(3), v0(3)		!Required simulation parameters
real(dp), intent(out) :: a, b, c, e			!Hyperbola geometrical parameters
real(dp), intent(out) :: xf, yf					!Last point coordinates
real(dp) :: x0, y0, vx0
real(dp) :: phi0, phif, dphi, phii, ri	!Plotting variables
real(dp) :: xi, yi, alpha, den, num
	!From input vectors to dummy variables
	x0 = r0(1)
	y0 = r0(2)
	vx0 = v0(1)

	a = 1/(vx0*vx0) !a = 1/(2*K0), K0 = 0.5*v0**2
	b = y0
	c = dsqrt(a*a + b*b)
	e = c/a

	write(oiu, "('*** THEORETICAL TRAJECTORY PARAMETERS ***')")
	write(oiu, "('Hyperbola geometric parameters')")
	write(oiu, "('a:', e12.4, '[au]')") a
	write(oiu, "('b:', e12.4, '[au]')") b
	write(oiu, "('c:', e12.4, '[au]')") c
	write(oiu, "('e:', e12.4, '[au]')") e

	!*****************************************************************************
	!Plotting the scattering trajectory, i.e. the Hyperbola's Left Branch
	!using the simulation initialization parameters

	!Compute the angle alpha. which limits the range for the left branch
	alpha = dacos(1/e)

	!Compute the angle coordinate of the initial simulation point.
	phif = datan2(y0,x0)
	if (y0 .lt. 0._dp) phif = phif + 2*PI

	!Substract PI from the resulting angle to get the angle correspoding
	!to the LAST POINT to be plotted using the polar coordinates hyperbola
	!equation. This is necessary because the way in which the left branch
	!is generated via this equation uses negative r values, i.e. -alpha<phi<alpha
	phif = phif - PI
	phi0 = -phif - 2*alpha
	dphi = (phif-phi0)/N

	!Cartesian coordinates trajectory equations:
			!ri = (b*b/a)/(1 - e*dcos(phii + alpha))
			!ri = (b*b/c)/(1/e - dcos(phii + alpha))
			!ri = (b*b)/(a - c*dcos(phii + alpha))
			!There is loss of significant figures on: a - c*dcos(phii + alpha)
			!for really small values of impact parameter, b
			!Using difference of squares trick to somewhat lessen the error
			!x - y = (x**2 - y**2)/(x + y)
	do i=0, N
		phii = phi0 + i*dphi
		num = a*a - c*c*dcos(phii + alpha)*dcos(phii + alpha)
		den = a + c*dcos(phii + alpha)
		ri = (b*b*den)/num
		xi = ri*dcos(phii)
		yi = ri*dsin(phii)
		write (ou, *) xi, yi
		if (i .eq. 0) then
			!Save coordinates for "final" point of the theoretical trajectory
			xf = xi
			yf = yi
		end if
	end do

	write (ou, *)
	write (ou, *)
end subroutine theoretical_trajectory

subroutine conserved_quantities(r, v, U, L)!(x, y, z, vx, vy, vz, L, K, V, U)
real(dp), intent(in) :: r(3), v(3)!x, y, z, vx, vy, vz
real(dp), intent(out) :: U, L
real(dp) :: Uk, Ue, Lx, Ly, Lz

	Uk = 0.5*(v(1)**2 + v(2)**2 + v(3)**2)!(vx*vx + vy*vy + vz*vz)
	Ue = 1/norm2(r)!(dsqrt(x*x + y*y + z*z))
	U = Uk + Ue

	Lx = r(2)*v(3) - r(3)*v(2)!  y*vz - z*vy
	Ly = r(3)*v(1) - r(1)*v(3)!  z*vx - x*vz
	Lz = r(1)*v(2) - r(2)*v(1)!  x*vy - y*vx
	L = dsqrt(Lx**2 + Ly**2 + Lz**2)
end subroutine conserved_quantities

!Acceleration for interaction between point particles for Velocity Verlet
!Algorithm using atomic units (au)
subroutine akP(rt, ri, a)
!rt: r target, position of scattering center, i.e. stationary particle
!ri: r incoming, position of moving/incident particle
!a: Acceleration of moving/incident particle
implicit none
real(dp), intent(in) :: rt(3), ri(3)
real(dp), intent(out) :: a(3)
real(dp) :: rs(3), r !rs: r separation vector and its magnitude
	rs = ri - rt
	r = norm2(rs)
	a = rs/(r**3)
end subroutine akP

subroutine vv_step(i, rt, t0, dt, t, r, v, a)
integer(i8), intent(in) :: i
real(dp), intent(in) :: rt(3)
real(dp), intent(in) :: t0, dt
real(dp), intent(inout) :: t, r(3), v(3), a(3)
	!Velocity Verlet algorithm
	t = t0 + i*dt
	r = r + v*dt + 0.5*a*dt*dt
	v = v + 0.5*a*dt
	call akP(rt, r, a)
	v = v + 0.5*a*dt
end subroutine vv_step

end program