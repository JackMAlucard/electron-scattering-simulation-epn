!===============================================================================
! Program    : rutherford_scattering_test_simulations
! Purpose    : This program is intended as a proof of concept to become familiar
!              with the Velocity Verlet algorithm, Rutherford Scattering, and
!              how the different geometric and physical parameters used impact
!              the resulting simulations.
! Description:
!   This program executes the following steps:
!     1. Read simulation input parameters from a file.
!     2. Compute the theoretical hyperbolic trajectory.
!     3. Compute the simulated trajectory using the Velocity Verlet algorithm.
!     4. Compute conserved quantities (energy and angular momentum).
!     5. Compute scattering angles and distnce of closest approach.
! Input      :
!   The following input parameters are read from the 'input.txt' file:
!     - Number of simulations and plotting points
!		  - Initial kinetic energy (K0) and horizontal distance (x0)
!     - Pairs of values per simulation: impact parameter (b) and time step (dt)
! Output     :
!   For each simulation,
!     - A '.dat' file with trajectory and conserved magnitudes data,
!     - An accompanying '_info.dat' file with detailed simulation metadata.
!   The '.dat' file contains two blocks of data separated by two blank lines.
!   The first block contains the theoretical trajectory like this:
!     x-coordinate, y-coordinate
!   The second block contains the simulated magnitudes as follows:
!     time, x-coordinate, y-coordinate, z-coordinate, total mechanical energy,
!       total mechanical energy percent deviation fron initial value,
!       angular momentum, angular momentum percent deviation fron initial value
!   The output values are all in atomic units. 
!===============================================================================
program rutherford_scattering_test_simulations
	
	use, intrinsic:: iso_fortran_env, only: stdin=>input_unit
	implicit none
	! Kind type parameters for increased real precision and integer length
	! Single precision reals, 6 digits, range 10**(-37) to 10**(37)-1; 32 bits
	integer, parameter :: sp = selected_real_kind(6, 37)
	! Double precision reals, 15 digits, range 10**(-307) to 10**(307)-1; 64 bits
	integer, parameter :: dp = selected_real_kind(15, 307)
	! Long length for integers, range -2**63 to 2**63-1; 64 bits
	integer, parameter :: i8 = selected_int_kind(18)
	! General constants
	real(dp), parameter :: PI = dacos(-1.d0)
	real(dp), parameter :: t0 = 0._dp
	
	! Input/output files, filenames, and format string
	character(len=*), parameter :: input_file = "input.txt", aux_file = "aux.txt"
	character(len=10) :: K0_char, x0_char, b_char, dt_char	! 10 characters max
	character(len=:), allocatable :: K0_char_trim, x0_char_trim
	character(len=:), allocatable :: b_char_trim, dt_char_trim
	character(len=:), allocatable :: output_file, info_output_file
	character(len=80) :: format_string
	
	! I/O units
	integer(i8), parameter :: input_values_unit = 11, input_chars_unit = 12
	integer(i8), parameter :: output_unit = 13, info_output_unit = 14
		
	! Fixed input parameters
	integer(i8) :: num_simulations, num_plot_ploints
	real(dp) :: K0, K0f, x0, x0f	! Initial kinetic energy and horizontal distance
	
	! Variable input parameters
	real(dp) :: b, dt							! Impact parameter and time step size
	
	! Theoretical trajectory computation
	real(dp) :: xf, yf
	
	! Simulated trajectory computation variables
	integer(i8) :: max_iterations
	real(dp) :: rt(3)						! Target electron position
	real(dp) :: r0(3), v0(3)		! Initial position and velocity
	real(dp) :: ti, ri(3)				! Time and position
	real(dp) :: vi(3)						! Velocity
	real(dp) :: ai(3)						! Acceleration
	
	! Conserved quantities computation
	real(dp) :: U0, L0, Ui, Li
		
	! Scattering angles computation and comparison
	real(dp) :: scat_angle_theoretical, scat_angle_plot, scat_angle_sim
	
	! Simulation timing
	logical :: estimated_time
	character(len=24) :: current_time
	real(sp) :: start_time, end_time, total_time
	
	! Approach to scattering center detection and measurement
	logical :: approaching_center
	real(dp) :: closest_distance_scattering_center, distance_scattering_center
	real(dp) :: closest_position_scattering_center(3)
	
	! Program structure variables
	integer(i8) :: i, j, k

	! Copy input file and open for numeric and char input
	call system ("cp "//input_file//" "//aux_file)
	open(unit=input_values_unit, file=input_file, status='unknown')
	open(unit=input_chars_unit, file=aux_file, status='unknown')

	! Skip the first 5 lines of the input file (header)
	do i = 1, 5
		read(input_values_unit, *)
		read(input_chars_unit, *)
	end do

	! Read input: number of simulations to be executed and points to be plotted
	read(input_values_unit, *) num_simulations
	read(input_chars_unit, *)
	read(input_values_unit, *) num_plot_ploints
	read(input_chars_unit, *)
	
	! Read input: fixed parameters
	read(input_values_unit, *) K0f				! Initial kinetic energy, K0 [keV]
	read(input_chars_unit, *) K0_char
	K0_char_trim = trim(K0_char)
	
	read(input_values_unit, *) x0f				! Initial horizontal distance, x0 [Å]
	read(input_chars_unit, *) x0_char
	x0_char_trim = trim(x0_char)

	! Skip the next 4 lines of the input file (intermediate section)
	do i = 1, 4
		read(input_values_unit, *)
		read(input_chars_unit, *)
	end do

	do k = 1, num_simulations
		! Assign fixed values
		K0 = K0f
		x0 = x0f

		! Read input: variable parameters
		read(input_values_unit, *) b				! Impact parameter, b [Å]
		read(input_chars_unit, *) b_char
		b_char_trim = trim(b_char)
		
		read(input_values_unit, *) dt				! Time step size, dt [aut]
		read(input_chars_unit, *) dt_char
		dt_char_trim = trim(dt_char)

		! Generate output file names
		output_file = 'K0_'//K0_char_trim//'keV_x0_'//x0_char_trim// &
			'Å_b_'//b_char_trim//'Å_dt_'//dt_char_trim//'aut'
		info_output_file = output_file//'_info.dat'
		output_file = output_file//'.dat'

		! Open output files
		open(unit=output_unit, file=output_file, status='unknown')
		open(unit=info_output_unit, file=info_output_file, status='unknown')

		! Simulation info to console and file
		write(info_output_unit, "('*** ELECTRON-ELECTRON &
			&RUTHERFORD SCATTERING ***')")
		write(info_output_unit, *)		! Blank line on Output Info File
		print "('SIMULATION ', i3, ' OUT OF ', i3)", k, num_simulations
		print*, 'K0: '//K0_char_trim//'[keV], x0: '//x0_char_trim// &
			'[Å], b: '//b_char_trim//'[Å], dt: '//dt_char_trim//'[aut]'

		! Initialize parameters and estimate number of time steps
		call parameter_initialization &
			(info_output_unit, num_plot_ploints, K0, x0, b, dt, r0, v0, &
			max_iterations)
		
		write(info_output_unit, *)		! Blank line on Output Info File
		
		! Theoretical trajectory
		call compute_theoretical_trajectory &
			(num_plot_ploints, output_unit, info_output_unit, r0, K0, xf, yf)
		
		write(info_output_unit, *)		! Blank line on Output Info File

		! Initialize simulation
		rt = 0

		ti = t0
		ri = r0
		vi = v0
		call acceleration_due_to_electron(ri, rt, ai)

		! Initialize conserved quantities
		call compute_conserved_quantities(r0, v0, U0, L0)
		Ui = U0
		Li = L0

		write(output_unit,*) ti, ri, Ui, 100*(dabs(Ui-U0)/U0), &
			Li, 100*(dabs(Li-L0)/L0)
		
		j = 1

		! Closest distance to the scattering center check
		approaching_center = .true.
		closest_distance_scattering_center = norm2(ri)
		closest_position_scattering_center = ri

		estimated_time = .true.
		call cpu_time(start_time)		! Start simulation timer

		do i = 1, max_iterations
			! Conditional logic used to plot only num_plot_ploints
			if ((mod(i,max_iterations/num_plot_ploints) .eq. 0) &
				.and. j .lt. num_plot_ploints) then
				! Update conserved quantities
				call compute_conserved_quantities(ri, vi, Ui, Li)

				! Write simulation values to file
				write(output_unit,*) ti, ri, Ui, 100*(dabs(Ui-U0)/U0), &
					Li, 100*(dabs(Li-L0)/L0)
				
				j = j + 1

				! Simulation execution time estimation (it runs only once)
				if (estimated_time) then
					! Compute iteration time after 1/num_plot_ploints of the simulation
					call cpu_time(end_time)
					total_time = end_time - start_time
					call cpu_time(start_time)
					call fdate(current_time)
					format_string = "('Estimated simulation time: ', f8.2, '[s]', &
						&f8.2, '[min]')"
					
					print "('Start time:                ', a)", current_time
					print format_string, total_time*num_plot_ploints, &
						total_time*num_plot_ploints/60
					
					write(info_output_unit, "('*** SIMULATION TIME ***')")
					write(info_output_unit, "('Start time:                ', a)") &
						current_time
					write(info_output_unit, format_string) total_time*num_plot_ploints, &
						total_time*num_plot_ploints/60
					
					estimated_time = .false.
				end if

			end if

			! Compute next time step using Velocity Verlet algorithm
			call velocity_verlet_step(i, rt, t0, dt, ti, ri, vi, ai)

			! Track approach to scattering center
			if (approaching_center) then
				distance_scattering_center = norm2(ri)
				
				if (distance_scattering_center .gt. &
					closest_distance_scattering_center) then
					approaching_center = .false.
				else
					closest_distance_scattering_center = distance_scattering_center
					closest_position_scattering_center = ri
				end if
				
			end if

		end do

		! Compute total simulation execution time
		call cpu_time(end_time)
		total_time = total_time + (end_time - start_time)

		call fdate(current_time)
		format_string = "('Total simulation time:     ', f8.2, '[s]', &
			&f8.2, '[min]')"
		
		print "('End time:                  ', a)", current_time
		print format_string, total_time, total_time/60
		
		write(info_output_unit, "('End time:                  ', a)") current_time
		write(info_output_unit, format_string) total_time, total_time/60

		write(info_output_unit, *)   ! Blank line on Output Info File
		write(info_output_unit, "('*** SCATTERING ANGLES COMPUTATION ***')")

		! Theoretical Scattering Angle
		scat_angle_theoretical = 2*datan(1/(2*K0*b))
		format_string = "('Theoretical Scattering Angle: ', f20.16, 'º')"
		write(info_output_unit, format_string) scat_angle_theoretical*180/PI
		
		! Plot Scattering Angle
		scat_angle_plot = datan2(yf,xf)
!		if (yf .lt. 0._dp) scat_angle_plot = 2*PI + scat_angle_plot
		format_string = "('Plot Scattering Angle:        ', f20.16, 'º')"
		write(info_output_unit, format_string) scat_angle_plot*180/PI
		
		! Simulation Scattering Angle
		scat_angle_sim = datan2(ri(2),ri(1))
!		if (ri(2) .lt. 0._dp) scat_angle_sim = 2*PI + scat_angle_sim
		format_string = "('Simulation Scattering Angle:  ', f20.16, 'º')"
		write(info_output_unit, format_string) scat_angle_sim*180/PI
		
		write(info_output_unit, *)   ! Blank line on Output Info File
		write(info_output_unit, "('*** SCATTERING ANGLES COMPARISON ***')")
		write(info_output_unit, "('Percent error between...')")

		! Comparison between simulation and plot scattering angles
		format_string = "('  simulation and plot angles:        ', d12.4, '%')"
		write(info_output_unit, format_string) &
			100*dabs(scat_angle_sim - scat_angle_plot)/scat_angle_plot
			
		! Comparison between simulation and theoretical scattering angles
		format_string = "('  simulation and theoretical angles: ', d12.4, '%')"
		write(info_output_unit, format_string) &
			100*dabs(scat_angle_sim - scat_angle_theoretical)/scat_angle_theoretical
			
		! Comparison between plot and theoretical scattering angle
		format_string = "('  plot and theoretical angles:       ', d12.4, '%')"
		write(info_output_unit, format_string) &
			100*dabs(scat_angle_plot - scat_angle_theoretical)/scat_angle_theoretical

		write(info_output_unit, *)		! Blank line on Output Info File

		write(info_output_unit, "('*** CLOSEST DISTANCE AND POSITION &
			&TO SCATTERING CENTER ***')")
		format_string = "('Closest distance:   ', d12.4, '[a0]')"
		write(info_output_unit, format_string) closest_distance_scattering_center
		format_string = "('Closest postion, x: ', d12.4, '[a0]')"
		write(info_output_unit, format_string) closest_position_scattering_center(1)
		format_string = "('Closest postion, y: ', d12.4, '[a0]')"
		write(info_output_unit, format_string) closest_position_scattering_center(2)

		! Close output files
		close(output_unit)
		close(info_output_unit)

		print*

	end do

	! Close input files opened for numeric and char input, delete duplicate file
	close(input_values_unit)
	close(input_chars_unit)
	call system ("rm "//aux_file)
	
	contains
	!=============================================================================
	! Subroutine: parameter_initialization
	! Purpose   : Initialize simulation parameters by reading and converting 
	!             values for initial conditions, estimating number of time steps.
	! Arguments :
	!   - integer(i8), intent(in) :: info_output_unit
	!       Unit number for the simulation information output file.
	!   - integer(i8), intent(inout) :: num_plot_ploints
	!       Number of points to be plotted. At most equal to number of time steps.
	!   - real(dp), intent(inout) :: K0
	!       Initial kinetic energy of the projectile. On input, in kiloelectron
	!       volts (keV). On output, converted to atomic units (Hartree, Eh).
	!   - real(dp), intent(inout) :: x0
	!       Initial horizontal distance between particles. On input, in angstroms
	!       (Å). On output, converted to atomic units (Bohr radius, a0).
	!   - real(dp), intent(inout) :: b
	!       Impact parameter. On input, in angstroms (Å). On output, converted
	!       to atomic units (Bohr radius, a0).
	!   - real(dp), intent(inout) :: dt
	!       Time step size in atomic units of time (aut).
	!   - real(dp), intent(out) :: r0(3)
	!       Initial position vector of the projectile electron (a0).
	!   - real(dp), intent(out) :: v0(3)
	!       Initial velocity vector of the projectile electron (a0/aut).
	!   - integer(i8), intent(out) :: max_iterations
	!       Estimated number of time steps in the simulation.
	!=============================================================================
	subroutine parameter_initialization &
		(info_output_unit, num_plot_ploints, K0, x0, b, dt, r0, v0, max_iterations)
		implicit none

		! Input/Output variables
		integer(i8), intent(in) :: info_output_unit
		integer(i8), intent(inout) :: num_plot_ploints
		real(dp), intent(inout) :: K0		! Initial kinetic energy (keV)
		real(dp), intent(inout) :: x0		! Initial horizontal distance (Å)
		real(dp), intent(inout) :: b		! Impact parameter (Å)
		real(dp), intent(inout) :: dt		! Time step size (aut)
		real(dp), intent(out) :: r0(3)	! Initial position vector (a0)
		real(dp), intent(out) :: v0(3)	! Initial velocity vector (a0/aut)
		integer(i8), intent(out) :: max_iterations

		! Local variables
		real(dp) :: conv_aux	! Auxiliary variable for unit conversion
		real(dp) :: v					! Velocity magnitude (a0/aut)
		real(dp) :: a					! Approximate constant acceleration (a0/aut^2)
		real(dp) :: tf				! Estimated final simulation time (aut)
		character(len=80) :: format_param_write_string

		! Output basic simulation configuration
		
		write(info_output_unit, "('*** SIMULATION PARAMETERS ***')")
		write(info_output_unit, "('Number of points to be plotted: ', i6)") &
			num_plot_ploints

		! Convert and print initial kinetic energy, K0
		conv_aux = K0
		K0 = K0*1.d3/27.21139	! Conversion from keV to Eh
		format_param_write_string = &
			"('Initial kinetic energy, K0:     ', d12.4, '[keV] =', d12.4, '[Eh]')"
		write(info_output_unit, format_param_write_string) conv_aux, K0

		! Convert and print initial horizontal distance, x0
		conv_aux = x0
		x0 = x0/0.5291772	! Conversion from Å to a0
		format_param_write_string = &
		"('Initial horizontal distance, x0:', d12.4, '[Å]   =', d12.4, '[a0]')"
		write(info_output_unit, format_param_write_string) conv_aux, x0

		! Convert and print impact parameter, b
		conv_aux = b
		b = b/0.5291772	! Conversion from Å to a0
		format_param_write_string = &
		"('Impact parameter, b:            ', d12.4, '[Å]   =', d12.4, '[a0]')"
		write(info_output_unit, format_param_write_string) conv_aux, b

		! Convert and print time step size, dt
		conv_aux = dt * 2.418d-17	! Convert from aut to seconds
		format_param_write_string = &
		"('Time step size, dt:             ', d12.4, '[aut] =', d12.4, '[s]')"
		write(info_output_unit, format_param_write_string) dt, conv_aux

		! Initialize position vector: projectile starts at (-x0, b, 0)
		r0 = (/-x0, b, 0._dp/)

		! Initialize velocity vector: motion along +x axis
		v0 = (/dsqrt(2*K0), 0._dp, 0._dp/)

		! Estimate final simulation time assuming constant acceleration
		a = 1.0_dp / norm2(r0)
		v = norm2(v0)
		tf = (-v + dsqrt(v**2 - 2*a*x0))/a
		tf = -2*tf

		! Calculate estimated number of simulation steps
		max_iterations = dint(tf/dt)

		! Adjust number of plot points if it exceeds maximum number of iterations
		if (max_iterations .lt. num_plot_ploints) num_plot_ploints = max_iterations

	end subroutine parameter_initialization

	!=============================================================================
	! Subroutine: compute_theoretical_trajectory
	! Purpose   : Compute the theoretical scattering trajectory of a projectile 
	!             particle moving towards a target particle of the same charge
	!             using the corresponding orbit equation: a hyperbola in polar
	!             coordinates. The geometric parameters of the hyperbola are
	!             computed from the physical parameters. The trajectory
	!             and geometric parameters are written to output files.
	!             
	! Arguments :
	!   - integer(i8), intent(in) :: num_plot_ploints
	!       Number of points to be plotted. At most equal to number of time steps.
	!   - integer(i8), intent(in) :: output_unit
	!       Unit number for the simulation output file.
	!   - integer(i8), intent(in) :: info_output_unit
	!       Unit number for the simulation information output file.
	!   - real(dp), intent(in) :: r0(3)
	!       Initial position vector of the projectile electron (a0).
	!   - real(dp), intent(in) :: K0
	!       Initial kinetic energy of the projectile (Hartree, Eh).
	!   - real(dp), intent(out) :: xf
	!       Final x-coordinate of the trajectory (a0).
	!   - real(dp), intent(out) :: yf
	!       Final y-coordinate of the trajectory (a0).
	!=============================================================================
	subroutine compute_theoretical_trajectory &
		(num_plot_ploints, output_unit, info_output_unit, r0, K0, xf, yf)
		implicit none

		! Input/Output variables
		integer(i8), intent(in) :: num_plot_ploints
		integer(i8), intent(in) :: output_unit, info_output_unit
		real(dp), intent(in) :: r0(3), K0				! Physical parameters
		real(dp), intent(out) :: xf, yf					! Last point coordinates

		! Local variables
		real(dp) :: l, e, a, b, c									! Hyperbola geometric parameters
		real(dp) :: x0, y0												! Initial position coordinates
		real(dp) :: phi0, phif, dphi, phii, ri		! Angular and radial coordinates
		real(dp) :: xi, yi												! Cartesian coordinates
		real(dp) :: theta0												! Rotation angle
		real(dp) :: delta, den, num								! Auxiliary variables
		integer :: i

		! Extract initial values
		x0 = r0(1)	! Initial horizontal distance
		y0 = r0(2)	! Impact parameter

		! Compute hyperbola geometric parameters
		l = 2*K0*(y0**2)									! Semi-latus rectum
		e = dsqrt(1 + 4*(K0**2)*(y0**2))	! Eccentricity
		a = 1/(2*K0)											! Semi-major axis
		b = y0														! Semi-minor axis
		c = dsqrt(a*a + b*b)							! Distance to the center

		! Print hyperbola parameters
		write(info_output_unit, "('*** THEORETICAL TRAJECTORY PARAMETERS ***')")
		write(info_output_unit, "('Hyperbola geometric parameters')")
		write(info_output_unit, "('l:', e12.4, '[a0]')") l
		write(info_output_unit, "('e:', e12.4, '[--]')") e
		write(info_output_unit, "('a:', e12.4, '[a0]')") a
		write(info_output_unit, "('b:', e12.4, '[a0]')") b
		write(info_output_unit, "('c:', e12.4, '[a0]')") c
		

		! Compute trajectory using hyperbola's polar equation (left branch only)
		delta = dabs(datan2(y0, x0) - PI)			! Small delta to avoid undefined ri
		phi0 = PI - 2*dacos(1/e) + delta			! Initial angle to start plotting
		phif = PI - delta											! Final angle
		dphi = (phif - phi0)/num_plot_ploints	! Angular step between points
		
		theta0 = PI - dacos(1/e)							! Hyperbola rotation angle

		do i = 0, num_plot_ploints
			phii = phi0 + i*dphi
			! Cartesian coordinates trajectory equation:
			! ri = -l/(1 - e*dcos(phii - theta0)) rewritten below to avoid loss of
			! significance for small values of b, using difference of squares trick
			num = -l*(1 + e*dcos(phii - theta0))
			den = 1 - (e*dcos(phii - theta0))**2
			ri = num/den

			xi = ri*dcos(phii)
			yi = ri*dsin(phii)

			write(output_unit, *) xi, yi

			if (i == 0) then
				! Save first point as theoretical final coordinates
				xf = xi
				yf = yi
			end if
			
		end do

		! Add separation in output file to plot using gnuplot
		write(output_unit, *)
		write(output_unit, *)

	end subroutine compute_theoretical_trajectory

	!=============================================================================
	! Subroutine: compute_conserved_quantities
	! Purpose   : Calculate conserved mechanical quantities for a Coulomb
	!              interaction: total energy (kinetic + potential) and the
	!              magnitude of the angular-momentum vector.
	! Arguments :
	!   - real(dp), intent(in)  :: r(3)
	!       Position vector of the projectile electron (a0).
	!   - real(dp), intent(in)  :: v(3)
	!       Velocity vector of the projectile electron (a0/aut).
	!   - real(dp), intent(out) :: U
	!       Total mechanical energy (Hartree, Eh).
	!   - real(dp), intent(out) :: L
	!       Magnitude of angular momentum (a0^2/aut, ħ).
	!=============================================================================
	subroutine compute_conserved_quantities(r, v, U, L)
		 implicit none

		 ! Input/Output variables
		 real(dp), intent(in) :: r(3)		! Position of the electron (a0)
		 real(dp), intent(in) :: v(3)		! Velocity of the electron (a0/aut)
		 real(dp), intent(out) :: U			! Total mechanical energy (Eh)
		 real(dp), intent(out) :: L			! Angular momentum magnitude (ħ)

		 ! Local variables
		 real(dp) :: Uk, Ue		! Energy components
		 real(dp) :: L_vec(3)	! Angular-momentum vector

		 ! Energy calculation
		 Uk = 0.5*(v(1)**2 + v(2)**2 + v(3)**2)	! Kinetic energy
		 Ue = 1/norm2(r)												! Electrostatic potential energy
		 U = Uk + Ue

		 ! Angular-momentum vector   L = r × v
		 L_vec(1) = r(2)*v(3) - r(3)*v(2)	! y*vz - z*vy
		 L_vec(2) = r(3)*v(1) - r(1)*v(3)	! z*vx - x*vz
		 L_vec(3) = r(1)*v(2) - r(2)*v(1)	! x*vy - y*vx
		 L = norm2(L_vec)

	end subroutine compute_conserved_quantities

	!=============================================================================
	! Subroutine: acceleration_due_to_electron
	! Purpose   : Calculate the acceleration vector experienced by a projectile
	!             electron due to the electrostatic interaction with a stationary
	!             target electron, based on their positions.
	! Arguments :
	!   - real(dp), intent(in) :: rp(3)
	!       Position vector of the projectile electron in atomic units (a0).
	!   - real(dp), intent(in) :: rt(3)
	!       Position vector of the target electron in atomic units (a0/aut).
	!   - real(dp), intent(out) :: a(3)
	!       Acceleration vector experienced by the projectile electron due to
	!       the target electron in atomic units (a0/aut^2).
	!=============================================================================
	subroutine acceleration_due_to_electron(rp, rt, a)
		implicit none

		! Input/Output variables
		real(dp), intent(in) :: rp(3)	! Position of the projectile electron (a0)
		real(dp), intent(in) :: rt(3)	! Position of the target electron (a0)
		real(dp), intent(out) :: a(3)	! Resulting acceleration vector (a0/aut^2)

		! Local variables
		real(dp) :: rs(3)	! Separation vector between the electrons (a0)
		real(dp) :: r			! Magnitude of the separation vector (a0)

		! Compute the separation vector between electrons and its magnitude
		rs = rp - rt
		r = norm2(rs)

		! Determine the acceleration using Coulomb's law in atomic units
		a = rs/(r**3)

	end subroutine acceleration_due_to_electron

	!=============================================================================
	! Subroutine: velocity_verlet_step
	! Purpose   : Advance the position, velocity, and acceleration of the moving
	!             particle by one Velocity-Verlet time step.
	! Arguments :
	!   - integer(i8), intent(in) :: i
	!       Current iteration index.
	!   - real(dp), intent(in) :: rt(3)
	!       Position of the stationary target charge (a0).
	!   - real(dp), intent(in) :: t0
	!       Initial simulation time (aut).
	!   - real(dp), intent(in) :: dt
	!       Time-step size (aut).
	!   - real(dp), intent(inout) :: t
	!       Current simulation time (aut).
	!   - real(dp), intent(inout) :: r(3)
	!       Position vector of the moving charge (a0).
	!   - real(dp), intent(inout) :: v(3)
	!       Velocity vector of the moving charge (a0 / aut).
	!   - real(dp), intent(inout) :: a(3)
	!       Acceleration acting on the moving charge (a0 / aut²).
	!=============================================================================
	subroutine velocity_verlet_step(i, rt, t0, dt, t, r, v, a)
		 implicit none

		 ! Input/Output variables
		 integer(i8), intent(in) :: i
		 real(dp), intent(in) :: rt(3), t0, dt
		 real(dp), intent(inout) :: t, r(3), v(3), a(3)

		 ! Time update
		 t = t0 + i * dt

		 ! Half-step velocity update
		 v = v + 0.5*a*dt
		 
		 ! Position update
		 r = r + v*dt
		 
		 ! Full-step acceleration update at new position
		 call acceleration_due_to_electron(r, rt, a)
		 
		 ! Second half-step velocity update
		v = v + 0.5*a*dt

	end subroutine velocity_verlet_step
	
end program rutherford_scattering_test_simulations