! Rutherford Scattering Test Simulations
! This program is intended as a proof of concept to become familiar with the
! Velocity Verlet algorithm, Rutherford Scattering, and how the different
! geometric and physical parameters used impact the resulting simulations.
program rutherfhord_scattering_test_simulations
	
	use, intrinsic:: iso_fortran_env, only: stdin=>input_unit
	implicit none
	! Kind type parameters for increased real precision and integer length
	! Single precision reals, 6 digits, range 10**(-37) to 10**(37)-1; 32 bits
	integer, parameter :: sp = selected_real_kind(6, 37)
	! Double precision reals, 15 digits, range 10**(-307) to 10**(307)-1; 64 bits
	integer, parameter :: dp = selected_real_kind(15, 307)
	! Long length for integers, range -2**63 to 2**63-1; 64 bits
	integer, parameter :: i8 = selected_int_kind(18)
	! General variables
	real(dp), parameter :: PI = dacos(-1.d0)
	real(dp), parameter :: t0 = 0._dp
	integer(i8) :: i, j!, option	!Program structure parameters
	integer(i8) :: N 						!Number of points to be plotted/simulated
	real(dp) :: K0, K0f, hd, hdf	!Fixed simulation parameters
	real(dp) :: b, dt							!Variable simulation parameters
	real(dp) :: rt(3)						!Target electron position
	real(dp) :: r0(3), v0(3)		!Simulation position and velocity
	integer(i8) :: T						!Number of iterations
	real(dp) :: a, c, e!, b			!Theoretical trajectory geometric parameters
	real(dp) :: xf, yf		!Last theoretical trajectory point generated
	real(dp) :: ti, ri(3)	!Time and space simulation variables
	real(dp) :: vi(3)			!Velocity simulation variables
	real(dp) :: ai(3) 		!Acceleration simulation variables
	logical :: estimated_time
	character(len=24) :: current_time
	real(sp) :: startT, endT, execTime			!Program timer variables
	integer(i8) :: NS, k
	character(len=*), parameter :: input_file = "test_input.in", aux_file = "aux.in"
	character(len=10) :: K0c, hdc, bc, dtc!AT MOST 10 CHARACTERS FOR INPUT VALUES!!!
	character(len=:), allocatable :: K0ct, hdct, bct, dtct
	character(len=:), allocatable :: output_file, output_file_info
	integer(i8), parameter :: inu = 11, icu = 12	!Input file as numbers and chars
	integer(i8), parameter :: ou = 13, oiu = 14	!Output file for data and info
	character(len=80) :: FMTS		!Format string
	
	!*******************************************************************************
	!Reading fixed parameters from input file
	!Open input file twice: unit 11 to use as numbers, unit 12 to use as chars
	call system ("cp "//input_file//" "//aux_file)
	open(unit=inu, file=input_file, status='unknown')
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
		
		estimated_time = .true.
		!Start timer RIGHT before the first iteration
		call cpu_time(startT)

		do i=1,T
			!Plotting only N points
			if ( (mod(i,T/N) .eq. 0) .and. j .lt. N ) then

				!Write values to file
				write(ou,*) ti, ri
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
		
	end do

	close(inu)
	close(icu)
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
		integer(i8), intent(in)    :: info_output_unit
		integer(i8), intent(inout) :: num_plot_ploints
		real(dp), intent(inout)    :: K0		! Initial kinetic energy (keV)
		real(dp), intent(inout)    :: x0		! Initial horizontal distance (Å)
		real(dp), intent(inout)    :: b			! Impact parameter (Å)
		real(dp), intent(inout)    :: dt		! Time step size (aut)
		real(dp), intent(out)      :: r0(3)	! Initial position vector (a0)
		real(dp), intent(out)      :: v0(3)	! Initial velocity vector (a0/aut)
		integer(i8), intent(out)   :: max_iterations

		! Local variables
		real(dp) :: conv_aux	! Auxiliary variable for unit conversion
		real(dp) :: v         ! Velocity magnitude (a0/aut)
		real(dp) :: a         ! Approximate constant acceleration (a0/aut^2)
		real(dp) :: tf        ! Estimated final simulation time (aut)
		character(len=80) :: format_param_write_string

		! Output basic simulation configuration
		write(info_output_unit, "('*** SIMULATION PARAMETERS ***')")
		write(info_output_unit, "('Number of points to be plotted, N: ', i6)") &
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
		"('Initial horizontal distance, hd:', d12.4, '[Å]   =', d12.4, '[a0]')"
		write(info_output_unit, format_param_write_string) conv_aux, x0

		! Convert and print impact parameter, b
		conv_aux = b
		b = b/0.5291772	! Conversion from Å to a0
		format_param_write_string = &
		"('Impact parameter, b:            ', d12.4, '[Å]   =', d12.4, '[a0]')"
		write(info_output_unit, format_param_write_string) conv_aux, b

		! Convert and print time step size, dt
		conv_aux = dt * 2.418d-17  ! Convert from aut to seconds
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
	
	
end program rutherfhord_scattering_test_simulations