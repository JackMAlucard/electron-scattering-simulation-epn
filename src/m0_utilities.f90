module m0_utilities
	implicit none

	! NUMERICAL STORAGE SIZE PARAMETERS FOR REAL AND INTEGER VALUES **************
	! Double precision real numbers,
	! 15 digits, range 10**(-307) to 10**(307)-1; 64 bits
	integer, parameter :: dp = selected_real_kind(15, 307)
	! Long length for integers, range -2**63 to 2**63-1; 64 bits
	integer, parameter :: i8 = selected_int_kind(18)

	! NUMERICAL CONSTANTS ********************************************************
	real(dp), parameter :: PI = dacos(-1._dp)
	real(dp), parameter :: INTERATOMIC_DIST_SIO2 = 1.77/0.5291772! a0, 1.77 Å
	real(dp), parameter :: INTERATOMIC_DIST_SIO2_INV = 1/INTERATOMIC_DIST_SIO2
	real(dp), parameter :: MATERIAL_HEIGHT_SIO2 = 0.5*INTERATOMIC_DIST_SIO2
	real(dp), parameter :: MEAN_FREE_PATH_SIO2 = 33.74649829! a0
	real(dp), parameter :: CROSS_SECTION_SIO2 = 1/MEAN_FREE_PATH_SIO2 ! Macro
	real(dp), parameter :: CELL_SCALE_FACTOR = 3
	real(dp), parameter :: CELL_LENGTH = CELL_SCALE_FACTOR*INTERATOMIC_DIST_SIO2
	real(dp), parameter :: CELL_LENGTH_INV = 1/CELL_LENGTH
	real(dp), parameter :: EFFECTIVE_DISTANCE = 30! a0
	real(dp), parameter :: MAX_EQUIVALENT_CHARGE = 180

	contains
	! UNIT CONVERSION SUBROUTINES ************************************************
	!=======================================================================
	! Subroutine: kev_to_atomic_energy_conversion
	! Purpose   : Convert energy from kilo-electron volts (keV) to atomic 
	!             units (Hartrees).
	! Arguments :
	!   - real(dp), intent(inout) :: energy
	!       The energy value to be converted. On input, it should be in 
	!       kilo-electron volts (keV). On output, it will be converted to 
	!       atomic units (Hartrees).
	!=======================================================================
	! Energy conversion from keV to hartree Eh (atomic)
	subroutine kev_to_atomic_energy_conversion(energy)
		implicit none
		real(dp), intent(inout) :: energy
		energy = energy*1.d3/27.21139
	end subroutine kev_to_atomic_energy_conversion
	!=======================================================================
	! Subroutine: si_to_atomic_distance_conversion
	! Purpose   : Convert distance from SI units (meters) to atomic units 
	!             (Bohr radius).
	! Arguments :
	!   - real(dp), intent(inout) :: distance
	!       The distance value to be converted. On input, it should be in 
	!       meters (SI units). On output, it will be converted to atomic 
	!       units (Bohr radius).
	!=======================================================================
	! Distance conversion from meters m (SI) to Bohr radii a0 (atomic)
	subroutine si_to_atomic_distance_conversion(distance)
		implicit none
		real(dp), intent(inout) :: distance
		distance = distance/0.5291772d-10
	end subroutine si_to_atomic_distance_conversion
	!=======================================================================
	! Subroutine: angstrom_to_atomic_distance_conversion
	! Purpose   : Convert distance from Angstroms to atomic units 
	!             (Bohr radius).
	! Arguments :
	!   - real(dp), intent(inout) :: distance
	!       The distance value to be converted. On input, it should be in 
	!       Angstroms. On output, it will be converted to atomic units 
	!       (Bohr radius).
	!=======================================================================
	! Distance conversion from angstroms to Bohr radii a0 (atomic)
	subroutine angstrom_to_atomic_distance_conversion(distance)
		implicit none
		real(dp), intent(inout) :: distance
		distance = distance/0.5291772
	end subroutine angstrom_to_atomic_distance_conversion
	!=======================================================================
	! Subroutine: atomic_to_si_distance_conversion
	! Purpose   : Convert distance from atomic units (Bohr radius) to 
	!             SI units (meters).
	! Arguments :
	!   - real(dp), intent(inout) :: distance
	!       The distance value to be converted. On input, it should be in 
	!       atomic units (Bohr radius). On output, it will be converted 
	!       to meters (SI units).
	!=======================================================================
	! Distance conversion from Bohr radii a0 (atomic) to meters (SI)
	subroutine atomic_to_si_distance_conversion(distance)
		implicit none
		real(dp), intent(inout) :: distance
		distance = distance*0.5291772d-10
	end subroutine atomic_to_si_distance_conversion

	! RANDOM NUMBER GENERATION SUBROUTINES****************************************
	!=======================================================================
	! Subroutine: random_std_uniform
	! Purpose   : Generate a random number following the standard uniform 
	!             distribution over the range (0,1] to avoid possible error 
	!             if 0 is generated.
!	/ Generate a random number uniformly distributed between 
!             0 and 1, using the standard uniform distribution.
	
	! Arguments :
	!   - real(dp), intent(out) :: u
	!       The generated random number following the standard uniform 
	!       distribution over the range (0,1].
!	The generated random number uniformly distributed between 
!       0 and 1.
	!=======================================================================
	! Subroutine which generates a random number following the standard uniform
	! distribution over the range (0,1] to avoid possible error if 0 is generated.
	! Adapted from: https://masuday.github.io/fortran_tutorial/random.html
	subroutine random_std_uniform(u)
		implicit none
		real(dp), intent(out) :: u
		real(dp) :: r
		call random_number(r)
		u = 1 - r
	end subroutine random_std_uniform
	!=======================================================================
	! Subroutine: random_std_normal
	! Purpose   : Generate two independent random numbers following the 
	!             standard normal distribution (mean = 0, standard deviation = 1) 
	!             using the Box-Muller transform.
	! Arguments :
	!   - real(dp), intent(out) :: z1
	!       The first generated random number following the standard 
	!       normal distribution.
	!   - real(dp), intent(out) :: z2
	!       The second generated random number following the standard 
	!       normal distribution.
	!=======================================================================
	! Subroutine which generates two random numbers following the standard 
	! normal distribution (mean = 0 and standard deviation = 1).
	! Adapted from: https://masuday.github.io/fortran_tutorial/random.html
	subroutine random_std_normal(z1, z2)
		implicit none
		real(dp), intent(out) :: z1, z2
		real(dp) :: u1, u2
		call random_std_uniform(u1)
		call random_std_uniform(u2)
		z1 = dsqrt(-2*dlog(u1))*dcos(2*PI*u2)
		z2 = dsqrt(-2*dlog(u1))*dsin(2*PI*u2)
	end subroutine random_std_normal
	!=======================================================================
	! Subroutine: random_normal
	! Purpose   : Generate two independent random numbers following a 
	!             normal distribution with specified mean and standard 
	!             deviation.
	! Arguments :
	!   - real(dp), intent(in) :: mu
	!       The mean of the normal distribution.
	!   - real(dp), intent(in) :: sigma
	!       The standard deviation of the normal distribution.
	!   - real(dp), intent(out) :: x1
	!       The first generated random number following the specified 
	!       normal distribution.
	!   - real(dp), intent(out) :: x2
	!       The second generated random number following the specified 
	!       normal distribution.
	!=======================================================================
	! Subroutine which generates two random number following the general
	! normal distribution with mean 'mu' and standard deviation 'sigma'.
	subroutine random_normal(mu, sigma, x1, x2)
		implicit none
		real(dp), intent(in) :: mu, sigma
		real(dp), intent(out) :: x1, x2
		real(dp) :: z1, z2
		call random_std_normal(z1, z2)
		x1 = mu + sigma*z1
		x2 = mu + sigma*z2
	end subroutine random_normal
	!=======================================================================
	! Subroutine: random_exponential
	! Purpose   : Generate a random number following an exponential 
	!             distribution with a specified rate parameter.
	! Arguments :
	!   - real(dp), intent(in) :: lambda
	!       The rate parameter (λ) of the exponential distribution.
	!   - real(dp), intent(out) :: x
	!       The generated random number following the exponential 
	!       distribution.
	!=======================================================================
	! Subroutine which generates a random number following the exponential
	! probability distribution with parameter 'lambda' > 0.
	! Adapted from:
	! https://www.eg.bucknell.edu/~xmeng/Course/CS6337/Note/master/node50.html
	subroutine random_exponential(lambda, x)
		implicit none
		real(dp), intent(in) :: lambda
		real(dp), intent(out) :: x
		real(dp) :: u
		call random_std_uniform(u)
		x = -dlog(u)/lambda
	end subroutine random_exponential
	
	! VECTOR TRANSFORMATION SUBROUTINES ******************************************
	!=======================================================================
	! Subroutine: vector_translation
	! Purpose   : Translate a 3-dimensional vector by adding a translation 
	!             vector to it.
	! Arguments :
	!   - real(dp), intent(in) :: translation_vector(3)
	!       The translation vector to be added.
	!   - real(dp), intent(inout) :: vector(3)
	!       The vector to be translated. On input, it contains the original 
	!       vector. On output, it contains the translated vector.
	!=======================================================================
	subroutine vector_translation(translation_vector, vector)
		implicit none
		real(dp), intent(in) :: translation_vector(3)
		real(dp), intent(inout) :: vector(3)
		vector = vector + translation_vector
	end subroutine vector_translation
	!=======================================================================
	! Subroutine: rotation_about_x_axis
	! Purpose   : Rotate a 3-dimensional vector about the x-axis by a given 
	!             angle.
	! Arguments :
	!   - real(dp), intent(in) :: theta
	!       The angle (in radians) by which the vector is to be rotated 
	!       about the x-axis.
	!   - real(dp), intent(inout) :: vector(3)
	!       The 3-dimensional vector to be rotated. On input, it contains 
	!       the original vector. On output, it contains the rotated vector.
	!=======================================================================
	subroutine rotation_about_x_axis(theta, vector)
		implicit none
		real(dp), intent(in) :: theta
		real(dp), intent(inout) :: vector(3)
		real(dp) :: rotation_matrix(3,3), rotated_vector(3)
		integer :: i
		rotation_matrix(1,:) = (/1._dp, 0._dp, 0._dp/)
		rotation_matrix(2,:) = (/0._dp, dcos(theta), -dsin(theta)/)
		rotation_matrix(3,:) = (/0._dp, dsin(theta), dcos(theta)/)
		rotated_vector = 0
		do i = 1, 3
			rotated_vector(:) = rotated_vector(:) + rotation_matrix(:,i)*vector(i)
		end do
		vector = rotated_vector
	end subroutine rotation_about_x_axis
	
	! INPUT AND OUTPUT FILE HANDLING SUBROUTINES *********************************
	!=======================================================================
	! Subroutine: read_input_parameters
	! Purpose   : Read input parameters from a file named 'input.txt'.
	! Arguments :
	!   - integer(i8), intent(out) :: num_electrons
	!       Number of electrons in the simulation.
	!   - real(dp), intent(out) :: beam_energy
	!       Energy of the electron beam in kilo-electron volts (keV).
	!   - real(dp), intent(out) :: energy_spread
	!       Energy spread of the electron beam in percentage.
	!   - real(dp), intent(out) :: spot_size_factor
	!       Factor determining the spot size of the electron beam.
	!   - real(dp), intent(out) :: beam_target_distance
	!       Distance from the electron beam source to the target material 
	!       surface in Angstroms (Å).
	!   - real(dp), intent(out) :: grazing_angle
	!       Grazing angle of the electron beam in degrees (º).
	!   - integer(i8), intent(out) :: material_boundaries(3)
	!       Indices of the material boundaries.
	!   - integer(i8), intent(out) :: num_plot_points
	!       Number of plot points for visualization.
	!   - real(dp), intent(out) :: dt
	!       Time step size in atomic units of time.
	!=======================================================================
	! Subroutine that reads the main simulation parameters from the 'input.txt' 
	! file. Parameter units are specified as comments in the 'input.txt' file.
	subroutine read_input_parameters &
		(num_electrons, beam_energy, energy_spread, spot_size_factor, &
		beam_target_distance, grazing_angle, material_boundaries, &
		beam_model_output_saving_enabled, num_plot_ploints, dt)
		implicit none
		integer(i8), intent(out) :: num_electrons, material_boundaries(3)
		real(dp), intent(out) :: beam_energy, energy_spread, beam_target_distance
		real(dp), intent(out) :: spot_size_factor
		real(dp), intent(out) :: grazing_angle, dt
		logical, intent(out) :: beam_model_output_saving_enabled
		integer(i8), intent(out) :: num_plot_ploints
		open(unit=42, file='input.txt', status='old', action='read')
			read(42, *) ! beam_energy (keV)
			read(42, *) beam_energy
			read(42, *) ! grazing_angle (º)
			read(42, *) grazing_angle
			read(42, *) ! num_electrons (int)
			read(42, *) num_electrons
			read(42, *) ! energy_spread (%)
			read(42, *) energy_spread
			read(42, *) ! spot_size_factor (int)
			read(42, *) spot_size_factor
			read(42, *) ! beam_target_distance (Å)
			read(42, *) beam_target_distance
			read(42, *) ! material_boundaries (int, int, int)
			read(42, *) material_boundaries
			read(42, *) ! num_plot_ploints (int)
			read(42, *) num_plot_ploints
			read(42, *) ! dt (atomic units of time), time step size
			read(42, *) dt
			read(42, *) ! Save beam model vectors as output? (boolean)
			read(42, *) beam_model_output_saving_enabled
		close(42)
	end subroutine read_input_parameters
	!=======================================================================
	! Subroutine: open_output_files
	! Purpose   : Open output files for writing electron trajectories, final
	!             positions of embedded and scattered electrons, scattering
	!             angles of scattered electrons, all final electron positions,
	!             time evolution of electrons' final state, and iteration 
	!             and total simulation times.
	! Arguments :
	!   - integer :: output_unit
	!       The unit number to be used as a base for opening the output files.
	!=======================================================================
	! Subroutine that opens the files that will store the output results of the
	! simulations.
	subroutine open_output_files(output_unit)
		implicit none
		integer :: output_unit
		! Open file to write electron trajectories
		open(unit=output_unit+1, file="electron_trajectories.dat", &
			status='replace', action='write')
		! Open file to write the final position of embedded electrons
		open(unit=output_unit+2, file='embedded_positions.dat', &
			status='replace', action='write')
		! Open file to write the final position of scattered electrons
		open(unit=output_unit+3, file='scattered_positions.dat', &
			status='replace', action='write')
		! Open file to write the scattered electrons scattering angles 
		open(unit=output_unit+4, file='scattering_angles.dat', &
			status='replace', action='write')
		! Open file to write all final electron positions
		open(unit=output_unit+5, file='final_electron_positions.dat', &
			status='replace', action='write')
		! Open file to write the time evolution of the electrons' final state
		open(unit=output_unit+6, file='final_state_evolution.dat', &
			status='replace', action='write')
		! Open file to write the iteration and total simulation times
		open(unit=output_unit+7, file='simulation_times.dat', &
			status='replace', action='write')
	end subroutine open_output_files
	!=======================================================================
	! Subroutine: write_simulation_status_information
	! Purpose   : Write simulation status information including the current 
	!             iteration count, number of electrons simulated, number 
	!             of electrons embedded, number of electrons scattered, 
	!             and total simulation time to a file named 
	!             'simulation_info.dat'.
	! Arguments :
	!   - integer :: output_unit
	!       The unit number used for writing the simulation status 
	!       information.
	!   - integer(i8) :: i
	!       The current iteration count.
	!   - integer(i8) :: num_electrons
	!       The total number of electrons to be simulated.
	!   - integer(i8) :: num_embedded
	!       The number of electrons embedded during the simulation.
	!   - integer(i8) :: num_scattered
	!       The number of electrons scattered during the simulation.
	!   - real(dp) :: total_time
	!       The total simulation time.
	!=======================================================================
	! Subroutine that writes the status of the simulation of after each 
	! trajectory computation is completed into a file.
	subroutine write_simulation_status_information &
		(output_unit, i, num_electrons, num_embedded, num_scattered, total_time)
		implicit none
		integer :: output_unit
		integer(i8) :: i, num_electrons
		integer(i8) :: num_embedded, num_scattered
		real(dp) :: total_time
		open(unit=output_unit, file='simulation_info.dat', &
			status='replace', action='write')
			write(output_unit, *) i, "out of", num_electrons, &
				"electron trajectories simulated"
			if (i .eq. num_electrons) write(output_unit, *) "SIMULATION COMPLETE!"
			write(output_unit, *) " Number of electrons embedded: ", num_embedded
			write(output_unit, *) " Number of electrons scattered:", num_scattered
			if (i .ne. num_electrons) then
				write(output_unit, *) "Simulation time thus far:", total_time
			else
				write(output_unit, *) "Total simulation time:", total_time
			end if
		close(output_unit)
	end subroutine write_simulation_status_information

end module m0_utilities