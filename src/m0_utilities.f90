module m0_utilities

	implicit none
	!=============================================================================
	! NUMERICAL CONSTANTS
	!=============================================================================

	! Kind type parameters for increased real precision and integer length
	! Double precision reals, 15 digits, range 10**(-307) to 10**(307)-1; 64 bits
	integer, parameter :: dp = selected_real_kind(15, 307)
	! Long length for integers, range -2**63 to 2**63-1; 64 bits
	integer, parameter :: i8 = selected_int_kind(18)

	! Constants related to physical properties and simulation parameters,
	! all initialized in the corresponding atomic units
	! Interatomic distance used in the SiO2 model, 1.77 Å, converted
	! to atomic units (Bohr radius, a0)
	real(dp), parameter :: INTERATOMIC_DIST_SIO2 = 1.77/0.5291772
	! Inverse of the interatomic distance in the SiO2 model (a0^-1)
	real(dp), parameter :: INTERATOMIC_DIST_SIO2_INV = 1/INTERATOMIC_DIST_SIO2
	! Vertical coordinate value corresponding to the height of the SiO2 material
	! top layer in the model (a0)
	real(dp), parameter :: MATERIAL_HEIGHT_SIO2 = 0.5*INTERATOMIC_DIST_SIO2
	! Approximated mean free path of electrons in the SiO2 model (a0)
	real(dp), parameter :: MEAN_FREE_PATH_SIO2 = 33.74649829
	! Approximated macroscopic cross section in the SiO2 model (a0)
	real(dp), parameter :: CROSS_SECTION_SIO2 = 1/MEAN_FREE_PATH_SIO2
	! Scaling factor for space partitioning cell size (real)
	real(dp), parameter :: CELL_SCALE_FACTOR = 3
	! Length of a space partitioning cell (a0)
	real(dp), parameter :: CELL_LENGTH = CELL_SCALE_FACTOR*INTERATOMIC_DIST_SIO2
	! Inverse of the space partitioning cell length (a0^-1)
	real(dp), parameter :: CELL_LENGTH_INV = 1/CELL_LENGTH
	! Effective distance for approximate interaction calculations used in
	! trajectory computation optimizations (a0)
	real(dp), parameter :: EFFECTIVE_DISTANCE = 30
	! Maximum charge value for super electrons, given as a real to be used in
	! mixed type division (e)
	real(dp), parameter :: MAX_SUPER_ELECTRON_CHARGE = 180
	! Value of the mathematical constant pi
	real(dp), parameter :: PI = dacos(-1._dp)

	contains
	!=============================================================================
	! UNIT CONVERSION SUBROUTINES
	!=============================================================================

	!=============================================================================
	! Subroutine: kev_to_atomic_energy_conversion
	! Purpose   : Convert energy from kiloelectron volts (keV) to atomic
	!             units (Hartree, Eh).
	! Arguments :
	!   - real(dp), intent(inout) :: energy
	!       The energy value to be converted. On input, it should be in
	!       kiloelectron volts (keV). On output, it will be converted to
	!       atomic units (Hartree, Eh).
	!=============================================================================
	subroutine kev_to_atomic_energy_conversion(energy)
		implicit none

		! Input/Output variable
		real(dp), intent(inout) :: energy

		! Convert energy from keV to atomic units (Hartree, Eh)
		energy = energy*1.d3/27.21139

	end subroutine kev_to_atomic_energy_conversion

	!=============================================================================
	! Subroutine: si_to_atomic_distance_conversion
	! Purpose   : Convert distance from SI units (meters, m) to atomic units
	!             (Bohr radius, a0).
	! Arguments :
	!   - real(dp), intent(inout) :: distance
	!       The distance value to be converted. On input, it should be in
	!       SI units (meters, m). On output, it will be converted to atomic
	!       units (Bohr radius, a0).
	!=============================================================================
	subroutine si_to_atomic_distance_conversion(distance)
		implicit none

		! Input/Output variable
		real(dp), intent(inout) :: distance

		! Convert distance from SI units (meters) to atomic units (Bohr radius, a0)
		distance = distance/0.5291772d-10

	end subroutine si_to_atomic_distance_conversion

	!=============================================================================
	! Subroutine: angstrom_to_atomic_distance_conversion
	! Purpose   : Convert distance from angstroms (Å) to atomic units
	!             (Bohr radius, a0).
	! Arguments :
	!   - real(dp), intent(inout) :: distance
	!       The distance value to be converted. On input, it should be in
	!       angstroms (Å). On output, it will be converted to atomic units
	!       (Bohr radius, a0).
	!=============================================================================
	subroutine angstrom_to_atomic_distance_conversion(distance)
		implicit none

		! Input/Output variable
		real(dp), intent(inout) :: distance

		! Convert distance from angstroms to atomic units (Bohr radius, a0)
		distance = distance/0.5291772

	end subroutine angstrom_to_atomic_distance_conversion

	!=============================================================================
	! Subroutine: atomic_to_si_distance_conversion
	! Purpose   : Convert distance from atomic units (Bohr radius, a0) to
	!             SI units (meters, m).
	! Arguments :
	!   - real(dp), intent(inout) :: distance
	!       The distance value to be converted. On input, it should be in
	!       atomic units (Bohr radius, a0). On output, it will be converted
	!       to SI units (meters, m).
	!=============================================================================
	subroutine atomic_to_si_distance_conversion(distance)
		implicit none

		! Input/Output variable
		real(dp), intent(inout) :: distance

		! Convert distance from atomic units (Bohr radius, a0) to SI units (meters)
		distance = distance*0.5291772d-10

	end subroutine atomic_to_si_distance_conversion

	!=============================================================================
	! RANDOM NUMBER GENERATION SUBROUTINES
	!=============================================================================

	!=============================================================================
	! Subroutine: random_std_uniform
	! Purpose   : Generate a random number following the standard uniform
	!             distribution over the range (0,1]. The range restriction
	!             helps prevent possible runtime-errors when using the
	!             resulting number.
	! Arguments :
	!   - real(dp), intent(out) :: u
	!       A random number following the standard uniform distribution
	!       over the range (0,1].
	! Adapted from :
	! https://masuday.github.io/fortran_tutorial/random.html
	!=============================================================================
	subroutine random_std_uniform(u)
		implicit none

		! Output variable
		real(dp), intent(out) :: u

		! Local variable
		real(dp) :: r

		! Generate a uniformly distributed random number in [0, 1),
		! and transform it so it is always on the range (0,1]
		call random_number(r)
		u = 1 - r

	end subroutine random_std_uniform

	!=============================================================================
	! Subroutine: random_std_normal
	! Purpose   : Generate two independent random numbers following the
	!             standard normal distribution (mean = 0, standard deviation = 1)
	!             using the Box-Muller transform.
	! Arguments :
	!   - real(dp), intent(out) :: z1
	!       The first independent random number generated following the
	!       standard normal distribution.
	!   - real(dp), intent(out) :: z2
	!       The second independent random number generated following the
	!       standard normal distribution.
	! Adapted from :
	! https://masuday.github.io/fortran_tutorial/random.html
	!=============================================================================
	subroutine random_std_normal(z1, z2)
		implicit none

		! Output variables
		real(dp), intent(out) :: z1, z2

		! Local variables
		real(dp) :: u1, u2

		! Generate two uniformly distributed random numbers in [0, 1)
		call random_std_uniform(u1)
		call random_std_uniform(u2)

		! Convert to standard normal distribution using Box-Muller transform
		z1 = dsqrt(-2*dlog(u1))*dcos(2*PI*u2)
		z2 = dsqrt(-2*dlog(u1))*dsin(2*PI*u2)

	end subroutine random_std_normal

	!=============================================================================
	! Subroutine: random_normal
	! Purpose   : Generate two independent random numbers following a
	!             general normal distribution with the specified
	!             parameters (mean = mu, standard deviation = sigma).
	! Arguments :
	!   - real(dp), intent(in) :: mu
	!       The mean of the normal distribution.
	!   - real(dp), intent(in) :: sigma
	!       The standard deviation of the normal distribution.
	!   - real(dp), intent(out) :: x1
	!       The first independent random number generated following the
	!       specified normal distribution.
	!   - real(dp), intent(out) :: x2
	!       The second independent random number generated following the
	!       specified normal distribution.
	! Adapted from :
	! https://masuday.github.io/fortran_tutorial/random.html
	!=============================================================================
	subroutine random_normal(mu, sigma, x1, x2)
		implicit none

		! Input/Output variables
		real(dp), intent(in) :: mu, sigma
		real(dp), intent(out) :: x1, x2

		! Local variables
		real(dp) :: z1, z2
		call random_std_normal(z1, z2)

		! Generate standard normal random numbers and scale to desired
		! mean and standard deviation
		x1 = mu + sigma*z1
		x2 = mu + sigma*z2

	end subroutine random_normal

	!=============================================================================
	! Subroutine: random_exponential
	! Purpose   : Generate a random number following an exponential
	!             distribution with a specified rate parameter (λ > 0).
	! Arguments :
	!   - real(dp), intent(in) :: lambda
	!       The rate parameter (λ) of the exponential distribution.
	!   - real(dp), intent(out) :: x
	!       A random number following the exponential distribution.
	! Adapted from :
	! https://www.eg.bucknell.edu/~xmeng/Course/CS6337/Note/master/node50.html
	!=============================================================================
	subroutine random_exponential(lambda, x)
		implicit none

		! Input/Output variables
		real(dp), intent(in) :: lambda
		real(dp), intent(out) :: x

		! Local variable
		real(dp) :: u

		! Generate a random number from an exponential distribution with
    ! parameter lambda
		call random_std_uniform(u)
		x = -dlog(u)/lambda

	end subroutine random_exponential

	!=============================================================================
	! VECTOR TRANSFORMATION SUBROUTINES
	!=============================================================================

	!=============================================================================
	! Subroutine: vector_translation
	! Purpose   : Translate a 3-dimensional vector by adding a translation
	!             vector to it.
	! Arguments :
	!   - real(dp), intent(in) :: translation_vector(3)
	!       The translation vector to be added.
	!   - real(dp), intent(inout) :: vector(3)
	!       The vector to be translated. On input, it contains the original
	!       vector. On output, it contains the translated vector.
	!=============================================================================
	subroutine vector_translation(translation_vector, vector)
		implicit none

		! Input/Output variables
		real(dp), intent(in) :: translation_vector(3)
		real(dp), intent(inout) :: vector(3)

		! Translate the vector by adding the translation_vector
		vector = vector + translation_vector

	end subroutine vector_translation

	!=============================================================================
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
	!=============================================================================
	subroutine rotation_about_x_axis(theta, vector)
		implicit none

		! Input/Output variables
		real(dp), intent(in) :: theta
		real(dp), intent(inout) :: vector(3)

		! Local variables
		real(dp) :: rotation_matrix(3,3), rotated_vector(3)
		integer :: i

		! Define the rotation matrix for a rotation about the x-axis by angle theta
		rotation_matrix(1,:) = (/1._dp, 0._dp, 0._dp/)
		rotation_matrix(2,:) = (/0._dp, dcos(theta), -dsin(theta)/)
		rotation_matrix(3,:) = (/0._dp, dsin(theta), dcos(theta)/)

		! Initialize rotated_vector to zero
		rotated_vector = 0

		! Initialize rotated_vector to zero
		do i = 1, 3
			rotated_vector(:) = rotated_vector(:) + rotation_matrix(:,i)*vector(i)
		end do

		! Update the input vector with the rotated result
		vector = rotated_vector

	end subroutine rotation_about_x_axis

	!=============================================================================
	! INPUT AND OUTPUT FILE HANDLING SUBROUTINES
	!=============================================================================

	!=============================================================================
	! Subroutine: read_input_parameters
	! Purpose   : Read input parameters from the 'input.txt' file.
	! Arguments :
	!   - integer(i8), intent(out) :: num_electrons
	!       Number of electrons in the beam for the simulation.
	!   - real(dp), intent(out) :: beam_energy
	!       Energy of the electron beam in kiloelectron volts (keV).
	!   - real(dp), intent(out) :: energy_spread
	!       Energy spread of the electron beam as a percentage (%).
	!   - real(dp), intent(out) :: spot_size_factor
	!       Scale factor determining the spot size of the electron beam.
	!   - real(dp), intent(out) :: beam_target_distance
	!       Distance from the electron beam source to the target material
	!       surface in angstroms (Å).
	!   - real(dp), intent(out) :: grazing_angle
	!       Grazing angle of the electron beam in degrees (º).
	!   - integer(i8), intent(out) :: material_boundaries(3)
	!       Indexes of the material model grid boundaries.
	!   - logical, intent(out) :: beam_model_output_saving_enabled
	!       Flag to enable saving beam model positions and velocities.
	!   - logical, intent(out) :: material_model_output_saving_enabled
	!       Flag to enable saving material model atom positions.
	!   - logical, intent(out) :: electron_trajectories_saving_enabled
	!       Flag to enable saving electron trajectories.
	!   - integer(i8), intent(out) :: num_plot_points
	!       Number of plot points for visualization.
	!   - real(dp), intent(out) :: dt
	!       Time step size in atomic units of time (aut).
	!=============================================================================
	subroutine read_input_parameters &
		(num_electrons, beam_energy, energy_spread, spot_size_factor, &
		beam_target_distance, grazing_angle, material_boundaries, &
		beam_model_output_saving_enabled, material_model_output_saving_enabled, &
		electron_trajectories_saving_enabled, num_plot_ploints, dt)
		implicit none

		! Output variables
		integer(i8), intent(out) :: num_electrons, material_boundaries(3)
		real(dp), intent(out) :: beam_energy, energy_spread, beam_target_distance
		real(dp), intent(out) :: spot_size_factor
		real(dp), intent(out) :: grazing_angle, dt
		logical, intent(out) :: beam_model_output_saving_enabled
		logical, intent(out) :: material_model_output_saving_enabled
		logical, intent(out) :: electron_trajectories_saving_enabled
		integer(i8), intent(out) :: num_plot_ploints

		! Open the input file
		open(unit=42, file='input.txt', status='old', action='read')

		! Read the simulation parameters. The extra read statement before
		! each actual parameter reading is there to omit the line with the
		! value description preceding every value in the 'input.txt' file
		read(42, *) ! Read beam energy in kiloelectronvolts (keV)
		read(42, *) beam_energy
		read(42, *) ! Read grazing angle in degrees (º)
		read(42, *) grazing_angle
		read(42, *) ! Read number of electrons in beam (integer)
		read(42, *) num_electrons
		read(42, *) ! Read beam energy spread as percentage (%)
		read(42, *) energy_spread
		read(42, *) ! Read spot size scale factor (real)
		read(42, *) spot_size_factor
		read(42, *) ! Read beam target distance in angstroms (Å)
		read(42, *) beam_target_distance
		read(42, *) ! Read material model grid index boundaries (integers)
		read(42, *) material_boundaries
		read(42, *) ! Read number of points to be plotted (integer)
		read(42, *) num_plot_ploints
		read(42, *) ! Read time step size in atomic units of time (aut)
		read(42, *) dt
		read(42, *) ! Read flag for beam model output saving (boolean)
		read(42, *) beam_model_output_saving_enabled
		read(42, *) ! Read flag for material model output saving (boolean)
		read(42, *) material_model_output_saving_enabled
		read(42, *) ! Read flag for electron trajectories saving (boolean)
		read(42, *) electron_trajectories_saving_enabled

		! Close the input file
		close(42)

	end subroutine read_input_parameters

	!=============================================================================
	! Subroutine: open_output_files
	! Purpose   : Open files for writing the following simulation outputs:
	!             electron trajectories (optionally), final positions of
	!             embedded and scattered electrons, scattering angles of
	!             scattered electrons, all final electron positions, time
	!             evolution of electrons' final state, and iteration
	!             and total simulation times.
	! Arguments :
	!   - integer, intent(in) :: output_unit
	!       The unit number to be used as a base for opening the output files.
	!   - logical, intent(in) :: electron_trajectories_saving_enabled
	!       Flag to determine if electron trajectories should be saved.
	!=============================================================================
	subroutine open_output_files &
		(output_unit, electron_trajectories_saving_enabled)
		implicit none

		! Input variables
		integer, intent(in) :: output_unit
		logical, intent(in) :: electron_trajectories_saving_enabled

		! Open file to save electron trajectories, if trajectory saving is enabled
		if (electron_trajectories_saving_enabled) then
			open(unit=output_unit+1, file="electron_trajectories.dat", &
				status='replace', action='write')
		end if

		! Open file to save the final position of embedded electrons
		open(unit=output_unit+2, file='embedded_positions.dat', &
			status='replace', action='write')

		! Open file to save the final position of scattered electrons
		open(unit=output_unit+3, file='scattered_positions.dat', &
			status='replace', action='write')

		! Open file to save the scattered electrons scattering angles
		open(unit=output_unit+4, file='scattering_angles.dat', &
			status='replace', action='write')

		! Open file to save all final electron positions
		open(unit=output_unit+5, file='final_electron_positions.dat', &
			status='replace', action='write')

		! Open file to save the electrons' final state evolution
		open(unit=output_unit+6, file='final_state_evolution.dat', &
			status='replace', action='write')

		! Open file to save each iteration and total simulation times
		open(unit=output_unit+7, file='simulation_times.dat', &
			status='replace', action='write')

	end subroutine open_output_files

	!=============================================================================
	! Subroutine: write_simulation_status_information
	! Purpose   : Write simulation status information up to the latest electron
	!             trajectory completed to a file named 'simulation_info.dat'.
	! Arguments :
	!   - integer, intent(in) :: output_unit
	!       The unit number used for writing the simulation status information.
	!   - integer(i8), intent(in) :: i
	!       The latest electron trajectory completed.
	!   - integer(i8), intent(in) :: num_electrons
	!       The total number of electron trajectories to be simulated.
	!   - integer(i8), intent(in) :: num_embedded
	!       Number of electrons embedded up to the latest trajectory completed.
	!   - integer(i8), intent(in) :: num_scattered
	!       Number of electrons scattered up to the latest trajectory completed.
	!   - real(dp) :: total_time
	!       Total simulation time up to the latest trajectory completed.
	!=============================================================================
	subroutine write_simulation_status_information &
		(output_unit, i, num_electrons, num_embedded, num_scattered, total_time)
		implicit none

		! Input variables
		integer, intent(in) :: output_unit
		integer(i8), intent(in) :: i, num_electrons
		integer(i8), intent(in) :: num_embedded, num_scattered
		real(dp), intent(in) :: total_time

		! Open file to save simulation status information
		open(unit=output_unit, file='simulation_info.dat', &
			status='replace', action='write')

		! Write number of electron trajectories completed so far
		write(output_unit, *) i, "out of", num_electrons, &
			"electron trajectories simulated"

		! Indicate completion if all electron trajectories are complete
		if (i .eq. num_electrons) write(output_unit, *) "SIMULATION COMPLETE!"

		! Write the number of electrons embedded and scattered so far
		write(output_unit, *) " Number of electrons embedded: ", num_embedded
		write(output_unit, *) " Number of electrons scattered:", num_scattered

		! Write simulation time based on progress so far
		if (i .ne. num_electrons) then
			write(output_unit, *) "Simulation time thus far:", total_time
		else
			write(output_unit, *) "Total simulation time:", total_time
		end if

		! Close the file
		close(output_unit)

	end subroutine write_simulation_status_information

end module m0_utilities