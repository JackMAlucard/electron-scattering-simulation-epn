module m1_electron_beam_model

	use m0_utilities, &
		only: dp, i8, PI, INTERATOMIC_DIST_SIO2, kev_to_atomic_energy_conversion, &
		angstrom_to_atomic_distance_conversion, random_normal, vector_translation, &
		rotation_about_x_axis

	implicit none
	contains
	!=============================================================================
	! Subroutine: electron_beam_parameters_unit_conversion
	! Purpose   : Convert electron beam parameters from their original units
	!             to atomic units for use in the simulation.
	! Arguments :
	!   - real(dp), intent(inout) :: beam_energy
	!       The energy of the electron beam. On input, it should be in
	!       kiloelectron volts (keV). On output, it will be converted to
	!       atomic units (Hartrees, Eh).
	!   - real(dp), intent(inout) :: beam_target_distance
	!       The distance from the electron beam source to the target
	!       material surface. On input, it should be in angstroms (Å).
	!       On output, it will be converted to atomic units (Bohr radius, a0).
	!   - real(dp), intent(inout) :: grazing_angle
	!       The grazing angle of the electron beam. On input, it should be
	!       in degrees. On output, it will be converted to radians.
	!=============================================================================
	subroutine electron_beam_parameters_unit_conversion &
		(beam_energy, beam_target_distance, grazing_angle)
		implicit none

		! Input/Output variables
		real(dp), intent(inout) :: beam_energy
		real(dp), intent(inout) :: beam_target_distance, grazing_angle

		! Convert units
		call kev_to_atomic_energy_conversion(beam_energy)
		call angstrom_to_atomic_distance_conversion(beam_target_distance)
		grazing_angle = grazing_angle*PI/180

	end subroutine electron_beam_parameters_unit_conversion

	!=============================================================================
	! Subroutine: set_up_electron_beam_model
	! Purpose   : Set up the electron beam model by initializing the
	!             positions, velocities, and accelerations for the electrons.
	!             - Positions are initially located on the xy-plane and follow
	!               a bivariate normal distribution symmetric around the origin.
	!               To achieve this, the x and y coordinates are independently
	!               drawn from a normal distribution with the same mean (mu = 0)
	!               and standard deviation (sigma), with the latter set in
	!               relation to the beam's spot size. In the model, the spot
	!               size is defined as the Full Width at Tenth Maximum (FWTM)
	!               of the distribution, spot size = 2*sqrt(2*ln(10))*sigma.
	!             - Velocities are initially directed along the +z axis, with
	!               their magnitudes determined by the kinetic energy
	!               distribution. The energy values follow a normal distribution
	!               with a mean equal to the primary beam energy and a standard
	!               deviation given by the energy spread.
	!             - The initial electron positions are translated along the
	!               negative direction of the z-axis. They are then rotated,
	!               along with the initial velocities, around the x-axis by the
	!               grazing angle. This transformation aligns the beam with the
	!               target material, positioning the electrons to simulate the
	!               beam's incidence on the target surface.
	!             - The positions and velocities are saved, before and after the
	!               transformations, if the output saving flag is enabled.
	! Arguments :
	!   - integer(i8), intent(in) :: num_electrons
	!       Number of electrons in the beam.
	!   - real(dp), intent(in) :: spot_size_factor
	!       Scale factor determining the spot size of the electron beam.
	!   - real(dp), intent(in) :: beam_energy
	!       Primary energy of the electron beam in atomic units (Hartrees, Eh).
	!   - real(dp), intent(in) :: energy_spread
	!       Energy spread of the electron beam as a percentage.
	!   - real(dp), intent(in) :: beam_target_distance
	!       Distance from the electron beam source to the target material
	!       surface in atomic units (Bohr radius, a0).
	!   - real(dp), intent(in) :: grazing_angle
	!       Grazing angle of the electron beam in radians.
	!   - logical, intent(in) :: output_saving_enabled
	!       Flag to determine if the beam model initial positions and
	!       velocities are saved.
	!   - real(dp), allocatable, intent(out) :: electron_positions(:,:)
	!       Array to store the initial positions of the electrons in atomic
	!       units (a0).
	!   - real(dp), allocatable, intent(out) :: electron_velocities(:,:)
	!       Array to store the initial velocities of the electrons in atomic
	!       units (a0/aut).
	!   - real(dp), allocatable, intent(out) :: electron_accelerations(:,:)
	!       Array to store the initial accelerations of the electrons in
	!       atomic units (a0/aut^2).
	!=============================================================================
  subroutine set_up_electron_beam_model &
		(num_electrons, spot_size_factor, beam_energy, energy_spread, &
		beam_target_distance, grazing_angle, output_saving_enabled, &
		electron_positions, electron_velocities, electron_accelerations)
		implicit none

		! Input/Output variables
		integer(i8), intent(in) :: num_electrons
		real(dp), intent(in) :: spot_size_factor
		real(dp), intent(in) :: beam_energy, energy_spread
		real(dp), intent(in) :: beam_target_distance, grazing_angle
		logical, intent(in) :: output_saving_enabled
		real(dp), allocatable, intent(out) :: electron_positions(:,:)
		real(dp), allocatable, intent(out) :: electron_velocities(:,:)
		real(dp), allocatable, intent(out) :: electron_accelerations(:,:)

		! Local variables
		real(dp) :: spot_size
		real(dp) :: positions_mu, positions_sigma, energy_mu, energy_sigma
		real(dp) :: x, y, z, vx, vy, vz
		real(dp) :: energy_1, energy_2
		real(dp) :: translation_vector(3)
		integer(i8) :: i

		! Allocate arrays for electron positions, velocities, and accelerations
		allocate(electron_positions(num_electrons,3))
		allocate(electron_velocities(num_electrons,3))
		allocate(electron_accelerations(num_electrons,3))

		! Set electron beam spot size
		spot_size = spot_size_factor*INTERATOMIC_DIST_SIO2

		! Set the mean and standard deviation for the normal distribution of
		! the initial positions
		positions_mu = 0
		positions_sigma = spot_size/(2*dsqrt(2*dlog(10._dp)))

		! Set the mean and standard deviation for the normal distribution of
		! kinetic energies
		energy_mu = beam_energy
		energy_sigma = beam_energy*(energy_spread/100)

		! Generate random positions and velocities for each electron
		do i = 1, num_electrons
			! Generate two independent random numbers following a normal distribution
			! for the x and y coordinates of the positions with z = 0
			call random_normal(positions_mu, positions_sigma, x, y)
			z = 0
			electron_positions(i,:) = (/x, y, z/)

			! Initialize velocities in the +z direction
			vx = 0
			vy = 0
			! Generate two random energy values every odd iteration;
			! compute velocity z-component using equation in atomic units
			if (mod(i,2) .ne. 0) then
				call random_normal(energy_mu, energy_sigma, energy_1, energy_2)
				vz = dsqrt(2*energy_1)
			else
				vz = dsqrt(2*energy_2)
			end if
			electron_velocities(i,:) = (/vx, vy, vz/)

			! Initialize accelerations to zero
			electron_accelerations(i,:) = 0
		end do

		! Save initial positions and velocities if output saving is enabled
		if (output_saving_enabled) then
			open(unit=42, file='beam_model.dat', status='replace', action='write')
			do i = 1, num_electrons
				write(42, *) electron_positions(i,:), electron_velocities(i,:)
			end do
			write(42, *)
			write(42, *)
		end if

		! Set translation vector for electron positions
		translation_vector = (/0._dp, 0._dp, -beam_target_distance/)

		! Translate and rotate electron positions, and rotate electron velocities
		do i = 1, num_electrons
			call vector_translation(translation_vector, electron_positions(i,:))
			call rotation_about_x_axis(grazing_angle, electron_positions(i,:))
			call rotation_about_x_axis(grazing_angle, electron_velocities(i,:))
		end do

		! Save transformed positions and velocities if output saving is enabled
		if (output_saving_enabled) then
			do i = 1, num_electrons
				write(42, *) electron_positions(i,:), electron_velocities(i,:)
			end do
			close(42)
		end if

	end subroutine set_up_electron_beam_model

end module m1_electron_beam_model