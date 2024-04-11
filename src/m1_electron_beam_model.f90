module m1_electron_beam_model
	use m0_utilities, &
		only: dp, i8, PI, INTERATOMIC_DIST_SIO2, kev_to_atomic_energy_conversion, &
		angstrom_to_atomic_distance_conversion, random_normal, vector_translation, &
		rotation_about_x_axis
	implicit none
	contains
	
	subroutine electron_beam_parameters_unit_conversion &
		(beam_energy, beam_target_distance, grazing_angle)
		implicit none
		real(dp), intent(inout) :: beam_energy
		real(dp), intent(inout) :: beam_target_distance, grazing_angle
		call kev_to_atomic_energy_conversion(beam_energy)
		call angstrom_to_atomic_distance_conversion(beam_target_distance)
		grazing_angle = grazing_angle*PI/180
	end subroutine electron_beam_parameters_unit_conversion

	! Subroutine which generates an array with N electron random starting positions
	! following a normal distribution whose Full Width at Tenth Maximum (FWTM) 
	! represents the spot size (beam diameter) of the beam.
	! The electrons are set up to have initial velocities pointing in 
	! the +z direction with values following a normal distribution of mean equal 
	! to the beam primary energy. The initial accelerations are null.
	! The electron beam is set up to have a grazing incidence angle psi and
	! located a distance d from the target dielectric material.
	! The subroutine parameters are the number of electrons num_electrons, 
	! the beam primary energy beam_energy (in keV), the beam spread (as a % of
	! the beam primary energy, the spot size (beam diameter) spot_size (in Å), 
	! the grazing incidence angle grazing_angle (in degrees), and distance 
	! to the target material distance_to_target (in Å).
	! The subroutine converts the magnitudes to the atomic system of units.
	! The vectors generated are in the beam reference system, the 'primed' system.
  subroutine setup_electron_beam_model &
		(num_electrons, spot_size_factor, beam_energy, energy_spread, &
		beam_target_distance, grazing_angle, electron_positions, &
		electron_velocities, electron_accelerations)
		implicit none
		integer(i8), intent(in) :: num_electrons
		integer, intent(in) :: spot_size_factor
		real(dp), intent(in) :: beam_energy, energy_spread
		real(dp), intent(in) :: beam_target_distance, grazing_angle
		real(dp), allocatable, intent(out) :: electron_positions(:,:)
		real(dp), allocatable, intent(out) :: electron_velocities(:,:)
		real(dp), allocatable, intent(out) :: electron_accelerations(:,:)
		real(dp) :: spot_size
		real(dp) :: positions_mu, positions_sigma, energy_mu, energy_sigma
		real(dp) :: x, y, z, vx, vy, vz
		real(dp) :: energy_1, energy_2
		real(dp) :: translation_vector(3)
		integer(i8) :: i
		! Allocating electron vectors
		allocate(electron_positions(num_electrons,3))
		allocate(electron_velocities(num_electrons,3))
		allocate(electron_accelerations(num_electrons,3))
		! Defining spot size
		spot_size = spot_size_factor*INTERATOMIC_DIST_SIO2
		! Electron positions follow a normal distribution with
		! Full Width at Tenth of Maximum (FWTM):
		! FWTM = 2*sqrt(2*ln(10))*sigma ~ 4.29193*sigma
		positions_mu = 0
		positions_sigma = spot_size/(2*dsqrt(2*dlog(10._dp)))
		! Electron energies follow a normal distribution with
		energy_mu = beam_energy
		energy_sigma = beam_energy*(energy_spread/100)
		! Generating electron arrays
		do i = 1, num_electrons
			! Generating positions
			! x and y coordinates
			call random_normal(positions_mu, positions_sigma, x, y)
			! z coordinate
			z = 0
			electron_positions(i,:) = (/x, y, z/)
			! Generating velocities
			vx = 0
			vy = 0
			! Since two numbers are generated in each call of random_normal
			! Using this if-else structure I avoid generating unnecessary numbers
			if (mod(i+1,2) .eq. 0) then
				call random_normal(energy_mu, energy_sigma, energy_1, energy_2)
				vz = dsqrt(2*energy_1)
			else
				vz = dsqrt(2*energy_2)
			end if
			electron_velocities(i,:) = (/vx, vy, vz/)
			! Initializing accelerations
			electron_accelerations(i,:) = 0
		end do
		! Transforming electron beam vectors
		translation_vector = (/0._dp, 0._dp, -beam_target_distance/)
		do i = 1, num_electrons
			! Initial positions translation
			call vector_translation(translation_vector, electron_positions(i,:))
!				electron_positions(i,:) = electron_positions(i,:) + T
			! Rotation around x-axis
			call rotation_about_x_axis(grazing_angle, electron_positions(i,:))
			call rotation_about_x_axis(grazing_angle, electron_velocities(i,:))
		end do
	end subroutine setup_electron_beam_model

end module m1_electron_beam_model