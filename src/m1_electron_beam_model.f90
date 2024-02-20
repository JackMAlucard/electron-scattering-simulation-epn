module m1_electron_beam_model
	implicit none
	use m0!, only: sp, dp, i8,...
	contains

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
	! subroutine electron_beam(e_pos, e_vel, e_acc, N_in, E_in, R_in, prob_dist_in)
	!subroutine electron_beam(N, E, R, eb_pos, eb_vel, eb_acc)
	!subroutine electron_beam(N, E, R, eb_r, eb_v, eb_a)
	!subroutine electron_beam(N, E, R, r_eb, v_eb, a_eb)
!  subroutine setup_electron_beam(num_electrons, spot_size, &
!                                 beam_energy, beam_spread, &
!                                 distance_to_target, grazing_angle, &
!                                 electron_positions, electron_velocities, &
!                                 electron_accelerations)
  subroutine setup_electron_beam_model(num_electrons, spot_size, &
                                       beam_energy, beam_spread, &
                                       distance_to_target, grazing_angle, &
                                       electron_positions, &
                                       electron_velocities, &
                                       electron_accelerations)
		implicit none
		integer(i8), intent(in) :: num_electrons
		real(dp), intent(in) :: spot_size, beam_energy, beam_spread
		real(dp), intent(in) :: distance_to_target, grazing_angle
		real(dp), allocatable, intent(out) :: electron_positions(:,:)
		real(dp), allocatable, intent(out) :: electron_velocities(:,:)
		real(dp), allocatable, intent(out) :: electron_accelerations(:,:)
		real(dp) :: positions_mu, positions_sigma, energy_mu, energy_sigma
		real(dp) :: x, y, z, vx, vy, vz
		real(dp) :: energy_1, energy_2
		real(dp) :: traslation_vector(3)
		integer(i8) :: i

		! Unit conversion
		call angstrom_to_atomic_distance_conversion(spot_size)
		call kev_to_atomic_energy_conversion(beam_energy)
		call angstrom_to_atomic_distance_conversion(distance_to_target)
		grazing_angle = grazing_angle*PI/180

		allocate(electron_positions(num_electrons,3))
		allocate(electron_velocities(num_electrons,3))
		allocate(electron_accelerations(num_electrons,3))
		
		! Electron positions follow a normal distribution with
		! Full Width at Tenth of Maximum (FWTM):
		! FWTM = 2*sqrt(2*ln(10))*sigma ~ 4.29193*sigma
		positions_mu = 0
		positions_sigma = spot_size/(2*dsqrt(2*dlog(10._dp)))
		
		! Electron energies follow a normal distribution with
		energy_mu = beam_energy
		energy_sigma = (energy_spread/100)*beam_energy

		do i = 1, num_electrons
			! Generating positions
			! x and y coordinates
			call random_normal(beam_mu, beam_sigma, x, y)
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
		translation_vector = (/0._dp, 0._dp, -distance_to_target/)
		
		do i = 1, num_electrons
			! Initial positions translation
			call vector_translation(traslation_vector, electron_positions(i,:))
!				electron_positions(i,:) = electron_positions(i,:) + T
			
			! Rotation around x-axis
			call rotation_about_x_axis(grazing_angle, electron_positions(i,:))
			call rotation_about_x_axis(grazing_angle, electron_velocities(i,:))
			
		end do

	end subroutine setup_electron_beam_model

end module m1_electron_beam_model