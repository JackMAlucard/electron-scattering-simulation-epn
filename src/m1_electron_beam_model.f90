module m1_electron_beam_model
	implicit none
	use m0!, only: sp, dp, i8,...
	contains
!INDENTING EVERYTHING INSIDE HERE PENDING...

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
  subroutine setup_electron_beam(num_electrons, spot_size, beam_energy &
                                 beam_spread, grazing_angle, distance_to_target &
                                 electron_positions, electron_velocities &
                                 electron_accelerations)
! subroutine setup_electron_beam_model(...)                                  
!subroutine electron_beam(N, E, R, eb_pos, eb_vel, eb_acc)
!subroutine electron_beam(N, E, R, eb_r, eb_v, eb_a)
!subroutine electron_beam(N, E, R, r_eb, v_eb, a_eb)
!real(dp), allocatable, intent(out) :: eb_pos(:,:), eb_vel(:,:), eb_acc(:,:)
integer(i8), intent(in) :: num_electrons
real(dp), intent(in) :: spot_size, beam_energy, beam_spread
real(dp), intent(in) :: grazing_angle, distance_to_target
real(dp), allocatable, intent(out) :: electron_positions(:,:)
real(dp), allocatable, intent(out) :: electron_velocities(:,:)
real(dp), allocatable, intent(out) :: electron_accelerations(:,:)
real(dp) :: beam_mu, beam_sigma
real(dp) :: x, y, z, vx, vy, vz
real(dp) :: Ei1, Ei2
integer(i8) :: i

	!Unit conversion
	call kev_to_atomic_energy_conversion(beam_energy)
	call angstrom_to_atomic_distance_conversion(spot_size)

	allocate(electron_positions(num_electrons,3))
	allocate(electron_velocities(num_electrons,3))
	allocate(electron_accelerations(num_electrons,3))
	
	!They follow the normal distribution with
	!Full Width at Tenth of Maximum (FWTM):
	!FWTM = 2*sqrt(2*ln(10))*sigma ~ 4.29193*sigma
	beam_mu = 0
	beam_sigma = spot_size/(2*dsqrt(2*dlog(10._dp)))

	do i = 1, N
		!Generate positions
		!x and y coordinates
		call random_normal(beam_mu, beam_sigma, x, y)
		!z coordinate
		z = 0
		electron_positions(i,:) = (/x, y, z/)

		!Generate velocities
		!For 10kev, the speed is ~ 0.19783585 c = 59350755.095 m/s
		!Or 58455098 m/s = 0.194850327 c (with Lorentz tranform?)
		vx = 0
		vy = 0

		!Since two numbers are generated in each call of random_normal
		!Using this if-else structure I avoid generating more numbers than necessary
		if (mod(i+1,2) .eq. 0) then
			! Electron energy values follow a normal distribution 
			! with a mean of E, and a mean deviation of 0.1*E
			call random_normal(beam_energy,(energy_spread/100)*beam_energy
			                   electron_energy_1, electron_energy_2)
			vz = dsqrt(2*electron_energy_1)
		else
			vz = dsqrt(2*electron_energy_2)
		end if
		electron_velocities(i,:) = (/vx, vy, vz/)

		!Initialize accelerations
		electron_accelerations(i,:) = 0
		
		!Transforming electron beam
		
		!Translation vector
		T = (/0._dp, 0._dp, -distance_to_target/)
		
		!Unit conversion
		grazing_angle = grazing_angle*PI/180
		call angstrom_to_atomic_distance_conversion(distance_to_target)
		
		do i = 1, num_electrons
			!Rotation around x-axis
			call rotation_about_x_axis(grazing_angle, electron_positions(i,:))
			call rotation_about_x_axis(grazing_angle, electron_velocities(i,:))

			!Initial positions translation
			electron_positions(i,:) = electron_positions(i,:) + T
		end do
		
	end do

end subroutine setup_electron_beam

end module m1_electron_beam_model