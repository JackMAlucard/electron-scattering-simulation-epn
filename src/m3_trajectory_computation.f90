module m3_trajectory_computation

	use m0_utilities, &
		only: dp, i8, PI, INTERATOMIC_DIST_SIO2_INV, MATERIAL_HEIGHT_SIO2, &
		CROSS_SECTION_SIO2, random_exponential

	implicit none
	contains
	!=============================================================================
	! Subroutine: number_of_iterations_estimation
	! Purpose   : Estimate an upper bound to the number of iterations required to
	!             simulate an electron trajectory given initial conditions and
	!             time step size. Adjust the number of points to be plotted
	!             accordingly if needed.
	!             - An electron trajectory is considered to be complete when
	!               the electron has traversed at least twice its initial
	!               distance to the target material.
	!             - The final time used in the upper bound estimation is seven
	!               times (a loose overestimation) the time it would take an
	!               electron to traverse twice the distance to the target
	!               material if it were not subjected to any acceleration during
	!               the trajectory.
	!             - The subroutine also adjusts the number of points to be
	!               plotted to not exceed the upper bound to the number of
	!               iterations set by the estimation if necessary.
	! Arguments :
	!   - real(dp), intent(in) :: r(3)
	!       Initial position vector of the electron in atomic units (a0).
	!   - real(dp), intent(in) :: v(3)
	!       Initial velocity vector of the electron in atomic units (a0/aut).
	!   - real(dp), intent(in) :: dt
	!       Time step size in atomic units of time (aut).
	!   - integer(i8), intent(out) :: max_iterations
	!       Estimated upper bound to the number of iterations required to
	!       simulate a complete electron trajectory for the simulation.
	!   - integer(i8), intent(inout) :: num_plot_ploints
	!       Number of points to be plotted for a complete electron trajectory.
	!=============================================================================
	subroutine number_of_iterations_estimation &
		(r, v, dt, max_iterations, num_plot_ploints)
		implicit none

		! Input/Output variables
		real(dp), intent(in) :: r(3), v(3), dt
		integer(i8), intent(out) :: max_iterations
		integer(i8), intent(inout) :: num_plot_ploints

		! Local variables
		real(dp) :: d0, v0, a0, tf

		! Compute electron's initial distance to target and initial speed
		d0 = norm2(r)
		v0 = norm2(v)

		! Set final time for the trajectory simulation as a loose overestimation
		tf = 7*(2*d0/v0)

		! Calculate maximum number of iterations based on final time and time step
		max_iterations = dint(tf/dt)

		! Adjust number of plot points if it exceeds maximum number of iterations
		if (max_iterations .lt. num_plot_ploints) num_plot_ploints = max_iterations

	end subroutine number_of_iterations_estimation

	!=============================================================================
	! Subroutine: get_nearest_atom_indices
	! Purpose   : Determine the indices of the nearest atom in the silica (SiO2)
	!             material grid relative to the position of a projectile electron,
	!             considering the material's boundary limits.
	! Arguments :
	!   - real(dp), intent(in) :: r(3)
	!       Position vector of the projectile electron in atomic units (a0).
	!   - integer(i8), intent(in) :: material_boundaries(3)
	!       Array containing each material grid index boundary.
	!   - integer(i8), intent(out) :: atom_indices(3)
	!       Array containing the indices of the nearest atom in the material grid.
	!=============================================================================
	subroutine get_nearest_atom_indices(r, material_boundaries, atom_indices)
		implicit none

		! Input/Output variables
		real(dp), intent(in) :: r(3)
		integer(i8), intent(in) :: material_boundaries(3)
		integer(i8), intent(out) :: atom_indices(3)

		! Local variables
		integer(i8) :: mbi, mbj, mbk, i, j, k

		! Extract material grid index boundaries
		mbi = material_boundaries(1)
		mbj = material_boundaries(2)
		mbk = material_boundaries(3)

		! Calculate the index in the x-direction and constrain within boundaries
		i = idnint(r(1)*INTERATOMIC_DIST_SIO2_INV)
		if (i .gt. mbi) i = mbi
		if (i .lt. -mbi) i = -mbi

		! Calculate the index in the y-direction and constrain within boundaries
		j = idnint(r(2)*INTERATOMIC_DIST_SIO2_INV)
		if (j .gt. 0) j = 0
		if (j .lt. -mbj) j = -mbj

		! Calculate the index in the z-direction and constrain within boundaries
		k = idnint(r(3)*INTERATOMIC_DIST_SIO2_INV)
		if (k .gt. mbk) k = mbk
		if (k .lt. -mbk) k = -mbk

		! Store the nearest atom indices
		atom_indices = (/i, j, k/)

	end subroutine get_nearest_atom_indices

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
		real(dp), intent(in) :: rt(3) ! Position of the target electron (a0)
		real(dp), intent(out) :: a(3) ! Resulting acceleration vector (a0/aut^2)

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
	! Subroutine: acceleration_due_to_atom
	! Purpose   : Calculate the acceleration vector experienced by a projectile
	!             electron due to the interaction with a stationary neutral target
	!             atom through the electrostatic Thomas-Fermi-Molière (TFM)
	!             potential based on their positions and the atom's atomic number.
	! Arguments :
	!   - real(dp), intent(in) :: rp(3)
	!       Position vector of the projectile electron in atomic units (a0).
	!   - real(dp), intent(in) :: rt(3)
	!       Position vector of the target atom in atomic units (a0/aut).
	!   - integer, intent(in) :: Z
	!       Atomic number of the stationary neutral atom.
	!   - real(dp), intent(in) :: cbrt_Z
	!       Cubic root of the atomic number of the stationary neutral atom.
	!   - real(dp), intent(out) :: a(3)
	!       Acceleration vector experienced by the projectile electron due to 
	!       the stationary target atom in atomic units (a0/aut^2).
	!=============================================================================
	subroutine acceleration_due_to_atom(rp, rt, Z, cbrt_Z, a)
		implicit none

		! Input/Output variables
		real(dp), intent(in) :: rp(3)		! Position of the projectile electron (a0)
    real(dp), intent(in) :: rt(3)		! Position of the target atom (a0)
    integer, intent(in) :: Z				! Atomic number of the target atom
    real(dp), intent(in) :: cbrt_Z	! Cubic root of the atomic number
    real(dp), intent(out) :: a(3)		! Resulting acceleration vector (a0/aut^2)
 
		! Local variables
		! Analytic potential coefficients
		real(dp), parameter :: a_coef(3) = (/0.10_dp, 0.55_dp, 0.35_dp/)
		real(dp), parameter :: b_coef(3) = (/6.00_dp, 1.20_dp, 0.30_dp/)
		! Inverse of the b0 coefficient
		real(dp), parameter :: b0_inv = ((2**7)/(3*PI)**2)**(1/3._dp)
		real(dp) :: rs(3)			! Separation vector between the electron and atom (a0)
		real(dp) :: r					! Magnitude of the separation vector (a0)
		real(dp) :: chi, psi	! Value of TFM screening function and its derivative
		real(dp) :: aux
		integer(i8) :: i

		! Compute the separation vector between electron and atom and its magnitude
		rs = rp - rt
		r = norm2(rs)

		! Initialize the TFM screening function (chi) and its derivative (psi)
		chi = 0
		psi = 0

		! Sum over the exponentials to compute chi and psi
		do i = 1, 3
			aux = a_coef(i)*dexp(-b_coef(i)*r*b0_inv*cbrt_Z)
			chi = chi + aux
			aux = aux*b_coef(i)
			psi = psi + aux
		end do

		! Compute the acceleration using the TFM potential formula in atomic units
		a = -(Z/r**3)*(chi + psi*b0_inv*cbrt_Z*r)*rs
	end subroutine acceleration_due_to_atom

	!=============================================================================
	! Subroutine: time_step
	! Purpose   : Perform a time step update for a projectile electron, including
	!             position, velocity, and acceleration updates, using the
	!             Velocity Verlet algorithm to simulate its trajectory.
	!             - The algorithm considers interactions with all embedded
	!               electrons, the nearest atom and its closest neighbors.
	!             - Acceleration due to embedded electrons is computed for every
	!               embedded electron present up to the current trajectory.
	!             - Acceleration due to neutral atoms is calculated only if the
	!               electron is within the material zone, taken to be below the
	!               height of the material's top layer.
	! Arguments :
	!   - integer(i8), intent(in) :: num_embedded
	!       Number of embedded electrons.
	!   - integer(i8), intent(in) :: material_boundaries(3)
	!       Material grid index boundaries.
	!   - real(dp), intent(in) :: embedded_positions(:,:)
	!       Positions of embedded electrons in atomic units (a0).
	!   - real(dp), intent(in) :: atom_positions(-material_boundaries(1)-1:, &
	!       -material_boundaries(2)-1:, -material_boundaries(3)-1:, :)
	!       Positions of material atoms in atomic units (a0).
	!   - integer, intent(in) :: atomic_numbers(-material_boundaries(1)-1:, &
	!       -material_boundaries(2)-1:, -material_boundaries(3)-1:)
	!       Atomic numbers of material atoms.
	!   - real(dp), intent(in) :: atomic_numbers_cbrt &
	!       (-material_boundaries(1)-1:, -material_boundaries(2)-1:, &
	!       -material_boundaries(3)-1:)
	!       Cubic roots of the atomic numbers of material atoms.
	!   - real(dp), intent(in) :: dt
	!       Time step size for the simulation in atomic units of time (aut)
	!   - real(dp), intent(inout) :: r(3)
	!       Position vector of the electron in atomic units (a0). Updated after
	!       the time step.
	!   - real(dp), intent(inout) :: v(3)
	!       Velocity vector of the electron in atomic units (a0/aut). Updated
	!       after the time step.
	!   - real(dp), intent(inout) :: a(3)
	!       Acceleration vector of the electron in atomic units (a0/aut^2).
	!       Updated after the time step.
	!=============================================================================
	subroutine time_step &
		(num_embedded, material_boundaries, embedded_positions, atom_positions, &
		atomic_numbers, atomic_numbers_cbrt, dt, r, v, a)
		implicit none

		! Input/Output variables
		integer(i8), intent(in) :: num_embedded, material_boundaries(3)
		real(dp), intent(in) :: embedded_positions(:,:)
		real(dp), intent(in) :: atom_positions(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:, :)
		integer, intent(in) :: atomic_numbers(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:)
		real(dp), intent(in) :: atomic_numbers_cbrt(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:)
		real(dp), intent(in) :: dt
		real(dp), intent(inout) :: r(3), v(3), a(3)

		! Local variables
		real(dp) :: rt(3), ak(3), cbrt_Z
		integer :: Z
		integer(i8) :: atom_indices(3)
		integer(i8) :: i, j, k

		! Half-step velocity update
		v = v + 0.5*a*dt

		! Position update
		r = r + v*dt

		! Full-step acceleration update
		! Reset acceleration for the current time step
		a = 0

		! Acceleration due to embedded electrons;
		! only computed if there is at least one embedded electron
		if (num_embedded .gt. 0) then
			do k = 1, num_embedded
				rt = embedded_positions(k,:)
				call acceleration_due_to_electron(r, rt, ak)
				a = a + ak
			end do
		end if

		! Acceleration due to nearest material atom and its closest neighbors
		if (r(2) .lt. MATERIAL_HEIGHT_SIO2) then
			call get_nearest_atom_indices(r, material_boundaries, atom_indices)
			do i = atom_indices(1) - 1, atom_indices(1) + 1
				do j = atom_indices(2) - 1, atom_indices(2) + 1
					do k = atom_indices(3) - 1, atom_indices(3) + 1
						! Only compute acceleration in positions with atoms
						Z = atomic_numbers(i,j,k)
						if (Z .ne. 0) then
							rt = atom_positions(i,j,k,:)
							cbrt_Z = atomic_numbers_cbrt(i,j,k)
							call acceleration_due_to_atom(r, rt, Z, cbrt_Z, ak)
							a = a + ak
						end if
					end do
				end do
			end do
		end if

		! Second half-step velocity update
		v = v + 0.5*a*dt

	end subroutine time_step

	!=============================================================================
	! Subroutine: compute_trajectory
	! Purpose   : Simulate the trajectory of a projectile electron aimed at the
	!             dielectric material considering interactions with embedded
	!             electrons and material atoms. The trajectory is computed using
	!             the Velocity Verlet algorithm. When the trajectory simulation
	!             is completed, the projectile electron can either get embedded
	!             in the material or be scattered.
	!             - The trajectory simulation ends when one of the following
	!               conditions is met:
	!               - The electron is inside the material and a collision with
	!                 an atomic electron occurs (embedded),
	!               - The electron is outside the material zone and its distance
	!                 to the origin exceeds its initial distance (scattered),
	!               - The electron is outside the material zone and the maximum
	!                 number of iterations is reached (scattered).
	!             - A collision with an atomic electron is simulated to occur
	!               when the distance the projectile electron travels inside
	!               the material exceeds a random distance drawn from an
	!               exponential distribution with parameter equal to the
	!               macroscopic cross section of the material.
	!             - The final positions of embedded and scattered electrons are
	!               stored in the respective arrays.
	!             - A specified number of trajectory points are saved to a file,
	!               if the electron trajectories saving flag is enabled.
	! Arguments :
	!   - logical, intent(in) :: electron_trajectories_saving_enabled
	!       Flag to determine if the electron trajectories are saved.
	!   - integer(i8), intent(in) :: num_plot_ploints
	!       Maximum number of points to be plotted.
	!   - integer(i8), intent(in) :: max_iterations
	!       Maximum number of iterations allowed for the trajectory computation.
	!   - integer, intent(in) :: output_unit
	!       Unit number for the output file.
	!   - integer(i8), intent(in) :: material_boundaries(3)
	!       Material grid index boundaries.
	!   - real(dp), intent(in) :: atom_positions(-material_boundaries(1)-1:, &
	!       -material_boundaries(2)-1:, -material_boundaries(3)-1:, :)
	!       Positions of material atoms in atomic units (a0).
	!   - integer, intent(in) :: atomic_numbers(-material_boundaries(1)-1:, &
	!       -material_boundaries(2)-1:, -material_boundaries(3)-1:)
	!       Atomic numbers of material atoms.
	!   - real(dp), intent(in) :: atomic_numbers_cbrt &
	!       (-material_boundaries(1)-1:, -material_boundaries(2)-1:, &
	!       -material_boundaries(3)-1:)
	!       Cubic roots of atomic numbers of material atoms.
	!   - real(dp), intent(in) :: dt
	!       Time step size for the simulation in atomic units of time (aut)
	!   - real(dp), intent(inout) :: r(3)
	!       Position vector of the electron in atomic units (a0). Updated during
	!       the simulation.
	!   - real(dp), intent(inout) :: v(3)
	!       Velocity vector of the electron in atomic units (a0/aut). Updated
	!       during the simulation.
	!   - real(dp), intent(inout) :: a(3)
	!       Acceleration vector of the electron in atomic units (a0/aut^2).
	!       Updated during the simulation.
	!   - integer(i8), intent(inout) :: num_embedded
	!       Number of embedded electrons. Updated at the end of the trajectory
	!       simulation if the electron gets embedded.
	!   - integer(i8), intent(inout) :: num_scattered
	!       Number of scattered electrons. Updated at the end of the trajectory
	!       simulation if the electron is scattered.
	!   - real(dp), intent(inout) :: embedded_positions(:,:)
	!       Positions of embedded electrons in atomic units (a0). Updated at the
	!       end of the trajectory simulation if the electron gets embedded.
	!   - real(dp), intent(inout) :: scattered_positions(:,:)
	!       Positions of scattered electrons in atomic units (a0). Updated at the
	!       end of the trajectory simulation if the electron is scattered.
	!=============================================================================
	subroutine compute_trajectory &
		(electron_trajectories_saving_enabled, num_plot_ploints, max_iterations, &
		output_unit, material_boundaries, atom_positions, atomic_numbers, &
		atomic_numbers_cbrt, dt, r, v, a, num_embedded, num_scattered, &
		embedded_positions, scattered_positions)
		implicit none

		! Input/Output variables
		logical, intent(in) :: electron_trajectories_saving_enabled
		integer(i8), intent(in) :: num_plot_ploints, max_iterations
		integer, intent(in) :: output_unit
		integer(i8), intent(in) :: material_boundaries(3)
		real(dp), intent(in) :: atom_positions(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:, :)
		integer, intent(in) :: atomic_numbers(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:)
		real(dp), intent(in) :: atomic_numbers_cbrt(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:)
		real(dp), intent(in) :: dt
		real(dp), intent(inout) :: r(3), v(3), a(3)
		integer(i8), intent(inout) :: num_embedded, num_scattered
		real(dp), intent(inout) :: embedded_positions(:,:), scattered_positions(:,:)

		! Local variables
		logical :: is_embedded, is_scattered, is_max_iteration
		logical :: in_material
		real(dp) :: initial_distance_to_target
		real(dp) :: distance_before_collision, distance_in_material
		real(dp) :: previous_position(3), actual_position(3), step_length
		real(dp) :: distance_to_target, t
		integer(i8) :: mbi, mbj, mbk
		integer(i8) :: i, j

		! Initialize end condition flags
		is_embedded = .false.
		is_scattered = .false.
		is_max_iteration = .false.

		! Initialize embedding end condition variable: inside material flag
		in_material = .false.

		! Initialize scattering end condition variables: distances to target
		initial_distance_to_target = norm2(r)
		distance_to_target = initial_distance_to_target

		! Intialize maximum iterations end condition variable: iteration counter
		i = 0

		! Initialize plot points counter
		j = 0
		
		! Main simulation loop
		do while (.not.(is_embedded .or. is_scattered .or. is_max_iteration))
			! Current time
			t = i*dt

			! Plot trajectory points, if conditions are met
			if (electron_trajectories_saving_enabled .and. &
				(mod(i,max_iterations/num_plot_ploints) .eq. 0) .and. &
				j .lt. num_plot_ploints) then
				! Write time and position to output file
				write(output_unit,*) t, r
				j = j + 1
			end if

			! Compute next time step using Velocity Verlet algorithm:
			! Update electron position r, velocity v, and acceleration a
			call time_step &
			(num_embedded, material_boundaries, embedded_positions, &
			atom_positions, atomic_numbers, atomic_numbers_cbrt, dt, r, v, a)

			! Update distance to target and iteration count
			distance_to_target = norm2(r)
			i = i + 1

			! Check for embedding end condition: electron is inside the material and
			! a collision with an atomic electron occurs
			if (r(2) .lt. MATERIAL_HEIGHT_SIO2) then
				! Only on first entry into material zone
				if (.not.(in_material)) then
					! Update inside material flag
					in_material = .true.

					! Generate travel distance before collision with atomic electron
					call random_exponential(CROSS_SECTION_SIO2, distance_before_collision)

					! Initialize distance traveled within material
					actual_position = r
					distance_in_material = 0
				end if

				! Update position and track distance traveled within material
				previous_position = actual_position
				actual_position = r
				step_length = norm2(actual_position - previous_position)
				distance_in_material = distance_in_material + step_length

				! Check if electron gets embedded due to collision with atomic electron
				if (distance_in_material .ge. distance_before_collision) then
					! Update embedded condition flag, embedded number and positions
					is_embedded = .true.
					num_embedded = num_embedded + 1
					embedded_positions(num_embedded,:) = r

					! Print end condition information to console
					print*, "Trajectory end --> Electron is embedded"
					print*, "Total iterations:", i
					print*, "Final electron position:", r
					print*, " Distance before collision:", distance_before_collision
					print*, " Distance traveled inside material:", distance_in_material
					print*, " Material boundaries:"
					mbi = material_boundaries(1)
					mbj = material_boundaries(2)
					mbk = material_boundaries(3)
					print*, "  x --> +/-", atom_positions(mbi,-mbj,mbk,1)
					print*, "  y -->    ", atom_positions(mbi,-mbj,mbk,2)
					print*, "  z --> +/-", atom_positions(mbi,-mbj,mbk,3)
					print*
				end if
			end if

			! Check for scattering end condition: electron is outside the material
			! and exceeds its initial distance to target
			if (distance_to_target .gt. initial_distance_to_target .and. &
			r(2) .gt. 0) then
				! Update scattered condition flag, scattered number and positions
				is_scattered = .true.
				num_scattered = num_scattered + 1
				scattered_positions(num_scattered,:) = r

				! Print end condition information to console
				print*, "Trajectory end --> Electron is scattered"
				print*, "Total iterations:", i
				print*, "Final electron position:", r
				print*, " Initial distance to target:", initial_distance_to_target
				print*, " Final distance to target:", distance_to_target
				print*
			end if

			! Check for scattering end condition: electron is outside the material
			! and maximum number of iterations is exceeded
			if (i .ge. max_iterations) then
				! Update maximum iterations and scattered conditions flag, scattered
				! number and positions
				is_max_iteration = .true.
				is_scattered = .true.
				num_scattered = num_scattered + 1
				scattered_positions(num_scattered,:) = r

				! Print end condition information to console
				print*, "Trajectory end --> Electron is scattered"
				print*, "Total iterations:", i
				print*, "Final electron position:", r
				print*, " Maximum number of iterations reached!"
				print*, " Maximum number of iterations:", max_iterations
				print*
			end if

		end do

	end subroutine compute_trajectory

	!=============================================================================
	! Subroutine: compute_scattering_angles
	! Purpose   : Compute the scattering angles alpha and beta of the last
	!             scattered electron based on its final position.
	!             - The alpha angle is an elevation angle, it measures the
	!               angle between the projection of the position vector (x, y, z)
	!               onto the xz-plane and the vector itself:
	!               alpha = arctan(y/sqrt(x^2 + z^2)) = atan2(y, sqrt(x^2 + z^2)).
	!             - The beta angle is an azimuthal angle, it measures the angle
	!               between the +z direction and the projection of the position
	!               vector (x, y, z) onto the xz-plane:
	!               beta = arctan(-x/z) = -arctan(x/z) = -atan2(x, z).
	! Arguments :
	!   - integer(i8), intent(in) :: num_scattered
	!       Number of scattered electrons.
	!   - real(dp), intent(in) :: scattered_positions(:,:)
	!       Array containing the positions of scattered electrons in atomic
	!       units (a0).
	!   - real(dp), intent(out) :: alpha
	!       Scattering angle alpha, elevation angle in degrees (º).
	!   - real(dp), intent(out) :: beta
	!       Scattering angle beta, azimuthal angle in degrees (º).
	!=============================================================================
	subroutine compute_scattering_angles &
		(num_scattered, scattered_positions, alpha, beta)
		implicit none

		! Input/Output variables
		integer(i8), intent(in) :: num_scattered
		real(dp), intent(in) :: scattered_positions(:,:)
		real(dp), intent(out) :: alpha, beta

		! Local variables
		real(dp) :: x, y, z, s

		! Extract coordinates of the last scattered electron
		x = scattered_positions(num_scattered,1)
		y = scattered_positions(num_scattered,2)
		z = scattered_positions(num_scattered,3)

		! Compute the magnitude of the projection vector on the xz-plane.
		s = dsqrt(x**2 + z**2)

		! Calculate scattering angles in degrees
		alpha = datan2(y, s)*180/PI
		beta = -datan2(x, z)*180/PI

	end subroutine compute_scattering_angles

end module m3_trajectory_computation