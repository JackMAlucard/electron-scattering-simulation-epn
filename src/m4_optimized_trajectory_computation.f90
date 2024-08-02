module m4_optimized_trajectory_computation

	use m0_utilities, &
		only: dp, i8, INTERATOMIC_DIST_SIO2_INV, CELL_SCALE_FACTOR, &
		CELL_LENGTH_INV, MAX_SUPER_ELECTRON_CHARGE, MATERIAL_HEIGHT_SIO2, &
		CROSS_SECTION_SIO2, EFFECTIVE_DISTANCE, random_exponential

	use m3_trajectory_computation, &
		only: time_step, acceleration_due_to_electron, acceleration_due_to_atom

	implicit none
	contains
	!=============================================================================
	! Subroutine: setup_cells_and_super_electrons
	! Purpose   : Partition the material zone into cells by initializing the
	!             super electron number, charges, and positions arrays.
	!             - The partition cells are used to merge embedded electrons
	!               into super electrons in each cell, which will in turn
	!               allow to compute an approximation of the time step updates
	!               in the Velocity Verlet algorithm used to to simulate a
	!               projectile electron's trajectory.
	!             - Each partition cell can only store up to a maximum amount of
	!               super electrons determined by the MAX_SUPER_ELECTRON_CHARGE
	!               numerical constant.
	!             - The partition cells indices boundaries are defined by the
	!               material grid index boundaries and the CELL_SCALE_FACTOR
	!               numerical constant.
	!             - Partition cells are defined by the indices i, j, and k.
	!               The range of each index is determined by the values in the
	!               partition_boundaries = (pbi, pbj, pbk) array as follows:
	!               - i ranges from -pbi to pbi skipping zero (2*pbi values),
	!               - j ranges from -pbj to -1 (pbj values),
	!               - k ranges from -pbk to pbk skipping zero (2*pbk values).
	! Arguments :
	!   - integer(i8), intent(in) :: num_electrons
	!       Total number of electrons in the simulation.
	!   - integer(i8), intent(in) :: material_boundaries(3)
	!       Material grid indices boundaries.
	!   - integer(i8), intent(out) :: partition_boundaries(3)
	!       Partition cells indices boundaries.
	!   - integer(i8), allocatable, intent(out) :: num_super_electrons(:,:,:)
	!       Number of super electrons in each cell.
	!   - integer, allocatable, intent(out) :: super_electron_charges(:,:,:,:)
	!       Charges of the super electrons in each cell in atomic units (e).
	!   - real(dp), allocatable, intent(out) :: super_electron_positions &
	!       (:,:,:,:,:)
	!       Positions of the super electrons in each cell in atomic units (a0).
	!=============================================================================
	subroutine setup_cells_and_super_electrons &
		(num_electrons, material_boundaries, partition_boundaries, &
		num_super_electrons, super_electron_charges, super_electron_positions)
		implicit none

		! Input/Output variables
		integer(i8), intent(in) :: num_electrons, material_boundaries(3)
		integer(i8), intent(out) :: partition_boundaries(3)
		integer(i8), allocatable, intent(out) :: num_super_electrons(:,:,:)
		integer, allocatable, intent(out) :: super_electron_charges(:,:,:,:)
		real(dp), allocatable, intent(out) :: super_electron_positions(:,:,:,:,:)

		! Local variables
		integer(i8) :: pbi, pbj, pbk
		integer :: max_num_super_electrons

		! Compute partition cells indices boundaries
		pbi = material_boundaries(1)/CELL_SCALE_FACTOR + 1
		pbj = material_boundaries(2)/CELL_SCALE_FACTOR + 1
		pbk = material_boundaries(3)/CELL_SCALE_FACTOR + 1
		partition_boundaries = (/pbi, pbj, pbk/)

		! Set maximum possible number of super electrons in a cell, considering
		! the extreme case in which all electrons are embedded in a single cell
		max_num_super_electrons = (num_electrons - 1)/MAX_SUPER_ELECTRON_CHARGE + 1

		! Allocate array for the number of super electrons in each cell
		allocate (num_super_electrons(-pbi:pbi,-pbj:-1,-pbk:pbk))
		! Initialize as zero, since there are no super electrons at the start
		num_super_electrons = 0

		! Allocate array for the charges of super electrons in each cell
		allocate &
		(super_electron_charges(-pbi:pbi,-pbj:-1,-pbk:pbk, &
			max_num_super_electrons))
		! Initialize as zero, since there are no super electrons at the start
		super_electron_charges = 0

		! Allocating array that stores the super electrons positions on all cells
		allocate &
		(super_electron_positions(-pbi:pbi,-pbj:-1,-pbk:pbk, &
			max_num_super_electrons,3))
		! Initialize as zero, the positions will be given as the mean of the
		! positions of the electrons that form the super electron
		super_electron_positions = 0

	end subroutine setup_cells_and_super_electrons

	!=============================================================================
	! Subroutine: get_cell_indexes
	! Purpose   : Determine the cell indices within the partitioned material
	!             zone corresponding to a given position considering the
	!             partition cells' boundary limits.
	! Arguments :
	!   - real(dp), intent(in) :: r(3)
	!       Position vector in atomic units (a0)
	!   - integer(i8), intent(in) :: partition_boundaries(3)
	!       Partition cells indices boundaries.
	!   - integer(i8), intent(out) :: cell_indexes(3)
	!       Computed cell indices corresponding to the position vector.
	!=============================================================================
	subroutine get_cell_indexes(r, partition_boundaries, cell_indexes)
		implicit none

		! Input/Output variables
		real(dp), intent(in) :: r(3)
		integer(i8), intent(in) :: partition_boundaries(3)
		integer(i8), intent(out) :: cell_indexes(3)

		! Local variables
		integer(i8) :: pbi, pbj, pbk
		integer(i8) :: i, j, k

		! Extract partition cells index boundaries
		pbi = partition_boundaries(1)
		pbj = partition_boundaries(2)
		pbk = partition_boundaries(3)

		! Calculate the index in the x-direction and constrain within boundaries
		if (r(1) .ge. 0) then
			i = dint(r(1)*CELL_LENGTH_INV) + 1
			if (i .gt. pbi) i = pbi
		else
			i = dint(r(1)*CELL_LENGTH_INV) - 1
			if (i .lt. -pbi) i = -pbi
		end if

		! Calculate the index in the y-direction and constrain within boundaries
		if (r(2) .ge. 0) then
			j = -1
		else
			j = dint(r(2)*CELL_LENGTH_INV) - 1
			if (j .lt. -pbj) j = -pbj
		end if

		! Calculate the index in the z-direction and constrain within boundaries
		if (r(3) .ge. 0) then
			k = dint(r(3)*CELL_LENGTH_INV) + 1
			if (k .gt. pbk) k = pbk
		else
			k = dint(r(3)*CELL_LENGTH_INV) - 1
			if (k .lt. -pbk) k = -pbk
		end if

		! Store the partition cell indices
		cell_indexes = (/i, j, k/)

	end subroutine get_cell_indexes

	!=============================================================================
	! Subroutine: update_super_electron_in_cell
	! Purpose   : Update the super electron arrays in a specific cell based on
	!             the position of an embedded electron.
	! Arguments :
	!   - real(dp), intent(in) :: r(3)
	!       Position vector of the embedded electron in atomic units (a0).
	!   - integer(i8), intent(in) :: partition_boundaries(3)
	!       Partition cells indices boundaries.
	!   - integer(i8), intent(inout) :: num_super_electrons &
	!       (-partition_boundaries(1):, -partition_boundaries(2):, &
	!       -partition_boundaries(3):)
	!       Number of super electrons in each cell.
	!   - integer, intent(inout) :: super_electron_charges &
	!       (-partition_boundaries(1):, -partition_boundaries(2):, &
	!       -partition_boundaries(3):, :)
	!       Charges of the super electrons in each cell in atomic units (e).
	!   - real(dp), intent(inout) :: super_electron_positions &
	!       (-partition_boundaries(1):, -partition_boundaries(2):, &
	!       -partition_boundaries(3):, :, :)
	!       Positions of the super electrons in each cell in atomic units (a0).
	!=============================================================================
	subroutine update_super_electron_in_cell &
		(r, partition_boundaries, num_super_electrons, &
		super_electron_charges, super_electron_positions)
		implicit none

		! Input/Output variables
		real(dp), intent(in) :: r(3)
		integer(i8), intent(in) :: partition_boundaries(3)
		integer(i8), intent(inout) :: num_super_electrons &
			(-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):)
		integer, intent(inout) :: super_electron_charges &
			(-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :)
		real(dp), intent(inout) :: super_electron_positions &
			(-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :, :)

		! Local variables
		integer(i8) :: cell_indexes(3), i, j, k, n
		integer :: super_electron_charge
		real(dp) :: super_electron_r(3)

		! Get embedded electron's cell indices based on its position r
		call get_cell_indexes(r, partition_boundaries, cell_indexes)
		i = cell_indexes(1)
		j = cell_indexes(2)
		k = cell_indexes(3)

		! Get number of super electrons in the cell
		n = num_super_electrons(i,j,k)

		! Cell status: There are NO super electrons in the cell yet
		if (n .eq. 0) then
			! Initialize the first super electron in the cell and update arrays
			n = 1
			num_super_electrons(i,j,k) = n
			super_electron_charges(i,j,k,n) = 1
			super_electron_positions(i,j,k,n,:) = r

		! Cell status: There is at least one super electron in the cell, and
		! the current super electron has NOT reached its maximum charge value
		else if (super_electron_charge .lt. MAX_SUPER_ELECTRON_CHARGE) then
			! Read current super electron charge and position
			super_electron_charge = super_electron_charges(i,j,k,n)
			super_electron_r = super_electron_positions(i,j,k,n,:)

			! Add extra charge and compute new mean position
			super_electron_charge = super_electron_charge + 1
			super_electron_r = &
			super_electron_r + (r - super_electron_r)/super_electron_charge

			! Update current super electron charge and position arrays
			super_electron_charges(i,j,k,n) = super_electron_charge
			super_electron_positions(i,j,k,n,:) = super_electron_r

		! Cell status: There is at least one super electron in the cell, and
		! the current super electron has reached its maximum charge value
		else if (super_electron_charge .eq. MAX_SUPER_ELECTRON_CHARGE) then
			! Update current super electron number
			n = n + 1

			! Initialize the next super electron in the cell and update arrays
			num_super_electrons(i,j,k) = n
			super_electron_charges(i,j,k,n) = 1
			super_electron_positions(i,j,k,n,:) = r

		end if

		! Print updated super electron info to console
		print*, "Super electron arrays updated!"
		print*, " Cell indexes:", i, j, k
		print*, " New super electron number:", n
		print*, " New super electron charge:", super_electron_charge
		print*, " New super electron position:", super_electron_r

	end subroutine update_super_electron_in_cell

	!=============================================================================
	! Subroutine: acceleration_due_to_super_electron
	! Purpose   : Calculate the acceleration vector experienced by a projectile
	!             electron due to the electrostatic interaction with a stationary
	!             target super electron, based on their positions.
	! Arguments :
	!   - real(dp), intent(in) :: rp(3)
	!       Position vector of the projectile electron in atomic units (a0).
	!   - real(dp), intent(in) :: rt(3)
	!       Position vector of the target super electron in atomic units (a0/aut).
	!   - integer, intent(in) :: Q
	!       Charge of the super electron in atomic units (e).
	!   - real(dp), intent(out) :: a(3)
	!       Acceleration vector experienced by the projectile electron due to
	!       the stationary target super electron in atomic units (a0/aut^2).
	!=============================================================================
	subroutine acceleration_due_to_super_electron(rp, rt, Q, a)
		implicit none

		! Input/Output variables
		real(dp), intent(in) :: rp(3)	! Position of the projectile electron (a0)
		real(dp), intent(in) :: rt(3)	! Position of the target super electron (a0)
		integer, intent(in) :: Q			! Charge of the target super electron (e)
		real(dp), intent(out) :: a(3)	! Resulting acceleration vector (a0/aut^2)

		! Local variables
		real(dp) :: rs(3)	! Separation vector between the projectile electron and
		                  ! target super electron (a0)
		real(dp) :: r 		! Magnitude of the separation vector (a0)

		! Compute the separation vector and its magnitude
		rs = rp - rt
		r = norm2(rs)

		! Determine the acceleration using Coulomb's law in atomic units
		a = (Q*rs)/(r**3)

	end subroutine acceleration_due_to_super_electron

	!=============================================================================
	! Subroutine: time_step_approximate
	! Purpose   : Perform an approximate time step for a projectile electron,
	!             including position, velocity, and acceleration approximate
	!             updates, using the Velocity Verlet algorithm to simulate
	!             its trajectory.
	!             - The algorithm only considers interactions with all embedded
	!               electrons through the corresponding super electrons.
	! Arguments :
	!   - integer(i8), intent(in) :: num_embedded
	!       Number of embedded electrons.
	!   - integer(i8), intent(in) :: partition_boundaries(3)
	!       Partition cells indices boundaries.
	!   - integer(i8), intent(inout) :: num_super_electrons &
	!       (-partition_boundaries(1):, -partition_boundaries(2):, &
	!       -partition_boundaries(3):)
	!       Number of super electrons in each cell.
	!   - integer, intent(inout) :: super_electron_charges &
	!       (-partition_boundaries(1):, -partition_boundaries(2):, &
	!       -partition_boundaries(3):, :)
	!       Charges of the super electrons in each cell in atomic units (e).
	!   - real(dp), intent(inout) :: super_electron_positions &
	!       (-partition_boundaries(1):, -partition_boundaries(2):, &
	!       -partition_boundaries(3):, :, :)
	!       Positions of the super electrons in each cell in atomic units (a0).
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
	subroutine time_step_approximate &
		(num_embedded, partition_boundaries, num_super_electrons, &
		super_electron_charges, super_electron_positions, dt, r, v, a)
		implicit none

		! Input/Output variables
		integer(i8), intent(in) :: num_embedded, partition_boundaries(3)
		integer(i8), intent(in) :: num_super_electrons(-partition_boundaries(1):, &
			-partition_boundaries(2):, -partition_boundaries(3):)
		integer, intent(in) :: super_electron_charges(-partition_boundaries(1):, &
			-partition_boundaries(2):, -partition_boundaries(3):, :)
		real(dp), intent(in) :: super_electron_positions &
		(-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :, :)
		real(dp), intent(in) :: dt
		real(dp), intent(inout) :: r(3), v(3), a(3)

		! Local variables
		real(dp) :: rt(3), ak(3)
		integer :: Q
		integer(i8) :: pbi, pbj, pbk, super_electron_num
		integer(i8) :: i, j, k, n

		! Half-step velocity update
		v = v + 0.5*a*dt

		! Position update
		r = r + v*dt

		! Full-step acceleration update
		! Reset acceleration for the current time step
		a = 0

		! Extract partition cells index boundaries
		pbi = partition_boundaries(1)
		pbj = partition_boundaries(2)
		pbk = partition_boundaries(3)

		! Acceleration due to super electrons;
		! only computed if there is at least one embedded electron
		if (num_embedded .gt. 0) then
			do i = -pbi, pbi
				! Omit array positions with i index 0 to reduce overhead;
				! since it is not a valid cell index it does not contain super electrons
				if (i .ne. 0) then

					do j = -1, -pbj, -1
						do k = -pbk, pbk
							! Get number of super electrons in the cell
							super_electron_num = num_super_electrons(i,j,k)

							! Only compute acceleration in cells with at least one super
							! electron; this condition also omits invalid cells with j index
							! 0 since they remain initialized with super_electron_num 0
							if (super_electron_num .gt. 0) then
								do n = 1, super_electron_num
									Q = super_electron_charges(i,j,k,n)
									rt = super_electron_positions(i,j,k,n,:)
									call acceleration_due_to_super_electron(r, rt, Q, ak)
									a = a + ak
								end do
							end if

						end do
					end do

				end if
			end do
		end if

		! Second half-step velocity update
		v = v + 0.5*a*dt

	end subroutine time_step_approximate

	!=============================================================================
	! Subroutine: compute_trajectory_optimized
	! Purpose   : Simulate the trajectory of a projectile electron aimed at the
	!             dielectric material considering interactions either with both
	!             the embedded electrons and material atoms, or only with the
	!             corresponding super electrons depending on the electron's
	!             distance to the target material in order to optimize total
	!             simulation time. The trajectory is computed using the Velocity
	!             Verlet algorithm. When the trajectory simulation is completed,
	!             the projectile electron can either get embedded in the material
	!             or be scattered.
	!             - Each trajectory is simulated using a combination of regular
	!               and approximated time step updates as implemented in the
	!               time_step and time_step_approximate subroutines, respectively.
	!             - Each time step subroutine considers the following interacions:
	!               - Embedded electrons and material atoms (time_step),
	!               - Super electrons corresponding to embedded electrons
	!                 (time_step_approximate).
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
	!   - integer(i8), intent(in) :: partition_boundaries(3)
	!       Partition cell indices boundaries.
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
	!   - integer(i8), intent(inout) :: num_super_electrons &
	!       (-partition_boundaries(1):, -partition_boundaries(2):, &
	!       -partition_boundaries(3):)
	!       Number of super electrons in each cell.
	!   - integer, intent(inout) :: super_electron_charges &
	!       (-partition_boundaries(1):, -partition_boundaries(2):, &
	!       -partition_boundaries(3):, :)
	!       Charges of the super electrons in each cell in atomic units (e).
	!   - real(dp), intent(inout) :: super_electron_positions &
	!       (-partition_boundaries(1):, -partition_boundaries(2):, &
	!       -partition_boundaries(3):, :, :)
	!       Positions of the super electrons in each cell in atomic units (a0).
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
	subroutine compute_trajectory_optimized &
		(electron_trajectories_saving_enabled, num_plot_ploints, max_iterations, &
		output_unit, material_boundaries, atom_positions, atomic_numbers, &
		atomic_numbers_cbrt, partition_boundaries, num_super_electrons, &
		super_electron_positions, super_electron_charges, dt, r, v, a, &
		num_embedded, num_scattered, embedded_positions, scattered_positions)
		implicit none

		! Input/Output variables
		logical, intent(in) :: electron_trajectories_saving_enabled
		integer(i8), intent(in) :: num_plot_ploints, max_iterations
		integer, intent(in) :: output_unit
		integer(i8), intent(in) :: material_boundaries(3), partition_boundaries(3)
		real(dp), intent(in) :: atom_positions(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:, :)
		integer, intent(in) :: atomic_numbers(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:)
		real(dp), intent(in) :: atomic_numbers_cbrt(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:)
		real(dp), intent(in) :: dt
		integer(i8), intent(inout) :: num_super_electrons &
			(-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):)
		integer, intent(inout) :: super_electron_charges &
			(-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :)
		real(dp), intent(inout) :: super_electron_positions &
			(-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :, :)
		real(dp), intent(inout) :: r(3), v(3), a(3)
		integer(i8), intent(inout) :: num_embedded, num_scattered
		real(dp), intent(inout) :: embedded_positions(:,:), scattered_positions(:,:)

		! Local variables
		real(dp) :: initial_distance_to_target
		logical :: is_embedded, is_scattered, is_max_iteration
		logical :: in_material
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

			! Compute next time step using Velocity Verlet algorithm;
			! Update electron position r, velocity v, and acceleration a;
			! Time step update: Approximated, projectile electron's distance to target
			! material is greater than EFFECTIVE_DISTANCE and not in material zone
			if (distance_to_target .gt. EFFECTIVE_DISTANCE .and. &
			.not. in_material) then
				call time_step_approximate &
				(num_embedded, partition_boundaries, num_super_electrons, &
				super_electron_charges, super_electron_positions, dt, r, v, a)

			! Time step update: Regular, projectile electron's distance to target
			! material is less than EFFECTIVE_DISTANCE or in material zone
			else
				call time_step &
				(num_embedded, material_boundaries, embedded_positions, &
				atom_positions, atomic_numbers, atomic_numbers_cbrt, dt, r, v, a)

			end if

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
					is_embedded = .true.
					num_embedded = num_embedded + 1
					embedded_positions(num_embedded,:) = r

					! Update embedded condition flag, embedded number and positions
					call update_super_electron_in_cell &
					(r, partition_boundaries, num_super_electrons, &
					super_electron_charges, super_electron_positions)

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

	end subroutine compute_trajectory_optimized

	!=============================================================================
	! Subroutine: write_final_super_electron_distribution
	! Purpose   : Write the final distribution of super electrons, including
	!             their charges and positions, to a file.
	! Arguments :
	!   - integer(i8), intent(in) :: partition_boundaries(3)
	!       Partition cells indices boundaries.
	!   - integer(i8), intent(inout) :: num_super_electrons &
	!       (-partition_boundaries(1):, -partition_boundaries(2):, &
	!       -partition_boundaries(3):)
	!       Number of super electrons in each cell.
	!   - integer, intent(inout) :: super_electron_charges &
	!       (-partition_boundaries(1):, -partition_boundaries(2):, &
	!       -partition_boundaries(3):, :)
	!       Charges of the super electrons in each cell in atomic units (e).
	!   - real(dp), intent(inout) :: super_electron_positions &
	!       (-partition_boundaries(1):, -partition_boundaries(2):, &
	!       -partition_boundaries(3):, :, :)
	!       Positions of the super electrons in each cell in atomic units (a0).
	!=============================================================================
	subroutine write_final_super_electron_distribution &
		(partition_boundaries, num_super_electrons, super_electron_charges, &
		super_electron_positions)
		implicit none

		! Input variables
		integer(i8), intent(in) :: partition_boundaries(3)
		integer(i8), intent(in) :: num_super_electrons &
			(-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):)
		integer, intent(in) :: super_electron_charges &
			(-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :)
		real(dp), intent(in) :: super_electron_positions &
			(-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :, :)

		! Local variables
		integer(i8) :: pbi, pbj, pbk
		integer(i8) :: i, j, k, n

		! Extract partition cells index boundaries
		pbi = partition_boundaries(1)
		pbj = partition_boundaries(2)
		pbk = partition_boundaries(3)

		! Open file to save the final super electron distribution
		open(unit=42, file='final_super_electron_distribution.dat', &
			status='replace', action='write')

			! Loop over all partition cells and write super electron data to the file
			do i = -pbi, pbi
				do j = -1, -pbj, -1
					do k = -pbk, pbk

						! Get super electron(s) charges and positions if there is at least
						! one super electron in the cell
						if (num_super_electrons(i,j,k) .gt. 0) then
							do n = 1, num_super_electrons(i,j,k)
								write(42, *) i, j, k, n, super_electron_charges(i,j,k,n), &
								super_electron_positions(i,j,k,n,:)
							end do
						end if

					end do
				end do
			end do

		! Close output file
		close(42)

	end subroutine write_final_super_electron_distribution

end module m4_optimized_trajectory_computation