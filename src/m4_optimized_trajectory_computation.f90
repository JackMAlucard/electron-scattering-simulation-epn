module m4_optimized_trajectory_computation
	use m0_utilities, &
		only: dp, i8, INTERATOMIC_DIST_SIO2_INV, CELL_SCALE_FACTOR, &
		CELL_LENGTH_INV, MAX_EQUIVALENT_CHARGE, MATERIAL_HEIGHT_SIO2, &
		CROSS_SECTION_SIO2, EFFECTIVE_DISTANCE, random_exponential
	use m3_trajectory_computation, &
		only: time_step, acceleration_due_to_electron, acceleration_due_to_atom
	implicit none
	contains
	!=======================================================================
	! Subroutine: setup_cells_and_super_electrons
	! Purpose   : Initialize the cell partitioning of the simulation space 
	!             and prepare storage for super electrons within each cell.
	! Arguments :
	!   - integer(i8), intent(in) :: num_electrons
	!       Total number of electrons in the simulation.
	!   - integer(i8), intent(in) :: material_boundaries(3)
	!       Boundaries of the material grid.
	!   - integer(i8), intent(out) :: partition_boundaries(3)
	!       Boundaries of the partitioned cells.
	!   - integer(i8), allocatable, intent(out) :: num_super_electrons(:,:,:)
	!       Number of super electrons in each cell.
	!   - integer, allocatable, intent(out) :: super_electron_charges(:,:,:,:)
	!       Charges of the super electrons in each cell.
	!   - real(dp), allocatable, intent(out) :: super_electron_positions(:,:,:,:,:)
	!       Positions of the super electrons in each cell.
	! Variables :
	!   - integer(i8) :: pbi, pbj, pbk
	!       Partition boundaries in each dimension.
	!   - integer :: max_super_electrons
	!       Maximum possible number of super electrons in a cell.
	!=======================================================================
	! Subroutine that defines the partition boundaries. It also allocates and 
	! initializes the arrays that store the super electrons information.
	subroutine setup_cells_and_super_electrons &
		(num_electrons, material_boundaries, partition_boundaries, &
		num_super_electrons, super_electron_charges, super_electron_positions)
		implicit none
		integer(i8), intent(in) :: num_electrons, material_boundaries(3)
		integer(i8), intent(out) :: partition_boundaries(3)
		integer(i8), allocatable, intent(out) :: num_super_electrons(:,:,:)
		integer, allocatable, intent(out) :: super_electron_charges(:,:,:,:)
		real(dp), allocatable, intent(out) :: super_electron_positions(:,:,:,:,:)
		integer(i8) :: pbi, pbj, pbk
		integer :: max_super_electrons
		! Computing partition boundaries
		pbi = material_boundaries(1)/CELL_SCALE_FACTOR + 1
		pbj = material_boundaries(2)/CELL_SCALE_FACTOR + 1
		pbk = material_boundaries(3)/CELL_SCALE_FACTOR + 1
		partition_boundaries = (/pbi, pbj, pbk/)
		! num_cells = (2*pbi)*(pbj)*(2*pbk)
		! Maximum possible number of super electrons on a cell, considering the
		! extreme case in which all of the electrons are embedded into a single cell
		max_super_electrons = (num_electrons - 1)/MAX_EQUIVALENT_CHARGE + 1
		! Allocating array that stores the number of super electrons on all cells
		allocate (num_super_electrons(-pbi:pbi,-pbj:-1,-pbk:pbk))
		! Initializing as zero, since there are no super electrons at the start
		num_super_electrons = 0
		! Allocating array that stores the super electrons charges on all cells
		allocate &
		(super_electron_charges(-pbi:pbi,-pbj:-1,-pbk:pbk,max_super_electrons))
		! Initializing as zero, since there are no super electrons at the start
		super_electron_charges = 0
		! Allocating array that stores the super electrons positions on all cells
		allocate &
		(super_electron_positions(-pbi:pbi,-pbj:-1,-pbk:pbk,max_super_electrons,3))
		! Initializing as zero, the positions will be given as the mean of the 
		! positions of the electrons that form the super electron
		super_electron_positions = 0
	end subroutine setup_cells_and_super_electrons
	!=======================================================================
	! Subroutine: get_cell_indexes
	! Purpose   : Determine the cell indexes corresponding to a given position 
	!             within the partitioned simulation space.
	! Arguments :
	!   - real(dp), intent(in) :: r(3)
	!       Position vector of the particle.
	!   - integer(i8), intent(in) :: partition_boundaries(3)
	!       Boundaries of the partitioned cells.
	!   - integer(i8), intent(out) :: cell_indexes(3)
	!       Computed cell indexes for the position vector.
	! Variables :
	!   - integer(i8) :: pbi, pbj, pbk
	!       Partition boundaries in each dimension.
	!   - integer(i8) :: i, j, k
	!       Cell indexes in each dimension.
	!=======================================================================
	! Subroutine that takes the position of an electron r(3), the inverse cell 
	! length (in a0) and determines the cell_indexes of the cubic cell to which 
	! the embedded electron belongs to in order to update the super electron in 
	! that cell. 
	! Electrons embedded after crossing the material boundaries will be included
	! in the cell at the boundary of the partitions.
	! The lowest possible index in the y-direction is -1 since there wouldn't be
	! any cells with index 0 in any direction. Since the arrays have those slots,
	! those MUST BE ignored, and as such were initialized and will remain zero.
	subroutine get_cell_indexes(r, partition_boundaries, cell_indexes)
		implicit none
		real(dp), intent(in) :: r(3)
		integer(i8), intent(in) :: partition_boundaries(3)
		integer(i8), intent(out) :: cell_indexes(3)
		integer(i8) :: pbi, pbj, pbk
		integer(i8) :: i, j, k
		! Getting partition cells boundaries
		pbi = partition_boundaries(1)
		pbj = partition_boundaries(2)
		pbk = partition_boundaries(3)
		if (r(1) .ge. 0) then
			i = dint(r(1)*CELL_LENGTH_INV) + 1
			if (i .gt. pbi) i = pbi
		else
			i = dint(r(1)*CELL_LENGTH_INV) - 1
			if (i .lt. -pbi) i = -pbi
		end if
		if (r(2) .ge. 0) then
			j = -1 ! The top cell includes a bit more space
		else
			j = dint(r(2)*CELL_LENGTH_INV) - 1
			if (j .lt. -pbj) j = -pbj
		end if
		if (r(3) .ge. 0) then
			k = dint(r(3)*CELL_LENGTH_INV) + 1
			if (k .gt. pbk) k = pbk
		else
			k = dint(r(3)*CELL_LENGTH_INV) - 1
			if (k .lt. -pbk) k = -pbk
		end if
		cell_indexes = (/i, j, k/)
	end subroutine get_cell_indexes
	!=======================================================================
	! Subroutine: update_super_electron_in_cell
	! Purpose   : Update the super electron properties in a specific cell based on 
	!             the position of an embedded electron.
	! Arguments :
	!   - real(dp), intent(in) :: r(3)
	!       Position vector of the embedded electron.
	!   - integer(i8), intent(in) :: partition_boundaries(3)
	!       Boundaries of the partitioned cells.
	!   - integer(i8), intent(inout) :: num_super_electrons(-partition_boundaries(1):, 
	!                     -partition_boundaries(2):, -partition_boundaries(3):)
	!       Array storing the number of super electrons in each cell.
	!   - integer, intent(inout) :: super_electron_charges(-partition_boundaries(1):, 
	!                     -partition_boundaries(2):, -partition_boundaries(3):, :)
	!       Array storing the charge of super electrons in each cell.
	!   - real(dp), intent(inout) :: super_electron_positions(-partition_boundaries(1):, 
	!                     -partition_boundaries(2):, -partition_boundaries(3):, :, :)
	!       Array storing the positions of super electrons in each cell.
	! Variables :
	!   - integer(i8) :: cell_indexes(3), i, j, k, n
	!       Indexes of the cell corresponding to the position vector and loop variables.
	!   - integer :: super_electron_charge
	!       Charge of the current super electron.
	!   - real(dp) :: super_electron_r(3)
	!       Position of the current super electron.
	!=======================================================================
	! Subroutine that updated the values of the super electron arrays in the
	! cell in which the projectile electron with position r gets embedded.
	subroutine update_super_electron_in_cell &
		(r, partition_boundaries, num_super_electrons, & 
		super_electron_charges, super_electron_positions)
		implicit none
		real(dp), intent(in) :: r(3)
		integer(i8), intent(in) :: partition_boundaries(3)
		integer(i8), intent(inout) :: num_super_electrons( &
			-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):)
		integer, intent(inout) :: super_electron_charges( &
			-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :)
		real(dp), intent(inout) :: super_electron_positions( &
			-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :, :)
		integer(i8) :: cell_indexes(3), i, j, k, n
		integer :: super_electron_charge
		real(dp) :: super_electron_r(3)
		! Getting embedded electron's cell based on its position r
		call get_cell_indexes(r, partition_boundaries, cell_indexes)
		i = cell_indexes(1)
		j = cell_indexes(2)
		k = cell_indexes(3)
		! Get number of super electrons in the cell
		n = num_super_electrons(i,j,k)
		! If there are no embedded electrons in the cell, yet,
		! i.e. no super electrons yet in the cell
		if (n .eq. 0) then
			! Updating super electrons arrays initializing the first one 
			n = 1
			num_super_electrons(i,j,k) = n
			super_electron_charges(i,j,k,n) = 1
			super_electron_positions(i,j,k,n,:) = r
		! If current super electron is not full of charge, add extra charge
		! and compute new super electron position
		else if (super_electron_charge .lt. MAX_EQUIVALENT_CHARGE) then
			! Read current super electron charge and position
			super_electron_charge = super_electron_charges(i,j,k,n)
			super_electron_r = super_electron_positions(i,j,k,n,:)
			! Compute new super electron charge and position
			super_electron_charge = super_electron_charge + 1
			super_electron_r = &
			super_electron_r + (r - super_electron_r)/super_electron_charge
			! Update super electron charge and position
			super_electron_charges(i,j,k,n) = super_electron_charge
			super_electron_positions(i,j,k,n,:) = super_electron_r
		! If current super electron is full of charge, initialize next one
		else if (super_electron_charge .eq. MAX_EQUIVALENT_CHARGE) then
			n = n + 1
			! Updating super electrons arrays initializing the next one
			num_super_electrons(i,j,k) = n
			super_electron_charges(i,j,k,n) = 1
			super_electron_positions(i,j,k,n,:) = r
		else
			!Everything else is an error
			print*, "ERROR UPDATING SUPER ELECTRON IN CELL"
		end if
		! Printing updated super electron info
		print*, "Super electron arrays updated!"
		print*, " Cell indexes:", i, j, k
		print*, " New super electron number:", n
		print*, " New super electron charge:", super_electron_charge
		print*, " New super electron position:", super_electron_r
	end subroutine update_super_electron_in_cell
	!=======================================================================
	! Subroutine : acceleration_due_to_super_electron
	! Purpose    : Calculate the acceleration on a particle due to a super electron.
	! Arguments  :
	!   - real(dp), intent(in) :: rp(3)
	!       Position vector of the particle experiencing the acceleration.
	!   - real(dp), intent(in) :: rt(3)
	!       Position vector of the super electron.
	!   - integer, intent(in) :: Q
	!       Charge of the super electron.
	!   - real(dp), intent(out) :: a(3)
	!       Resultant acceleration vector due to the super electron.
	! Variables  :
	!   - real(dp) :: rs(3)
	!       Separation vector between the particle and the super electron.
	!   - real(dp) :: r
	!       Magnitude of the separation vector.
	!=======================================================================
	! Acceleration of a projectile electron interacting with an embedded 
	! super electron of charge q via the electrostatic Coulomb potential. 
	! The input and output variables are in atomic units (au).
	! rp: r projectile, position of the incident electron (a0)
	! rt: r target, equivalent position of the stationary electron (a0)
	! a: Acceleration of the incident electron (aua)
	subroutine acceleration_due_to_super_electron(rp, rt, Q, a)
		implicit none
		real(dp), intent(in) :: rp(3), rt(3)
		integer, intent(in) :: Q
		real(dp), intent(out) :: a(3)
		real(dp) :: rs(3), r !rs, r: separation vector and its magnitude
		rs = rp - rt
		r = norm2(rs)
		a = (Q*rs)/(r**3)
	end subroutine acceleration_due_to_super_electron
	!=======================================================================
	! Subroutine : time_step_approximate
	! ORIGINAL:
	! Purpose    : Perform a time step update for a projectile electron in the far zone,
	!              including position, velocity, and acceleration updates.
	! UPDATED AS APPROXIMATE:
	! Purpose    : Perform a single time step update for the position, velocity,
	!              and acceleration of a projectile electron using an approximate
	!              method. This subroutine computes the influence of embedded
	!              super electrons on the projectile electron's acceleration.
	! Arguments  :
	!   - integer(i8), intent(in) :: num_embedded
	!       Number of embedded electrons.
	!   - integer(i8), intent(in) :: partition_boundaries(3)
	!       Boundaries of the partition cells.
	!   - integer(i8), intent(in) :: num_super_electrons(-partition_boundaries(1):, &
	!         -partition_boundaries(2):, -partition_boundaries(3):)
	!       Number of super electrons in each cell.
	!   - integer, intent(in) :: super_electron_charges(-partition_boundaries(1):, &
	!         -partition_boundaries(2):, -partition_boundaries(3):, :)
	!       Charges of the super electrons in each cell.
	!   - real(dp), intent(in) :: super_electron_positions( &
	!         -partition_boundaries(1):, -partition_boundaries(2):, &
	!         -partition_boundaries(3):, :, :)
	!       Positions of the super electrons in each cell.
	!   - real(dp), intent(in) :: dt
	!       Time step for the simulation.
	!   - real(dp), intent(inout) :: r(3)
	!       Position vector of the projectile electron.
	!   - real(dp), intent(inout) :: v(3)
	!       Velocity vector of the projectile electron.
	!   - real(dp), intent(inout) :: a(3)
	!       Acceleration vector of the projectile electron.
	! Variables  :
	!   - real(dp) :: rt(3)
	!       Position vector of the target super electron.
	!   - real(dp) :: ak(3)
	!       Acceleration vector due to from a single super electron.
	!   - integer :: Q
	!       Charge of the super electron.
	!   - integer(i8) :: pbi, pbj, pbk
	!       Boundaries of the partition cells.
	!   - integer(i8) :: super_electron_num
	!       Number of super electrons in the current/in a specific cell.
	!   - integer(i8) :: i, j, k, n
	!       Loop counters for traversing through the partition cells and super electrons.
	!=======================================================================
	! TO BE EDITED
	! Velocity Verlet Step for Far Field (bulk electrons)
	! Subroutine which computes a time step of the Velocity Verlet algorithm used to
	! compute the trajectory of a projectile electron interacting with Nb equivalent
	! embedded electrons characterized by the b_neq, b_ec, and b_er arrays. 
	subroutine time_step_approximate &
		(num_embedded, partition_boundaries, num_super_electrons, &
		super_electron_charges, super_electron_positions, dt, r, v, a)
		implicit none
		integer(i8), intent(in) :: num_embedded, partition_boundaries(3)
		integer(i8), intent(in) :: num_super_electrons(-partition_boundaries(1):, &
			-partition_boundaries(2):, -partition_boundaries(3):)
		integer, intent(in) :: super_electron_charges(-partition_boundaries(1):, &
			-partition_boundaries(2):, -partition_boundaries(3):, :)
		real(dp), intent(in) :: super_electron_positions( &
			-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :, :)
		real(dp), intent(in) :: dt
		real(dp), intent(inout) :: r(3), v(3), a(3)
		real(dp) :: rt(3), ak(3)
		integer :: Q
		integer(i8) :: pbi, pbj, pbk, super_electron_num
		integer(i8) :: i, j, k, n
		! Half-step velocity update
		v = v + 0.5*a*dt
		! Position update
		r = r + v*dt
		! Full-step acceleration update
		a = 0
		! Getting partition cells boundaries
		pbi = partition_boundaries(1)
		pbj = partition_boundaries(2)
		pbk = partition_boundaries(3)
		! Acceleration due to embedded super electrons
		! ONLY computed if there is AT LEAST one embedded electron
		if (num_embedded .gt. 0) then
			do i = -pbi, pbi
				! To avoid checking for unused cells with index 0
				if (i .ne. 0) then
					do j = -1, -pbj, -1
						do k = -pbk, pbk
							! Getting number of super electrons in the cell
							super_electron_num = num_super_electrons(i,j,k)
							! This also avoids unused cells, their super_electron_num is 0
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
	!=======================================================================
	! Subroutine : compute_trajectory_optimized
	! Purpose    : Compute the trajectory of a projectile electron using an optimized
	!              algorithm, including updates for position, velocity, and acceleration.
	!              The subroutine tracks the electron's path and determines if it gets
	!              embedded in the material, scattered, or reaches the maximum number of iterations.
	! Arguments  :
	!   - integer(i8), intent(in) :: num_plot_ploints
	!       Number of points to plot.
	!   - integer(i8), intent(in) :: max_iterations
	!       Maximum number of iterations for the trajectory computation.
	!   - integer, intent(in) :: output_unit
	!       Output unit for writing the trajectory data.
	!   - integer(i8), intent(in) :: material_boundaries(3)
	!       Boundaries of the material cells.
	!   - integer(i8), intent(in) :: partition_boundaries(3)
	!       Boundaries of the partition cells.
	!   - real(dp), intent(in) :: atom_positions(-material_boundaries(1)-1:, &
	!         -material_boundaries(2)-1:, -material_boundaries(3)-1:, :)
	!       Positions of the atoms within the material cells.
	!   - integer, intent(in) :: atom_charges(-material_boundaries(1)-1:, &
	!         -material_boundaries(2)-1:, -material_boundaries(3)-1:)
	!       Charges of the atoms within the material cells.
	!   - real(dp), intent(in) :: atom_charges_cbrt(-material_boundaries(1)-1:, &
	!         -material_boundaries(2)-1:, -material_boundaries(3)-1:)
	!       Cube root of the charges of the atoms within the material cells.
	!   - real(dp), intent(in) :: dt
	!       Time step for the simulation.
	!   - integer(i8), intent(inout) :: num_super_electrons( &
	!         -partition_boundaries(1):, -partition_boundaries(2):, &
	!         -partition_boundaries(3):)
	!       Number of super electrons in each cell.
	!   - integer, intent(inout) :: super_electron_charges( &
	!         -partition_boundaries(1):, -partition_boundaries(2):, &
	!         -partition_boundaries(3):, :)
	!       Charges of the super electrons in each cell.
	!   - real(dp), intent(inout) :: super_electron_positions( &
	!         -partition_boundaries(1):, -partition_boundaries(2):, &
	!         -partition_boundaries(3):, :, :)
	!       Positions of the super electrons in each cell.
	!   - real(dp), intent(inout) :: r(3)
	!       Position vector of the projectile electron.
	!   - real(dp), intent(inout) :: v(3)
	!       Velocity vector of the projectile electron.
	!   - real(dp), intent(inout) :: a(3)
	!       Acceleration vector of the projectile electron.
	!   - integer(i8), intent(inout) :: num_embedded
	!       Number of embedded electrons.
	!   - integer(i8), intent(inout) :: num_scattered
	!       Number of scattered electrons.
	!   - real(dp), intent(inout) :: embedded_positions(:,:)
	!       Positions of the embedded electrons.
	!   - real(dp), intent(inout) :: scattered_positions(:,:)
	!       Positions of the scattered electrons.
	! Variables  :
	!   - real(dp) :: initial_distance_to_target
	!       Initial distance from the projectile electron to the target.
	!   - logical :: is_embedded, is_scattered, is_max_iteration
	!       Flags indicating the end conditions of the trajectory computation.
	!   - logical :: in_material
	!       Flag indicating if the electron is inside the material.
	!   - real(dp) :: distance_before_collision, distance_in_material
	!       Distances related to the embedding process.
	!   - real(dp) :: previous_position(3), actual_position(3), step_length
	!       Variables for tracking the electron's position and step length.
	!   - real(dp) :: distance_to_target, t
	!       Variables for tracking the distance to the target and the current time step.
	!   - integer(i8) :: mbi, mbj, mbk
	!       Boundaries of the material cells.
	!   - integer(i8) :: i, j
	!       Loop counters.
	!=======================================================================
! TO BE EDITED, AND OMG, WHAT A CHORE IT LOOKS LIKE IT WILL BE
!Subroutine which computes the trajectory of the i-th projectile electron, out
!of N electrons in the electron beam, and plots P points in unit file number ou.
!This subroutine uses the Mixed Method, i.e. it computes the acceleration due to
!equivalent embedded electrons in the Far Field Regimen when the distance of the
!i-th projectile is greater than the Far Field effective distance, and the 
!acceleration due to each invidivual electron otherwise, as well as the
!acceleration due to the material atoms when it is very close to the material.
!The projectile interacts with a set of Ne embedded electrons out of a maximum
!of N total possible embedded electrons and a **SCC** bulk of 
!N_a=(2*Nx + 1)*(Ny + 1)*(2*Nz + 1) neutral atoms with atomic number Z.
!The beam initial positions, velocities, and accelerations are stored in the 
!(N,3) arrarys r, v and a respectively.
!The embedded and scattered electrons final positions are stored in the (N,3)
!arrays e_emb and e_sct respectively.
!The simulation is ended either when the electron distance dk to the origin of 
!coordinates (which approximates the center of the bulk) is greater than that 
!initial distance d0, or when the number of iterations k is equal or greater 
!than the estimated number of iterations T.
!The simulation is also ended if the electron gets embedded on the bulk, at
!which point the number of embedded electrons N_e is updated/increased.
!The information for the Mix Method simulation is stored in the additional 
!arrays that have been previously defined in the add_electron_to_bin
!subroutine.
	subroutine compute_trajectory_optimized &
		(electron_trajectories_saving_enabled, num_plot_ploints, max_iterations, &
		output_unit, material_boundaries, atom_positions, atom_charges, & 
		atom_charges_cbrt, partition_boundaries, num_super_electrons, &
		super_electron_positions, super_electron_charges, dt, r, v, a, &
		num_embedded, num_scattered, embedded_positions, scattered_positions)
		implicit none
		logical, intent(in) :: electron_trajectories_saving_enabled
		integer(i8), intent(in) :: num_plot_ploints, max_iterations
		integer, intent(in) :: output_unit
		integer(i8), intent(in) :: material_boundaries(3), partition_boundaries(3)
		real(dp), intent(in) :: atom_positions(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:, :)
		integer, intent(in) :: atom_charges(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:)
		real(dp), intent(in) :: atom_charges_cbrt(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:)
		real(dp), intent(in) :: dt
		integer(i8), intent(inout) :: num_super_electrons( &
			-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):)
		integer, intent(inout) :: super_electron_charges( &
			-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :)
		real(dp), intent(inout) :: super_electron_positions( &
			-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :, :)
		real(dp), intent(inout) :: r(3), v(3), a(3)
		integer(i8), intent(inout) :: num_embedded, num_scattered
		real(dp), intent(inout) :: embedded_positions(:,:), scattered_positions(:,:)
		real(dp) :: initial_distance_to_target
		logical :: is_embedded, is_scattered, is_max_iteration
		logical :: in_material
		real(dp) :: distance_before_collision, distance_in_material
		real(dp) :: previous_position(3), actual_position(3), step_length
		real(dp) :: distance_to_target, t
		integer(i8) :: mbi, mbj, mbk
		integer(i8) :: i, j
		! Initializing end conditions as false
		is_embedded = .false.
		is_scattered = .false.
		is_max_iteration = .false.
		! Initializing variables related to end conditions
		i = 0
		in_material = .false.
		initial_distance_to_target = norm2(r)
		distance_to_target = initial_distance_to_target
		! Initialize points to be plotted counter
		j = 0
		do while (.not.(is_embedded .or. is_scattered .or. is_max_iteration))
			t = i*dt ! t0 = 0, just to print to file, not needed for computations
			! Plotting only P points
			if (electron_trajectories_saving_enabled .and. &
				(mod(i,max_iterations/num_plot_ploints) .eq. 0) .and. &
				j .lt. num_plot_ploints) then
				! Write values to file
				write(output_unit,*) t, r
				j = j + 1
			end if
			! Computing next Velocity Verlet time step integration depending on the
			! projectile electron distance to target material center in relation to 
			! the EFFECTIVE_DISTANCE and whether or not it's in the material
			if (distance_to_target .gt. EFFECTIVE_DISTANCE .and. &
			.not. in_material) then
				! Computing next Velocity Verlet time step integration using far zone
				! approximation
				call time_step_approximate &
				(num_embedded, partition_boundaries, num_super_electrons, &
				super_electron_charges, super_electron_positions, dt, r, v, a)
			else
				! Computing next Velocity Verlet time step integration using near zone
				! approximation
				call time_step &
				(num_embedded, material_boundaries, embedded_positions, &
				atom_positions, atom_charges, atom_charges_cbrt, dt, r, v, a)
			end if
			! Updating variables related to end conditions for distance and iterations
			distance_to_target = norm2(r)
			i = i + 1
			! Check if any of the end conditions are met
			! End condition for embedded electrons
			! The electron must be below material level
			if (r(2) .lt. MATERIAL_HEIGHT_SIO2) then
				! ONLY the first time it enters the material
				! Generate random number following exponential distribution
				! Initialize actual electron position (inside material)
				! Initialize electron distance travelled inside material
				if (.not.(in_material)) then
					in_material = .true.
					call random_exponential(CROSS_SECTION_SIO2, distance_before_collision)
					actual_position = r
					distance_in_material = 0
				end if
				! Updating variables related to end condition for embedding due to 
				! random atomic electron collision 
				previous_position = actual_position
				actual_position = r
				step_length = norm2(actual_position - previous_position)
				distance_in_material = distance_in_material + step_length
				! Did the electron got embedded due to collision with atomic electron?
				if (distance_in_material .ge. distance_before_collision) then
					is_embedded = .true.
					num_embedded = num_embedded + 1
					embedded_positions(num_embedded,:) = r
					! Updating super electron arrays based on new embedded electron
					call update_super_electron_in_cell &
					(r, partition_boundaries, num_super_electrons, & 
					super_electron_charges, super_electron_positions)
					print*, "Trajectory end --> Electron is embedded"
					print*, "Total iterations:", i
					print*, "Final electron position:", r
					print*, " Distance before collision:", distance_before_collision
					print*, " Distance travelled inside material:", distance_in_material
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
			! End condition for total distance travelled
			if (distance_to_target .gt. initial_distance_to_target .and. &
			r(2) .gt. 0) then
				is_scattered = .true.
				num_scattered = num_scattered + 1
				scattered_positions(num_scattered,:) = r
				print*, "Trajectory end --> Electron is scattered"
				print*, "Total iterations:", i
				print*, "Final electron position:", r
				print*, " Initial distance to target:", initial_distance_to_target
				print*, " Final distance to target:", distance_to_target
				print*
			end if
			! End condition for maximum number of iterations
			if (i .ge. max_iterations) then
				is_max_iteration = .true.
				print*, "Trajectory end --> Maximum number of iterations reached"
				print*, "Total iterations:", i
				print*, "Final electron position:", r
				print*, " Maximum number of iterations:", max_iterations
				print*
			end if
		end do
	end subroutine compute_trajectory_optimized
	!=======================================================================
	! Subroutine : write_final_super_electron_distribution
	! Purpose    : Write the final distribution of super electrons, including
	!              their charges and positions, to a file.
	! Arguments  :
	!   - integer(i8), intent(in) :: partition_boundaries(3)
	!       Boundaries of the partition cells.
	!   - integer(i8), intent(in) :: num_super_electrons( &
	!         -partition_boundaries(1):, -partition_boundaries(2):, &
	!         -partition_boundaries(3):)
	!       Number of super electrons in each cell.
	!   - integer, intent(in) :: super_electron_charges( &
	!         -partition_boundaries(1):, -partition_boundaries(2):, &
	!         -partition_boundaries(3):, :)
	!       Charges of the super electrons in each cell.
	!   - real(dp), intent(in) :: super_electron_positions( &
	!         -partition_boundaries(1):, -partition_boundaries(2):, &
	!         -partition_boundaries(3):, :, :)
	!       Positions of the super electrons in each cell.
	! Variables  :
	!   - integer(i8) :: pbi, pbj, pbk
	!       Boundaries of the partition cells.
	!   - integer(i8) :: i, j, k, n
	!       Loop counters.
	!=======================================================================
	subroutine write_final_super_electron_distribution &
		(partition_boundaries, num_super_electrons, super_electron_charges, &
		super_electron_positions)
		implicit none
		integer(i8), intent(in) :: partition_boundaries(3)
		integer(i8), intent(in) :: num_super_electrons( &
			-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):)
		integer, intent(in) :: super_electron_charges( &
			-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :)
		real(dp), intent(in) :: super_electron_positions( &
			-partition_boundaries(1):, -partition_boundaries(2):, &
			-partition_boundaries(3):, :, :)
		integer(i8) :: pbi, pbj, pbk
		integer(i8) :: i, j, k, n
		! Getting partition cells boundaries
		pbi = partition_boundaries(1)
		pbj = partition_boundaries(2)
		pbk = partition_boundaries(3)
		open(unit=42, file='final_super_electron_distribution.dat', &
			status='replace', action='write')
			do i = -pbi, pbi
				do j = -1, -pbj, -1
					do k = -pbk, pbk
						! If there are more than zero super electrons in the i, j, k cell,
						! get the super electron(s) charges and positions
						if (num_super_electrons(i,j,k) .gt. 0) then
							do n = 1, num_super_electrons(i,j,k)
								write(42, *) i, j, k, n, super_electron_charges(i,j,k,n), &
								super_electron_positions(i,j,k,n,:)
							end do
						end if
					end do
				end do
			end do
		close(42)
	end subroutine write_final_super_electron_distribution
end module m4_optimized_trajectory_computation