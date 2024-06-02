module m4_optimized_trajectory_computation
	use m0_utilities, &
		only: dp, i8, INTERATOMIC_DIST_SIO2_INV, CELL_SCALE_FACTOR, &
		CELL_LENGTH_INV, MAX_EQUIVALENT_CHARGE, MATERIAL_HEIGHT_SIO2, &
		CROSS_SECTION_SIO2, EFFECTIVE_DISTANCE, random_exponential
	use m3_trajectory_computation, &
		only: acceleration_due_to_electron, acceleration_due_to_atom
	implicit none
	contains

	! Subroutine that takes the position of a projectile electron in 
	! the material zone and outputs the indexes to a nearby material atom.
	subroutine get_nearby_atom_indexes(r, material_boundaries, atom_indexes)
		implicit none
		real(dp), intent(in) :: r(3)	! Projectile electron position
		integer(i8), intent(in) :: material_boundaries(3)
		integer(i8), intent(out) :: atom_indexes(3)
		integer(i8) :: mbi, mbj, mbk, i, j, k
		! Getting material grid boundaries
		mbi = material_boundaries(1)
		mbj = material_boundaries(2)
		mbk = material_boundaries(3)
		if (r(1) .ge. 0) then
			i = dint(r(1)*INTERATOMIC_DIST_SIO2_INV) + 1
			if (i .gt. mbi) i = mbi
		else
			i = dint(r(1)*INTERATOMIC_DIST_SIO2_INV) - 1
			if (i .lt. -mbi) i = -mbi
		end if
		if (r(2) .ge. 0) then
			j = 0
		else
			j = dint(r(2)*INTERATOMIC_DIST_SIO2_INV) - 1
			if (j .lt. -mbj) j = -mbj
		end if
		if (r(3) .ge. 0) then
			k = dint(r(3)*INTERATOMIC_DIST_SIO2_INV) + 1
			if (k .gt. mbk) k = mbk
		else
			k = dint(r(3)*INTERATOMIC_DIST_SIO2_INV) - 1
			if (k .lt. -mbk) k = -mbk
		end if
		atom_indexes = (/i, j, k/)
	end subroutine get_nearby_atom_indexes
	
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
	
	! TO BE EDITED
	! Subroutine which computes a time step of the Velocity Verlet algorithm
	! used to compute the trajectory of a projectile electron interacting with N_e
	! embedded electrons out of N_et total possible embedded electrons and a **SCC**
	! bulk of N_a=(2*Nx + 1)*(Ny + 1)*(2*Nz + 1) neutral atoms with atomic number Z.
	! The position of the embedded electrons is stored in the array e_emb(N_et,3).
	! Acceleration due to embedded electrons is computed for EVERY embedded electron
	! and ONLY if there is AT LEAST one embedded electron.
	! Acceleration due to neutral atoms is ONLY computed if near the material zone, 
	! i.e. its height is less than half the interatomic distance and lower.
	! This subroutine only computes the acceleration due to nearby neutral atoms,
	! in comparison with the original vv_step, which computed the acceleration from
	! all neutral atoms in the bulk. This is a good approximation and significantly 
	! speeds up the simulation.
	! This subroutine is used to compute the steps when near the material zone,
	! which is why it is called Near Field (NF).
	subroutine time_step_near_zone &
		(num_embedded, material_boundaries, embedded_positions, atom_positions, &
		atom_charges, atom_charges_cbrt, dt, r, v, a)
		implicit none
		integer(i8), intent(in) :: num_embedded, material_boundaries(3)
		real(dp), intent(in) :: embedded_positions(:,:)
		real(dp), intent(in) :: atom_positions(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:, :)
		integer, intent(in) :: atom_charges(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:)
		real(dp), intent(in) :: atom_charges_cbrt(-material_boundaries(1)-1:, &
			-material_boundaries(2)-1:, -material_boundaries(3)-1:)
		real(dp), intent(in) :: dt
		real(dp), intent(inout) :: r(3), v(3), a(3)
		real(dp) :: rt(3), ak(3), cbrt_Z
		integer :: Z
		integer(i8) :: atom_indexes(3)
		integer(i8) :: i, j, k
		! Half-step velocity update
		v = v + 0.5*a*dt
		! Position update
		r = r + v*dt
		! Full-step acceleration update
		a = 0
		! Acceleration due to embedded electrons
		! ONLY computed if there is AT LEAST one embedded electron
		if (num_embedded .gt. 0) then
			do k = 1, num_embedded
				rt = embedded_positions(k,:)
				call acceleration_due_to_electron(r, rt, ak)
				a = a + ak
			end do
		end if
		! Acceleration due to material atoms
		! ONLY computed if near the material zone, i.e. if the electron 
		! projectile position y-coordinate is less than the material's height
		if (r(2) .lt. MATERIAL_HEIGHT_SIO2) then
			call get_nearby_atom_indexes(r, material_boundaries, atom_indexes)
			do i = atom_indexes(1) - 1, atom_indexes(1) + 1
				do j = atom_indexes(2) - 1, atom_indexes(2) + 1
					do k = atom_indexes(3) - 1, atom_indexes(3) + 1
						!Only compute acceleration in positions with atoms
						Z = atom_charges(i,j,k)
						if (Z .ne. 0) then
							rt = atom_positions(i,j,k,:)
							cbrt_Z = atom_charges_cbrt(i,j,k)
							call acceleration_due_to_atom(r, rt, Z, cbrt_Z, ak)
							a = a + ak
						end if
					end do
				end do
			end do
		end if
		! Second half-step velocity update
		v = v + 0.5*a*dt
	end subroutine time_step_near_zone

	! TO BE EDITED
	! Velocity Verlet Step for Far Field (bulk electrons)
	! Subroutine which computes a time step of the Velocity Verlet algorithm used to
	! compute the trajectory of a projectile electron interacting with Nb equivalent
	! embedded electrons characterized by the b_neq, b_ec, and b_er arrays. 
	subroutine time_step_far_zone &
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
	end subroutine time_step_far_zone

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
		(num_plot_ploints, max_iterations, output_unit, material_boundaries, & 
		atom_positions, atom_charges, atom_charges_cbrt, partition_boundaries, & 
		num_super_electrons, super_electron_positions, super_electron_charges, & 
		dt, r, v, a, num_embedded, num_scattered, embedded_positions, & 
		scattered_positions)
		implicit none
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
			if ((mod(i,max_iterations/num_plot_ploints) .eq. 0) .and. &
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
				call time_step_far_zone &
				(num_embedded, partition_boundaries, num_super_electrons, &
				super_electron_charges, super_electron_positions, dt, r, v, a)
			else
				! Computing next Velocity Verlet time step integration using near zone
				! approximation
				call time_step_near_zone &
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

end module m4_optimized_trajectory_computation