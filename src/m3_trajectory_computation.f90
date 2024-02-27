module m3_trajectory_computation
	implicit none
	use m0!, only ...
	contains

	! Acceleration of a projectile Electron interacting with a stationary target 
	! Electron via the electrostatic Coulomb potential. 
	! The input and output variables are in the atomic unit system (au).
	! rp: r projectile, position of the incident electron (a0).
	! rt: r target, position of the stationary electron (a0).
	! a: Acceleration of the incident electron (...)
	subroutine acceleration_due_to_electron(rp, rt, a)
		implicit none
		real(dp), intent(in) :: rp(3), rt(3)
		real(dp), intent(out) :: a(3)
		real(dp) :: rs(3), r !rs, r: separation vector and its magnitude
		rs = rp - rt
		r = norm2(rs)
		a = rs/(r**3)
	end subroutine acceleration_due_to_electron

	! Acceleration of a projectile Electron interacting with a stationary neutral 
	! Atom of atomic number Z via the electrostatic Thomas-Fermi potential.
	! The input and output are in atomic units system (au).
	! rp: r projectile, position of the incident electron (a0).
	! rt: r target, position of the stationary neutral atom (a0).
	! Z: Atomic number of the stationary neutral atom.
	! cbrt_Z: Cube root of the atomic number of the stationary neutral atom.
	! a: Acceleration of the incident electron (...)
	subroutine acceleration_due_to_atom(rp, rt, Z, cbrt_Z, a)
		implicit none
		real(dp), intent(in) :: rp(3), rt(3)
		integer, intent(in) :: Z, cbrt_Z
		real(dp), intent(out) :: a(3)
		real(dp), parameter :: a_coef(3) = (/0.10_dp, 0.55_dp, 0.35_dp/)
		real(dp), parameter :: b_coef(3) = (/6.00_dp, 1.20_dp, 0.30_dp/)
		real(dp), parameter :: b0_inv = ((2**7)/(3*PI)**2)**(1/3._dp)!b0_inv = 1/b0
		!real(dp), parameter :: b0_inv = 1.1295078101832_dp !b0_inv = 1/b0
		real(dp) :: rs(3), r !rs, r: separation vector and its magnitude
		real(dp) :: chi, psi, aux
		integer(i8) :: i
		! Electron-atom relative position
		rs = rp - rt
		r = norm2(rs)
		! Computing chi and psi
		chi = 0
		psi = 0
		do i = 1, 3
			aux = a_coef(i)*dexp(-b_coef(i)*r*b0_inv*cbrt_Z)
			chi = chi + aux
			aux = aux*b_coef(i)
			psi = psi + aux
		end do
		a = -(Z/r**3)*(chi + psi*b0_inv*cbrt_Z*r)*rs ! CHECK THAT IT IS CORRECT
	!	a = -Z*( (chi/(r**3)) + ((psi*b0_inv*cbrt_Z)/(r**2)) )*rs
	end subroutine acceleration_due_to_atom
	
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
	! CONSIDER GIVING THIS ONE AN ADJECTIVE
	! subroutine time_step(N_e, N_et, emb, Nx, Ny, Nz, atoms, Z, d, dt, r, v, a)
	subroutine time_step &
		(num_embedded, material_boundaries, embedded_positions, atom_positions, &
		atom_charges, dt, r, v, a)
		implicit none
		integer(i8), intent(in) :: num_embedded, material_boundaries(3)
		real(dp), intent(in) :: embedded_positions(:,:), atom_positions(:,:,:,:)
		integer, intent(in) :: , atom_charges(:,:,:)
		real(dp), intent(in) :: dt
		real(dp), intent(inout) :: r(3), v(3), a(3)
		real(dp) :: rt(3), ak(3)
		integer :: Z
		integer(i8) :: Nx, Ny, Nz
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
		! Getting material grid boundaries
		Nx = material_boundaries(1)
		Ny = material_boundaries(2)
		Nz = material_boundaries(3)
		! Acceleration due to material atoms
		! ONLY computed if near the material zone, i.e. if the electron 
		! projectile position y-coordinate is less than the material's height
		if (r(2) .lt. MATERIAL_HEIGHT_SIO2) then
			do i = -Nx, Nx
				do j = -Ny, 0, -1
					do k = -Nz, Nz
						!Only compute acceleration in positions with atoms
						Z = atom_charges(i,j,k)
						if (Z .ne. 0) then
							rt = atom_positions(i,j,k,:)
							call acceleration_due_to_atom(r, rt, Z, ak)
							a = a + ak
						end if
					end do
				end do
			end do
		end if
		! Second half-step velocity update
		v = v + 0.5*a*dt
	end subroutine time_step

	! Subroutine which computes the trajectory of the i-th projectile electron, out
	! of N electrons in the electron beam, and plots P points in unit file number ou.
	! This subroutine uses the Near Field Regimen, i.e. computes the field of each
	! invidivual electron regardless of distance via the vv_step subroutine.
	! The projectile interacts with a set of N_e embedded electrons out of a maximum
	! of N total possible embedded electrons and a **SCC** bulk of 
	! N_a=(2*Nx + 1)*(Ny + 1)*(2*Nz + 1) neutral atoms with atomic number Z.
	! The beam initial positions, velocities, and accelerations are stored in the 
	! (N,3) arrarys r, v and a respectively.
	! The embedded and scattered electrons final positions are stored in the (N,3)
	! arrays e_emb and e_sct respectively. 
	! The simulation is ended either when the electron distance d to the origin of 
	! coordinates (which approximates the center of the bulk) is greater than the 
	! initial distance d0, or when the number of iterations k is equal or greater 
	! than the estimated number of iterations T (how is this estimated?)
	! The simulation is also ended if the electron gets embedded on the bulk, at
	! which point the number of embedded electrons N_e is updated/increased.
	! I NEED to explain WHY these conditions for the end of the trajectory are set.
	subroutine compute_trajectory &
		(num_plot_ploints, max_iterations, output_unit, &
		material_boundaries, atom_positions, atom_charges, &
		initial_distance_to_target, dt, r, v, a &
		num_embedded, num_scattered, embedded_positions, scattered_positions)
		implicit none
		integer(i8), intent(in) :: num_plot_ploints, max_iterations
		integer(i8), intent(in) :: output_unit, material_boundaries(3)
		real(dp), intent(in) :: atom_positions(:,:,:,)
		integer, intent(in) :: atom_charges(:,:,:)
		real(dp), intent(in) :: initial_distance_to_target, dt
		real(dp), intent(inout) :: r(3), v(3), a(3)
		integer(i8), intent(inout) :: num_embedded, num_scattered
		real(dp), intent(inout) :: embedded_positions(:,:), scattered_positions(:,:)
		logical :: is_embedded, is_scattered, is_max_iteration
		logical :: in_material
		real(dp) :: distance_before_collision, distance_in_material
		real(dp) :: previous_position(3), actual_position(3), step_length
		real(dp) :: distance_to_target, t
		integer(i8) :: Nx, Ny, Nz
		integer(i8) :: i, j
		! Initializing end conditions as false
		is_embedded = .false.
		is_scattered = .false.
		is_max_iteration = .false.
		! Initializing variables related to end conditions
		i = 0
		in_material = .false.
		distance_to_target = initial_distance_to_target
		! Initialize points to be plotted counter
		j = 0
		do while (.not.(is_embedded .or. is_scattered .or. is_iterations))
			t = i*dt ! t0 = 0, just to print to file, not needed for computations
			! Plotting only P points
			if ((mod(k,max_iterations/num_plot_ploints) .eq. 0) .and. &
			j .lt. num_plot_ploints) then
				! Write values to file
				write(output_unit,*) t, r
				j = j + 1
			end if
			! Computing next Velocity Verlet time step integration
			call time_step &
			(num_embedded, material_boundaries, embedded_positions, &
			atom_positions, atom_charges, dt, r, v, a)
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
					print*, "Trajectory end --> Electron is embedded"
					print*, "Total iterations:", i
					print*, "Final electron position:", r
					print*, " Distance before collision:", distance_before_collision
					print*, " Distance travelled inside material:", distance_in_material
					print*, " Material boundaries:"
					Nx = material_boundaries(1)
					Ny = material_boundaries(2)
					Nz = material_boundaries(3)
					print*, "  x --> +/-", atom_positions(Nx,0,0,1)
					print*, "  y -->    ", atom_positions(0,-Ny,0,2)
					print*, "  z --> +/-", atom_positions(0,0,Nz,3)
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
	end subroutine compute_trajectory

end module m3_trajectory_computation