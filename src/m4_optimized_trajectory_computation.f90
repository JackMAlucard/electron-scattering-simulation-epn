module m4_optimized_trajectory_computation
	implicit none
	use m0
	use m3!Specify which subroutines?
	contains

	! Subroutine that takes the position of a projectile electron in 
	! the material zone and outputs the indices to a nearby material atom.
	subroutine get_nearby_atom_indices(r, material_boundaries, atom_indices)
		implicit none
		real(dp), intent(in) :: r(3)	! Projectile electron position
		integer(i8), intent(in) :: material_boundaries
		integer(i8), intent(out) :: atom_indices(3)
		integer(i8) :: Nx, Ny, Nz, i, j, k
		! Getting material grid boundaries
		Nx = material_boundaries(1)
		Ny = material_boundaries(2)
		Nz = material_boundaries(3)
		if (r(1) .ge. 0) then
			i = dint(r(1)*INTERATOMIC_DIST_SIO2_INV) + 1
			if (i .gt. Nx) i = Nx
		else
			i = dint(r(1)*INTERATOMIC_DIST_SIO2_INV) - 1
			if (i .lt. -Nx) i = -Nx
		end if
		if (r(2) .ge. 0) then
			j = 0
		else
			j = dint(r(2)*INTERATOMIC_DIST_SIO2_INV) - 1
			if (j .lt. -Ny) j = -Ny
		end if
		if (r(3) .ge. 0) then
			k = dint(r(3)*INTERATOMIC_DIST_SIO2_INV) + 1
			if (k .gt. Nz) k = Nz
		else
			k = dint(r(3)*INTERATOMIC_DIST_SIO2_INV) - 1
			if (k .lt. -Nz) k = -Nz
		end if
		atom_indices = (/i, j, k/)
	end subroutine get_nearby_atom_indices
	
	! Subroutine that defines the cells indices, effectively partitioning the
	! material zone into these cells. It also allocates and initializes the 
	! arrays that store the super electrons information for all cells.
	subroutine setup_cells_and_super_electrons &
		(num_electrons, material_boundaries, cells_boundaries, &
		num_super_electrons_in_cell, super_electron_charges, & 
		super_electron_positions)
		implicit none
		integer(i8), intent(in) :: num_electrons, material_boundaries(3)
		integer(i8), intent(out) :: cells_boundaries(3)
		integer(i8), intent(out) :: num_super_electrons_in_cell(:,:,:)
		integer(i8), intent(out) :: super_electron_charges(:,:,:,:)
		real(dp), intent(out) :: super_electron_positions(:,:,:,:,:)
		integer(i8) :: NCx, NCy, NCz
		integer :: max_super_electrons
		! Computing partition cells boundaries
		NCx = material_boundaries(1)/CELL_SCALE_FACTOR + 1
		NCy = material_boundaries(2)/CELL_SCALE_FACTOR + 1
		NCz = material_boundaries(3)/CELL_SCALE_FACTOR + 1
		cells_boundaries = (/NCx, NCy, NCz/)
		! num_cells = (2*NCx)*(NCy)*(2*NCz)
		! Maximum possible number of super electrons on a cell, considering the
		! extreme case in which all of the electrons are embedded into a single cell
		max_super_electrons = (num_electrons - 1)/MAX_EQUIVALENT_CHARGE + 1
		! Allocating array that stores the number of super electrons on all cells
		allocate (num_super_electrons_in_cell(-NCx:NCx,-NCy:-1,-NCz:NCz))
		! Initializing as zero, since there are no super electrons at the start
		num_super_electrons_in_cell = 0
		! Allocating array that stores the super electrons charges on all cells
		allocate &
		(super_electron_charges(-NCx:NCx,-NCy:-1,-NCz:NCz,max_super_electrons))
		! Initializing as zero, since there are no super electrons at the start
		super_electron_charges = 0
		! Allocating array that stores the super electrons positions on all cells
		allocate &
		(super_electron_positions(-NCx:NCx,-NCy:-1,-NCz:NCz,max_super_electrons,3))
		! Initializing as zero, the positions will be given as the mean of the 
		! positions of the electrons that form the super electron
		super_electron_positions = 0
	end subroutine setup_cells_and_super_electrons

	! Subroutine that takes the position of an electron r(3), the inverse cell 
	! length (in a0) and determines the cell_indices of the cubic cell to which 
	! the embedded electron belongs to in order to update the super electron in 
	! that cell. 
	! Electrons embedded after crossing the material boundaries will be included
	! in the cell at the boundary of the partitions.
	! The lowest possible index in the y-direction is -1 since there wouldn't be
	! any cells with index 0 in any direction. Since the arrays have those slots,
	! those MUST BE ignored, and as such were initialized and will remain zero.
	subroutine get_cell_indices(r, cells_boundaries, cell_indices)
		implicit none
		real(dp), intent(in) :: r(3)
		integer(i8), intent(in) :: cells_boundaries(3)
		integer(i8), intent(out) :: cell_indices(3)
		integer(i8) :: NCx, NCy, NCz
		integer(i8) :: i, j, k
		! Getting partition cells boundaries
		NCx = cells_boundaries(1)
		NCy = cells_boundaries(2)
		NCz = cells_boundaries(3)
		if (r(1) .ge. 0) then
			i = dint(r(1)*CELL_LENGTH_INV) + 1
			if (i .gt. NCx) i = NCx
		else
			i = dint(r(1)*CELL_LENGTH_INV) - 1
			if (i .lt. -NCx) i = -NCx
		end if
		if (r(2) .ge. 0) then
			j = -1 ! The top cell includes a bit more space
		else
			j = dint(r(2)*CELL_LENGTH_INV) - 1
			if (j .lt. -NCy) j = -NCy
		end if
		if (r(3) .ge. 0) then
			k = dint(r(3)*CELL_LENGTH_INV) + 1
			if (k .gt. NCz) k = NCz
		else
			k = dint(r(3)*CELL_LENGTH_INV) - 1
			if (k .lt. -NCz) k = -NCz
		end if
		cell_indices = (/i, j, k/)
	end subroutine get_cell_indices
	
	! Subroutine that updated the values of the super electron arrays in the
	! cell in which the projectile electron with position r gets embedded.
	subroutine update_super_electron_in_cell &
		(r, cells_boundaries, num_super_electrons, & 
		super_electron_charges, super_electron_positions)
		implicit none
		real(dp), intent(in) :: r(3)
		integer(i8), intent(in) :: cells_boundaries(3)
		integer(i8), intent(inout) :: num_super_electrons(:,:,:)
		integer, intent(inout) :: super_electron_charges(:,:,:,:)
		real(dp), intent(inout) :: super_electron_positions(:,:,:,:,:)
		integer(i8) :: cell_indices(3), i, j, k
		integer(i8) :: super_electron_num
		integer :: super_electron_charge
		real(dp) :: super_electron_r
		! Getting embedded electron's cell based on its position r
		call get_cell_indices(r, cells_boundaries, cell_indices)
		i = cell_indices(1)
		j = cell_indices(2)
		k = cell_indices(3)		
		! Getting number of super electrons in the cell
		super_electron_num = num_super_electrons(i,j,k)
		! Getting current super electron charge
		super_electron_charge = super_electron_charges(i,j,k,n)
		! Getting current super electron position
		super_electron_r = super_electron_positions(i,j,k,n,:)
		! If there are no embedded electrons in the cell, yet,
		! i.e. no super electrons yet in the cell
		if (super_electron_num .eq. 0) then
			! Updating super electrons arrays initializing the first one 
			num_super_electrons(i,j,k) = 1
			super_electron_charges(i,j,k,super_electron_num) = 1
			super_electron_positions(i,j,k,super_electron_num,:) = r
		! If current super electron is not full of charge, add extra charge
		! and compute new super electron position
		else if (super_electron_charge .lt. MAX_EQUIVALENT_CHARGE) then
			super_electron_charge = super_electron_charge + 1
			super_electron_r = &
			super_electron_r + (r - super_electron_r)/super_electron_charge
			! Updating super electron arrays after adding a new embedded electron
			super_electron_charges(i,j,k,super_electron_num) = super_electron_charge
			super_electron_positions(i,j,k,super_electron_num,:) = super_electron_r
		! If current super electron is full of charge, initialize next one
		else if (super_electron_charge .eq. MAX_EQUIVALENT_CHARGE) then
			super_electron_num = super_electron_num + 1
			! Updating super electrons arrays initializing the next one
			num_super_electrons(i,j,k) = super_electron_num
			super_electron_charges(i,j,k,super_electron_num) = 1
			super_electron_positions(i,j,k,super_electron_num,:) = r
		else
			!Everything else is an error (TB REMOVED)
			print*, "ERROR"
			print*, "Cell indices:", i,j,k
			print*, "Super electron number:", super_electron_num
			print*, "Super electron charge:", super_electron_charge
			print*, "Super electron position:", super_electron_r
		end if
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
		atom_charges, dt, r, v, a)
		implicit none
		integer(i8), intent(in) :: num_embedded, material_boundaries(3)
		real(dp), intent(in) :: embedded_positions(:,:), atom_positions(:,:,:,:)
		integer, intent(in) :: , atom_charges(:,:,:)
		real(dp), intent(in) :: dt
		real(dp), intent(inout) :: r(3), v(3), a(3)
		real(dp) :: rt(3), ak(3)
		integer :: Z
		integer(i8) :: atom_indices(3)
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
			call get_nearby_atom_indices(r, material_boundaries, atom_indices)
			do i = atom_indices(1) - 1, atom_indices(1) + 1
				do j = atom_indices(2) - 1, atom_indices(2) + 1
					do k = atom_indices(3) - 1, atom_indices(3) + 1
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
	end subroutine time_step_near_zone

	! TO BE EDITED
	! Velocity Verlet Step for Far Field (bulk electrons)
	! Subroutine which computes a time step of the Velocity Verlet algorithm used to
	! compute the trajectory of a projectile electron interacting with Nb equivalent
	! embedded electrons characterized by the b_neq, b_ec, and b_er arrays. 
	subroutine time_step_far_zone &
		(num_embedded, cells_boundaries, num_super_electrons, &
		super_electron_charges, super_electron_positions, material_boundaries, &
		atom_positions, atom_charges, dt, r, v, a)
		implicit none
		integer(i8), intent(in) :: num_embedded, num_super_electrons(:,:,:)
		integer(i8), intent(in) :: cells_boundaries(3), material_boundaries(3)
		integer, intent(in) :: super_electron_charges(:,:,:,:), atom_charges(:,:,:)
		real(dp), intent(in) :: super_electron_positions(:,:,:,:,:)
		real(dp), intent(in) :: atom_positions(:,:,:,:)
		real(dp), intent(in) :: dt
		real(dp), intent(inout) :: r(3), v(3), a(3)
		real(dp) :: rt(3), ak(3)
		integer :: Q, Z
		integer(i8) :: NCx, NCy, NCz, super_electron_num
		integer(i8) :: i, j, k, n
		! Half-step velocity update
		v = v + 0.5*a*dt
		! Position update
		r = r + v*dt
		! Full-step acceleration update
		a = 0
		! Getting partition cells boundaries
		NCx = cells_boundaries(1)
		NCy = cells_boundaries(2)
		NCz = cells_boundaries(3)
		! Acceleration due to embedded super electrons
		! ONLY computed if there is AT LEAST one embedded electron
		if (num_embedded .gt. 0) then
			do i = -NCx, NCx
				! To avoid checking for unused cells with index 0
				if (i .ne. 0) then
					do j = -1,-NCy, -1
						do k = -NCz, NCz
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
		! IT DOESN'T LOOK LIKE THIS NEXT PART IS NECESSARY, MUST TEST!!!!
		! Acceleration due to material atoms
		! ONLY computed if near the material zone, i.e. if the electron 
		! projectile position y-coordinate is less than the material's height
		if (r(2) .lt. MATERIAL_HEIGHT_SIO2) then
			call get_nearby_atom_indices(r, material_boundaries, atom_indices)
			do i = atom_indices(1) - 1, atom_indices(1) + 1
				do j = atom_indices(2) - 1, atom_indices(2) + 1
					do k = atom_indices(3) - 1, atom_indices(3) + 1
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
	end subroutine time_step_far_zone

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
subroutine trajectory_opt(i,N,P,ou,Ne,Ns,emb,sct,Nx,Ny,Nz,atoms,Z,Nxb,Nyb,Nzb,max_eq,max_neq,b_neq,b_ec,b_er,d0,FFed,T,d,dt,l,r,v,a)
!REGULAR TRAJECTORY VARIABLES
integer(i8), intent(in) :: i, N, P, ou, Nx, Ny, Nz, T
integer, intent(in) :: Z(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1)
real(dp), intent(in) :: atoms(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1,3)
real(dp), intent(in) :: d0, d, dt, l
!MIXED METHOD VARIABLES, WHICH USE VV_STEP_FF
integer(i8), intent(in) :: Nxb, Nyb, Nzb !Bins boundary indices
integer(i8), intent(in) :: max_eq		!Maximum equivalent charge, 180 from tests (max 182)
integer(i8), intent(in) :: max_neq	!Maximum number of equivalent electrons per bin, N_beam/max_eq+1
integer(i8), intent(inout) :: b_neq(-Nxb:Nxb,-Nyb:-1,-Nzb:Nzb)!Current equivalent electron array
integer(i8), intent(inout) :: b_ec(-Nxb:Nxb,-Nyb:-1,-Nzb:Nzb,0:max_neq)!Equivalent charges array
real(dp), intent(inout) :: b_er(-Nxb:Nxb,-Nyb:-1,-Nzb:Nzb,0:max_neq,3)!Equivalent positions array
real(dp), intent(in) :: FFed	!Far Field effective distance
integer(i8), intent(inout) :: Ne, Ns
real(dp), intent(inout) :: emb(N,3), sct(N,3)
real(dp), intent(inout) :: r(N,3), v(N,3), a(N,3)
real(dp) :: dk
integer(i8) :: j, k
!Internal variables for equivalent electron bin sorting
integer(i8) :: ib,jb,kb
integer :: bin_size
real(dp) :: d_bin, d_bin_inv
logical :: end_embedded, end_scattered, end_iterations	!Trajectory end conditions
logical :: in_material
real(dp) :: r_mat_i(3), r_mat_f(3), dr_mat(3), d_mat, d_embedded
real(dp) :: tk
	!Setting bin size (distance) and computing its inverse for sort_electron
	!subroutine (THIS SHOULD GO ON MAIN)
	bin_size = 3	!MUST CHECK IF IT IS THE SAME VALUE IN main
	d_bin = bin_size*d
	d_bin_inv = 1/d_bin
	!Initialize end conditions as false
	end_embedded = .false.
	end_scattered = .false.
	end_iterations = .false.
	!Initialize variables related to end conditions
	k = 0
	dk = d0
	in_material = .false.
	r_mat_f = 0
	d_mat = 0
	!Initialize points to be plotted counter
	j = 0
	do while ( .not.(end_embedded .or. end_scattered .or. end_iterations) )
		tk = k*dt !t0 = 0, just to print to file, not needed for computations
		!Plotting only P points
		if ( (mod(k,T/P) .eq. 0) .and. j .lt. P ) then
			!Write values to file
			write(ou,*) tk, r(i,:)	!I'LL HAVE TO READ FROM FILE TO COMPARE ON TEST!
			j = j + 1
		end if
		!Computing next Velocity Verlet step
		!Check if the electron distance to center is greater than the far field effective distance ffed (a0)
		if (dk .gt. FFed .and. r(i,2) .gt. 0._dp) then
			!Computing next Velocity Verlet step using Far Field Approximation
			call vv_step_FF(Ne,Nxb,Nyb,Nzb,max_eq,max_neq,b_neq,b_ec,b_er,Nx,Ny,Nz,atoms,Z,d,dt,r(i,:), v(i,:), a(i,:))
		else
			!Computing next Velocity Verlet step for Near Field
			call vv_step(Ne, N, emb, Nx, Ny, Nz, atoms, Z, d, dt, r(i,:), v(i,:), a(i,:))
		end if
		!Updating variables related to end conditions for distance and iterations
		dk = norm2(r(i,:))
		k = k + 1
		!Check if any of the end conditions are met
		!End condition for embedded electrons
		!The electron must be below material level
		if (r(i,2) .lt. 0._dp) then
			!Updating variables related to random embedding due to attenuation
			r_mat_i = r_mat_f
			r_mat_f = r(i,:)
			dr_mat = r_mat_f - r_mat_i
			d_mat = d_mat + norm2(dr_mat)
			!ONLY the first time it enters the material
			!Generate random number following exponential distribution
			!Save first electron position inside material
			!Initialize first electron distance travelled inside material
			!This overwrites the update in variables related to random embedding  
			if (.not.(in_material)) then
				in_material = .true.
				call random_exponential(d_embedded,l)
				r_mat_i = r(i,:)
				d_mat = 0
			end if
			!Did the electron got embedded due to attenuation?
			if (d_mat .ge. d_embedded) then
				end_embedded = .true.
				print*, "End condition for embedded electron due to attenuation is met"
				Ne = Ne + 1
				emb(Ne,:) = r(i,:)
				!Update equivalent electron bins
				!Sort embedded electron into their respective equivalent electron bin
				call sort_electron(emb(Ne,:),d_bin_inv,Nxb,Nyb,Nzb,ib,jb,kb)
				!Add embedded electron into each respective bin
				call add_electron_to_bin(ib,jb,kb,emb(Ne,:),Nxb,Nyb,Nzb,max_eq,max_neq,b_neq,b_ec,b_er)
				!Print embedded electron info message to console
				print*, "Final electron position:", r(i,:)
				print*, "Attenuation distance (random exp):", d_embedded
				print*, "Distance travelled inside material:", d_mat
				print*, "x_scc: +/-", atoms(Nx,0,0,1)
				print*, "y_scc:", atoms(0,-Ny,0,2)
				print*, "z_scc: +/-", atoms(0,0,Nz,3)
				print*, "k", k
				!Print equivalent electron info message to console
				print*, "Equivalent electron info:"
				print*, "i,j,k:", ib,jb,kb
				print*, "neq:", b_neq(ib,jb,kb)
				print*, "ec:", b_ec(ib,jb,kb,b_neq(ib,jb,kb))
				print*, "er:", b_er(ib,jb,kb,b_neq(ib,jb,kb),:)
				print*
			end if
		end if
		!End condition for distance
		if (dk .gt. d0 .and. r(i,2) .gt. 0._dp) then
			end_scattered = .true.
			Ns = Ns + 1
			sct(Ns,:) = r(i,:)
			print*, "End condition for scattered electron is met"
			print*, "d0:", d0
			print*, "dk:", dk
			print*, "k", k
		end if
		!End condition for iterations
		if (k .ge. T) then
			end_iterations = .true.
			print*, "End condition for iterations is met"
			print*, "Final electron position:", r(i,:)
			print*, "T:", T
			print*, "k", k
		end if
	end do
end subroutine trajectory_opt

end module m4_optimized_trajectory_computation