module m4
use m0
use m3!Specify which subroutines?
implicit none

contains

!Subroutine to sort electron inside the material zone
subroutine sort_in_material(r,d,Nx,Ny,Nz,i,j,k)
real(dp), intent(in) :: r(3), d	!Electron position and interatomic distance
integer(i8), intent(in) :: Nx, Ny, Nz	!Material bound indices
integer(i8), intent(out) :: i, j, k
	if (r(1) .ge. 0._dp) then
		i = dint(r(1)/d) + 1
		if (i .gt. Nx) i = Nx
	else
		i = dint(r(1)/d) - 1
		if (i .lt. -Nx) i = -Nx
	end if
	if (r(2) .ge. 0._dp) then
		j = 0
	else
		j = dint(r(2)/d) - 1
		if (j .lt. -Ny) j = -Ny
	end if
	if (r(3) .ge. 0._dp) then
		k = dint(r(3)/d) + 1
		if (k .gt. Nz) k = Nz
	else
		k = dint(r(3)/d) - 1
		if (k .lt. -Nz) k = -Nz
	end if
end subroutine sort_in_material

!Subroutine that takes the position of an electron r(3), the inverse of the 
!bin "interatomic" distance d_bin_inv in a0 and determines the identifying indices i,j,k 
!of the cubic bin to which the electron belongs to in order to group electrons 
!together into an equivalent electron with charge greater or equal than one.
!Included failsafe that avoids exceeding the bin region boundary indices.
subroutine sort_electron(r,d_bin_inv,Nxb,Nyb,Nzb,i,j,k)
real(dp), intent(in) :: r(3), d_bin_inv	!Electron position and interatomic distance
integer(i8), intent(in) :: Nxb, Nyb, Nzb
integer(i8), intent(out) :: i, j, k
	if (r(1) .ge. 0._dp) then
		i = dint(r(1)*d_bin_inv) + 1
		if (i .gt. Nxb) i = Nxb
	else
		i = dint(r(1)*d_bin_inv) - 1
		if (i .lt. -Nxb) i = -Nxb
	end if
	if (r(2) .ge. 0._dp) then
		j = dint(r(2)*d_bin_inv) + 1
		if (j .gt. -1) j = -1
	else
		j = dint(r(2)*d_bin_inv) - 1
		if (j .lt. -Nyb) j = -Nyb
	end if
	if (r(3) .ge. 0._dp) then
		k = dint(r(3)*d_bin_inv) + 1
		if (k .gt. Nzb) k = Nzb
	else
		k = dint(r(3)*d_bin_inv) - 1
		if (k .lt. -Nzb) k = -Nzb
	end if
end subroutine sort_electron

!Subroutine that takes an electron position and equivalent electron bin indices,
!and updates the equivalent electron bin arrays correspondingly.
!There is no failsafe in the off case that the number of charges for a bin 
!exceeds the value of max_neq, however it may not be necessary to add one by 
!definition.
subroutine add_electron_to_bin(i,j,k,r,Nxb,Nyb,Nzb,max_eq,max_neq,b_neq,b_ec,b_er)
integer(i8), intent(in) :: i, j, k	!Electron bin indices
real(dp), intent(in) :: r(3)	!Electron position
integer(i8), intent(in) :: Nxb, Nyb, Nzb !Bins boundary indices
integer(i8), intent(in) :: max_eq		!Maximum equivalent charge, 180 from tests (max 182)
integer(i8), intent(in) :: max_neq	!Maximum number of equivalent electrons per bin, N_beam/max_eq+1
integer(i8), intent(inout) :: b_neq(-Nxb:Nxb,-Nyb:-1,-Nzb:Nzb)!Current equivalent electron array
integer(i8), intent(inout) :: b_ec(-Nxb:Nxb,-Nyb:-1,-Nzb:Nzb,0:max_neq)!Equivalent charges array
real(dp), intent(inout) :: b_er(-Nxb:Nxb,-Nyb:-1,-Nzb:Nzb,0:max_neq,3)!Equivalent positions array
integer(i8) :: neq, ec
real(dp) :: er(3)
	!Read number of equivalent electrons in the ijk bin
	neq = b_neq(i,j,k)
	!Read charge of current equivalent electron in the ijk bin
	ec = b_ec(i,j,k,neq)
	!Read equivalent position of equivalent electron in the ijk bin
	er = b_er(i,j,k,neq,:)
	!If current equivalent electron is not full of charge, add new charge
	!And compute new equivalent position
	if ((ec .ge. 0) .and. (ec .lt. max_eq)) then
		ec = ec + 1
		b_ec(i,j,k,neq) = ec
		er = er + (r-er)/ec!Running mean computed so as to minimize precision errors
		b_er(i,j,k,neq,:) = er
	!If current equivalent electron is full, move to next equivalent electron
	!And add new charge there
	else if (ec .eq. max_eq) then
		neq = neq + 1
		b_neq(i,j,k) = neq
		b_ec(i,j,k,neq) = 1
		b_er(i,j,k,neq,:) = r
	else
	!Everything else is an error
		print*, "ERROR"
		print*, "i,j,k:", i,j,k
		print*, "neq:", neq
		print*, "ec:", ec
		print*, "er", er
	end if
end subroutine add_electron_to_bin

!Acceleration of a projectile Electron interacting with a stationary target 
!equivalent electron of charge q via the electrostatic Coulomb potential. 
!The input and output variables are in atomic units (au).
!rp: r projectile, position of the incident electron (a0)
!rt: r target, equivalent position of the stationary electron (a0)
!a: Acceleration of the incident electron (aua)
subroutine akEE_FF(q, rp, rt, a)
integer(i8), intent(in) :: q
real(dp), intent(in) :: rp(3), rt(3)
real(dp), intent(out) :: a(3)
real(dp) :: rs(3), r !rs: r separation vector and its magnitude
	rs = rp - rt
	r = norm2(rs)
	a = (q*rs)/(r**3)
end subroutine akEE_FF

!Subroutine which computes a time step of the Velocity Verlet algorithm
!used to compute the trajectory of a projectile electron interacting with N_e
!embedded electrons out of N_et total possible embedded electrons and a **SCC**
!bulk of N_a=(2*Nx + 1)*(Ny + 1)*(2*Nz + 1) neutral atoms with atomic number Z.
!The position of the embedded electrons is stored in the array e_emb(N_et,3).
!Acceleration due to embedded electrons is computed for EVERY embedded electron
!and ONLY if there is AT LEAST one embedded electron.
!Acceleration due to neutral atoms is ONLY computed if near the material zone, 
!i.e. its height is less than half the interatomic distance and lower.
!This subroutine only computes the acceleration due to nearby neutral atoms,
!in comparison with the original vv_step, which computed the acceleration from
!all neutral atoms in the bulk. This is a good approximation and significantly 
!speeds up the simulation.
!This subroutine is used to compute the steps when near the material zone,
!which is why it is called Near Field (NF).
subroutine vv_step_NF(N_e, N_et, emb, Nx, Ny, Nz, atoms, Z, d, dt, r, v, a)
integer(i8), intent(in) :: N_e, N_et, Nx, Ny, Nz
integer, intent(in) :: Z(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1)
real(dp), intent(in):: emb(N_et,3), atoms(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1,3)
real(dp), intent(in) :: d, dt	!Interatomic distance and time step size
real(dp), intent(inout) :: r(3), v(3), a(3)
real(dp) :: ak(3)
real(dp) :: material_height
integer(i8) :: im, jm, km
integer(i8) :: i, j, k
	!Velocity Verlet step
	!It is not necessary to compute the time t = t0 + i*dt
	r = r + v*dt + 0.5*a*dt*dt
	v = v + 0.5*a*dt
	!Computing acceleration after time step
	a = 0
	!Acceleration due to embedded electrons,
	!ONLY computed if there is AT LEAST one embedded electron
	if (N_e .gt. 0) then
		do k = 1, N_e
			call akEE(r, emb(k,:), ak)
			a = a + ak
		end do
	end if
	!Acceleration due to neutral atoms
	!ONLY computed if near the material zone, i.e. 
	material_height = atoms(0,0,0,2) + 0.5*d
	if (r(2) .lt. material_height) then
		call sort_in_material(r,d,Nx,Ny,Nz,im,jm,km)
		do i = im-1, im+1 
			do j = jm-1, jm+1
				do k = km-1, km+1
					!Only compute acceleration with positions with atoms
					if (Z(i, j, k) .ne. 0) then
						call akEA(r, atoms(i, j, k, :), Z(i, j, k), ak)
						a = a + ak
					end if
				end do
			end do
		end do
	end if
	v = v + 0.5*a*dt
end subroutine vv_step_NF

!Velocity Verlet Step for Far Field (bulk electrons)
!Subroutine which computes a time step of the Velocity Verlet algorithm used to
!compute the trajectory of a projectile electron interacting with Nb equivalent
!embedded electrons characterized by the b_neq, b_ec, and b_er arrays. 
subroutine vv_step_FF(Nb,Nxb,Nyb,Nzb,max_eq,max_neq,b_neq,b_ec,b_er,Nx,Ny,Nz,atoms,Z,d,dt,r,v,a)
integer(i8), intent(in) :: Nb !Number of embedded electrons
integer(i8), intent(in) :: Nxb, Nyb, Nzb !Bins boundary indices
integer(i8), intent(in) :: max_eq		!Maximum equivalent charge, 180 from tests (max 182)
integer(i8), intent(in) :: max_neq	!Maximum number of equivalent electrons per bin, N_beam/max_eq+1
integer(i8), intent(in) :: b_neq(-Nxb:Nxb,-Nyb:-1,-Nzb:Nzb)!Current equivalent electron array
integer(i8), intent(in) :: b_ec(-Nxb:Nxb,-Nyb:-1,-Nzb:Nzb,0:max_neq)!Equivalent charges array
real(dp), intent(in) :: b_er(-Nxb:Nxb,-Nyb:-1,-Nzb:Nzb,0:max_neq,3)!Equivalent positions array
integer(i8), intent(in) :: Nx, Ny, Nz
real(dp), intent(in):: atoms(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1,3)
integer, intent(in) :: Z(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1)
real(dp), intent(in) :: d, dt
real(dp), intent(inout) :: r(3), v(3), a(3)
real(dp) :: ak(3)
integer(i8) :: i, j, k
integer(i8) :: im, jm, km
integer(i8) :: neq, ec
real(dp) :: er(3)
real(dp) :: material_height
	!Velocity Verlet step
	!It is not necessary to compute the time t = t0 + i*dt
	r = r + v*dt + 0.5*a*dt*dt
	v = v + 0.5*a*dt
	!Computing acceleration after time step
	a = 0
	!Acceleration due to embedded electrons,
	!ONLY computed if there is AT LEAST one embedded electron
	if (Nb .gt. 0) then
		do i = -Nxb, Nxb
			if (i .ne. 0) then	!To avoid checking for unused bins with index 0 
				do j = -1,-Nyb, -1
					do k = -Nzb, Nzb
						if (b_neq(i,j,k) .gt. 0) then
							do neq = 1, b_neq(i,j,k)
								ec = b_ec(i,j,k,neq)
								er = b_er(i,j,k,neq,:)
								call akEE_FF(ec, r, er, ak)
								a = a + ak
							end do
						end if
					end do
				end do
			end if
		end do
	end if
	!Acceleration due to neutral atoms
	!ONLY computed if near the material zone, i.e. 
	material_height = atoms(0,0,0,2) + 0.5*d
	if (r(2) .lt. material_height) then
	call sort_in_material(r,d,Nx,Ny,Nz,im,jm,km)
		do i = im-1, im+1 
			do j = jm-1, jm+1
				do k = km-1, km+1
					!Only compute acceleration from positions with atoms
					if (Z(i, j, k) .ne. 0) then
						call akEA(r, atoms(i, j, k, :), Z(i, j, k), ak)
						a = a + ak
					end if
				end do
			end do
		end do
	end if
	v = v + 0.5*a*dt
end subroutine vv_step_FF

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

end module