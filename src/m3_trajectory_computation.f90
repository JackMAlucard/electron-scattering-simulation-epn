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
	!subroutine akEE(rp, rt, a)
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
	!subroutine akEA(rp, rt, Z, a)
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
		!Electron-atom relative position
		rs = rp - rt
		r = norm2(rs)
		!Computing chi and psi
		chi = 0
		psi = 0
		do i = 1, 3
			aux = a_coef(i)*dexp(-b_coef(i)*r*b0_inv*cbrt_Z)
			chi = chi + aux
			aux = aux*b_coef(i)
			psi = psi + aux
		end do
		a = -(Z/r**3)*(chi + psi*b0_inv*cbrt_Z*r)*rs
	!	a = -Z*( (chi/(r**3)) + ((psi*b0_inv*cbrt_Z)/(r**2)) )*rs
	end subroutine acceleration_due_to_atom

!Subroutine which computes a time step of the Velocity Verlet algorithm
!used to compute the trajectory of a projectile electron interacting with N_e
!embedded electrons out of N_et total possible embedded electrons and a **SCC**
!bulk of N_a=(2*Nx + 1)*(Ny + 1)*(2*Nz + 1) neutral atoms with atomic number Z.
!The position of the embedded electrons is stored in the array e_emb(N_et,3).
!Acceleration due to embedded electrons is computed for EVERY embedded electron
!and ONLY if there is AT LEAST one embedded electron.
!Acceleration due to neutral atoms is ONLY computed if near the material zone, 
!i.e. its height is less than half the interatomic distance and lower.
subroutine vv_step(N_e, N_et, emb, Nx, Ny, Nz, atoms, Z, d, dt, r, v, a)
integer(i8), intent(in) :: N_e, N_et, Nx, Ny, Nz
integer, intent(in) :: Z(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1)
real(dp), intent(in):: emb(N_et,3), atoms(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1,3)
real(dp), intent(in) :: d, dt	!Interatomic distance and time step size
real(dp), intent(inout) :: r(3), v(3), a(3)
real(dp) :: ak(3)
real(dp) :: material_height
!integer(i8) :: im, jm, km
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
		do i = -Nx, Nx 
			do j = -Ny, 0, -1
				do k = -Ny, Ny
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
end subroutine vv_step

!Subroutine which computes the trajectory of the i-th projectile electron, out
!of N electrons in the electron beam, and plots P points in unit file number ou.
!This subroutine uses the Near Field Regimen, i.e. computes the field of each
!invidivual electron regardless of distance via the vv_step subroutine.
!The projectile interacts with a set of N_e embedded electrons out of a maximum
!of N total possible embedded electrons and a **SCC** bulk of 
!N_a=(2*Nx + 1)*(Ny + 1)*(2*Nz + 1) neutral atoms with atomic number Z.
!The beam initial positions, velocities, and accelerations are stored in the 
!(N,3) arrarys r, v and a respectively.
!The embedded and scattered electrons final positions are stored in the (N,3)
!arrays e_emb and e_sct respectively. 
!The simulation is ended either when the electron distance d to the origin of 
!coordinates (which approximates the center of the bulk) is greater than the 
!initial distance d0, or when the number of iterations k is equal or greater 
!than the estimated number of iterations T (how is this estimated?)
!The simulation is also ended if the electron gets embedded on the bulk, at
!which point the number of embedded electrons N_e is updated/increased.
!I NEED to explain WHY these conditions for the end of the trajectory are set.
subroutine trajectory(i, N, P, ou, Ne, Ns, emb, sct, Nx, Ny, Nz, atoms, Z, d0, T, d, dt, l, r, v, a)
integer(i8), intent(in) :: i, N, P, ou, Nx, Ny, Nz, T
integer, intent(in) :: Z(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1)
real(dp), intent(in) :: atoms(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1,3)
real(dp), intent(in) :: d0, d, dt, l
integer(i8), intent(inout) :: Ne, Ns
real(dp), intent(inout) :: emb(N,3), sct(N,3)
real(dp), intent(inout) :: r(N,3), v(N,3), a(N,3)
real(dp) :: dk
integer(i8) :: j, k
logical :: end_embedded, end_scattered, end_iterations	!Trajectory end conditions
logical :: in_material
real(dp) :: r_mat_i(3), r_mat_f(3), dr_mat(3), d_mat, d_embedded
real(dp) :: tk
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
			write(ou,*) tk, r(i,:)
			j = j + 1
		end if
		!Computing next Velocity Verlet step
		call vv_step(Ne, N, emb, Nx, Ny, Nz, atoms, Z, d, dt, r(i,:), v(i,:), a(i,:))
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
				print*, "Final electron position:", r(i,:)
				print*, "Attenuation distance (random exp):", d_embedded
				print*, "Distance travelled inside material:", d_mat
				print*, "x_scc: +/-", atoms(Nx,0,0,1)
				print*, "y_scc:", atoms(0,-Ny,0,2)
				print*, "z_scc: +/-", atoms(0,0,Nz,3)
				print*, "k", k
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
end subroutine trajectory

end module m3_trajectory_computation