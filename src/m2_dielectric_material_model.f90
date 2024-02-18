module m2_dielectric_material_model
	implicit none
	use m0!, only ...
	contains

	! Subroutine that computes integer grid boundaries from rectangular cuboid
	! region of space dimensions
	! subroutine length_to_integer_boundaries(length_boundaries, integer_boundaries)
	! ...

	! Information to be printed from this subroutine
	!	print*, "Nx", Nx, "Ny", Ny, "Nz", Nz
	!	print*, "Lx", Nx*d, "Ly", Ny*d, "Lz", Nz*d
	!	N = (2*Nx + 1)*(Ny + 1)*(2*Nz + 1)
	!	print*, "Total atom spaces:", N

	! Subroutine that sets up a simple model for SiO2
	! This subroutine takes as input three integers stored in the 
	! "grid_boundaries" array and generates the "grid_points" and "grid_charges"
	! arrays as a result.
	! It uses an approximated interatomic distance d = 1.77 (in Å) taken as the
	! sum of the covalent radius of both silicon (1.11 Å) and oxigen (0.66 Å).
	! The "grid_points" array store the coordinates of each grid point located
	! in the region of space where the material is present in the model.
	! The material extends between the spacial coordinates 
	! x = (-Nx*d, Nx*d),  
	! y = (-Ny, 0), and
	! z = (-Nz*d, Nz*d)
	! ... (PENDING STILL, REVIEW NEXT COMMENT LINES)
	! This subroutine generates an array with the i, j, and k indexes of the atoms
	! in the cubic cell, as well as the x, y, z, position of that atom, and 
	! its charges.
	! The bulk occupies the space in the surface reference system as follows:
	! -Lx/2 < x < Lx/2 (width of the rectangular cuboid)
	!    0  < y < -Ly  (height/depth of the rectangular cuboid)
	! -Lz/2 < z < Lz/2 (length of the rectangular cuboid)
	! The subroutine takes the integer indices bounds Nx, Ny, Nz corresponding to
	! the dimensions of the rectangular cuboid as well as the array 'atoms' of
	! dimensions (2*Nx+1, Ny+1, 2*Nz+1, 3) as inout variables. The first three 
	! indices label the atom via its x, y, and z integer position labels, the last 
	! index corresponds to the position coordinates.
	! The subroutine also transforms the magnitudes to au.
	subroutine setup_simple_silica_model(grid_boundaries, grid_points, &
																			 grid_charges)
		integer(i8), intent(in) :: grid_boundaries(3)
		real(dp), allocatable, intent(out):: grid_points(:,:,:,:)
		integer, allocatable, intent(out):: grid_charges(:,:,:)
		real(dp), parameter :: interatomic_distance = 1.77!Å
		real(dp) :: x, y, z
		integer(i8) :: Nx, Ny, Nz
		integer(i8) :: i, j, k
				
		!Unit conversion
		call angstrom_to_atomic_distance_conversion(interatomic_distance)
		
		Nx = grid_boundaries(1)
		Ny = grid_boundaries(2)
		Nz = grid_boundaries(3)
		
		! Grid positions and charges indices
		!|i| = 2*Nx + 1	--> i = -Nx, ..., 0, ..., Nx
		!|j| = Ny + 1		--> j = 0, ..., -Ny
		!|k| = 2*Nz + 1	--> k = -Nz, ..., 0, ..., Nz
		allocate(grid_points(-Nx:Nx,-Ny:0,-Nz:Nz,3))
		allocate(grid_charges(-Nx:Nx,-Ny:0,-Nz:Nz))
		
		! Originally, both arrays had an extra slot after both boundaries
		! I don't remember if this functionality is really necessary,
		! since the for loop does not reach those indices here,
		! however, I did initially fill both arrays with zeroes, 
		! so, I'll keep these lines for the time being here in case they
		! really are necessary after all.
		!	allocate(ssm_r(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1, 3))
		!	allocate(ssm_Z(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1))
		!	ssm_r = 0
		!	ssm_Z = 0

		do i = -Nx, Nx
			x = i*interatomic_distance
			do j = 0, -Ny, -1
				y = j*interatomic_distance
				do k = -Nz, Nz
					z = k*interatomic_distance
					grid_points(i, j, k, :) = (/x, y, z/)
					!Silicon positions, i and k indices are both even 
					if (mod(abs(i),2) .eq. 0 .and. mod(abs(k),2) .eq. 0) then
						grid_charges(i,j,k) = 14
					!Empty positions, i and k indices are both odd
					else if (mod(abs(i),2) .eq. 1 .and. mod(abs(k),2) .eq. 1) then
						grid_charges(i,j,k) = 0
					!Oxygen positions, all other combinations
					else
						grid_charges(i,j,k) = 8
					end if
				end do
			end do
		end do

	end subroutine simple_silica_model

end module m2_dielectric_material_model