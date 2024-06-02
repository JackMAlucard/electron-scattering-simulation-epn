module m2_dielectric_material_model
	use m0_utilities, &
		only: dp, i8, INTERATOMIC_DIST_SIO2
	implicit none
	contains
	!=======================================================================
	! Subroutine: setup_simple_silica_model
	! Purpose   : Setup a simple silica model by generating atom positions 
	!             and charges within the specified material boundaries.
	! Arguments :
	!   - integer(i8), intent(in) :: material_boundaries(3)
	!       Array containing the material grid boundaries in three dimensions.
	!   - real(dp), allocatable, intent(out) :: atom_positions(:,:,:,:)
	!       Array to store the positions of the atoms.
	!   - integer, allocatable, intent(out) :: atom_charges(:,:,:)
	!       Array to store the charges of the atoms.
	!   - real(dp), allocatable, intent(out) :: atom_charges_cbrt(:,:,:)
	!       Array to store the cubic root of the charges of the atoms.
	!=======================================================================
	! Subroutine that sets up a simple model for SiO2
	! This subroutine takes as input three integers stored in the 
	! "grid_boundaries" array and generates the "grid_points" and "grid_charges"
	! arrays as a result.
	! It uses an approximated interatomic distance d = 1.77 (in Å) taken as the
	! sum of the covalent radius of both silicon (1.11 Å) and oxigen (0.66 Å).
	! The "grid_points" array store the coordinates of each grid point located
	! in the region of space where the material is present in the model.
	! The material extends between the spacial coordinates 
	! x = (-mbi*d, mbi*d),  
	! y = (-mbj, 0), and
	! z = (-mbk*d, mbk*d)
	! ... (PENDING STILL, REVIEW NEXT COMMENT LINES)
	! This subroutine generates an array with the i, j, and k indexes of the atoms
	! in the cubic cell, as well as the x, y, z, position of that atom, and 
	! its charges.
	! The bulk occupies the space in the surface reference system as follows:
	! -Lx/2 < x < Lx/2 (width of the rectangular cuboid)
	!    0  < y < -Ly  (height/depth of the rectangular cuboid)
	! -Lz/2 < z < Lz/2 (length of the rectangular cuboid)
	! The subroutine takes the integer indexes bounds mbi, mbj, mbk corresponding to
	! the dimensions of the rectangular cuboid as well as the array 'atoms' of
	! dimensions (2*mbi+1, mbj+1, 2*mbk+1, 3) as inout variables. The first three 
	! indexes label the atom via its x, y, and z integer position labels, the last 
	! index corresponds to the position coordinates.
	! The subroutine also transforms the magnitudes to au.
	subroutine setup_simple_silica_model &
		(material_boundaries, atom_positions, atom_charges, atom_charges_cbrt)
		implicit none
		integer(i8), intent(in) :: material_boundaries(3)
		real(dp), allocatable, intent(out) :: atom_positions(:,:,:,:)
		integer, allocatable, intent(out) :: atom_charges(:,:,:)
		real(dp), allocatable, intent(out) :: atom_charges_cbrt(:,:,:)
		real(dp) :: x, y, z
		integer(i8) :: mbi, mbj, mbk
		integer(i8) :: i, j, k
		! Getting material grid boundaries
		mbi = material_boundaries(1)
		mbj = material_boundaries(2)
		mbk = material_boundaries(3)
		! Grid points indexes
		! |i| = 2*mbi + 1	--> i = -mbi, ..., 0, ..., mbi
		! |j| = mbj + 1		--> j = 0, ..., -mbj
		! |k| = 2*mbk + 1	--> k = -mbk, ..., 0, ..., mbk
		! The arrays are allocated with an extra position added
		! to both sides of their boundaries, which won't store atoms
		! which is why they are initialized and left with zero values
		allocate(atom_positions(-mbi-1:mbi+1, -mbj-1:1, -mbk-1:mbk+1, 3))
		allocate(atom_charges(-mbi-1:mbi+1, -mbj-1:1, -mbk-1:mbk+1))
		allocate(atom_charges_cbrt(-mbi-1:mbi+1, -mbj-1:1, -mbk-1:mbk+1))
		atom_positions = 0
		atom_charges = 0
		atom_charges_cbrt = 0
		do i = -mbi, mbi
			x = i*INTERATOMIC_DIST_SIO2
			do j = 0, -mbj, -1
				y = j*INTERATOMIC_DIST_SIO2
				do k = -mbk, mbk
					z = k*INTERATOMIC_DIST_SIO2
					atom_positions(i,j,k,:) = (/x, y, z/)
					! Silicon positions, i and k indexes are both even 
					if (mod(abs(i),2) .eq. 0 .and. mod(abs(k),2) .eq. 0) then
						atom_charges(i,j,k) = 14
					! Empty positions, i and k indexes are both odd
					else if (mod(abs(i),2) .eq. 1 .and. mod(abs(k),2) .eq. 1) then
						atom_charges(i,j,k) = 0
					! Oxygen positions, all other combinations
					else
						atom_charges(i,j,k) = 8
					end if
				end do
			end do
		end do
		atom_charges_cbrt = atom_charges**(1/3._dp)
	end subroutine setup_simple_silica_model

end module m2_dielectric_material_model