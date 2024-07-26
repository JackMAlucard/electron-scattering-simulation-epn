module m2_dielectric_material_model

	use m0_utilities, &
		only: dp, i8, INTERATOMIC_DIST_SIO2

	implicit none
	contains
	!=============================================================================
	! Subroutine: setup_simple_silica_model
	! Purpose   : Set up a simple silica (SiO2) model by generating atom
	!             positions, charges, and cubic root of charges on a 3D grid
	!             within the specified material grid index boundaries.
	!             - The model represents silicon atoms as adjacent to four
	!               oxygen atoms within each plane of the material, with these
	!               planes stacked vertically in the y-axis.
	!             - Material grid points are defined by the indices i, j, and k.
	!               The range of each index is determined by the values in the
	!               material_boundaries = (mbi, mbj, mbk) array as follows:
	!               - i ranges from -mbi to mbi (2*mbi + 1 values),
	!               - j ranges from -mbj to 0 (mbj + 1 values),
	!               - k ranges from -mbk to mbk (2*mbk + 1 values).
	!             - The interatomic distance between different coplanar atoms
	!               and between planes is set as 1.77 Å (converted to atomic
	!               units). This is the sum of the covalent radii of both
	!               silicon (1.11 Å) and oxigen (0.66 Å).
	! Arguments :
	!   - integer(i8), intent(in) :: material_boundaries(3)
	!       Array containing each material grid index boundary.
	!   - logical, intent(in) :: output_saving_enabled
	!       Flag to determine if the material model atom positions are saved.
	!   - real(dp), allocatable, intent(out) :: atom_positions(:,:,:,:)
	!       Array to store the positions of the atoms.
	!   - integer, allocatable, intent(out) :: atom_charges(:,:,:)
	!       Array to store the charges of the atoms.
	!   - real(dp), allocatable, intent(out) :: atom_charges_cbrt(:,:,:)
	!       Array to store the cubic root of the charges of the atoms.
	!=============================================================================
	subroutine setup_simple_silica_model &
		(material_boundaries, output_saving_enabled, atom_positions, &
		atom_charges, atom_charges_cbrt)
		implicit none

		! Input/Output variables
		integer(i8), intent(in) :: material_boundaries(3)
		logical, intent(in) :: output_saving_enabled
		real(dp), allocatable, intent(out) :: atom_positions(:,:,:,:)
		integer, allocatable, intent(out) :: atom_charges(:,:,:)
		real(dp), allocatable, intent(out) :: atom_charges_cbrt(:,:,:)
		
		! Local variables
		real(dp) :: x, y, z
		integer(i8) :: mbi, mbj, mbk
		integer(i8) :: i, j, k, n
		
		! Extract material grid index boundaries
		mbi = material_boundaries(1)
		mbj = material_boundaries(2)
		mbk = material_boundaries(3)
		
		! Allocate arrays for atom positions, atom charges, and cubic root of 
		! atom charges with each index boundary extended by one in each direction;
    ! outermost positions will not store atoms
		allocate(atom_positions(-mbi-1:mbi+1, -mbj-1:1, -mbk-1:mbk+1, 3))
		allocate(atom_charges(-mbi-1:mbi+1, -mbj-1:1, -mbk-1:mbk+1))
		allocate(atom_charges_cbrt(-mbi-1:mbi+1, -mbj-1:1, -mbk-1:mbk+1))

		! Initialize arrays to zero
		atom_positions = 0
		atom_charges = 0
		atom_charges_cbrt = 0

		! Set atom positions and charges based on grid position
		do i = -mbi, mbi
			x = i*INTERATOMIC_DIST_SIO2
			do j = 0, -mbj, -1
				y = j*INTERATOMIC_DIST_SIO2
				do k = -mbk, mbk
					z = k*INTERATOMIC_DIST_SIO2
					atom_positions(i,j,k,:) = (/x, y, z/)
					if (mod(abs(i),2) .eq. 0 .and. mod(abs(k),2) .eq. 0) then
						! Assign silicon atom charge where both i and k are even
						atom_charges(i,j,k) = 14
					else if (mod(abs(i),2) .eq. 1 .and. mod(abs(k),2) .eq. 1) then
						! Assign zero charge to empty positions where both i and k are odd
						atom_charges(i,j,k) = 0
					else
						! Assign oxygen atom charge to all other positions
						atom_charges(i,j,k) = 8
					end if
				end do
			end do
		end do

		! Compute the cubic root of the atom charges
		atom_charges_cbrt = atom_charges**(1/3._dp)

		! Save silicon and oxygen atom positions if output saving is enabled
		if (output_saving_enabled) then
			open(unit=42, file='material_model.dat', status='replace', action='write')
			do n = 1, 2
				do i = -mbi, mbi
					do j = 0, -mbj, -1
						do k = -mbk, mbk
							if (n .eq. 1 .and. atom_charges(i,j,k) .eq. 14) then
								! Write silicon atom positions in the first iteration
								write(42, *) atom_positions(i,j,k,:)
							else if (n .eq. 2 .and. atom_charges(i,j,k) .eq. 8) then
								! Write oxygen atom position in the second iteration
								write(42, *) atom_positions(i,j,k,:)
							end if
						end do
					end do
				end do
				write(42, *)
				write(42, *)
			end do
			close(42)
		end if

	end subroutine setup_simple_silica_model

end module m2_dielectric_material_model