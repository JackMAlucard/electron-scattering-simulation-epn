module m2_dielectric_material_model

	use m0_utilities, &
		only: dp, i8, INTERATOMIC_DIST_SIO2

	implicit none
	contains
	!=============================================================================
	! Subroutine: setup_simple_silica_model
	! Purpose   : Set up a simple silica (SiO2) model by generating atom
	!             positions, atomic numbers, and cubic roots of atomic
	!             numbers on a 3D grid within the specified material grid
	!             index boundaries.
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
	!             - The arrays used to store the positions, atomic numbers and
	!               cubic roots of atomic numbers are allocated with one extra
	!               index in each direction for all dimensions (i, j, and k).
	!               This additional space allows for simpler handling of boundary
	!               conditions when traversing the arrays for calculations.
	! Arguments :
	!   - integer(i8), intent(in) :: material_boundaries(3)
	!       Array containing each material grid index boundary.
	!   - logical, intent(in) :: output_saving_enabled
	!       Flag to determine if the material model atom positions are saved.
	!   - real(dp), allocatable, intent(out) :: atom_positions(:,:,:,:)
	!       Array to store the positions of the atoms with each coordinate in
	!       atomic units (a0).
	!   - integer, allocatable, intent(out) :: atomic_numbers(:,:,:)
	!       Array to store the atomic numbers of the atoms.
	!   - real(dp), allocatable, intent(out) :: atomic_numbers_cbrt(:,:,:)
	!       Array to store the cubic roots of the atomic numbers of the atoms.
	!=============================================================================
	subroutine setup_simple_silica_model &
		(material_boundaries, output_saving_enabled, atom_positions, &
		atomic_numbers, atomic_numbers_cbrt)
		implicit none

		! Input/Output variables
		integer(i8), intent(in) :: material_boundaries(3)
		logical, intent(in) :: output_saving_enabled
		real(dp), allocatable, intent(out) :: atom_positions(:,:,:,:)
		integer, allocatable, intent(out) :: atomic_numbers(:,:,:)
		real(dp), allocatable, intent(out) :: atomic_numbers_cbrt(:,:,:)

		! Local variables
		real(dp) :: x, y, z
		integer(i8) :: mbi, mbj, mbk
		integer(i8) :: i, j, k, n

		! Extract material grid index boundaries
		mbi = material_boundaries(1)
		mbj = material_boundaries(2)
		mbk = material_boundaries(3)

		! Allocate arrays for atom positions, atomic numbers, and cubic root of 
		! atomic numbers with one extra index in each direction for all dimensions;
		! outermost positions will not store atoms
		allocate(atom_positions(-mbi-1:mbi+1, -mbj-1:1, -mbk-1:mbk+1, 3))
		allocate(atomic_numbers(-mbi-1:mbi+1, -mbj-1:1, -mbk-1:mbk+1))
		allocate(atomic_numbers_cbrt(-mbi-1:mbi+1, -mbj-1:1, -mbk-1:mbk+1))

		! Initialize arrays to zero
		atom_positions = 0
		atomic_numbers = 0
		atomic_numbers_cbrt = 0

		! Set atom positions and atomic numbers based on grid position
		do i = -mbi, mbi
			x = i*INTERATOMIC_DIST_SIO2
			do j = 0, -mbj, -1
				y = j*INTERATOMIC_DIST_SIO2
				do k = -mbk, mbk
					z = k*INTERATOMIC_DIST_SIO2
					atom_positions(i,j,k,:) = (/x, y, z/)
					if (mod(abs(i),2) .eq. 0 .and. mod(abs(k),2) .eq. 0) then
						! Assign silicon atomic number where both i and k are even
						atomic_numbers(i,j,k) = 14
					else if (mod(abs(i),2) .eq. 1 .and. mod(abs(k),2) .eq. 1) then
						! Assign zero atomic number to indices where both i and k are odd
						atomic_numbers(i,j,k) = 0
					else
						! Assign oxygen atomic number to all other positions
						atomic_numbers(i,j,k) = 8
					end if
				end do
			end do
		end do

		! Compute the cubic root of the atomic numbers
		atomic_numbers_cbrt = atomic_numbers**(1/3._dp)

		! Save silicon and oxygen atom positions if output saving is enabled
		if (output_saving_enabled) then
			open(unit=42, file='material_model.dat', status='replace', action='write')
			do n = 1, 2
				do i = -mbi, mbi
					do j = 0, -mbj, -1
						do k = -mbk, mbk
							if (n .eq. 1 .and. atomic_numbers(i,j,k) .eq. 14) then
								! Write silicon atom positions in the first iteration
								write(42, *) atom_positions(i,j,k,:)
							else if (n .eq. 2 .and. atomic_numbers(i,j,k) .eq. 8) then
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