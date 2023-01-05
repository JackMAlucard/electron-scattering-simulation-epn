module MS2
use kind_parameters, only: sp, dp, i8
use units!, only: dinAUfromSI
use MS1, only: random_uniform
implicit none

contains

!This subroutine generates an array with the i, j, and k indexes of the atoms
!in the cubic cell, as well as the x, y, z, position of that atom.
!it takes as input parameters the interatomic distance d (in Å), and
!the dimensions of the rectangular cuboid of the space where the atoms are
!located Lx, Ly, Lz in (in Å), width, height, and length respectively.
!The bulk occupies the space in the surface reference system as follows:
!-Lx/2 < x < Lx/2 (width of the rectangular cuboid)
!   0  < y < -Ly  (height/depth of the rectangular cuboid)
!-Lz/2 < z < Lz/2 (length of the rectangular cuboid)
!The subroutine outputs the integer indices bounds Nx, Ny, Nz corresponding to
!the dimensions of the rectangular cuboid as well as the array 'atoms' of
!dimensions (2*Nx+1, Ny+1, 2*Nz+1, 3). The first three indices label the atom
!via its x, y, and z integer position labels, the last index labels the atom
!position coordinates
!The subroutine also transforms the magnitudes to au.
!All the input parameters are optional, and in case none in inserted, the
!standard values are:
!d = 2.35 Å
!Lx = 4*2.35 Å
!Ly = 2.35 Å
!Lz = 4*2.35 Å
subroutine simple_cubic_crystal(atoms, Nx, Ny, Nz, N, d_in, Lx_in, Ly_in, Lz_in)
real(dp), allocatable, intent(out):: atoms(:,:,:,:)	!Indexes and positions
integer(i8), intent(out) :: Nx, Ny, Nz, N	!Indexes bounds and number of atoms
real(dp), intent(in), optional :: d_in				!Interatomic distance
real(dp), intent(in), optional :: Lx_in, Ly_in, Lz_in	!Cuboid dimensions
integer(i8) :: i, j, k
real(dp) :: x, y, z
real(dp) :: d, Lx, Ly, Lz
	
	!Standard values
	if (present(d_in)) then
		d = d_in
	else
		d = 2.35
	end if
	if (present(Lx_in)) then
		Lx = Lx_in
	else
		Lx = 3*2.35
	end if
	if (present(Ly_in)) then
		Ly = Ly_in
	else
		Ly = 2.35
	end if
	if (present(Lz_in)) then
		Lz = Lz_in
	else
		Lz = 3*2.35
	end if
	
	!Unit conversion
	d = dinAU(d)
	Lx = dinAU(Lx)
	Ly = dinAU(Ly)
	Lz = dinAU(Lz)

	Nx = dint((Lx/2)/d) + 1
	Ny = dint(Ly/d) + 1
	Nz = dint((Lz/2)/d) + 1
	print*, "Lx/2", Lx/2, "Ly", Ly, "Lz/2", Lz/2
	print*, "Nx", Nx, "Ny", Ny, "Nz", Nz
	print*, "x_Nx", Nx*d, "y_Ny", Ny*d, "z_Nz", Nz*d

	!|i| = 2*Nx + 1	: -Nx, ..., 0, ..., Nx
	!|j| = Ny + 1		: 0, ..., Ny
	!|k| = 2*Nz + 1	: -Nz, ..., 0, ..., Nz
	N = (2*Nx + 1)*(Ny + 1)*(2*Nz + 1)
	print*, "Total atoms:", N
	allocate(atoms(-Nx:Nx, -Ny:0, -Nz:Nz, 3))

	do i = -Nx, Nx
		x = i*d
		do j = 0, -Ny, -1
			y = j*d
			do k = -Nz, Nz
				z = k*d
				atoms(i, j, k, :) = (/x, y, z/)
			end do
		end do
	end do

end subroutine

!Subroutine that simulates a simple model of SiO2
!This subroutine generates an array with the i, j, and k indexes of the atoms
!in the cubic cell, as well as the x, y, z, position of that atom, and 
!its charges.
!It takes as input parameters the interatomic distance d=1.77 (in Å), taken as 
!the sum of the covalent radius of both silicon (1.11 A) and oxigen (0.66 A)
!The dimensions of the rectangular cuboid of the space where the atoms are
!located Lx, Ly, Lz in (in Å), width, height, and length respectively.
!NOT USED AS INPUT NOR OUTPUT, JUST AS A CONCEPT
!The bulk occupies the space in the surface reference system as follows:
!-Lx/2 < x < Lx/2 (width of the rectangular cuboid)
!   0  < y < -Ly  (height/depth of the rectangular cuboid)
!-Lz/2 < z < Lz/2 (length of the rectangular cuboid)
!The subroutine takes the integer indices bounds Nx, Ny, Nz corresponding to
!the dimensions of the rectangular cuboid as well as the array 'atoms' of
!dimensions (2*Nx+1, Ny+1, 2*Nz+1, 3) as inout variables. The first three 
!indices label the atom via its x, y, and z integer position labels, the last 
!index corresponds to the position coordinates.
!The subroutine also transforms the magnitudes to au.
subroutine simple_silica_model(ssm_r, ssm_Z, Nx, Ny, Nz, d_in)
real(dp), allocatable, intent(out):: ssm_r(:,:,:,:)	!Indexes and positions array
integer, allocatable, intent(out):: ssm_Z(:,:,:)	!Indexes and charge array
integer(i8), intent(inout) :: Nx, Ny, Nz	!Indexes bounds and number of atoms
real(dp), intent(in) :: d_in	!Interatomic distance, sum of two covalent radius
integer(i8) :: i, j, k, N
real(dp) :: x, y, z
real(dp) :: d
	
	!Unit conversion
	d = dinAU(d_in)

	print*, "Nx", Nx, "Ny", Ny, "Nz", Nz
	print*, "Lx", Nx*d, "Ly", Ny*d, "Lz", Nz*d

	!|i| = 2*Nx + 1	: -Nx, ..., 0, ..., Nx
	!|j| = Ny + 1		: 0, ..., Ny
	!|k| = 2*Nz + 1	: -Nz, ..., 0, ..., Nz
	N = (2*Nx + 1)*(Ny + 1)*(2*Nz + 1)
	print*, "Total atom spaces:", N
	allocate(ssm_r(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1, 3))
	allocate(ssm_Z(-Nx-1:Nx+1, -Ny-1:1, -Nz-1:Nz+1))
	ssm_r = 0
	ssm_Z = 0

	do i = -Nx, Nx
		x = i*d
		do j = 0, -Ny, -1
			y = j*d
			do k = -Nz, Nz
				z = k*d
				ssm_r(i, j, k, :) = (/x, y, z/)
				!Silicon positions, i and z indices are both even 
				if (mod(abs(i),2) .eq. 0 .and. mod(abs(k),2) .eq. 0) then
					ssm_Z(i,j,k) = 14
				!Empty positions, i and z indices are both odd
				else if (mod(abs(i),2) .eq. 1 .and. mod(abs(k),2) .eq. 1) then
					ssm_Z(i,j,k) = 0
				!Oxygen positions, all other combinations
				else
					ssm_Z(i,j,k) = 8
				end if
			end do
		end do
	end do

end subroutine simple_silica_model

!This subroutine generates an array with the positions of N electrons randomly
!distributed uniformly in a cuboid of width Lx, height Ly, and length Lz (in Å).
!The electrons occupies the space in the surface reference system as follows:
!-Lx/2 < x < Lx/2 (width of the rectangular cuboid)
!   0  < y < -Ly  (height/depth of the rectangular cuboid)
!-Lz/2 < z < Lz/2 (length of the rectangular cuboid)
subroutine random_embedded_electrons(rnd_elect, N, d, Lx, Ly, Lz)
real(dp), allocatable, intent(out) :: rnd_elect(:,:)
integer(i8), intent(in) :: N
real(dp), intent(in) :: d, Lx, Ly, Lz
real(dp) :: Lxi, Lyi, Lzi
integer(i8) :: i, Nx, Ny, Nz
	
	Nx = dint((Lx/2)/d) + 1
	Ny = dint(Ly/d) + 1
	Nz = dint((Lz/2)/d) + 1
	
	allocate(rnd_elect(N,3))
	
	do i=1, N
		call random_uniform(Lxi, -Nx*d, Nx*d)
		call random_uniform(Lyi, -Ny*d, 0._dp)
		call random_uniform(Lzi, -Nz*d, Nz*d)
		rnd_elect(i,:) = (/Lxi, Lyi, Lzi/)
	end do

end subroutine random_embedded_electrons

end module