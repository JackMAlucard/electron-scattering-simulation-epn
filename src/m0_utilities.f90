module m0_utilities
	implicit none

	! NUMERICAL STORAGE SIZE PARAMETERS FOR REAL AND INTEGER VALUES***************
	! Double precision real numbers,
	! 15 digits, range 10**(-307) to 10**(307)-1; 64 bits
	integer, parameter :: dp = selected_real_kind(15, 307)
	! Long length for integers, range -2**63 to 2**63-1; 64 bits
	integer, parameter :: i8 = selected_int_kind(18)

	! NUMERICAL CONSTANTS*********************************************************
	real(dp), parameter :: PI = dacos(-1._dp)

	contains
	! UNIT CONVERSION SUBROUTINES*************************************************
	! Energy conversion from keV to hartree Eh (atomic)
	!	real(dp) function EinAU(EinkeV)
	subroutine kev_to_atomic_energy_conversion(energy)
		implicit none
		real(dp), intent(inout) :: energy
		energy = energy*1.d3/27.21139
	end subroutine kev_to_atomic_energy_conversion
	
	! Distance conversion from meters m (SI) to Bohr radii a0 (atomic)
	! real(dp) function dinAUfromSI(dinSI)
	subroutine si_to_atomic_distance_conversion(distance)
		implicit none
		real(dp), intent(inout) :: distance
		distance = distance/0.5291772d-10
	end subroutine si_to_atomic_distance_conversion

	! Distance conversion from angstroms to Bohr radii a0 (atomic)
!	real(dp) function dinAU(dinangstroms)
	subroutine angstrom_to_atomic_distance_conversion(distance)
		implicit none
		real(dp), intent(inout) :: distance
		distance = distance/0.5291772
	end subroutine angstrom_to_atomic_distance_conversion

	! Distance conversion from Bohr radii a0 (atomic) to meters (SI)
!	real(dp) function dinSI(dinAU)
	subroutine atomic_to_si_distance_conversion(distance)
		implicit none
		real(dp), intent(inout) :: distance
		distance = distance*0.5291772d-10
	end subroutine atomic_to_si_distance_conversion

	! RANDOM NUMBER GENERATION SUBROUTINES****************************************
	! Subroutine which generates a random number following the standard uniform
	! distribution over the range (0,1] to avoid possible error if 0 is generated.
	! Adapted from: https://masuday.github.io/fortran_tutorial/random.html
	subroutine random_std_uniform(u)
		implicit none
		real(dp), intent(out) :: u
		real(dp) :: r
		call random_number(r)
		u = 1 - r
	end subroutine random_std_uniform

	! Subroutine which generates two random numbers following the standard 
	! normal distribution (mean = 0 and standard deviation = 1).
	! Adapted from: https://masuday.github.io/fortran_tutorial/random.html
	subroutine random_std_normal(z1, z2)
		implicit none
		real(dp), intent(out) :: z1, z2
		real(dp) :: u1, u2
		call random_std_uniform(u1)
		call random_std_uniform(u2)
		z1 = dsqrt(-2*dlog(u1))*dcos(2*PI*u2)
		z2 = dsqrt(-2*dlog(u1))*dsin(2*PI*u2)
	end subroutine random_std_normal

	! Subroutine which generates two random number following the general
	! normal distribution with mean 'mu' and standard deviation 'sigma'.
	subroutine random_normal(mu, sigma, x1, x2)
		implicit none
		real(dp), intent(in) :: mu, sigma
		real(dp), intent(out) :: x1, x2
		real(dp) :: z1, z2
		call random_stdnormal(z1, z2)
		x1 = mu + sigma*z1
		x2 = mu + sigma*z2
	end subroutine random_normal

	! Subroutine which generates a random number following the exponential
	! probability distribution with parameter 'lambda' > 0.
	! Adapted from:
	! https://www.eg.bucknell.edu/~xmeng/Course/CS6337/Note/master/node50.html
	subroutine random_exponential(lambda, x)
		implicit none
		real(dp), intent(in) :: lambda
		real(dp), intent(out) :: x
		real(dp) :: r
		call random_std_uniform(u)
		x = -dlog(u)/lambda
	end subroutine random_exponential

	!VECTOR TRANSFORMATIONS SUBROUTINES*******************************************
	subroutine vector_translation(traslation_vector, vector)
		implicit none
		real(dp), intent(in) :: translation_vector(3)
		real(dp), intent(inout) :: vector
		vector = vector + translation_vector
	end subroutine vector_translation
	
!	subroutine matrix_vector_product(M, V)
	subroutine rotation_about_x_axis(theta, vector)
		implicit none
		real(dp), intent(in) :: angle
		real(dp), intent(inout) :: vector(3)
		real(dp) :: R(3,3), aux_vector(3)
		integer :: i
		
		!Rotation matrix
		R(1,:) = (/1._dp, 0._dp, 0._dp/)
		R(2,:) = (/0._dp, dcos(theta), -dsin(theta)/)
		R(3,:) = (/0._dp, dsin(theta), dcos(theta)/)
		
		aux_vector = 0

		do i = 1, 3
			aux_vector(:) = aux_vector(:) + Rx(:,i)*vector(i)
		end do
		
		vector = aux_vector

	end subroutine rotation_about_x_axis

end module m0_utilities