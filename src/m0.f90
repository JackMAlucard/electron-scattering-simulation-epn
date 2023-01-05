module m0
implicit none

!NUMERICAL STORAGE SIZE PARAMETERS FOR REAL AND INTEGER VALUES******************
!Single precision real numbers, 6 digits, range 10**(-37) to 10**(37)-1; 32 bits
integer, parameter :: sp = selected_real_kind(6, 37)
!Double precision real numbers, 15 digits, range 10**(-307) to 10**(307)-1; 64 bits
integer, parameter :: dp = selected_real_kind(15, 307)
!Quadruple precision real numbers, 33 digits, range 10**(-4931) to 10**(4931)-1; 128 bits
integer, parameter :: qp = selected_real_kind(33, 4931)

!Char length for integers, range -2**7 to 2**7-1; 8 bits
integer, parameter :: i1 = selected_int_kind(2)
!Short length for integers, range -2**15 to 2**15-1; 16 bits
integer, parameter :: i2 = selected_int_kind(4)
!Length of default integers, range -2**31 to 2**31-1; 32 bits
integer, parameter :: i4 = selected_int_kind(9)
!Long length for integers, range -2**63 to 2**63-1; 64 bits
integer, parameter :: i8 = selected_int_kind(18)

real(dp), parameter :: PI = dacos(-1._dp)

contains
!UNIT CONVERSION SUBROUTINES****************************************************
!Energy conversion from keV to Eh (au)
real(dp) function EinAU(EinkeV)
implicit none
real(dp), intent(in) :: EinkeV
	EinAU = EinkeV*1.d3/27.21139
end function EinAU

!Distance conversion from angstroms to Bohr radii a0 (au)
real(dp) function dinAU(dinangstroms)
implicit none
real(dp), intent(in) :: dinangstroms
	dinAU = dinangstroms/0.5291772
end function dinAU

!Distance conversion from Bohr radii a0 (au) to meters (SI)
real(dp) function dinSI(dinAU)
implicit none
real(dp), intent(in) :: dinAU
	dinSI = dinAU*0.5291772d-10
end function dinSI

!Distance conversion from meters (SI) to Bohr radii a0 (au)
real(dp) function dinAUfromSI(dinSI)
implicit none
real(dp), intent(in) :: dinSI
	dinAUfromSI = dinSI/0.5291772d-10
end function dinAUfromSI

!PRN GENERATION SUBROUTINES*****************************************************
!Subroutine which generates a pseudorandom number from the standard uniform
!distribution over the range (0,1].
!Taken from https://masuday.github.io/fortran_tutorial/random.html
subroutine random_stduniform(u)
real(dp), intent(out) :: u
real(dp) :: r
	call RANDOM_NUMBER(r)
	u = 1 - r
end subroutine random_stduniform

!Subroutine which generates a pseudorandom number from the uniform
!distribution over the range [a,b), assuming a<b.
!Modified from https://masuday.github.io/fortran_tutorial/random.html
subroutine random_uniform(x,a,b)
real(dp), intent(out) :: x
real(dp), intent(in) :: a, b
real(dp) :: r
	call RANDOM_NUMBER(r)
	x = (b-a)*r + a
end subroutine random_uniform

!Subroutine which generates two pseudorandom number from the normal
!distribution (mean=0 and standard deviation=1).
!Modified from https://masuday.github.io/fortran_tutorial/random.html
subroutine random_stdnormal(x1,x2)
real(dp), intent(out) :: x1, x2
real(dp) :: u1, u2
	call random_stduniform(u1)
	call random_stduniform(u2)
	x1 = dsqrt(-2*dlog(u1))*dcos(2*PI*u2)
	x2 = dsqrt(-2*dlog(u1))*dsin(2*PI*u2)
end subroutine random_stdnormal

!Subroutine which generates two pseudorandom number from the general
!normal distribution (mean=mu and standard deviation=sigma).
subroutine random_normal(x1,x2,mu,sigma)
real(dp), intent(out) :: x1, x2
real(dp), intent(in) :: mu, sigma
real(dp) :: u1, u2
	call random_stdnormal(u1,u2)
	x1 = u1*sigma + mu
	x2 = u2*sigma + mu
end subroutine random_normal

!Subroutine which generates a pseudorandom number following the exponential
!probability distribution with rate parameter lambda > 0 and mean 1/lambda.
!Taken from https://www.eg.bucknell.edu/~xmeng/Course/CS6337/Note/master/node50.html
subroutine random_exponential(x,lambda)
real(dp), intent(out) :: x
real(dp), intent(in) :: lambda
real(dp) :: r
	call random_stduniform(r)
	x = -dlog(r)/lambda
end subroutine random_exponential

!MATRIX PRODUCT SUBRUTINE*******************************************************
subroutine matrix_vector_product(M, V)
real(dp), intent(in) :: M(:,:)
real(dp), intent(inout) :: V(:)
real(dp), allocatable :: V_aux(:)
integer :: N, i, j
	N = size(V)
	allocate(V_aux(N))
	V_aux = 0._dp

	do i=1, N
		do j=1, N
			V_aux(i) = V_aux(i) + M(i,j)*V(j)
		end do
	end do

	V = V_aux
	deallocate(V_aux)
end subroutine matrix_vector_product

end module m0