module MS1
use kind_parameters, only: sp, dp, i8
use units!, only: EinAU, dinAUfromSI
implicit none
real(dp), parameter :: PI = dacos(-1._dp)

!subroutine Main_Subroutine_1()

!end Main_Subroutine_1

contains

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

!Subroutine which generates an array with N electron random starting positions
!inside a circle of radius R moving in the +z' direction with energy E
!The z'-velocities are distributed around the corresponding value of primary
!energy E. The initial accelerations are null.
!The subroutine parameters are the number of electrons N_in, the beam primary
!energy E_in (in keV), the beam radius R_in (in Å), and one out of four types
!of probability distributions for the electron positions:
!1: Normal distribution, R=FWTM*sigma
!2: Normal distribution, R=FWHM*sigma
!3: Non-uniform distribution
!4: Uniform distribution
!The subroutine also transforms the magnitudes to au.
!All the input parameters are optional, and in case none in inserted, the
!standard values are:
!N = 10000 electrons
!E = 1 keV
!R = 400 Å = 4e-5 mm
!prob_dist = 1, normal distribution, R=FWTM*sigma
!The vectors generated are in the beam reference system, the 'primed' system.
subroutine electron_beam(e_pos, e_vel, e_acc, N_in, E_in, R_in, prob_dist_in)
real(dp), allocatable, intent(out) :: e_pos(:,:), e_vel(:,:), e_acc(:,:)
integer(i8), intent(in), optional :: N_in
real(dp), intent(in), optional :: E_in, R_in
integer, intent(in), optional :: prob_dist_in
integer(i8) :: N
real(dp) :: E, R
integer :: prob_dist
real(dp) :: ri, phii, xi, yi, zi, vxi, vyi, vzi!, axi, ayi, azi
real(dp) :: Ei1, Ei2
integer(i8) :: i

	!Standard values
	if (present(N_in)) then
		N = N_in
	else
		N = 10000
	end if
	if (present(E_in)) then
		E = E_in
	else
		E = 1!keV
	end if
	if (present(R_in)) then
		R = R_in
	else
		R = 400!Å = 4e-5 mm
	end if
	if (present(prob_dist_in)) then
		prob_dist = prob_dist_in
	else
		prob_dist = 1!Normal distribution, R=FWTM*sigma
	end if

	!Unit conversion
	E = EinAU(E)
	R = dinAU(R)

	allocate(e_pos(N,3))
	allocate(e_vel(N,3))
	allocate(e_acc(N,3))

	do i = 1, N
		!Generate positions
		!x' and y' coordinates
		select case (prob_dist)
			case (1)
			!Normal distribution, R=FWTM*sigma
			!Full Width at Tenth of Maximum (FWTM) = 2*sqrt(2*ln(10))*sigma ~ 4.29193*sigma
				call random_normal(xi,zi,0._dp,R/4.29193_dp)
			case (2)
			!Normal distribution, R=FWHM*sigma
			!Full Width at Half Maximum (FWHM) = 2*sqrt(2*ln(2))*sigma ~ 2.35482*sigma
				call random_normal(xi,zi,0._dp,R/2.35482_dp)
			!Non-uniform distribution
			case (3)
				call random_uniform(ri, 0._dp, R)
				call random_uniform(phii, 0._dp, 2*PI)
				xi = ri*dcos(phii)
				zi = ri*dsin(phii)
			!Uniform distribution
			case (4)
				call random_stduniform(ri)
				ri = dsqrt(ri)*R
				call random_uniform(phii, 0._dp, 2*PI)
				xi = ri*dcos(phii)
				zi = ri*dsin(phii)
		end select
		!z' coordinate
		yi = 0
		e_pos(i,:) = (/xi, yi, zi/)

		!Generate velocities
		!For 10kev, the speed is ~ 0.19783585 c = 59350755.095 m/s
		!Or 58455098 m/s = 0.194850327 c (with Lorentz tranform?)
		vxi = 0
		vzi = 0

		!Since two numbers are generated in each call of random_normal
		!Using this if-else structure I avoid generating more numbers than necessary
		if (mod(1+i,2) .eq. 0) then
			!Energy values wuth a mean of E, and a mean deviation of 0.1*E
			call random_normal(Ei1,Ei2,E,0.1*E)
			vyi = -dsqrt(2*Ei1)
		else
			vyi = -dsqrt(2*Ei2)
		end if
		e_vel(i,:) = (/vxi, vyi, vzi/)

		!Fill accelerations
		e_acc(i,:) = (/0._dp, 0._dp, 0._dp/)
	end do

end subroutine electron_beam


!Subroutine which transforms the positions and velocities to the surface
!reference frame according to the incident angle theta (in degrees) and the
!distance d (in Å) between the initial electron beam position and the target.
!The subroutine also transforms the magnitudes to au.
!All the input parameters are optional, and in case none is inserted, the
!standard values are:
!theta = 45º
!d = 1.d4 Å = 0.001 mm
subroutine trnsf_to_surf_ref_sys(e_pos, e_vel, e_acc, theta_in, d_in)
real(dp), intent(inout) :: e_pos(:,:), e_vel(:,:), e_acc(:,:)
real(dp), intent(in), optional :: theta_in, d_in
real(dp) :: theta, d
real(dp) :: T(3), R(3,3)
integer(i8) :: N, i

	!Standard values
	if (present(theta_in)) then
		theta = theta_in
	else
		theta = 45!º
	end if
	if (present(d_in)) then
		d = d_in
	else
		d = 1.d4!Å = 0.001 mm
	end if

	!Angle in radians
	theta = theta*PI/180
	print*, "theta (rad)", theta

	d = dinAU(d)
	print*, "d (a0)", d

	!Translation vector
	T = (/0._dp, d*dcos(theta), -d*dsin(theta)/)
	print*, T

	!Rotation matrix
	R(1,:) = (/1._dp, 0._dp, 0._dp/)
	R(2,:) = (/0._dp, dcos(theta), dsin(theta)/)
	R(3,:) = (/0._dp, -dsin(theta), dcos(theta)/)

	!Number of electrons/vectors in array
	N = size(e_pos)/3

	do i=1, N
		!Rotation around x-axis (surface reference system)
		call matrix_vector_product(R, e_pos(i,:))
		call matrix_vector_product(R, e_vel(i,:))

		!Position translation in zy-plane (surface reference system)
		e_pos(i,:) = e_pos(i,:) + T
	end do
end subroutine


end module