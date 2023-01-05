module m1
use m0!, only: sp, dp, i8
implicit none


contains

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