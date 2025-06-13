program scattering_angles_progression
	implicit none
	integer, parameter :: dp = selected_real_kind(15, 307)
	! THIS PROGRAM ONLY WORKS FOR REGULAR SIMULATIONS, NOT FOR OPTIMIZED
	! The program will be run from the results/main folder 
	! It will scan all sea.dat files in the subfolders main-f1, ..., main-f9, main-g1, ..., main-g9.
	! And it will save the results in a corresponding sea_progression.dat file in said folder
	! IMPORTANT: I will need a different program to copy all sea_progression.dat files from each
	! subfolder to another single desired path, adding the suffix -fx or -gx.
	! The progression will be measured as follows:
	! The for loop will go from 1 to 5000, it will only consider scattering angle values for the
	! mean where they are present, that is, simulations where electrons were embedded will not
	! be used in the mean. The range for the mean will be a variable to see which one fits best,
	! but it will start with 50. 
	!
	! SCATTERING ANGLES FILE VALUES
	real(dp) :: elevation_angle_alpha, azimuthal_angle_beta
	integer :: aux_int
	! ARRAYS TO STORE VALUES
	real(dp), allocatable :: alpha_array(:), beta_array(:)
	! MEAN VALUES
	integer :: half_interval, scattered_mean_computation_threshold
	real(dp) :: alpha_mean, beta_abs_mean
	integer :: num_scattered_in_interval
	! FOR TESTING
	integer :: num_scattered_alpha, num_scattered_beta
	! FILE NAMES
	character(2) :: suffix ! e.g. f1, g9
	character(*), parameter :: suffix_letters(2) = (/'f', 'g'/)
	character(*), parameter :: suffix_numbers(9) = (/'1', '2', '3', '4', '5', '6', '7', '8', '9'/)
	character(27) :: output_filename
	character(106) target_folder
	character(10) :: target_filename
	! Number of rows on each file
	integer, parameter :: num_rows(18) = &
		(/3990, 3642, 3408, &
			4097, 3809, 3509, &
			4154, 3867, 3615, &
			3079, 2677, 2180, &
			3203, 2878, 2306, &
			3381, 2968, 2515/)
	! OTHER VARIABLES
	integer :: i, j, k, n
	integer :: copy_results
	
	! For each angle value, I'll have a 1D array of dimensions (1 + 25:5000 + 25)
	! The matrix will be initialized to -1._dp
	! The first column will store the actual angle values, for when electrons are scattered
	! Leaving the positions where there are no scattered electrons as -1 (less than 0)
	! The elevation angle is strictly positive already, 
	! and for the azimuthal angle, the absolute value will be used
	! The mean will be computed for each point considering a variable interval (with a test value of 50),
	! Considering only the number of scattered electrons in the interval
	
	! Read options for half interval length and whether or not to copy results  
	print "('Input half interval length (int, e.g. 250): ', $)"
	read*, half_interval
	print*
	print*
	
	print*, "Copy results to thesis Graphics folder?"
	print*, "1. Yes"
	print*, "2. No"
	print "('Choose: ', $)"
	read*, copy_results
	print*
	
	! Allocate angle value arrays
	allocate(alpha_array(1 - half_interval:5000 + half_interval))
	allocate(beta_array(1 - half_interval:5000 + half_interval))
	
	! Loop over all data files
	do i = 1, 2
		do j = 1, 9
			! Open file
			suffix = suffix_letters(i)//suffix_numbers(j)
			open(unit=10, file='main-'//suffix//'/sea.dat', status='unknown')
			
			! Initialize angle arrays to -1
			alpha_array = -1
			beta_array = -1
			
			! Load data from file to array
			do k = 1, num_rows(j + (i-1)*9)
				read(10,*) elevation_angle_alpha, azimuthal_angle_beta, aux_int
				
				! Only store values where electrons are scattered
				alpha_array(aux_int) = elevation_angle_alpha
				beta_array(aux_int) = dabs(azimuthal_angle_beta)
				
			end do
			
			close(10)
			
			! Compute mean values and save to file
			output_filename = 'main-'//suffix//'/sea_progression.dat'
			open(unit=10, file=output_filename, status='unknown')
			
			do k = 1, 5000
				! Initialize scattering angle progression mean values
				alpha_mean = 0
				beta_abs_mean = 0
!				num_scattered_in_interval = 0
				! FOR TESTING
				num_scattered_alpha = 0
				num_scattered_beta = 0
				
				do n = k - half_interval, k + half_interval
				 ! SIMPLIFY IF TESTING GOES WELL
					if (alpha_array(n) .gt. 0) then
						alpha_mean = alpha_mean + alpha_array(n)
						num_scattered_alpha = num_scattered_alpha + 1
					end if
					
					if (beta_array(n) .gt. 0) then
						beta_abs_mean = beta_abs_mean + beta_array(n)
						num_scattered_beta = num_scattered_beta + 1
					end if
					
				end do

				! FOR TESTING
				if (num_scattered_alpha .ne. num_scattered_beta) then
					print*, 'ERROR: num_scattered_alpha .ne. num_scattered_beta'
					print*, 'main-'//suffix//'/sea.dat'
					print*, '  k', k
					print*, '  num_scattered_alpha', num_scattered_alpha
					print*, '  num_scattered_beta', num_scattered_beta
				end if
				
				! The means is computed as different than zero only if there are more
				! scattered electrons than a specific value
				scattered_mean_computation_threshold = 10
!				scattered_mean_computation_threshold = half_interval/100
!				scattered_mean_computation_threshold = half_interval/200
				
				if (num_scattered_alpha .ge. scattered_mean_computation_threshold) then
					alpha_mean = alpha_mean/num_scattered_alpha
				else
					! THIS CAN BE IMPROVED SO AS TO TAKE THE PREVIOUS MEAN VALUE, 
					! TO SHOW PROGRESSION, WHICH IS THE IDEA, AND NOT SET IT BACK TO ZERO
					alpha_mean = 0
				end if
				
				if (num_scattered_alpha .ge. scattered_mean_computation_threshold) then
					beta_abs_mean = beta_abs_mean/num_scattered_beta
				else
					! THIS CAN BE IMPROVED SO AS TO TAKE THE PREVIOUS MEAN VALUE, 
					! TO SHOW PROGRESSION, WHICH IS THE IDEA, AND NOT SET IT BACK TO ZERO
					beta_abs_mean = 0
				end if

				write(10,*) alpha_mean, beta_abs_mean, k
				
			end do
			
			close(10)
			
			! Copy files, if option is toggled
			if (copy_results .eq. 1) then
				target_folder = '../../../epn-graduation-project-thesis/&
				&Graphics/scattering_simulation_plots/scattering_angles_progression/'
				target_filename = 'sap-'//suffix//'.dat'
				
!				print*, target_folder
!				print*, target_folder//target_filename
				call system ("cp "//output_filename//" "//target_folder//target_filename)
				
			end if
			
		end do
	end do
	
end program scattering_angles_progression