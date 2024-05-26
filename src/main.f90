program main
	use m0_utilities, &
		only: i8, dp, read_input_parameters, open_output_files, &
		write_simulation_status_information
	use m1_electron_beam_model, &
		only: electron_beam_parameters_unit_conversion, setup_electron_beam_model
	use m2_dielectric_material_model, &
		only: setup_simple_silica_model
	use m3_trajectory_computation, &
		only: number_of_iterations_estimation, compute_trajectory, &
		compute_scattering_angles
	implicit none
	! m1_electron_beam_model module variables ************************************
	integer(i8) :: num_electrons
	real(dp) :: spot_size_factor
	real(dp) :: beam_energy, energy_spread
	real(dp) :: beam_target_distance, grazing_angle
	real(dp), allocatable :: electron_positions(:,:), electron_velocities(:,:)
	real(dp), allocatable :: electron_accelerations(:,:)
	! m2_dielectric_material_model module variables ******************************
	integer(i8) :: material_boundaries(3)
	real(dp), allocatable :: atom_positions(:,:,:,:)
	integer, allocatable :: atom_charges(:,:,:)
	real(dp), allocatable :: atom_charges_cbrt(:,:,:)
	! m3_trajectory_computation module variables *********************************
	integer(i8) :: num_plot_ploints, max_iterations
	real(dp) :: dt, r(3), v(3), a(3)
	integer(i8) :: num_embedded, num_scattered
	real(dp), allocatable :: embedded_positions(:,:), scattered_positions(:,:)
	real(dp) :: alpha, beta
	! Additional local variables *************************************************
	integer(i8) :: i, j, k
	real(dp) :: start_time, end_time, total_time
	integer(i8) :: previous_num_embedded, previous_num_scattered
	integer, parameter :: output_unit = 7
	integer :: seed_size
	integer, allocatable :: seed(:)
	! TEST, GETTING TO THE BOTTOM OF IT
!	real(dp) :: rt(3), cbrt_Z
!	integer :: Z
	integer(i8) :: mbi, mbj, mbk, ii, jj, kk
	
	! SEEDING RANDOM NUMBERS GENERATION FOR REPRODUCIBILITY **********************
	call random_seed(size=seed_size)
	allocate(seed(seed_size))
	seed = 23466992! Arbitrary seed for all elements
	call random_seed(put=seed)
	deallocate(seed)
	
	! READING INPUT PARAMETERS FROM FILE *****************************************
	call read_input_parameters &
		(num_electrons, beam_energy, energy_spread, spot_size_factor, &
		beam_target_distance, grazing_angle, material_boundaries, &
		num_plot_ploints, dt)
	
	! SETTING UP ELECTRON BEAM MODEL**********************************************
	! Converting input parameters to the appropriate units
	call electron_beam_parameters_unit_conversion &
		(beam_energy, beam_target_distance, grazing_angle)
	! Initializing electron beam model variables and arrays
	call setup_electron_beam_model &
		(num_electrons, spot_size_factor, beam_energy, energy_spread, &
		beam_target_distance, grazing_angle, electron_positions, &
		electron_velocities, electron_accelerations)
	
	! SETTING UP DIELECTRIC MATERIAL MODEL ***************************************
	call setup_simple_silica_model &
		(material_boundaries, atom_positions, atom_charges, atom_charges_cbrt)
		
	! COMPUTING ELECTRON TRAJECTORIES ********************************************
	! Allocating and initializing embedded and scattered positions arrays
	allocate(embedded_positions(num_electrons,3))
	allocate(scattered_positions(num_electrons,3))
	embedded_positions = 0
	scattered_positions = 0
	! Initializing number of embedded and scattered electrons
	num_embedded = 0
	num_scattered = 0
	previous_num_embedded = 0
	previous_num_scattered = 0
	! Opening files to store the output results of the simulation
	call open_output_files(output_unit)
	! Getting simulation start time
	call cpu_time(start_time)
	! Computing the trajectories of each electron in the beam
	do i = 1, num_electrons
		print*, "Simulating", i, "out of", num_electrons, "electron trajectories"
		! Getting current electron initial conditions
		r = electron_positions(i,:)
		v = electron_velocities(i,:)
		a = electron_accelerations(i,:)
		! Loose estimation of maximum number of iterations
		call number_of_iterations_estimation &
			(r, v, dt, max_iterations, num_plot_ploints)
		! TEST, GETTING TO THE BOTTOM OF IT
	mbi = material_boundaries(1)
	mbj = material_boundaries(2)
	mbk = material_boundaries(3)
	print*, 'TEST PERSONALIZED INDEXES BOUNDARIES!'
		print*, 'CHARGES BOUNDARIES', lbound(atom_charges, 1), ubound(atom_charges, 1), &
		lbound(atom_charges, 2), ubound(atom_charges, 2), &
		lbound(atom_charges, 3), ubound(atom_charges, 3)
		print*, 'CHARGES_CBRT BOUNDARIES', lbound(atom_charges_cbrt, 1), &
		ubound(atom_charges_cbrt, 1), lbound(atom_charges_cbrt, 2), &
		ubound(atom_charges_cbrt, 2), lbound(atom_charges_cbrt, 3), &
		ubound(atom_charges_cbrt, 3)
		print*, 'POSITIONS BOUNDARIES', lbound(atom_positions, 1), &
		ubound(atom_positions, 1), lbound(atom_positions, 2), &
		ubound(atom_positions, 2), lbound(atom_positions, 3), ubound(atom_positions, 3)
!	do ii = -mbi, mbi
!!			print*, "TEST i LOOP", ii
!				do jj = 0, -mbj, -1
!!				print*, " TEST j LOOP", jj
!					do kk = -mbk, mbk
!!					print*, " TEST k LOOP", kk
!						!Only compute acceleration in positions with atoms
!						Z = atom_charges(ii,jj,kk)
!!						print*, "  TEST charge Z", Z, Z .ne. 0
!						cbrt_Z = atom_charges_cbrt(ii,jj,kk)
!!						print*, "TEST charge_cbrt", cbrt_Z
!						rt = atom_positions(ii,jj,kk,:)
!!						print*, "  TEST position rt", rt
!						if (Z .eq. 0) then
!!							print*, ii, jj, kk
!						end if
!						print*, ii, jj, kk, Z, cbrt_Z, rt
!					end do
!				end do
!			end do
!	print*, "TEST SUCCESS!"
	read(*,*)
		! Computing current electron trajectory
		call compute_trajectory &
			(num_plot_ploints, max_iterations, output_unit+1, material_boundaries, & 
			atom_positions, atom_charges, atom_charges_cbrt, dt, r, v, a, &
			num_embedded, num_scattered, embedded_positions, scattered_positions)
		! Writing empty lines for separation between trajectories in output file
		write(output_unit+1, *)
		write(output_unit+1, *)
		! Updating number of embedded and scattered electrons and writing to files
		! Case: Electron is embedded
		if (num_embedded .gt. previous_num_embedded) then
			previous_num_embedded = num_embedded
			write(output_unit+2, *) i, embedded_positions(num_embedded,:)
		! Case: Electron is scattered
		else if (num_scattered .gt. previous_num_scattered) then
			previous_num_scattered = num_scattered
			write(output_unit+3, *) i, scattered_positions(num_scattered,:)
			! Computing and writing scattering angles to file
			call compute_scattering_angles &
				(num_scattered, scattered_positions, alpha, beta)
			write(output_unit+4, *) i, alpha, beta
		end if
		! Writing final electron position to file 
		write(output_unit+5, *) i, r
		! Writing electrons' accumulated final states to file
		write(output_unit+6, *) i, num_embedded, num_scattered, &
			real(num_embedded,dp)/i, real(num_scattered,dp)/i
		! Computing simulation time up to current electron
		call cpu_time(end_time)
		total_time = end_time - start_time
		! Printing simulation information on console
		print*, "  Number of electrons embedded: ", num_embedded
		print*, "  Number of electrons scattered:", num_scattered
		print*, "Simulation time thus far:", total_time
		! Writing simulation status information to output file
		call write_simulation_status_information &
			(output_unit, i, num_electrons, num_embedded, num_scattered, total_time)
	end do
	! Closing files that stored the output results of the simulation
	do i = 1, 6
		close(output_unit+i)
	end do
	! Computing total simulation time
	call cpu_time(end_time)
	total_time = end_time - start_time
	print*, "SIMULATION COMPLETE!"
	print*, "  Number of electrons embedded: ", num_embedded
	print*, "  Number of electrons scattered:", num_scattered
	print*, "Total simulation time:", total_time
	! Writing final simulation status information to output file
	call write_simulation_status_information &
		(output_unit, i, num_electrons, num_embedded, num_scattered, total_time)
	! Deallocating allocatable arrays used
	deallocate(electron_positions)
	deallocate(electron_velocities)
	deallocate(electron_accelerations)
	deallocate(atom_positions)
	deallocate(atom_charges)
	deallocate(atom_charges_cbrt)
	deallocate(embedded_positions)
	deallocate(scattered_positions)

end program main