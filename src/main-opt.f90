program main
	use m0_utilities, &
		only: i8, dp, read_input_parameters, open_output_files, &
		write_simulation_status_information
	use m1_electron_beam_model, &
		only: electron_beam_parameters_unit_conversion, setup_electron_beam_model
	use m2_dielectric_material_model, &
		only: setup_simple_silica_model
	use m3_trajectory_computation, &
		only: number_of_iterations_estimation, compute_scattering_angles
	use m4_optimized_trajectory_computation, &
		only: setup_cells_and_super_electrons, compute_trajectory_optimized, &
		write_final_super_electron_distribution
	implicit none
	! m1_electron_beam_model module variables ************************************
	integer(i8) :: num_electrons
	real(dp) :: spot_size_factor
	real(dp) :: beam_energy, energy_spread
	real(dp) :: beam_target_distance, grazing_angle
	logical :: beam_model_output_saving_enabled
	logical :: material_model_output_saving_enabled
	real(dp), allocatable :: electron_positions(:,:), electron_velocities(:,:)
	real(dp), allocatable :: electron_accelerations(:,:)
	! m2_dielectric_material_model module variables ******************************
	integer(i8) :: material_boundaries(3)
	real(dp), allocatable :: atom_positions(:,:,:,:)
	integer, allocatable :: atom_charges(:,:,:)
	real(dp), allocatable :: atom_charges_cbrt(:,:,:)
	! m3_trajectory_computation module variables *********************************
	integer(i8) :: max_iterations
	real(dp) :: alpha, beta
	! m4_optimized_trajectory_computation variables ******************************
	integer(i8) :: num_plot_ploints, partition_boundaries(3)
	real(dp) :: dt, r(3), v(3), a(3)
	integer(i8) :: num_embedded, num_scattered
	integer(i8), allocatable :: num_super_electrons(:,:,:)
	integer, allocatable :: super_electron_charges(:,:,:,:)
	real(dp), allocatable :: super_electron_positions(:,:,:,:,:)
	real(dp), allocatable :: embedded_positions(:,:), scattered_positions(:,:)
	! Additional local variables *************************************************
	integer(i8) :: i, j, k
	real(dp) :: start_time, end_time, total_time
	real(dp) :: iteration_start_time, iteration_end_time, iteration_time
	integer(i8) :: previous_num_embedded, previous_num_scattered
	integer, parameter :: output_unit = 7
	integer :: seed_size
	integer, allocatable :: seed(:)
	
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
		beam_model_output_saving_enabled, material_model_output_saving_enabled, &
		num_plot_ploints, dt)
	
	! SETTING UP ELECTRON BEAM MODEL**********************************************
	! Converting input parameters to the appropriate units
	call electron_beam_parameters_unit_conversion &
		(beam_energy, beam_target_distance, grazing_angle)
	! Initializing electron beam model variables and arrays
	call setup_electron_beam_model &
		(num_electrons, spot_size_factor, beam_energy, energy_spread, &
		beam_target_distance, grazing_angle, beam_model_output_saving_enabled, &
		electron_positions, electron_velocities, electron_accelerations)
	
	! SETTING UP DIELECTRIC MATERIAL MODEL ***************************************
	call setup_simple_silica_model &
		(material_boundaries, material_model_output_saving_enabled, &
		atom_positions, atom_charges, atom_charges_cbrt)
	
	! COMPUTING OPTIMIZED ELECTRON TRAJECTORIES **********************************
	! Setting up partition cells and super electron arrays
	call setup_cells_and_super_electrons &
		(num_electrons, material_boundaries, partition_boundaries, &
		num_super_electrons, super_electron_charges, super_electron_positions)
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
	!! Getting total simulation start time
	call cpu_time(start_time)
	! Computing the trajectories of each electron in the beam
	do i = 1, num_electrons
		print*, "Simulating", i, "out of", num_electrons, &
		"electron trajectories (optimized)"
		! Computing current iteration simulation time
		call cpu_time(iteration_start_time)
		! Getting current electron initial conditions
		r = electron_positions(i,:)
		v = electron_velocities(i,:)
		a = electron_accelerations(i,:)
		! Loose estimation of maximum number of iterations
		call number_of_iterations_estimation &
			(r, v, dt, max_iterations, num_plot_ploints)
		! Computing current electron trajectory (optimized)
		call compute_trajectory_optimized &
			(num_plot_ploints, max_iterations, output_unit+1, material_boundaries, & 
			atom_positions, atom_charges, atom_charges_cbrt, partition_boundaries, & 
			num_super_electrons, super_electron_positions, super_electron_charges, & 
			dt, r, v, a, num_embedded, num_scattered, embedded_positions, & 
			scattered_positions)
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
		! Computing current iteration time and simulation time thus far
		call cpu_time(iteration_end_time)
		iteration_time = iteration_end_time - iteration_start_time
		total_time = iteration_end_time - start_time
		! Writing iteration and total simulation times to file
		write(output_unit+7, *) i, iteration_time, total_time
		! Printing simulation information on console
		print*, "  Number of electrons embedded: ", num_embedded
		print*, "  Number of electrons scattered:", num_scattered
		print*, "Iteration time:", iteration_time
		print*, "Simulation time thus far:", total_time
		! Writing simulation status information to output file
		call write_simulation_status_information &
			(output_unit, i, num_electrons, num_embedded, num_scattered, total_time)
	end do
	! Closing files that stored the output results of the simulation
	do i = 1, 7
		close(output_unit+i)
	end do
	! Computing total simulation time
	call cpu_time(end_time)
	total_time = end_time - start_time
	print*, "SIMULATION COMPLETE!"
	print*, "  Number of electrons embedded: ", num_embedded
	print*, "  Number of electrons scattered:", num_scattered
	print*, "Total simulation time:", total_time
	! Writing final super electron distribution magnitudes to output file
	call write_final_super_electron_distribution &
		(partition_boundaries, num_super_electrons, super_electron_charges, &
		super_electron_positions)
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
	deallocate(num_super_electrons)
	deallocate(super_electron_positions)
	deallocate(super_electron_charges)

end program