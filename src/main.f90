program main

	use m0_utilities, &
		only: i8, dp, read_input_parameters, open_output_files, &
		write_simulation_status_information

	use m1_electron_beam_model, &
		only: electron_beam_parameters_unit_conversion, set_up_electron_beam_model

	use m2_dielectric_material_model, &
		only: set_up_simple_silica_model

	use m3_trajectory_computation, &
		only: number_of_iterations_estimation, compute_trajectory, &
		compute_scattering_angles

	implicit none
	!=============================================================================
	! MAIN PROGRAM VARIABLES
	!=============================================================================

	! Electron beam model variables
	integer(i8) :: num_electrons
	real(dp) :: spot_size_factor
	real(dp) :: beam_energy, energy_spread
	real(dp) :: beam_target_distance, grazing_angle
	real(dp), allocatable :: electron_positions(:,:), electron_velocities(:,:)
	real(dp), allocatable :: electron_accelerations(:,:)
	logical :: beam_model_output_saving_enabled

	! Dielectric material model variables
	integer(i8) :: material_boundaries(3)
	real(dp), allocatable :: atom_positions(:,:,:,:)
	integer, allocatable :: atomic_numbers(:,:,:)
	real(dp), allocatable :: atomic_numbers_cbrt(:,:,:)
	logical :: material_model_output_saving_enabled

	! Electron trajectories' computation variables
	integer(i8) :: num_plot_ploints, max_iterations
	real(dp) :: dt, r(3), v(3), a(3)
	integer(i8) :: num_embedded, num_scattered
	real(dp), allocatable :: embedded_positions(:,:), scattered_positions(:,:)
	real(dp) :: alpha, beta
	integer, parameter :: output_unit = 42
	logical :: electron_trajectories_saving_enabled

	! Additional variables
	integer(i8) :: i, j, k
	real(dp) :: start_time, end_time, total_time
	real(dp) :: iteration_start_time, iteration_end_time, iteration_time
	integer(i8) :: previous_num_embedded, previous_num_scattered
	integer :: seed_size
	integer, allocatable :: seed(:)

	!=============================================================================
	! RANDOM NUMBERS SEEDS AND SIMULATION INPUT PARAMETERS
	!=============================================================================

	! Seed random number generation for reproducibility
	call random_seed(size=seed_size)
	allocate(seed(seed_size))
	seed = 23466992 ! Arbitrary seed for all elements
	call random_seed(put=seed)
	deallocate(seed)
	! Adapted from :
	! https://masuday.github.io/fortran_tutorial/random.html

	! Read simulation input parameters from 'input.txt' file
	call read_input_parameters &
		(num_electrons, beam_energy, energy_spread, spot_size_factor, &
		beam_target_distance, grazing_angle, material_boundaries, &
		beam_model_output_saving_enabled, material_model_output_saving_enabled, &
		electron_trajectories_saving_enabled, num_plot_ploints, dt)

	!=============================================================================
	! SET UP ELECTRON BEAM AND DIELECTRIC MATERIAL MODELS
	!=============================================================================

	! Convert electron beam input parameters to atomic units and radians
	call electron_beam_parameters_unit_conversion &
		(beam_energy, beam_target_distance, grazing_angle)

	! Initialize the positions, velocities, and accelerations for the electrons
	! in the electron beam model
	call set_up_electron_beam_model &
		(num_electrons, spot_size_factor, beam_energy, energy_spread, &
		beam_target_distance, grazing_angle, beam_model_output_saving_enabled, &
		electron_positions, electron_velocities, electron_accelerations)

	! Initialize atom positions, atomic numbers, and cubic root of atomic numbers
	! in the dielectric material model
	call set_up_simple_silica_model &
		(material_boundaries, material_model_output_saving_enabled, &
		atom_positions, atomic_numbers, atomic_numbers_cbrt)

	!=============================================================================
	! COMPUTE ELECTRON TRAJECTORIES
	!=============================================================================

	! Allocate and initialize embedded and scattered positions arrays
	allocate(embedded_positions(num_electrons,3))
	allocate(scattered_positions(num_electrons,3))
	embedded_positions = 0
	scattered_positions = 0

	! Initialize number of embedded and scattered electrons
	num_embedded = 0
	num_scattered = 0
	previous_num_embedded = 0
	previous_num_scattered = 0

	! Open files to save the simulation outputs
	call open_output_files (output_unit, electron_trajectories_saving_enabled)

	! Set total simulation start time
	call cpu_time(start_time)

	! Simulate the trajectory of each electron in the beam
	do i = 1, num_electrons
		print*, "Simulating", i, "out of", num_electrons, "electron trajectories"

		! Set current iteration simulation time
		call cpu_time(iteration_start_time)

		! Get current electron initial conditions
		r = electron_positions(i,:)
		v = electron_velocities(i,:)
		a = electron_accelerations(i,:)

		! Estimate an upper bound to the number of iterations required to simulate
		! a complete electron trajectory
		call number_of_iterations_estimation &
			(r, v, dt, max_iterations, num_plot_ploints)

		! Simulate current electron trajectory
		call compute_trajectory &
			(electron_trajectories_saving_enabled, num_plot_ploints, max_iterations, &
			output_unit+1, material_boundaries, atom_positions, atomic_numbers, &
			atomic_numbers_cbrt, dt, r, v, a, num_embedded, num_scattered, &
			embedded_positions, scattered_positions)

		! Write empty lines for separation between trajectories in output file,
		! if electron trajectories saving is enabled
		if (electron_trajectories_saving_enabled) then
			write(output_unit+1, *)
			write(output_unit+1, *)
		end if

		! Update embedded and scattered electrons' numbers and save to output files
		! Simulation result: Electron is embedded
		if (num_embedded .gt. previous_num_embedded) then
			previous_num_embedded = num_embedded
			write(output_unit+2, *) i, embedded_positions(num_embedded,:)

		! Simulation result: Electron is scattered
		else if (num_scattered .gt. previous_num_scattered) then
			previous_num_scattered = num_scattered
			write(output_unit+3, *) i, scattered_positions(num_scattered,:)

			! Compute scattering angles and save to output file
			call compute_scattering_angles &
				(num_scattered, scattered_positions, alpha, beta)
			write(output_unit+4, *) i, alpha, beta

		end if

		! Save final electron position to output file
		write(output_unit+5, *) i, r

		! Update evolution of the electrons' final state and save to output file
		write(output_unit+6, *) i, num_embedded, num_scattered, &
			real(num_embedded,dp)/i, real(num_scattered,dp)/i

		! Compute current iteration time and simulation time up to current iteration
		call cpu_time(iteration_end_time)
		iteration_time = iteration_end_time - iteration_start_time
		total_time = iteration_end_time - start_time

		! Save current iteration and total simulation times
		write(output_unit+7, *) i, iteration_time, total_time

		! Print simulation information to console
		print*, "  Number of electrons embedded: ", num_embedded
		print*, "  Number of electrons scattered:", num_scattered
		print*, "Iteration time:", iteration_time
		print*, "Simulation time thus far:", total_time

		! Save simulation status information
		call write_simulation_status_information &
			(output_unit, i, num_electrons, num_embedded, num_scattered, total_time)
	end do

	! Close output files
	do i = 1, 7
		close(output_unit+i)
	end do

	! Compute total simulation time
	call cpu_time(end_time)
	total_time = end_time - start_time

	! Print simulation end information to console
	print*, "SIMULATION COMPLETE!"
	print*, "  Number of electrons embedded: ", num_embedded
	print*, "  Number of electrons scattered:", num_scattered
	print*, "Total simulation time:", total_time

	! Save final simulation status information to output file
	call write_simulation_status_information &
		(output_unit, i, num_electrons, num_embedded, num_scattered, total_time)

	! Deallocate allocatable arrays used in the simulation
	deallocate(electron_positions)
	deallocate(electron_velocities)
	deallocate(electron_accelerations)
	deallocate(atom_positions)
	deallocate(atomic_numbers)
	deallocate(atomic_numbers_cbrt)
	deallocate(embedded_positions)
	deallocate(scattered_positions)

end program main