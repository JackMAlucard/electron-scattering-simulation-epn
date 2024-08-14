!===============================================================================
! Program    : main_electron_scattering
! Purpose    : Simulate the behavior of an electron beam interacting with a
!              dielectric material by computing the trajectories of each
!              electron in the beam sequentially.
! Description:
!   This program executes the following main steps:
!     1. Set random seeds for simulation reproducibility.
!     2. Read simulation input parameters from input file.
!     3. Set up the electron beam model and dielectric material model.
!     4. Compute the trajectories of individual electrons sequentially.
!     5. Determine whether each electron gets embedded or is scattered.
!     6. Compute scattering angles for scattered electrons.
!     7. Save simulation results to output files for further analysis.
! Input      :
!   The following input parameters are read from the 'input.txt' file:
!     - beam_energy
!         Beam energy in kiloelectronvolts (keV), real
!		  - grazing_angle
!         Grazing angle in degrees (º), real
!     - num_electrons
!         Number of electrons in the beam, integer
!     - energy_spread
!         Beam energy spread as a percentage of the beam energy (%), real
!     - spot_size_factor
!         Beam spot size scale factor, real
!     - beam_target_distance
!         Beam-target material distance in angstroms (Å), real
!     - material_boundaries
!         Material model grid index boundaries, three integers
!     - num_plot_ploints
!         Number of points to be plotted for each electron trajectory for
!         when electron trajectories saving flag is enabled, integer
!     - dt
!         Time step size in atomic units of time (aut), real
!     - beam_model_output_saving_enabled
!         Flag for beam model output saving, boolean
!		  - material_model_output_saving_enabled
!         Flag for material model output saving, boolean
!		  - electron_trajectories_saving_enabled
!         Flag for electron trajectories saving, boolean
! Output     :
!   The output values corresponding to simulated magnitudes are all in atomic
!   units. The output files and the format of each of their lines are:
!     - electron_trajectories.dat (optional)
!         Electron trajectories, if trajectory saving is enabled. Each
!         trajectory in the file is separated by two blank lines.
!         Output format:
!           time, x-coordinate, y-coordinate, z-coordinate
!     - embedded_positions.dat
!         Final positions of embedded electrons. Output format:
!           electron index, x-coordinate, y-coordinate, z-coordinate
!     - scattered_positions.dat
!         Final positions of scattered electrons. Output format:
!           electron index, x-coordinate, y-coordinate, z-coordinate
!     - scattering_angles.dat
!         Scattering angles of scattered electrons. Output format:
!           electron index, elevation angle alpha, azimuthal angle beta
!     - final_electron_positions.dat
!         Final positions of all electrons. Output format:
!           electron index, x-coordinate, y-coordinate, z-coordinate
!     - final_state_evolution.dat
!         Evolution of the final state of the electrons in the beam.
!         Output format:
!           electron index, number of electrons embedded up to current electron,
!           number of electrons scattered up to current electron,
!           fraction of electrons embedded up to current electron,
!           fraction of electrons scattered up to current electron
!     - simulation_times.dat
!         Simulation time of each electron trajectory and total simulation time.
!         Output format:
!           electron index, current electron trajectory simulation time,
!           total simulation time up to current electron
!     - simulation_info.dat
!         Simulation status information. It follows a specific text format.
!===============================================================================
program main_electron_scattering

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
	! PROGRAM VARIABLES DECLARATION
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

end program main_electron_scattering