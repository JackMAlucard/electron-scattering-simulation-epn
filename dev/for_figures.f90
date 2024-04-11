!####################    AFTER call setup_electron_beam_model of refactored main   ############################

! ***************  Storing values to file FOR TESTING  ***********************
	open(unit=output_unit+10, file='haz.dat', status='replace', action='write')
		! Initialize beam model with beam_target_distance = 0 and grazing_angle = 0
		! call...
		do i = 1, num_electrons
			r = electron_positions(i,:)
			v = electron_velocities(i,:)
			a = electron_accelerations(i,:)
			write(output_unit+10, *) r, v, a, 0
		end do
		write(output_unit+10, *)
		write(output_unit+10, *)
		! Initialize beam model with beam_target_distance and grazing_angle to use
		! call...
		do i = 1, num_electrons
			r = electron_positions(i,:)
			v = electron_velocities(i,:)
			a = electron_accelerations(i,:)
			write(output_unit+10, *) r, v, a, 0
		end do
	close(output_unit+10)
	! UNCOMMENT/COMMENT TO SEE/HIDE PLOTS
!	call system('gnuplot "haz.gp"')
!	call system('gnuplot "vel.gp"')
	print*, "Testing dummy input"
	read(*,*)
	! ***************  Storing values to file FOR TESTING  ***********************
	

!####################    AFTER call setup_simple_silica_model of refactored main   ############################

! ***************  Storing values to file FOR TESTING  ***********************
	! Storing Si and O positions in different files
	! call...
	open(unit=output_unit+20, file='material_Si_positions.dat', status='replace', action='write')
	open(unit=output_unit+30, file='material_O_positions.dat', status='replace', action='write')
		mbi = material_boundaries(1)
		mbj = material_boundaries(2)
		mbk = material_boundaries(3)
		do i = -mbi, mbi
			do j = 0, -mbj, -1
				do k = -mbk, mbk
					if (atom_charges(i,j,k) .eq. 14) write(output_unit+20, *) atom_positions(i,j,k,:)
					if (atom_charges(i,j,k) .eq. 8) write(output_unit+30, *) atom_positions(i,j,k,:)
				end do
			end do
		end do
	close(output_unit+20)
	close(output_unit+30)
	! PLOT THAT ONLY SHOWS SSM STRUCTURE
!	call system('gnuplot "ssm.gp"')
	! INFO PRINT
	print*, "Interatomic distance (a0):", INTERATOMIC_DIST_SIO2
	print*, "Material length boundaries (10):"
	print*, "x: +/-", material_boundaries(1)*INTERATOMIC_DIST_SIO2
	print*, "y:    ", material_boundaries(2)*INTERATOMIC_DIST_SIO2
	print*, "z: +/-", material_boundaries(3)*INTERATOMIC_DIST_SIO2
	print*, "Macroscopic cross section SiO2 (a0):", CROSS_SECTION_SIO2
	print*, "Mean free path SiO2 (a0):", MEAN_FREE_PATH_SIO2
	! ***************  Storing values to file FOR TESTING  ***********************
	
	
!####################    AFTER closing all files, just before deallocaring arrays in refactored main   ############################

! ***************  Storing values to file FOR TESTING  ***********************
	!PLOT THAT SHOWS SIMULATED EMBEDDED, SCATTERED ELECTRONS AND SCC STRUCTURE 	
	!call system('gnuplot "particles.gp"')
	
	!PLOT THAT SHOWS ELECTRON TRAJECTORIES 	
	!call system('gnuplot "beam_trajectories_zx.gp"')
	!call system('gnuplot "beam_trajectories_zy.gp"')
	! ***************  Storing values to file FOR TESTING  ***********************
	

