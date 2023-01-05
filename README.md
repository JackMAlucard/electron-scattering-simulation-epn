# Simulation-of-Elastic-Scattering-Processes-of-1-keV-to-10-keV-electrons-in-dielectric-materials

## Main Electron Trajectories Simulation

The source code for the main simulation can be found in the 'src' folder.

### Code Structure

The simulation code consists of five subroutines and the main program.

- Module 0 (m0): Contains numerical storage size parameters definition for real and integer variables, as well as subroutines for Pseudo-Random Numbers (PRNs) generation, unit conversion, and vector-matrix product.   
- Module 1 (m1): Contains subroutines used to generate and locate the electron beam.
- Module 2 (m2): Contains subroutines used to generate the dielectric material model atom locations and charge.
- Module 3 (m3): Contains subroutines used to: compute the electron acceleration due to other electrons and neutral atoms,; the position, velocity, and acceleration for each time step using the Velocity Verlet Algorithm; and finally the complete electron trajectory. 
- Module 4 (m4): Contains subroutines used to optimize the computation of electron trajectories.

### Compilation

To run a simulation:

1. Compile each of the modules with the command:
	`gfortran -c m0.f90 m1.f90 m2.f90 m3.f90 m4.f90`
2. Compile the main program with the command:
	`gfortran m0.o m1.o m2.o m3.o m4.o main.f90 -o main.o`
3. Run the main program executable generated 'main.o' with the command:
	`./main.o`

At the moment, the simulation parameters are set within the 'main.f90' program. This file needs to be modified in order to change them if necessary.