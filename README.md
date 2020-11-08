# CompMod
Computer Modelling Projects - Year 3

This program will simulate an N-body astronomical interaction through Newtonian gravity,
using the Velocity Verlet time integration algorithm to update values of position and velocity over time.
It will describe the Solar System, including the Sun, Pluto, Earth’s moon and the Halley’s Comet.

                  -------------------------------------------------------------

Units : 
Mass : Kilogram
Distance : AU (Astronomical Units) ( 1 AU = 1.496E+11)
Time : days (1 day = 86400 s)

                  -------------------------------------------------------------

Instructions
 - Unzip the folder and find the files
	- ParticleManyBody.py : Contains the main execution of the program
	- Astro_Body.py  : contains the class for the astronomical body objects
	- param.txt : Input file containing the simulation parameters in the order 

	  timestep
	  numstep
	  start time

	- particle.txt : Input file containing the planet names positions,velocities and masses(label,x,y,z,vx,vy,vz,mass)
	  (Taken from http://ssd.jpl.nasa.gov/horizons.cgi on date '2458911.500000000 = A.D. 2020-Mar-03 00:00:00.0000 TDB')
	- traj.xyz : Output trajectory file for the VMD simulation
	- Observables.txt : Output file with the apoapsis,periapses and orbital periods for all bodies except Sun

 - Run the program from the terminal using the command
	$ python ParticleManyBody.py particle.txt param.txt traj.xyz
 - A plot of total energy of the system against time will appear and the observables and trajectory will be written into respective files
 - Run the VMD simulation from the terminal using the command
	$ vmd traj.xyz
 - Manipulate the representations within vmd as required

                  -------------------------------------------------------------

NOTE: The program needs to run for atleast 80,000 days of simulation time to produce a value for orbital period of Pluto.
      The product of timestep and numstep must be greater than 80,000 days.
	
				
