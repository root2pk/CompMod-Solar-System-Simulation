"""
PROJECT A : N- Body Astronomical Simulation
This program will simulate an N-body astronomical interaction through Newtonian gravity
using the Velocity Verlet time integration algorithm to update values of position and velocity over time.
It will describe the Solar System, including the Sun, Pluto, Earth’s moon and the Halley’s Comet.

The bodies interact under Newtoniaan Gravity F = Gm1m2/r^2

The program produces a plot of the Total Energy of the system against timestep
It also produces two output files, one describing the trajectories of the bodies over time,
and the other with the apoapsis,peripasis and orbital periods of the bodies.

Units:
Mass : Kilogram
Distance : AU (Astronomical Units) ( 1 AU = 1.496E+11)
Time : days (1 day = 86400 s)

"""

import sys
import math
import numpy as np
import matplotlib.pyplot as pyplot
from Astro_Body import Astro_Body

def apoapsis(dist_list):
    """
        Function to return the value of apoapsis of a particular body

        :param dist_list: List of distances of body from centre of attraction
    """

    return max(dist_list)           # Maximum distance away from centre of attraction

def peripasis(dist_list):
    """
        Function to return the value of periapsis of a particular body

        :param dist_list: List of distances of body from centre of rotation
    """
    return min(dist_list)          # Minimum distance away from centre of attraction

def orbit_period(dist_list,time_list):
    """
        Function to return the value of orbital period of a particular body

        :param dist_list: List of distances of body from centre of rotation
        :param time_list: Corresponding time list of body
    """

    # If the list ascends first
    if dist_list[1] > dist_list[0]:

        # Find first peak in the list
        for i in range(len(dist_list)):
            if dist_list[i+1]<dist_list[i]:
                peak = i
                break
        # Find first trough after first peak in the list
        for i in range(peak,len(dist_list)):
            if dist_list[i+1]>dist_list[i]:
                trough = i
                break

    # If the list descends first
    else:

        # Find first trough in the list
        for i in range(len(dist_list)):
            if dist_list[i+1]>dist_list[i]:
                trough = i
                break
        # Find first peak after first trough in the list
        for i in range(trough,len(dist_list)):
            if dist_list[i+1]<dist_list[i]:
                peak = i
                break

    # Orbital period is twice the time difference betwen crest and trough
    return 2*abs(time_list[trough] - time_list[peak])

###                 MAIN FUNCTION         ###

def main():
    # Read name of output file from command line
    if len(sys.argv)!=4:
        print("Wrong number of arguments.")
        print("Usage: " + sys.argv[0] + " <input file 1>" + "<input file 2>")
        quit()
    else:
        infile1_name = sys.argv[1]
        infile2_name = sys.argv[2]
        outfile_name = sys.argv[3]

    # Open output file
    outfile = open(outfile_name, "w")              # Trajectory file
    outfile2 = open("Observables.txt","w")         # Observables file

    # Open input files
    infile1 = open(infile1_name, "r")              # particle details file
    infile2 = open(infile2_name, "r")              # simulation parameters file


    # Set up simulation parameters
    lines = infile2.readlines()
    dt = float(lines[0])                           # Time step dt
    numstep = int(lines[1])                        # number of steps
    time = float(lines[2])                         # start time

    # Check if there are enough days to calculate Pluto orbit (largest orbit)
    if (numstep*dt) < 80000:
        print("Not enough number of steps to calculate orbital period of all bodies")
        print("Increase numstep or timestep(dt)")
        print("numstep*dt must be greater than 80000")
        quit()


    # Set up astronomical body parameters
    P = Astro_Body.read_from_file(infile1)

    # Centre of Mass correction
    Astro_Body.com_correction(P)

    # Counting the number of bodies being simulated
    no = len(P)

    # Find iniitial conditions
    KE = 0
    # adding Individual kintetic energies
    for i in range(len(P)):
        KE = KE + Astro_Body.kinetic_energy(P[i])
    # Total eenrgy of the system
    energy =  Astro_Body.potential_energy(P) + KE

    # Initialise data lists for plotting later
    time_list = [time]
    energy_list = [energy]
    force = []
    dist_list = np.zeros((no,numstep+1))   # Array to store distance values for observable calculations

    # Finding out the index numbers of the Moon, Earth and Sun (For observable calculation)
    for i in range(no):
        if P[i].label == "Moon":
            ind_Moon = i
        if P[i].label == "Earth":
            ind_Earth = i
        if P[i].label == "Sun":
            ind_Sun = i

    # Finding the initial values of distance from centre of attraction for all bodies
    for i in range(no):
        # For moon, check orbit around Earth
        if i == ind_Moon:
            pos = P[i].position - P[ind_Earth].position
            dist = np.linalg.norm(pos)
            dist_list[ind_Moon,0] = dist
        # For others, check orbit around Sun
        else:
            pos = P[i].position - P[ind_Sun].position
            dist = np.linalg.norm(pos)
            dist_list[i,0] = dist


    ##########      SIMULATION       ###########


    # Calculation of Initial Forces

    F_matrix = Astro_Body.compute_force(P)
    for i in range(no):
        F = np.array([0,0,0])
        for j in range(no):
            F = F + np.array([F_matrix[i,j,0],F_matrix[i,j,1],F_matrix[i,j,2]])
        force.append(F)

    # Start the time integration loop

    # For loop to iterate over timesteps
    for n in range(numstep):

        Tot_KE = 0.0
        force_new = []
        #For loop to iterate over N astronomical bodies
        # Update Body Position
        for i in range(no):
            P[i].leap_position(dt, force[i])

        # Update forces
        F_matrix = Astro_Body.compute_force(P)
        for i in range(no):
            F = np.array([0,0,0])
            for j in range(no):
                F = F + np.array([F_matrix[i,j,0],F_matrix[i,j,1],F_matrix[i,j,2]])
            force_new.append(F)

        # Update body velocity by averaging
        # current and new forces
        for i in range(no):
            P[i].leap_velocity(dt, 0.5*(force[i]+force_new[i]))


        for i in range(no):
            force[i] = force_new[i]                   # Update new force value
            Tot_KE += Astro_Body.kinetic_energy(P[i]) # Find total kinetic energy of system

            # Add to distance lists for each body
            # Moon condition
            if i == ind_Moon:
                pos = P[i].position - P[ind_Earth].position
                dist = np.linalg.norm(pos)
                dist_list[ind_Moon,n+1] = dist
            # Others condition
            else:
                pos = P[i].position - P[ind_Sun].position
                dist = np.linalg.norm(pos)
                dist_list[i,n+1] = dist


        # Increase time
        time +=dt

        # Total Energy of the system
        energy = Tot_KE + Astro_Body.potential_energy(P)

        # Append information to data lists
        time_list.append(time)
        energy_list.append(energy)

        # Writing out trajectory file for VMD
        outfile.write(str(no) + "\n")
        outfile.write("Point = %d\n" % (j + 1))
        for i in range(no):
            outfile.write(str(P[i]) + "\n")


    #############   POST - SIMULATION  #####################

    # Write Observable values into file

    for i in range(no):
        # Ensure we're not calculating for the Sun
        if i != ind_Sun:
            apo = apoapsis(dist_list[i,:])                    # Apoapsis
            peri = peripasis(dist_list[i,:])                  # Periapsis
            OP = orbit_period(dist_list[i,:],time_list)       # Orbital Period
            #outfile2.write(str(P[i].label)+" - Apo-apsis : "+str(apo)+" AU, Peri-apsis : "+str(peri)+" AU, Orbital Period : "+str(OP)+" days ("+str(OP/365)+" years)\n")
            outfile2.write(str(P[i].label)+" - Apo-apsis : {0:.4f} AU, Peri-apsis : {1:.4f} AU, Orbital Period : {2:.2f} days ({3:.5f} years)\n".format(apo,peri,OP,OP/365))
    # Close all files
    infile1.close()
    infile2.close()
    outfile.close()
    outfile2.close()

    # Display average value of Total Energy of the System

    print("Average value of Total Energy of System : " + str(sum(energy_list)/len(energy_list)) + " kg*(AU^2)/(days^2)")

    # Plot total system energy to screen
    pyplot.title('Total Energy of the system vs. time')
    pyplot.xlabel('Time (days)')
    pyplot.ylabel('Energy kg*(AU^2)/(days^2)')
    pyplot.plot(time_list, energy_list)
    pyplot.show()


# Execute main method, but only when directly invoked
if __name__ == "__main__":
    main()
