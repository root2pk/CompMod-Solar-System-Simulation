# Universal Gravitational constant( in AU^3 kg^-1 days^-2)
G = -1.48818E-34

import numpy as np
class Astro_Body(object):
    """
    Class to describe 3D particles.

    Properties:
    label(string) - Name of the object body
    position(array) - position vector along the x,y,z axis
    velocity(array) - velocity vector along the x,y,z axis
    mass(float) - body mass

    Methods:
    * Formatted output
    * Kinetic energy
    * First-order velocity update
    * Second order position updates

    Static methods:
    * Vector separation
    * Read From File
    * Force Calculation
    * Potential Energy Calculation
    * Centre of Mass Motion Correction

    """

    def __init__(self, lab, pos, vel, mass):
        """
        Initialise a Particle1D instance

        :param pos: position as numpy array
        :param vel: velocity as numpy array
        :param mass: mass as float
        """
        self.label = lab
        self.position = pos
        self.velocity = vel
        self.mass = mass


    def __str__(self):
        """
        Define output format.
        For particle p=(2.0, 0.5, 1.0) this will print as
        "x = 2.0, v = 0.5, m = 1.0"
        """
        return self.label + " " + str(self.position[0]) + " " + str(self.position[1]) + " " + str(self.position[2])


    def kinetic_energy(self):
        """
        Return kinetic energy as
        1/2*mass*vel^2
        """
        return 0.5*self.mass*(np.linalg.norm(self.velocity))**2


    # Time integration methods
    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)

        :param dt: timestep as float
        :param force: force on particle as float
        """
        self.velocity += dt*(force/self.mass)


    def leap_position(self, dt, force):
        """
        Second-order position update,
        r(t+dt) = r(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force as float
        """
        self.position += dt*self.velocity + 0.5*(dt**2)*(force/self.mass)


    #Static Methods
    @staticmethod

    def vect_sep(p1, p2):
        """
        Function to return the vector separation between two particles

        :param p1: position of body p1
        :param p2: position of body p2
        """
        return p2.position-p1.position


    def read_from_file(infile):
        """
        A static method to create a particle from a file. The method reads through the lines of the file.
        The form of the file lines should be:
        "label,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,mass"

        :param: A file that python can read through
        :return: Astro_Body object
        """
        P = []
        infile.seek(0)
        lines = infile.readlines()

        for line in lines:
            token = line.split(",")

            label = str(token[0])
            x = float(token[1])
            y = float(token[2])
            z = float(token[3])
            vx = float(token[4])
            vy = float(token[5])
            vz = float(token[6])
            m = float(token[7])

            body = Astro_Body(label,np.array([x,y,z]),np.array([vx,vy,vz]),m)
            P.append(body)

        return P


    def compute_force(P):
        """
            Method for computing the force matrix for the system
            at any point of time.

            Returns a 3D array of the form F[i,j,k], where F[i,j,:] is the
            force acting on a body i due to body j. k varies between 0,1 and 2;
            indicating the x,y and z components.
            Uses the idea that F[i,j,:] = -F[j,i,:],
            and F[i,i,:] = 0

            :param P: list of Astro_Body objects
        """

        # Number of bodies
        no = len(P)
        # Initialising the matrix
        F_matrix = np.zeros((no,no,3))

        for i in range(no):
            for j in range(no):
                # Regular computation
                if i<j:
                    vect_r = Astro_Body.vect_sep(P[j],P[i])
                    r = np.linalg.norm(vect_r)
                    F = ((G*P[i].mass*P[j].mass)/r**3)*vect_r
                    F_matrix[i,j,:] = F

                # Force acting on a body due to itself
                elif i==j:
                    F_matrix[i,j,:] = [0,0,0]

                # Simplifying Condition
                else:
                    F_matrix[i,j,:] = -F_matrix[j,i,:]

        return F_matrix

    def potential_energy(P):

        """
            Method used to calculate the total potential energy of the system

            This method calculates a unuique potential energy interaction between each body
            and returns the scalar sum.

            :param P: List of Astro_Body objects
        """

        # Initalising the total potential energy = 0
        PE = 0
        # Number of bodies
        no = len(P)

        for i in range(no):
            for j in range(no):
                # To eliminate double counting
                if i<j:
                    r = np.linalg.norm(Astro_Body.vect_sep(P[i],P[j]))
                    PE += (G*P[i].mass*P[j].mass)/r

        return PE

    def com_correction(P):
        """
            Method to correct the values for initial kinetic energy for
            list of Astro_Body objects by reducing the velocity of the
            centre of mass of the system from all bodies.

            :param P: List of Astro_Body objects

        """
        # Number of bodies
        no = len(P)
        # Initialising momentum to be 0
        momentum = np.array([0,0,0])
        # Initialising sum of masses to be 0
        mass_sum = 0.0
        for i in range(no):

            momentum = momentum + P[i].mass*P[i].velocity
            mass_sum += P[i].mass

        # Velocuty of centre of mass
        vcom = momentum/mass_sum
        # Subtracting the velocity of the centre of mass from each body
        for i in range(no):
            P[i].velocity = P[i].velocity - vcom
