import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import sys

from parameters import *
from utils import *
    
def sample_diffused_pdf(v):
    # get a random step (in rho & z) 
    DT = physics_parameters["DT"]
    DL = physics_parameters["DL"]

    dt = sim_parameters["dt"]

    rho = st.norm.rvs(scale = np.sqrt(2*8*DL*dt))
    z = st.norm.rvs(loc = v*dt, scale = np.sqrt(2*8*DT*dt))
    return rho, z

def sample_azimuthal_pdf():
    theta = st.uniform.rvs(scale = 2*np.pi)
    return theta

class Efield:
    def __init__(self, transv, longit):
        # TODO load Efield from map, get value from interp
        self.field = None
        self.transv = transv
        self.longit = longit
    def value(self, pos):
        return np.array([self.transv, 0., 0.5 + self.longit]) # nominal = 0.5 kV/cm
        
class charge:
    def __init__(self, pos_i):
        # TODO
        self.pos_i = pos_i # store initial position
        self.pos = pos_i # position vector which changes over time
        self.fate = None # is the electron still in play?
        self.history = []
        self.arrivalT = 0
    def drift(self, this_E):
        # for each timestep, sample the next position from the bulk drift PDF:
        # n(rho, z, t) = (n0/(4*pi*DT*t*sqrt(4*pi*DL*t)))*exp(-(z - v*t)^2/(4*DL*t) - lambda*v*t)*exp(-rho^2/(4*DT*t))
        #
        # z is drift direction (direction of E field)
        #
        # repeat until z position is at the cathode
        while self.fate == None:
            localEfield = this_E.value(self.pos)

            v_drift = physics_parameters["v"](localEfield)

            drift_direction = -norm(localEfield)
            perp_direction1, perp_direction2 = get_perpendicular_vectors(drift_direction)
            
            # get the coordinates for the displacement
            rho, z = sample_diffused_pdf(sim_parameters["dt"], v_drift)
            phi = sample_azimuthal_pdf()

            dx = z*drift_direction + rho*(np.cos(phi)*perp_direction1 + np.sin(phi)*perp_direction2)
            
            self.history.append(self.pos)

            self.pos = self.pos + dx

            self.arrivalT += sim_parameters["dt"]
            
            # check for the particle to finish
            if self.pos[2] <= 0:
                self.fate = 1 # fate of 1 means the electron made it to the anode

def sample_from_cathode_target(z0 = 50):
    # TODO
    # get a random position on the cathode target (8mm diameter on the cathode plane)
    # z = detector_parameters["cathode position"]
    z = z0
    R = detector_parameters["target radius"]

    rho = np.power(st.uniform.rvs()*np.power(R, 2), 0.5)
    phi = st.uniform.rvs(scale = 2*np.pi)
    
    x = rho*np.cos(phi)
    y = rho*np.sin(phi) 

    pos = np.array([x, y, z])
    return pos

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', type = str,
                        help = 'where to save the destinations')
    parser.add_argument('-t', '--transverse', type = float,
                        help = 'amount of transverse (x-direction) drift field to add')
    parser.add_argument('-l', '--longitudinal', type = float,
                        help = 'amount of longitudinal (z-direction) drift field to add')
    parser.add_argument('-N', '--N', type = int,
                        help = 'number of photoelectrons to drift')
    parser.add_argument('-z', '--z0', type = float,
                        help = 'nominal distance from cathode to anode')
    args = parser.parse_args()

    outFile = args.output
    
    thisEfield = Efield(args.transverse, args.longitudinal)
    
    # Npe = int(st.norm.rvs(loc = physics_parameters["npe"],
    #                       scale = physics_parameters["npe_sigma"]))

    # print (Npe)

    # Npe = 5000
    Npe = args.N
    
    arrivalTimes = []
    finalXs = []
    finalYs = []
    finalZs = []

    initXs = []
    initYs = []
    initZs = []
    
    fig3d = plt.figure()
    ax = fig3d.add_subplot(projection = '3d')
    
    for i in range(Npe):
        starting_position = sample_from_cathode_target(z0 = args.z0)
        this_charge = charge(starting_position)

        this_charge.drift(thisEfield)

        x, y, z = np.array(this_charge.history).T
        
        ax.plot(x, y, z)

        arrivalTimes.append(this_charge.arrivalT)
        finalXs.append(this_charge.pos[0])
        finalYs.append(this_charge.pos[1])
        finalZs.append(this_charge.pos[2])

        initXs.append(this_charge.pos_i[0])
        initYs.append(this_charge.pos_i[1])
        initZs.append(this_charge.pos_i[2])

        # print ( this_charge.arrivalT )


    ax.set_xlabel(r'y (cm)')
    ax.set_ylabel(r'z (cm)')
    ax.set_zlabel(r'x (cm)')
    
    # plt.figure()
    # plt.hist(arrivalTimes)
    # print (arrivalTimes)
    # print (np.min(arrivalTimes),
    #        np.mean(arrivalTimes),
    #        np.max(arrivalTimes))

    np.save(outFile,
            np.array([finalXs,
                      finalYs,
                      finalZs,
                      arrivalTimes]))

    # exp_wid = np.sqrt(4*physics_parameters["DT"]*np.mean(arrivalTimes))
    # print ("expected width", exp_wid)
    # exp_z_wid = 2*np.sqrt(physics_parameters["DL"]*np.mean(arrivalTimes))
    # print ("expected (timing) width", exp_z_wid)
    
    # plt.figure()
    # plt.hist(np.array(finalXs) - np.array(initXs),
    #          density = True)
    # plt.plot(np.linspace(-3e-1, 3.e-1, 1000),
    #          st.norm.pdf(np.linspace(-3.e-1, 3.e-1, 1000),
    #                      loc = 0,
    #                      scale = exp_wid))

    # plt.figure()
    # print (np.min(finalZs), np.mean(finalZs), np.max(finalZs))
    # plt.hist(np.array(finalZs),
    #          density = True)
    # plt.plot(np.linspace(-3e-1, 3.e-1, 1000),
    #          st.norm.pdf(np.linspace(-3.e-1, 3.e-1, 1000),
    #                      loc = 0,
    #                      scale = exp_z_wid))
    plt.show()
