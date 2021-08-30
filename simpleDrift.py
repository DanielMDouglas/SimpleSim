import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import sys

from utils import *
from parameters import *

def sample_diffused_pdf(v):
    """
    get the length of the step in the axial (z) and radial (rho) directions.
    These directions are defined in terms of the local E field
    """
    # get a random step (in rho & z) 
    DT = physics_parameters["DT"]
    DL = physics_parameters["DL"]

    dt = sim_parameters["dt"]

    rho = st.norm.rvs(scale = np.sqrt(2*8*DL*dt))
    z = st.norm.rvs(loc = v*dt, scale = np.sqrt(2*8*DT*dt))
    return rho, z

def sample_azimuthal_pdf():
    """
    get the azimuthal (where the axial is the drift direction) direction for a drift step.
    There is no preferred direction, so the distribution is flat
    """
    theta = st.uniform.rvs(scale = 2*np.pi)
    return theta

def sample_from_cathode_target():
    """ 
    get a random position on the cathode target (8mm diameter on the cathode plane)
    """
    z = detector_parameters["cathode position"]
    R = detector_parameters["target radius"]

    rho = np.power(st.uniform.rvs()*np.power(R, 2), 0.5)
    phi = st.uniform.rvs(scale = 2*np.pi)
    x = rho*np.cos(phi)
    y = rho*np.sin(phi) 

    pos = np.array([x, y, z])
    return pos

def sample_from_hemisphere():
    """
    return a random unit vector where the z component is negative
    """
    az = 2*np.pi*st.uniform.rvs() # the distribution of azimuthal angles is flat
    pol = np.arccos(-st.uniform.rvs()) # the distribution of polar angles is weighted for even coverage per solid angle

    x = np.cos(az)*np.sin(pol)
    y = np.sin(az)*np.sin(pol)
    z = np.cos(pol)

    return np.array([x, y, z])
    
class Efield:
    def __init__(self, transv, longit):
        # TODO load Efield from map, get value from interp
        self.field = None
        self.transv = transv
        self.longit = longit
    def value(self, pos):
        return np.array([self.transv, 0., detector_parameters['nominal field'] + self.longit])
        
class charge:
    def __init__(self, pos_i):
        # TODO
        self.pos_i = pos_i # store initial position
        self.pos = pos_i # position vector which changes over time
        self.fate = None # is the electron still in play?
        self.history = []
        self.arrivalT = 0
    def drift(self, this_E, use_bulk = False):
        ## TODO: add bulk drift method, pass a flag to use bulk vs. brownian drift
        
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
            rho, z = sample_diffused_pdf(v_drift)
            phi = sample_azimuthal_pdf()

            dx = z*drift_direction + rho*(np.cos(phi)*perp_direction1 + np.sin(phi)*perp_direction2)
            
            self.history.append(self.pos)

            self.pos = self.pos + dx

            self.arrivalT += sim_parameters["dt"]
            
            # check for the particle to finish
            if self.pos[2] <= 0:
                self.fate = 1 # fate == 1 means the electron made it to the anode

class tracklet:
    def __init__(self, thisSource = 'Cs137'):
        self.pos = sample_from_cathode_target() # z of the cathode, x, y, sampled within source shape
        self.dir = sample_from_hemisphere() # az is flat in [0, 2pi], zenith is cos(zen) in [0, -pi/2]

        self.Ei = emission_spectrum(thisSource) # Initial energy [MeV]
        
        self.length = betaRange(self.Ei) # track length [cm]

        # dEdx = betadEdx(Ei)
        self.dEdx = self.Ei/self.length # deposited energy density [MeV/cm]

        # dQdx = recomb(dEdx)
        self.Qtot = self.Ei*physics_parameters["R"]/physics_parameters["w"] # total ionized charge [e]
        self.dQdx = self.dEdx*physics_parameters["R"]/physics_parameters["w"] # ioniization density [e/cm]

        print(self.Ei, self.length, self.dEdx, self.Qtot, self.dQdx)

    def generate_charge(self):
        """
        from the tracklet shape defined in the initializer, generate a charge
        """
        dist = st.uniform.rvs(0, self.length)

        thisPos = self.pos + self.dir*dist
        
        return charge(thisPos)
        

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', type = str,
                        default = 'driftHits.npy',
                        help = 'where to save the destinations')
    parser.add_argument('-t', '--transverse', type = float,
                        default = 0,
                        help = 'amount of transverse (x-direction) drift field to add')
    parser.add_argument('-l', '--longitudinal', type = float,
                        default = 0,
                        help = 'amount of longitudinal (z-direction) drift field to add')
    parser.add_argument('-N', '--N', type = int,
                        default = int(1e4),
                        help = 'number of photoelectrons to drift')
    parser.add_argument('-z', '--z0', type = float,
                        default = 50,
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

    thisTracklet = tracklet()
    nChargeBundles = int(thisTracklet.Qtot/sim_parameters["scalingF"])
    # for i in range(Npe):
    print ("drifting ", nChargeBundles, " discrete charge bundles")
    for i in range(nChargeBundles):
        
        this_charge = thisTracklet.generate_charge()

        print ("charge ", i, " originates at ", this_charge.pos)

        # starting_position = sample_from_cathode_target(z0 = args.z0)
        # this_charge = charge(starting_position)

        this_charge.drift(thisEfield)
        
        x, y, z = np.array(this_charge.history).T
        
        arrivalTimes.append(this_charge.arrivalT)
        finalXs.append(this_charge.pos[0])
        finalYs.append(this_charge.pos[1])
        finalZs.append(this_charge.pos[2])

        initXs.append(this_charge.pos_i[0])
        initYs.append(this_charge.pos_i[1])
        initZs.append(this_charge.pos_i[2])

        print ("charge ", i, " terminates at ", this_charge.pos)

    np.save(outFile,
            np.array([finalXs,
                      finalYs,
                      finalZs,
                      arrivalTimes]))
