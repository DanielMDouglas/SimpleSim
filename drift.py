import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import sys

import h5py
from scipy.interpolate import RegularGridInterpolator

from utils import *
from parameters import *
from eventRecord import *

def import_hdf5(file):
    f = h5py.File(file,'r')
    a_group_key = list(f.keys())[0]
    d = list(f[a_group_key])
    return np.array(d)

def sample_diffused_pdf(v):
    """
    get the length of the step in the axial (z) and radial (rho) directions.
    These directions are defined in terms of the local E field
    """
    # get a random step (in rho & z)
    DT = physics_parameters["DT"]
    DL = physics_parameters["DL"]

    dt = sim_parameters["dt"]

    # rho = st.norm.rvs(scale = np.sqrt(2*8*DL*dt))
    # x is the drift direction
    x = st.norm.rvs(loc=v*dt, scale=np.sqrt(4*DL*dt))
    y = st.norm.rvs(scale=np.sqrt(4*DT*dt))
    z = st.norm.rvs(scale=np.sqrt(4*DT*dt))
    return x, y, z


def sample_azimuthal_pdf():
    """
    get the azimuthal (where the axial is the drift direction) direction for a drift step.
    There is no preferred direction, so the distribution is flat
    """
    theta = st.uniform.rvs(scale=2*np.pi)
    return theta


def sample_from_cathode_target(yC = None, zC = None):
    """ 
    get a random position on the cathode target (8mm diameter on the cathode plane)
    """
    x = -1.0
    # x = detector_parameters["cathode position"] 
    R = detector_parameters["target radius"]

    rho = np.power(st.uniform.rvs()*np.power(R, 2), 0.5)
    phi = st.uniform.rvs(scale=2*np.pi)

    if yC == None:
        yC = 0
    if zC == None:
        zC = 0
    
    z = rho*np.cos(phi) + zC
    y = rho*np.sin(phi) + yC

    pos = np.array([x, y, z])
    return pos


def sample_from_hemisphere():
    """
    """
    az = 2*np.pi*st.uniform.rvs()  # the distribution of azimuthal angles is flat
    # the distribution of polar angles is weighted for even coverage per solid angle
    pol = np.arccos(-st.uniform.rvs())

    # Beam z, Zenith y, x drift
    x = np.sin(az)*np.sin(pol)
    y = np.cos(pol)
    z = np.cos(az)*np.sin(pol)

    # Standard spherical coords
    # x = np.cos(az)*np.sin(pol)
    # y = np.sin(az)*np.sin(pol)
    # z = np.cos(pol)

    return np.array([x, y, z])


def sample_from_face():
    """
    Return a position from a rectangle near the upstream face 
    """

    # Make sampling face larger than actual detector upstream face to ensure 
    # that some rock muons come in from the sides. 
    x0 = np.random.uniform(-10, 40, 1)[0] #(0,30) detector bound
    y0 = np.random.uniform(-25, 25, 1)[0] #(-15,15) detector bound
    z0 = -15.1  # 0.1 cm behind the upstream face (located at -15 cm)

    return np.array([x0, y0, z0])

def sample_from_bounding_sphere():
    """
    return a random position on a sphere that 
    """
    az = 2*np.pi*st.uniform.rvs()  # the distribution of azimuthal angles is flat
    # the distribution of polar angles is weighted for even coverage per solid angle
    pol = np.arccos(1 - 2*st.uniform.rvs())
    r = sim_parameters["generation sphere radius"]
    center = sim_parameters["generation sphere center"]

    # Beam z, Zenith y, x drift
    x = r*np.sin(az)*np.sin(pol)
    y = r*np.cos(pol)
    z = r*np.cos(az)*np.sin(pol)

    # Standard spherical coords
    # x = r*np.cos(az)*np.sin(pol)
    # y = r*np.sin(az)*np.sin(pol)
    # z = r*np.cos(pol)

    return np.array([x, y, z]) + center


atm_spect = np.loadtxt('atm_flux.dat')


def sample_from_CR_spectrum():
    ind = np.random.choice(atm_spect.shape[0])
    E, zen = atm_spect[ind]
    return 1.e3*E, zen


beam_spect = np.loadtxt('beam_flux.dat')


def sample_from_beam_spectrum():
    ind = np.random.choice(beam_spect.shape[0])
    E, zen, az = beam_spect[ind]
    return 1.e3*E, zen, az



class Efield:
    def __init__(self):

        f_root_read = TFile(f"yifan/Efield_TH3_raw.root","READ")
        # f_root_read = TFile(f"Efield_TH3_raw.root","READ")
        # Define input TH3F's (efield)
        self.Ex = f_root_read.Get("th3_Ex_raw")
        self.Ey = f_root_read.Get("th3_Ey_raw")
        self.Ez = f_root_read.Get("th3_Ez_raw")

        self.Ex.SetDirectory(ROOT.nullptr)
        self.Ey.SetDirectory(ROOT.nullptr)
        self.Ez.SetDirectory(ROOT.nullptr)

        # self.Ex.GetXaxis()
        # GetNbins()
        # GetBinLowEdge(n)
        # GetBinCenter(n)

        # print( self.Ex )  

    def value(self, pos):

        #From Dan (John's coord system)
        #x is [-0.3225, 0.3225], y is [-0.625, 0.625], z is [-0.3, 0.3]

        x,y,z = pos/100 # Convert cm to m
        print(x,y,z)
        # print(self.Ex)
        # print( self.Ex.Interpolate(x,y,z) )
        # Calculate at some position (flip z and x between our and John's frameworks)
        E1 = self.Ex.Interpolate(z,y,x)
        E2 = self.Ey.Interpolate(z,y,x)
        E3 = self.Ez.Interpolate(z,y,x)
        E_v = np.array([E3, E2, E1]) #kV/m
        E_v = E_v/100
        # print(E_v)
        return E_v #in kV/cm

# class Efield:
#     def __init__(self, transv, longit, pert_files):

#         ex_map, ey_map, ez_map = pert_files
#         # ex_map = 'ex.hdf5'
#         # ey_map = 'ey.hdf5'
#         # ez_map = 'ez.hdf5'
#         f_x = import_hdf5(ex_map)    
#         f_y = import_hdf5(ey_map)
#         f_z = import_hdf5(ez_map)

#         #Number of grid points in the map (same for each component of the electric field)
#         nx, ny, nz = f_x.shape

#         #Set up the same grid points for the Python interpolator 
#         x_range = detector_parameters['detector bounds'][0]
#         y_range = detector_parameters['detector bounds'][1]
#         z_range = detector_parameters['detector bounds'][2] 

#         x = np.linspace(x_range[0], x_range[1], nx)
#         y = np.linspace(y_range[0], y_range[1], ny)
#         z = np.linspace(z_range[0], z_range[1], nz)

#         #Interpolate 
#         self.dE_x = RegularGridInterpolator((x, y, z), f_x)
#         self.dE_y = RegularGridInterpolator((x, y, z), f_y)
#         self.dE_z = RegularGridInterpolator((x, y, z), f_z)


#         self.field = None
#         self.transv = transv
#         self.longit = longit

#     def value(self, pos):
#         # x drift
#         # flat component
#         flatField = np.array( [detector_parameters['nominal field'] + self.longit, 0., self.transv] )

#         # add a perturbation

#         # a ball of charge
#         # Q = 10.
#         # R = 5.
#         # c = detector_parameters['detector center']
#         # displ = pos - c
#         # dist = mag(displ)

#         # if dist < R:
#         #     pertField = Q/np.power(R, 3)*displ
#         # else:
#         #     pertField = Q/np.power(dist, 3)*displ

#         #Taken from loaded, interpolated field 
#         pertField = np.array( [self.dE_x(pos), self.dE_y(pos), self.dE_z(pos)] ).flatten()
        
#         return flatField + pertField
#         # return flatField


class charge:
    def __init__(self, pos_i):
        # TODO
        self.pos_i = np.array(pos_i)  # store initial position
        self.pos = np.array(pos_i)  # position vector which changes over time
        self.fate = None  # is the electron still in play?
        self.history = []
        self.arrivalT = 0
        self.weight = sim_parameters["scalingF"]

    def drift(self, this_E, drift_model='randomWalk'):
        self.history.append(self.pos)

        if not is_in_tpc(self.pos):
            self.fate = 2  # fate == 2 means the electron was generated outside of the tpc volume

        # for each timestep, sample the next position from the bulk drift PDF:
        # n(rho, z, t) = (n0/(4*pi*DT*t*sqrt(4*pi*DL*t)))*exp(-(z - v*t)^2/(4*DL*t) - lambda*v*t)*exp(-rho^2/(4*DT*t))
        #
        # z is drift direction (direction of E field)
        #
        # repeat until z position is at the cathode
        if drift_model == 'bulk':
            # a simplified drift model that runs in a single step
            # assumes the field is constant everywhere
            # and the angle of deflection is not too great

            while self.fate == None:

                localEfield = this_E.value(self.pos)

                v_drift = physics_parameters["v"](localEfield)
                v = mag(v_drift)

                t0 = self.pos[0]/v

                x = 0
                z = st.norm.rvs(loc=self.pos[1], scale=np.sqrt(
                    4*physics_parameters["DT"]*t0))
                y = st.norm.rvs(loc=self.pos[2], scale=np.sqrt(
                    4*physics_parameters["DT"]*t0))
                t = st.norm.rvs(loc=t0, scale=np.sqrt(
                    4*physics_parameters["DL"]*t0)/v)

                self.pos = np.array([x, y, z])
                self.history.append(self.pos)
                self.arrivalT = t

                self.fate = 1

            self.weight *= np.exp(-self.arrivalT/physics_parameters['lt'])

        elif drift_model == 'randomWalk':
            while self.fate == None:
                localEfield = this_E.value(self.pos)
                # print("Field")
                # print(localEfield)

                v_drift = physics_parameters["v"](localEfield)

                # drift_direction = norm(localEfield) #For John's solution 
                drift_direction = -norm(localEfield) #For our solution (drift direction opposite)
                perp_direction1, perp_direction2 = get_perpendicular_vectors(
                    drift_direction)

                # get the coordinates for the displacement
                dx, dy, dz = sample_diffused_pdf(v_drift)
                phi = sample_azimuthal_pdf()

                dl = dx*drift_direction + dz*perp_direction1 + dy*perp_direction2

                self.pos = self.pos + dl
                self.history.append(self.pos)

                self.arrivalT += sim_parameters["dt"]

                # check for the particle drifted to the anode
                if self.pos[0] <= -29: # For John's solution
                # if self.pos[0] <= 0: # For our work (we flipped direction of drift)
                    self.fate = 1  # fate == 1 means the electron made it to the anode
                    self.weight *= np.exp(-self.arrivalT/physics_parameters['lt'])

        else:
            raise ValueError(drift_model + " is not a valid drift model!")


class radSource:
    def __init__(self, thisSource='Cs137'):
        self.source = thisSource

    def generate_tracklet(self):
        # z of the cathode, x, y, sampled within source shape
        pos = sample_from_cathode_target()
        # az is flat in [0, 2pi], zenith is cos(zen) in [0, -pi/2]
        dir = sample_from_hemisphere()
        Ei = emission_spectrum(self.source)  # Initial energy [MeV]

        length = physics_parameters["e-Ar range"](Ei)  # track length [cm]

        # dEdx = betadEdx(Ei)
        dEdx = Ei/length  # deposited energy density [MeV/cm]

        # dQdx = recomb(dEdx)
        # total ionized charge [e]
        Qtot = Ei*physics_parameters["R"]/physics_parameters["w"]
        # ioniization density [e/cm]
        dQdx = dEdx*physics_parameters["R"]/physics_parameters["w"]

        return tracklet(pos, dir, Ei, length)


class tracklet:
    def __init__(self, pos, dir, edep, length):
        self.pos = pos
        self.dir = dir

        self.edep = edep

        self.length = length

        # dEdx = betadEdx(Ei)
        self.dEdx = self.edep/self.length  # deposited energy density [MeV/cm]

        # dQdx = recomb(dEdx)
        # total ionized charge [e]
        self.Qtot = self.edep*physics_parameters["R"]/physics_parameters["w"]
        # ioniization density [e/cm]
        self.dQdx = self.dEdx*physics_parameters["R"]/physics_parameters["w"]

    def generate_charge(self):
        """
        from the tracklet shape defined in the initializer, generate a charge
        """
        dist = st.uniform.rvs(0, self.length)

        thisPos = self.pos + self.dir*dist

        return charge(thisPos)


class muonTrack:
    def __init__(self, pos, dir, Ei):
        self.pos = pos
        self.dir = dir
        self.Ei = Ei

        self.segment_length = sim_parameters['dx']

        # track length [cm]
        self.length = physics_parameters["mu-Ar range"](self.Ei)

        self.dE = [] # deposited energy per segment [MeV]
        Eseg = self.Ei
        while Eseg >= 0:
            dEdx = physics_parameters['mu-Ar dEdx'](Eseg)
            dE = dEdx*self.segment_length
            self.dE.append(dE)
            Eseg -= dE

    def generate_segments(self):
        """
        from the track shape defined in the initializer, generate a charge
        """

        nSegments = int(self.length/self.segment_length)

        trackletList = []

        for i, dE in enumerate(self.dE):
            pos = self.pos + self.dir*i*self.segment_length
            dir = self.dir
            edep = 1

            trackletList.append(tracklet(pos, dir, dE, self.segment_length))

        return trackletList

class cosmicRayTrack:
    def __init__(self, condition):
        self.throw_pos_dir()
        while not condition(self):
            self.throw_pos_dir()
        
        self.track = muonTrack(self.pos, self.dir, self.Ei)

    def throw_pos_dir(self):
        self.pos = sample_from_bounding_sphere()
        self.Ei, zen = sample_from_CR_spectrum()
        az = 2*np.pi*st.uniform.rvs()

        # z beam, y zenith, x drift
        self.dir = np.array([np.sin(az)*np.sin(zen),
                             np.cos(zen),
                             np.cos(az)*np.sin(zen)])
class rockMuonTrack:
    def __init__(self, condition):
        self.throw_pos_dir()
        while not condition(self):
            self.throw_pos_dir()

        self.track = muonTrack(self.pos, self.dir, self.Ei)

    def throw_pos_dir(self):
        self.pos = sample_from_face()
        self.Ei, zen, az = sample_from_beam_spectrum()

        self.dir = np.array([np.sin(az)*np.sin(zen),
                             np.cos(zen),
                             np.cos(az)*np.sin(zen)])

class isotropicMuonTrack:
    def __init__(self, condition):
        self.throw_pos_dir()
        # if the track is not poniting inwards, try again
        while not condition(self):
            self.throw_pos_dir()

        self.track = muonTrack(self.pos, self.dir, self.Ei)

    def throw_pos_dir(self):
        self.pos = sample_from_bounding_sphere()
        self.Ei = 1000

        az = 2*np.pi*st.uniform.rvs() 
        zen = np.arccos(1 - 2*st.uniform.rvs())

        # z beam, y zenith, x drift
        self.dir = np.array([np.sin(az)*np.sin(zen),
                             np.cos(zen),
                             np.cos(az)*np.sin(zen)])


if __name__ == '__main__':
    import argparse
    # from ROOT import TFile, TH2F, TH3F, TCanvas, TLegend, gStyle, TGaxis

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output', type=str,
                        default='driftHits.npy',
                        help='where to save the destinations')
    # parser.add_argument('-t', '--transverse', type=float,
    #                     default=0,
    #                     help='amount of transverse (z-direction) drift field to add')
    # parser.add_argument('-l', '--longitudinal', type=float,
    #                     default=0,
    #                     help='amount of longitudinal (x-direction) drift field to add')

    # parser.add_argument('-p', '--pertFields', nargs='+',
    #                     default=['ex.hdf5','ey.hdf5','ez.hdf5'],
    #                     help='names of 3 files with field perturbations')

    parser.add_argument('-N', '--N', type=int,
                        default=int(1e2),
                        help='number of photoelectrons to drift')
    parser.add_argument('-z', '--z0', type=float,
                        default=15, #default is 30
                        help='nominal distance from cathode to anode')
    parser.add_argument('-g', '--generator', type=str,
                        default='targetArray', #cosmic
                        help='the generator to use for building charge clouds')
    parser.add_argument('-d', '--drift', type=str,
                        default='randomWalk',
                        help='which drift model to use')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='set verbosity')

    args = parser.parse_args()

    outFile = args.output

    thisEfield = Efield()
    # thisEfield = Efield(args.transverse, args.longitudinal, args.pertFields)
    # thisEfield = Efield(args.transverse, args.longitudinal)

    Npe = args.N

    thisEventRecord = eventRecord()

    trackConditions = lambda thisTrack: is_inwards(thisTrack) & \
                                        is_anode_crosser(thisTrack) & \
                                        is_cathode_crosser(thisTrack)

    arrivalTimes = []
    finalXs = []
    finalYs = []
    finalZs = []

    initXs = []
    initYs = []
    initZs = []

    weights = []

    if args.generator == 'radSource':
        thisSource = radSource()
        thisTracklet = thisSource.generate_tracklet()
        nChargeBundles = int(thisTracklet.Qtot/sim_parameters["scalingF"])
        charges = [thisTracklet.generate_charge()
                   for i in range(nChargeBundles)]

    elif args.generator == 'cosmic':
        thisTrack = cosmicRayTrack(trackConditions).track

        thisEventRecord.pos = thisTrack.pos
        thisEventRecord.dir = thisTrack.dir
        thisEventRecord.length = thisTrack.length
        
        theseTracklets = thisTrack.generate_segments()
        charges = []
        for thisTracklet in theseTracklets:
            nChargeBundles = int(thisTracklet.Qtot/sim_parameters["scalingF"])
            for i in range(nChargeBundles):
                charges.append(thisTracklet.generate_charge())
        nChargeBundles = len(charges)
                
    elif args.generator == 'rock':
        thisTrack = rockMuonTrack(trackConditions).track

        thisEventRecord.pos = thisTrack.pos
        thisEventRecord.dir = thisTrack.dir
        thisEventRecord.length = thisTrack.length

        theseTracklets = thisTrack.generate_segments()
        charges = []
        for thisTracklet in theseTracklets:
            nChargeBundles = int(thisTracklet.Qtot/sim_parameters["scalingF"])
            for i in range(nChargeBundles):
                charges.append(thisTracklet.generate_charge())
        nChargeBundles = len(charges)

    elif args.generator == 'iso':
        trackGenerator = isotropicMuonTrack(trackConditions)

        thisTrack = trackGenerator.track

        thisEventRecord.pos = thisTrack.pos
        thisEventRecord.dir = thisTrack.dir
        thisEventRecord.length = thisTrack.length

        theseTracklets = thisTrack.generate_segments()
        charges = []
        for thisTracklet in theseTracklets:
            nChargeBundles = int(thisTracklet.Qtot/sim_parameters["scalingF"])
            for i in range(nChargeBundles):
                charges.append(thisTracklet.generate_charge())
        nChargeBundles = len(charges)
                
    elif args.generator == 'point':
        nChargeBundles = args.N
        charges = [charge([0, 0, detector_parameters["cathode position"]])
                   for i in range(nChargeBundles)]

    elif args.generator == 'disc':
        nChargeBundles = args.N
        charges = [charge(sample_from_cathode_target())
                   for i in range(nChargeBundles)]

    elif args.generator == 'targetArray':

        nTargetsPerRow = 8
        nTargetsPerCol = 8

        wallMargin = 5. # cm
        bounds = detector_parameters['detector bounds']
        
        # x is [-0.3225, 0.3225], y is [-0.625, 0.625], z is [-0.3, 0.3] in John's coords. 
        # We need to switch x and z for our simulation. 
        targetY = np.linspace(bounds[1][0] + wallMargin, bounds[1][1] - wallMargin, nTargetsPerRow)
        targetZ = np.linspace(bounds[2][0] + wallMargin, bounds[2][1] - wallMargin, nTargetsPerCol)
        # The x and z coordinates are switched between John's solution and our simulation. 
        # targetZ = np.linspace(bounds[2][0] + wallMargin, bounds[1][1] - wallMargin, nTargetsPerCol)

        targetLocs = [(y, z) for y in targetY for z in targetZ]
        
        nChargeBundles = args.N

        charges = [charge(sample_from_cathode_target(yC = y, zC = z))
                   for i in range(nChargeBundles)
                   for y, z in targetLocs]

        nChargeBundles *= len(targetLocs)

    else:
        raise ValueError("the specified generator is not a recognized option!")

    # for i in range(Npe):
    if args.verbose:
        print("drifting " + str(nChargeBundles) + " discrete charge bundles")
    for i, this_charge in enumerate(charges):

        # if args.verbose:
        #     print("charge " + str(i) + " originates at " + str(this_charge.pos))

        this_charge.drift(thisEfield, args.drift)

        x, y, z = np.array(this_charge.history).T

        if this_charge.fate == 1:
            arrivalTimes.append(this_charge.arrivalT)
            finalXs.append(this_charge.pos[0])
            finalYs.append(this_charge.pos[1])
            finalZs.append(this_charge.pos[2])

            initXs.append(this_charge.pos_i[0])
            initYs.append(this_charge.pos_i[1])
            initZs.append(this_charge.pos_i[2])

            weights.append(this_charge.weight)

            # if args.verbose:
            #     print("charge " + str(i) +
            #           " terminates at " + str(this_charge.pos))

    if args.verbose:
        print("writing record of " + str(len(arrivalTimes)) + " charges")

    n_charges = len(arrivalTimes)

    if n_charges > 0: #Check if the muon sliced through the volume 
        thisEventRecord.inside = True

    thisEventRecord.chargeMap = np.array([finalXs,
                                        finalYs,
                                        finalZs,
                                        arrivalTimes])

    thisEventRecord.QdepMap = np.array([initXs,
                                        initYs,
                                        initZs,
                                        arrivalTimes])

    thisEventRecord.weights = weights

    np.save(outFile, np.array([thisEventRecord]))


