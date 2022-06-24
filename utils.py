import numpy as np
import scipy.stats as st
from numpy import linalg as LA
from array import array
import ROOT
from ROOT import TFile, TH2F, TH3F, TCanvas, TLegend, gStyle, TGaxis
import math

def mag(vect):
    return np.sqrt(np.sum(np.power(vect, 2)))


def norm(vect):
    return vect/mag(vect)


def get_perpendicular_vectors(vect):
    # TODO
    # perp1 = np.array([1, 0, 0])
    # perp2 = np.array([0, 1, 0])
    perp1 = np.random.randn(3)
    perp1 -= perp1.dot(vect)*vect
    perp1 /= np.linalg.norm(perp1)

    perp2 = np.cross(vect, perp1)

    return perp1, perp2


def is_in_tpc(pos):

    from parameters import detector_parameters

    bounds = detector_parameters['detector bounds']

    is_in_x_bounds = (pos[0] >= bounds[0][0]) & (pos[0] <= bounds[0][1])
    is_in_y_bounds = (pos[1] >= bounds[1][0]) & (pos[1] <= bounds[1][1])
    is_in_z_bounds = (pos[2] >= bounds[2][0]) & (pos[2] <= bounds[2][1])

    return is_in_x_bounds and is_in_y_bounds and is_in_z_bounds

def is_inwards(thisTrack):
    """
    Check to see that the track is pointing towards the 
    detector center
    """
    from parameters import detector_parameters

    return np.dot(norm(thisTrack.pos - detector_parameters['detector center']), thisTrack.dir) < 0

def is_anode_crosser(thisTrack):
    """
    Check to see if a given track (in slope-intercept form)
    will cross the anode detector boundary defined in `parameters.py`
    (this boundary is the one on the low side of the x-axis)
    """

    pos, direc = thisTrack.pos, thisTrack.dir

    from parameters import detector_parameters
    
    bounds = detector_parameters['detector bounds']

    # find the extent along the track where the track crosses the anod plane
    extent = (bounds[0][0] - pos[0])/direc[0]

    # and the coordinates in y, z at that point
    yInt = (pos + extent*direc)[1]
    zInt = (pos + extent*direc)[2]

    # are those points on the detector face?
    if yInt > bounds[1][0] and yInt < bounds[1][1]:
        if zInt > bounds[2][0] and zInt < bounds[2][1]:
            return True

    return False

def is_cathode_crosser(thisTrack):
    """
    Check to see if a given track (in slope-intercept form)
    will cross the cathode detector boundary defined in `parameters.py`
    (this boundary is the one on the high side of the x-axis)
    """

    pos, direc = thisTrack.pos, thisTrack.dir
    
    from parameters import detector_parameters

    bounds = detector_parameters['detector bounds']

    # find the extent along the track where the track crosses the anod plane
    extent = (bounds[0][1] - pos[0])/direc[0]

    # and the coordinates in y, z at that point
    yInt = (pos + extent*direc)[1]
    zInt = (pos + extent*direc)[2]
    
    # are those points on the detector face?
    if yInt > bounds[1][0] and yInt < bounds[1][1]:
        if zInt > bounds[2][0] and zInt < bounds[2][1]:
            return True

    return False
