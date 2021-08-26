import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import os
import sys

from parameters import *
from utils import *

pixelPitch = 0.4434
nPix = 16

xBins = np.linspace(-pixelPitch*(nPix - 0.5), pixelPitch*(nPix - 0.5), 2*nPix)
yBins = np.linspace(-pixelPitch*(nPix - 0.5), pixelPitch*(nPix - 0.5), 2*nPix)
xBinCenters = 0.5*(xBins[:-1] + xBins[1:])
yBinCenters = 0.5*(yBins[:-1] + yBins[1:])
X, Y = np.meshgrid(xBinCenters, yBinCenters)

tBins = 0.2*np.arange(2000)

class pixel:
    def __init__(self, arrivalTimes):
        bins, self.binnedT = np.histogram(arrivalTimes, bins = tBins)
        # self.binnedT *= sim_parameters["scalingF"]
        self.arrivalTimes = arrivalTimes
    def find_hits(self):
        accumulatedCharge = np.cumsum(self.arrivalTimes)
        hitTime = np.interpolate(detector 

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type = str,
                        default = 'driftHits.npy',
                        help = 'where read the destinations from')
    args = parser.parse_args()

    finalLocs = np.load(args.input)

    x, y, z, t = finalLocs

    counts, xBinEdges, yBinEdges = np.histogram2d(x, y, bins = [xBins, yBins])

    # plt.hist2d(x, y)
    # plt.show()

    # print (np.digitize(np.array([x, y]), bins = np.array([xBins, yBins])))
    # TODO: find a way to assign arrival times to pixels based on arrival position

    
