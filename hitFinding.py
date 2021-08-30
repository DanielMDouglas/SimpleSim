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

tBins = np.linspace(290, 330, 41)
tBinCenters = 0.5*(tBins[:-1] + tBins[1:])

class pixel:
    def __init__(self, (xCenter, yCenter), arrivalTimes):
        boundaries = []
        self.xCenter = xCenter
        self.yCenter = yCenter
        # arrival times is a histogram (tBins, nChargesArrived)
        self.arrivalTimes = arrivalTimes
        self.hits = []
    def find_hits(self):
        thisT = 0
        thisQ = 0
        scaledThreshold = detector_parameters["pixel threshold"]/sim_parameters["scalingF"]

        if np.any(self.arrivalTimes):
            print (self.arrivalTimes)
            plt.step(tBinCenters, self.arrivalTimes)

            accumulatedCharge = np.cumsum(self.arrivalTimes)

            print (accumulatedCharge)
            
            # hitTime
            # hitTime = np.interp(scaledThreshold, accumulatedCharge, tBinCenters)

            print (hitTime)

            plt.step(tBinCenters, accumulatedCharge)

            digitizedCharge = accumulatedCharge(hitTime)
            # plt.axvline(x = hitTime)
            # plt.axhline(y = scaledThreshold)
            
            plt.show()

        # if it crosses the threshold, when (at what time tick does the qAccum > threshold)?
        # go + 2us
        # get the value of qAccum
            
        # integrate arriving charge
        # when integral hits threshold, seek + 2 us
        # yield a "hit" with clock timestamp + integrated charge up to then
        self.hits.append(hit(thisT, thisQ))
        # subtract that charge
        # repeat the integration from there until the end of the window
        return self.hits
        
        
class hit:
    def __init__(self, t, q):
        self.t = t
        self.q = q
        
class hitMap:
    def __init__(self, hits):
        self.hits = hits
        
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type = str,
                        default = 'driftHits.npy',
                        help = 'where read the destinations from')
    args = parser.parse_args()

    finalLocs = np.load(args.input)

    x, y, z, t = finalLocs


    # TODO: find LArPix clock rate!
    # for now, we are assuming 100 MHz = 10 ns
    
    print (tBins)
    # counts, xBinEdges, yBinEdges = np.histogram2d(x, y, bins = [xBins, yBins])
    # bin in x, y, t
    sample = np.array([x, y, t]).T
    print (sample.shape)
    # counts, xBinEdges, yBinEdges, tBinEdges = np.histogramdd(sample, bins = [xBins, yBins, tBins])
    counts, edges = np.histogramdd(sample, bins = [xBins, yBins, tBins])
    print (counts.shape)

    pixels = []
    for i, xCenter in enumerate(xBinCenters):
        for j, yCenter in enumerate(yBinCenters):
            tHist = counts[i, j, :]
            # create a new pixel with the right center (according to the bin centers)
            thisPixel = pixel((xCenter, yCenter), tHist)
            # assign that pixel the corresponding t histogram (column of counts)

            thisPixel.find_hits() # do integration from t = 0 right and finds thresho...

            pixels.append(thisPixel)

    
    # print (np.digitize(np.array([x, y]), bins = np.array([xBins, yBins])))
    # TODO: find a way to assign arrival times to pixels based on arrival position

    
