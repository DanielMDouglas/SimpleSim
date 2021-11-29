import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import os
import sys

from parameters import *
from utils import *

# Beam z, Zenith y, x drift
def form_hits(finalLocs):
    """
    Return array of hits (x center, y center, t) 
    """
    x, y, z, t = finalLocs
    sample = np.array([z, y, t]).T

    thresh = detector_parameters["pixel threshold"] # e's
    pixelPitch = 0.4434
    nPix = 34

    zBins = np.linspace(-pixelPitch*(nPix - 0.5), pixelPitch*(nPix - 0.5), 2*nPix)
    yBins = np.linspace(-pixelPitch*(nPix - 0.5), pixelPitch*(nPix - 0.5), 2*nPix)
    zBinCenters = 0.5*(zBins[:-1] + zBins[1:])
    yBinCenters = 0.5*(yBins[:-1] + yBins[1:])
    # X, Y = np.meshgrid(xBinCenters, yBinCenters)

    t_start, t_end = int( round(min(t)) )-1, int( round(max(t)) ) + 1 #us
    t_range_us = ( (t_end - t_start) )
    n_tbins_200ns = t_range_us * 5
    tBins = np.linspace( t_start, t_end, n_tbins_200ns )
    # tBins = np.linspace(290, 330, 41)
    tBinCenters = 0.5*(tBins[:-1] + tBins[1:])

    counts, edges = np.histogramdd(sample, bins = [zBins, yBins, tBins])

    H = []
    for i, zCenter in enumerate(zBinCenters):
        for j, yCenter in enumerate(yBinCenters):
            tHist = counts[i, j, :] # At some pixel, get the counts at each time tick (us)
            counts_accum = 0
            counts_registered = 0
            collection_clock_ticks = 0
            for k, count in enumerate(tHist): 
                counts_accum += count
                if counts_accum > thresh:
                    collection_clock_ticks += 1
                if collection_clock_ticks == 8:
                    counts_collected = counts_accum 
                    H.append( [zCenter, yCenter, t_start + 0.2*k, 1e2*counts_collected] )
                if collection_clock_ticks == 11: 
                    counts_accum = 0
                    collection_clock_ticks = 0
                    
    return np.array( H )



        
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type = str,
                        default = 'driftHits.npy',
                        help = 'where read the destinations from')
    parser.add_argument('-o', '--output', type=str,
                        default='pixels.npy',
                        help='where to save')
    args = parser.parse_args()
    outFile = args.output
    finalLocs = np.load(args.input)

    H = form_hits(finalLocs) # Pixels are in the y,z plane
    # Beam z, Zenith y, x drift

    np.save(outFile, H)