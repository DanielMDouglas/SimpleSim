import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import os
import sys

from parameters import *
from utils import *

# set up binning
# right now, bins are set up such that x = y = 0 (target center)
# is aligned with the center of a pixel
pixelPitch = 0.4434
nPix = 16

xBins = np.linspace(-pixelPitch*(nPix - 0.5), pixelPitch*(nPix - 0.5), 2*nPix)
yBins = np.linspace(-pixelPitch*(nPix - 0.5), pixelPitch*(nPix - 0.5), 2*nPix)
xBinCenters = 0.5*(xBins[:-1] + xBins[1:])
yBinCenters = 0.5*(yBins[:-1] + yBins[1:])
X, Y = np.meshgrid(xBinCenters, yBinCenters)

tBins = np.linspace(0, 500, 501) # assume 2 us timing resolution, starting at flash time
tBinCenters = 0.5*(tBins[:-1] + tBins[1:])

nomDriftV = physics_parameters['v'](detector_parameters["nominal field"])
noiseLevel = detector_parameters["noise level"]
hitThreshold = noiseLevel + 5*np.sqrt(noiseLevel) # hit threshold is a 5 sigma (assuming noise) count level  

class targetLocHypothesis:
    def __init__(self, xLocs, yLocs, zLocs, tLocs):
        self.reconstruct_flash(xLocs, yLocs, zLocs, tLocs)

    def reconstruct_flash(self, xLocs, yLocs, zLocs, tLocs):
        nSel = len(xLocs)
    
        counts, xBinEdges, yBinEdges = np.histogram2d(xLocs, yLocs, bins = [xBins, yBins])
    
        noise = st.poisson.rvs(noiseLevel, size = counts.shape)

        readout = counts + noise
        
        self.xCenter = np.sum(X*readout.T)/nSel
        self.yCenter = np.sum(Y*readout.T)/nSel

        self.xSpread = np.sqrt(np.sum(np.power(X, 2)*readout.T)/nSel - np.power(self.xCenter, 2))
        self.ySpread = np.sqrt(np.sum(np.power(Y, 2)*readout.T)/nSel - np.power(self.yCenter, 2))
        
        tCounts, tBinEdges = np.histogram(tSel, bins = tBins)
    
        tNoise = st.poisson.rvs(noiseLevel, size = tCounts.shape)

        tReadout = tCounts + tNoise
        tReadout = np.where(tReadout > hitThreshold, tReadout - noiseLevel*np.ones_like(tReadout), 0)

        self.tCenter = np.sum(tBinCenters*tReadout)/nSel
        self.tSpread = np.sqrt(np.sum(np.power(tBinCenters, 2)*tReadout)/nSel - np.power(self.tCenter, 2))
            
            
if __name__ == '__main__':

    driftDir = 'driftHits'

    print ("nominal drift velocity", nomDriftV)

    nTrials = 1000
    
    # specify some files for reconstructing from
    # 
    files = {'500': {'fileList': ['driftHits/hits_500Vcm.npy'],
                     'label': '500 V/cm',
                     'transv': 0,
                     'longit': 0},
    }

    transv = []
    longit = []
    recoLoc = []
    recoUnc = []

    tLoc = []
    tUnc = []

    zLoc = []
    zUnc = []
    
    for thisE, data in files.items():
        xObs = np.array([])
        yObs = np.array([])
        zObs = np.array([])
        tObs = np.array([])
        
        fileList = data['fileList']
        label = data['label']
        for inFile in fileList:
            thisX, thisY, thisZ, thisT = np.load(inFile)
            xObs = np.concatenate((xObs, thisX))
            yObs = np.concatenate((yObs, thisY))
            zObs = np.concatenate((zObs, thisZ))
            tObs = np.concatenate((tObs, thisT))

        res = []
        acc = []

        tRes = []
        tAcc = []

        zRes = []
        zAcc = []
        
        # nSels = [50, 100, 500, 1000, 5000, 10000]
        nSels = np.logspace(1, 3, 4, dtype = int)
        for nSel in nSels:
    
            measured_x_centers = []
            measured_y_centers = []
            measured_t_centers = []
            measured_z_centers = []

            for i in range(nTrials):
                selected_indices = np.random.choice(xObs.size, nSel, replace  = False)
                xSel = xObs[selected_indices]
                ySel = yObs[selected_indices]
                zSel = zObs[selected_indices]
                tSel = tObs[selected_indices]

                reco = targetLocHypothesis(xSel, ySel, zSel, tSel)

                measured_x_centers.append(reco.xCenter)
                measured_y_centers.append(reco.yCenter)
                measured_t_centers.append(reco.tCenter)
                measured_z_centers.append(reco.tCenter*nomDriftV)
                
            acc.append(np.mean(measured_x_centers))
            res.append(np.std(measured_x_centers))

            tAcc.append(np.mean(measured_t_centers))
            zAcc.append(np.mean(measured_z_centers))

            tRes.append(np.std(measured_t_centers))
            zRes.append(np.std(measured_z_centers))

            if nSel == 10000:
                transv.append(data['transv'])
                longit.append(data['longit'])
                recoLoc.append(np.mean(measured_x_centers))
                recoUnc.append(np.std(measured_x_centers))

                tLoc.append(np.mean(measured_t_centers))
                tUnc.append(np.std(measured_t_centers))

                zLoc.append(np.mean(measured_z_centers))
                zUnc.append(np.std(measured_z_centers))


        print (label)
        print (np.mean(xObs), np.std(xObs))
        print (np.mean(yObs), np.std(yObs))
        print (np.mean(zObs), np.std(zObs))

        print(np.mean(measured_x_centers))
        print(np.std(measured_x_centers))
