import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import os
import sys

from parameters import *
from utils import *

class targetLocHypothesis:
    def __init__(self):
        self.set_center(0, 0)
        
    def set_center(self, x, y):
        self.center = np.array([x, y])
        
    def likelihood(self, x, y):
        probs = []
        for xi, yi in zip(x, y):
            if self.is_inside(xi, yi):
                probs.append(1.)
            else:
                probs.append(0.)
        return np.prod(probs)

    def frac_inside(self, x, y):
        inside = 0
        Ntot = float(len(x))
        for xi, yi in zip(x, y):
            if self.is_inside(xi, yi):
                inside += 1
        return float(inside)/Ntot
    
    def is_inside(self, x, y):
        dist_from_center = mag(np.array([x, y]) - self.center)
        return dist_from_center <= detector_parameters["target radius"]
            
if __name__ == '__main__':

    xresFig = plt.figure()
    xresAx = xresFig.gca()

    xaccFig = plt.figure()
    xaccAx = xaccFig.gca()

    tdistFig = plt.figure()
    tdistAx = tdistFig.gca()

    tResFig = plt.figure()
    tResAx = tResFig.gca()

    zResFig = plt.figure()
    zResAx = zResFig.gca()
    
    tAccFig = plt.figure()
    tAccAx = tAccFig.gca()

    zAccFig = plt.figure()
    zAccAx = zAccFig.gca()

    driftDir = 'driftHits'

    nomDriftV = physics_parameters['v'](0.5)
    print ("nominal drift velocity", nomDriftV)

    pixelPitch = 0.4434
    nPix = 16
    # nPix = 8
    xBins = np.linspace(-pixelPitch*(nPix - 0.5), pixelPitch*(nPix - 0.5), 2*nPix)
    yBins = np.linspace(-pixelPitch*(nPix - 0.5), pixelPitch*(nPix - 0.5), 2*nPix)
    xBinCenters = 0.5*(xBins[:-1] + xBins[1:])
    yBinCenters = 0.5*(yBins[:-1] + yBins[1:])
    X, Y = np.meshgrid(xBinCenters, yBinCenters)

    # tBins = np.linspace(0, 500, 251) # assume 2 us timing resolution, starting at flash time
    tBins = np.linspace(0, 500, 501) # assume 2 us timing resolution, starting at flash time
    tBinCenters = 0.5*(tBins[:-1] + tBins[1:])
    
    noiseLevel = 200
    hitThreshold = noiseLevel + 5*np.sqrt(noiseLevel) # hit threshold is a 5 sigma (assuming noise) count level  
    # hitThreshold = 0

    print (xBins)
    print (xBinCenters)

    nTrials = 1000
    
    # files = {500: {'fileList': ['driftHits/hits_500Vcm_'+str(i)+".npy"
    #                             for i in range(1, 9)],
    #                'label': '500 V/cm',
    #                'transv': 0},
    #          # 475: {'fileList': ['driftHits/hits_475Vcm_'+str(i)+".npy"
    #          #                    for i in range(1, 9)],
    #          #       'label': '475 V/cm'},
    #          # 525: {'fileList': ['driftHits/hits_525Vcm_'+str(i)+".npy"
    #          #                    for i in range(1, 9)],
    #          #       'label': '525 V/cm'},
    #          # 1000: {'fileList': ['driftHits/hits_1000Vcm_'+str(i)+".npy"
    #          #                     for i in range(1, 13)],
    #          #        'label': '1000 V/cm'},
    #          '500 + 0.5 transv': {'fileList': ['driftHits/hits_500Vcm_0-5Vcmtrans_'+str(i)+".npy"
    #                                          for i in range(1, 9)],
    #                             'label': '500 V/cm + 0.5 V/cm transverse',
    #                               'transv': 0.5},
    #          '500 + 1 transv': {'fileList': ['driftHits/hits_500Vcm_1Vcmtrans_'+str(i)+".npy"
    #                                           for i in range(1, 9)],
    #                             'label': '500 V/cm + 1 V/cm transverse',
    #                             'transv': 1},
    #          '500 + 1.5 transv': {'fileList': ['driftHits/hits_500Vcm_1-5Vcmtrans_'+str(i)+".npy"
    #                                            for i in range(1, 9)],
    #                               'label': '500 V/cm + 1.5 V/cm transverse',
    #                               'transv': 1.5},
    #          '500 + 2 transv': {'fileList': ['driftHits/hits_500Vcm_2Vcmtrans_'+str(i)+".npy"
    #                                           for i in range(1, 9)],
    #                             'label': '500 V/cm + 2 V/cm transverse',
    #                             'transv': 2},
    #          '500 + 2.5 transv': {'fileList': ['driftHits/hits_500Vcm_2-5Vcmtrans_'+str(i)+".npy"
    #                                           for i in range(1, 9)],
    #                             'label': '500 V/cm + 2.5 V/cm transverse',
    #                               'transv': 2.5},
    #          '500 + 3 transv': {'fileList': ['driftHits/hits_500Vcm_3Vcmtrans_'+str(i)+".npy"
    #                                           for i in range(1, 9)],
    #                             'label': '500 V/cm + 3 V/cm transverse',
    #                             'transv': 3},
    #          '500 + 3.5 transv': {'fileList': ['driftHits/hits_500Vcm_3-5Vcmtrans_'+str(i)+".npy"
    #                                            for i in range(1, 9)],
    #                               'label': '500 V/cm + 3.5 V/cm transverse',
    #                               'transv': 3.5},
    #          '500 + 4 transv': {'fileList': ['driftHits/hits_500Vcm_4Vcmtrans_'+str(i)+".npy"
    #                                          for i in range(1, 9)],
    #                             'label': '500 V/cm + 4 V/cm transverse',
    #                             'transv': 4},
    #          '500 + 4.5 transv': {'fileList': ['driftHits/hits_500Vcm_4-5Vcmtrans_'+str(i)+".npy"
    #                                            for i in range(1, 9)],
    #                               'label': '500 V/cm + 4.5 V/cm transverse',
    #                               'transv': 4.5},
    #          '500 + 1% transv': {'fileList': ['driftHits/hits_500Vcm_5Vcmtrans_'+str(i)+".npy"
    #                                           for i in range(1, 9)],
    #                              'label': '500 V/cm + 5 V/cm transverse',
    #                              'transv': 5},
    #          # '500 + 5% transv': {'fileList': ['driftHits/hits_500Vcm_25Vcmtrans_'+str(i)+".npy"
    #          #                                  for i in range(1, 9)],
    #          #                     'label': '500 V/cm + 25 V/cm transverse',
    #          #                     'transv': 25},
    # }
    files = {'475': {'fileList': ['driftHits/hits_475Vcm_'+str(i)+".npy"
                                  for i in range(1, 9)],
                     'label': '475 V/cm',
                     'transv': 0,
                     'longit': -25},
            '485': {'fileList': ['driftHits/hits_485Vcm_'+str(i)+".npy"
                                 for i in range(1, 5)],
                    'label': '485 V/cm',
                    'transv': 0,
                    'longit': -15},
            '490': {'fileList': ['driftHits/hits_490Vcm_'+str(i)+".npy"
                                 for i in range(1, 5)],
                    'label': '490 V/cm',
                    'transv': 0,
                    'longit': -10},
             '495': {'fileList': ['driftHits/hits_495Vcm_'+str(i)+".npy"
                                 for i in range(1, 5)],
                    'label': '495 V/cm',
                    'transv': 0,
                    'longit': -5},
             '500': {'fileList': ['driftHits/hits_500Vcm_'+str(i)+".npy"
                                  for i in range(1, 9)],
                     'label': '500 V/cm',
                     'transv': 0,
                   'longit': 0},
             '505': {'fileList': ['driftHits/hits_505Vcm_'+str(i)+".npy"
                                 for i in range(1, 5)],
                     'label': '505 V/cm',
                     'transv': 0,
                     'longit': 5},
             '510': {'fileList': ['driftHits/hits_510Vcm_'+str(i)+".npy"
                                  for i in range(1, 5)],
                     'label': '510 V/cm',
                     'transv': 0,
                     'longit': 10},
             '515': {'fileList': ['driftHits/hits_515Vcm_'+str(i)+".npy"
                                  for i in range(1, 5)],
                     'label': '515 V/cm',
                     'transv': 0,
                     'longit': 15},
            '525': {'fileList': ['driftHits/hits_525Vcm_'+str(i)+".npy"
                                 for i in range(1, 9)],
                    'label': '525 V/cm',
                    'transv': 0,
                    'longit': 25},

             # 1000: {'fileList': ['driftHits/hits_1000Vcm_'+str(i)+".npy"
             #                     for i in range(1, 13)],
             #        'label': '1000 V/cm',
             #        'transv': 0,
             #        'longit': 500},
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
        nSels = np.logspace(1, 4, 13, dtype = int)
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

                # fig = plt.figure()
                
                # plt.hist2d(xSel, ySel, bins = [xBins, yBins])
                # plt.xlabel(r'x (cm)')
                # plt.ylabel(r'y (cm)')
                # plt.colorbar()
                # plt.show()
                counts, xBinEdges, yBinEdges = np.histogram2d(xSel, ySel, bins = [xBins, yBins])

                noise = st.poisson.rvs(noiseLevel, size = counts.shape)
                # noise = np.zeros_like(counts)

                readout = counts + noise

                # readout = readout[readout > hitThreshold]
                readout = np.where(readout > hitThreshold, readout - noiseLevel*np.ones_like(readout), 0)

                measured_x_centers.append(np.sum(X*readout.T)/nSel)
                measured_y_centers.append(np.sum(Y*readout.T)/nSel)

                tCounts, tBinEdges = np.histogram(tSel, bins = tBins)

                tNoise = st.poisson.rvs(noiseLevel, size = tCounts.shape)
                
                tReadout = tCounts + tNoise
                tReadout = np.where(tReadout > hitThreshold, tReadout - noiseLevel*np.ones_like(tReadout), 0)
                tCenter = np.sum(tBinCenters*tReadout)/nSel
                
                measured_t_centers.append(tCenter)
                measured_z_centers.append(tCenter*nomDriftV)
                
            # plt.figure()

            # plt.hist2d(measured_x_centers,
            #            measured_y_centers)

            # print ("x measurement width", np.std(measured_x_centers))
            # plt.show()
            acc.append(np.mean(measured_x_centers))
            res.append(np.std(measured_x_centers))

            tAcc.append(np.mean(measured_t_centers))
            zAcc.append(np.mean(measured_z_centers))

            tRes.append(np.std(measured_t_centers))
            zRes.append(np.std(measured_z_centers))

            if nSel == 10000:
                if label in ['475 V/cm', '500 V/cm', '525 V/cm']:
                    tdistAx.hist(tObs, label = label, histtype = 'step', bins = tBins, density = True)

                transv.append(data['transv'])
                longit.append(data['longit'])
                recoLoc.append(np.mean(measured_x_centers))
                recoUnc.append(np.std(measured_x_centers))

                tLoc.append(np.mean(measured_t_centers))
                tUnc.append(np.std(measured_t_centers))

                zLoc.append(np.mean(measured_z_centers))
                zUnc.append(np.std(measured_z_centers))

                # plt.figure()
                # plt.imshow(counts.T)
                # plt.title('Signal')
                # plt.xlabel(r'x (cm)')
                # plt.ylabel(r'y (cm)')
                # plt.colorbar()

                # plt.figure()
                # plt.imshow(noise.T)
                # plt.title('Noise')
                # plt.xlabel(r'x (cm)')
                # plt.ylabel(r'y (cm)')
                # plt.colorbar()

                plt.figure()
                plt.imshow(readout.T)
                # plt.title('Signal + Noise')
                plt.xlabel(r'z (cm)')
                plt.ylabel(r'y (cm)')
                cb = plt.colorbar()
                cb.set_label(r'bkg. subtracted hits')
                plt.show()

                # plt.figure()

                # plt.title(label)
                # plt.hist2d(xSel,
                #            ySel,
                #            bins = [xBins, yBins])

                # print (label)
                # print ( X ) 
                # print ( counts.T ) 
                # print ( X*counts.T ) 
                # print ( np.sum(X*counts.T) ) 
                # print ( np.sum(X*counts.T)/nSel ) 
                # plt.show()


        print (label)
        print (np.mean(xObs), np.std(xObs))
        print (np.mean(yObs), np.std(yObs))
        print (np.mean(zObs), np.std(zObs))

        print(np.mean(measured_x_centers))
        print(np.std(measured_x_centers))
        
        xresAx.scatter(nSels, res, label = label)
        xaccAx.scatter(nSels, np.abs(acc), label = label)

        tResAx.scatter(nSels, tRes, label = label)
        tAccAx.scatter(nSels, tAcc, label = label)

        zResAx.scatter(nSels, zRes, label = label)
        zAccAx.scatter(nSels, zAcc, label = label)

        # xaccAx.errorbar(nSels,
        #                 np.abs(acc),
        #                 yerr = res,
        #                 label = label)

    # x resolution (n_elec vs sigma)
    xresAx.legend()
    xresAx.loglog()
    xresAx.set_xlabel(r'$N_{\mathrm{electrons}}$')
    xresAx.set_ylabel(r'$\sigma(x_{\mathrm{reco}})$')

    # x accuracy (n_elec vs mean position) 
    xaccAx.legend()
    xaccAx.semilogx()
    # xaccAx.loglog()
    xaccAx.set_xlabel(r'$N_{\mathrm{electrons}}$')
    xaccAx.set_ylabel(r'$| \; \mu(x_{\mathrm{reco}}) \; |$')

    # arrival time distribution  
    tdistAx.legend()
    tdistAx.set_xlabel(r'Arrival time ($\mu$s)')
    # tdistAx.set_ylabel(r'')

    # t resolution (n_elec vs timing sigma)  
    tResAx.legend()
    tResAx.loglog()
    tResAx.set_xlabel(r'$N_{\mathrm{electrons}}$')
    tResAx.set_ylabel(r'$\sigma(t)$')

    # t accuracy (n_elec vs timing mean)  
    tAccAx.legend()
    tAccAx.semilogx()
    # tAccAx.loglog()
    tAccAx.set_xlabel(r'$N_{\mathrm{electrons}}$')
    tAccAx.set_ylabel(r'$\mu(t)$')

    # z resolution (n_elec vs. sigma(z))
    zResAx.legend()
    zResAx.loglog()
    zResAx.set_xlabel(r'$N_{\mathrm{electrons}}$')
    zResAx.set_ylabel(r'$\sigma(z_{\mathrm{reco}})$')

    # z accuracy (n_elec vs. z mean)
    zAccAx.legend()
    zAccAx.semilogx()
    # tAccAx.loglog()
    zAccAx.set_xlabel(r'$N_{\mathrm{electrons}}$')
    zAccAx.set_ylabel(r'$| \; \mu(z_{\mathrm{reco}}) \; |$')

    
    fig = plt.figure()
    ax = fig.gca()
    # ax.scatter(transv, np.abs(recoLoc),
    #            # yerr = recoUnc,
    #            # ls = 'none',
    # )
    # ax.errorbar(transv, np.abs(recoLoc),
    #             yerr = recoUnc,
    #             ls = 'none',
    #             marker = '+',
    # )
    ax.errorbar(longit, np.abs(recoLoc),
                yerr = recoUnc,
                ls = 'none',
                marker = '+',
    )
    # ax.set_xlabel(r'Transverse E field (V/cm)')
    ax.set_xlabel(r'Longitudinal E field excess (V/cm)')
    ax.set_ylabel(r'Reconstructed Position (cm)')

    fig = plt.figure()
    ax = fig.gca()
    ax.errorbar(longit, np.abs(tLoc),
                yerr = tUnc,
                ls = 'none',
                marker = '+',
    )
    ax.set_xlabel(r'Longitudinal E field excess (V/cm)')
    ax.set_ylabel(r'Arrival time ($\mu$s)')
    
    fig = plt.figure()
    ax = fig.gca()
    ax.errorbar(longit, np.abs(zLoc),
                yerr = zUnc,
                ls = 'none',
                marker = '+',
    )
    ax.set_xlabel(r'Longitudinal E field excess (V/cm)')
    ax.set_ylabel(r'Reconstructed Cathode distance (cm)')

    plt.show() 
