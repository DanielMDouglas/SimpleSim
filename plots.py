import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from eventRecord import *
from parameters import *

# Pretty fonts for figures 
mpl.rc('text', usetex = True)
mpl.rc('font', family='SignPainter')


def draw_boundaries(ax):
#     from parameters import detector_parameters
#     bounds = [[-15, 15], [0, 30], [-15, 15]]
    bounds = detector_parameters['detector bounds']

    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][0], bounds[1][0]],
            [bounds[2][0], bounds[2][0]], color='black', ls='--')
    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][1], bounds[1][1]],
            [bounds[2][0], bounds[2][0]], color='black', ls='--')
    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][0], bounds[1][0]],
            [bounds[2][1], bounds[2][1]], color='black', ls='--')
    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][1], bounds[1][1]],
            [bounds[2][1], bounds[2][1]], color='black', ls='--')

    ax.plot([bounds[0][0], bounds[0][0]],
            [bounds[1][0], bounds[1][1]],
            [bounds[2][0], bounds[2][0]], color='black', ls='--')
    ax.plot([bounds[0][1], bounds[0][1]],
            [bounds[1][0], bounds[1][1]],
            [bounds[2][0], bounds[2][0]], color='black', ls='--')
    ax.plot([bounds[0][0], bounds[0][0]],
            [bounds[1][0], bounds[1][1]],
            [bounds[2][1], bounds[2][1]], color='black', ls='--')
    ax.plot([bounds[0][1], bounds[0][1]],
            [bounds[1][0], bounds[1][1]],
            [bounds[2][1], bounds[2][1]], color='black', ls='--')

    ax.plot([bounds[0][0], bounds[0][0]],
            [bounds[1][0], bounds[1][0]],
            [bounds[2][0], bounds[2][1]], color='black', ls='--')
    ax.plot([bounds[0][1], bounds[0][1]],
            [bounds[1][0], bounds[1][0]],
            [bounds[2][0], bounds[2][1]], color='black', ls='--')
    ax.plot([bounds[0][0], bounds[0][0]],
            [bounds[1][1], bounds[1][1]],
            [bounds[2][0], bounds[2][1]], color='black', ls='--')
    ax.plot([bounds[0][1], bounds[0][1]],
            [bounds[1][1], bounds[1][1]],
            [bounds[2][0], bounds[2][1]], color='black', ls='--')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('input', nargs='+',
                        type=str,
                        help='data to show')

    args = parser.parse_args()
    thisRecord = np.load(args.input[0], allow_pickle=True)[0]

    x, y, z, t = thisRecord.chargeMap
    # x, y, z, t = thisRecord.QdepMap

    fig = plt.figure()
    # ax = fig.add_subplot(projection = '3d')

    ax = Axes3D(fig)
    ax.scatter(x, y, z, c=t, s = 1)
    #     ax.scatter(x, y, z, c = t)

    # also plot the deposited charge positions
    x, y, z, t = thisRecord.QdepMap
    ax.scatter(x, y, z, c='gray', s = 1, marker = '+')

    
    pos = thisRecord.pos
    dir = thisRecord.dir
    length = thisRecord.length

    p0 = pos
    p1 = pos + length*dir

    plt.plot([p0[0], p1[0]], [p0[1], p1[1]], [p0[2], p1[2]], ls = '--', c = 'r')

    # target_radius = 0.2
    # tspace = np.linspace(0, 2*np.pi, 1000)
    # ax.plot(target_radius*np.cos(tspace), target_radius*np.sin(tspace), 50, color = 'black', ls = '--')

    draw_boundaries(ax)

    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')

    ax.set_xlabel(r'Drift Direction [cm]',fontsize = 15)
    ax.set_ylabel(r'Zenith Direciton [cm]',fontsize = 15)
    ax.set_zlabel(r'Beam Direction [cm]',fontsize = 15)
    plt.tick_params(axis='both', which='both', labelsize = 15, direction = 'in')

    bounds = detector_parameters['detector bounds']
    margin = 2.
    ax.set_xlim(bounds[0][0] - margin, bounds[0][1] + margin)
    ax.set_ylim(bounds[1][0] - margin, bounds[1][1] + margin)
    ax.set_zlim(bounds[2][0] - margin, bounds[2][1] + margin)

    plt.show()

    # pixelPitch = 0.4434
    # nPix = 3
    # xBins = np.linspace(-pixelPitch*(nPix - 0.5), pixelPitch*(nPix - 0.5), 2*nPix)
    # yBins = np.linspace(-pixelPitch*(nPix - 0.5), pixelPitch*(nPix - 0.5), 2*nPix)
    # xBinCenters = 0.5*(xBins[:-1] + xBins[1:])
    # yBinCenters = 0.5*(yBins[:-1] + yBins[1:])
    # X, Y = np.meshgrid(xBinCenters, yBinCenters)

    # tBins = np.linspace(290, 330, 41)
    # tBinCenters = 0.5*(tBins[:-1] + tBins[1:])

    # # counts, edges = np.histogramdd(, bins = [xBins, yBins, tBins])

    # # print (counts.shape)

    # import matplotlib as mpl

    # plt.hist2d(x, y, bins = [xBins, yBins], norm = mpl.colors.LogNorm())
    # cb = plt.colorbar(label = r'raw charge per pad [e]')
    # plt.xlabel(r'x [cm]')
    # plt.ylabel(r'y [cm]')
    # plt.show()
