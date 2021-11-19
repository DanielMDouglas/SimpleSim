import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def draw_boundaries(ax):
    from parameters import detector_parameters
    bounds = detector_parameters['detector bounds']
    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][0], bounds[1][0]],
            [bounds[2][0], bounds[2][0]] , color = 'black', ls = '--')
    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][1], bounds[1][1]],
            [bounds[2][0], bounds[2][0]] , color = 'black', ls = '--')
    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][0], bounds[1][0]],
            [bounds[2][1], bounds[2][1]] , color = 'black', ls = '--')
    ax.plot([bounds[0][0], bounds[0][1]],
            [bounds[1][1], bounds[1][1]],
            [bounds[2][1], bounds[2][1]] , color = 'black', ls = '--')

    ax.plot([bounds[0][0], bounds[0][0]],
            [bounds[1][0], bounds[1][1]],
            [bounds[2][0], bounds[2][0]] , color = 'black', ls = '--')
    ax.plot([bounds[0][1], bounds[0][1]],
            [bounds[1][0], bounds[1][1]],
            [bounds[2][0], bounds[2][0]] , color = 'black', ls = '--')
    ax.plot([bounds[0][0], bounds[0][0]],
            [bounds[1][0], bounds[1][1]],
            [bounds[2][1], bounds[2][1]] , color = 'black', ls = '--')
    ax.plot([bounds[0][1], bounds[0][1]],
            [bounds[1][0], bounds[1][1]],
            [bounds[2][1], bounds[2][1]] , color = 'black', ls = '--')

    ax.plot([bounds[0][0], bounds[0][0]],
            [bounds[1][0], bounds[1][0]],
            [bounds[2][0], bounds[2][1]] , color = 'black', ls = '--')
    ax.plot([bounds[0][1], bounds[0][1]],
            [bounds[1][0], bounds[1][0]],
            [bounds[2][0], bounds[2][1]] , color = 'black', ls = '--')
    ax.plot([bounds[0][0], bounds[0][0]],
            [bounds[1][1], bounds[1][1]],
            [bounds[2][0], bounds[2][1]] , color = 'black', ls = '--')
    ax.plot([bounds[0][1], bounds[0][1]],
            [bounds[1][1], bounds[1][1]],
            [bounds[2][0], bounds[2][1]] , color = 'black', ls = '--')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument('input', nargs = '+',
                        type = str,
                        help = 'data to show')

    args = parser.parse_args()
    data = np.load(args.input[0])

    print (data.shape)

    x, y, z, t = data

    fig = plt.figure()
    # ax = fig.add_subplot(projection = '3d')
    ax = Axes3D(fig)

    ax.scatter(x, y, z, c = t)

    # target_radius = 0.2
    # tspace = np.linspace(0, 2*np.pi, 1000)
    # ax.plot(target_radius*np.cos(tspace), target_radius*np.sin(tspace), 50, color = 'black', ls = '--')

    draw_boundaries(ax)
    
    ax.set_xlabel(r'x [cm]')
    ax.set_ylabel(r'y [cm]')
    ax.set_zlabel(r'z [cm]')
    
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

   # Testing testing




