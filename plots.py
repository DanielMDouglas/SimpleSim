import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Pretty fonts for figures 
mpl.rc('text', usetex = True)
mpl.rc('font', family='SignPainter')


def draw_boundaries(ax):
#     from parameters import detector_parameters
    bounds = [[-15, 15], [0, 30], [-15, 15]]

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
    data = np.load(args.input[0])

#     print(data.shape)

    x, y, z, t = data

    fig = plt.figure()
    # ax = fig.add_subplot(projection = '3d')

    ax = Axes3D(fig)
    ax.scatter(z, x, y, c=t, s = 1)
    #     ax.scatter(x, y, z, c = t)

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

    ax.set_xlabel(r'Beam Direction [cm]',fontsize = 15)
    ax.set_ylabel(r'Drift Direction [cm]',fontsize = 15)
    ax.set_zlabel(r'Zenith Direciton [cm]',fontsize = 15)
    plt.tick_params(axis='both', which='both', labelsize = 15, direction = 'in')

#     ax.set_xlabel(r'x [cm]')
#     ax.set_ylabel(r'y [cm]')
#     ax.set_zlabel(r'z [cm]')

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
