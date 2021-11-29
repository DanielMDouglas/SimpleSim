import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from utils import *

# Pretty fonts for figures 
mpl.rc('text', usetex = True)
mpl.rc('font', family='SignPainter')

pixelPitch = 0.4434
nPix = 34

zBins = np.linspace(-pixelPitch*(nPix - 0.5), pixelPitch*(nPix - 0.5), 2*nPix)
yBins = np.linspace(-pixelPitch*(nPix - 0.5), pixelPitch*(nPix - 0.5), 2*nPix)
zBinCenters = 0.5*(zBins[:-1] + zBins[1:])
yBinCenters = 0.5*(yBins[:-1] + yBins[1:])

if __name__ == '__main__':
        import argparse
        parser = argparse.ArgumentParser()

        parser.add_argument('input', nargs='+',
                                type=str,
                                help='data to show')

        args = parser.parse_args()
        pixels = np.load(args.input[0])

        z, y, t, q = pixels.T
        v = 0.1544 #cm/us for E = 0.5 kV/cm

        # Format expected by PCA function 
        d = np.array([z,y,v*t]).T
        zs,ys = get_pca(d) #Get PCA (two points in z,y defined by principal axis)

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_aspect('equal')
        plt.scatter(z, y, c = v*t, s = 10, cmap='viridis', marker='s')

        #Plot PCA with dotted black line 
        plt.plot(zs,ys, c = 'black', linestyle = 'dashed', linewidth = 2)

        cbar= plt.colorbar()
        cbar.set_label("X [cm]",fontsize = 20)
        cbar.ax.tick_params(labelsize=20)
        plt.xlabel("Z [cm]",fontsize = 20)
        plt.ylabel("Y [cm]",fontsize = 20)
        plt.xlim([-15, 15])
        plt.ylim([-15, 15])
        plt.grid(True)
        plt.tick_params(axis='both', which='both', labelsize = 20, direction = 'in')

        # Major ticks every 20, minor ticks every 5
        major_ticks = np.arange(-15, 15, 5)
        # Set the minor tick marks to match the pixel binning
        minor_ticks = zBins # same as yBins

        ax.set_xticks(major_ticks)
        ax.set_xticks(minor_ticks, minor=True)
        ax.set_yticks(major_ticks)
        ax.set_yticks(minor_ticks, minor=True)

        # And a corresponding grid
        ax.grid(which='both')

        # Or if you want different settings for the grids:
        ax.grid(which='minor', alpha=0.2)
        ax.grid(which='major', alpha=0.5)

        # ax.set_aspect('equal')

        plt.show()


        
