import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Pretty fonts for figures 
mpl.rc('text', usetex = True)
mpl.rc('font', family='SignPainter')

if __name__ == '__main__':
        import argparse
        parser = argparse.ArgumentParser()

        parser.add_argument('input', nargs='+',
                                type=str,
                                help='data to show')

        args = parser.parse_args()
        thisRecord = np.load(args.input[0], allow_pickle=True)[0]
        x, y, z, t = thisRecord.chargeMap
        v = 0.1544 #cm/us for E = 0.5 kV/cm

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_aspect('equal')
        plt.scatter(z, y, c = v*t, s = 1, cmap='viridis')
        cbar= plt.colorbar()
        cbar.set_label("X [cm]",fontsize = 20)
        cbar.ax.tick_params(labelsize=20)
        plt.xlabel("Z [cm]",fontsize = 20)
        plt.ylabel("Y [cm]",fontsize = 20)
        plt.xlim([-15, 15])
        plt.ylim([-15, 15])
        plt.grid(True)
        plt.tick_params(axis='both', which='both', labelsize = 20, direction = 'in')
        plt.show()


        
