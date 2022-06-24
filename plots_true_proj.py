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
        # x, y, z, t = thisRecord.chargeMap
        v = 0.1544 #cm/us for E = 0.5 kV/cm

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(1, 1, 1)
        ax.set_aspect('equal')
        # plt.scatter(z, y, c = v*t, s = 1, cmap='viridis')

        x, y, z, t = thisRecord.QdepMap
        plt.scatter(z, y, c='gray', s = 2, marker = 'o', alpha = 0.1) 
        x, y, z, t = thisRecord.chargeMap
        plt.scatter(z, y, c='blue', s = 2, marker = 'o', alpha = 0.1) 

        # cbar= plt.colorbar()
        # cbar.set_label("X [cm]",fontsize = 20)
        # cbar.ax.tick_params(labelsize=20)
        plt.xlabel("X [cm]",fontsize = 20)
        plt.ylabel("Y [cm]",fontsize = 20)
        plt.xlim([-30, 30])
        plt.ylim([-60, 60])
        plt.grid(True)
        plt.tick_params(axis='both', which='both', labelsize = 20, direction = 'in')

        xticks = np.arange(-30, 40, 10)
        yticks = np.arange(-60, 70, 10)

        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
  
        plt.show()