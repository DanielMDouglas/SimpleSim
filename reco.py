import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
import os
import sys

from parameters import *
v = physics_parameters["v_nominal"] #drift velocity 

# Will not use this for Cathode-Anode crossers. There, just get direction with two points. 
def get_pca( X ): #Expects array of z, y, x

    j = 0
    w_X_list = []
    w = [1]*len(X) #Equal weights
    while j < len(X):
        w_X_list.append( 1.0*X[j] )
        j += 1
    mu = np.sum( w_X_list, axis = 0 )/np.sum(w)

    X = X - mu
    W = np.diag(w)
    C = (X.T).dot(W).dot(X)/np.sum(w) #Weighted covariance matrix
    eigenvalues, eigenvectors = LA.eig(C)
    eigenvectors = -eigenvectors.T # The columns (not the rows) of LA.eig(C) are the eigenvectors. 

    # Sort the e-vals and e-vecs by the e-val (from max to min).
    a = eigenvalues.argsort()[::-1] # Returns the indices that would sort the array (max to min).
    eigenvalues = eigenvalues[a]
    eigenvectors = eigenvectors[a] # Column (e-vec) i corresponds to the ith e-val.    
    
    vec = eigenvectors[0] #Primary Axis
    L = np.sqrt(eigenvalues[0]) #Using this as a proxy for mid-length (no good reason)
    v = vec*2*L  #vector from midpoint in direction of PCA
    p1, p2 = mu - v, mu + v
    
    # Two points along PCA 
    z1,y1,x1 = p1
    z2,y2,x2 = p2
    
    # For now only save z,y components (plane projections)
    return [z1, z2], [y1, y2], [x1, x2]


if __name__ == '__main__':
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('-i', '--input', type = str,
                        default = 'pixels.npy',
                        help = 'where read the destinations from')
        parser.add_argument('-o', '--output', type=str,
                        default='evt.npy',
                        help='where to save the destinations')
                                
        args = parser.parse_args()
        outFile = args.output

        oldRecord = np.load(args.input, allow_pickle=True)[0]
        pixels = oldRecord.hitMap
        z, y, t, q = pixels.T
        d = np.array([z,y,v*t]).T

        newRecord = oldRecord
        newRecord.pointsPCA = get_pca(d) #Get PCA (two points in z,y defined by principal axis)

        np.save(outFile, np.array([newRecord]))
