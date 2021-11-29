import numpy as np
import scipy.stats as st
from numpy import linalg as LA


def mag(vect):
    return np.sqrt(np.sum(np.power(vect, 2)))


def norm(vect):
    return vect/mag(vect)


def get_perpendicular_vectors(vect):
    # TODO
    # perp1 = np.array([1, 0, 0])
    # perp2 = np.array([0, 1, 0])
    perp1 = np.random.randn(3)
    perp1 -= perp1.dot(vect)*vect
    perp1 /= np.linalg.norm(perp1)

    perp2 = np.cross(vect, perp1)

    return perp1, perp2



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
    return [z1, z2], [y1, y2]
