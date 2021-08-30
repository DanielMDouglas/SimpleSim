import numpy as np
import scipy.stats as st
    
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
