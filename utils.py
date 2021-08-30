import numpy as np

def sample_diffused_pdf(v):
    """
    get the length of the step in the axial (z) and radial (rho) directions.
    These directions are defined in terms of the local E field
    """
    # get a random step (in rho & z) 
    DT = physics_parameters["DT"]
    DL = physics_parameters["DL"]

    dt = sim_parameters["dt"]

    rho = st.norm.rvs(scale = np.sqrt(2*8*DL*dt))
    z = st.norm.rvs(loc = v*dt, scale = np.sqrt(2*8*DT*dt))
    return rho, z

def sample_azimuthal_pdf():
    """
    get the azimuthal (where the axial is the drift direction) direction for a drift step.
    There is no preferred direction, so the distribution is flat
    """
    theta = st.uniform.rvs(scale = 2*np.pi)
    return theta

def sample_from_cathode_target():
    """ 
    get a random position on the cathode target (8mm diameter on the cathode plane)
    """
    z = detector_parameters["cathode position"]
    R = detector_parameters["target radius"]

    rho = np.power(st.uniform.rvs()*np.power(R, 2), 0.5)
    phi = st.uniform.rvs(scale = 2*np.pi)
    
    x = rho*np.cos(phi)
    y = rho*np.sin(phi) 

    pos = np.array([x, y, z])
    return pos

def sample_from_hemisphere():
    """
    return a random unit vector where the z component is negative
    """
    az = 2*np.pi*st.uniform.rvs() # the distribution of azimuthal angles is flat
    pol = np.arccos(-st.uniform.rvs()) # the distribution of polar angles is weighted for even coverage per solid angle

    x = np.cos(az)*np.sin(pol)
    y = np.sin(az)*np.sin(pol)
    z = np.cos(pol)

    return np.array([x, y, z])
    
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
