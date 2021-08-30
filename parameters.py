from utils import *

def drift_model(E, temperature = 89):
    magE = mag(E)
    
    if magE > 4.0:
        magE = 4.0

    if ( temperature < 87.0 ) or ( temperature > 94.0 ):
        print ( "DriftModel Warning!: Temperature value of ")
        print ( temperature )
        print ( "K is outside of range covered by drift velocity" )
        print ( " parameterization.  Returned value may not be correct" )

    tShift = -87.203 + temperature
    xFit = 0.0938163 - 0.0052563*tShift - 0.0001470*tShift*tShift
    uFit = 5.18406 + 0.01448*tShift - 0.003497*tShift*tShift - 0.000516*tShift*tShift*tShift
    
    # Icarus Parameter set
    # Use as default
    P1 = -0.04640 # K^-1
    P2 = 0.01712 # K^-1
    P3 = 1.88125 # (kV/cm)^-1
    P4 = 0.99408 # kV/cm
    P5 = 0.01172 # (kV/cm)^-P6
    P6 = 4.20214
    T0 = 105.749 # K

    # Walkowiak Parameter set
    P1W = -0.01481; # K^-1
    P2W = 0.0075; # K^-1
    P3W = 0.141; # (kV/cm)^-1
    P4W = 12.4; # kV/cm
    P5W = 1.627; # (kV/cm)^-P6
    P6W = 0.317;
    T0W = 90.371; # K

    # From Craig Thorne . . . currently not documented
    # smooth transition from linear at small fields to
    # icarus fit at most fields to Walkowiak at very high fields
    if magE < xFit:
        vd = magE*uFit
    elif magE < 0.619:
        vd = ((P1*(temperature-T0)+1)
	      *(P3*magE*np.log(1+P4/magE) + P5*np.power(magE, P6))
	      +P2*(temperature-T0))
    elif magE < 0.699:
        vd = 12.5*(magE - 0.619)*((P1W*(temperature-T0W)+1)
			          *(P3W*magE*np.log(1+P4W/magE) + P5W*np.power(magE, P6W))
			          +P2W*(temperature-T0W))
        + 12.5*(0.699 - magE)*((P1*(temperature-T0)+1)
			     *(P3*magE*np.log(1+P4/magE) + P5*np.power(magE, P6))
			     +P2*(temperature-T0))
    else:
        vd = ((P1W*(temperature-T0W)+1)
	      *(P3W*magE*np.log(1+P4W/magE) + P5W*np.power(magE, P6W))
	      +P2W*(temperature-T0W))

    vd /= 10

    return vd; # cm/us

def betadEdx(E):
    # do a thing return dEdx
    #
    return 1 # cm
    
def betaRange(E):
    """ this function returns a length of a tracklet from a given initial E """
    # read file
    # interpolate
    # return range
    return 10 # cm

def emission_spectrum(thisSource = 'Cs137'):
    if thisSource == 'Cs137':
        return 0.512 # [MeV]
    else:
        raise Exception ("This is not a valid source!")
    
# need to add recombination and work
physics_parameters = {"DT": 8.8e-6,         # transverse diffusion,       cm * cm / us
                      "DL": 4.0e-6,         # longitudinal diffusion,     cm * cm / us
                      "v":  drift_model,    # drift velocity (function),  cm / us
                      "lt": 10.e3,          # lifetime,                   us
                      "npe": 100,           # number of photoelectrons
                      "npe_sigma": 10,      # error on number of pe
                      "dEdx": betadEdx,
                      "R": 0.66,            # recombination factor
                      "w": 23.6e-6}         # ionization w.f. MeV
                      
sim_parameters = {"dt": 5.e-1,              # time step for integrating the electron's path,   us
                  "scalingF": 100}         # electrons per charge bundle

detector_parameters = {"cathode position": 50, # distance from cathode to anode,     cm
                       "target radius": 0.2,   # radius of the cathode target,       cm
                       "noise level": 200,     # ENC (electrons)
                       "nominal field": 0.5,   # nominal field strength,             kV / cm
                       "pixel threshold": 4e2 }  # threshold for pixel hit ,           # e-
