import batman
import numpy as np
from lightkurve import LightCurve

###function to fetch batman lc according to user params (single transit, noise
#   type, pl/star params, etc.)
#   inject noise here

###function to instantiate a batman transit
def batman_transit(period, rp, a, u=[0.4804, 0.1867], t0=0., inc=90., ecc=0.,
                   w=90., limb_dark='quadratic', cadence=0.01):
    """
    Function to instantiate a simulated 'LightCurve' object using the BATMAN
    package developed by Laura Kreidberg. Please see BATMAN documentation for a
    more in-depth description of each parameter.

    Parameters
    ----------
    period : float
       Orbital period in days of the simulated planet.
    rp : float
       Radius of planet in units of stellar radii.
    a : float
       Semi-major axis of planet orbit in units of stellar radii.
    u : array
       Limb darkening coefficients. Defaults to a G2V star in the Kepler 
       bandpass.
    t0 : float
       Mid-transit time.
    inc : float
       Orbital inclination in degrees.
    ecc : float
       Eccentricity of planet orbit. May cause simulation to run slowly if set 
       to a nonzero eccentricity and noise is added using a separate function.
    w : float
       Longitude of periastron in degrees.
    limb_dark : str
       Limb darkening model.
    cadence : float
       Cadence of data points in minutes.

    Returns
    -------
    lc : 'LightCurve' object
       The light curve of the simulated planet transit for one full orbital 
       period.
    """
    params = batman.TransitParams()
    params.t0 = t0
    params.per = period
    params.rp = rp
    params.a = a
    params.inc = inc
    params.ecc = ecc
    params.w = w
    params.u = u
    params.limb_dark = str(limb_dark)

    t = np.linspace(-period/2, period/2, int((period*24*60)//cadence))

    m = batman.TransitModel(params, t)
    flux = m.light_curve(params)

    lc = LightCurve(t, flux, flux_err=None)

    lc.params = params
    
    return lc

    
###function to extend transit to arbitrary length
