import batman
import numpy as np
from lightkurve import LightCurve
from pathlib import Path

from .utils import name_to_tic
from .fetch_lc import gather_lc

###function to fetch batman lc according to user params (single transit, noise
#   type, pl/star params, etc.)
#   inject noise here
def full_batlc(period, rp, a, noise='obs', lc=None, t0=None, sectors='all',
               **kwargs):
    """
    Function to retrive a full simulated BATMAN light curve complete with 
    injected noise.

    !!Allow custom noise model!!
    !!Allow injection into custom light curve w/ lc arg instead of tic!!
    !!Update name processing once name processing wrapper function written!!
    !!subtract 1 from simulated lc then add it to real lc!!

    Parameters
    ----------
    """
    if not lc and noise == 'obs':
        #randomly choose lc to inject into from preselected list
        path = Path(__file__).parent / "./files/lc_list.csv"
        lc_list = np.genfromtxt(path, delimiter=',')
        lc_list = np.delete(lc_list, 0, 0)[:, 0]

        lc = int(np.random.choice(lc_list))
        tic = lc
        print(tic)

    if isinstance(lc, str):
        #perform name processing and set lc to TIC ID
        lc = name_to_tic(str(lc))
        tic = lc
        
    if isinstance(lc, int):
        #inject into chosen TIC ID using gather_lc
        tic = lc
        lc, sectors = gather_lc(lc, sectors=sectors, return_sectors=True)
        
    if lc is not None:
        #inject into chosen lc
        time = lc.time
        flux = lc.flux
        flux_err = None
        if hasattr(lc, 'flux_err'):
            flux_err = lc.flux_err

        if not t0:
            t0_range = time[time < (time[0] + period)]
            t0 = np.random.choice(t0_range)
        else:
            t0 = t0
            
        model = batman_transit(period=period, rp=rp, a=a, t0=t0, time=time,
                               **kwargs)

        flux = flux + model.flux - 1

        lc = LightCurve(time, flux, flux_err=flux_err)
        
    else:
        lc = batman_transit(period=period, rp=rp, a=a, **kwargs)

    #add sim params to dict and add as attr to lc
        
    return lc
        #return clean BATMAN light curve w/o noise
        
def batman_transit(period, rp, a, u=[0.4804, 0.1867], t0=0., inc=90., ecc=0.,
                   w=90., limb_dark='quadratic', cadence=0.01, length=None,
                   time=None):
    """
    Function to generate a simulated 'LightCurve' object using the BATMAN
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
       Mid-transit time of the first transit.
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
    length : float or None
       Length of output light curve in days. If set to None, defaults to one 
       orbital period in length centered on the first transit.
    time : array
       Array of time points during which the light curve will be simulated.

    Returns
    -------
    lc : 'LightCurve' object
       The light curve of the simulated planet transit for one full orbital 
       period. Light curve will have attribute 'params' that is an object
       containing all of the input simulated parameters.
    """
    if length is not None and time is not None:
        raise ValueError('Please only specify either length or time')
    
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

    if not length:
        length = period

    if time.all() is None:
        t = np.linspace(-period/2, length-(period/2),
                        int((period*24*60)//cadence))
    else:
        t = time
        
    m = batman.TransitModel(params, t)
    flux = m.light_curve(params)

    lc = LightCurve(t, flux, flux_err=None)

    lc.params = params
    
    return lc
