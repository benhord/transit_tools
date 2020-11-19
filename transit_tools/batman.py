import batman
import numpy as np
from lightkurve import LightCurve
from pathlib import Path

from .utils import name_to_tic
from .fetch_lc import gather_lc

###function to fetch batman lc according to user params (single transit, noise
#   type, pl/star params, etc.)
#   inject noise here
def full_batlc(period, rp, a, noise='obs', inlc=None, t0=None, sectors='all',
               cadence='2min', **kwargs):
    """
    Function to retrive a full simulated BATMAN light curve complete with 
    injected noise.

    !!Allow custom noise model!!
    !!Update name processing once name processing wrapper function written!!
    !!Utilize actual stellar params somehow?!!
    !!Search for known planets in the real system and pass those on as well!!
    !!Write docstrings!!
    !!Pass which TIC ID was used!!
    !!Allow for generation of gaussian noise!!
    !!Allow custom lc for inlc!!
    !!Raise error for noise='obs' and inlc=None!!

    Parameters
    ----------
    period : float
       Orbital period to be simulated. Should be in days.
    rp : float
       Radius of planet in units of host star radii. Essentially Rp/Rs.
    a : float
       Semi-major axis of the orbit in units of host star radii.
    noise : str or None
       Noise to be included in simulation. Current options are 'obs' for 
       injection into pre-existing light curve or None for a simulated light
       curve without noise. Gaussian noise is expected in the future.
    inlc : None or int
       TIC ID of light curve for simulations to be injected into. Can be set to
       None for a light curve to be randomly chosen from a pre-selected list in
       transit_tools/files/lc_list.csv if noise='obs'.
    t0 : float or None
       Mid-transit time of the first transit. If None, t0 will be randomly 
       generated to be somewhere within 1 period of the start of the 
       observation. Must be in the same units as period.
    sectors : str or array
       Sectors to be considered for retrieving light curves specified by inlc.
       If 'all', all available light curves for selected inlc will be 
       retrieved, otherwise only those contained in the user-provided array will
       be retrieved.
    cadence : str
       Cadence of TESS light curve to retrieve. Options are '2min', 'ffi_ml', or
       'eleanor'. This method is passed to transit_tools.fetch_lc.gather_lc and
       follows those hierarchy rules if a light curve at the given cadence
       cannot be found.
    kwargs
       Additional arguments to be passed to transit_tools.batman_transit.
    """
    tic = None
    
    if not inlc and noise == 'obs':
        #randomly choose lc to inject into from preselected list
        path = Path(__file__).parent / "./files/lc_list.csv"
        lc_list = np.genfromtxt(path, delimiter=',')
        lc_list = np.delete(lc_list, 0, 0)[:, 0]

        inlc = int(np.random.choice(lc_list))
        tic = inlc

    if isinstance(inlc, str):
        #perform name processing and set lc to TIC ID
        inlc = name_to_tic(str(inlc))
        tic = inlc
        
    if isinstance(inlc, int):
        #inject into chosen TIC ID using gather_lc
        tic = inlc
        inlc, sectors = gather_lc(inlc, sectors=sectors, return_sectors=True,
                                  method=cadence)
        
    if inlc is not None:
        #inject into chosen lc
        time = inlc.time
        flux = inlc.flux
        flux_err = None
        if hasattr(inlc, 'flux_err'):
            flux_err = inlc.flux_err

        if not t0:
            t0_range = time[time < (time[0] + period)]
            t0 = np.random.choice(t0_range)
        else:
            t0 = t0
            
        model = batman_transit(period=period, rp=rp, a=a, t0=t0, time=time,
                               **kwargs)

        params = model.params
        
        flux = flux + model.flux - 1

        inlc = LightCurve(time, flux, flux_err=flux_err)

        inlc.sectors = sectors
        inlc.params = params
        inlc.tic = tic
        
    elif noise == None:
        inlc = batman_transit(period=period, rp=rp, a=a, t0=t0, **kwargs)
        inlc.tic = None
        inlc.sectors = None
        
    inlc.known_pls = {'orbital_period' : inlc.params.per,
                      't0' : inlc.params.t0,
                      'rprs' : inlc.params.rp, 'a' : inlc.params.a,
                      'inc' : inlc.params.inc, 'ecc' : inlc.params.ecc,
                      'w' : inlc.params.w, 'u' : inlc.params.u,
                      'limb_dark' : inlc.params.limb_dark,
                      'pl_name' : 'batman'}
        
    return inlc
        
def batman_transit(period, rp, a, u=[0.4804, 0.1867], t0=0., inc=90., ecc=0.,
                   w=90., limb_dark='quadratic', cadence=0.01, length=None,
                   time=None, **kwargs):
    """
    Function to generate a simulated 'LightCurve' object using the BATMAN
    package developed by Laura Kreidberg. Please see BATMAN documentation for a
    more in-depth description of each parameter.

    !!Fix length argument. Make it remove parts of light curve from boths sides
      of transit!!
    !!Integrate supersample_factor and exp_time arguments better since they're
      likely going to be common for TESS. Maybe automatically do it for user!!

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
    kwargs
       Additional arguments to be passed to batman.TransitModel. For TESS
       light curves, it is recommended to specify the supersample_factor and
       exp_time arguments for best results.

    Returns
    -------
    lc : 'LightCurve' object
       The light curve of the simulated planet transit for one full orbital 
       period. Light curve will have attribute 'params' that is an object
       containing all of the input simulated parameters.
    """
    if length is not None and time is not None:
        raise ValueError('Please only specify either length or time')

    if t0 is None:
        t0 = 0.
    
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

    if time is None:
        t = np.linspace(-period/2, length-(period/2),
                        int((period*24*60)//cadence))
    else:
        t = time
        
    m = batman.TransitModel(params, t, **kwargs)
    flux = m.light_curve(params)

    lc = LightCurve(t, flux, flux_err=None)

    lc.params = params
    
    return lc
