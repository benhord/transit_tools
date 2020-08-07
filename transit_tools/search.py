import numpy as np
import transitleastsquares as tls
import astropy
import math
from lightkurve import LightCurve

from .utils import rms

#tls_default #do these really need to be broken up? no, combine

#tls_graze

#tls_box
#def tls_search(lc, time=None, flux=None, dy=None, shape='default'):
def tls_search(*args, tic=None, shape='default', star_params=None,
               rms_calc=True, norm_val=1., clean_lc=False, starparams_out=False,
               del_dur=1., verbose=False, nthreads=6, **kwargs):
    """
    Function to perform a search for periodic signals using the Transit Least 
    Squares (TLS) algorithm developed by Hippke & Heller 2018. While slower than
    Box Least Squares, the transit shape used in the search is more realistic.

    Parameters
    ----------
    *args : 'LightCurve' object or multiple numpy array arguments
       If the len of *args = 1, then the argument is assumed to be a lightkurve
       'LightCurve' object with at least two columns, labeled 'time' and 'flux',
       respectively, with an optional 'flux_err' column as the third column. If 
       the len of *args > 1, then it is assumed the user is passing time, flux, 
       and flux_err (optional), respectively. These columns or arguments should 
       be arrays of equal length.
    tic : int or None
       TIC ID of the source that the light curve comes from. This will be used 
       to query the TIC for the stellar parameters of the system. May be set to
       None if a full dictionary of stellar params are provided to the 
       star_params keyword.
    shape : str
       The shape used by TLS to search for periodic signals. The user may 
       specify 'default', 'grazing', or 'box'. See Hippke & Heller 2018 for an
       in-depth description of these shapes.
    star_params : dict or None
       A dictionary containing stellar parameters to be used in the TLS search.
       The dictionary can contain an array of limb-darkening 
       parameters, stellar radius, lower radius error, upper radius error, 
       stellar mass, lower mass error, and upper mass error labeled 'ab', 
       'rstar', 'rlow', 'rhigh', 'mstar', 'mlow', and 'mhigh', respectively. The
       error values are the errors themselves and not the upper and lower values
       for each of the parameters. A partial list may be included, but in this
       case, the TIC must also be given.
    rms_calc : bool
       A flag to denote whether the root mean square error will be applied in 
       the case that error values are not provided.
    norm_val : float
       Value that the light curve is normalized to. Default is 1. Only 1 or 0 
       are valid normalizations for TLS.
    clean_lc : bool
       Flag to indicate whether or not to output a cleaned lightcurve with the
       recovered periodic signal masked out. Results in an additional expected
       output.
    starparams_out : bool
       Flag to indicate whether or not to output the dictionary of stellar
       parameters used in the TLS search. Results in an additional expected
       output.
    del_dur : float
       How many durations worth of data points should be excluded from cleaned
       light curve centered on the transit center. Default is 1. Values < 1 will
       result in some in-transit points remaining whild values > 1 will remove
       some points outside the transit.
    verbose : bool
       Flag to have function print more while it runs.
    nthreads : int
       Number of threads to be used for running the signal search. Many times,
       cores have the capability to run multiple threads, so be sure to check
       your machine to optimize this parameter.
    kwargs
       Optional arguments passed to the transitleastsquares.power function.

    Returns
    -------
    results : dict
       Results of the TLS fit. See TLS documentation for the contents of this
       dictionary and descriptions of each element.
    cleaned_lc : 'LightCurve' object, optional
       A light curve with the transits masked out based on the results of the 
       TLS search.
    """
    dy = None

    #processing inputs
    if not tic and (not star_params or len(star_params) != 7):
        raise ValueError('Either tic or full star_params dictionary must' +
                         ' be given!')

    if len(args) == 1:
        lc = args[0]
        time = lc.time
        flux = lc.flux
        try:
            dy = lc.flux_err
        except:
            print('No flux errors provided')

    elif len(args) > 1:
        time = args[0]
        flux = args[1]
        if len(args) == 3:
            dy = args[2]
        else:
            print('No flux errors provided')
            
    if rms_calc == True and not isinstance(dy, np.ndarray):
        dy = np.ones(len(flux)) * rms(flux, norm_val=norm_val)
        print('RMS will be used for errors')

    #get catalog info and/or parse user-provided values
    if tic and (not star_params or len(star_params) != 7):
        print(f'Gathering stellar params that were not provided...', end='\r')

        ab, R_star, R_star_lowe, R_star_highe, M_star, M_star_lowe, M_star_highe = tls.catalog_info(TIC_ID=int(tic))
        cat = {'ab' : ab, 'rstar' : R_star, 'rlow' : R_star_lowe,
               'rhigh' : R_star_highe, 'mstar' : M_star, 'mlow' : M_star_lowe,
               'mhigh' : M_star_highe}
        
        if not star_params:
            star_params = {}

        missing = list(set(cat) - set(star_params))
        star_params.update({k: cat[k] for k in missing})
        
        print('Gathering stellar params that were not provided... Done!')

    #quality control for stellar params
    dc = star_params
        
    dc['rstar'] = 1.0 if math.isnan(dc['rstar']) else dc['rstar']
    dc['mstar'] = 1.0 if math.isnan(dc['mstar']) else dc['mstar']
    dc['mlow'] = 0.1 if math.isnan(dc['mlow']) else dc['mlow']
    dc['mhigh'] = 0.1 if math.isnan(dc['mhigh']) else dc['mhigh']
    dc['rlow'] = 0.1 if math.isnan(dc['rlow']) else dc['rlow']
    dc['rhigh'] = 0.1 if math.isnan(dc['rhigh']) else dc['rhigh']

    rmax = dc['rstar'] + dc['rhigh']
    rmin = dc['rstar'] - dc['rlow']
    mmax = dc['mstar'] + dc['mhigh']
    mmin = dc['mstar'] - dc['mlow']

    if verbose:
        print('Stellar params used:')
        for i in list(dc.keys()):
            print(str(i) + ' = ' + str(dc[i]))
        print('(defaults are solar and 0.1 for errors)')
            
    #beginning TLS search
    print('Searching using TLS using %s shape...' % shape)
    model = tls.transitleastsquares(t=time, y=flux, dy=dy)
    results = model.power(R_star=dc['rstar'], R_star_min=rmin, R_star_max=rmax,
                          M_star=dc['mstar'], M_star_min=mmin, M_star_max=mmax,
                          u=dc['ab'], transit_template=shape,
                          use_threads=nthreads, **kwargs)

    #cleaning light curve if lc_clean flag
    if clean_lc:
        intransit = tls.transit_mask(time, results.period,
                                     del_dur * results.duration, results.T0)
        time2 = time[~intransit]
        flux2 = flux[~intransit]
        time2, flux2 = tls.cleaned_array(time2, flux2)

        lc_clean = LightCurve(time=time2, flux=flux2, flux_err=dy[~intransit])

    if clean_lc and not starparams_out:
        return results, lc_clean
    elif not clean_lc and starparams_out:
        return results, dc
    elif clean_lc and starparams_out:
        return results, lc_clean, dc
    else:
        return results
    
#bls from astropy

#river plot tls

#river plot bls
