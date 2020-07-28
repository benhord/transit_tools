import numpy as np
from transitleastsquares import transitleastsquares, catalog_info, cleaned_array
import astropy

#tls_default #do these really need to be broken up? no, combine

#tls_graze

#tls_box
#def tls_search(lc, time=None, flux=None, dy=None, shape='default'):
def tls_search(*args, shape='default', star_params=None):
    """
    Function to perform a search for periodic signals using the Transit Least 
    Squares (TLS) algorithm developed by Hippke & Heller 2018. While slower than
    Box Least Squares, the transit shape used in the search is more realistic.

    !!Include option to have function calculate RMS for light curve!!
    !!Include option to have user provide stellar params for search!!
    !!Output both search values and cleaned light curve for future runs!!

    Parameters
    ----------
    *args : 'LightCurve' object or multiple numpy array arguments
       If the len of *args = 1, then the argument is assumed to be a lightkurve
       'LightCurve' object with at least two columns, labeled 'time' and 'flux',
       respectively, with an optional 'flux_err' column as the third column. If 
       the len of *args > 1, then it is assumed the user is passing time, flux, 
       and flux_err (optional), respectively.
    shape : str
       The shape used by TLS to search for periodic signals. The user may 
       specify 'default', 'grazing', or 'box'. See Hippke & Heller 2018 for an
       in-depth description of these shapes.
    star_params : dict or None
       A dictionary containing 

    Returns
    -------
    params : dict
    cleaned_lc : 'LightCurve' object
    """
    if len(args) == 1:
        lc = args[0]
        time = lc.time
        flux = lc.flux
        try:
            flux_err = lc.flux_err
        except:
            flux_err = None #placeholder for RMS
            print('No errors provided')

    elif len(args) > 1:
        time = args[0]
        flux = args[1]
        try:
            flux_err = args[2]
        except:
            flux_err = None #placeholder for RMS
            print('No errors provided')

    #get catalog info or parse user-provided values
            
    if flux_err:
        model = transitleastsquares(time=time, flux=flux, dy=flux_err,
                                    shape=shape)
    else:
        model = transitleastsquares(time=time, flux=flux, shape=shape)
        
    print('Searching using TLS using %s shape...' % shape)
    #check to make sure lc is lc object and if not, use time, flux, dy
    
#bls from astropy

#river plot tls

#river plot bls
