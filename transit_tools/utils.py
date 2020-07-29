import numpy as np
from astroquery.mast import Catalogs, Observations

#function to import catalog info for source
   #stellar information
   #RA, DEC
   #whatever is in the TIC
   #GAIA position/other info?

   #break into separate functions to provide greater flexibility?
   #for any stellar information, put any relevant info into both a self.star
   #  attribute but also a self.star_params_tls attribute

#function to import observation information
   #observation information based on self.method keyword (eg ccd, sector, etc.)

#function to find and report data on planets already-known in the system

#function to update stellar params of lightcurve object (do automatically?,
#   keyword to update in the fetch command?)

#common name processing so that a TIC isn't the only thing that can be passed
   #maybe first have general parser that gives canonical name then second stage
   # that queries TIC for canonical name
def name_to_tic(name):
    """
    Function to convert common name to TIC ID. Queries the MAST for TIC entry
    nearest to known position for common name.

    Parameters
    ----------
    name : str
       Common name to be converted to TIC.

    Returns
    -------
    tic : int
       TIC ID of closest match to input name's position from TIC on MAST.
    """
    if not isinstance(name, str):
        raise ValueError('Name must be a string.')

    cat = Catalogs.query_object(name, radius=0.02, catalog="TIC")
    tic = int(cat[0]['ID'])
    
    return tic

def tic_to_name(tic):
    """
    Function to determine the common name of a TIC ID, if it has one.

    Parameters
    ----------
    tic : int
       The TIC ID of the object for which the common name is desired.

    Returns
    -------
    name : str
       The common name of the input TIC ID.
    """
    if not isinstance(tic, int):
        raise ValueError('TIC must be an integer')

    cat = Catalogs.query_criteria(catalog="TIC", ID=int(tic))
    ra = cat[0]['ra']
    dec = cat[0]['dec']
    
    results = Simbad.query_region(coord.SkyCoord(ra, dec, unit=(u.deg, u.deg)),
                                  radius='0d0m5s')

    name = str(results[0]['MAIN_ID'].decode('utf-8'))
    
    return str(name)

#calculate rms of a data set
def rms(data, norm_val=1.):
    """
    Calculates the Root Mean Square of the provided data.

    Parameters
    ----------
    data : numpy array
       Data for which the root mean square will be calculated
    norm_val : float or numpy array
       The value(s) that the data array is normalized to or the value of the
       model that the data values are being compared against if RMS is non-
       uniform.

    Returns
    -------
    rms : float
       The calculated root mean square of the data.
    """
    rms = np.sqrt(np.sum((data - norm_val) ** 2) / len(data))

    return rms

#method to display best fit planet parameters from search
def search_summary(results, routine='tls'):
    """
    Function to display periodic signal search results in a user-friendly 
    format.

    !!Update to print BLS parameters as well!!
    !!Upgrade to generic print statements agnostic of routine. (Just search for
    keywords and print whichever are available!!
    !!Consider moving to search.py!!
    !!Propogate uncertainties to derived quantities where TLS doesn't!!
    !!Format sig figs better!!

    Parameters
    ----------
    results : dict
       Dictionary containing results from the signal search
    routine : str
       String defining which search method was used. Important for formatting
       the outputs.
    """
    if routine == 'tls':
        print('Period = %.5f +/- %.5f' %
              (results.period, results.period_uncertainty))
        print('t0 = %.5f' % results.T0)
        print('Duration = %.5f' % results.duration)
        print('Avg depth = %.5f +/- %.5f' %
              (results.depth_mean[0], results.depth_mean[1]))
        print('SDE = %.5f' % results.SDE)
        print('FAP = %.5e' % results.FAP)
        print('Odd-Even difference = %.5f' % results.odd_even_mismatch)
        print('SNR = %.5f' % results.snr)
        print('Rp/Rs = %.5f' % results.rp_rs)
        print('Transit count = %s' % results.transit_count)
        print('')

    if routine == 'bls':
        print('Still working on this functionality! Please be patient.')
