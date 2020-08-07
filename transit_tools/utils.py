import numpy as np
from astroquery.mast import Catalogs, Observations, Tesscut
from astroquery.simbad import Simbad
from urllib import request
import astropy.coordinates as coord
import astropy.units as u
import operator

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
def tessobs_info(tic=None, ra=None, dec=None):
    """
    Function to retrieve observation information for objects observed by TESS.

    !!Update to include exp time, pixel location, other observation-specific
      quantities!!

    Parameters
    ----------
    tic : int or None
       TIC ID of target to be queried. Must not be None if ra and dec are None.
    ra : float or None
       RA of target to be queried. Must not be None if tic is None.
    dec : float or None
       Dec of target to be queried. Must not be None if tic is None.

    Returns
    -------
    info : dict
       Dictionary continaing TESS observation info.
    """
    if not tic and not ra and not dec:
        raise ValueError('Please provide either a TIC ID or both RA and Dec')

    if not ra or not dec:
        cat = Catalogs.query_criteria(catalog="TIC", ID=int(tic))
        ra = cat[0]['ra']
        dec = cat[0]['dec']

    coords = coord.SkyCoord(ra, dec, unit='deg')
    sector_table = Tesscut.get_sectors(coordinates=coords)
    
    if len(sector_table) == 0:
        print('Object not observed by TESS')

    sec_name = []
    sec = []
    cam = []
    ccd = []

    for i in range(len(sector_table)):
        sec_name.append(sector_table[i]['sectorName'])
        sec.append(sector_table[i]['sector'])
        cam.append(sector_table[i]['camera'])
        ccd.append(sector_table[i]['ccd'])

    info = {'sectorName' : sec_name, 'sector' : sec, 'camera' : cam,
            'ccd' : ccd}

    return info

def coord_to_tic(ra, dec):
    """
    Function to convert input RA and Dec coordinates to the nearest TIC ID from 
    the TESS Input Catalog (TIC).

    Parameters
    ----------
    ra : float
       The RA of the target source.
    dec : float
       The Dec of the target source.

    Returns
    -------
    tic : int
       TIC ID of the source nearest to the input RA and Dec from the TIC.
    """
    cat = Catalogs.query_region((str(ra) + ' ' + str(dec)), catalog="TIC")
    tic = int(cat[0]['ID'])
    
    return tic

def known_pls(name=None, ra=None, dec=None, verbose=False):
    """
    A function to gather information on any known planets in a given system. 
    Queries Simbad for objects and queries the Gaia catalog if RA/Dec are not
    provided.

    Parameters
    ----------
    name : str
       Common name of the system being checked. Optional if RA/Dec are provided.
    ra : float
       RA in decimal degrees. Optional if name is provided. If provided, Dec is
       also required.
    Dec : float
       Dec in decimal degrees. Optional if name is provided. If provided, RA is
       also required.
    verbose : bool
       Flag to determine whether some of the parameters of the known planets in
       the system are printed.

    Returns
    -------
    pl_info : list of dicts
       List containing a dictionary of all known planet parameters for each 
       planet in the queried system.
    """
    if not name and not ra and not dec:
        raise ValueError('Either name or both RA & Dec must be provided')

    if not ra and not dec:
        results = Catalogs.query_object(str(name), radius=0.02, catalog="Gaia")
        ra = results[0]['ra']
        dec = results[0]['dec']

    cat = Simbad.query_region(coord.SkyCoord(ra, dec, unit=(u.deg, u.deg)),
                              radius='0d0m5s')

    print(cat)
    
    pls = len(cat) - 1
    if verbose:
        print(str(pls) + ' known planets found in system')
        print('Gathering info for each planet...')

    pl_info = []
        
    if pls > 0:
        alphastr = 'abcdefghijklmnopqrstuvwxyz'
        null = None
        true = True
        false = False

        for i in range(pls):
            #pl = str(cat[i+1]['MAIN_ID'].decode('utf-8'))
            pl = (str(sorted(cat, key=operator.itemgetter('MAIN_ID'))[i+1]
                      ['MAIN_ID'].decode('utf-8')))
            urlname = (pl[:-1] + "%20" + pl[-1])

            print(urlname)
            
            link = ("https://exo.mast.stsci.edu/api/v0.1/exoplanets/" +
                    urlname + "/properties")
            info = request.urlopen(link).read().decode('utf-8')
            info = eval(info)

            try:
                pl_info.append(info[0])
                query_fault = False
            except:
                query_fault = True
                pl_info.append(str(pl))
                
        if verbose and not query_fault:
            pl = pl_info
            for n in range(pls):
                print(pl[n]['canonical_name'])
                print('Period = %.5f +/- %.5f %s' %
                      (pl[n]['orbital_period'],
                       pl[n]['orbital_period_upper'],
                       pl[n]['orbital_period_unit']))
                if pl[n]['transit_time']:
                    print('t0 = %.5f +2457000 BTJD +/- %.5f' %
                          (pl[n]['transit_time']-56999.5,
                           pl[n]['transit_time_upper']))
                else:
                    print('t0 = None')
                print('Duration = %.5f +/- %.5f %s' %
                      (float(pl[n]['transit_duration'] or 0),
                       float(pl[n]['transit_duration_upper'] or 0),
                       str(pl[n]['transit_duration_unit'] or 'None')))
                print('Transit depth = %.5f +/- %.5f' %
                      (float(pl[n]['transit_depth'] or 0),
                       float(pl[n]['transit_depth_upper'] or 0)))
                print('Orbital distance = %.5f +/- %.5f %s' %
                      (float(pl[n]['orbital_distance'] or 0),
                       float(pl[n]['orbital_distance_upper'] or 0),
                       str(pl[n]['orbital_distance_unit'] or 'None')))
                print('Rp/Rs = %.5f +/- %.5f' %
                      (float(pl[n]['Rp/Rs'] or 0),
                       float(pl[n]['Rp/Rs_upper'] or 0)))
                print('Radius = %.5f +/- %.5f %s' %
                      (float(pl[n]['Rp'] or 0), float(pl[n]['Rp_upper'] or 0),
                       str(pl[n]['Rp_unit'] or 'None')))
                print('Mass = %.5f +/- %.5f %s' %
                      (float(pl[n]['Mp'] or 0), float(pl[n]['Mp_upper'] or 0),
                       str(pl[n]['Mp_unit'] or 'None')))
                print('Disposition: %s' % pl[n]['disposition'])
                print('')

        elif verbose and query_fault:
            print('Known planet found but parameters were not found in MAST' +
                  'for some reason.')
                
    return pl_info

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

def tic_to_name(tic, ra=None, dec=None):
    """
    Function to determine the common name of a TIC ID or given RA/Dec position, 
    if it has one. Queries the MAST and Simbad to gather this information.

    Parameters
    ----------
    tic : int
       The TIC ID of the object for which the common name is desired.
    ra : float
       The RA in decimal degrees. Optional with Dec to circumvent querying MAST.
    dec : float
       The Dec in decimal degrees. Optional with TA to circumvent querying MAST.

    Returns
    -------
    name : str
       The common name of the input TIC ID.
    """
    if not isinstance(tic, int):
        raise ValueError('TIC must be an integer')

    if not ra and not dec:
        cat = Catalogs.query_criteria(catalog="TIC", ID=int(tic))
        ra = cat[0]['ra']
        dec = cat[0]['dec']
    
    results = Simbad.query_region(coord.SkyCoord(ra, dec, unit=(u.deg, u.deg)),
                                  radius='0d0m5s')

    name = str(results[0]['MAIN_ID'].decode('utf-8'))
    
    return str(name)

##def name_processing(): #wrap all name processing into one fn, search for TOIs
  #intake coords as well.

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
        print('Period = %.5f +/- %.5f d' %
              (results.period, results.period_uncertainty))
        print('t0 = %.5f BTJD' % results.T0)
        print('Duration = %.5f d' % results.duration)
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

#def overlap(range1, range2):
    """
    Function to determine if two value ranges overlap.

    !!Revisit later. Not urgent right now!!

    Parameters
    ----------
    range1 : list
       2 element list with the lower and upper bounds of the first range, 
       respectively.
    range2 : list
       2 element list with the lower and upper bound of the second range,
       respectively.

    Returns
    -------
    intersect : bool
       Boolean indicating whether the ranges overlap or not.
    """

def fold(time, flux, period, flux_err=None, midpoint=None):
    """
    Folds a timeseries based on a given period with the option to provide a 
    midpoint around which to fold.

    !!Doesn't work with current version. Needs update!!

    Parameters
    ----------
    time : array
       The time of the timeseries that is to be folded. Must be identical in
       length to the flux array.
    flux : array
       The flux of the timeseries that is to be folded. Must be identical in 
       length to the time array.
    period : float
       The period with which the timeseries will be folded. Must be given in the
       same units as the time array.
    flux_err : array or None
       The corresponding error values for the timeseries. Must be identical in
       length to the time and flux arrays. If not None, an additional output
       will be expected.
    midpoint : float or None
       The point around which the timeseries will be folded. If set to None, the
       first element in the time array will be used.

    Returns
    -------
    folded_time : array
       Array with the folded time values.
    folded_flux : array
       Array with the folded flux values.
    folded_flux_err : array
       Array with the folded flux_err values.
    """
    if not midpoint:
        midpoint = time[0]

    epoch = np.floor((time - midpoint) / period)
    folded_time = (time - midpoint) / period - epoch

    sortIndi = np.argsort(folded_time)
    folded_time = folded_time[sortIndi]
    folded_flux = np.roll(flux[sortIndi], len(flux) // 2)
    if flux_err is None:
        pass
    else:
        folded_flux_err = np.roll(flux_err[sortIndi], len(flux_err) // 2)

    if flux_err is None:
        return folded_time, folded_flux
    else:
        return folded_time, folded_flux, folded_flux_err
