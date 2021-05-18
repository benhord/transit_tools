import numpy as np
from astroquery.mast import Catalogs, Observations, Tesscut
from astroquery.simbad import Simbad
from urllib import request
import astropy.coordinates as coord
import astropy.units as u
import operator
import sys
import http.client as httplib
from urllib.parse import quote as urlencode
import json
from astroquery.nasa_exoplanet_archive import NasaExoplanetArchive as nea
from astropy.time import Time

#function to import catalog info for source
   #stellar information
   #RA, DEC
   #whatever is in the TIC
   #GAIA position/other info?

   #break into separate functions to provide greater flexibility?
   #for any stellar information, put any relevant info into both a self.star
   #  attribute but also a self.star_params_tls attribute
def catalog_info(tic=None, ra=None, dec=None, cat='all', out_cat=False):
    """
    Function to fetch catalog info about the target system from various 
    catalogs. The GAIA catalog takes precedence over the TIC in the values
    where they overlap, such as RA and Dec.

    !!Change tic to full name processing!!
    !!Change to iterable catalog queries based on user input!!
    !!Allow for addition of new keywords rather than updating!!

    Parameters
    ----------
    tic : int or None
       TIC ID of target object. If None, ra and dec must both not be None.
    ra : float or None
       RA of target object. If None, tic must not be None.
    dec : float or None
       Dec of target object. If None, tic must not be None.
    cat : string
       Catalog(s) to be queried. Currently, only 'tic' and 'gaia' are supported.
    out_cat : bool
       Flag determining whether the catagories searched will be output as an 
       array of strings. If True, an additional output will be expected.

    Returns
    -------
    info : dict
       Dictionary of stellar parameters for target object.
    catalogs : array, optional
       Array containing strings denoting which catalogs were queried for stellar
       parameters.
    """
    if tic is None and (ra is None or dec is None):
        raise ValueError('Please enter either a TIC ID or RA/Dec pair')

    catalogs = []
    info = None
    
    if tic and (ra is None or dec is None):
        query = Catalogs.query_object(('TIC ' + str(tic)), catalog='TIC')
        info = [dict(zip(query.colnames, row)) for row in query][0]
        
        ra = info['ra']
        dec = info['dec']

        if cat != 'tic' and cat != 'TIC' and cat != 'all':
            info = None
        else:
            catalogs.append('tic')
        
    if info is None and (cat == 'tic' or cat == 'TIC' or cat == 'all'):
        query = Catalogs.query_object((str(ra) + ' ' + str(dec)), catalog='TIC')
        info = [dict(zip(query.colnames, row)) for row in query][0]
        
        catalogs.append('tic')

    if cat == 'gaia' or cat == 'GAIA' or cat == 'all':
        query = Catalogs.query_object((str(ra) + ' ' + str(dec)),
                                      catalog='gaia')
        if info is not None:
            info.update([dict(zip(query.colnames, row)) for row in query][0])
        else:
            info = [dict(zip(query.colnames, row)) for row in query][0]

        catalogs.append('gaia')

    if out_cat:
        return info, catalogs
    else:
        return info

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

def known_pls(name=None, ra=None, dec=None, radius=5.0, table='exoplanets',
              values='all', verbose=False):
    """
    A function to gather information on any known planets in a given system. 
    Queries the NASA Exoplanet Archive for objects and their known parameters.

    !!Reduce number of columns queried with each iteration for shorter runtime!!
    !!Allow to search for planets that are not confirmed on the archive!!
    !!Allow user to pass values through as list to be queried for more specific
         query values!!

    Parameters
    ----------
    name : str
       Common name of the system being checked. Optional if RA/Dec are provided.
    ra : float
       RA in decimal degrees. Optional if name is provided. If provided, Dec is
       also required.
    dec : float
       Dec in decimal degrees. Optional if name is provided. If provided, RA is
       also required.
    radius : float
       Radius in arcseconds around which the provided RA and Dec will be searched
       for planets.
    table : str
       Specifies the table to search for planet parameters. See documentation on
       the Exoplanet Archive website for a full list of possible tables and their
       contents. Default is the 'exoplanets' table, which is the default for the
       Exoplanet Archive.
    values : str
       Specifies how many values are collected. Current supported options are
       'minimum' and 'all'.
    verbose : bool
       Flag to determine whether some of the parameters of the known planets in
       the system are printed.

    Returns
    -------
    info : list of dicts
       List containing a dictionary of all known planet parameters for each 
       planet in the queried system.
    """
    if not name and (not ra or not dec):
        raise ValueError('Either name or both RA & Dec must be provided')

    if values == 'minimum':
        select = ('pl_name, pl_orbper, pl_orbpererr1, pl_orbpererr2, ' +
                  'pl_tranmid, pl_tranmiderr1, pl_tranmiderr2, pl_trandur, ' +
                  'pl_trandurerr1, pl_trandurerr2, pl_trandep, pl_trandeperr1,' +
                  ' pl_trandeperr2, pl_ratdor, pl_ratdorerr1, pl_ratdorerr2, ' +
                  'pl_ratror, pl_ratrorerr1, pl_ratrorerr2, pl_radj, ' +
                  'pl_radjerr1, pl_radjerr2, pl_bmassj, pl_bmassjerr1, ' +
                  'pl_bmassjerr2, pl_hostname')
    else:
        select = '*'
        
    
    if name is not None:
        results = nea.query_object(str(name), table=table, select=select)

    elif ra is not None and dec is not None:
        results = nea.query_region(
            coordinates=coord.SkyCoord(ra * u.deg, dec * u.deg),
            radius = radius * u.arcsec,
            table=table,
            select=select
            )
    
    pls = len(results)
    if verbose:
        print(str(pls) + ' known planets found in system')
        print('Gathering info for each planet...')
        
    if pls > 0:
        names = results.colnames
        info = [dict(zip(names, row)) for row in results]

        for i in range(pls):
            info[i]['t0'] = (Time(val=info[i]['pl_tranmid'].value,
                                  format='jd').to_value(format='mjd') - 56999.5)
        
        query_fault = False
    else:
        query_fault = True
                
    if verbose and not query_fault:
        for i in range(pls):
            try:
                pl = info[i]
                print(pl['pl_name'])
                print('Period = %.5f +%.5f %.5f %s' %
                      (pl['pl_orbper'].value,
                       pl['pl_orbpererr1'].value,
                       pl['pl_orbpererr2'].value,
                       pl['pl_orbper'].unit))
                print('t0 = %.5f +2457000 BTJD +%.5f %.5f' %
                      (pl['t0'],
                       pl['pl_tranmiderr1'].value,
                       pl['pl_tranmiderr2'].value))
                print('Duration = %.5f +%.5f %.5f %s' %
                      (pl['pl_trandur'].value,
                       pl['pl_trandurerr1'].value,
                       pl['pl_trandurerr2'].value,
                       pl['pl_trandur'].unit))
                print('Transit depth = %.5f +%.5f -%.5f' %
                      (pl['pl_trandep'].value,
                       pl['pl_trandeperr1'].value,
                       pl['pl_trandeperr2'].value))
                print('a/Rs = %.5f +%.5f %.5f' %
                      (pl['pl_ratdor'],
                       pl['pl_ratdorerr1'],
                       pl['pl_ratdorerr2']))
                print('Rp/Rs = %.5f +%.5f %.5f' %
                      (pl['pl_ratror'],
                       pl['pl_ratrorerr1'],
                       pl['pl_ratrorerr2']))
                print('Radius = %.5f +%.5f %.5f %s' %
                      (pl['pl_radj'].value,
                       pl['pl_radjerr1'].value,
                       pl['pl_radjerr2'].value,
                       pl['pl_radj'].unit))
                print('Mass = %.5f +%.5f %.5f %s' %
                      (pl['pl_bmassj'].value,
                       pl['pl_bmassjerr1'].value,
                       pl['pl_bmassjerr2'].value,
                       pl['pl_bmassj'].unit))
                print('')
            except:
                continue

    elif verbose and query_fault:
        print('Known planet found but parameters were not found in Exoplanet ' +
              'Archive for some reason.')
        info = None
                
    return info

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

    !!Keysort so planet doesn't come first?!!

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

def mastQuery(request):
    """
    Function to make queries to the MAST easier. Helper function. See MAST site
    for more details.

    Parameters
    ----------
    request : dict
       Dictionary of parameters that constitute the query to the MAST.

    Returns
    -------
    head : str
       String containing header information.
    content : str
       String containing desired information.
    """
    server='mast.stsci.edu'

    version = ".".join(map(str, sys.version_info[:3]))

    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain", "User-agent":"python-requests/"+version}

    requestString = json.dumps(request)
    requestString = urlencode(requestString)

    conn = httplib.HTTPSConnection(server)

    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    conn.close()

    return head, content

def canonical_name(name):
    """
    Function to obtain the canonical name of an object from the MAST.
    
    Parameters
    ----------
    name : str
       Name of object to be queried for canonical name.
    
    Results
    -------
    canonical_name : str
       Canonical name of the object according to the MAST.
    """
    request = {'service' : 'Mast.Name.Lookup', 'params' :
               {'input' : str(name), 'format' : 'json'},}

    headers, outString = mastQuery(request)

    outData = json.loads(outString)

    #return outData['resolvedCoordinate'][0]['canonicalName']
    return outData

#function to convert semi-major axis into stellar radii units (or other units)
#  given stellar radius and what units semi-major axis is in
