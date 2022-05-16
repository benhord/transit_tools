import numpy as np
from astroquery.mast import Observations, Tesscut
from astropy.coordinates import SkyCoord
from lightkurve import TessLightCurveFile, LightCurve
import lightkurve as lk
import pandas as pd
#import eleanor
import os
import pickle

from .utils import rms, tessobs_info, coord_to_tic

#misc to-do: add lc.cadence attribute for each method that gathers the cadnece
# from the header (or infers if header is unavailable) to be passed as part of
# lightkurve object (could be gathered by helper function later)

def gather_lc(coords=None, tic=None, name=None, cadence='shortest', ffi_only=False, method='2min', sectors='all', eleanor_flag=False,
              return_method=False,
              return_sectors=False, obsinfo=None, verbose=False, **kwargs):
    #!!Add ability to support common names as inputs in absence of TIC!!
    #!!Add ability for sector cuts in FFI light curves!!
    #!!Add obsinfo keyword to pass obsinfo if it exists!!
    #!!Add keyword to track 'author' and 'cadence' keywords for each sector?!!
    """
    Function to gather the light curve of a given TIC using the specified 
    method. Currently, 2 minute SPOC pipeline light curves, machine learning FFI
    light curves, and eleanor light curves are supported.

    Parameters
    ----------
    coords : tuple or array
       Coords of input
    name : string
       common name
    ffi_only : bool
       whether to default to 2min or have everything be ffis
    tic : int
       TESS Input Catalog ID for desired target. At this time, common names are
       not accepted input, only TIC IDs.
    method : str
       The method with which the light curve will be acquired. Options are 
       '2min', 'ffi_ml', and 'eleanor'.
    sectors : str or list or numpy array
       List of sectors to be included in the fetching of the light curve. If
       'all' or None is passed, all available light curves will be fetched. 
       Thresholds can be passed according to the valid syntax of the method
       specified.
    return_method : bool
       A flag to indicate whether the method used to gather the light curve 
       will be returned. Useful if the chosen method is not known or expected to
       return a valid light curve. An additional output will be expected.
    return_sectors : bool
       A flag to indicate the sectors that the light curve was recovered from.
       An additional output will be expected. Currently in progress for most
       methods.
    obsinfo : dict
       A dictionary of observational information that can be passed to make some 
       processes run faster if inlc is specified. Assumes output format of
       transit_tools.tessobs_info command.
    **kwargs
       Additional arguments to be passed to the selected method fetching 
       function.

    Returns
    -------
    lc : 'LightCurve'
       Light curve containing light curves from all sectors contained within the
       query of sectors.
    method : str, optional
       The method with which the output light curve was gathered.
    sectors : numpy array, optional
       The sectors that the light curve was gathered from.
    """
    if coords is None and tic is None and name is None:
        raise ValueError("You must specify either (RA, Dec), a TIC ID, or " +
                         "the name of the target!")

    #Fetching the TESS sectors that the target was observed in
    if coords is not None:
        coordinates = SkyCoord(coords[0], coords[1], unit="deg")
        sector_table = Tesscut.get_sectors(coordinates=coordinates)
    else:
        if tic is not None:
            name = "TIC " + str(tic)
            
        sector_table = Tesscut.get_sectors(objectname=str(name))

    if len(sector_table) == 0:
        raise ValueError('No valid sectors found for specified target!')
        
    #Comprehension of sector inputs
    if sectors != 'all' and sectors != None:
        if not isinstance(sectors, list):
            sectors = list(sectors)
        
        secs = list(set(sector_table['sector'].value) & set(sectors))
    else:
        secs = sector_table['sector'].value

    sec_clipboard = secs

    #Trying to fetch light curves from MAST
    if not eleanor_flag: #and not custom_only?
        try:
            lc, mast_sec = get_mastlc(tic=tic, name=name, coords=coords,
                                      sectors=sec_clipboard, out_sec=True,
                                      cadence=cadence, **kwargs)
            sec_clipboard = [x for x in sec_clipboard if x not in mast_sec]
        except:
            print('No light curves found on the MAST. Trying another method...')

    #Trying to generate light curves from eleanor

    ####
    #ELEANOR FUNCITON NEEDS TO BE DOUBLE-CHECKED!
    ####
    if eleanor_flag == True and len(sec_clipboard) > 0:
        try:
            lc, el_sec = get_eleanor(tic=tic, sectors=secs, out_sec=True,
                                     **kwargs)
            sec_clipboard = [x for x in sec_clipboard if x not in el_sec]
        except:
            print('eleanor could not find light curves for these sectors.' +
                  ' Trying another method...')

    #Trying to generate light curves via FFI cutouts
    if len(sec_clipboard) > 0:
        try:
            lc, get_ffilc = ffi_cutout()
            sec_clipboard = [x for x in sec_clipboard if x not in cut_sec]
        except:
            print('Issue with fetching FFI cutouts!')

    #Printing the sectors (if any) that did not have light curves
    if len(sec_clipboard) > 0:
        print('Sectors ' + str(sec_clipboard) + ' did not have light curves!!')



            
    #Loop through sectors to fetch light curves
    #if ffi_only is not None:
    #    try:
    #        lc = get_mastlc()
            #if return_sectors:
            #    lc, sectors = get_2minlc(tic, secs, out_sec=return_sectors,
            #                             **kwargs)
            #else:
            #    lc = get_2minlc(tic, secs, out_sec=return_sectors, **kwargs)
    #    except:
    #        print('No TESS 2 minute light curves found! Trying FFIs...')
    #        method = 'ffi_ml'
            
    #if method == 'ffi_ml':
    #    try:
    #        if return_sectors:
    #            lc, sectors = get_mlffi(tic=tic, sectors=secs,
    #                                    out_sec=return_sectors, **kwargs)
    #        else:
    #            lc = get_mlffi(tic=tic, sectors=secs, out_sec=return_sectors,
    #                           **kwargs)
    #    except:
    #        print('No ML light curves found locally. Trying with eleanor...')
    #        method = 'eleanor'
        
    if method == 'eleanor':
        try:
            if return_sectors:
                lc, sectors = get_eleanor(tic=tic, sectors=secs,
                                          out_sec=return_sectors, **kwargs)
            else:
                lc = get_eleanor(tic=tic, sectors=secs,
                                 out_sec=return_sectors, **kwargs)
        except:
            raise ValueError('No light curves found for the specified sectors!')

    if return_method and not return_sectors:
        return lc, method
    elif not return_method and return_sectors:
        return lc, secs
    elif return_method and return_sectors:
        return lc, method, secs
    elif not return_method and not return_sectors:
        return lc
            
def get_mastlc(name=None, coords=None, tic=None, sectors='all',
               author=['SPOC', 'TESS-SPOC'], cadence='shortest', out_sec=False,
               strict_cadence=False):
    #!!!Need to do some comprehension if multiple authors are specified so
    #   that their list order also denotes their priority order!!!
    #!!!mix_cadence keyword!!!
    #!!!ffi_only or no_ffi flags!!
    #!!!Add ability to bin sectors to a given cadence!!!
    """
    Function to retrieve 2 minute cadence TESS lightcurve for a given TIC ID 
    and given sectors. Returns a combined lightcurve in the form of a
    lightkurve object. If light curves from multiple sectors are combined, each
    light curve is individually normalized prior to combining.

    Parameters
    ----------
    name : str
       Name of the target requested.
    coords : tuple, 2-element list, or SkyCoord coordinates
       RA and Dec of the target being queried. Can be in decimal or sexigesimal.
    tic : integer
       TESS Input Catalog ID for desired target. At this time, common names are
       not accepted input, only TIC IDs.
    sectors : list, numpy array, or 'all'
       List of desired sectors to include when fetching the SPOC-processed light
       curve. If 'all' is specified, all available light curves
       will be downloaded.
    author : str or list
       Specifies the author of the light curves. Useful for retrieving light
       curves created by non-SPOC entities that have been posted on MAST.
       Default ['SPOC', 'TESS-SPOC'].
    cadence : str or int
       Specifies the desired cadence of the light curves. Options are 'shortest'
       for the shortest cadence in each sector, '2min' or 120 for 2 minute
       cadence, and '20sec' or 20 for 20 second cadence. Default is 'shortest'.
    out_sec : bool
       A flag to determine whether the sectors from which light curves were
       downloaded are included as an output. If True, command will provide two
       outputs, the light curve object and a numpy array of sectors, in that 
       order.
    strict_cadence : bool
       A flag to determine whether or not to be strict with cadence 
       requirements. If set to True, only the specified cadence(s) will be
       considered and light curves will not be fetched with other cadences if
       the specified ones are not available.

    Returns
    -------
    lc : 'LightCurve'
       Combined light curve of all available light curves at TESS 2 minute 
       cadence for specified TIC ID.
    secs : numpy array, optional
       List of sectors from which light curve was gathered.
    """
    if coords is None and tic is None and name is None:
        raise ValueError("Valid coordinates, a TIC ID, or a name must be " +
                         "specified!")

    if cadence == '2min' or cadence == 120 or cadence == '2minutes':
        exptime = 120 #add def_authors here for short or long cadence flags?
        #add flag for only short cadence and no ffis
    elif cadence == '20sec' or cadence == 20 or cadence == '20seconds':
        exptime = 20
    elif cadence == '30min' or cadence == 1800 or cadence == '30minutes':
        exptime = 1800
    elif cadence == '10min' or cadence == 600 or cadence == '10minutes':
        exptime = 600

    if tic is not None:
        name = "TIC " + str(tic)

    if name is not None:
        search_result = lk.search_lightcurve(name, author=author)
    else:
        search_result = lk.search_lightcurve((coords[0]+" "+coords[1]),
                                             author=author)

    if len(search_result) == 0:
        raise ValueError('No valid sectors found for target with specified ' +
                         'inputs!')

    secs = list(map(int, [row[12:] for row in search_result.table['mission']]))

    #remove unwanted sectors
    if sectors != 'all':
        if not isinstance(sectors, list):
            sectors = list(sectors)
        search_result = search_result[:][[(element in sectors) for element in
                                          secs]]
        if len(search_result) == 0:
            raise ValueError('None of the specified sectors found for the ' +
                             'specified target!')
    
    #remove duplicate sectors based on shortest cadence or specify exptime
    if cadence != 'shortest' and strict_cadence is True:
        search_result = search_result[search_result.table['exptime'] == exptime]

    df = pd.DataFrame([[int(row['mission'][12:]), int(row['exptime'])] for
                       row in search_result.table],
                      columns=['sector', 'exptime'])

    keep = np.zeros(len(df['sector']), dtype=bool)
    keep[df.loc[df.groupby('sector').exptime.idxmin()].index] = True

    search_result = search_result[keep]

    if len(search_result) == 0:
        raise ValueError('No light curves found at the specified cadence!')

    secs = list(map(int, [row[12:] for row in search_result.table['mission']]))
    
    #download light curves
    lc_col = search_result.download_all()
    lc = lc_col.stitch()

    if out_sec:
        return lc, secs
    else:
        return lc

#search for ML FFIs
def get_mlffi(tic=None, ra=None, dec=None, sectors='all',
              flux_type='corr', out_sec=False):
    """
    For use on tesseract only. Fetches FFI light curves made by Brian Powell.

    Parameters
    ----------
    tic : int or None
       TIC ID of target. If None, both ra and dec must not be None.
    ra : float or None
       RA of target. If None, tic must not be None. If not None, dec must also
       not be None.
    dec : float or None
       Dec of target. If None, tic must not be None. If not None, ra must also 
       not be None.
    sectors : str or list or array
       The desired sectors to make the light curve from. May be set to 'all' to
       use all available light curves.
    flux_type : str
       The type of correction applied to the light curve. Options are 'raw', 
       'corr', and 'pca'.
    out_sec : bool
       Flag determining whether or not the sectors that the light curve was
       generated from are output. If True, an additional output will be 
       expected.

    Returns
    -------
    lc : 'LightCurve' object
       The combined light curve.
    sectors : array
       The sectors from which the output light curve was generated.
    """
    if not tic and not ra and not dec:
        raise ValueError('Please provide valid input for either tic or ra/dec')

    if sectors == 'all':
        if not ra and not dec:
            info = tessobs_info(tic=tic)
        else:
            info = tessobs_info(ra=ra, dec=dec)

        sectors = list(set(info['sector']))
        
    if not isinstance(sectors, list):
        raise ValueError('Sectors must be a list!')

    if not tic:
        tic = coord_to_tic(ra, dec)

    camera_arr = []
    chip_arr = []
    Tmag_arr = []

    sects = sectors

    if not os.path.isdir('/data/tessraid/bppowel1/'):
        raise ValueError('Not on tesseract')
    
    for i in range(len(sects)):
        path = '/data/tessraid/bppowel1/tesslcs_sector_'+str(sects[i])+'_104'
        lc_files = []
        
        for (dirpath, dirnames, filenames) in os.walk(path):
            for name in filenames:
                lc_files.append(os.path.join(dirpath, name))

        lc_files = [t for t in lc_files if 'tesslc_' in t]
        tics = [int(t.split('.')[0].split('_')[-1]) for t in lc_files]
        
        path = [s for s in lc_files if str(tic) in s]
        
        try:
            fp = open(str(path[0]), 'rb')
        except:
            print('Target not found in Sector %s' % sects[i])
            sects.remove(sects[i]) #trouble line, may cause silent error
            continue
        
        data = pickle.load(fp)
        fp.close()

        Tmag_arr.append(data[2])
        camera_arr.append(data[4])
        chip_arr.append(data[5])
        
        time = data[6]
        flux_err = data[10]

        if flux_type == 'corr':
            flux = data[8]
        elif flux_type == 'raw':
            flux = data[7]
        elif flux_type == 'pca':
            flux = data[9]

        if i == 0:
            lc = LightCurve(time, flux, flux_err=flux_err)
        else:
            sec_lc = LightCurve(time, flux, flux_err=flux_err)
            lc.append(sec_lc)
            
    lc = lc.normalize()
    
    lc.Tmag = Tmag_arr
    lc.camera = camera_arr
    lc.chip = chip_arr
    
    if out_sec:
        return lc, sects
    else:
        return lc

#eleanor
def get_eleanor(sectors='all', tic=None, coords=None, out_sec=False, height=15,
                width=15, bkg_size=31, do_psf=False, do_pca=False,
                out_flux='corr_flux', norm=True, errorcalc=True,
                qual_flag=True):
    #!!Add more docstrings for all keywords!!
    #!!Add common name processing instead of just tic!!
    """
    Function to get a light curve from the TESS full frame images (FFIs) using
    the Python package eleanor.

    Parameters
    ----------
    sectors : str or array or list
       The sectors that eleanor will use to produce the light curve. If set to 
       'all', then all available sectors will be used in the light curve
       production.
    tic : int or None
       TIC ID for the object that a light curve is desired for. If set to None,
       coords must have a valid input.
    coords : tuple of floats
       The RA and Dec of the object that a light curve is desired for. Must be
       of the form (RA, Dec) in decimal degrees. If set to None, the tic 
       argument cannot be None.
    out_sec : bool
       Flag controlling whether an array containing the sectors used to extract
       the light curve will be output. If True, an additional output will be 
       expected.
    height : int
       Height in pixels of the postage stamp with which to extract the light 
       curve.
    width : int
       Height in pixels of the postage stamp with which to extract the light
       curve.
    bkg_size : int
       Background size to be considered for the background subtraction from the
       light curve.
    do_psf : bool
       Flag to determine whether a PSF-corrected light curve will be generated
       as an additional option to the corrected light curve.
    do_pca : bool
       Flag to deteremine whether a PCA-corrected light curve will be generated
       as an additional option to the corrected light curve.
    out_flux : str
       Which of the light curves to output. Options are 'corr_flux', 'psf_flux',
       and 'pca_flux'. Only one may be selected. If either 'psf_flux' or 
       'pca_flux' are selected, the do_psf and do_pca flags must be set to True,
       respectively.
    norm : bool
       Flag determining whether the light curve will be normalized prior to 
       output.
    errorcalc : bool
       Flag determining whether the RMS errors will be calculated for the light
       curve.
    qual_flag : bool
       Flag determining whether the timestamps with bad quality flags will be 
       excluded automatically.

    Returns
    -------
    lc : 'LightCurve' object
       The combined light curve from each sector for the coordinates or TIC ID 
       requested.
    sectors : array
       Optional output array containing the sectors that the combined light 
       curve was extracted from.
    """
    import eleanor
    
    if tic is None and coords is None:
        raise ValueError('Please make sure either tic or coords have valid ' +
                         'input')
    
    if tic:
        star = eleanor.multi_sectors(tic=tic, sectors=sectors)
    else:
        star = eleanor.multi_sectors(coords=coords, sectors=sectors)
    
    secs = []
    data = []
    
    for s in star:
        datum = eleanor.TargetData(s, height=height, width=width,
                                   bkg_size=bkg_size, do_psf=do_psf,
                                   do_pca=do_pca)
        data.append(datum)
        
        sec = s.sector
        secs.append(sec)

    for i in range(len(data)):
        q = data[i].quality == 0
        time = data[i].time
            
        if out_flux == 'corr_flux':
            flux = data[i].corr_flux
        elif out_flux == 'pca_flux':
            flux = data[i].pca_flux
        elif out_flux == 'psf_flux':
            flux = data[i].psf_flux

        if qual_flag:
            time = time[q]
            flux = flux[q]
            
        if norm:
            flux = flux/np.median(flux)

        flux_err = None
        if errorcalc:
            flux_err = np.ones(len(flux)) * rms(flux)
            
        if i == 0:
            lc = LightCurve(time, flux, flux_err=flux_err)
        else:
            sec_lc = LightCurve(time, flux, flux_err=flux_err)
            lc = lc.append(lc)

    if out_sec:
        return lc, secs
    else:
        return lc

def get_ffilc(name=None, tic=None, coords=None, sectors='all'):
    """
    """
    if name is None and tic is None and coords is None:
        raise ValueError('Please specify a name, TIC ID, or coordinates!')

    if name is None and tic is not None:
        name = 'TIC ' + str(tic)

    if coords is not None:
        if isinstance(coords, tuple) or isinstance(coords, list):
            coords = SkyCoord(coords[0], coords[1], unit='deg')
        
    if sectors == 'all' and coords is not None:
        sector_table = Tesscut.get_sectors(coordinates=coord)
        sectors = list(map(int, [row[6:10] for row in
                                 sector_table['sectorName']]))
        #tesscut add more here to parse sectors?
    elif sectors == 'all' and name is not None:
        sector_table = Tesscut.get_sectors(objectname=name)
        sectors = list(map(int, [row[6:10] for row in
                                 sector_table['sectorName']]))



        
    
    if name is None and tic is not None:
        name = "TIC " + str(tic)
        search_result = lk.search_tesscut(name, sector=sectors)
    else:
        print('hi')
