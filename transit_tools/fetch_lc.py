import numpy as np
from astroquery.mast import Observations
from lightkurve import TessLightCurveFile, LightCurve
import eleanor
import os
import pickle

from .utils import rms, tessobs_info, coord_to_tic

#misc to-do: add lc.cadence attribute for each method that gathers the cadnece
# from the header (or infers if header is unavailable) to be passed as part of
# lightkurve object (could be gathered by helper function later)

def gather_lc(tic, method='2min', sectors='all', return_method=False,
              return_sectors=False, **kwargs):
    """
    Function to gather the light curve of a given TIC using the specified 
    method. Currently, 2 minute SPOC pipeline light curves, machine learning FFI
    light curves, and eleanor light curves are supported.

    !!Add ability to support common names as inputs in absence of TIC!!

    Parameters
    ----------
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
    if sectors == None: sectors = 'all'
    
    if method == '2min':
        try:
            if return_sectors:
                lc, sectors = get_2minlc(tic, sectors, out_sec=return_sectors,
                                         **kwargs)
            else:
                lc = get_2minlc(tic, sectors, out_sec=return_sectors, **kwargs)
        except:
            print('No TESS 2 minute light curves found! Trying FFIs...')
            method = 'ffi_ml'
            
    if method == 'ffi_ml':
        try:
            if return_sectors:
                lc, sectors = get_mlffi(tic, sectors, out_sec=return_sectors,
                                        **kwargs)
            else:
                lc = get_mlffi(tic, sectors, out_sec=return_sectors, **kwargs)
        except:
            print('No ML light curves found locally. Trying with eleanor...')
            method = 'eleanor'
            
    if method == 'eleanor':
        try:
            if return_sectors:
                lc, sectors = get_eleanor(tic=tic, sectors=sectors,
                                          out_sec=return_sectors, **kwargs)
            else:
                lc = get_eleanor(tic=tic, sectors=sectors,
                                 out_sec=return_sectors, **kwargs)
        except:
            raise ValueError('No light curves found for the specified sectors!')

    if return_method and not return_sectors:
        return lc, method
    elif not return_method and return_sectors:
        return lc, sectors
    elif return_method and return_sectors:
        return lc, method, sectors
    elif not return_method and not return_sectors:
        return lc
            
def get_2minlc(tic, sectors='all', thresh=None, out_sec=False):
    """
    Function to retrieve 2 minute cadence TESS lightcurve for a given TIC ID 
    and given sectors. Returns a combined lightcurve in the form of a
    lightkurve object. If light curves from multiple sectors are combined, each
    light curve is individually normalized prior to combining.

    Parameters
    ----------
    tic : integer
       TESS Input Catalog ID for desired target. At this time, common names are
       not accepted input, only TIC IDs.
    sectors : list, numpy array or 'all'
       List of desired sectors to include when fetching the SPOC-processed light
       curve. If 'all' is specified, all available 2 minute PDCSAP light curves
       will be downloaded.
    thresh : string of form AA##
       Threshold to specify a range of sectors without knowing the specific 
       sectors that contain the target. 'AA' should be either 'gt' or 'lt' for 
       'greater than' and 'less than', respectively. ## is the threshold sector.
       The sector number specified in the threshold will not be included in the
       sector query. EX: 'lt13' will return all 2 minute light curves of the 
       target from sectors prior to, but not including, sector 13.
    out_sec : boolean
       A flag to determine whether the sectors from which light curves were
       downloaded are included as an output. If True, command will provide two
       outputs, the light curve object and a numpy array of sectors, in that 
       order.

    Returns
    -------
    lc : 'LightCurve'
       Combined light curve of all available light curves at TESS 2 minute 
       cadence for specified TIC ID.
    secs : numpy array, optional
       List of sectors from which light curve was gathered.
    """

    obsTable = Observations.query_criteria(dataproduct_type=['timeseries'], 
                                           target_name=tic,
                                           obs_collection='TESS')

    lc_str = "lc.fits"
    good_ind = []
    secs = []

    if thresh: #get sectors from gt or lt range
        for i in range(len(obsTable['dataURL'])):
            if (thresh[:2] == 'lt' and
                str(obsTable['dataURL'][i]).find(lc_str) > 0 and
                int(obsTable['dataURL'][i][37:41]) < int(thresh[2:])):
                good_ind.append(i)
            if (thresh[:2] == 'gt' and
                str(obsTable['dataURL'][i]).find(lc_str) > 0 and
                int(obsTable['dataURL'][i][37:41]) > int(thresh[2:])):
                good_ind.append(i)

    elif sectors != 'all': #get just specified sectors
        for sec in sectors:
            sec_str = "s" + str(sec).zfill(4)
            for i in range(len(obsTable['dataURL'])):
                if (str(obsTable['dataURL'][i]).find(lc_str) > 0 and
                    str(obsTable['dataURL'][i]).find(sec_str) > 0):
                    good_ind.append(i)
                    
    elif sectors == 'all': #get all sectors from list
        for i in range(len(obsTable['dataURL'])):
            if str(obsTable['dataURL'][i]).find(lc_str) > 0:
                good_ind.append(i)

    if len(good_ind) == 0:
        #print("WARNING: NO VALID SECTORS IN GIVEN RANGE!!")
        raise ValueError("ERROR: No valid sectors in given range!! " +
                         "Try running the command again with a different " +
                         "sector range.")
                
    obsTable = obsTable[good_ind] #returns only good indices for table
    
    #get lightcurve from MAST
    id = str(tic).zfill(16)
    for i in range(len(good_ind)):
        sec_str = "s" + str(obsTable['dataURL'][i][37:41])
        secs.append(int(obsTable['dataURL'][i][37:41]))
        lc_loc = ("https://archive.stsci.edu/missions/tess/tid/" + str(sec_str)
                  + "/" + id[:4] + "/" + id[4:8] + "/" + id[8:12] + "/" +
                  id[12:16] + "/" + str(obsTable['dataURL'][i])[18:])

        if i == 0:
            lc = (TessLightCurveFile(lc_loc).PDCSAP_FLUX.normalize().
                  remove_nans())
            
        else:
            lc_new = (TessLightCurveFile(lc_loc).PDCSAP_FLUX.normalize().
                      remove_nans())
            lc = lc.append(lc_new)

    if out_sec:
        return lc, np.array(secs)
    else:
        return lc

#search for ML FFIs
def get_mlffi(tic=None, ra=None, dec=None, sectors='all',
              flux_type='corr_flux', out_sec=False):
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
    
    for i in range(len(sectors)):
        path = '/data/tessraid/bppowel1/tesslcs_sector_'+str(sectors[i])+'_104'
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
            print('Target not found in Sector %s' % sectors[i])
            sectors.remove(sectors[i])
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

    lc.Tmag = Tmag_arr
    lc.camera = camera_arr
    lc.chip = chip_arr

    if out_sec:
        return lc, sectors
    else:
        return lc

#eleanor
def get_eleanor(sectors='all', tic=None, coords=None, out_sec=False, height=15,
                width=15, bkg_size=31, do_psf=False, do_pca=False,
                out_flux='corr_flux', norm=True, errorcalc=False):
    """
    Function to get a light curve from the TESS full frame images (FFIs) using
    the Python package eleanor.

    !!Add more docustrings for all keywords!!
    !!Add common name processing instead of just tic!!

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

    Returns
    -------
    lc : 'LightCurve' object
       The combined light curve from each sector for the coordinates or TIC ID 
       requested.
    sectors : array
       Optional output array containing the sectors that the combined light 
       curve was extracted from.
    """
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
        time = data[i].time[q]

        if out_flux == 'corr_flux':
            flux = data[i].corr_flux[q]
        elif out_flux == 'pca_flux':
            flux = data[i].pca_flux[q]
        elif out_flux == 'psf_flux':
            flux = data[i].psf_flux[q]

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
