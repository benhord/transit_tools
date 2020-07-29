import numpy as np
from astroquery.mast import Observations
from lightkurve import TessLightCurveFile

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
            lc, sectors = get_2minlc(tic, sectors, out_sec=True, **kwargs)
        except:
            print('No TESS 2 minute light curves found! Trying FFIs...')
            method = 'ffi_ml'
            
    elif method == 'ffi_ml':
        try:
            lc, sectors = get_mlffi(tic, sectors, out_sec=True, **kwargs)
        except:
            print('No ML light curves found locally. Trying with eleanor...')
            method = 'eleanor'
            
    elif method == 'eleanor':
        try:
            lc, sectors = get_eleanor(tic, sectors, out_sec=True, **kwargs)
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
                
    obsTable=obsTable[good_ind] #returns only good indices for table
    
    #get lightcurve from MAST
    id = str(tic).zfill(16)
    for i in range(len(good_ind)):
        sec_str = "s" + str(obsTable['dataURL'][i][37:41])
        secs.append(int(obsTable['dataURL'][i][37:41]))
        lc_loc=("https://archive.stsci.edu/missions/tess/tid/" + str(sec_str) +
                "/" + id[:4] + "/" + id[4:8] + "/" + id[8:12] + "/" +
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

#eleanor
