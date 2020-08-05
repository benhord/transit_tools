import numpy as np

#global light curve processing function
def process_lc(lc, flatten_length=None, binsize=1, remove_outliers=None,
               remove_nans=True, mask_noise=False, **kwargs):
    """
    Function to process a light curve and return it after processing. The order
    in which processing occurs is removing NaNs, removing noise based on sector,
    binning, removing outliers manually, and finally flattening, which includes
    not only flattening using the scipy.signal.savgol_filter function, but also
    iterative outlier clipping. See 'lightkurve' documentation for further
    information on the capabilities of these functions.

    !!Update so processing occurs in order of keywords listed!!

    Parameters
    ----------
    lc : LightCurve object
       The light curve to be processed.
    flatten_length : odd int or None
       The window length to be passed to scipy.signal.savgol_filter for the
       flattening. If None, no flattening will be performed.
    binsize : int
       The number of points per bin that will be combined.
    remove_outliers : int or float or None
       If set to a number, it represents the number of sigma above which to 
       remove outliers from. If set to None, no outliers will be removed 
       separately from those already removed by the flattening.
    remove_nans : bool
       Flag to determine if NaNs will be removed from the light curve.
    mask_noise : bool
       Flag to determine whether parts of the light curve will be removed 
       based on how noisy those parts of each sector are on average and 
       depending on the type of light curve.
    kwargs
       Additional arguments to be passed to the 'flatten' method of the 
       LightCurve object.

    Returns
    -------
    lc : LightCurve object
       Processed light curve. Also now has the attribute lc.trend if flattening
       was performed. This attribute contains the trend of the light curve
       prior to flattening.
    """
    if remove_nans:
        lc = lc.remove_nans()
    
    if mask_noise:
        print('Feature in progress')
    
    if binsize:
        lc = lc.bin(binsize)

    if remove_outliers:
        lc = lc.remove_outliers(sigma=remove_outliers)
        
    if flatten_length:
        full_mask = np.zeros(len(lc.time), dtype=bool)
        
        if hasattr(lc, 'known_pls') and isinstance(lc.known_pls[0], dict):
            #mask based off of known planets when flattening
            
            for i in range(len(lc.known_pls)):
                period = lc.known_pls[i]['orbital_period']
                try:
                    t0 = lc.known_pls[i]['transit_time']
                    duration = lc.known_pls[i]['transit_duration']
                except:
                    print('Known planet %s not masked during flattening' %
                          lc.known_pls[i]['canonical_name'])
                    
                mask = np.ones(len(lc.time), dtype=bool)
                start = (lc.time[0] + (1 - ((lc.time[0] - t0) / period % 1)) *
                         period)
                for i in np.arange(o, len(lc.time) // period):
                    mask[(lc.time < (start + i * period + duration)) &
                         (lc.time > (start + i * period - duration))] = 0

                full_mask = [any(tup) for tup in zip(full_mask, mask)]
                
        lc, trend = lc.flatten(window_length=flatten_length, return_trend=True,
                               mask=full_mask, **kwargs)

        lc.trend = trend
    
    return lc
    

#light curve update function (could this just be self = lc in main.py?)
