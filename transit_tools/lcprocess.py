import numpy as np

#global light curve processing function
def process_lc(lc, flatten_length=None, binsize=1, remove_outliers=None,
               remove_nans=True, mask_noise=False, mask_plinds=None, **kwargs):
    """
    Function to process a light curve and return it after processing. The order
    in which processing occurs is removing NaNs, removing noise based on sector,
    binning, removing outliers manually, and finally flattening, which includes
    not only flattening using the scipy.signal.savgol_filter function, but also
    iterative outlier clipping. See 'lightkurve' documentation for further
    information on the capabilities of these functions.

    !!Update so processing occurs in order of keywords listed!!
    !!Update to let user specify which planets from lc.known_pls will be 
         included in mask for flattening (maybe by passing arr of indices?)!!

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
    mask_plinds : list or array
       Series of indices in known_pls property of the light curve that will be
       masked out before flattening. If None, all planets in known_pls will be
       masked out before flattening.
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
        lc = lc.bin(time_bin_size=binsize)

    if remove_outliers:
        lc = lc.remove_outliers(sigma=remove_outliers)
        
    if flatten_length:
        
        #if len(lc.known_pls) > 0 and isinstance(lc.known_pls[0], dict):
        if lc.known_pls is not None and isinstance(lc.known_pls[0], dict):
            #mask based off of known planets when flattening
            full_mask = np.zeros(len(lc.time), dtype=bool)

            if mask_plinds is None:
                mask_plinds = np.arange(0, len(lc.known_pls), 1)
            
            for i in mask_plinds:
                period = lc.known_pls[i]['pl_orbper'].value
                try:
                    t0 = lc.known_pls[i]['t0']# - 56999.5 #BTJD
                    duration = lc.known_pls[i]['pl_trandur'].value
                except:
                    print('Known planet %s not masked during flattening' %
                          lc.known_pls[i]['pl_name'])
                    
                mask = np.ones(len(lc.time), dtype=bool)
                start = (lc.time[0] + (1 - ((lc.time[0] - t0) / period % 1)) *
                         period)
                for i in np.arange(0, len(lc.time) // period):
                    mask[(lc.time < (start + i * period + duration)) &
                         (lc.time > (start + i * period - duration))] = 0

                full_mask = np.array([any(tup) for tup in zip(full_mask, mask)])

        else:
            full_mask = np.ones(len(lc.time), dtype=bool)
                
        lc, trend = lc.flatten(window_length=flatten_length, return_trend=True,
                               mask=~full_mask, **kwargs)

        lc.trend = trend
    
    return lc

#median absolute deviation (MAD) function

#mask function
def mask(lc, upper=None, lower=None, out_mask=True):
    """
    Function to mask out sections of the user-provided light curve.

    !!Add interpolation option to fill in gaps with noise!!
    !!Add option to input a mask that can then be altered with further
      removals from light curve!!
    
    Parameters
    ----------
    lc : 'LightCurve' object
       The light curve to be masked.
    upper : float
       Upper limit timestamp of the masked region. If lower is None, all 
       timestamps before this value will be masked out. This value is not
       inclusive and will remain in the light curve.
    lower : float
       Lower limit timestamp of the masked region. If upper is None, all 
       timestamps after this value will be masked out. This value is not 
       inclusive and will remain in the light curve.
    out_mask : bool
       Flag to determine whether the mask generated will be output. If True, an
       additional output will be expected.

    Returns
    -------
    masked_lc : 'LightCurve' object
       The light curve with the specified portions masked out.
    mask : array, optional
       A boolean array of identical length to the time array from the input 
       light curve with masked portions set to False and everything else set to 
       True.
    """
    if not isinstance(lc, object):
        raise ValueError('Please make sure that the input light curve is a ' +
                         'LightCurve object')
        
    time = lc.time

    if upper is not None and lower is not None:
        mask = ((time < lower) | (time > upper))
    
    elif upper is not None and lower is None:
        mask = (time > upper)
    elif upper is None and lower is not None:
        mask = (time < lower)

    masked_lc = lc[mask]
        
    if out_mask:
        return masked_lc, mask
    else:
        return lc

def gen_transitmask(lc, period, t0, duration):
    """
    Generates boolean array for masking out periodic signals in a light curve.

    !!Need to allow user to pass in mask and add to it!!
    !!Need to allow for period and t0 to be list and iterate through!!
    """
    mask = np.zeros_like(lc.time, dtype=bool)
    
    for i in range(int((max(lc.time.value) - min(lc.time.value)) // period + 1)):
        mask[(lc.time.value > (t0 + (period * i) - (duration / 2))) & (lc.time.value < (t0 + (period * i) + (duration / 2)))] = True

    return mask
