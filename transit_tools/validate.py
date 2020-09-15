import numpy as np
import lightkurve as lk

import triceratops.triceratops as tr
from transit_tools.utils import tessobs_info

# display starfield and star table
#    once tpf/aperture downloaded, save as obj in main so
#       that it doesn't have to be reloaded each time
def get_starfield(tic, sectors=None, aperture=None, cadence='2min',
                  depth=None, target_out=False, show_plots=True):
    """
    Function to fetch and display the stars around a TESS target.

    !!Finish docstrings!!
    !!Allow for option to save plots and tables!!
    !!Allow for other cadences!!

    Parameters
    ----------
    tic : int
       TIC ID for the desired target source.
    sectors : array or list or None
       TESS sectors in which the target was observed. If None, all
       valid observed sectors will be collected and used.
    aperture : numpy array or None
       Array of aperture arrays with pixel coordinates of the form
       [col, row] for each pixel included in the aperture for each
       TESS sector. If None, the apertures will attempt to be 
    """
    if sectors is None:
        sectors = tessobs_info(tic=tic)['sector']
    
    target = tr.target(ID=tic, sectors=sectors)
    
    if aperture is None and cadence == '2min':
        aperture = []
        
        for i in range(len(sectors)):
            res = lk.search_targetpixelfile(('TIC' + str(tic)),
                                            mission='TESS',
                                            sector=int(sectors[i]))
            tpf = res.download(quality_bitmask='default')
            mask = []
            
            for j in range(len(tpf.pipeline_mask)):
                for k in range(len(tpf.pipeline_mask[0])):
                    if tpf.pipeline_mask[j][k]:
                        mask.append(
                            [(tpf.column + k), (tpf.row + j)]
                        )

            mask = np.array(mask)
            aperture.append(mask)

        aperture = np.array(aperture)

    if show_plots and aperture is not None and len(aperture) != 0:
        for i in range(len(sectors)):
            target.plot_field(sector=sectors[i],
                              ap_pixels=aperture[i])
    elif show_plots and (aperture is None or len(aperture) == 0):
        for i in range(len(sectors)):
            target.plot_field(sector=sectors[i])

    if depth is not None and aperture is not None and len(aperture) != 0:
        target.calc_depths(tdepth=depth, all_ap_pixels=aperture)

    print(target.stars)
    
    if target_out:
        return target

# calcfpp_tr
#    save as obj bc later func will multi this? (no multi
#       func, maybe save)
def calcfpp_tr(lc=None, *args, period=None, t0=None, depth=None,
               sectors=None, binsize=1, folded=True, target_in=None,
               target_out=False):
    """
    Function to calculate the FPP for a signal using TRICERATOPS.

    !!Make it so that single value can be passed for flux_err!!
    !!Add bin_num to specify number of bins across LC instead of 
      points per bin!!
    !!Remove NaNs!!

    Parameters
    ----------
    lc : object
       Object that contains a folded light curve with time in units
       of days from midtransit and flux normalized to 1. At 
       minimum, this must contain time and flux attributes, but may 
       contain flux_err as well.
    *args
       Time, flux, and optional flux_err arguments. Arguments passed
       must be arrays containing time from midtransit in days, flux
       normalized to 1, and optionally flux_err for each data point.
       lc must be None for args to be specified.
    period : float
       Orbital period of the signal in days. Not required if input 
       light curve is already folded.
    depth : float
       Depth of midtransit. Not required if target_in is not None.
    sectors : array or list
       Array or list of sectors in which this signal was observed.
       Not required if target_in is not None.
    binsize : int
       Number of points per bin for binning the folded light curve.
    folded : bool
       Flag indicating if the light curve input is folded or not.
    target_in : object
       Input target object if TRICERATOPS has been run previously or
       the get_starfield plot has been run. Default is None.
    target_out : bool
       Flag to indicate whether the target object used in this 
       function should be returned. If True, an additional output
       will be expected.

    Returns
    -------
    target : object
       Optional output containing the target object used in the 
       TRICERATOPS run in this function.
    """
    if lc is None:
        if len(args) == 1:
            raise ValueError('Please specify both time and flux')
        else:
            time = args[0]
            flux = args[1]
            if len(args) == 3:
                flux_err = flux_err

        

    #if not folded:
        #fold




    
# display table, plots for outcome

# full run (no, save for main, also somewhat included in above func)

# function to convert various aperture types (e.g. eleanor) to a
#    valid TRICERATOPS format

# all VESPA commands (LOTS of work/thinking needed for this)
