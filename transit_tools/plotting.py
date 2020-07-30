import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

##function to save all diagnostic plots as combined png

##function to generate vetting sheet
def tls_vetsheet(lc, results=0, show=True, save=False, savename=None):
    """
    Function to plot a vetting sheet after running the 
    transit_tools.signal_search function. Output is formatted to fit onto a
    standard 8"x11" sheet of paper.

    !!Add functionality to plot without having to run search?!!

    Parameters
    ----------
    lc : 'lightcurve' object
       Input transit_tools 'lightcurve' object that has had a signal search
       performed on it.
    results : int
       Index of results attribute array of provided 'lightcurve' object. 
       Indicates which set of results to plot if signal_search produced more
       than one set of output results. Can be set to -1 to display the most
       recent signal run that did not meet the significance threshold.
    show : bool or str
       Flag to determine whether plots will be displayed or not. Must be set to
       False or 'both' for output matplotlib object to be expected.
    save : bool
       Flag to determine whether the plots will be saved as a PNG.
    savename : str or None
       File name for plots to be saved as if save is set to True.

    Returns
    -------
    plots : matplotlib object
       Output matplotlib plot object. Optional if show is set to False or 
       'both'.
    """
    if results == -1:
        res = lc.bad_search[0]
    else:
        res = lc.results[results]
    
    #setting up figure
    fig = plt.figure(figsize=(8, 10))
    gs = fig.add_gridspec(5, 2)

    #phase-folded light curve with transit model
    ax = fig.add_subplot(gs[0, 1])


    #raw light curve
    ax = fig.add_subplot(gs[1, :])


    #processed light curve with transit model
    ax = fig.add_subplot(gs[2, :])


    #TLS periodogram
    ax = fig.add_subplot(gs[3, :])


    #secondary eclipse
    ax = fig.add_subplot(gs[4, 0])


    #Odd-Even comparison
    ax = fig.add_subplot(gs[4, 1])
    
    

    if show:
        plt.show()

    if show == 'both':
        return plots
    

##search-specific plots
