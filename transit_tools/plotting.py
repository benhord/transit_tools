import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from transitleastsquares import transit_mask
import numpy as np
from lightkurve import LightCurve

import transit_tools.constants as c

##function to save all diagnostic plots as combined png

##function to generate vetting sheet
def tls_vetsheet(lc, results=0, show=True, save=False, savename='vetsheet.png'):
    """
    Function to plot a vetting sheet after running the 
    transit_tools.signal_search function. Output is formatted to fit onto a
    standard 8"x11" sheet of paper.

    !!Add functionality to plot without having to run search?!!
    !!Add vertical lines for expected transit times on unprocessed LC!!
    !!Change x-axis phase values for odd-even plot!!
    !!Make processed and unprocessed x-axes line up, especially when LC is 
      masked. Maybe include grayed-out parts that were masked!!

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
    if results == -1 and len(lc.results) > 0:
        time = lc.cleanlc[-1].time
        flux = lc.cleanlc[-1].flux
        flux_err = lc.cleanlc[-1].flux_err
    elif len(lc.results) > 0:
        res = lc.results[results]

    if results == -1:
        res = lc.bad_search[0]

    if results == 0 or len(lc.results) == 0:
        time = lc.time
        flux = lc.flux
        flux_err = lc.flux_err
    else:
        time = lc.cleanlc[results].time
        flux = lc.cleanlc[results].flux
        flux_err = lc.cleanlc[results].flux_err
    
    #setting up figure
    fig = plt.figure(figsize=(8, 9.2))
    gs = fig.add_gridspec(5, 2)

    #phase-folded light curve with transit model
    ax = fig.add_subplot(gs[0, 1])

    fold = LightCurve(res.folded_phase, res.folded_y).bin(20)
    
    ax.plot(res.model_folded_phase, res.model_folded_model, color='red',
            alpha=0.7)
    ax.scatter(res.folded_phase, res.folded_y, color='blue', s=0.2, alpha=0.5,
            zorder=2)
    ax.scatter(fold.time, fold.flux, s=2., c='k')
    
    ax.set_xlim(0.45, 0.55)
    ax.set_xlabel('Phase')
    ax.set_ylabel('Relative Flux')

    #raw light curve
    ax = fig.add_subplot(gs[1, :])

    ax.scatter(lc.raw_lc.time, lc.raw_lc.flux, s=0.2, c='b', alpha=0.5)
    
    ax.set_xlabel('Time [BTJD]')
    ax.set_ylabel('Raw Relative Flux')
    ax.set_xlim(lc.raw_lc.time.min(), lc.raw_lc.time.max())

    if hasattr(lc, 'trend'):
        plt.plot(lc.trend.time, lc.trend.flux, c='g')

    #processed light curve with transit model
    ax = fig.add_subplot(gs[2, :])

    transit_time = transit_mask(time, res.period, res.duration, res.T0)
    time_notrans = time[~transit_time]
    flux_notrans = flux[~transit_time]
    flux_err_notrans = flux_err[~transit_time]

    ax.scatter(time[transit_time], flux[transit_time], color='red', s=0.2,
               zorder=0)
    ax.scatter(time[~transit_time], flux[~transit_time], color='blue',
               alpha=0.5, s=0.2, zorder=0)
    ax.plot(res.model_lightcurve_time, res.model_lightcurve_model, alpha=0.5,
            color='red', zorder=1)
    ax.set_xlim(time.min(), time.max())
    ax.set_ylim(flux.min()-0.001, flux.max()+0.001)
    ax.set_xlabel('Time [BTJD]')
    ax.set_ylabel('Relative Flux')
    
    #TLS periodogram
    ax = fig.add_subplot(gs[3, :])
    
    ax.axvline(res.period, alpha=0.4, lw=3)
    ax.set_xlim(np.min(res.periods), np.max(res.periods))
    for n in range(2, 10):
        ax.axvline(n * res.period, alpha=0.4, lw=1, linestyle='dashed')
        ax.axvline(res.period / n, alpha=0.4, lw=1, linestyle='dashed')

    ax.set_ylabel('SDE')
    ax.set_xlabel('Period [days]')
    ax.plot(res.periods, res.power, color='k', lw=0.5)
    ax.set_xlim(0, max(res.periods))

    #secondary eclipse
    ax = fig.add_subplot(gs[4, 0])

    ax.plot(res.model_folded_phase-0.5,
            np.roll(res.model_folded_model, len(res.model_folded_model)//2),
            color='red')
    ax.scatter(res.folded_phase-0.5,
               np.roll(res.folded_y, len(res.folded_y)//2), color='blue',
               s=0.2, alpha=0.5, zorder=2)
    ax.set_xlim(-0.05, 0.05)
    ax.set_xlabel('Phase')
    ax.set_ylabel('Relative Flux')

    ax.axvspan(-(res.duration/2/res.period), (res.duration/2/res.period),
               alpha=0.3, color='orange', label='Transit Duration')

    ax.legend(loc=2, fontsize='x-small')
    
    #Odd-Even comparison
    ax = fig.add_subplot(gs[4, 1])

    oe = LightCurve(time, flux, flux_err)

    oe_odd = oe.fold(2*res.period, res.T0)
    oe_even = oe.fold(2*res.period, res.T0+res.period)

    ax.scatter(oe_odd.time+0.5, oe_odd.flux, s=0.2, c='b', label='Odd')
    ax.scatter(oe_even.time+0.5, oe_even.flux, s=0.2, c='r', label='Even')
    
    ax.set_xlim(0.475, 0.525)
    ax.set_xlabel('Phase')
    ax.set_ylabel('Relative Flux')
    ax.set_xticks([0.48, 0.49, 0.50, 0.51, 0.52])
    ax.set_xticklabels(['0.46', '0.48', '0.50', '0.52', '0.54'])
    ax.legend(loc=2, fontsize='x-small')
    
    #plot summary text
    fig.text(0.12, 0.97, s=(str(lc.name) + '  (TIC ' + str(lc.tic) + ')'),
             fontweight='bold')
    fig.text(0.04, 0.95,
             s=(r'P = %.5f +/- %.5f d, $t_{0}$ = %.5f BTJD' %
                (res.period, res.period_uncertainty, res.T0)))
    fig.text(0.04, 0.93,
             s=(r'$T_{dur}$ = %.5f d,  %s/%s transits with data' %
                (res.duration, res.distinct_transit_count, res.transit_count)))
    fig.text(0.04, 0.91,
             s=('SDE = %.2f,  SNR = %.2f,  FAP = %.3e' %
                (res.SDE, res.snr, res.FAP)))
    fig.text(0.04, 0.89,
             s=(r'$R_{P}$/$R_{*}$ = %.4f,  $R_{P}$ = %.3f $R_{\bigoplus}$ = %.3f $R_{Jup}$' %
                (res.rp_rs,
                 (lc.star_params_tls['rstar']*res.rp_rs*c.Rsolar_m)/c.Rearth_m,
                 (lc.star_params_tls['rstar']*res.rp_rs*c.Rsolar_m)/c.Rjup_m)))
    fig.text(0.04, 0.87,
             s=(r'odd/even mismatch = %.2f $\sigma$,  $\delta$ = %.4f' %
                (res.odd_even_mismatch, 1-res.depth)))
    fig.text(0.04, 0.85,
             s=(r'$R_{*}$ = %.2f (+%.2f, -%.2f) $R_{\bigodot}$'
                % (lc.star_params_tls['rstar'], lc.star_params_tls['rhigh'],
                   lc.star_params_tls['rlow'])))
    fig.text(0.04, 0.83, s=(r'$M_{*}$ = %.2f (+%.2f, -%.2f) $M_{\bigodot}$' %
                            (lc.star_params_tls['mstar'],
                             lc.star_params_tls['mhigh'],
                             lc.star_params_tls['mlow'])))
    if lc.star_params is not None and lc.star_params['Tmag'] is not None:
        fig.text(0.04, 0.81, s=('Tmag = %.2f' % (lc.star_params['Tmag'])))
    
    plt.tight_layout()
    
    if show:
        plt.show()

    if save:
        plt.savefig(savename)
    

##search-specific plots
