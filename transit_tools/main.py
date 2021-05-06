import numpy
import pandas as pd
from lightkurve.lightcurve import *
import matplotlib.pyplot as plt

from .fetch_lc import gather_lc
from .search import *
from .utils import *
from .plotting import *
from .lcprocess import *
from .batman import *
from .constants import *

class lightcurve(LightCurve):
    #pass mission-specific classes to preserve some speciality?
    
    """
    Initialization function for all of transit_tools that defines the light 
    curve to be used. This class wraps `lightkurve.LightCurve` and as such 
    retains all functionality of that class with some additions. The light 
    curve can either be user-defined with a `lightkurve.LightCurve` object 
    or arrays of time, flux, and flux error or can be generated via a target 
    name.

    Parameters
    ----------
    lc : `lightkurve.LightCurve` or None
        `lightkurve.LightCurve` object. Mutually exclusive with ``*args`` 
        and ``obj`` inputs. Default ``None``.
    *args : 
        Up to three arguments specifying the time, flux, and flux error of the 
        light curve. Mutually exclusive with ``lc`` and ``obj`` inputs.
    obj : str or None
        Name of target to search MAST for a light curve. Mutually exclusive 
        with ``lc`` and ``*args`` inputs. Default ``None``.
    method : str
        Method for acquiring a specified light curve through a search of MAST. 
        Options are currently ``'2min'``, ``'eleanor'``, and ``'batman'``. If 
        ``'batman'`` is provided, user will also need to pass ``**kwargs`` to 
        the ``transit_tools.batman`` function. Default is ``'2min'``.
    sector : list or None
        List of sectors to include in final assembled light curve. If ``None`` 
        is specified, all possible sectors will be included. Default ``None``.
    mission : str
        Current placeholder. Will eventually allow user to specify mission from 
        which light curve will be assembled. At this point, only TESS is 
        supported.
    find_knownpls : bool
        Flag to control whether a query to the NASA Exoplanet Archive is 
        searched for any known planets in the target system. Default ``True``.
    values : list or str
        Quantities to retrieve from the NASA Exoplanet Archive if 
        ``find_knownpls`` is True. If set to ``'all'``, all available 
        quantities will be retrieved. Default ``'all'``.
    **kwargs : dict
        Additional arguments to pass to `lighturve.LightCurve` initializer.
    
    Attributes
    ----------
    time : `np.array`
        The timestamps of the light curve.
    flux : `np.array`
        The flux at each timestamp in the light curve.
    flux_err : `np.array`
        The flux error for each flux in the light curve.
    method : str or None
        The method used to generate the light curve.
    sector : list
        The list of sectors the generated light curve represents.
    name : str or None
        The name of the target from the light curve query.
    tic : int
        TIC ID of target, if it has one.
    known_pls : `np.array` of dicts or None
        The known planets in the system, if any.
    star_params : dict
        A dictionary containing known star params if a query for such has been
        performed.
    star_params_tls : dict
        The star params used for Transit Least Squares searches, if such a 
        search has been run.
    raw_lc : object
        Nested `lightcurve` object with `time`, `flux`, and `flux_err` 
        attributes of the original input lightcurve prior to processing or 
        masking.
    trend : `np.array`
        Trend of light curve based on the `process` method, if detrending was
        performed.
    mask_arr : `np.array`
        Boolean array of equal length to the light curve mapping the masked data
        points.
    routine : str
        The shape used by the most recent Transit Least Squares search on the
        light curve.
    sde_thresh : float
        The Signal Detection Efficiency threshold for a significant detection
        specified in the most recent Transit Least Squares search.
    results : `np.array` of dicts
        Output results of each iteration of transit searches from the most 
        recent call of the `signal_search` method.
    cleanlc : list of objects
        Light curves with significant periodic signals from the most recent call
        of the `signal_search` method subtracted.
    bad_search : dict
        Results from the most recent call of `signal_search` that did not meet
        the specified `sde_thresh` value.
    """
    
    def __init__(self, lc=None, *args, obj=None, method="2min", sector=None,
                 mission='TESS', find_knownpls=True, values='all', **kwargs):
        #Have 2min, ffi_ml (which points to Brian's lcs), and eleanor.
        #This assumes you're looking for TESS lcs, can be generalized later.
        #Defaults to looking for 2min lc using get_2minlc script from
        #vespa_wrapper but falls back to ffi_ml then to eleanor in this order
        #if 2min is specified with some explanatory messages. Use property
        #decorators to control which sectors can be used and the cadence.

        #!!Make conversion to name robust for other missions!!!
        #!!Enable IDs other than TIC to be used and cross-referenced to find
        #  the correct ID for the requested mission!!
        #!!Enable passing custom light curves as separate method option!!
        #!!**kwargs to pass to gather_lc command!!
        #!!Enable RA/Dec search for light curves!!
        #!!If lc is not none or args are given, allow for obj to still search
        #  for known_pls and stellar info!!

        self.sector = sector
        self.method = method
        self.star_params = None
        
        if lc is not None or len(args) != 0:
            self.method = 'custom'
        
        if isinstance(obj, str):
            try:
                self.tic = name_to_tic(obj)
                self.name = str(obj)
            except:
                raise ValueError('For some reason, the name provided was not ' +
                                 'found on the MAST. Please try again.')

        elif isinstance(obj, int):
            try:
                self.tic = int(obj)
                self.name = tic_to_name(obj)
            except:
                print('For some reason, the name provided was not found on ' +
                      'the MAST. Proceeding with just the TIC.')
                self.name='Null'
        
        @property
        def sector(self):
            """'The TESS sector(s) being used'"""
            return self.sector

        @sector.setter
        def sector(self, sector):
            if not isinstance(sector, list) and (value != None):
                values = np.array(sector)
            self.sector = sector

        @property
        def method(self):
            """'The method for acquiring the TESS light curve'"""
            return self.method
        
        @method.setter
        def method(self, method):
            if not method in ['2min', 'ffi_ml', 'eleanor']:
                raise ValueError('Please specify a supported light curve type.')
            self.method = method

        @property
        def mission(self):
            "The mission that collected the desired light curve."
            return self.mission

        @mission.setter
        def mission(self, mission):
            if not mission in ['Tess', 'tess', 'TESS']:
                raise ValueError('Specified mission not currently supported!')
            self.mission = mission
            
        if self.method != 'custom' and self.method != 'batman':
            self.lc, self.method, self.sector = gather_lc(
                self.tic,
                method=str(self.method),
                sectors=self.sector,
                return_method=True,
                return_sectors=True,
                **kwargs
            )
            self.id = self.tic

        if self.method == 'batman' or self.method == 'BATMAN' or self.method == 'Batman':
            self.id = 'batman'
            self.name = self.id
            self.star_params_tls = default_star_params
            
            lc = full_batlc(sectors=self.sector, **kwargs)
            self.lc = lc
            self.tic = lc.tic
            if lc.tic is not None:
                cat = catalog_info(tic=lc.tic)
                self.known_pls = known_pls(ra=cat['ra'], dec=cat['dec'], values=values)
                self.known_pls.append(lc.known_pls)
                #allow for sim params to be passed in same form as real
            else: #add step in case inlc is not None but known_pls=None
                #print('inlc not known')
                self.known_pls = [lc.known_pls]
            self.sector = lc.sectors
            
            self.time = lc.time
            self.flux = lc.flux
            self.flux_err = lc.flux_err

            #update known pls w/ simulated params
            
        if self.method == 'custom':
            self.id = 'custom'
            self.tic = None
            self.name = self.id
            self.star_params_tls = default_star_params
            
            if lc is not None:
                self.lc = lc
            elif len(args) == 1:
                raise ValueError('Please make sure you provide both time ' +
                                 'and flux or an object name to query')
            elif len(args) > 1:
                self.time = args[0]
                self.flux = args[1]
                self.flux_err = None
                
                if len(args) == 3:
                    self.flux_err = args[2]
                    
        super().__init__(time=self.lc.time, flux=self.lc.flux,
                         flux_err=self.lc.flux_err)

        if find_knownpls and self.method != 'custom' and self.method !='batman':
            try:
                if self.name != 'Null':
                    self.known_pls = known_pls(name=self.name, values=values)
                else:
                    self.known_pls = known_pls(name=('TIC ' + str(self.tic)),
                                               values=values)
            except:
                self.known_pls = None
            
        
    ###Method to process light curve
    def process(self, **kwargs):
        """
        Method to process light curve for further analysis. Must be performed on
        a `lightkurve.LightCurve` object.

        Parameters
        ----------
        **kwargs : dict
            Arguments to be passed to the ``process_lc``.
        """
        lc = process_lc(self, **kwargs)

        self._updatelc(lc)

    ###method to update lc
    def _updatelc(self, lc):
        """
        Method to update the light curve that will be used if manual
        alteration is performed that the user wishes to save without completely 
        overriding the original light curve.

        Parameters
        ----------
        lc : `lightkurve.LightCurve` or `transit_tools.lightcurve`
            The light curve that will become the main, updated light curve.
        """
        #check for self.flux_err attribute
        if not hasattr(self, 'raw_lc'):
            raw_lc = LightCurve(self.time, self.flux, flux_err=self.flux_err)
            self.raw_lc = raw_lc

        self.time = lc.time
        self.flux = lc.flux
        #check if lc.flux_err exists
        self.flux_err = lc.flux_err

        if hasattr(lc, 'trend'):
            self.trend = lc.trend

    def reset(self):
        """
        Method to bring the raw, original light curve back as the main, 
        working light curve. Useful if the user is testing different processing
        methods and does not wish to reload a 'lightcurve' instance after each
        one.
        """
        if not hasattr(self, 'raw_lc'):
            print('Working light curve already raw light curve')

        else:
            self.time = self.raw_lc.time
            self.flux = self.raw_lc.flux
            self.flux_err = self.raw_lc.flux_err

            delattr(self, 'raw_lc')
            if hasattr(self, 'mask_arr'):
                delattr(self, 'mask_arr')

    def mask(self, **kwargs):
        """
        Method to mask out parts of the light curve. Wraps the lcprocess.mask
        function. Can be performed iteratively to mask out multiple areas or
        ranges of times in the light curve.

        Parameters
        ----------
        **kwargs : dict
            Arguments to be passed to `transit_tools.lcprocess.mask`, with the 
            exception of the out_mask argument.
        """
        if not hasattr(self, 'raw_lc'):
            raw_lc = LightCurve(self.time, self.flux, self.flux_err)
            self.raw_lc = raw_lc

        lc, mask_arr = mask(lc=self, out_mask=True, **kwargs)

        self._updatelc(lc)

        if hasattr(self, 'mask_arr'):
            self.mask_arr = [all(tup) for tup in zip(self.mask_arr, mask_arr)]
        else:
            self.mask_arr = mask_arr
            
    ###Method for user-provided stellar params (utils.py), input is a dict
       #these will be used as default and override any other gathered params
    
    ###Method to run TLS
    ###Method to run BLS
       ###combine with TLS into single search function
    def signal_search(self, routine='tls', plot_live=False, max_runs=5, sde=7.0,
                      exact=False, **kwargs):
        #!!Update to allow for search of known planets first and rejection of 
        #  signal until known signal is found. Quit after X trials!!
        #!!Update incomplete docustring!!
        #!!If periods too close, increase del_dur and run again or just raise
        #  warning and run again without logging a significant detection!!
        #!!Allow to run set number of iterations. Just set sde to 0?!!
        #!!Add option to optimize between searches to use exoplanet to subtract
        #  out transit model!!
        """
        Method to search the light curve for periodic signals using a variety of
        search algorithms, including Transit Least Squares and Box Least Squares
        from a number of packages. 
        
        Parameters
        ----------
        routine : str
            The desired routine to be used for finding periodic signals in the
            light curve. Options are 'tls' for Transit Least Squares and 'bls'
            for Box Least Squares.
        plot_live : bool
            Placeholder argument for the option to view vetting plots as runs 
            are finished for a more interactive transit search.
        max_runs : int
            The maximum number of runs allowed. This will be the number of runs
            for BLS runs.
        sde : float
            The threshold for the Source Detection Efficiency to be used to 
            determine whether a signal is significant or not in TLS runs.
        exact : bool
            Flag to indicate that the exact number of iterations specified in 
            max_runs will be performed.
        """
        self.routine = routine
        self.sde_thresh = sde

        if not hasattr(self, 'raw_lc'):
            self.raw_lc = LightCurve(self.time, self.flux,
                                     flux_err=self.flux_err)
        
        @property
        def routine(self):
            """The routine for searching for periodic signals"""
            return self.routine
        
        @routine.setter
        def routine(self, routine):
            if not routine in ['tls', 'TLS', 'bls', 'BLS']:
                raise ValueError('Please specify a supported routine type.')
            self.routine = routine
            
        self.results = np.array([])
        self.cleanlc = []
        self.bad_search = np.array([])
        run = 0
        
        #start iteration loop here and enclose both TLS and BLS code in it
        while run < max_runs:
            print('Run ' + str(run + 1) + ' for source ' + str(self.id))

            if len(self.cleanlc) == 0:
                time = self.time
                flux = self.flux
                flux_err = self.flux_err
            else:
                time = self.cleanlc[run-1].time
                flux = self.cleanlc[run-1].flux
                flux_err = self.cleanlc[run-1].flux_err
                
            if self.routine == 'tls' or self.routine == 'TLS':
                if not hasattr(self, 'star_params_tls'):
                    self.star_params_tls = None
                    
                #change so tic or star params don't have to be passed for sims
                    
                results_i, cleanlc_i, self.star_params_tls = tls_search(
                    time,
                    flux,
                    flux_err,
                    tic=self.tic,
                    star_params=self.star_params_tls,
                    clean_lc=True,
                    starparams_out=True,
                    **kwargs
                )
                
                if results_i.SDE < sde:
                    print('No further significant signals found. ' + str(run) +
                          ' total signals recovered')
                    self.bad_search = np.append(self.bad_search, results_i)
                    break
                
                self.results = np.append(self.results, results_i)
                self.cleanlc.append(LightCurve(cleanlc_i.time, cleanlc_i.flux,
                                               cleanlc_i.flux_err))
                
            if self.routine == 'bls' or self.routine == 'BLS':
                print('BLS is not implemented yet. Please be patient!')

                tmp_lc = LightCurve(time, flux, flux_err)
                
                results_i, cleanlc_i = bls_search(tmp_lc, clean_lc=True,
                                                  **kwargs)

                self.results = np.append(self.results, results_i)
                self.cleanlc.append(cleanlc_i)

            if plot_live:
                print('Please be patient. This feature is being worked on!')

            run += 1

    #Method to save signal_search results as new csv or to append to existing
            
    ###Method to plot vetting sheet based off of "search_method" property flag
    #which will be tls_def, tls_grz, tls_box, or bls
       #requires signal_search to run first, spits out error otherwise.
       ##Individual method to save all diagnostic plots and other methods to
       #   view each individually.

    def vetsheet(self, pls='all', save=False, **kwargs):
        #!!Maybe plot failed run if no significant signals are found?!!

        """
        Method to plot the vetting sheet for a given set of signal_search
        results.

        Parameters
        ----------
        pls : int or str
            signal_search iteration to display results for. If set to `'all'`, 
            all significant results will be displayed in separate windows. If 
            set to -1, the most recent set of results that did not meet the 
            significance threshold will be displayed.
        save : bool
            Flag to determine whether the plot is saved or not. Save filename 
            is specified and passed as a part of `**kwargs**`.
        **kwargs : dict
           Additional arguments to be passed to the 
           `transit_tools.plotting.tls_vetsheet` command.
        """
        if not hasattr(self, 'routine'):
            raise ValueError('Please run signal_search first')
        
        if pls == 'all':
            results = range(len(self.results))
        elif pls >= 0:
            results = [pls]

        if pls == -1:
            tls_vetsheet(self, results=-1, save=save, **kwargs)
        else:
            for i in results:
                tls_vetsheet(self, results=i, save=save, **kwargs)
                #add argument to combine pngs if necessary
                #if save:
                    #combine pngs

    #def saveplot(self, filename='summary.png'):
    #    """
    #    Function to save the vet sheet and other diagnostic plots.

    #    Parameters
    #    ----------
    #    filename : str
    #       Filename to save 
    #    """
    #    plt.savefig(filename)
        
    #Method to highlight where transits would be expected on current display
    #   axis based on known_pls. Different colors for different planets. Option
    #   to output expected transit times (into new attribute?) for each.

    #def save_fullplots(self): #flag to both display and save
        #output vetting sheet
        #vet = tls_vetsheet()
        #output river plot
        #others?
        #combine

    ###Method to run DAVE

    ###Method to run VESPA?

    ###Method to generate river plot? Implemented in plot method?

    ###method to implement REBOUND
    #   prioritizes exoplanet values over initial fit values but will take
    #   initial fit values if exoplanet wasn't run or explicity told to

    ###method to implement exoplanet

    ###method to implement BATMAN? (likely separate tool, or initialize object
    #   as simulated BATMAN light curve)
    #   Make whatever was simulated the self.known_pls parameter

    #gather stellar catalog info
    def get_starparams(self, cat='all'):
        #!!Check if star_params have already been collected or not!!
        #!!Add full docstrings!!
        #!!Add verbose option to call other function to print nicely!!
        
        """
        Method to gather stellar parameters and save them as the star_params
        and stellar_params_tls attributes.

        Parameters
        ----------
        cat : str
            Specify which calatog to retrieve stellar parameters from. If 
            ``'all'`` is entered, values will be gathered from the Gaia DR2 and 
            TESS Input Catalog, with the Gaia DR2 taking precedence
        **kwargs : dict
            Additional arguments to be passed to the ``'utils.catalog_info'`` 
            command.
        """
        #check to see if star_params already exist and if so, just skip to
        #   printing or something
        
        info, catalogs = catalog_info(tic=self.tic, cat=cat, out_cat=True)

        self.star_params = info
        self.catalogs = catalogs

        #need to check each value for NaNs or bad values?
        self.star_params_tls = {'mstar' : info['mass'], 'mlow' : info['e_mass'],
                                'mhigh' : info['e_mass'], 'rstar' : info['rad'],
                                'rlow' : info['e_rad'], 'rhigh' : info['e_rad']}

    #update star params command to allow user-input star params
        
    #print formatted star_params
        

    #print formatted search summary
    def searchsum(self):
        #!!Include indication if any result matches up with known planet!!
        """
        Method to display periodic signal search results in a user-friendly
        format. Wraps search_summary from utils.py.
        """
        if not hasattr(self, 'results'):
            raise ValueError('Please run a signal search first.')

        print('Search results for Source ' + str(self.tic))
        print('---------------------------------------')
        print(str(len(self.results)) + ' significant signals found')
        print('Used routine ' + str(self.routine))
        print('')
        
        for i in range(len(self.results)):
            if  ['tls', 'TLS'].count(self.routine):

                




                
                print('Signal ' + str(i+1))
                
                search_summary(self.results[i], routine='tls')

                #print('Disposition: ' + str(disp))
                
            if ['bls', 'BLS'].count(self.routine):
                search_summary(self.results[i], routine='bls')

    ##Method to update known_pls attribute with user-provided values
    def update_pls(self, name=None, per=None, t0=None, dur=None, params=None,
                   append=False):
        #!!Not finished. Revisit. Appending might not be working!!
        """
        Method to update the `lightcurve.known_pls` attribute with user-defined 
        values or a provided dictionary.

        Parameters
        ----------
        name : str or list or None
            Name(s) of the planets.
        per : float,` np.array` or None
            Period of the planet(s). Will expect a period in days.
        t0 : float, `np.array` or None
            Mid-transit time(s) of the first transit(s). Will expect a value in 
            MJD.
        dur : float, `np.array` or None
            Duration of transit(s). Will expect a value in days.
        params : dict, list or None
            Dictionary of additional parameters or all orbital parameters. Will
            expect keys similar to those output by 
            `transit_tools.utils.known_pls`.
        append : bool
            Flag whether or not to append provided planet parameters as a new 
            entry in the `lightcurve.known_pls` attribute or to overwrite the 
            current entry with the user-provided values.
        """
        if name is None and per is None and t0 is None and dur is None and params is not None:
            if not append:
                self.known_pls = params
            elif isinstance(params, list):
                for i in range(len(params)):
                    self.known_pls.append(params[i])
            else:
                self.known_pls.append(params)
                
        elif isinstance(any(name, per, t0, dur), list):
            conv = lambda i : i or ''
            [name, per, t0, dur] = [conv(i) for i in [name, per, t0, dur]]
            #turn all non-lists into lists
            if not isinstance(name, list): name = [name]
            if not isinstance(per, list): per = [per]
            if not isinstance(t0, list): t0 = [t0]
            if not isinstance(dur, list): dur = [dur]
            
            entries = max([len(i) for i in [name, per, t0, dur]])

            name += [''] * (entries - len(name))
            per += [''] * (entries - len(per))
            t0 += [''] * (entries - len(t0))
            dur += [''] * (entries - len(dur))
            
            #check if params is None and break up into separate dicts or just
            #   append each individually like above
            if params is None:
                pls = []
                
                for i in range(entries):
                    dic = {}
                    dic['canonical_name'] = name[i]
                    dic['orbital_period'] = per[i]
                    dic['transit_time'] = t0[i]
                    dic['transit_duration'] = dur[i]

                    pls.append(dic)

                self.known_pls = pls
                    
            
        #elif they're just floats or ints:

        else:
            pass
