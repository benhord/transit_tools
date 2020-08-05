import numpy
import pandas as pd
from lightkurve.lightcurve import *
import matplotlib.pyplot as plt

from .fetch_lc import gather_lc
from .search import *
from .utils import *
from .plotting import *
from .lcprocess import *

class lightcurve(LightCurve):
    """Description
    of
    class

    #pass mission-specific classes to preserve some speciality?

    Parameters
    ----------
    param : description
    """
    
    def __init__(self, obj, method="2min", sector=None, mission='TESS',
                 find_knownpls=True):
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
        
        if isinstance(obj, str):
            try:
                self.tic = name_to_tic(obj)
                self.name = str(obj)
            except:
                raise ValueError('For some reason, the name provided was not ' +
                                 'found on the MAST. Please try again.')

        elif isinstance(obj, int):
            try:
                self.name = tic_to_name(obj)
                self.tic = int(obj)
            except:
                print('For some reason, the name provided was not found on ' +
                      'the MAST. Proceeding with just the TIC.')
        
        self.sector = sector
        self.method = method
        self.star_params = None
        
        @property
        def sector(self):
            """The TESS sector(s) being used"""
            return self.sector

        @sector.setter
        def sector(self, sector):
            if not isinstance(sector, list) and (value != None):
                values = np.array(sector)
            self.sector = sector

        @property
        def method(self):
            """The method for acquiring the TESS light curve"""
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

        self.lc, self.method, self.sector = gather_lc(self.tic,
                                                      method=str(self.method),
                                                      sectors=self.sector,
                                                      return_method=True,
                                                      return_sectors=True)

        super().__init__(time=self.lc.time, flux=self.lc.flux,
                         flux_err=self.lc.flux_err)

        if find_knownpls:
            try:
                self.known_pls = known_pls(self.name)
            except:
                self.known_pls = None
            
        
    ###Method to process light curve
    def process(self, **kwargs):
        """
        Method to process light curve for further analysis. Must be performed on
        a LightCurve object.

        Parameters
        ----------
        kwargs
           Arguments to be passed to the process_lc.
        """
        lc = process_lc(self, **kwargs)

        self.updatelc(lc)

    #method to update lc
    def updatelc(self, lc):
        """
        Function to update the light curve that will be used if manual
        alteration is performed that the user wishes to save without completely 
        overriding the original light curve.

        Parameters
        ----------
        lc : 'LightCurve' object
        """
        #check for self.flux_err attribute
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
        Function to bring the raw, original light curve back as the main, 
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

    ###Method for user-provided stellar params (utils.py), input is a dict
       #these will be used as default and override any other gathered params
    
    ###Method to run TLS
    ###Method to run BLS
       ###combine with TLS into single search function
    def signal_search(self, routine='tls', plot_live=False, max_runs=5, sde=7.0,
                      **kwargs):
        """
        Method to search the light curve for periodic signals using a variety of
        search algorithms, including Transit Least Squares and Box Least Squares
        from a number of packages. 

        !!Update to allow for search of known planets first and rejection of 
        signal until known signal is found. Quit after X trials!!
        !!Update incomplete docustring!!

        Parameters
        ----------
        routine : str
           The desired routine to be used for finding periodic signals in the
           light curve. Options are 'tls' for Transit Least Squares and 'bls'
           for Box Least Squares.
        """
        self.routine = routine
        self.sde = sde

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
            print('Run ' + str(run + 1) + ' for source ' + str(self.tic))
            
            if self.routine == 'tls' or self.routine == 'TLS':
                if not hasattr(self, 'star_params_tls'):
                    self.star_params_tls = None

                if len(self.cleanlc) == 0:
                    time = self.time
                    flux = self.flux
                    flux_err = self.flux_err
                else:
                    time = self.cleanlc[run-1].time
                    flux = self.cleanlc[run-1].flux
                    flux_err = self.cleanlc[run-1].flux_err
                    
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
                
            if routine == 'bls' or routine == 'BLS':
                print('BLS is not implemented yet. Please be patient!')

            if plot_live:
                print('Please be patient. This feature is being worked on!')

            run += 1

    #Method to save signal_search results as new csv or to append to existing
            
    ###Method to plot vetting sheet based off of "search_method" property flag
    #which will be tls_def, tls_grz, tls_box, or bls
       #requires signal_search to run first, spits out error otherwise.
       ##Individual method to save all diagnostic plots and other methods to
       #   view each individually.

    def vetsheet(self, pls='all', **kwargs):
        """
        Function to plot the vetting sheet for a given set of signal_search
        results.

        Parameters
        ----------
        pls : int or str
           signal_search iteration to display results for. If set to 'all', all
           significant results will be displayed in separate windows. If set to
           -1, the most recent set of results that did not meet the significance
           threshold will be displayed.
        kwargs
           Additional arguments to be passed to the tls_vetsheet command.
        """
        if not hasattr(self, 'routine'):
            raise ValueError('Please run signal_search first')
        
        if pls == 'all':
            results = range(len(self.results))
        elif pls >= 0:
            results = [pls]

        if pls == -1:
            tls_vetsheet(self, results=-1, **kwargs)
        else:
            for i in results:
                tls_vetsheet(self, results=i, **kwargs)

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

    #method to implement exoplanet

    #method to implement BATMAN? (likely separate tool, or initialize object as
    #   simulated BATMAN light curve)
    #   Make whatever was simulated the self.known_pls parameter

    #print formatted catalog info

    #print formatted search summary
    def searchsum(self):
        """
        Function to display periodic signal search results in a user-friendly
        format. Wraps search_summary from utils.py.

        !!Include indication if any result matches up with known planet!!
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

            if ['bls', 'BLS'].count(self.routine):
                search_summary(self.results[i], routine='bls')

    ##Method to update known_pls attribute with user-provided values
