import numpy
import pandas as pd
from lightkurve.lightcurve import *

from .fetch_lc import gather_lc
from .search import *

class lightcurve(LightCurve):
    """Description
    of
    class

    #pass mission-specific classes to preserve some speciality?

    Parameters
    ----------
    param : description
    """
    
    def __init__(self, tic, method="2min", sector=None, mission='TESS'):
        #Have 2min, ffi_ml (which points to Brian's lcs), and eleanor.
        #This assumes you're looking for TESS lcs, can be generalized later.
        #Defaults to looking for 2min lc using get_2minlc script from
        #vespa_wrapper but falls back to ffi_ml then to eleanor in this order
        #if 2min is specified with some explanatory messages. Use property
        #decorators to control which sectors can be used and the cadence.
        self.tic = int(tic)
        self.sector = sector
        self.method = method
        
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

        
    ###Method to process light curve

    ###Method for user-provided stellar params (helper.py), input is a dict
       #these will be used as default and override any other gathered params
    
    ###Method to run TLS
    ###Method to run BLS
       ###combine with TLS into single search function
    def signal_search(self, routine='tls', plot_live=False):
        """
        Method to search the light curve for periodic signals using a variety of
        search algorithms, including Transit Least Squares and Box Least Squares
        from a number of packages. 

        !!Update to allow for search of known planets first and rejection of 
        signal until known signal is found. Quit after X trials!!

        Parameters
        ----------
        routine : str
        """
        print('Searching for periodic signal...')
        #search until signal is no longer significant. Dumps results into arrays
        #to account for finding multiple planets

    ###Method to plot vetting sheet based off of "search_method" property flag
    #which will be tls_def, tls_grz, tls_box, or bls
       #requires signal_search to run first, spits out error otherwise.
       #make it customizeable so it doesn't print a ton of stuff each time?
       #ie allow the user to say 'not XYZ' or 'include XYZ' beyond defaults.

    ###Method to run DAVE

    ###Method to run VESPA?

    ###Method to generate river plot? Implemented in plot method?

    #method to implement REBOUND

    #method to implement exoplanet

    #print formatted catalog info
