import numpy
import pandas as pd
from lightkurve import LightCurve

from transit_tools.fetch_lc import ffi_ml, ffi_eleanor, 2min_pdcsap

class lightcurve(LightCurve):
    """Description
    of
    class

    Parameters
    ----------
    param : description
    """
    
    def __init__(self, tic, method="2min", sector=None):
        #Have 2min, ffi_ml (which points to Brian's lcs), and eleanor.
        #This assumes you're looking for TESS lcs, can be generalized later.
        #Defaults to looking for 2min lc using get_2minlc script from
        #vespa_wrapper but falls back to ffi_ml then to eleanor in this order
        #if 2min is specified with some explanatory messages. Use property
        #decorators to control which sectors can be used and the cadence.

        self.tic = int(tic)

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
            return self.method = method
            
    #Method to process light curve
        
    #Method to run TLS

    #Method to run BLS

    #Method to plot vetting sheet based off of "search_method" property flag
    #which will be tls_def, tls_grz, tls_box, or bls

    #Method to run DAVE

    #Method to run VESPA?

    #Method to generate river plot?

    #method to implement REBOUND

    #method to implement exoplanet
