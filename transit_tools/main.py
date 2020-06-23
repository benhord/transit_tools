import numpy
import pandas as pd
from lightkurve import LightCurve

class lightcurve(LightCurve):
    """Description
    of
    class

    Parameters
    ----------
    param : description
    """
    
    def __init__(self, tic, method="2min"):
        #Have 2min, ffi_ml (which points to Brian's lcs), and eleanor.
        #This assumes you're looking for TESS lcs, can be generalized later.
        #Defaults to looking for 2min lc using get_2minlc script from
        #vespa_wrapper but falls back to ffi_ml then to eleanor in this order
        #if 2min is specified with some explanatory messages.

    #Method to process light curve
        
    #Method to run TLS

    #Method to run BLS

    #Method to plot vetting sheet based off of "search_method" property flag
    #which will be tls_def, tls_grz, tls_box, or bls

    #Method to run DAVE

    #Method to run VESPA?

    #Method to generate river plot?
