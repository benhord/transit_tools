import numpy as np
from astroquery.mast import Catalogs, Observations

#function to import catalog info for source
   #stellar information
   #RA, DEC
   #whatever is in the TIC
   #GAIA position/other info?

   #break into separate functions to provide greater flexibility?
   #for any stellar information, put any relevant info into both a self.star
   #  attribute but also a self.star_params_tls attribute

#function to import observation information
   #observation information based on self.method keyword (eg ccd, sector, etc.)

#function to find and report data on planets already-known in the system

#function to update stellar params of lightcurve object

#common name processing so that a TIC isn't the only thing that can be passed

#calculate rms of a data set
def rms(data, norm_val=1.):
    """
    Calculates the Root Mean Square of the provided data.

    Parameters
    ----------
    data : numpy array
       Data for which the root mean square will be calculated
    norm_val : float or numpy array
       The value(s) that the data array is normalized to or the value of the
       model that the data values are being compared against if RMS is non-
       uniform.

    Returns
    -------
    rms : float
       The calculated root mean square of the data.
    """
    rms = np.sqrt(np.sum((data - norm_val) ** 2) / len(data))

    return rms
