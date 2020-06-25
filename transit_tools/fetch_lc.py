import os
import pickle
import matplotlib.pyplot as plt
from astroquery.mast import Observations

sector = 13
tic = 254113311

def 2min_pdcsap(tic, sector, thresh=None):

    #Query MAST for all instances of observations associated with TIC
    obsTable = Observations.query_criteria(dataproduct_type=['timeseries'], 
                                           target_name=tic,
                                           obs_collection='TESS')

    lc_str = "lc.fits"
    good_ind = []

    if thresh: #get sectors from gt or lt range
        for i in range(len(obsTable['dataURL'])):
            if (thresh[:2] == 'lt' and
                str(obsTable['dataURL'][i]).find(lc_str) > 0 and
                int(obsTable['dataURL'][i][37:41]) < int(thresh[2:])):
                good_ind.append(i)
            if (thresh[:2] == 'gt' and
                str(obsTable['dataURL'][i]).find(lc_str) > 0 and
                int(obsTable['dataURL'][i][37:41]) > int(thresh[2:])):
                good_ind.append(i)
    
    elif sectors != 'all': #get just specified sectors
        for sec in sectors:
            sec_str = "s" + str(sec).zfill(4)
            for i in range(len(obsTable['dataURL'])):
                if (str(obsTable['dataURL'][i]).find(lc_str) > 0 and
                    str(obsTable['dataURL'][i]).find(sec_str) > 0):
                    good_ind.append(i)

    elif sectors == 'all': #get all sectors from list
        for i in range(len(obsTable['dataURL'])):
            if str(obsTable['dataURL'][i]).find(lc_str) > 0:
                good_ind.append(i)

    if len(good_ind) == 0:
        #print("WARNING: NO VALID SECTORS IN GIVEN RANGE!!")
        raise ValueError("ERROR: No valid sectors in given range!! " +
                         "Try running the command again with a different " +
                         "sector range.")
                
    obsTable=obsTable[good_ind] #returns only good indices for table

    #get lightcurve from MAST
    id = str(tic).zfill(16)
    for i in range(len(good_ind)):
        sec_str = "s" + str(obsTable['dataURL'][i][37:41])
        lc_loc=("https://archive.stsci.edu/missions/tess/tid/" + str(sec_str) +
                "/" + id[:4] + "/" + id[4:8] + "/" + id[8:12] + "/" +
                id[12:16] + "/" + str(obsTable['dataURL'][i])[18:])

        if i == 0:
            lc = (TessLightCurveFile(lc_loc).PDCSAP_FLUX.normalize().
                  remove_nans())
            
        else:
            lc_new = (TessLightCurveFile(lc_loc).PDCSAP_FLUX.normalize().
                  remove_nans())
            lc = lc.append(lc_new)

    return lc
    
def ffi_ml(tic, sector): #edit to iterate through sectors

    # Gather all light curve file paths and TICs for specified sector
    light_curve_files = []
    path='/data/tessraid/bppowel1/tesslcs_sector_'+str(sector)+'_104'

    for (dirpath,dirnames,filenames) in os.walk(path):
        for name in filenames:
            light_curve_files.append(os.path.join(dirpath,name))

    light_curve_files=[t for t in light_curve_files if 'tesslc_' in t]
    tics=[int(t.split('.')[0].split('_')[-1]) for t in light_curve_files]

    # Find path to individual light curve
    path = [s for s in light_curve_files if str(tic) in s]
    print("Light curve found at %s" % path[0])

    # Import light curve from pickle file
    fp = open(str(path[0]), 'rb')
    data = pickle.load(fp)
    fp.close()

    ra = data[1]
    dec = data[2]
    Tmag = data[3]
    camera = data[4]
    chip = data[5]
    time = data[6]
    raw_flux = data[7]
    corr_flux = data[8] #use this one
    pca_flux = data[9]
    flux_err = data[10]

    q = data[11] == 0
    time = time[q]
    raw_flux = raw_flux[q]
    corr_flux = corr_flux[q]
    pca_flux = pca_flux[q]
    
    return time, raw_flux, corr_flux, pca_flux, flux_err, ra, dec, Tmag, camera, chip
    
#plt.scatter(time, corr_flux)
#plt.show()
