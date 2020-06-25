import os
import pickle
import matplotlib.pyplot as plt

sector = 13
tic = 254113311

def fetch_ffi_ml(tic, sector):

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
