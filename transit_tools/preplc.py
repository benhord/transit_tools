import numpy as np

def gather_lc(method, sectors):
    if method == '2min':
        #try and except to change to next method?
        try:
            lc = get_2minlc(sectors)
        except:
            print('No TESS 2 minute light curves found! Trying FFIs...')
            method = 'ffi_ml'
            
    elif method == 'ffi_ml':
        try:
            lc = get_mlffi(sectors)
        except:
            print('No ML light curves found locally. Trying with eleanor...')
            method = 'eleanor'
            
    elif method == 'eleanor':
        try:
            lc = get_eleanor(sectors)
        except:
            raise ValueError('No light curves found for the specified sectors!')
            
#get_2minlc

#search for ML FFIs

#eleanor
