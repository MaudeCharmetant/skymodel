import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from astropy.constants import k_B, h, c
import sz_tools as sz
import multiprocessing as mu
from joblib import Parallel, delayed

k_B = k_B.si.value
h = h.si.value
c = c.si.value

def mbb_spec_log(freq_log, A, T, beta, f_0 = 857e9, K_CMB=False, K_RJ=False):
    '''Computes the flux density of a modified black body in MJy/sr. 
    Parameters
    ----------
    freq: float or float array
        Frequency in Hz.
    A_dust: float
        Amplitude of the modified black body.
    T: float
        Dust temperature.
    beta: float
        Spectral index.
    f_0: float
        Pivot frequency in Hz.
    K_CMB: bool, optional
        The output will be given in units of K_CMB. Default: False
    K_RJ: bool, optional
        The output will be given in units of K_RJ. Default: False
    Returns
    -------
    dust: float or float array
        Flux density of a modified black body.
    '''

    #dust_log = np.log(A) + (3.+beta)*np.log(freq/f_0) + np.log(np.exp(h*f_0/k_B/T)-1) - np.log(np.exp(h*freq/k_B/T)-1)
    dust_log = np.log(A) + (3.+beta)*(freq_log-np.log(f_0)) + np.log(np.exp(h*f_0/k_B/T)-1) - np.log(np.exp(h*np.exp(freq_log)/k_B/T)-1)

    if K_CMB is True:
        dust_log = convert_units(freq, dust_log, mjy2cmb=True)
    elif K_RJ is True:
        dust_log = convert_units(freq, dust_log, mjy2rj=True)

    return(dust_log)

### read cita CIB maps

nside = 4096

npix = hp.pixelfunc.nside2npix(nside)
cib = np.zeros((5,npix), dtype = np.float32)

path = "/vol/arc3/data1/sz/CITA/newest/"

cib[0,:] = hp.fitsfunc.read_map(path + "cib_nu0143.fits", dtype = np.float64)
cib[1,:] = hp.fitsfunc.read_map(path + "cib_nu0217.fits", dtype = np.float64)
cib[2,:] = hp.fitsfunc.read_map(path + "cib_nu0353.fits", dtype = np.float64)
cib[3,:] = hp.fitsfunc.read_map(path + "cib_nu0545.fits", dtype = np.float64)
cib[4,:] = hp.fitsfunc.read_map(path + "cib_nu0857.fits", dtype = np.float64)


### define model and regression code

def fit_cib(data,i):
	freq = np.array([143,217,353,545,857])*1e9
	try:
		#p_opt, cov = curve_fit(lambda f, T, b: sz.mbb_spec(f, data[2], T, b, f_0 = 353e9), freq, data, p0=(10,1.5), sigma = np.ones(5)/100., bounds = [(5,0.5),(300,6)])
		p_opt, cov = curve_fit(lambda f, T, b: mbb_spec_log(f, data[2], T, b, f_0 = 353e9), np.log(freq), np.log(data), p0=(10,1.5), sigma = np.ones(5)/100., bounds = [(5,0.5),(300,6)])
		T_Dust = p_opt[0]
		beta = p_opt[1]
	except:
		T_Dust = -1
		beta = - 1

	#f.value += 1
	#print(p_opt)

	return([i, T_Dust, beta])


### fit the cita maps pixel by pixel using multiprocessing

threads = 128

dust_maps = Parallel(verbose=1, n_jobs=threads)(delayed(fit_cib)(cib[:,i], i) for i in np.arange(npix))
dust_maps = np.transpose(np.array(dust_maps))


### store obtained results as fits files

out_path = "/vol/arc3/data1/sz/CCATp_sky_model/workspace_maude/cita_components/"

A_Dust = np.float32(cib[2,:])
T_Dust = np.float32(dust_maps[1,:])
beta = np.float32(dust_maps[2,:])

hp.fitsfunc.write_map(out_path + "CITA_CIB_A_DUST_4096_log.fits", A_Dust, overwrite=True, dtype = np.float32)
hp.fitsfunc.write_map(out_path + "CITA_CIB_T_DUST_4096_log.fits", T_Dust, overwrite=True, dtype = np.float32)
hp.fitsfunc.write_map(out_path + "CITA_CIB_beta_DUST_4096_log.fits", beta, overwrite=True, dtype = np.float32)
