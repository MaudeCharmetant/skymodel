import numpy as np 
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

k_B = cst.k_B.si.value
h = cst.h.si.value
c = cst.c.si.value
T_CMB = 2.7255

def pix_reso(nside,arcmin=True): 
    
            
    """
    Compute the resolution of Healpy pixels in function of the nside of the map.  

    Parameters
    ----------
    
    nside : int
        Nside, number of separation of the healpy blocks. 
    arcmin : bool
        if True return the resolution in arcmin, otherwise in radians. 
        
    Returns
    -------
    float
        The resolution of one healpy pixel. 

    """
    
    #Number of pixels : 
    N_pix = hp.pixelfunc.nside2npix(nside)
    
    #Compute the unit : 
    if arcmin == True :        
        multi=(180*60/np.pi)        
    else:        
        multi=1
        
    #Compute the resolution :     
    reso = np.sqrt(((4*np.pi)/N_pix)*multi**2)
    
    #Feedback user : 
    print('The pixel resolution has been computed')

    return reso
  
  
def WN_map(nside,data_path,file_out,noise,unit_noise,arcmin,units,pictures_path): 

        
    """
    Function which create a White noise map for a given noise/arcmin or noise/radians.  

    Parameters
    ----------
    
    nside : int
        Nside, number of separation of the healpy blocks. 
    data_path : str
        Path where the datas will be stored.
    file_out : str
        Name of the White noise healpy map.
    noise : float 
        noise value desired in any units by radians or arcmin.
    unit_noise : float 
        resolution of the noise, for exemple 1' or 1 radians. 
    arcmin : bool 
        if true mean that the noise is given in /arcmin. 
    units : str 
        units of the noise 
    pictures_path : str 
        path where the pictures will be stored.
        
    Returns
    -------
    array
        Array contaning the Variarion of intensity produced by tSZ over the fequencies. 

    """
    
    #Compute the average noise in each pixel : 
    noise = (noise * unit_noise)/pix_reso(nside,arcmin)
        
    #Create White-noise map :
    Npix = hp.pixelfunc.nside2npix(nside) #Compute the number of pixels
    m =  np.random.normal(0,noise,Npix) #Random normal distribution centered over the desired noise
    
    #Display and save the distribution: 
    plt.hist(m,200,normed=True) #Display histogram of the distribution
    plt.title('$Noise$ $distribution$')
    plt.xlabel('$Noise$ $in$ '+ units)
    plt.ylabel('$Distribution$')
    plt.savefig(pictures_path + 'WN_PDF_'+str(int(noise))+units+'.png')

    #Display map and save map: 
    WN_map = hp.write_map(data_path + file_out,m, overwrite=True) #Save the distribution into a healpy map
    hp.mollview(m, title="White Noise map : n= "+ str(int(noise)) + units, norm='hist',unit=units)
    plt.savefig(pictures_path + 'WN_map_'+str(int(noise))+units+'.png') 
    plt.show()
    
    return WN_map 
  
  
def sample_sphere_uniform(n, mask = None, radec = True):
	'''Draws uniformly sampled tuples of coordinates on the sphere. 
	All-sky masks in the healpix format can be applied, in which case 
	masked areas will be excluded.
	Parameters
	----------
	n: int
		Number of data points to be drawn
	mask: float array, optional
		All-sky healpix mask. If a mask is used data points will 
		only be drawn in areas that are not masked. Default: None
	radec: bool
		Determines the coordinate system of the output. If True, 
		equatorial coordinates will be returned, i.e. RA, DEC (fk5). 
		If False, galactic coordinates are returned. Default: True
	Returns
	-------
	phi, theta
		longitude and latitude of sampled points in equatorial or 
		galactic coordinate system
	'''
    
	if mask is None:
		phi = 360 * np.random.random(n)
		theta = np.arccos(2*np.random.random(n) - 1)*180/np.pi - 90
	else:
		nside = hp.get_nside(mask)

		phi = np.zeros(n)
		theta = np.zeros(n)

		i = int(0)
		while i < n:
			phi_guess = 360 * np.random.random()
			theta_guess = np.arccos(2*np.random.random() - 1)*180/np.pi - 90

			index = hp.ang2pix(nside, phi_guess, theta_guess, lonlat = True)
			
			if mask[index] != 0:
				phi[i] = phi_guess
				theta[i] = theta_guess
				i += 1

	if radec is True:
		print("output will be fk5 coordinates (RA, DEC) for the equinox J2000")
		c = SkyCoord(phi, theta, frame='galactic', unit='deg')
		c_fk5 = c.transform_to('fk5')
		phi = c_fk5.ra.degree
		theta = c_fk5.dec.degree
	else:
		print("output will be galactic coordinates (glon, glat)")

	return(phi, theta)
  
    
    def convert_units(freq, values, cmb2mjy=False, mjy2cmb=False, rj2mjy=False, mjy2rj=False, cmb2rj=False, rj2cmb=False):
	'''Convert observed signal at given frequencies to different units. 
	Parameters
	----------
	freq: float or float array
		Frequency in Hz.
	values: float or float array
		Measured signal.
	cmb2mjy: bool, optional
		If True, the input is assumed to be K_CMB, the output will be MJy/sr.
		Default: False
	mjy2cmb: bool, optional
		If True, the input is assumed to be MJy/sr, the output will be K_CMB.
		Default: False
	rj2mjy: bool, optional
		If True, the input is assumed to be K_RJ, the output will be MJy/sr.
		Default: False
	mjy2rj: bool, optional
		If True, the input is assumed to be MJy/sr, the output will be K_RJ.
		Default: False
	cmb2rj: bool, optional
		If True, the input is assumed to be K_CMB, the output will be K_RJ.
		Default: False
	rj2cmb: bool, optional
		If True, the input is assumed to be K_RJ, the output will be K_CMB.
		Default: False
	Returns
	-------
	converted_signal: float or float array
		Converted signal
	'''

	x = h * freq / k_B / T_CMB
    
	if cmb2mjy is True:
		conversion = 1e20 * 2*k_B**3.*T_CMB**2. / (h * c)**2. * x**4. * np.exp(x) / (np.exp(x)-1)**2.
	elif mjy2cmb is True:
		conversion = 1/(1e20 * 2*k_B**3.*T_CMB**2. / (h * c)**2. * x**4. * np.exp(x) / (np.exp(x)-1)**2.)
	elif rj2mjy is True:
		conversion = 1e20 * 2*freq**2.*k_B/c**2.
	elif mjy2rj is True:
		conversion = 1/(1e20 * 2*freq**2.*k_B/c**2.)
	elif cmb2rj is True:
		conversion = (k_B*T_CMB/h)**2. * x**4. * np.exp(x) / (np.exp(x)-1)**2. / freq**2.
	elif rj2cmb is True:
		conversion = 1/((k_B*T_CMB/h)**2. * x**4. * np.exp(x) / (np.exp(x)-1)**2. / freq**2.)        
	else:
		print("Not sure which units are given and what should be returned.")

	converted_signal = conversion * values

	return(converted_signal)


def generate_cib(freq, nside_out=4096, beam_FWHM = None, template = "SO", unit = "mjy"):
    '''Computes an all-sky CIB map at a given frequency and nside based on . 

    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz.
    nside_out: float
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    beam_FWHM: bool, optional
        If set, the output will be convolved with a gaussian. The FWHM of the Gaussian
        in units of arcmin is given by the provided value. Default: None
    template: bool, optional
        Determines the all-sky foregrounds templates to be used to build the sky model.
        If 'SO' is chosen, the Simons Observatory sky model provided by Colin Hill and 
        based on the simulations by Sehgal et al. (2010) is used. If 'CITA' is chosen,
        the used templates will be based on the WebSky Extragalactic CMB Mocks provided 
        by CITA. Default: 'SO'
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'mjy'.

    Returns
    -------
    cib: float array
        Healpix all-sky map of the CIB mission that the specified frequency.
    '''

    #Load all-sky parameter value maps	
    if template == "SO":
        path = "/vol/arc3/data1/sz/CCATp_sky_model/workspace_jens/so_components/"
        A = hp.fitsfunc.read_map(path + "SO_CIB_A_DUST_4096.fits", dtype = np.float32)    
        T = hp.fitsfunc.read_map(path + "SO_CIB_T_DUST_4096.fits", dtype = np.float32)
        beta = hp.fitsfunc.read_map(path + "SO_CIB_beta_DUST_4096.fits", dtype = np.float32)
        f_0 = 353e9
    elif template == "CITA":
        path = "/vol/arc3/data1/sz/CCATp_sky_model/workspace_jens/cita_components/"
        A = hp.fitsfunc.read_map(path + "CITA_CIB_A_DUST_4096.fits", dtype = np.float32)    
        T = hp.fitsfunc.read_map(path + "CITA_CIB_T_DUST_4096.fits", dtype = np.float32)
        beta = hp.fitsfunc.read_map(path + "CITA_CIB_beta_DUST_4096.fits", dtype = np.float32)
        f_0 = 353e9
    else:
        print("Waring: Unknown template requested! Output will be based on SO sky model")

    #Compute CIB brightness at given frequency
    cib = A * (freq/f_0)**(3.+beta) * (np.exp(h*f_0/k_B/T)-1) / (np.exp(h*freq/k_B/T)-1)
    del A, T, beta
    
    #Re-bin map if necessary
    if hp.get_nside(cib) != nside_out:
        cib = hp.pixelfunc.ud_grade(cib, nside_out = nside_out)

    #Smooth map if necessary
    if beam_FWHM is not None:
        print("begin smoothing")
        cib = hp.sphtfunc.smoothing(cib, iter = 0, lmax = 2*nside-1, fwhm = beam_FWHM/60*np.pi/180)

    #Convert units if necessary
    if unit == "mjy":
        None
    elif unit == "cmb":
        cib = convert_units(freq, cib, mjy2cmb=True)
    elif unit == "rj":
        cib = convert_units(freq, cib, mjy2rj=True)
    else:
        print("Waring: Unknown unit! Output will be in MJy/sr")

    #Return output
    return(np.float32(cib))
