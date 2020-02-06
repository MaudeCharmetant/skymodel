import numpy as np 
import healpy as hp
from astropy.io import fits
from astropy import constants as cst
import pysm
from tqdm import tqdm 
from pysm.nominal import models

k_B = cst.k_B.si.value
h = cst.h.si.value
c = cst.c.si.value
T_CMB = 2.7255


def convert_units(freq, values, cmb2mjy = False, mjy2cmb = False, rj2mjy = False, mjy2rj = False, cmb2rj = False, rj2cmb = False):

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


def px_size(nside,arcmin=True):   
            
    """
    Compute the size of Healpy pixels in function of the nside of the map.  

    Parameters
    ----------
    
    nside : int
        Nside, number of separation of the healpy blocks. 
    arcmin : bool
        if True return the size in arcmin, otherwise in radians. 
        
    Returns
    -------
    float
        The size of one healpy pixel. 
    """
    
    #Number of pixels : 
    N_pix = hp.pixelfunc.nside2npix(nside)
    
    #Compute the unit : 
    if arcmin == True :        
        multi=(180*60/np.pi)        
    else:        
        multi=1
        
    #Compute the resolution :     
    size = np.sqrt(((4*np.pi)/N_pix)*multi**2)

    return(size)


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


def project_maps(allsky_map, RA, DEC, map_size = 10, pixel_size = 0.4):

	'''Creates gnomic projections of HEALPix all-sky maps.

	Parameters
	----------

	allsky_map: float array
		numpy array containing a healpy all-sky map with a valid nside
	RA: float or float array array, optional
		Right acention of objects, fk5 coordinates are required.
	DEC: float or float array, optional
		Declination of objects, fk5 coordinates are required.
	map_size: float, optional
		Size of the desired projected map in degree, map will be square. Default: 10
	pixel_size: float, optional
		Pixel size of the desired projected map in arcmin. Default: 0.4

	Returns
	-------

	maps: float array
		Array containing the projected maps
	'''

	RA = np.asarray(RA)
	DEC = np.asarray(DEC)

	n = len(RA)
	npix = int(map_size*60 / pixel_size)

	maps = np.zeros((n, npix, npix), dtype=np.float32)

	for i in np.arange(n):
			progress(i,n, percentage=True)
			maps[i,:,:] = hp.visufunc.gnomview(allsky_map, coord = ['G', 'C'], rot=[RA[i],DEC[i]], reso = pixel_size, xsize = npix, return_projected_map = True, no_plot=True) 

	return(maps)
    

def simulate_gal_foregrounds(freq, components = "all", nside_out = 4096, lmax = None, beam_FWHM = None, unit = "cmb"):

    '''Computes an all-sky galactic foregrounds noise map at a given frequency and nside using 
    the Python Sky model (PySM, Thorne et al. 2017), which is build from Planck data. 

    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz. 
    components: list of stings
        List of components to be included in the galactic foreground model. Possible 
        components are "gal_synchrotron", "gal_dust", "gal_freefree", and "gal_ame".
        Alternatively, all components are included if the variable is set to "all".
        Default: "all"	
    nside_out: float
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float
        Maximum value of the multipolemoment at which the atmospheric power spectrum
        wil be computed. Default: 3*nside_out-1    
    beam_FWHM: bool, optional
        If set, the output will be convolved with a gaussian. The FWHM of the Gaussian
        in units of arcmin is given by the provided value. Default: None
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.

    Returns
    -------
    foregrounds: float array
        Healpix all-sky map of the galactic foregrounds at the specified frequency.
    '''

    if components == "all":
        components = ["gal_synchrotron", "gal_dust", "gal_freefree", "gal_ame"]	

    #Define foreground model
    sky_config = {
        'synchrotron' : models("s1", 512),
        'dust' : models("d1", 512),
        'freefree' : models("f1", 512),
        'ame' : models("a1", 512),
    }

    #Initialise Sky 
    sky = pysm.Sky(sky_config)

    #Compute component maps
    foregrounds = np.zeros(hp.nside2npix(512))
    rj2mjy_factor = convert_units(freq, 1e-6, rj2mjy = True)
    
    if "gal_synchrotron" in components:
        foregrounds += sky.synchrotron(freq/1e9)[0,:] * rj2mjy_factor

    if "gal_dust" in components:
        foregrounds += sky.dust(freq/1e9)[0,:] * rj2mjy_factor
          
    if "gal_freefree" in components:
        foregrounds += sky.freefree(freq/1e9)[0,:] * rj2mjy_factor
        
    if "gal_ame" in components:
        foregrounds += sky.ame(freq/1e9)[0,:] * rj2mjy_factor

    #Define smoothing kernal, upgrade, and smooth map
    fwhm = 10
    if beam_FWHM is not None:
        if beam_FWHM > fwhm:
            fwhm = np.sqrt(fwhm*2 + beam_FWHM**2)
            
    foregrounds = hp.pixelfunc.ud_grade(foregrounds, nside_out = nside_out)
    foregrounds = hp.sphtfunc.smoothing(foregrounds, iter = 0, fwhm = fwhm/60*np.pi/180, lmax = lmax)
    
    #Convert units if necessary
    if unit == "mjy":
        None
    elif unit == "cmb":
        foregrounds = convert_units(freq, foregrounds, mjy2cmb=True)
    elif unit == "rj":
        foregrounds = convert_units(freq, foregrounds, mjy2rj=True)
    else:
        print("Waring: Unknown unit! Output will be in K_CMB")
        
    return(foregrounds)


def simulate_cib(freq, nside_out = 4096, beam_FWHM = None, template = "SO", unit = "cmb"):

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
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.

    Returns
    -------
    cib: float array
        Healpix all-sky map of the CIB mission that the specified frequency.
    '''

    #Load all-sky parameter value maps
    if template != "SO" and template != "CITA":
        print("Waring: Unknown template requested! Output will be based on SO sky model")
        template = "SO"

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
        print("Waring: Unknown unit! Output will be in K_CMB")

    #Return output
    return(np.float32(cib))


def simulate_radio_ps(freq, nside_out = 4096, lmax = None, beam_FWHM = None, template = "SO", unit = "cmb"):

    '''Computes an all-sky radio point source map at a given frequency and nside based on 
    the simulations provided by Sehgal et al. (2010), which have been recalibrated by the
    SO collaboration. Du to the complex SEDs of the radio sources, bilinear interpolation
    is applied for input frequencies between 27 and 353 GHz. For higher frequencies, a null
    map is returned.

    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz.
    nside_out: float
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float
        Maximum value of the multipolemoment at which the atmospheric power spectrum
        wil be computed. Default: 3*nside_out-1            
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
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.

    Returns
    -------
    radio_ps: float array
        Healpix all-sky map of the radio point source emission at the specified frequency.
    '''

    if lmax is None:
        lmax = int(3*nside_out-1)
    
    #Define path to data and frequencies
    path = "/vol/arc3/data1/sz/SO_sky_model/sky_maps/"    
    nu = np.array([27,30,39,44,70,93,100,143,145,217,225,280,353])*1e9
    nu_names = ['027','030','039','044','070','093','100','143','145','217','225','280','353']

    #Read data files
    npix = hp.pixelfunc.nside2npix(4096)
    data = np.zeros((len(nu), npix), dtype = np.float32)    

    for i in np.arange(len(nu)):
        file_name = path + nu_names[i] + "_rad_pts_healpix_nopell_Nside4096_DeltaT_uK_fluxcut148_7mJy_lininterp.fits"
        data[i,:] = hp.fitsfunc.read_map(file_name, dtype = np.float32) * convert_units(nu[i], 1e-6, cmb2mjy=True)

    #Interpolate data points
    radio_ps = np.zeros(hp.nside2npix(4096))

    if template != "SO" and template != "CITA":
        print("Waring: Unknown template requested! Output will be based on SO sky model")
        template = "SO"	
	
    if template == "SO":
        
        if freq > 353e9:
            print("Warning: Input frequency lies beyoned the 353 GHz. Since higher frequencies are not constraint by simulations, the data will be 0.")
        else:
            for i in tqdm(np.arange(hp.nside2npix(4096))):
                radio_ps[i] = np.interp(freq, nu, data[:,i])
                
    elif template == "CITA":
        print("Warning: No radio PS template provided by the CITA simulations")

    #Re-bin map if necessary
    if hp.get_nside(radio_ps) != nside_out:
        radio_ps = hp.pixelfunc.ud_grade(radio_ps, nside_out = nside_out)

    #Smooth map if necessary
    if beam_FWHM is not None:
        print("begin smoothing")
        radio_ps = hp.sphtfunc.smoothing(radio_ps, iter = 0, lmax = lmax, fwhm = beam_FWHM/60*np.pi/180)

    #Convert units if necessary
    if unit == "mjy":
        None
    elif unit == "cmb":
        radio_ps = convert_units(freq, radio_ps, mjy2cmb=True)
    elif unit == "rj":
        radio_ps = convert_units(freq, radio_ps, mjy2rj=True)
    else:
        print("Waring: Unknown unit! Output will be in K_CMB")

    #Return output
    return(np.float32(radio_ps))


def simulate_cmb(freq, cl_file = None, lensed = True, nside_out = 4096, lmax = None, beam_FWHM = None, template = "SO", unit = "cmb"): 
         
    """
    Function that computes a CMB map from a power spectrum.

    Parameters
    ----------

    freq: float or float array
        Frequency of the output map in Hz.
    cl_file : str or array 
        Name of the .dat contaning the values of the power spectrum given by CAMB.
        Or array containing the power spectrum to generate random maps in Kelvin.
    lensed : bool 
    	if True select the lensed CMB map when possible. This is only possible for 'CITA' and 'Sehgal.'
    nside_out: float
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float
        Maximum value of the multipolemoment at which the atmospheric power spectrum
        wil be computed. Default: 3*nside_out-1            
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
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.
    
    data_path : str
        Path were the data of the maps are stored and we the cutout are going to be stored. 
    data_save : str 
    	Path were the datas will be saved.
        
    Returns
    -------
    array
        Contaning the CMB map. 

    """
                
    if cl_file is not None:
    
        #Load the datas : 
        data = np.loadtxt(cl_file)
        ell = data[:,0]
        TT = data[:,1] * 2 * np.pi / ell / (ell+1)  # Take only the first column, which is the temperature T

        #Because the monopole and dipole are not in the data of the power spectrum given by CAMB we need to add them back
        #They need to be 0 because we usually remove them no to influence our studies. 
        TT_final = np.insert(TT,[0],[0,0])

        #Compute the CMB map from the power spectrum :
        CMB = hp.sphtfunc.synfast(TT_final, nside_out, lmax = lmax)
 
    else: 

        if template != "SO" and template != "CITA" and template != "Sehgal":
            print("Waring: Unknown template requested! Output will be based on SO sky model")
            template = "SO"			
		
        if template == 'CITA':
        
            if lensed == True: 
        
                CMB = hp.read_map('/vol/arc3/data1/sz/CCATp_sky_model/workspace_maude/'+' CMB_lensed_CITA_mK', dtype = np.float32)

            else: 
            
                CMB = hp.read_map('/vol/arc3/data1/sz/CCATp_sky_model/workspace_maude/'+' CMB_unlensed_CITA_mK', dtype = np.float32)     
    
        if template == 'SO': 

            data_path = '/vol/arc3/data1/sz/SO_sky_model/CMB_SZ_maps/'
            file_name = 'Sehgalsimparams_healpix_4096_KappaeffLSStoCMBfullsky_phi_SimLens_Tsynfastnopell_fast_lmax8000_nside4096_interp2.5_method1_1_lensed_map.fits'

        if template == 'Sehgal':

            data_path = '/vol/arc3/data1/sz/Sehgal/'

            if lensed == True:

                file_name = '030_lensedcmb_healpix.fits'

            else: 

                file_name = '030_unlensedcmb_healpix.fits'

        CMB = hp.read_map(data_path + file_name, dtype = np.float32)

    #Re-bin map if necessary
    if hp.get_nside(CMB) != nside_out:
        radio_ps = hp.pixelfunc.ud_grade(CMB, nside_out = nside_out)

    #Smooth map if necessary
    if beam_FWHM is not None:
        print("begin smoothing")
        CMB = hp.sphtfunc.smoothing(CMB, iter = 0, lmax = lmax, fwhm = beam_FWHM/60*np.pi/180)

    #Convert units if necessary
    if unit == "mjy":
        CMB = convert_units(freq, CMB/1e6, cmb2mjy=True)
    elif unit == "cmb":
        CMB /= 1e6
    elif unit == "rj":
        CMB = convert_units(freq, CMB/1e6, cmb2rj=True)
    else:
        print("Waring: Unknown unit! Output will be in K_CMB")   

    return(CMB)


def simulate_tSZ(freq, nside_out = 4096, lmax = None, beam_FWHM = None, template = "SO", unit = "cmb"):
    
    """
    Function which compute tSZ maps at different frequencies and different nside. 
    
    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz.
    nside_out: float
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float
        Maximum value of the multipolemoment at which the atmospheric power spectrum
        wil be computed. Default: 3*nside_out-1            
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
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.
        
    Returns
    -------
    tSZ: float array
        Healpix allsky map contaning the tSZ map at a given frequency. 
    """

    #Read data
    if template != "SO" and template != "CITA" and template != "Sehgal":
        print("Waring: Unknown template requested! Output will be based on SO sky model")
        template = "SO"

    if template == 'SO':
        
        data_path='/vol/arc3/data1/sz/SO_sky_model/CMB_SZ_maps/'
        file_in ='tSZ_skymap_healpix_nopell_Nside4096_y_tSZrescale0p75.fits'
        y_map = hp.read_map(data_path + file_in, dtype = np.float32)

        
    if template == 'CITA': 
        
        data_path='/vol/arc3/data1/sz/CITA/'
        file_in = 'tsz.fits'
        y_map = hp.read_map(data_path + file_in, dtype = np.float32)        

        
    if template == 'Sehgal': 
        
        data_path='/vol/arc3/data1/sz/Sehgal/'
        file_in='030_tsz_healpix.fits'
        
        #Create Compton y-map : 
        tSZ_30GHz = hp.read_map(data_path + file_in, dtype = np.float32) * convert_units(freq, 1e-6, mjy2cmb=True)
           
        x_nu = np.array((h*30e9)/(k_B*T_CMB))    
        tSZ_SED = ((x_nu*(np.exp(x_nu)+1)/(np.exp(x_nu)-1))-4)

        y_map = tSZ_30GHz / tSZ_SED / T_CMB

    #Get tSZ at different frequencies : 
    x_nu = np.array((h*freq)/(k_B*T_CMB))    
    t_SZ_SED = ((x_nu*(np.exp(x_nu)+1)/(np.exp(x_nu)-1))-4)

    tSZ = (y_map * t_SZ_SED * T_CMB)
        
    #Re-bin map if necessary
    if hp.get_nside(tSZ) != nside_out:
        tSZ = hp.pixelfunc.ud_grade(tSZ, nside_out = nside_out)

    #Smooth map if necessary
    if beam_FWHM is not None:
        print("begin smoothing")
        tSZ = hp.sphtfunc.smoothing(tSZ, iter = 0, lmax = lmax, fwhm = beam_FWHM/60*np.pi/180)

    #Convert units if necessary
    if unit == "mjy":
        tSZ = convert_units(freq, tSZ, cmb2mjy=True)
    elif unit == "cmb":
        None
    elif unit == "rj":
        tSZ = convert_units(freq, tSZ, cmb2rj=True)
    else:
        print("Waring: Unknown unit! Output will be in K_CMB")

    return(tSZ)
  

def simulate_kSZ(freq, nside_out = 4096, lmax = None, beam_FWHM = None, template = "SO", unit = "cmb"):

    """
    Function which computes kSZ maps at different frequencies and different nside. 
    
    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz.
    nside_out: float
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float
        Maximum value of the multipolemoment at which the atmospheric power spectrum
        wil be computed. Default: 3*nside_out-1            
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
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.
        
    Returns
    -------
    kSZ: float array
        Healpy all-sky map contaning the kSZ map at a given frequency. 
    """

    #Read data
    if template != "SO" and template != "CITA" and template != "Sehgal":
        print("Waring: Unknown template requested! Output will be based on SO sky model")
        template = "SO"

    if template == 'SO':
        
        data_path='/vol/arc3/data1/sz/SO_sky_model/CMB_SZ_maps/'
        file_in = '148_ksz_healpix_nopell_Nside4096_DeltaT_uK.fits'
        kSZ = hp.read_map(data_path + file_in, dtype = np.float32)*1e-6
        
    if template == 'CITA':
        
        data_path = '/vol/arc3/data1/sz/CITA/'
        file_in = 'ksz.fits'
        kSZ = hp.read_map(data_path + file_in, dtype = np.float32)*1e-6
        
    if template == 'Sehgal' : 
              
        data_path='/vol/arc3/data1/sz/Sehgal/'
        file_in = '030_ksz_healpix.fits'  

        #Create compton-y map : 
        kSZ = hp.read_map(data_path + file_in, dtype = np.float32) * convert_units(freq, 1e-6, mjy2cmb=True)
 
    #Re-bin map if necessary
    if hp.get_nside(kSZ) != nside_out:
        kSZ = hp.pixelfunc.ud_grade(kSZ, nside_out = nside_out)

    #Smooth map if necessary
    if beam_FWHM is not None:
        print("begin smoothing")
        kSZ = hp.sphtfunc.smoothing(kSZ, iter = 0, lmax = lmax, fwhm = beam_FWHM/60*np.pi/180)

    #Convert units if necessary
    if unit == "mjy":
        kSZ = convert_units(freq, kSZ, cmb2mjy=True)
    elif unit == "cmb":
        None
    elif unit == "rj":
        kSZ = convert_units(freq, kSZ, cmb2rj=True)
    else:
        print("Waring: Unknown unit! Output will be in K_CMB")

    return(kSZ)


def simulate_white_noise(freq, noise_level, nside_out = 4096, unit_noise = 1, arcmin = True, unit = 'cmb'): 
        
    """
    Function which create a White noise map for a given noise/arcmin or noise/radians. 
    By default the code expect the noise level to be given in microK_CMB/arcmin

    Parameters
    ----------

    freq: float or float array
        Frequency of the output map in Hz.
    noise_level : float 
        noise level desired in any units by radians or arcmin.
    nside_out: float
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    unit_noise : float 
        resolution of the noise, for exemple 1' or 1 radians. 
    arcmin : bool 
        if true mean that the noise is given in /arcmin. 
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'mjy'.
        
    Returns
    -------
    array
        Array contaning the Variarion of intensity produced by tSZ over the fequencies. 

    """
    
    #Compute the average noise in each pixel : 
    sigma_noise = (noise_level * unit_noise)/px_size(nside_out,arcmin)
        
    #Create White-noise map :
    npix = hp.pixelfunc.nside2npix(nside_out) #Compute the number of pixels
    noise_map =  np.random.normal(0, sigma_noise, npix)*1e-6 #Random normal distribution centered over the desired noise

    #Convert units if necessary
    if unit == "cmb":
        None
    elif unit == "mjy":
        noise_map = convert_units(freq, noise_map, cmb2mjy=True)
    elif unit == "rj":
        noise_map = convert_units(freq, noise_map, cmb2rj=True)
    else:
        print("Waring: Unknown unit! Output will be in K_CMB")
    
    return(np.float32(noise_map))


def simulate_atmosphere(freq, nside_out = 4096, lmax = None, beam_FWHM = None, unit = "cmb"):

    '''Computes an all-sky atmospheric noise map at a given frequency and nside based on 
    the SO noise model presented by the SO Collaboration (2019) and using the model parameters
    provided by Choi et al. (2019). 

    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz. Must be a valid SO or CCAT-prime central
        band frequency, i.e. 27, 29, 93, 145, 225, 279, 220, 280, 350, 405, or 860 GHz.
    nside_out: float
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float
        Maximum value of the multipolemoment at which the atmospheric power spectrum
        wil be computed. Default: 3*nside_out-1    
    beam_FWHM: bool, optional
        If set, the output will be convolved with a gaussian. The FWHM of the Gaussian
        in units of arcmin is given by the provided value. Default: None
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.

    Returns
    -------
    noise_map: float array
        Healpix all-sky map of the atmospheric noise at the specified frequency.
    '''

    if lmax is None:
        lmax = 3*nside_out-1

    #Define frequencies and noise characteristics of the SO and CCAT-prime
    nu = np.array([27, 29, 93, 145, 225, 279, 220, 280, 350, 405, 860])*1e9
    N_white = np.array([2.351167E-04, 6.414462E-05, 2.975014E-06, 3.458125E-06, 2.011171E-05, 1.272230E-04, 1.8e-5, 6.4e-5, 9.3e-4, 1.2e-2, 2.8e4])
    N_red = np.array([9.248376E-06, 2.236786E-06, 2.622661E-05, 3.298953E-04, 1.462620E-02, 8.579506E-02, 1.6e-2, 1.1e-1, 2.7e0, 1.7e1, 6.1e6])
    ell_knee = 1000 
    alpha_knee = -3.5
    
    #Check if input frequency is a valid SO or CCAT-prime central band frequency
    index = freq == nu
    if np.sum(index) == 0:

        print("Warning: Input frequency is not a valid SO or CCAT-prime band. Try 27, 29, 93, 145, 225, 279, 220, 280, 350, 405, or 860 GHz")
        return(None)

    else:

        #Compute power spectrum
        ell = np.linspace(1,lmax,lmax)
        Cl = N_red[index] * (ell/ell_knee)**alpha_knee + N_white[index]

        #Create all-sky map
        noise_map = hp.sphtfunc.synfast(Cl, nside_out, lmax=lmax)/1e6

        #Convert units if necessary
        if unit == "mjy":
            noise_map = convert_units(freq, noise_map, cmb2mjy=True)
        elif unit == "cmb":
            None
        elif unit == "rj":
            noise_map = convert_units(freq, noise_map, cmb2rj=True)
        else:
            print("Waring: Unknown unit! Output will be in K_CMB")

        #Return output
        return(np.float32(noise_map))


def ccatp_sky_model(freq, sensitivity = None, components = "all", red_noise = False, cl_file = None, lensed = True, out_file = None, nside_out = 4096, lmax = None, beam_FWHM = None, template = "SO", unit = "cmb"):

    '''Computes an all-sky map of the simulated microwave sky at the specified frequency. 
    The CCAT-prime sky model uses all-sky templates in the healpix format and includes the
    most important galactic foregrounds, extragalactic backgrounds, the tSZ and kSZ effects 
    of galaxy clusters, and instrumental white/red noise. The individual components can be
    selectively switched on and of.

    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz.
    sensitivity: float
        sensitivity of the instrument in units of micro K_CMB - arcmin. Default: None
    components: string or list of stings
        List containing the names of the components to be included in the sky model.
	    Any combination of the following individual components is possible: "gal_synchrotron", 
	    "gal_dust", "gal_freefree", "gal_ame", "cib", "radio_ps", "cmb", "tsz", "ksz". All 
	    components are used if components is set to "all". Default: "all"
    red_noise: bool
        If True, a realistic white noise + red noise atmospheric model is added to the
	sky model in case the imput frequency is a valid SO or CCAT-prime central
        band frequency, i.e. 27, 29, 93, 145, 225, 279, 220, 280, 350, 405, or 860 GHz.
	    Default: False
    cl_file: string
        name of a file containing a CMB power spectrum computed e.g. with CAMP. The
	    first column of the file has to correspond to ell, the second column to 
	    Cl ell(ell+1)/2pi. If set, a random realization of the CMB based on the provided
	    powerspectrum will be added to the data.
    lensed: bool
        If True, lensed SO, Sehgal and CITA CMB maps will be used. Default: True 
    out_file: string
        If set, the results will be written as a healpy .fits file of the given name. 
    nside_out: float
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float
        Maximum value of the multipolemoment. Default: 3*nside_out-1            
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
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.

    Returns
    -------
    allsky_map: float array
        Healpix all-sky map of the microwave sky at the specified frequency and nside.
    '''

    if lmax is None:
        lmax = int(3*nside_out-1)

    if components == "all":
        components = ["gal_synchrotron", "gal_dust", "gal_freefree", "gal_ame", "cib", "radio_ps", "cmb", "tsz", "ksz"]

    #build sky model from the individual components
    allsky_map = np.zeros(hp.nside2npix(4096))

    if ("gal_synchrotron" in components) or ("gal_dust" in components) or ("gal_freefree" in components) or ("gal_ame" in components):
        print("Computing Galactic foregounds...")
        allsky_map += simulate_gal_foregrounds(freq, components, nside_out = nside_out, lmax = lmax, beam_FWHM = None, unit = unit)
        print("Galactic foregounds complete.")

    if "cib" in components:
        print("Computing CIB...")
        allsky_map += simulate_cib(freq, nside_out = nside_out, beam_FWHM = None, template = template, unit = unit)
        print("CIB complete.")

    if "radio_ps" in components:
        print("Computing radio PS...")
        allsky_map += simulate_radio_ps(freq, nside_out = nside_out, lmax = lmax, beam_FWHM = None, template = template, unit = unit)
        print("Radio PS complete.")

    if "cmb" in components:
        print("Computing CMB...")
        allsky_map += simulate_cmb(freq, cl_file = cl_file, lensed = lensed, nside_out = nside_out, lmax = lmax, beam_FWHM = None, template = template, unit = unit)
        print("CMB complete.")

    if "tsz" in components:
        print("Computing tSZ...")
        allsky_map += simulate_tSZ(freq, nside_out = nside_out, lmax = lmax, beam_FWHM = None, template = template, unit = unit)
        print("tSZ complete.")

    if "ksz" in components:
        print("Computing kSZ...")
        allsky_map += simulate_kSZ(freq, nside_out = nside_out, lmax = lmax, beam_FWHM = None, template = template, unit = unit)
        print("kSZ complete.")

    #Smooth map if necessary
    if beam_FWHM is not None:
        print("Begin smoothing...")
        allsky_map = hp.sphtfunc.smoothing(allsky_map, iter = 1, lmax = lmax, fwhm = beam_FWHM/60*np.pi/180)
        print("Smoothing complete.")

    #Add noise if necessary
    if red_noise is True:
        if sensitivity is True:
            print("Warning! You request to apply both red noise + white noise and white noise. The White noise sensitivity parameter will be overwritten to None. If this is not what you want than check the settings and re-run.")
            sensitivity = None
        print("Computing red noise...")        
        allsky_map += simulate_atmosphere(freq, nside_out = nside_out, lmax = lmax, beam_FWHM = None, unit = unit)
        print("Red noise complete.")

    if sensitivity is not None and red_noise is False:
        print("Computing white noise...")        
        allsky_map += simulate_white_noise(freq, sensitivity, nside_out = nside_out, unit = unit)
        print("White noise complete.")

    if out_file is not None:
        print("The result will be written to " + out_file)
        hp.fitsfunc.write_map(out_file, allsky_map, overwrite = True, dtype = np.float32)
        allsky_map = None

    print("Done!")

    return(allsky_map)
