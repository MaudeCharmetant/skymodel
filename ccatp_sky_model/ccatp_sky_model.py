import numpy as np 
import healpy as hp
from astropy.io import fits, ascii
from astropy import constants as cst
from astropy.coordinates import SkyCoord
from astropy import units as u
import pysm
from tqdm import tqdm 
from pysm.nominal import models
import os.path

os_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'masks/')
data_path = '/vol/arc3/data1/sz/CCATp_sky_model/templates/'

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
		print('Not sure which units are given and what should be returned.')

	converted_signal = conversion * values

	return(converted_signal)


def px_size(nside,arcmin=True):   
            
    '''
    Computes the size of Healpy pixels in function of the nside of the map.  

    Parameters
    ----------
    
    nside: int
        Nside, number of separation of the healpy blocks. 
    arcmin: bool, optional
        if True return the size in arcmin, otherwise in radians. 
        
    Returns
    -------
    size: float
        The size of one healpy pixel. 
    '''
    
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
		only be drawn in areas that are not masked. If mask is set
		to 'advACT', 'SPT', 'Dust', or 'NVSS', the respective 
		survey masks will be used. Default: None
	radec: bool, optional
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
		if mask == 'advACT':
			mask = hp.read_map(os_path + 'Adv_ACT_survey_mask.fits', dtype = np.int16)  
		elif mask == 'SPT':
			mask = hp.read_map(os_path + 'SPT-SZ_survey_mask.fits', dtype = np.int16)
		elif mask == 'Dust':
			mask = hp.read_map(os_path + 'galactic_dust_mask.fits', dtype = np.int16)
		elif mask == 'NVSS':
			mask = hp.read_map(os_path + 'galactic_dust_nvss_mask.fits', dtype = np.int16)
	
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
		print('output will be fk5 coordinates (RA, DEC) for the equinox J2000')
		c = SkyCoord(phi, theta, frame='galactic', unit='deg')
		c_fk5 = c.transform_to('fk5')
		phi = c_fk5.ra.degree
		theta = c_fk5.dec.degree
	else:
		print('output will be galactic coordinates (glon, glat)')

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
			maps[i,:,:] = hp.visufunc.gnomview(allsky_map, coord = ['G', 'C'], rot=[RA[i],DEC[i]], reso = pixel_size, xsize = npix, return_projected_map = True, no_plot=True) 

	return(maps)


def sigmoid_filter(l_0, d_l, lmax):
    
    '''Computes a sigmoid-shaped filer high-pass window in spherical harmonic space. 
    
    Parameters
    ----------
    l_0: float
        Determines the spatial scale above which the sky signal will be filtered.
    d_l: float
        Determines the slope of the widow.
    lmax: float
        Maximum value of ell to be used
        
    Returns
    -------
    window: float array
        Computed window in spherical harmonic space
    '''    
    
    ell = np.arange(lmax)
    window = 1/(1+np.exp((-ell+l_0)/d_l))
    
    return(window)


def return_mask(survey, nside_out = 256, coord = 'G'):

    '''Returns the specified all-sky survey map. 
    
    Parameters
    ----------
    survey: sting
        Defines which survey mask will be returned. The options are 'advACT', 'SPT',
        'Dust', and 'NVSS'. 
    nside_out: float, optional
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 256
    coord: sting, optional
        Defines the coordinate system of the output mask. 'G' --> Galactic, 
        'E' --> Ecliptic, 'C' --> Equatorial. Default: 'G'
    
    Returns
    -------
    mask: float array
        Healpix all-sky mask.
    '''    
    
    #read mask
    if survey == 'advACT':
        mask = hp.read_map(os_path + 'Adv_ACT_survey_mask.fits', dtype = np.int16)  
    elif survey == 'SPT':
        mask = hp.read_map(os_path + 'SPT-SZ_survey_mask.fits', dtype = np.int16)
    elif survey == 'Dust':
        mask = hp.read_map(os_path + 'galactic_dust_mask.fits', dtype = np.int16)
    elif survey == 'NVSS':
        mask = hp.read_map(os_path + 'galactic_dust_nvss_mask.fits', dtype = np.int16)

    #change coordinate system if necessary
    if (coord == 'E') or (coord == 'C'):
        mask = hp.ud_grade(mask, nside_out = 2048)
        r = hp.Rotator(coord = ('G', coord))
        mask = hp.Rotator.rotate_map_pixel(r, mask)
        mask = hp.ud_grade(mask, nside_out = 256)
        mask[mask != 0] = 1

    #Re-bin map if necessary
    if nside_out != hp.get_nside(mask):
        mask = hp.ud_grade(mask, nside_out = nside_out)
        
    return(np.int16(mask))


def simulate_gal_foregrounds(freq, components = 'all', nside_out = 4096, lmax = None, beam_FWHM = None, intrinsic_FWHM = 10, unit = 'cmb'):

    '''Computes an all-sky galactic foregrounds noise map at a given frequency and nside using 
    the Python Sky model (PySM, Thorne et al. 2017), which is build from Planck data. 

    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz. 
    components: list of stings, optional
        List of components to be included in the galactic foreground model. Possible 
        components are 'gal_synchrotron', 'gal_dust', 'gal_freefree', and 'gal_ame'.
        Alternatively, all components are included if the variable is set to 'all'.
        Default: 'all'	
    nside_out: float, optional
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float, optional
        Maximum value of the multipolemoment at which the atmospheric power spectrum
        wil be computed. Default: 3*nside_out-1    
    beam_FWHM: bool, optional
        If set, the output will be convolved with a gaussian. The FWHM of the Gaussian
        in units of arcmin is given by the provided value. Default: None
    intrinsic_FWHM: float, optional
        Determines the with of a gaussian that is always applied to the Galactic foreground
        maps after they have been upgraded to a higher nside. Default: 10
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.

    Returns
    -------
    foregrounds: float array
        Healpix all-sky map of the galactic foregrounds at the specified frequency.
    '''

    if lmax is None:
        lmax = int(3*nside_out-1)

    if components == 'all':
        components = ['gal_synchrotron', 'gal_dust', 'gal_freefree', 'gal_ame']	

    #Define foreground model
    sky_config = {
        'synchrotron' : models('s1', 512),
        'dust' : models('d1', 512),
        'freefree' : models('f1', 512),
        'ame' : models('a1', 512),
    }

    #Initialise Sky 
    sky = pysm.Sky(sky_config)

    #Compute component maps
    foregrounds = np.zeros(hp.nside2npix(512))
    rj2mjy_factor = convert_units(freq, 1e-6, rj2mjy = True)
    
    if 'gal_synchrotron' in components:
        foregrounds += sky.synchrotron(freq/1e9)[0,:] * rj2mjy_factor

    if 'gal_dust' in components:
        foregrounds += sky.dust(freq/1e9)[0,:] * rj2mjy_factor
          
    if 'gal_freefree' in components:
        foregrounds += sky.freefree(freq/1e9)[0,:] * rj2mjy_factor
        
    if 'gal_ame' in components:
        foregrounds += sky.ame(freq/1e9)[0,:] * rj2mjy_factor

    #Define smoothing kernal, upgrade, and smooth map
    fwhm = intrinsic_FWHM
    if beam_FWHM is not None:
        if beam_FWHM > intrinsic_FWHM:
            fwhm = np.sqrt(intrinsic_FWHM**2 + beam_FWHM**2)
	    
    foregrounds = hp.pixelfunc.ud_grade(foregrounds, nside_out = nside_out)
    foregrounds = hp.sphtfunc.smoothing(foregrounds, iter = 0, fwhm = fwhm/60*np.pi/180, lmax = lmax)
    
    #Convert units if necessary
    if unit == 'mjy':
        None
    elif unit == 'cmb':
        foregrounds = convert_units(freq, foregrounds, mjy2cmb=True)
    elif unit == 'rj':
        foregrounds = convert_units(freq, foregrounds, mjy2rj=True)
    else:
        foregrounds = convert_units(freq, foregrounds, mjy2cmb=True)
        print('Waring: Unknown unit! Output will be in K_CMB')
        
    return(foregrounds)


def simulate_cib(freq, nside_out = 4096, lmax = None, beam_FWHM = None, template = 'WebSky', unit = 'cmb'):

    '''Computes an all-sky CIB map at a given frequency and nside based on . 

    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz. If freq=-1 return a frequency independant CIB. 
    nside_out: float, optional
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float, optional
        Maximum value of the multipolemoment at which the atmospheric power spectrum
        wil be computed. Default: 3*nside_out-1  	
    beam_FWHM: bool, optional
        If set, the output will be convolved with a gaussian. The FWHM of the Gaussian
        in units of arcmin is given by the provided value. Default: None
    template: bool, optional
        Determines the all-sky foregrounds templates to be used to build the sky model.
        If 'Sehgal' is chosen, simulations by Sehgal et al. (2010) are used.
        If 'SO' is chosen, the Simons Observatory sky model provided by Colin Hill and 
        based on the simulations by Sehgal et al. (2010) is used. If 'WebSky' is chosen,
        the used templates will be based on the WebSky Extragalactic CMB Mocks provided 
        by CITA. If 'SO_reproduced' is chosen, the SO sky model is reproduced directly 
        from the Sehgal et al. (2010) data. Default: 'WebSky'
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.

    Returns
    -------
    cib: float array
        Healpix all-sky map of the CIB mission that the specified frequency.
    '''

    if lmax is None:
        lmax = int(3*nside_out-1)

    if freq == -1:
        freq = 353e9
        y_CIB = True
    else :	
        y_CIB = False
	
    #Load all-sky parameter value maps
    if template != 'SO' and template != 'WebSky' and template != 'Sehgal' and template != 'SO_reproduced':
        print('Waring: Unknown template requested! Output will be based on WebSky sky model')
        template = 'WebSky'

    if template == 'SO':
        A = hp.fitsfunc.read_map(data_path + 'CIB/SO_CIB_A_DUST_4096.fits', dtype = np.float32)    
        T = hp.fitsfunc.read_map(data_path + 'CIB/SO_CIB_T_DUST_4096.fits', dtype = np.float32)
        beta = hp.fitsfunc.read_map(data_path + 'CIB/SO_CIB_beta_DUST_4096.fits', dtype = np.float32)
        f_0 = 353e9

    elif (template == 'Sehgal') or (template == 'SO_reproduced'):
        A = hp.fitsfunc.read_map(data_path + 'CIB/Sehgal_CIB_A_DUST_8192.fits', dtype = np.float32)
        T = hp.fitsfunc.read_map(data_path + 'CIB/Sehgal_CIB_T_DUST_8192.fits', dtype = np.float32)
        beta = hp.fitsfunc.read_map(data_path + 'CIB/Sehgal_CIB_beta_DUST_8192.fits', dtype = np.float32)
        f_0 = 350e9

        if template == 'SO_reproduced':
            A *= 0.75
	
    elif template == 'WebSky':
        A = hp.fitsfunc.read_map(data_path + 'CIB/CITA_CIB_A_DUST_4096.fits', dtype = np.float32)    
        T = hp.fitsfunc.read_map(data_path + 'CIB/CITA_CIB_T_DUST_4096.fits', dtype = np.float32)
        beta = hp.fitsfunc.read_map(data_path + 'CIB/CITA_CIB_beta_DUST_4096.fits', dtype = np.float32)
        f_0 = 353e9

    #Compute CIB brightness at given frequency
    cib = A * (freq/f_0)**(3.+beta) * (np.exp(h*f_0/k_B/T)-1) / (np.exp(h*freq/k_B/T)-1)
    x_nu = np.array((h*freq)/(k_B*T_CMB))   
    tSZ_SED = ((x_nu*(np.exp(x_nu)+1)/(np.exp(x_nu)-1))-4)

    if y_CIB != True: 
        cib = cib
    else:
        cib =  cib / tSZ_SED 
	
    del A, T, beta
    
    #Re-bin map if necessary
    if hp.get_nside(cib) != nside_out:
        cib = hp.pixelfunc.ud_grade(cib, nside_out = nside_out)

    #Smooth map if necessary
    if beam_FWHM is not None:
        print('begin smoothing')
        cib = hp.sphtfunc.smoothing(cib, iter = 0, lmax = lmax, fwhm = beam_FWHM/60*np.pi/180)

    #Convert units if necessary
    #Get the frequency independent y-map : 
    if y_CIB != True:	  
        #Convert units if necessary
        if unit == 'mjy':
            None
        elif unit == 'cmb':
            cib = convert_units(freq, cib, mjy2cmb=True)
        elif unit == 'rj':
            cib = convert_units(freq, cib, mjy2rj=True)
        else:
            cib = convert_units(freq, cib, mjy2cmb=True)
            print('Waring: Unknown unit! Output will be in K_CMB')
    else: 
            cib = convert_units(freq, cib, mjy2cmb=True)	

    #Return output
    return(np.float32(cib))


def simulate_radio_ps(freq, nside_out = 4096, lmax = None, beam_FWHM = None, template = 'WebSky', unit = 'cmb'):

    '''Computes an all-sky radio point source map at a given frequency and nside based on 
    the simulations provided by Sehgal et al. (2010), which have been recalibrated by the
    SO collaboration. The original simulations by Sehgal et al. are modelled by a curved 
    power law and can be extrapolated to frequencies beyoned 350 GHz, while the sources in 
    the SO simulations show more complex spectra due to ringing artifacts caused by the 
    processing of the maps with alter_alm(). Due to these complex SEDs, bilinear interpolation
    is applied for input frequencies between 27 and 353 GHz. For higher frequencies, a null
    map is returned.

    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz.
    nside_out: float, optional
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float, optional
        Maximum value of the multipolemoment at which the atmospheric power spectrum
        wil be computed. Default: 3*nside_out-1            
    beam_FWHM: bool, optional
        If set, the output will be convolved with a gaussian. The FWHM of the Gaussian
        in units of arcmin is given by the provided value. Default: None
    template: bool, optional
        Determines the all-sky foregrounds templates to be used to build the sky model.
        If 'Sehgal' is chosen, simulations by Sehgal et al. (2010) are used.
        If 'SO' is chosen, the Simons Observatory sky model provided by Colin Hill and 
        based on the simulations by Sehgal et al. (2010) is used. If 'WebSky' is chosen,
        the used templates will be based on the WebSky Extragalactic CMB Mocks provided 
        by CITA. If 'SO_reproduced' is chosen, the SO sky model is reproduced directly 
        from the Sehgal et al. (2010) data. Default: 'WebSky'
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
    
    if template != 'SO' and template != 'WebSky' and template != 'Sehgal' and template != 'SO_reproduced':
        print('Waring: Unknown template requested! Output will be based on SO sky model')
        template = 'WebSky'	
	
    if template == 'SO':

        npix = hp.pixelfunc.nside2npix(4096)
        radio_ps = np.zeros(npix)
        
        #Define frequencies
        nu = np.array([27,30,39,44,70,93,100,143,145,217,225,280,353])*1e9
        nu_names = ['027','030','039','044','070','093','100','143','145','217','225','280','353']

        #Interpolate data points
        if freq > np.max(nu):
            print('Warning: Input frequency lies beyoned the 353 GHz. Since higher frequencies are not constraint by simulations, the data will be 0.')
        else:
		
            #Read data files
            data = np.zeros((len(nu), npix), dtype = np.float32)

            for i in np.arange(len(nu)):
                file_name = data_path + 'radio_ps/' + nu_names[i] + '_rad_pts_healpix_nopell_Nside4096_DeltaT_uK_fluxcut148_7mJy_lininterp.fits'
                data[i,:] = hp.fitsfunc.read_map(file_name, dtype = np.float32) * convert_units(nu[i], 1e-6, cmb2mjy=True)		
		
            for i in tqdm(np.arange(npix)):
                radio_ps[i] = np.interp(freq, nu, data[:,i])
		
            del data		

    elif (template == 'Sehgal') or (template == 'SO_reproduced'):

        if template == 'SO_reproduced':

            file_name = data_path + 'radio_ps/148_rad_pts_healpix.fits'
            map_148 = hp.fitsfunc.read_map(file_name, dtype = np.float32)/1e6
            
            solid_angle = 4*np.pi / hp.nside2npix(8192)        
            threshold = 7e-3 / solid_angle / 1e6 #defines 7 mJy source cut at 148 GHz in units of MJy/sr for nside = 8192
            flagged = map_148 > threshold
            del map_148

        #Read data files
        file_name = data_path + 'radio_ps/030_rad_pts_healpix.fits'
        map_30 = hp.fitsfunc.read_map(file_name, dtype = np.float32)/1e6

        file_name = data_path + 'radio_ps/Sehgal_radio_ps_spectral_index_8192.fits'
        spectral_index = hp.fitsfunc.read_map(file_name, dtype = np.float32)

        file_name = data_path + 'radio_ps/Sehgal_radio_ps_f_crit_8192.fits'
        crit_freq = hp.fitsfunc.read_map(file_name, dtype = np.float32)

        #Extrapolate fluxes
        radio_ps = map_30*(freq/30e9)**spectral_index * np.exp(-freq/crit_freq)

        if template == 'SO_reproduced':
            radio_ps[flagged] = 0

        del map_30, spectral_index, crit_freq	
				
    elif template == 'WebSky':
        print('Warning: No radio PS template provided by the WebSky simulations')
        npix = hp.pixelfunc.nside2npix(4096)
        radio_ps = np.zeros(npix)

    #Re-bin map if necessary
    if hp.get_nside(radio_ps) != nside_out:
        radio_ps = hp.pixelfunc.ud_grade(radio_ps, nside_out = nside_out)

    #Smooth map if necessary
    if beam_FWHM is not None:
        print('begin smoothing')
        radio_ps = hp.sphtfunc.smoothing(radio_ps, iter = 0, lmax = lmax, fwhm = beam_FWHM/60*np.pi/180)

    #Convert units if necessary
    if unit == 'mjy':
        None
    elif unit == 'cmb':
        radio_ps = convert_units(freq, radio_ps, mjy2cmb=True)
    elif unit == 'rj':
        radio_ps = convert_units(freq, radio_ps, mjy2rj=True)
    else:
        radio_ps = convert_units(freq, radio_ps, mjy2cmb=True)
        print('Waring: Unknown unit! Output will be in K_CMB')

    #Return output
    return(np.float32(radio_ps))


def simulate_cmb(freq, cl_file = None, lensed = True, nside_out = 4096, lmax = None, beam_FWHM = None, template = 'WebSky', unit = 'cmb'): 
         
    '''
    Function that computes a CMB map from a power spectrum.

    Parameters
    ----------

    freq: float or float array
        Frequency of the output map in Hz. 
    cl_file: str or array, optional 
        Name of the .dat contaning the values of the power spectrum given by CAMB.
        Or array containing the power spectrum to generate random maps in Kelvin.
        Default: None	
    lensed: bool, optional
    	If True select the lensed CMB map when possible. This is only possible for 'WebSky' and 'Sehgal.'
        Default: True
    nside_out: float, optional
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float, optional
        Maximum value of the multipolemoment at which the atmospheric power spectrum
        wil be computed. Default: 3*nside_out-1            
    beam_FWHM: bool, optional
        If set, the output will be convolved with a gaussian. The FWHM of the Gaussian
        in units of arcmin is given by the provided value. Default: None
    template: bool, optional
        Determines the all-sky foregrounds templates to be used to build the sky model.
        If 'Sehgal' is chosen, simulations by Sehgal et al. (2010) are used.
        If 'SO' is chosen, the Simons Observatory sky model provided by Colin Hill and 
        based on the simulations by Sehgal et al. (2010) is used. If 'WebSky' is chosen,
        the used templates will be based on the WebSky Extragalactic CMB Mocks provided 
        by CITA. If 'SO_reproduced' is chosen, the SO sky model is reproduced directly 
        from the Sehgal et al. (2010) data. Default: 'WebSky'
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.
        
    Returns
    -------
    float array
        Healpix allsky map contaning the CMB map at a given frequency. 

    '''
                
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

        if template != 'SO' and template != 'WebSky' and template != 'Sehgal' and template != 'SO_reproduced':
            print('Waring: Unknown template requested! Output will be based on WebSky sky model')
            template = 'WebSky'			
		
        if template == 'WebSky':
        	
            if lensed == True:
			
                file_name = 'CMB/CMB_lensed_CITA_mK.fits'

            else:
		
                file_name = 'CMB/CMB_unlensed_CITA_mK.fits'
            
            CMB = hp.read_map(data_path + file_name, dtype = np.float32)/1e6     
    
        elif template == 'SO': 

            file_name = 'CMB/Sehgalsimparams_healpix_4096_KappaeffLSStoCMBfullsky_phi_SimLens_Tsynfastnopell_fast_lmax8000_nside4096_interp2.5_method1_1_lensed_map.fits'

            CMB = hp.read_map(data_path + file_name, dtype = np.float32)/1e6 		
		
        elif (template == 'Sehgal') or (template == 'SO_reproduced'):
			
            if lensed == True:

                file_name = 'CMB/030_lensedcmb_healpix.fits'

            else: 

                file_name = 'CMB/030_unlensedcmb_healpix.fits'

            CMB = hp.read_map(data_path + file_name, dtype = np.float32) * convert_units(30e9, 1e-6, mjy2cmb=True)

    #Re-bin map if necessary
    if hp.get_nside(CMB) != nside_out:
        CMB = hp.pixelfunc.ud_grade(CMB, nside_out = nside_out)

    #Smooth map if necessary
    if beam_FWHM is not None:
        print('begin smoothing')
        CMB = hp.sphtfunc.smoothing(CMB, iter = 0, lmax = lmax, fwhm = beam_FWHM/60*np.pi/180)	
	
    #Convert units if necessary : 
    if unit == 'mjy':
        CMB = convert_units(freq, CMB, cmb2mjy=True)
    elif unit == 'cmb':
        None
    elif unit == 'rj':
        CMB = convert_units(freq, CMB, cmb2rj=True)
    else:	
        print('Waring: Unknown unit! Output will be in K_CMB')   

    return(CMB)


def simulate_tSZ(freq, nside_out = 4096, lmax = None, beam_FWHM = None, template = 'WebSky', unit = 'cmb'):
    
    '''
    Function which compute tSZ maps at different frequencies and different nside. 
    
    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz. If freq = -1 a unitless y-map is returned. 
    nside_out: float, optional
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float, optional
        Maximum value of the multipolemoment at which the atmospheric power spectrum
        wil be computed. Default: 3*nside_out-1            
    beam_FWHM: bool, optional
        If set, the output will be convolved with a gaussian. The FWHM of the Gaussian
        in units of arcmin is given by the provided value. Default: None
    template: bool, optional
        Determines the all-sky foregrounds templates to be used to build the sky model.
        If 'Sehgal' is chosen, simulations by Sehgal et al. (2010) are used.
        If 'SO' is chosen, the Simons Observatory sky model provided by Colin Hill and 
        based on the simulations by Sehgal et al. (2010) is used. If 'WebSky' is chosen,
        the used templates will be based on the WebSky Extragalactic CMB Mocks provided 
        by CITA. If 'SO_reproduced' is chosen, the SO sky model is reproduced directly 
        from the Sehgal et al. (2010) data. Default: 'WebSky'
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.
        
    Returns
    -------
    tSZ: float array
        Healpix allsky map contaning the tSZ map at a given frequency. 
    '''

    #Read data

    if template != 'SO' and template != 'WebSky' and template != 'Sehgal' and template != 'SO_reproduced':
        print('Waring: Unknown template requested! Output will be based on WebSky sky model')
        template = 'WebSky'

    if template == 'SO':
        
        file_in ='tSZ/tSZ_skymap_healpix_nopell_Nside4096_y_tSZrescale0p75.fits'
        y_map = hp.read_map(data_path + file_in, dtype = np.float32)

        
    elif template == 'WebSky': 
        
        file_in = 'tSZ/tsz.fits'
        y_map = hp.read_map(data_path + file_in, dtype = np.float32)        
        
    elif (template == 'Sehgal') or (template == 'SO_reproduced'): 
        
        file_in='tSZ/030_tsz_healpix.fits'
        
        #Create Compton y-map : 
        tSZ_30GHz = hp.read_map(data_path + file_in, dtype = np.float32) * convert_units(30e9, 1e-6, mjy2cmb=True)
           
        x_nu = np.array((h*30e9)/(k_B*T_CMB))    
        tSZ_SED = ((x_nu*(np.exp(x_nu)+1)/(np.exp(x_nu)-1))-4)

        y_map = tSZ_30GHz / tSZ_SED / T_CMB

        if template == 'SO_reproduced':
            y_map *= 0.75

    #Get tSZ at different frequencies : 
    x_nu = np.array((h*freq)/(k_B*T_CMB))    
    t_SZ_SED = ((x_nu*(np.exp(x_nu)+1)/(np.exp(x_nu)-1))-4)

    if freq != -1: 
        tSZ = (y_map * t_SZ_SED * T_CMB)
    else:
        tSZ =  y_map
        
    #Re-bin map if necessary
    if hp.get_nside(tSZ) != nside_out:
        tSZ = hp.pixelfunc.ud_grade(tSZ, nside_out = nside_out)

    #Smooth map if necessary
    if beam_FWHM is not None:
        print('begin smoothing')
        tSZ = hp.sphtfunc.smoothing(tSZ, iter = 0, lmax = lmax, fwhm = beam_FWHM/60*np.pi/180)

    if freq != -1: 	
        #Convert units if necessary
        if unit == 'mjy':
            tSZ = convert_units(freq, tSZ, cmb2mjy=True)
        elif unit == 'cmb':
            None
        elif unit == 'rj':
            tSZ = convert_units(freq, tSZ, cmb2rj=True)
        else:
            print('Waring: Unknown unit! Output will be in K_CMB')

    return(tSZ)
  

def simulate_kSZ(freq, nside_out = 4096, lmax = None, beam_FWHM = None, template = 'WebSky', unit = 'cmb'):

    '''
    Function which computes kSZ maps at different frequencies and different nside. 
    
    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz. 
    nside_out: float, optional
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float, optional
        Maximum value of the multipolemoment at which the atmospheric power spectrum
        wil be computed. Default: 3*nside_out-1            
    beam_FWHM: bool, optional
        If set, the output will be convolved with a gaussian. The FWHM of the Gaussian
        in units of arcmin is given by the provided value. Default: None
    template: bool, optional
        Determines the all-sky foregrounds templates to be used to build the sky model.
        If 'Sehgal' is chosen, simulations by Sehgal et al. (2010) are used.
        If 'SO' is chosen, the Simons Observatory sky model provided by Colin Hill and 
        based on the simulations by Sehgal et al. (2010) is used. If 'WebSky' is chosen,
        the used templates will be based on the WebSky Extragalactic CMB Mocks provided 
        by CITA. If 'SO_reproduced' is chosen, the SO sky model is reproduced directly 
        from the Sehgal et al. (2010) data. Default: 'WebSky'
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.
        
    Returns
    -------
    kSZ: float array
        Healpy all-sky map contaning the kSZ map at a given frequency. 
    '''

    #Read data
    if template != 'SO' and template != 'WebSky' and template != 'Sehgal' and template != 'SO_reproduced':
        print('Waring: Unknown template requested! Output will be based on SO sky model')
        template = 'WebSky'

    if template == 'SO':
        
        file_in = 'kSZ/148_ksz_healpix_nopell_Nside4096_DeltaT_uK.fits'
        kSZ = hp.read_map(data_path + file_in, dtype = np.float32)/1e6
        
    elif template == 'WebSky':
        
        file_in = 'kSZ/ksz.fits'
        kSZ = hp.read_map(data_path + file_in, dtype = np.float32)/1e6
        
    elif (template == 'Sehgal') or (template == 'SO_reproduced'):
              
        file_in = 'kSZ/030_ksz_healpix.fits'  
        kSZ = hp.read_map(data_path + file_in, dtype = np.float32) * convert_units(30e9, 1e-6, mjy2cmb=True)
 
    #Re-bin map if necessary
    if hp.get_nside(kSZ) != nside_out:
        kSZ = hp.pixelfunc.ud_grade(kSZ, nside_out = nside_out)

    #Smooth map if necessary
    if beam_FWHM is not None:
        print('begin smoothing')
        kSZ = hp.sphtfunc.smoothing(kSZ, iter = 0, lmax = lmax, fwhm = beam_FWHM/60*np.pi/180)
	
    #Convert units if necessary
    if unit == 'mjy':
        kSZ = convert_units(freq, kSZ, cmb2mjy=True)
    elif unit == 'cmb':
        None
    elif unit == 'rj':
        kSZ = convert_units(freq, kSZ, cmb2rj=True)
    else:
        print('Waring: Unknown unit! Output will be in K_CMB')

    return(kSZ)


def simulate_white_noise(freq, noise_level, nside_out = 4096, unit_noise = 1, arcmin = True, unit = 'cmb'): 
        
    '''
    Function which create a White noise map for a given noise/arcmin or noise/radians. 
    By default the code expect the noise level to be given in microK_CMB/arcmin

    Parameters
    ----------

    freq: float or float array
        Frequency of the output map in Hz.If freq=0 we get a frequency independent map. 
    noise_level: float, optional 
        noise level desired in any units of micro K_CMB by radians or arcmin.
    nside_out: float, optional
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    unit_noise: float, optional 
        resolution of the noise, for exemple 1' or 1 radians. Default: 1
    arcmin: bool, optional 
        if true mean that the noise is given in /arcmin. Default: True
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'mjy'.
        
    Returns
    -------
    array
        Array contaning the white noise map at a given frequency and resolution. 

    '''
    
    #Compute the average noise in each pixel : 
    sigma_noise = (noise_level * unit_noise)/px_size(nside_out,arcmin)
        
    #Create White-noise map :
    npix = hp.pixelfunc.nside2npix(nside_out) #Compute the number of pixels
    noise_map =  np.random.normal(0, sigma_noise, npix)*1e-6 #Random normal distribution centered over the desired noise
	
    #Convert units if necessary
    if unit == 'cmb':
        None
    elif unit == 'mjy':
        noise_map = convert_units(freq, noise_map, cmb2mjy=True)
    elif unit == 'rj':
        noise_map = convert_units(freq, noise_map, cmb2rj=True)
    else:
        print('Waring: Unknown unit! Output will be in K_CMB')
    
    return(np.float32(noise_map))


def simulate_atmosphere(freq, nside_out = 4096, lmax = None, beam_FWHM = None, unit = 'cmb', no_white = False):

    '''Computes an all-sky atmospheric noise map at a given frequency and nside based on 
    the SO noise model presented by the SO Collaboration (2019) and using the model parameters
    provided by Choi et al. (2019). 

    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz. Must be a valid SO or CCAT-prime central
        band frequency, i.e. 27, 39, 93, 145, 225, 279, 220, 280, 350, 405, or 860 GHz.
    nside_out: float, optional
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float, optional
        Maximum value of the multipolemoment at which the atmospheric power spectrum
        wil be computed. Default: 3*nside_out-1    
    beam_FWHM: bool, optional
        If set, the output will be convolved with a gaussian. The FWHM of the Gaussian
        in units of arcmin is given by the provided value. Default: None
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.
    no_white: bool, optional
        If True, only the red noise component is simulated. Default: False

    Returns
    -------
    noise_map: float array
        Healpix all-sky map of the atmospheric noise at the specified frequency.
    '''

    if lmax is None:
        lmax = 3*nside_out-1

    #Define frequencies and noise characteristics of the SO and CCAT-prime.
    #The values of N_white and N_red for the SO correspond to the baseline sensitivities and are obtained using v. 3.1.1. of the SO noise model. 
    nu = np.array([27, 39, 93, 145, 225, 279, 220, 280, 350, 405, 860])*1e9	
    N_white = np.array([4.31747108e-04, 1.07936777e-04, 5.46429934e-06, 8.41194778e-06, 4.21628035e-05, 2.42857748e-04, 1.8e-5, 6.4e-5, 9.3e-4, 1.2e-2, 2.8e4])
    N_red = np.array([7.88002445e-06, 1.90584033e-06, 2.23462248e-05, 2.81085357e-04, 1.24621704e-02, 7.31011726e-02, 1.6e-2, 1.1e-1, 2.7e0, 1.7e1, 6.1e6])
    ell_knee = 1000 
    alpha_knee = -3.5
    
    #Check if input frequency is a valid SO or CCAT-prime central band frequency
    index = freq == nu
    if np.sum(index) == 0:

        print('Warning: Input frequency is not a valid SO or CCAT-prime band. Try 27, 39, 93, 145, 225, 279, 220, 280, 350, 405, or 860 GHz')
        return(None)

    else:

        #Compute power spectrum
        ell = np.linspace(1,lmax,lmax)
        Cl = N_red[index] * (ell/ell_knee)**alpha_knee 
 
        if no_white is True:
            Cl += N_white[index]

        #Create all-sky map
        noise_map = hp.sphtfunc.synfast(Cl, nside_out, lmax=lmax)/1e6

        #Convert units if necessary
        if unit == 'mjy':
            noise_map = convert_units(freq, noise_map, cmb2mjy=True)
        elif unit == 'cmb':
            None
        elif unit == 'rj':
            noise_map = convert_units(freq, noise_map, cmb2rj=True)
        else:
            print('Waring: Unknown unit! Output will be in K_CMB')

        #Return output
        return(np.float32(noise_map))


def simulate_iras_ps(freq, nside_out = 4096, beam_FWHM = None, unit = 'cmb'):

    '''Computes a map that contains all ~250,000 point sources from the IRAS PS catalog.
    The measured IRAS FIR spectra have been fit with modified blackbodies and are extrapolated
    to the mm/sum-mm regime. This function is not part of the default CCAT-p sky model and 
    serves only for allow to reproduce the forecast results obtained with previous sky models  

    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz. Must be a valid SO or CCAT-prime central
        band frequency, i.e. 27, 39, 93, 145, 225, 279, 220, 280, 350, 405, or 860 GHz.
    nside_out: float, optional
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096   
    beam_FWHM: bool, optional
        If set, the output will be convolved with a gaussian. The FWHM of the Gaussian
        in units of arcmin is given by the provided value. Default: None
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.

    Returns
    -------
    ps_model: float array
        Healpix all-sky map containing the IRAS point sources.
    '''

    if beam_FWHM is None:
        print('Warning: beam is not allowed to be None or 0. beam_FWHM will be set to 1 arcmin')
        beam_FWHM = 1

    #read data
    data = ascii.read(data_path + 'catalogs/IRAS_PSC_fit_results.txt')
    RA = np.array(data['RA'])
    DEC = np.array(data['DEC'])
    A = np.array(data['A'])
    T = np.array(data['T'])
    n = len(RA)

    #compute spectra
    flux = A*freq**1.3 * 2 * h * freq**3 / c**2 / (np.exp(h*freq/k_B/T)-1)
    sigma = beam_FWHM/60*np.pi/180 / (2*np.sqrt(2*np.log(2)))
    amplitude_ps = flux / (2*np.pi*sigma**2) / 1e6

    #compute positions
    coord = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='fk5')
    
    gl = coord.galactic.l.value
    gb = coord.galactic.b.value

    phi = gl * np.pi/180.
    theta = (90-gb) * np.pi/180.

    vector = hp.pixelfunc.ang2vec(theta, phi)

    #build map
    ps_model = np.zeros(hp.pixelfunc.nside2npix(nside_out), dtype=np.float32)
    for i in tqdm(np.arange(n)):
        if amplitude_ps[i] > 0:
            index = hp.query_disc(nside_out, vector[i,:], 5*sigma, inclusive = True)
            vec = hp.pixelfunc.pix2vec(nside_out, index)
            distances = hp.rotator.angdist(vector[i,:], vec)
            ps_model[index] += amplitude_ps[i]*np.exp(-0.5*(distances**2/sigma**2))

    #Convert units if necessary
    if unit == 'cmb':
        ps_model = convert_units(freq, ps_model, mjy2cmb=True)
    elif unit == 'mjy':
        None
    elif unit == 'rj':
        ps_model = convert_units(freq, ps_model, mjy2rj=True)
    else:
        ps_model = convert_units(freq, ps_model, mjy2cmb=True)
        print('Waring: Unknown unit! Output will be in K_CMB')
		   
    return(ps_model)


def simulate_nvss_ps(freq, nside_out = 4096, beam_FWHM = None, unit = 'cmb'):

    '''Computes a map that contains all ~1,800,000 point sources from the NVSS PS catalog 
    (Condon et al. 1998). The measured 1.4 GHz fluxes are extrapolated to the mm/sub-mm 
    regime my assigning each source a random spectral index drawn from N(-0.5, 0.1). This 
    function is not part of the default CCAT-p sky model and serves only for allow to 
    reproduce the forecast results obtained with previous sky models  

    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz. Must be a valid SO or CCAT-prime central
        band frequency, i.e. 27, 39, 93, 145, 225, 279, 220, 280, 350, 405, or 860 GHz.
    nside_out: float, optional
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096    
    beam_FWHM: bool, optional
        If set, the output will be convolved with a gaussian. The FWHM of the Gaussian
        in units of arcmin is given by the provided value. Default: None
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'cmb'.

    Returns
    -------
    ps_model: float array
        Healpix all-sky map containing the NVSS point sources.
    '''

    if beam_FWHM is None:
        print('Warning: beam is not allowed to be None or 0. beam_FWHM will be set to 1 arcmin')
        beam_FWHM = 1

    #read data
    data = ascii.read(data_path + 'catalogs/NVSS_ps_results.txt')
    RA = np.array(data['RA'])
    DEC = np.array(data['DEC'])
    flux = np.array(data['flux@1.4GHz'])
    alpha = np.array(data['alpha'])
    n = len(RA)

    #compute spectra
    flux = flux*(freq/1.4e9)**(alpha)
    sigma = beam_FWHM/60*np.pi/180 / (2*np.sqrt(2*np.log(2)))
    amplitude_ps = flux / (2*np.pi*sigma**2) / 1e6

    #compute positions
    coord = SkyCoord(ra=RA*u.degree, dec=DEC*u.degree, frame='fk5')

    gl = coord.galactic.l.value
    gb = coord.galactic.b.value

    phi = gl * np.pi/180.
    theta = (90-gb) * np.pi/180.

    vector = hp.pixelfunc.ang2vec(theta, phi)

    #build map
    ps_model = np.zeros(hp.pixelfunc.nside2npix(nside), dtype=np.float32)
    for i in tqdm(np.arange(n)):
        if amplitude_ps[i] > 0:
            index = hp.query_disc(nside, vector[i,:], 5*sigma, inclusive = True)
            vec = hp.pixelfunc.pix2vec(nside, index)
            distances = hp.rotator.angdist(vector[i,:], vec)
            ps_model[index] += amplitude_ps[i]*np.exp(-0.5*(distances**2/sigma**2))

    #Convert units if necessary
    if unit == 'cmb':
        ps_model = convert_units(freq, ps_model, mjy2cmb=True)
    elif unit == 'mjy':
        None
    elif unit == 'rj':
        ps_model = convert_units(freq, ps_model, mjy2rj=True)
    else:
        ps_model = convert_units(freq, ps_model, mjy2cmb=True)
        print('Waring: Unknown unit! Output will be in K_CMB')
		   
    return(ps_model)


def ccatp_sky_model(freq, sensitivity = None, components = 'all', red_noise = False, cl_file = None, lensed = True, out_file = None, nside_out = 4096, lmax = None, beam_FWHM = None, template = 'WebSky', unit = 'cmb'):

    '''Computes an all-sky map of the simulated microwave sky at the specified frequency. 
    The CCAT-prime sky model uses all-sky templates in the healpix format and includes the
    most important galactic foregrounds, extragalactic backgrounds, the tSZ and kSZ effects 
    of galaxy clusters, and instrumental white/red noise. The individual components can be
    selectively switched on and of.

    Parameters
    ----------
    freq: float or float array
        Frequency of the output map in Hz. If freq =-1 mean that we get frequency independent maps. 
    sensitivity: float, optional
        sensitivity of the instrument in units of micro K_CMB - arcmin. Default: None
    components: string or list of stings, optional
        List containing the names of the components to be included in the sky model.
        Any combination of the following individual components is possible: 'gal_synchrotron', 
        'gal_dust', 'gal_freefree', 'gal_ame', 'cib', 'radio_ps', 'cmb', 'tsz', 'ksz'. All 
        components are used if components is set to 'all'. Default: 'all'
    red_noise: bool, optional
        If True, a realistic white noise + red noise atmospheric model is added to the
	sky model in case the imput frequency is a valid SO or CCAT-prime central
        band frequency, i.e. 27, 39, 93, 145, 225, 279, 220, 280, 350, 405, or 860 GHz.
        Default: False
    cl_file: string, optional
        name of a file containing a CMB power spectrum computed e.g. with CAMP. The
        first column of the file has to correspond to ell, the second column to 
        Cl ell(ell+1)/2pi. If set, a random realization of the CMB based on the provided
        powerspectrum will be added to the data. Default: None
    lensed: bool, optional
        If True, lensed SO, Sehgal and WebSky CMB maps will be used. Default: True 
    out_file: string, optional
        If set, the results will be written as a healpy .fits file of the given name. 
        Default: None	
    nside_out: float, optional
        Healpix nside parameter of the output map. Must be a valid value for nside.
        Default: 4096
    lmax: float, optional
        Maximum value of the multipolemoment. Default: 3*nside_out-1            
    beam_FWHM: bool, optional
        If set, the output will be convolved with a gaussian. The FWHM of the Gaussian
        in units of arcmin is given by the provided value. Default: None
    template: bool, optional
        Determines the all-sky foregrounds templates to be used to build the sky model.
        If 'Sehgal' is chosen, simulations by Sehgal et al. (2010) are used.
        If 'SO' is chosen, the Simons Observatory sky model provided by Colin Hill and 
        based on the simulations by Sehgal et al. (2010) is used. If 'WebSky' is chosen,
        the used templates will be based on the WebSky Extragalactic CMB Mocks provided 
        by CITA. If 'SO_reproduced' is chosen, the SO sky model is reproduced directly 
        from the Sehgal et al. (2010) data. Default: 'WebSky'
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

    if components == 'all':
        components = ['gal_synchrotron', 'gal_dust', 'gal_freefree', 'gal_ame', 'cib', 'radio_ps', 'cmb', 'tsz', 'ksz']

    #build sky model from the individual components
    allsky_map = np.zeros(hp.nside2npix(nside_out))

    if ('gal_synchrotron' in components) or ('gal_dust' in components) or ('gal_freefree' in components) or ('gal_ame' in components):
        print('Computing Galactic foregounds...')
        allsky_map += simulate_gal_foregrounds(freq, components, nside_out = nside_out, lmax = lmax, beam_FWHM = None, unit = unit)
        print('Galactic foregounds complete.')

    if 'cib' in components:
        print('Computing CIB...')
        allsky_map += simulate_cib(freq, nside_out = nside_out, lmax = lmax, beam_FWHM = None, template = template, unit = unit)
        print('CIB complete.')

    if 'radio_ps' in components:
        print('Computing radio PS...')
        allsky_map += simulate_radio_ps(freq, nside_out = nside_out, lmax = lmax, beam_FWHM = None, template = template, unit = unit)
        print('Radio PS complete.')

    if 'cmb' in components:
        print('Computing CMB...')
        allsky_map += simulate_cmb(freq, cl_file = cl_file, lensed = lensed, nside_out = nside_out, lmax = lmax, beam_FWHM = None, template = template, unit = unit)
        print('CMB complete.')

    if 'tsz' in components:
        print('Computing tSZ...')
        allsky_map += simulate_tSZ(freq, nside_out = nside_out, lmax = lmax, beam_FWHM = None, template = template, unit = unit)
        print('tSZ complete.')

    if 'ksz' in components:
        print('Computing kSZ...')
        allsky_map += simulate_kSZ(freq, nside_out = nside_out, lmax = lmax, beam_FWHM = None, template = template, unit = unit)
        print('kSZ complete.')

    #Smooth map if necessary
    if beam_FWHM is not None:
        print('Begin smoothing...')
        allsky_map = hp.sphtfunc.smoothing(allsky_map, iter = 1, lmax = lmax, fwhm = beam_FWHM/60*np.pi/180)
        print('Smoothing complete.')

    #Add noise if necessary
    if red_noise is True:
        if sensitivity is True:
            print('Warning! You request to apply both red noise + white noise and white noise. The White noise sensitivity parameter will be overwritten to None. If this is not what you want than check the settings and re-run.')
            sensitivity = None
        print('Computing red noise...')        
        allsky_map += simulate_atmosphere(freq, nside_out = nside_out, lmax = lmax, beam_FWHM = None, unit = unit)
        print('Red noise complete.')

    if sensitivity is not None and red_noise is False:
        print('Computing white noise...')        
        allsky_map += simulate_white_noise(freq, sensitivity, nside_out = nside_out, unit = unit)
        print('White noise complete.')

    if out_file is not None:
        print('The result will be written to ' + out_file)
        hp.fitsfunc.write_map(out_file, allsky_map, overwrite = True, dtype = np.float32)
        allsky_map = None

    print('Done!')

    return(allsky_map)
