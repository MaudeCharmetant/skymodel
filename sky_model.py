import numpy as np 
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as cst
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.7255)

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

  
def random_map(PS,nside,lmax):  
    
    """
    Function generating random realization of Healpix map from a given power spectrum. 

    Parameters
    ----------
     
    PS : array
        Array containing the power spectrum that will be used to generate the map.
    nside : int 
        Number of seperation of a healpix pixel. 
    lmax : int 
        maximum l that can be reach, by default lmax=3*nside-1.
        
    Returns
    -------
    array
        Array containing the map generated from the power spectrum.

    """
        
    #Generate display and save the map generated from the Power spectrum : 
    random_map = hp.sphtfunc.synfast(PS, nside, lmax=lmax, mmax=lmax, alm=False, pol=False, pixwin=False, fwhm=0.0, sigma=None, new=False, verbose=True)
    
    #Feedback user : 
    print('The random map have been generated')

    
    return random_map 

def udgrade_NSIDE(maps,nside):
    
    """
    Function which upgrade or degrade the NSIDE of an HEALPix map. 

    Parameters
    ----------
    
    maps : array
        Array containing the map we want to ud_grade.
    nside : int 
        Value of the new resolution we want to apply to the healpix map to upgrade or degrade it.
        
    Returns
    -------
    array
        Array contaning the upgraded map.

    """
    
    #Upgarde or degrade the map : 
    ud_map = hp.pixelfunc.ud_grade(map_in=maps, nside_out=nside, pess=False, order_in='RING', order_out=None, power=None, dtype=None)
    
    #Feedback to the operator : 
    print('the map have been ud_grade to '+str(nside))
    
    return ud_map

def alm2map_CITA(data_path,file_name, nside,lmax):
    
    """
    Function which create a map out of a file containing the alm under a tuple format. 

    Parameters
    ----------
    
    data_path : str
        Path were the data of the maps are stored and we the cutout are going to be stored. 
    file_name : str
        Name of the file containing the alms. 
    name_final_maps : str
        Name under which the maps image and fits file are going to be saved.    
    nside : int 
        Number of cut in the the healpix originals pixels. It is link to the final number of pixels in the map, N_pix = 
        12 x nside^2 
    title_map : str 
        Title that will be displayed on the image of the map. 
    unit_map : str 
        Units of the maps that are displayed, for exemple K_CMB, MJy/sr. This is going to be use for figure titles and 
        names of files.         
    pictures_path : str
        Path where we are going to save the pictures.   
    lmax : int 
        Maximum l that can be reach. By default lmax = 3*nside-1.
        
    Returns
    -------
    array
        Containing the map at a given frequency. 

    """
    
    #Open and read the file : 
    file = data_path + file_name
    hdul = fits.open(file)
    data = hdul[1].data

    #Get each parts of the alms : 
    j = complex(0,1) #Create the pure imaginary i 
    alm_1 = data['real'][:] #Take the full column of data under the name 'real'
    alm_2 = data['imag'][:]
    alm_T2 = np.array(alm_2,dtype=complex) #Make alm_2 a complex array 
    alm_T = alm_1 + alm_T2*j 
    
    #Define the l : 
    ell = np.arange(lmax) #Array going from 1 to lmax
    ellfactor = ell*(ell+1)/(2.*np.pi) #Array containing all the values of the factor used in CMB science  

    #Display map reconstruct from the alms : 
    map_T = hp.alm2map(alm_T, nside, pol=False, inplace=False) #Make a map out of the alm
      
    #Feeedback operator : 
    print('Map have been computed from the alm') 
    
    return map_T

def Simulate_CMB(data_path,file_in,nside,maps_unit,lmax,types,unit_out,freq,lensed,nside_out): 
    
        
    """
    Function which compute a CMB map from a power spectrum.

    Parameters
    ----------
    
    data_path : str
        Path were the data of the maps are stored and we the cutout are going to be stored. 
    file_in : str or array 
        Name of the .dat contaning the values of the power spectrum given by CAMB.
        Or array containing the power spectrum to generate random maps. 
    nside : int
        Resolution of the map. The power spectrum will go until l = 3 * nside -1.
    maps_unit : str 
        Unit in which the original map or power spectrum is coming in.       
    lmax : int 
        Maximum Multipole desired or that can be reach by the data, in which case lmax = 3*nside-1
    types : str 
    	type of CMB map you want to produce. 1. random 2. from CAMB 3. from CITA 4. from SO 5. from Sehgal.
    unit_out : str 
    	Unit in which you want to get your CMB map in. Can be : 'K', 'mK', 'Jysr', 'MJysr', 'RJ'
    freq : float 
    	frequency in which you want to produce a CMB map. Has to be given on Hz. 
    lensed : bool 
    	if True select the lensed CMB map when possible. This is only possible for 'CITA' and 'Sehgal.'
    nside_out : int 
    	If you wish to change the nside of the CMB map, this is avalaible for 'CITA', 'Sehgal' and 'SO'.
        
    Returns
    -------
    array
        Contaning the CMB map. 

    """
    
    if types == 'random':
        
        #Compute and convert the random map : 
        if maps_unit == 'mK' :            
        
            CMB = random_map(PS=file_in,nside=nside,lmax=lmax)
            
            if unit_out == 'K':
            
                CMB = CMB*10**-6
            
            if unit_out =='mK':
                
                CMB=CMB

            if unit_out == 'MJysr':  
                
                CMB = CMB*10**-6
            
                CMB = convert_units(freq=freq, values=CMB, cmb2mjy=True, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)
                
            if unit_out == 'Jysr': 
                
                CMB = CMB*10**-6
            
                CMB = convert_units(freq=freq, values=CMB, cmb2mjy=True, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)
                
                CMB = CMB*10**6
            
            if unit_out == 'RJ': 
        
                CMB = CMB*10**-6
            
                CMB = convert_units(freq=freq, values=CMB, cmb2mjy=False, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=True, rj2cmb=False)   
        
        #Compute and convert the random map : 
        if maps_unit =='K': 
            
            if unit_out == 'K':
            
                CMB = CMB
            
            if unit_out =='mK':
                
                CMB=CMB*10**6

            if unit_out == 'MJysr':  
            
                CMB = convert_units(freq=freq, values=CMB, cmb2mjy=True, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)
                
            if unit_out == 'Jysr': 
            
                CMB = convert_units(freq=freq, values=CMB, cmb2mjy=True, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)
                
                CMB = CMB*10**6
            
            if unit_out == 'RJ': 
            
                CMB = convert_units(freq=freq, values=CMB, cmb2mjy=False, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=True, rj2cmb=False)
                
    if types == 'CAMB':
    
        #Load the datas : 
        data = np.loadtxt(data_path + file_name)
        TT_data = data[:,1] # Take only the first column, which is the temperature T

        #Vairables : 
        ell = np.arange(2,lmax) #Array going from 2 to lmax
        ellfactor = ell*(ell+1)/(2.*np.pi) #Array containing all the values of the factor used in CMB science 

        TT = TT_data[ell] / ellfactor #The file given by CAMB is the power spectrum multiplied by this factor 
        #Because the monopole and dipole are not in the data of the power spectrum given by CAMB we need to add them back
        #They need to be 0 because we usually remove them no to influence our studies. 
        TT_1 = np.insert(TT,0,0)
        TT_final = np.insert(TT_1,0,0)

        fb.file2FITS(data=TT_final,dtype=np.float32,data_path=data_path,name_save='TT_CMB_CAMB_'+str(nside), overwrite=True)

        #Compute the CMB map from the power spectrum :
        CMB = hp.sphtfunc.synfast(TT_final, nside, lmax=lmax, mmax=lmax, alm=False, pol=False, pixwin=False)
        
        #Convert the map from micro to ... : 
        if unit_out == 'K':
            
            CMB = CMB*10**-6
            
        if unit_out =='mK':
                
            CMB=CMB

        if unit_out == 'MJysr':  
                
            CMB = CMB*10**-6
            
            CMB = convert_units(freq=freq, values=CMB, cmb2mjy=True, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)
                
        if unit_out == 'Jysr': 
                
            CMB = CMB*10**-6
            
            CMB = convert_units(freq=freq, values=CMB, cmb2mjy=True, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)
                
            CMB = CMB*10**6
            
        if unit_out == 'RJ': 
        
            CMB = CMB*10**-6
            
            CMB = convert_units(freq=freq, values=CMB, cmb2mjy=False, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=True, rj2cmb=False) 
        
    if types == 'CITA' or types == 'Sehgal' or types =='SO' : 
        
        if types == 'CITA':
        
            if lensed == True: 
        
                CMB = alm2map_CITA(data_path='/vol/arc3/data1/sz/CITA/',file_name='lensed_alm.fits', nside=4096,
                               lmax=4096*3-1)

            else: 
            
                CMB = alm2map_CITA(data_path='/vol/arc3/data1/sz/CITA/',file_name='unlensed_alm.fits', nside=4096,
                               lmax=4096*3-1)      
    
        if types == 'SO': 

            data_path = '/vol/arc3/data1/sz/SO_sky_model/CMB_SZ_maps/'
            file_name = 'Sehgalsimparams_healpix_4096_KappaeffLSStoCMBfullsky_phi_SimLens_Tsynfastnopell_fast_lmax8000_nside4096_interp2.5_method1_1_lensed_map.fits'
            CMB = hp.read_map(data_path + file_name)


        if types == 'Sehgal':

            data_path = '/vol/arc3/data1/sz/Sehgal/'
            data_save = '/vol/arc3/data1/sz/CCATp_sky_model/workspace_maude/Sehgal/'

            if lensed == True:

                file_name = '030_lensedcmb_healpix.fits'
                CMB = hp.read_map(data_path + file_name)

            else: 

                file_name = '030_unlensedcmb_healpix.fits'
                CMB = hp.read_map(data_path + file_name)


            #Create y_kSZ map : 
            conv = conv_vector(30)
            CMB = (CMB / conv)*T_CMB*10**6


        if unit_out == 'K':

            CMB = CMB*10**-6

        if unit_out =='mK':

            CMB=CMB

        if unit_out == 'MJysr':  

            CMB = CMB*10**-6

            CMB = convert_units(freq=freq, values=CMB, cmb2mjy=True, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                              cmb2rj=False, rj2cmb=False)

        if unit_out == 'Jysr': 

            CMB = CMB*10**-6

            CMB = convert_units(freq=freq, values=CMB, cmb2mjy=True, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                              cmb2rj=False, rj2cmb=False)

            CMB = CMB*10**6

        if unit_out == 'RJ': 

            CMB = CMB*10**-6

            CMB = convert_units(freq=freq, values=CMB, cmb2mjy=False, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                              cmb2rj=True, rj2cmb=False) 

        if nside_out < nside or nside_out > nside: 

            CMB = udgrade_NSIDE(maps=CMB, nside=nside_out)

        else: 

            CMB = CMB    
   
    return CMB



def D_I_tSZ(x,y):
    
    """
    Function which compute the tSZ spectral shape. 

    Parameters
    ----------
    
    x : array
        Frequency range over which the tSZ spectral shape will be computed. 
    y : float
        Value of the Compton-y parameter assumed here. 
        
    Returns
    -------
    array
        Array contaning the Variarion of intensity produced by tSZ over the fequencies. 

    """
    
    #Compute Delta I : 
    I_0 = (2*(cst.k_B.value*T_CMB)**3)/(cst.h.value*cst.c.value)**2    
    x_nu = np.array((cst.h.value*x)/(cst.k_B.value*T_CMB))    
    Delta_I = np.array(I_0*y*(x_nu**4)*(np.exp(x_nu)/(np.exp(x_nu)-1)**2)*((x_nu*(np.exp(x_nu)+1)/(np.exp(x_nu)-1))-4))
    
    #Give feedback to the operator : 
    print("Delta I as been computed ")
    
    return   Delta_I

def mixing_vector(frequency): 
    
    """
    Function which compute the multiplying vector to transform a y-map into a tSZ(f). 

    Parameters
    ----------
    
    frequency : array
        Array or single number containing the frequencies we want to compute the mixing vector for. 
        
    Returns
    -------
    array
        Array contaning the multiplying vector. 

    """    
    
    #Initilisation : 
    freq = np.arange(0,1000)*10**9  
    mix_vect = []
    Delta_I = D_I_tSZ(freq,1)

    #For each frequency channel, compute Delta_I : 
    mix_vect.append(Delta_I[int(frequency)]*(10**20))
    
        
    #Give feeback to the operator :    
    print('The mixing vector is : ',mix_vect)
   
    return mix_vect

def simulate_tSZ(simu,freq,unit_out,rescale,nside,nside_out):

    
    """
    Function which compute tSZ maps at different frequencies and different nside. 
    
    Parameters
    ----------
    
    simu : str
        Name of the simulation we want a tSZ map of. Possibles choices : 'SO','CITA','Sehgal'
    freq : float 
        frequency in which you want to produce a tSZ map. Has to be given on Hz.
    unit_out : str 
        Unit in which you want to get your tSZ map in. Can be : 'K', 'mK', 'Jysr', 'MJysr', 'RJ'
    rescale : bool 
        if True, in the case of 'SO' divide by the rescaling factor that was applied to the datas.  
    nside : int
        Original nside of the simulations, SO:4096, CITA:4096, Sehgal:8192
    nside_out : int 
        If you wish to change the nside of the tSZ map, this is avalaible for 'CITA', 'Sehgal' and 'SO'.
        
    Returns
    -------
    array
        Contaning the tSZ map at a given frequency. 
    """
    
    if simu == 'SO':
        
        #Fixed datas : 
        data_path='/vol/arc3/data1/sz/SO_sky_model/CMB_SZ_maps/'
        data_save = '/vol/arc3/data1/sz/CCATp_sky_model/workspace_maude/SO/'
        pictures_path = '/vol/arc3/data1/sz/SO_sky_model/pictures/'
        file_in ='tSZ_skymap_healpix_nopell_Nside4096_y_tSZrescale0p75.fits'
        
        tSZ = hp.read_map(data_path + file_in)
        tSZ = tSZ * 2.726e6
        
    if simu == 'CITA': 
        
        #Fixed datas :         
        data_path='/vol/arc3/data1/sz/CITA/'
        data_save = '/vol/arc3/data1/sz/CCATp_sky_model/workspace_maude/CITA/'
        pictures_path = '/vol/arc3/data1/sz/CITA/pictures/'
        file_in = 'tsz.fits'
        
        tSZ = hp.read_map(data_path + file_in)        
        tSZ = tSZ * 2.726e6
        
    if simu == 'Sehgal': 
        
        #Fixed datas :         
        data_path='/vol/arc3/data1/sz/Sehgal/'
        data_save = '/vol/arc3/data1/sz/CCATp_sky_model/workspace_maude/Sehgal/'
        pictures_path = '/vol/arc3/data1/sz/Sehgal/pictures/'
        file_in='030_tsz_healpix.fits'
        
        #Create compton-y map : 
        tSZ_freq = hp.read_map(data_path + file_in)
           
        multiplier = mixing_vector(30)
        tSZ = (tSZ_freq / multiplier) #Compton-y map 
         

    #Get tSZ at different frequencies :           
    multiplier = mixing_vector(freq*10**-9)
    tSZ = (tSZ * multiplier)
        
    if simu =='SO' and rescale == 'True': 
            
        tSZ_map = tSZ / 0.75 #Factor used to bring Sehgal(2010) values of tSZ in agreement with Planck,ACT,SPT
        
    else:
            
        tSZ_map = tSZ
            
    if unit_out == 'mK':
            
        tSZ_map = tSZ_map
            
    if unit_out == 'K':
            
        tSZ_map = tSZ_map*10**-6

    if unit_out == 'MJysr':  
                
        tSZ_map = tSZ_map*10**-6
            
        tSZ_map = convert_units(freq=freq, values=tSZ_map, cmb2mjy=True, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)
                
    if unit_out == 'Jysr': 
                
        tSZ_map = tSZ_map*10**-6
            
        tSZ_map = convert_units(freq=freq, values=tSZ_map, cmb2mjy=True, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)
                
        tSZ_map = tSZ_map*10**6
            
    if unit_out == 'RJ': 
        
        tSZ_map = tSZ_map*10**-6
            
        tSZ_map = convert_units(freq=freq, values=tSZ_map, cmb2mjy=False, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=True, rj2cmb=False) 
            
    if nside_out < nside or nside_out > nside: 
            
        tSZ_map = udgrade_NSIDE(maps=tSZ_map, nside=nside_out)
        
    else: 
            
        tSZ_map = tSZ_map
            
    hp.mollview(map=tSZ_map, coord=None, nest=False,title='',norm='hist', xsize=2000,return_projected_map=True)
        
   
    return tSZ_map


def DT_kSZ(x,y):

    """
    Function which compute the kSZ spectral shape.

    Parameters
    ----------
    
    x : str
        Path were the data of the maps are stored and we the cutout are going to be stored. 
    y : str
        Name of the data file. 
        
    Returns
    -------
    str
        Tell us where the function stored the datas and images.

    """
    
    I_0 = 2*(k_B*T_CMB)**3/(h*c)**2*1e26
    x_nu = (h*x)/(k_B*T_CMB)
    Delta_T = ((np.exp(x_nu)-1)**2)/(I_0*x_nu**4*np.exp(x_nu))
    
    print("Delta T as been computed ")
    
    return   Delta_T

def conv_vector(frequency): 
	
    """
    Function which compute the conversion from kSZ to ykSZ.

    Parameters
    ----------
    
    frequency : float
        Frequency of the original kSZ map.
        
    Returns
    -------
    array
        Contanining the conversion factor.

    """
    
    freq = np.arange(0,1000)*10**9  
    Delta_T = DT_kSZ(freq,1)
    conv_vect = []
    
    conv_vect.append(1/Delta_T[int(frequency)])
    print('The conversion vector is : ',conv_vect)

    return conv_vect

def simulate_kSZ(simu,freq,unit_out,nside,nside_out):

    """
    Function which compute kSZ maps at different frequencies and different nside. 
    
    Parameters
    ----------
    
    simu : str
        Name of the simulation we want a kSZ map of. Possibles choices : 'SO','CITA','Sehgal'
    freq : float 
        frequency in which you want to produce a kSZ map. Has to be given on Hz.
    unit_out : str 
        Unit in which you want to get your kSZ map in. Can be : 'K', 'mK', 'Jysr', 'MJysr', 'RJ'
    nside : int
        Original nside of the simulations, SO:4096, CITA:4096, Sehgal:8192
    nside_out : int 
        If you wish to change the nside of the kSZ map, this is avalaible for 'CITA', 'Sehgal' and 'SO'.
        
    Returns
    -------
    array
        Contaning the kSZ map at a given frequency. 
    """
        
    if simu == 'SO':
        
        #Fixed datas : 
        data_path='/vol/arc3/data1/sz/SO_sky_model/CMB_SZ_maps/'
        data_save = '/vol/arc3/data1/sz/CCATp_sky_model/workspace_maude/SO/'
        pictures_path = '/vol/arc3/data1/sz/SO_sky_model/pictures/'
        file_in = '148_ksz_healpix_nopell_Nside4096_DeltaT_uK.fits'
        
        ykSZ_map = hp.read_map(data_path +file_in)
        
    if simu == 'CITA':
        
        #Fixed datas : 
        data_path = '/vol/arc3/data1/sz/CITA/'
        file_in = 'ksz.fits'
        
        ykSZ_map = hp.read_map(data_path +file_in)
        
    if simu == 'Sehgal' : 
        
        #Fixed datas :         
        data_path='/vol/arc3/data1/sz/Sehgal/'
        data_save = '/vol/arc3/data1/sz/CCATp_sky_model/workspace_maude/Sehgal/'
        pictures_path = '/vol/arc3/data1/sz/Sehgal/pictures/'
        file_in = '030_ksz_healpix.fits'    
        
        kSZ_map = hp.read_map(data_path + file_in)
        
        #Create y_kSZ map : 
        conv = conv_vector(30)
        ykSZ_map = (kSZ_map / conv)*T_CMB*10**6
        
    
    if unit_out == 'mK': 
        
        ykSZ_map = ykSZ_map 
    
    if unit_out == 'K':
            
        ykSZ_map = ykSZ_map*10**-6


    if unit_out == 'MJysr':  
                
        ykSZ_map = ykSZ_map*10**-6
            
        ykSZ_map = convert_units(freq=freq, values=ykSZ_map, cmb2mjy=True, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)
                
    if unit_out == 'Jysr': 
                
        ykSZ_map = ykSZ_map*10**-6
            
        ykSZ_map = convert_units(freq=freq, values=ykSZ_map, cmb2mjy=True, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)
                
        ykSZ_map = ykSZ_map*10**6
            
    if unit_out == 'RJ': 
        
        ykSZ_map = ykSZ_map*10**-6
            
        ykSZ_map = convert_units(freq=freq, values=ykSZ_map, cmb2mjy=False, mjy2cmb=False, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=True, rj2cmb=False) 
            
    if nside_out < nside or nside_out > nside: 
            
        ykSZ_map = udgrade_NSIDE(maps=ykSZ_map, nside=nside_out)
            
    else:
            
        ykSZ_map =ykSZ_map
            
    hp.mollview(map=ykSZ_map, coord=None, nest=False,title='',norm='hist', xsize=2000,return_projected_map=True)
                 
                 
    return ykSZ_map

def generate_atmosphere(freq, nside_out=4096, lmax = None, beam_FWHM = None, unit = "mjy"):
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
        'rj' --> K_RJ (brightness temperature). Default: 'mjy'.

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
            noise_map = sz.convert_units(freq, noise_map, cmb2mjy=True)
        elif unit == "cmb":
            None
        elif unit == "rj":
            noise_map = sz.convert_units(freq, noise_map, cmb2rj=True)
        else:
            print("Waring: Unknown unit! Output will be in MJy/sr")

        #Return output
        return(np.float32(noise_map))

def generate_atmosphere(freq, components, nside_out=4096, lmax = None, beam_FWHM = None, unit = "mjy"):
    '''Computes an all-sky galactic foregrounds noise map at a given frequency and nside using 
    the Python Sky model (PySM, Thorne et al. 2017), which is build from Planck data. 

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
    unit: bool, optional
        Determines the units of the output map. The available units are 'mjy' --> MJy/sr
        (specific intensity), 'cmb' --> K_CMB (thermodynamic temperature), and 
        'rj' --> K_RJ (brightness temperature). Default: 'mjy'.

    Returns
    -------
    foregrounds: float array
        Healpix all-sky map of the galactic foregrounds at the specified frequency.
    '''

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
    rj2mjy_factor = sz.convert_units(freq, 1e-6, rj2mjy = True)
    
    if np.sum("dust" == components) == 1:
        foregrounds += sky.dust(freq/1e9)[0,:] * rj2mjy_factor
        
    if np.sum("synchrotron" == components) == 1:
        foregrounds += sky.synchrotron(freq/1e9)[0,:] * rj2mjy_factor
        
    if np.sum("freefree" == components) == 1:
        foregrounds += sky.freefree(freq/1e9)[0,:] * rj2mjy_factor
        
    if np.sum("ame" == components) == 1:
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
        foregrounds = sz.convert_units(freq, foregrounds, mjy2cmb=True)
    elif unit == "rj":
        foregrounds = sz.convert_units(freq, foregrounds, mjy2rj=True)
    else:
        print("Waring: Unknown unit! Output will be in MJy/sr")
        
    return(foregrounds)
