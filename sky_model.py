import numpy as np 
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from astropy.io import fits
from astropy import constants as cst

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
	
	h=cst.h.value 
	k_B = cst.k_B.value
	c = cst.c.value
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
Markdown
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
        
    if types == 'CITA':
        
        if lensed == True: 
        
            CMB = alm2map_CITA(data_path='/vol/arc3/data1/sz/CITA/',file_name='lensed_alm.fits', nside=4096,
                               lmax=4096*3-1)

        else: 
            
            CMB = alm2map_CITA(data_path='/vol/arc3/data1/sz/CITA/',file_name='unlensed_alm.fits', nside=4096,
                               lmax=4096*3-1)      
            
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
		CMB=CMB
            
    if types == 'SO': 
        
        data_path = '/vol/arc3/data1/sz/SO_sky_model/CMB_SZ_maps/'
        file_name = 'Sehgalsimparams_healpix_4096_KappaeffLSStoCMBfullsky_phi_SimLens_Tsynfastnopell_fast_lmax8000_nside4096_interp2.5_method1_1_lensed_map.fits'
        CMB = hp.read_map(data_path + file_name)
        
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
		CMB=CMB
    
    if types == 'Sehgal':
        
        data_path = '/vol/arc3/data1/sz/Sehgal/'
        data_save = '/vol/arc3/data1/sz/CCATp_sky_model/workspace_maude/Sehgal/'
        
        if lensed == True:
            
            file_name = '030_lensedcmb_healpix.fits'
            CMB = hp.read_map(data_path + file_name)
            
        else: 
            
            file_name = '030_unlensedcmb_healpix.fits'
            CMB = hp.read_map(data_path + file_name)
                
        if unit_out == 'Jysr': 
                
            CMB = CMB 
                
        if unit_out == 'MJysr': 
                
            CMB = CMB * 10**-6
                
        if unit_out == 'K':
            
            CMB = CMB * 10**-6
                
            CMB = convert_units(freq=freq, values=CMB, cmb2mjy=False, mjy2cmb=True, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)
                
        if unit_out == 'mK':
                
            CMB = CMB * 10**-6
                
            CMB = convert_units(freq=freq, values=CMB, cmb2mjy=False, mjy2cmb=True, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)  
                
            CMB = CMB * 10**6
            
        if unit_out == 'RJ': 
        
            CMB = CMB*10**-6
            
            CMB = convert_units(freq=freq, values=CMB, cmb2mjy=False, mjy2cmb=False, rj2mjy=False, mjy2rj=True, 
                              cmb2rj=True, rj2cmb=False)              
            
            
        if nside_out < nside or nside_out > nside: 
            
            CMB = udgrade_NSIDE(maps=CMB, nside=nside_out)
            
        else: 
            
            CMB = CMB
        
    return CMB

def simulate_kSZ(simu,freq,unit_out,nside,nside_out):
    
    if simu == 'SO':
        
        #Fixed datas : 
        data_path='/vol/arc3/data1/sz/SO_sky_model/CMB_SZ_maps/'
        data_save = '/vol/arc3/data1/sz/CCATp_sky_model/workspace_maude/SO/'
        pictures_path = '/vol/arc3/data1/sz/SO_sky_model/pictures/'
        file_in = data_path + '148_ksz_healpix_nopell_Nside4096_DeltaT_uK.fits'
        
        ykSZ_map = hp.read_map(file_in,nest=False,partial=False)
    
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
            
    if simu == 'Sehgal' : 
        
        #Fixed datas :         
        data_path='/vol/arc3/data1/sz/Sehgal/'
        data_save = '/vol/arc3/data1/sz/CCATp_sky_model/workspace_maude/Sehgal/'
        pictures_path = '/vol/arc3/data1/sz/Sehgal/pictures/'
        file_in = '030_ksz_healpix.fits'
        
    
        ykSZ_map = hp.read_map(data_path + file_in)
                
        if unit_out == 'Jysr': 
                
            ykSZ_map = ykSZ_map
                
        if unit_out == 'MJysr': 
                
            ykSZ_map = ykSZ_map * 10**-6
                
        if unit_out == 'K':
            
            ykSZ_map = ykSZ_map * 10**-6
                
            ykSZ_map = convert_units(freq=freq, values=ykSZ_map, cmb2mjy=False, mjy2cmb=True, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)
                
        if unit_out == 'mK':
                
            ykSZ_map = ykSZ_map * 10**-6
                
            ykSZ_map = convert_units(freq=freq, values=ykSZ_map, cmb2mjy=False, mjy2cmb=True, rj2mjy=False, mjy2rj=False, 
                          cmb2rj=False, rj2cmb=False)  
                
            ykSZ_map = ykSZ_map * 10**6
            
        if unit_out == 'RJ': 
        
            ykSZ_map = ykSZ_map*10**-6
            
            ykSZ_map = convert_units(freq=freq, values=ykSZ_map, cmb2mjy=False, mjy2cmb=False, rj2mjy=False, mjy2rj=True, 
                              cmb2rj=True, rj2cmb=False)              
            
            
        if nside_out < nside or nside_out > nside: 
            
            ykSZ_map = udgrade_NSIDE(maps=ykSZ_map, nside=nside_out)
            
        else: 
            
            ykSZ_map = ykSZ_map
            
        hp.mollview(map=ykSZ_map, coord=None, nest=False,title='',norm='hist', xsize=2000,return_projected_map=True)
        
    
    if simu == 'CITA':
        
        #Fixed datas : 
        data_path = '/vol/arc3/data1/sz/CITA/'
        file_in = 'ksz.fits'
        
        ykSZ_map = hp.read_map(data_path + file_in)
        
        if unit_out == 'K':
            
            ykSZ_map = ykSZ_map*10**-6
            
        if unit_out =='mK':
                
            ykSZ_map=ykSZ_map

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
            
            ykSZ_map = ykSZ_map
            
        hp.mollview(map=ykSZ_map, coord=None, nest=False,title='',norm='hist', xsize=2000,return_projected_map=True)
    
    
    return ykSZ_map
