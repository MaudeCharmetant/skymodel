import numpy as np 
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
from astropy.io import fits

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


def CMB_make_CAMB(data_path,file_name,nside,name_writemap,title_map,maps_unit,pictures_path,name_imfile,lmax,
                  name_imPS): 
    
        
    """
    Function which compute a CMB map from a power spectrum.

    Parameters
    ----------
    
    data_path : str
        Path were the data of the maps are stored and we the cutout are going to be stored. 
    file_name : str
        Name of the .dat contaning the values of the power spectrum given by CAMB.
    nside : int
        Resolution of the map. The power spectrum will go until l = 3 * nside -1.
    name_writemap : str 
        Name of the FITS file under which the map is going to be saved. 
    title_map : str 
        Title displayed on the image of the Healpix map.
    maps_unit : str 
        Unit displayed on the image of the Healpix map.         
    pictures_path : str 
        Path where we are going to save the pictures. 
    name_imfile : str 
        Name under which the image of the map is going to be saved. 
    lmax : int 
        Maximum Multipole desired or that can be reach by the data, in which case lmax = 3*nside-1
    name_imPS : str 
        Name given to the image containing the power spectrum of your CMB map. 
        
    Returns
    -------
    str
        Tell us where the function stored the images. 

    """
    
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
    hp.write_map(data_path + name_writemap,CMB,overwrite=True) # Wrtie the map as a FITS file
    hp.mollview(map=CMB,coord=None, nest=False, title=title_map, unit=maps_unit, norm='hist')
    plt.savefig(pictures_path + name_imfile  + '.png') #Save the image of the filtered map 
    plt.show()
    
    #Display the Power spectrum
    fig, ax = plt.subplots()
    ax.plot(ell,TT_final[ell]*ellfactor,'-') #Display the power spectrum multiply by the ell factor
    ax.set_title('$CMB$ $power$ $spectrum$',fontsize=28)
    ax.set_xlabel('$l$',fontsize=28)
    ax.set_ylabel('$l(l+1)/2\pi C_{l}^{TT}$',fontsize=28) 
    plt.xticks(fontsize=28)
    plt.yticks(fontsize=28)
    plt.rcParams['figure.figsize'] = [15, 10]
    plt.savefig(pictures_path + name_imPS + '.png') #Save the image of the power spectrum
    plt.show()
    
    #Feedback operator : 
    print('Images of the maps and Power spectrum saved at : '+pictures_path+' FITS files saves in : '+data_path)
    
    
    return CMB
