from numpy import *
from astropy.constants import G, M_sun, au
from astropy import units as u
from astropy import constants as const
from astropy.convolution import convolve, Gaussian2DKernel

def width_convolution_kernel(v, v_source):
    '''
    This function returns the width of the Gaussian convolution kernel
    that that we can take into account the envelope of the source.
    '''
    a = -2#slope in pixels/(km s^-1) 
    b = 10#Position at v = 0 in pixels.
    return a*abs(v-v_source) + b
    
def ellipse_mask(im, xc, yc, maj_ax, min_ax):
    #Function that takes the 2d array im and returns circular mask for this 
    #image. The circle is centered at (xc,yc) and has a major axis of maj_ax and 
    #a minor axis of min_ax. 
    #All in pixels.
    #Dimensions of the input array.        
    (nx,ny) = im.shape
    #Creating an coordinate grid.
    y,x = mgrid[0:nx,0:ny]
    
    mask = ((x - xc)**2/(maj_ax)**2 + (y - yc)**2/(min_ax)**2 < 1)

    return mask
    
def Gaussian_Smoothing_Kernel(im, xc, yc, sigma_x, sigma_y, rotation):
    #This function returns a Gaussian smoothing kernel the size of the 2D input
    #array. The center of the Gaussian is at (xc, yc), the sigma's in the x and
    #y direction are given by sigma_x and sigma_y (all in pixels), the rotation 
    #of the Gaussian is given by rotation in degrees.
      
    #Dimensions of the input array.        
    (nx,ny) = im.shape
    #Creating an coordinate grid.
    y,x = mgrid[0:nx,0:ny]
    #Change rotation angle to radians.
    rot = radians(rotation)

    #Lets first rotate the coordinates.
    xr=x*cos(rot)-y*sin(rot) 
    yr=x*sin(rot)+y*cos(rot)
    #Rotate the center
    xcr = xc*cos(rot)-yc*sin(rot) 
    ycr = xc*sin(rot)+yc*cos(rot)
    
    #Creating the 2D Gaussian.
    temp = exp(-(((xr-xcr)/sigma_x)**2 +((yr-ycr)/sigma_y)**2)/2)
    #Return the normalized Gaussian.
    return temp/sum(temp)

def Convolve_Datacube(Datacube, Sigma_x, Sigma_y, Rotation, Incl, Velocities, Vlsr, Process_Number):
    '''
    This function convolves a given datacube with an Gaussian Smoothing 
    kernel. 
    Datacube = the 3D array which contains the data that needs to be 
    convolved.
    Sigma_x = the width of the Gaussian convolution kernel in x 
    direction.
    Sigma_y = same as x but now in the y direction.
    Rotation = the rotation of the convolution kernel.
    Incl = inclination of the source.
    Velocities = an array containing the velocities of the datacube's 
    channels. 
    Vlsr = the velocity of the source. 
    Process_Number = the number of this process. Used for reconstructing
    the datacube later.
    '''
    #An elliptical smoothing kernel.
    #Smoothing_Kernel_temp = ellipse_mask(zeros((2*round(B_maj_pix_FWHM) + 5, 2*round(B_maj_pix_FWHM) + 5)), (2*round(B_maj_pix_FWHM) + 5)/2 - 1, (2*round(B_min_pix_FWHM) + 5)/2 - 1, B_maj_pix_FWHM, B_min_pix_FWHM)
    #Normalizing the smoothing kernel. 
    #Smoothing_Kernel = Smoothing_Kernel_temp/sum(Smoothing_Kernel_temp)
    
    #Creating a Gaussian smoothing kernel. 
    smoothing_kernel = Gaussian_Smoothing_Kernel(zeros((51,51)), 25, 25, Sigma_x, Sigma_y, 90 + Rotation)
    
    #Creating the array where the results will be saved. 
    datacube_convolved = zeros_like(Datacube)
    
    #The shape of the datacube.
    nz, ny, nx = Datacube.shape
    
    #Now we convolve over all the channels. 
    i = 0
    while i < nz:
        #If the source is very inclined, almost edge on, we need an 
        #extra convolution to make the thindisk a little bit thicker. 
        if Incl > 84*pi/180.:
            Width_Kernel = width_convolution_kernel(Velocities[i], Vlsr)
            Smoothing_Kernel_Disk = Gaussian_Smoothing_Kernel(zeros((51,51)), 25, 25, Width_Kernel, Width_Kernel, 0)
            Datacube[i,:,:] = convolve(Datacube[i,:,:], Smoothing_Kernel_Disk)
        
        datacube_convolved[i,:,:] = convolve(Datacube[i,:,:], smoothing_kernel)
        i+= 1
        
    return Process_Number, datacube_convolved
 
def Thindisk(mstar, rc, size, int0, fwhm, Boxsize, channel1, channel2, channel3, channel4, incl, material_v):
    '''
    This function creates an thindisk datacube based on given variables and some
    predetermined constants.
    mstar = mass of the central object in solar masses.
    rc = the centrifugal radius in AU.
    size = the size of the total object in AU.
    int0 = the highest intensity of the source in JY/BEAM.
    fwhm = the fwhm of the intensity distribution in JY/BEAM.
    Boxsize = the width of the box that zooms in on the center in pixels.
    channel1, channel2 = between (and including) these channels the data is convolved
    same for channel3 and channel4. 
    '''
    #Definition of some constants.
    
    #Some source info.
    dist = 140#pc
    vlsr = 5.95#km s^-1
    
    pa = 273#Degrees
    
    #The center of the continuum. 
    xcen = 261.696696133
    ycen = 246.99512412

    #Info about the line.   
    linewidth = 0.1

    #Info about the beam.
    B_maj = 1.651514404350E-04#deg
    B_min = 1.320303728183E-04#deg
    BPA = 7.119091033936E+00#deg

    #Info about the datacube.
    NAXIS   = 4                                                
    NAXIS1  = 512                                             
    NAXIS2  = 512                                               
    NAXIS3  = 100                                          
    NAXIS4  = 1
                                          
    CRPIX1  = 2.570000000000E+02#deg                                               
    CDELT1  = -1.388888888889E-05 #deg                                              
    CRVAL1  = 6.997458333333E+01#deg  
                                            
    CRPIX2  = 2.570000000000E+02#deg                                                
    CDELT2  = 1.388888888889E-05#deg                                              
    CRVAL2  = 2.605277777778E+01#deg 
    
    begin_freq = 2.195603600000E+11#Hz
    delta_freq = -1.171799244995E+05#Hz
    begin_pos = 1.000000000000E+00 - 1#pix
    rest_freq = 2.195603600000E+11#Hz

    #The speed of light is.
    c = const.c.value/1000#km s^-1

    #Calculating the begin velocity.
    begin_v = c * (square(rest_freq) - square(begin_freq - delta_freq*begin_pos)) / ( square(rest_freq) + square(begin_freq - delta_freq*begin_pos))
    #Now we calculate the delta v
    begin_v_plus_one = c * (square(rest_freq) - square(begin_freq - delta_freq*(begin_pos + 1))) / ( square(rest_freq) + square(begin_freq - delta_freq*(begin_pos + 1)))
    delta_v = begin_v - begin_v_plus_one
    
    CRPIX3  = 1.00000000000E+00#km s^-1                                              
    CDELT3  = delta_v#km s^-1   
    CRVAL3  = begin_v#km s^-1     
                                            
    CRPIX4  = 1.00000000000E+00                                           
    CDELT4  = 1.00000000000E+00                                              
    CRVAL4  = 1.00000000000E+00 
                                           
    LWIDTH  = delta_v#km s^-1                                             
    LSTEP   = delta_v#km s^-1                                         
    LSTART  = begin_v#km s^-1 

    #---------------------------------------------------------------------------
    #The variables that will be used in the code get their values. 
    v_begin = LSTART
    chanwidth = LSTEP
    nchan = NAXIS3
    npix = NAXIS1

    #Compute the initial grid 
    #New version.

    ra_offset = (arange(npix) - xcen)*abs(CDELT1)*u.deg.to(u.arcsec)
    dec_offset = (arange(npix) - ycen)*abs(CDELT2)*u.deg.to(u.arcsec)

    ra_grid, dec_grid = meshgrid(ra_offset, dec_offset)
    veloc = arange(nchan)*chanwidth + v_begin#Added 29-9-2015

    # Convert the projected coordinates to cylindrical coordinates in the plane of the disk

    pa *= pi / 180. # degrees -> rad
    x_grid = (ra_grid * cos(-pa) - dec_grid * sin(-pa)) # disk major axis
    y_grid = (ra_grid * sin(-pa) + dec_grid * cos(-pa)) # disk minor axis
    theta = zeros((npix, npix)) # theta = 0 along the l.o.s.
    r = zeros((npix, npix))
    mask = y_grid != 0
    incl *= pi / 180. # degrees -> rad
    theta[mask] = 2 * arctan((y_grid[mask] / cos(incl)) \
                             / (x_grid[mask] \
                                + sqrt(x_grid[mask]**2 + (y_grid[mask] / cos(incl))**2))) - pi/2
    r = sqrt(x_grid**2 + (y_grid / cos(incl))**2) * dist # AU

    # Compute the line peak intensity

    peakint = zeros((npix, npix))

    sigma = fwhm / (2 * sqrt(2 * log(2)))
    peakint = int0 * exp(-r / (2 * sigma**2) / dist)

    if size:
        peakint[r > size] = 0.

    # Compute the disk velocity 

    vr = zeros((npix, npix))
    vtheta = zeros((npix, npix))
    mask = r >= rc
    if material_v == 'infalling':
        vr[mask] = sqrt(2 * G * mstar * M_sun / (r[mask] * au) - (G * mstar * M_sun * rc * au) / (r[mask] * au)**2)
    elif material_v == 'rotating':
        vtheta[mask] = sqrt(G * mstar * M_sun * rc * au) / (r[mask] * au)
    elif material_v == 'rotating+infalling' or material_v == 'fixed_size':
        vr[mask] = sqrt(2 * G * mstar * M_sun / (r[mask] * au) - (G * mstar * M_sun * rc * au) / (r[mask] * au)**2)
        vtheta[mask] = sqrt(G * mstar * M_sun * rc * au) / (r[mask] * au)
    elif material_v == 'rotating_only_in_rc':
        vr = vr
        vtheta = vtheta
    else:
        print 'I do not know the velocity profile of the material.'
    mask = (r < rc) * (r != 0)
    if material_v == 'infalling':
        vr[mask] = sqrt(2 * G * mstar * M_sun / (r[mask] * au) - (G * mstar * M_sun * rc * au) / (r[mask] * au)**2)
    elif material_v == 'rotating' or material_v == 'rotating_only_in_rc' or material_v == 'rotating+infalling' or material_v == 'fixed_size':
        vtheta[mask] = sqrt((G * mstar * M_sun) / (r[mask] * au)) # assume Keplerian rotation within rc
    else:
        print 'I do not know the velocity profile of the material.'
    
    # Compute its projection along the line of sight

    vproj = zeros((npix, npix))
    mask = r !=0
    vproj[mask] = (cos(theta[mask])*vr[mask] + sin(theta[mask])*vtheta[mask])*sin(incl)
    vproj *= 1e-3 # m/s -> km/s
    vproj += vlsr

    # Compute the synthetic datacube
    
    mask = r != 0
    sigma = zeros((npix, npix)) # local linewidth (FWHM/sqrt(8*ln(2))

    sigma = float(linewidth) # km/s
    intensity = zeros((nchan, npix, npix))
    intensity = peakint[newaxis,:,:] * exp(-(veloc[:,newaxis,newaxis]- vproj)**2 / (2 * sigma**2)) # Fixme: singularity at (0,0) !
    
    #---------------------------------------------------------------------------
    #We convolve with the beam, added since 29-9-2015 
    #---------------------------------------------------------------------------
    
    #Changes the FWHM to sigma.
    def FWHM_to_Sigma(FWHM):
        return FWHM/(2.*sqrt(2.*log(2.)))

    #First we calculate the smoothing kernel. 
    #So we need to know the major and minor axis of the beam in pixels.
    B_maj_pix_FWHM = 0.5*B_maj/abs(CDELT1)
    B_min_pix_FWHM = 0.5*B_min/abs(CDELT2)

    #Converting it to sigma's. 
    B_maj_pix_sigma = FWHM_to_Sigma(B_maj_pix_FWHM)   
    B_min_pix_sigma = FWHM_to_Sigma(B_min_pix_FWHM) 

    #The datacube which includes the thindisk data with some channels convolved
    #is created here. At the same time it is sliced in the correct size.
    intensity_convolved = intensity[:, NAXIS1/2-1 - Boxsize/2.:NAXIS1/2-1 + Boxsize/2., NAXIS2/2-1 - Boxsize/2.:NAXIS2/2-1 + Boxsize/2.]

    #Lets create an mask to select the channels we want to convolve. 
    mask = zeros(nchan, dtype = bool)
    mask[channel1:(channel2+1)] = True
    mask[channel3:(channel4+1)] = True
    
    #Selecting the channels that will be convolved. 
    channels_to_be_convolved = intensity_convolved[mask,:,:]
    veloc_to_be_convolved = veloc[mask]
    
    #Now we divide the intensity datacube in different parts so we can 
    #divide the work over multiple processes.
    import multiprocessing as mp
    
    #The amount of processors in the computer.
    N_processors = mp.cpu_count()
    #We want to use a maximum of 4 processors.
    if N_processors > 4:
        N_processors = 4
    
    #Creating the arrays containing the begin and end channels for the 
    #splitting.
    begin_channel = zeros(N_processors)
    end_channel = zeros(N_processors)
    
    #The number of channels that each of the processors will roughly 
    #receive. 
    nz, ny, nx = channels_to_be_convolved.shape
    number_of_channels = nz/N_processors
    
    for i in range(0,N_processors):
        begin_channel[i] = number_of_channels*i 
        end_channel[i] = number_of_channels*i + number_of_channels
        #To make sure we do all channels, we need this if statement. 
        if i == N_processors - 1:
            end_channel[i] = nz 

    pool = mp.Pool(processes = N_processors)
    
    #Convolving the datacube using all the processors on the computer.
    results = [pool.apply_async(Convolve_Datacube, args=(channels_to_be_convolved[begin_channel[x]:end_channel[x],:,:], B_maj_pix_sigma, B_min_pix_sigma, BPA, incl, veloc_to_be_convolved[begin_channel[x]:end_channel[x]], vlsr, x)) for x in range(0,N_processors)]

    #results only contains the pointers, so now we have to go and 'get' 
    #the arrays.
    temp = [p.get() for p in results]
    
    #Sorting the results so make sure we return them in the correct order.
    temp.sort()

    #Selecting the output.
    output = [r[1] for r in temp]
    #And putting it together in to one datacube. 
    channels_convolved = concatenate((output),axis = 0)
    
    #putting the convolved channels back in the datacube. 
    intensity_convolved[mask] = channels_convolved
    
    #Return the result. 
    return intensity_convolved

if __name__ == '__main__':
    Test_Datacube = Thindisk(1, 500, 1000, 5, 10, 100, 15, 40, 50, 65, 80, 'rotating+infalling')

    import matplotlib.pyplot as plt

    Nz,Ny,Nx = Test_Datacube.shape

    for i in range(0, Nz):
        plt.imshow(Test_Datacube[i,:,:], origin = 'lower')
        plt.colorbar()
        plt.title('i = ' + str(i))
        plt.show()
