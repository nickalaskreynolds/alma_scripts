# thindisk.py --- Compute a synthetic datacube for thin disk

#Changes made to code since 4-12-2015 by Steven Bos.
#Include:
#The channel velocities are now based on the velocity of the first channel,
#which has to be given. 
#For incl > 85 degrees the code also convolves with a Gaussian beam, its
#width dependent on the velocity, so make the disk thicker. 
#The intensity datacube is convolved with the beam of the telescope. 
#Now also uses multiple processors for the convolving. 

import argparse
import ConfigParser
from numpy import *
from astropy.constants import G, M_sun, au
from astropy import units as u
from astropy.io import fits
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
        
def main(File):
    # Read the parameters file
    params = ConfigParser.ConfigParser()
    params.read(File)
    #---------------------------------------------------------------------------
    #Info about the source.
    mstar = params.getfloat("disk", "mstar")
    rc = params.getfloat("disk", "rc")
    size = params.getfloat("disk", "size")
    incl = params.getfloat("disk", "incl")
    pa = params.getfloat("disk", "pa")
    dist = params.getint("disk", "dist")
    lineint = params.get("line", "intensity").split(",")
    linewidth = params.get("line", "width")
    frequency = params.getfloat("line", "frequency")
    vlsr = params.getfloat("line", "vlsr")
    #The center of the disk in pixels.
    xcen = params.getint("disk", "xcen")
    ycen = params.getint("disk", "ycen")
    #Info for the beam.
    B_maj = params.getfloat("cube", "Bmaj")
    B_min = params.getfloat("cube", "Bmin")
    BPA = params.getfloat("cube", "BPA")

    #Info for the header.
    NAXIS   = params.getint("cube", "NAXIS")                                               
    NAXIS1  = params.getint("cube", "NAXIS1")                                              
    NAXIS2  = params.getint("cube", "NAXIS2")                                                 
    NAXIS3  = params.getint("cube", "NAXIS3")                                             
    NAXIS4  = params.getint("cube", "NAXIS4")                                                             
                                                      
    CRPIX1  = params.getfloat("cube", "CRPIX1")                                        
    CDELT1  = params.getfloat("cube", "CDELT1")                                                     
    CRVAL1  = params.getfloat("cube", "CRVAL1")                                                    
                                                         
    CRPIX2  = params.getfloat("cube", "CRPIX2")                                                      
    CDELT2  = params.getfloat("cube", "CDELT2")                                                    
    CRVAL2  = params.getfloat("cube", "CRVAL2")                                                     
                                                        
    CRPIX3  = params.getfloat("cube", "CRPIX3")                                                    
    CDELT3  = params.getfloat("cube", "CDELT3")                                                       
    CRVAL3  = params.getfloat("cube", "CRVAL3")                                                    
                                                        
    CRPIX4  = params.getfloat("cube", "CRPIX4")                                                 
    CDELT4  = params.getfloat("cube", "CDELT4")                                                    
    CRVAL4  = params.getfloat("cube", "CRVAL4")                                                  
                                                                                            
    LWIDTH  = params.getfloat("cube", "LWIDTH")                                                   
    LSTEP   = params.getfloat("cube", "LSTEP")                                                
    LSTART  = params.getfloat("cube", "LSTART")    
  
    #Output name.
    fitsname = params.get("output", "name")
    #---------------------------------------------------------------------------
    #Some important variables. 
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
    if lineint[0] == "powerlaw":
        int_r1, r1, int_expn = map(lambda x: float(x), lineint[1:4]) 
        mask = r != 0
        peakint[mask] = int_r1 * pow(r[mask] / dist / r1, int_expn)
    elif lineint[0] == "gaussian":
        int0, fwhm = map(lambda x: float(x), lineint[1:3])
        sigma = fwhm / (2 * sqrt(2 * log(2)))
        peakint = int0 * exp(-r / (2 * sigma**2) / dist)
    else:
        int_ring, r1, r2 = map(lambda x: float(x), lineint[1:4]) # ring
        mask = ((r / dist) > r1) * ((r / dist) < r2)
        peakint[mask] = int_ring
    if size:
        peakint[r > size] = 0.

    # Compute the disk velocity 
    vr = zeros((npix, npix))
    vtheta = zeros((npix, npix))
    mask = r >= rc
    vr[mask] = sqrt(2 * G * mstar * M_sun / (r[mask] * au) - (G * mstar * M_sun * rc * au) / (r[mask] * au)**2)
    vtheta[mask] = sqrt(G * mstar * M_sun * rc * au) / (r[mask] * au)
    mask = (r < rc) * (r != 0)
    vtheta[mask] = sqrt((G * mstar * M_sun) / (r[mask] * au)) # assume Keplerian rotation within rc
    
    # Compute its projection along the line of sight
    vproj = zeros((npix, npix))
    mask = r !=0
    
    vproj[mask] = (cos(theta[mask])*vr[mask] + sin(theta[mask])*vtheta[mask])*sin(incl)
    vproj *= 1e-3 # m/s -> km/s
    vproj += vlsr
    
    '''
    vproj_1 = zeros((npix, npix))
    vproj_2 = zeros((npix, npix))
    
    #First we make an projection of all radial velocities.
    vproj_1[mask] = cos(theta[mask])*vr[mask]*sin(incl)*1e-3 + vlsr
    #Then we make an projection of all tangential velocities.
    vproj_2[mask] = sin(theta[mask])*vtheta[mask]*sin(incl)*1e-3 + vlsr
    '''
    
    # Compute the synthetic datacube
    mask = r != 0
    sigma = zeros((npix, npix)) # local linewidth (FWHM/sqrt(8*ln(2))
    if "*vkep" in linewidth:
        linewidth = float(linewidth.split("*vkep", 2)[0])
        sigma[mask] = linewidth * sqrt(G * mstar * M_sun / (r[mask] * au)) # fraction of the Keplerian velocity
        sigma *= 1e-3 # m/s -> km/s
    else:
        sigma = float(linewidth) # km/s
    
    intensity = zeros((nchan, npix, npix))
    intensity = peakint[newaxis,:,:] * exp(-(veloc[:,newaxis,newaxis]- vproj)**2 / (2 * sigma**2)) # Fixme: singularity at (0,0) !
    '''
    intensity_1 = zeros((nchan, npix, npix))
    intensity_2 = zeros((nchan, npix, npix))
    
    #We create separate datacubes for both the radial and tangential 
    #datacubes.
    intensity_1 = peakint[newaxis,:,:] * exp(-(veloc[:,newaxis,newaxis]- vproj_1)**2 / (2 * sigma**2)) # Fixme: singularity at (0,0) !
    intensity_2 = peakint[newaxis,:,:] * exp(-(veloc[:,newaxis,newaxis]- vproj_2)**2 / (2 * sigma**2)) # Fixme: singularity at (0,0) !
    
    #We then add the two datacubes.
    intensity = intensity_1 + intensity_2
    '''
    #---------------------------------------------------------------------------
    #Test plotting
    #---------------------------------------------------------------------------
    '''
    import matplotlib.pyplot as plt
    
    plt.figure(1)
    plt.imshow(ra_grid, origin = 'lower')
    plt.title('RA grid')
    plt.colorbar()


    plt.figure(2)
    plt.imshow(dec_grid, origin = 'lower')
    plt.title('DEC grid')
    plt.colorbar()

    plt.figure(3)
    plt.imshow(x_grid, origin = 'lower')
    plt.title('x grid')
    plt.colorbar()

    plt.figure(4)
    plt.imshow(y_grid, origin = 'lower')
    plt.title('y grid')
    plt.colorbar()

    plt.figure(5)
    plt.imshow(r, origin = 'lower')
    plt.title('r')
    plt.colorbar()

    plt.figure(6)
    plt.imshow(theta/pi, origin = 'lower')
    plt.title('theta')
    plt.colorbar()

    plt.figure(7)
    plt.imshow(peakint, origin = 'lower')
    plt.title('Intensity')
    plt.colorbar()


    plt.figure(8)
    plt.imshow(vr, origin = 'lower')
    plt.title('radial velocity')
    plt.colorbar()

    plt.figure(9)
    plt.imshow(vtheta, origin = 'lower')
    plt.title('tangetial velocity')
    plt.colorbar()

    plt.figure(10)
    plt.imshow(vproj, origin = 'lower')
    plt.title('projected velocity')
    plt.colorbar()
    plt.show()
    '''
    #---------------------------------------------------------------------------
    
    print('Saving the thin disk model without convolving.')
    #First we save the non convolved disk.    
    hdu = fits.PrimaryHDU(intensity)
    hdr = hdu.header
    hdr["CRPIX1"] = CRPIX1
    hdr["CDELT1"] = CDELT1
    hdr["CRVAL1"] = CRVAL1
    hdr["CTYPE1"] = "RA---SIN"
    hdr["CRPIX2"] = CRPIX2
    hdr["CDELT2"] = CDELT2
    hdr["CRVAL2"] = CRVAL2
    hdr["CTYPE2"] = "DEC--SIN"
    hdr["CRPIX3"] = CRPIX3
    hdr["CDELT3"] = CDELT3
    hdr["CRVAL3"] = CRVAL3
    hdr["CTYPE3"] = "VELO-LSR"
    hdr["CRPIX4"] = CRPIX4
    hdr["CDELT4"] = CDELT4
    hdr["CRVAL4"] = CRVAL4
    hdr["CTYPE4"] = "STOKES" 
    hdr["CUNIT1"] = "DEG"
    hdr["CUNIT2"] = "DEG"
    hdr["CUNIT3"] = "M/S"
    hdr["BUNIT"] = "JY/BEAM"
    hdr["RESTFREQ"] = frequency*1e6
    hdr["BMAJ"] = B_maj
    hdr["BMIN"] = B_min
    hdr["BPA"] = BPA
    hdr["LWIDTH"] = LWIDTH
    hdr["LSTEP"] = LSTEP
    hdr["LSTART"] = LSTART
    hdr["NAXIS"] = NAXIS
    hdr["NAXIS1"] = NAXIS1
    hdr["NAXIS2"] = NAXIS2
    hdr["NAXIS3"] = NAXIS3
    hdr["NAXIS4"] = NAXIS4 
    hdu.writeto("%s.fits" % fitsname, clobber = True)

    #---------------------------------------------------------------------------
    #First we convolve with the beam, added since 29-9-2015 
    #---------------------------------------------------------------------------
    def FWHM_to_Sigma(FWHM):
        #Changes the FWHM to sigma.
        return FWHM/(2.*sqrt(2.*log(2.)))

    #First we calculate the smoothing kernel. 
    #So we need to know the major and minor axis of the beam in pixels.
    B_maj_pix_FWHM = 0.5*B_maj/abs(CDELT1)
    B_min_pix_FWHM = 0.5*B_min/abs(CDELT2)

    #Converting it to sigma's. 
    B_maj_pix_sigma = FWHM_to_Sigma(B_maj_pix_FWHM)   
    B_min_pix_sigma = FWHM_to_Sigma(B_min_pix_FWHM) 

    #Now we divide the intensity datacube in different parts so we can 
    #divide the work over multiple processes.
    import multiprocessing as mp
    
    #The amount of processors in the computer.
    N_processors = mp.cpu_count()
    
    #Creating the arrays containing the begin and end channels for the 
    #splitting.
    begin_channel = zeros(N_processors)
    end_channel = zeros(N_processors)
    
    #The number of channels that each of the processors will roughly 
    #receive. 
    number_of_channels = nchan/N_processors
    
    #Here we decide which channels go to which processor.
    for i in range(0,N_processors):
        begin_channel[i] = number_of_channels*i 
        end_channel[i] = number_of_channels*i + number_of_channels
        #To make sure we do all channels, we need this if statement. 
        #We basically give all but one of the processors the same amount 
        #of channels to convolve. The rest is given to one processor.  
        if i == N_processors - 1:
            end_channel[i] = nchan 

    pool = mp.Pool(processes = N_processors)
    
    #Convolving the datacube using all the processors on the computer.
    print 'Starting the convolving.' 
    results = [pool.apply_async(Convolve_Datacube, args=(intensity[begin_channel[x]:end_channel[x],:,:], B_maj_pix_sigma, B_min_pix_sigma, BPA, incl, veloc[begin_channel[x]:end_channel[x]], vlsr, x)) for x in range(0,N_processors)]

    #results only contains the pointers, so now we have to go and 'get' 
    #the arrays.
    temp = [p.get() for p in results]
    
    #Sorting the results so make sure we return them in the correct order.
    temp.sort()

    #Selecting the output.
    output = [r[1] for r in temp]
    #And putting it together in to one datacube. 
    intensity_convolved = concatenate((output),axis = 0)
    
    #---------------------------------------------------------------------------
    # Export the datacube to FITS

    print('Saving the thin disk model with convolving.')    
    hdu = fits.PrimaryHDU(intensity_convolved)
    hdr = hdu.header
    hdr["CRPIX1"] = CRPIX1
    hdr["CDELT1"] = CDELT1
    hdr["CRVAL1"] = CRVAL1
    hdr["CTYPE1"] = "RA---SIN"
    hdr["CRPIX2"] = CRPIX2
    hdr["CDELT2"] = CDELT2
    hdr["CRVAL2"] = CRVAL2
    hdr["CTYPE2"] = "DEC--SIN"
    hdr["CRPIX3"] = CRPIX3
    hdr["CDELT3"] = CDELT3
    hdr["CRVAL3"] = CRVAL3
    hdr["CTYPE3"] = "VELO-LSR"
    hdr["CRPIX4"] = CRPIX4
    hdr["CDELT4"] = CDELT4
    hdr["CRVAL4"] = CRVAL4
    hdr["CTYPE4"] = "STOKES" 
    hdr["CUNIT1"] = "DEG"
    hdr["CUNIT2"] = "DEG"
    hdr["CUNIT3"] = "KM/S"
    hdr["BUNIT"] = "JY/BEAM"
    hdr["RESTFREQ"] = frequency*1e6
    hdr["BMAJ"] = B_maj
    hdr["BMIN"] = B_min
    hdr["BPA"] = BPA
    hdr["LWIDTH"] = LWIDTH
    hdr["LSTEP"] = LSTEP
    hdr["LSTART"] = LSTART
    hdr["NAXIS"] = NAXIS
    hdr["NAXIS1"] = NAXIS1
    hdr["NAXIS2"] = NAXIS2
    hdr["NAXIS3"] = NAXIS3
    hdr["NAXIS4"] = NAXIS4 
    hdu.writeto("%s_convolved_with_beam.fits" % fitsname, clobber = True)

if __name__ == '__main__':
    main('Test_1.ini')
