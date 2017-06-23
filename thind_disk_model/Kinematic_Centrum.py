import numpy as np
from astropy import units as u
from astropy.io import fits
from Peak_PV_diagram import *
from astropy.coordinates import SkyCoord

def Degrees_to_Pixels(header, x_deg, y_deg, x_deg_err, y_deg_err):
    #Getting the necessary information from the header.
    CRPIX1 = header["CRPIX1"] 
    CDELT1 = header["CDELT1"] 
    CRVAL1 = header["CRVAL1"] 

    CRPIX2 = header["CRPIX2"]
    CDELT2 = header["CDELT2"] 
    CRVAL2 = header["CRVAL2"]

    #Calculating the difference between value defined in the header and 
    #x_center_deg 
    delta_x = x_deg - CRVAL1
    delta_y = y_deg - CRVAL2 

    #Now calculating the x and y coordinates in pixels.
    x_pix = CRPIX1 + delta_x/CDELT1
    y_pix = CRPIX2 + delta_y/CDELT2

    #Now calculating the errors on the x and y coordinates in degrees.
    x_pix_err = x_deg_err/np.abs(CDELT1)
    y_pix_err = y_deg_err/np.abs(CDELT2) 

    return x_pix, y_pix, x_pix_err, y_pix_err
    #End function

#Using the information in the header, this function transforms pixels 
#coordinates and their errors to degrees. 
def Pixels_to_Degrees(header, x_pix, y_pix, x_pix_err, y_pix_err):
    #Getting the necessary information from the header.
    CRPIX1 = header["CRPIX1"] 
    CDELT1 = header["CDELT1"] 
    CRVAL1 = header["CRVAL1"] 

    CRPIX2 = header["CRPIX2"]
    CDELT2 = header["CDELT2"] 
    CRVAL2 = header["CRVAL2"]  

    #Calculating the difference between value defined in the header and 
    #x_center_pix 
    delta_x = x_pix - CRPIX1
    delta_y = y_pix - CRPIX2  

    #Now calculating the x and y coordinates in degrees.
    x_deg = CRVAL1 + delta_x*CDELT1
    y_deg = CRVAL2 + delta_y*CDELT2

    #Now calculating the errors on the x and y coordinates in degrees.
    x_deg_err = x_pix_err*np.abs(CDELT1)
    y_deg_err = y_pix_err*np.abs(CDELT2)

    return x_deg, y_deg, x_deg_err, y_deg_err
    #End function

#This function determines the centrum of the M0 map (.fits file) using Gaussian
#fitting.  
#The following outputs can be chosen:
#deg = degrees
#pix = pixels 
#hmsdms = hmsdms
#The output always is (x, y, x_err, y_err)
#If Print_Output = True, all the coordinates will be printed. 
def Find_Centrum_M0_Map(M0_map, Output = 'deg', Print_Output = False):
    #Open the fits file. 
    hdu = fits.open(M0_map)[0]
    #Putting the data and header in seperate files.
    Header = hdu.header
    M0_Data = hdu.data
    
    #If the continuum is used it can have more dimensions than the 
    #program can handle. So we have to fix that.
    if len(M0_Data.shape) == 4:
		M0_Data = M0_Data[0,0,:,:]

    #---------------------------------------------------------------------------
    #Some useful information from the header

    #Pixel dimensions in degrees.
    PixelWidth_RA_deg = Header['CDELT1']
    PixelWidth_DEC_deg = Header['CDELT2']

    #Beam size
    Beam_x_deg = Header['BMAJ']
    Beam_y_deg = Header['BMIN']

    #Converting it to pixels
    Beam_x_pix = Beam_x_deg/np.abs(PixelWidth_RA_deg)
    Beam_y_pix = Beam_y_deg/PixelWidth_DEC_deg
    #---------------------------------------------------------------------------

    #First we going to find the center by Gaussian fitting routines. Then we 
    #will convert the resulting pixel coordinates to the WCS.

    #We will not use the entire picture for fitting, only the part around the 
    #central peak. 
    #Selecting the data that will be used for the fitting.   
    y_position_max, x_position_max  = np.unravel_index(M0_Data.argmax(), M0_Data.shape)

    print 'Coordinates max total object in pixels:'
    print 'x: ' + str(x_position_max)
    print 'y: ' + str(y_position_max)


    box_size = 25

    M0_Data[:y_position_max-box_size, :] = 0
    M0_Data[y_position_max+box_size:, :] = 0
    M0_Data[:,:x_position_max-box_size] = 0
    M0_Data[:, x_position_max+box_size:] = 0

    #Fitting the Gaussian to the M0 image.
    parameters_Gaussian_object = FitGauss2D(M0_Data)[0]
    x_pixel_centrum = parameters_Gaussian_object[2]
    y_pixel_centrum = parameters_Gaussian_object[1] 
    
    #Determing the errors on the coordinates of the centrum.
    x_pixel_err_centrum, y_pixel_err_centrum = Gaussian_Center_Position_Standard_Deviation(np.var(M0_Data), parameters_Gaussian_object[0], parameters_Gaussian_object[4], parameters_Gaussian_object[3], Beam_x_pix, Beam_y_pix)

    if Print_Output:
        print 'Coordinates center total object in pixels:'
        print 'x: ' + str(x_pixel_centrum)
        print 'y: ' + str(y_pixel_centrum)
        print 'The errors in the position determination:'
        print 'x error: ' + str(x_pixel_err_centrum)
        print 'y error: ' + str(y_pixel_err_centrum)
        print ' '

    #Now we convert the pixel coordinates to degrees.
    x_deg_centrum, y_deg_centrum, x_deg_err_centrum,  y_deg_err_centrum = Pixels_to_Degrees(Header, x_pixel_centrum, y_pixel_centrum, x_pixel_err_centrum, y_pixel_err_centrum)

    if Print_Output:
            print 'Coordinates center total object in degrees:'
            print 'x: ' + str(x_deg_centrum)
            print 'y: ' + str(y_deg_centrum)
            print 'The errors in the position determination:'
            print 'x error: ' + str(x_deg_err_centrum)
            print 'y error: ' + str(y_deg_err_centrum)
            print ' '

    #Now we convert the coordinates from degrees to hours:minutes:seconds and 
    #degrees:arcminutes:arcseconds

    coords = SkyCoord(ra=x_deg_centrum*u.degree, dec=y_deg_centrum*u.degree)
    coords_err = SkyCoord(ra=x_deg_err_centrum*u.degree, dec=y_deg_err_centrum*u.degree)
    if Print_Output:
        print 'Coordinates center total object in hms and dms:'
        print coords.to_string('hmsdms')
        print 'The errors in the position determination in hms and dms:'
        print coords_err.to_string('hmsdms')
        print ' '

    if Output == 'pix':
        return x_pixel_centrum, y_pixel_centrum, x_pixel_err_centrum, y_pixel_err_centrum

    if Output == 'deg':
        return x_deg_centrum, y_deg_centrum, x_deg_err_centrum, y_deg_err_centrum 

    if Output == 'hmsdms':
        return coords, coords_err
    #End function

if __name__ == '__main__':
    print ' '
    print 'C18O'
    directory = 'C18O_1/Moments_Images/'
    File = 'Moment_0_map_L1165-C18O_1.fits'
    test = directory + File 
    Find_Centrum_M0_Map(test, Print_Output = True)
    print '13CO'
    directory = '13CO_1/Moments_Images/'
    File = 'Moment_0_map_L1165-13CO_1.fits'
    test = directory + File 
    Find_Centrum_M0_Map(test, Print_Output = True)
    print 'HCO+'
    directory = 'HCO+_1/Moments_Images/'
    File = 'Moment_0_map_L1165-HCO+_1.fits'
    test = directory + File 
    Find_Centrum_M0_Map(test, Print_Output = True)
    print 'Optimization 1'
    directory = 'Optimizations/Optimization_1/Moments_Images/'
    File = 'Moment_0_map_Optimization_1_convolved_with_beam.fits'
    test = directory + File 
    Find_Centrum_M0_Map(test, Print_Output = True)
    print 'test 1'
    directory = 'Thin_Disk_Model_1/Moments_Images/'
    File = 'Moment_0_map_test_1_convolved_with_beam.fits'
    test = directory + File 
    Find_Centrum_M0_Map(test, Print_Output = True)

'''
--------------------------------------------------------------------------------
Unused Code:
--------------------------------------------------------------------------------

        mask_center = circle_mask(M0_Data, x_pixel_centrum, y_pixel_centrum, 2)# np.sqrt(x_pixel_err_centrum**2 + y_pixel_err_centrum**2)) 

        plt.imshow(M0_Data*~mask_center, origin = 'lower')
        plt.colorbar()
        plt.show()

    def circle_mask(im, xc, yc, rcirc):
        ny, nx = im.shape
        y,x = np.mgrid[0:nx,0:ny]
        r = np.sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc))
        return ( (r < rcirc))

    #Creating a circle mask around the desired coordinates.
    mask = circle_mask(M0_Data, x_position_max, y_position_max, box_size)

    #Applying the mask to the data.
    M0_Data *= mask

'''
