import numpy as np
import pyfits as pf
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as const
from astropy.coordinates import SkyCoord
from astropy.io import fits
from wcsaxes import WCS

#Make the plots nice
import Nice_Plots_2
Nice_Plots_2.set_style()
    
def Save_Important_Header_Data(header):
    '''
    This functions gets an fits file header and saves the information that is 
    important for other actions. The Name list contains all the names of the data
    to saved, while Data contains the actual data.
    header = [header] 
    '''
    #The header data we want to save.
    Names = ['NAXIS', 'NAXIS1', 'NAXIS2', 'CRPIX1', 'CDELT1', 'CRVAL1',
                      'CTYPE1', 'CRPIX2', 'CDELT2', 'CRVAL2', 'CTYPE2',
                      'BUNIT', 'BMAJ', 'BMIN', 'BPA']
    #Creating the list where we will save the data.
    Data = []
    #Saving the data. 
    i = 0
    while i < len(Names):
        Data.append(header[Names[i]])
        i += 1
        #End while

    #As we only have two axis left, we have to change NAXIS.
    Data[0] = 2

    return Names, Data

def Create_Continuum_Contour(Header_Data_Cube):
    File = 'L1527_Cont_2.fits'
    Directory = 'Cont_2/'
    
    #Opening the file.
    hdu = fits.open(Directory + File)[0]
    #Putting the data and header in seperate files.
    header = hdu.header
    Data = hdu.data
    Data = Data[0,0,:,:]
    
    #Calculating the std of the image.
    Std_Data = np.std(Data[0:100,0:100])
    Contour_Levels = np.array([5, 20, 40, 60, 80, 100])*Std_Data#[5,15, 25, 35, 45]
        
    #Header information.
    CRVAL_RA_Cont = header['CRVAL1']                                             
    CDELT_RA_Cont = header['CDELT1']                                        
    CRPIX_RA_Cont = header['CRPIX1']  
    
    CRVAL_DEC_Cont = header['CRVAL2']                                             
    CDELT_DEC_Cont = header['CDELT2']                                        
    CRPIX_DEC_Cont = header['CRPIX2']
    N2_Cont = header['NAXIS2']  
    
    CRVAL_RA_Cube = Header_Data_Cube['CRVAL1']                                             
    CDELT_RA_Cube = Header_Data_Cube['CDELT1']                                        
    CRPIX_RA_Cube = Header_Data_Cube['CRPIX1']  
    
    CRVAL_DEC_Cube = Header_Data_Cube['CRVAL2']                                             
    CDELT_DEC_Cube = Header_Data_Cube['CDELT2']                                        
    CRPIX_DEC_Cube = Header_Data_Cube['CRPIX2']
    N2_Cube = Header_Data_Cube['NAXIS2'] 
  
    #For simplicity we assume that the continuum always has a higher 
    #resolution than the data. Meaning we can bin the continuum data.
    #We also assume CDELT_DEC_Cont == CDELT_RA_Cont.
    if CDELT_DEC_Cont != CDELT_DEC_Cube or CDELT_RA_Cont != CDELT_RA_Cube:
        print 'CDELTs are not the same in the continuum and the data!'
        print '\nCDELT DEC cont: ' + str(CDELT_DEC_Cont)
        print 'CDELT DEC data: ' + str(CDELT_DEC_Cube)
        print '\nCDELT RA cont: ' + str(CDELT_RA_Cont)
        print 'CDELT RA data: ' + str(CDELT_RA_Cube)
        print '\nRescaling...'
        
        from numpy.lib.stride_tricks import as_strided
        #Found this function on the internet.
        def strided_rescale(g, bin_fac):
            strided = as_strided(g,
                shape=(g.shape[0]//bin_fac, g.shape[1]//bin_fac, bin_fac, bin_fac),
                strides=((g.strides[0]*bin_fac, g.strides[1]*bin_fac)+g.strides))
            return strided.mean(axis=-1).mean(axis=-1) 
        
        BinFactor = np.round(CDELT_DEC_Cube/CDELT_DEC_Cont,0)
        #Binning the data.
        Data = strided_rescale(Data, BinFactor)
        print 'Bin factor is: ' + str(BinFactor)
        
        #Rescaling the DELTs
        CDELT_DEC_Cont = CDELT_DEC_Cont*BinFactor
        CDELT_RA_Cont = CDELT_RA_Cont*BinFactor

        #Rescaling the pixel values.
        CRPIX_DEC_Cont = CRPIX_DEC_Cont/BinFactor
        CRPIX_RA_Cont = CRPIX_RA_Cont/BinFactor
        
        #Also the amount of pixels along axis 2.
        N2_Cont = N2_Cont/BinFactor
            
    #If the continuum is smaller pixel wise than the data we have to pad
    #instead of crop. 
    if N2_Cont < N2_Cube:
        print 'Padding the continuum data.'
        #First we fill the final image with zeros. 
        Final_Data = np.zeros((N2_Cube, N2_Cube))

        #We also take into account different coordinates for the 
        #reference pixels.
        Delta_DEC = np.round(CRPIX_DEC_Cube - (CRPIX_DEC_Cont - (CRVAL_DEC_Cont - CRVAL_DEC_Cube)/CDELT_DEC_Cube), 0)
        Delta_RA = np.round(CRPIX_RA_Cube - (CRPIX_RA_Cont - (CRVAL_RA_Cont - CRVAL_RA_Cube)/CDELT_RA_Cube), 0)
        offset = 0
        #Putting the data in the new bigger file.
        Final_Data[Delta_DEC:Delta_DEC+N2_Cont, Delta_RA + offset:Delta_RA+N2_Cont + offset] = Data
    #Now we crop, if necessary. 
    else:
        print 'Cropping the continuum data.'
        #The degrees where the edges of the data cube are.
        Deg_Left_Cube = CRVAL_RA_Cube + (0-CRPIX_RA_Cube)*CDELT_RA_Cube
        Deg_Right_Cube = CRVAL_RA_Cube + (N2_Cube-CRPIX_RA_Cube)*CDELT_RA_Cube
        Deg_Down_Cube = CRVAL_DEC_Cube + (0-CRPIX_DEC_Cube)*CDELT_DEC_Cube
        Deg_Up_Cube = CRVAL_DEC_Cube + (N2_Cube-CRPIX_DEC_Cube)*CDELT_DEC_Cube

        #The pixel coordinates in the cont image.
        Pix_Left_Cont = (Deg_Left_Cube - CRVAL_RA_Cont)/CDELT_RA_Cont + CRPIX_RA_Cont
        Pix_Right_Cont = (Deg_Right_Cube - CRVAL_RA_Cont)/CDELT_RA_Cont + CRPIX_RA_Cont
        Pix_Up_Cont = (Deg_Up_Cube - CRVAL_DEC_Cont)/CDELT_DEC_Cont + CRPIX_DEC_Cont
        Pix_Down_Cont = (Deg_Down_Cube - CRVAL_DEC_Cont)/CDELT_DEC_Cont + CRPIX_DEC_Cont

        #Cropping
        Final_Data = Data[Pix_Down_Cont:Pix_Up_Cont, Pix_Left_Cont:Pix_Right_Cont]

    return Final_Data, Contour_Levels

def Image_Moments(File, Directory, begin, end, Object, Line, Zoom = False, boxsize = 100, vmin_M1 = 0, vmax_M1 = 5, vmin_M2 = 0, vmax_M2 = 5, Thindisk = False):
    '''
    This function creates several image moments of the given .fits file in the 
    given directory. The following variables must be given:
    File: File name of the fits file, including .fits.
    Directory: Directory where the file is located.
    begin: Begin channel, not all channels contain usefull data and this is the 
    first one the programm should look at.
    end: The last channel it should look at, remember the first channel is 0 in the
    array.
    Object: the object name
    Line: in which molecular line (13CO, C18O, etc) it has been observed.
    Zoom: if true it zooms in on the center of the image with total pixel width of
    boxsize. 
    Thindisk: true for data that has been generated with the thindisk model, 
    because the data structure is slightly different.  
    '''
    print(' ')  
    print('Creating the image moments for ' + File[:-5])
    #---------------------------------------------------------------------------
    #Getting the data from the fits file.
    #---------------------------------------------------------------------------
    #Read the data.
    hdu = fits.open(Directory + File)[0]
        
#    print hdu.header.shape
#    print hdu.data.shape
    #For the wcs plotting we have to change the amount of axis.
    hdu.header['NAXIS'] = 2

    #Also we need to remove information about the third and fourth axis.
    #We have one exception, for some reason the thindisk model does not include 
    #NAXIS4 but does include CUNIT3. So we have to account for that.
    if Thindisk:
        hd=['NAXIS3','CTYPE3','CRVAL3','CDELT3','CRPIX3','CTYPE4','CRVAL4','CDELT4','CRPIX4', 'CUNIT3']
    else:
        hd=['NAXIS3','NAXIS4','CTYPE3','CRVAL3','CDELT3','CRPIX3','CTYPE4','CRVAL4','CDELT4','CRPIX4', 'PC03_01', 'PC04_01', 'PC03_02', 'PC04_02', 'PC01_03', 'PC02_03', 'PC03_03', 'PC04_03', 'PC01_04', 'PC02_04', 'PC03_04', 'PC04_04', 'CUNIT3', 'CUNIT4']
    for i in range(len(hd)):
        hdu.header.remove(hd[i])
        
    #print hdu.header
    #Now wcs knows how our coordinate system works.
    wcs = WCS(hdu.header)
    
    #We reload the data, because for the rest of our programm we need the 
    #information about the third and fourth axis.
    hdu = fits.open(Directory + File)[0]
    #Putting the data and header in seperate files.
    header = hdu.header
    Data = hdu.data

    #Saving the header data we want to keep. 
    Header_saved_Names, Header_saved_Data = Save_Important_Header_Data(header)

    #Size in pixels of each channel.
    Nx = header['NAXIS1']
    Ny = header['NAXIS2']

    #Determining the begin_v, delta_v and the amount of images N.
    N = header['NAXIS3']
    if header['CTYPE3']  == 'VELO-LSR':
        print 'True'
        begin_v = header['LSTART']
        delta_v = header['LWIDTH']
    elif header['CTYPE3']  == 'FREQ':
        #Reading the data in frequencies. We have to bring this to velocities. 
        begin_freq = header['CRVAL3']
        delta_freq = header['CDELT3']
        begin_pos = header['CRPIX3'] - 1
        rest_freq = header['RESTFRQ']

        #The speed of light is.
        c = const.c.value/1000#km s^-1

        #Calculating the begin velocity.
        begin_v = c * (np.square(rest_freq) - np.square(begin_freq - delta_freq*begin_pos)) / ( np.square(rest_freq) + np.square(begin_freq - delta_freq*begin_pos))
        #Now we calculate the delta v
        begin_v_plus_one = c * (np.square(rest_freq) - np.square(begin_freq - delta_freq*(begin_pos + 1))) / ( np.square(rest_freq) + np.square(begin_freq - delta_freq*(begin_pos + 1)))
        delta_v = np.round(begin_v - begin_v_plus_one, 2)

    #The used units.
    Image_Units = header['BUNIT']

    #Pixelsizes
    x_pixel_size = header['CDELT1']
    y_pixel_size = header['CDELT2']

    #Creating the data for beam size in the plots.
    major_axis = header['BMAJ']
    minor_axis = header['BMIN']
    
    print 'Beam properties: '
    print 'major axis: ' + str(major_axis*3600)
    print 'minor axis: ' + str(minor_axis*3600) 
    
    #The position angle of the beam. 
    Beam_Position_angle_deg = header['BPA']
    #Changing the position angle to radians and correcting for the different 
    #definition of the angle. 
    Beam_Position_angle_rad = np.radians(Beam_Position_angle_deg + 90)
    theta = np.arange(0.0, 360.0, 1.0)*np.pi/180.0
    #Also adjusting for the pixel size. 
    x_beam_temp = 0.5*major_axis/x_pixel_size*np.cos(theta)
    y_beam_temp = 0.5*minor_axis/y_pixel_size*np.sin(theta)
    #Now we rotate using the rotation matrix. 
    x_beam = np.cos(Beam_Position_angle_rad)*x_beam_temp - np.sin(Beam_Position_angle_rad)*y_beam_temp
    y_beam = np.sin(Beam_Position_angle_rad)*x_beam_temp + np.cos(Beam_Position_angle_rad)*y_beam_temp

    #The position of the center.
    if Zoom:
        x_beam += (Nx - boxsize)/2 + boxsize/10
        y_beam += (Ny - boxsize)/2 + boxsize/10    
    else:
        x_beam += Nx/10
        y_beam += Ny/10 

    #---------------------------------------------------------------------------
    #Calculating the moments
    #---------------------------------------------------------------------------
    
    if not Thindisk:
        Data = Data[0,:,:,:]
        
    mean_image = np.mean(Data[begin:end,:,:], 0)
    median_image = np.median(Data[begin:end,:,:], 0)
    max_image = np.max(Data[begin:end,:,:], 0)
    min_image = np.min(Data[begin:end,:,:], 0)
    sum_image = np.sum(Data[begin:end,:,:], 0)
    integrated_image = sum_image*delta_v

    #Now calculate the integrated image that will be used for the other
    #(velocity maps, velocity dispersion map) calculations.
    integrated_image_calculations = 0*mean_image
    limit = 3#Sigma
    i=begin
    while i < end + 1:
        std_data = np.std(Data[i,:,:])        
        Temp = Data[i,:,:]*1
        Temp[Temp < limit*std_data] = 0  
        integrated_image_calculations += Temp*delta_v
        i += 1

    #Now calculate the velocity field/intensity weighted coordinate.
    temp_image = 0*mean_image
    i=begin
    while i < end + 1:
        std_data = np.std(Data[i,:,:])        
        Temp = Data[i,:,:]*1
        Temp[Temp < limit*std_data] = 0  
        temp_image += Temp*(begin_v+i*delta_v)*delta_v
        i += 1
    velocity_field_image = np.divide(temp_image, integrated_image_calculations)

    #Now calculate the the intensity weighted dispersion of the coordinate/
    #veloctity dispersion fields.
    temp_image = 0*mean_image
    i = begin
    while i < end + 1: 
        std_data = np.std(Data[i,:,:])        
        Temp = Data[i,:,:]*1
        Temp[Temp < limit*std_data] = 0
        temp_image += Temp*(begin_v+i*delta_v - velocity_field_image)**2*delta_v
        i += 1
    velocity_dispersion_field_image = np.sqrt(np.divide(temp_image, np.abs(integrated_image_calculations)))

    #Now calculate the standard deviation about the mean of the spectrum.
    temp_image = 0*mean_image
    i = begin
    while i < end + 1:    
        temp_image += (Data[i,:,:] - mean_image)**2
        i += 1
    standard_deviation_image = np.sqrt(temp_image/(N-1))

    #Now calculate the root mean square of the spectrum.
    temp_image = 0*mean_image
    i = begin
    while i < end + 1:    
        temp_image += (Data[i,:,:])**2
        i += 1
    root_mean_square_image = np.sqrt(temp_image/N)

    #Now calculate the absolute mean deviation of the spectrum.
    temp_image = 0*mean_image
    i = begin
    while i < end + 1:    
        temp_image += abs(Data[i,:,:] - mean_image)
        i += 1
    absolute_mean_deviation_image = temp_image/N

    #---------------------------------------------------------------------------
    #Contour lines
    #---------------------------------------------------------------------------

    #For the max images we want contours based on percentages of the max. 
    Max_Value = np.max(max_image)
    Max_Contour_Levels = np.arange(0.1, 1, 0.1)*Max_Value 

    #For the min images we want contours based on percentages of the max. 
    Min_Value = np.min(min_image)
    Min_Contour_Levels = np.arange(0.1, 1, 0.1)*Min_Value

    #Contour for mean
    std_mean = np.std(mean_image) 
    Mean_Contour_Levels = np.arange(3, 11, 2)*std_mean 

    #Contour for median
    std_median = np.std(median_image) 
    Median_Contour_Levels = np.arange(3, 11, 2)*std_median 

    #Contour for sum
    std_sum = np.std(sum_image) 
    Sum_Contour_Levels = np.arange(3, 11, 2)*std_sum

    #Contour for integrated
    std_integrated = np.std(integrated_image) 
    Integrated_Contour_Levels = np.arange(3, 11, 2)*std_integrated 

    Shape_Images = np.shape(integrated_image)
    #---------------------------------------------------------------------------
    #Plotting
    #---------------------------------------------------------------------------
    #Creating the directory for saving the images, if it does not already exists.
    save_directory = Directory + 'Moments_Images/' 
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
    
    print('Plotting the images.')    
    #---------------------------------------------------------------------------
    #Plotting the mean image.
    #---------------------------------------------------------------------------
    '''
    print 'Plotting the mean image.'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = wcs)  
    cax = ax.imshow(mean_image)
    ax.contour(mean_image, Mean_Contour_Levels, colors = 'white')
    ax.fill(x_beam, y_beam, alpha=1, facecolor='white', edgecolor='white', linewidth=1, zorder=2)

    ax.invert_yaxis()

    cbar = fig.colorbar(cax)
    cbar.set_label(Image_Units)

    lon = ax.coords['ra']
    lat = ax.coords['dec']

    ax.set_title('Mean image of ' + Object + ' in ' + Line)
    lon.set_axislabel('J2000 Right Ascension')
    lat.set_axislabel('J2000 Declination')

    lon.set_major_formatter('hh:mm:ss')
    lat.set_major_formatter('dd:mm:ss')

    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    lat.set_minor_frequency(10)

    ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
    
    if Zoom:
        ax.set_xlim((Nx - boxsize)/2, (Nx + boxsize)/2)
        ax.set_ylim((Ny - boxsize)/2, (Ny + boxsize)/2)

    plt.savefig(save_directory + 'mean_image_' + File[:-5] + '.jpg',bbox_inches='tight')
    plt.savefig(save_directory + 'mean_image_' + File[:-5] + '.pdf',bbox_inches='tight')
    plt.close()

    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(mean_image)
    HDR = HDU.header
    #Saving the correct header data. 
    i = 0
    while i < len(Header_saved_Names):
        HDR[Header_saved_Names[i]] = Header_saved_Data[i]
        i += 1
        #End while
    HDU.writeto(save_directory + 'mean_image_' + File[:-5] + '.fits', clobber = True)
    #---------------------------------------------------------------------------
    #Plotting the median image.
    #---------------------------------------------------------------------------
    print 'Plotting the median image.'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = wcs)  
    cax = ax.imshow(median_image)
    ax.contour(median_image, Median_Contour_Levels, colors = 'white')
    ax.fill(x_beam, y_beam, alpha=1, facecolor='white', edgecolor='white', linewidth=1, zorder=2)

    ax.invert_yaxis()

    cbar = fig.colorbar(cax)
    cbar.set_label(Image_Units)

    lon = ax.coords['ra']
    lat = ax.coords['dec']

    ax.set_title('Median image of ' + Object + ' in ' + Line)
    lon.set_axislabel('J2000 Right Ascension')
    lat.set_axislabel('J2000 Declination')

    lon.set_major_formatter('hh:mm:ss')
    lat.set_major_formatter('dd:mm:ss')

    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    lat.set_minor_frequency(10)

    ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
    
    if Zoom:
        ax.set_xlim((Nx - boxsize)/2, (Nx + boxsize)/2)
        ax.set_ylim((Ny - boxsize)/2, (Ny + boxsize)/2)

    plt.savefig(save_directory + 'median_image_' + File[:-5] + '.jpg',bbox_inches='tight')
    plt.savefig(save_directory + 'median_image_' + File[:-5] + '.pdf',bbox_inches='tight')
    plt.close()
    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(median_image)
    HDR = HDU.header
    #Saving the correct header data. 
    i = 0
    while i < len(Header_saved_Names):
        HDR[Header_saved_Names[i]] = Header_saved_Data[i]
        i += 1
        #End while
    HDU.writeto(save_directory + 'median_image_' + File[:-5] + '.fits', clobber = True)
    #---------------------------------------------------------------------------
    #Plotting the max image.
    #---------------------------------------------------------------------------
    print 'Plotting the max image.'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = wcs)  
    cax = ax.imshow(max_image)
    ax.contour(max_image, Max_Contour_Levels, colors = 'white')
    ax.fill(x_beam, y_beam, alpha=1, facecolor='white', edgecolor='white', linewidth=1, zorder=2)

    ax.invert_yaxis()

    cbar = fig.colorbar(cax)
    cbar.set_label(Image_Units)

    lon = ax.coords['ra']
    lat = ax.coords['dec']

    ax.set_title('Max image of ' + Object + ' in ' + Line)
    lon.set_axislabel('J2000 Right Ascension')
    lat.set_axislabel('J2000 Declination')

    lon.set_major_formatter('hh:mm:ss')
    lat.set_major_formatter('dd:mm:ss')

    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    lat.set_minor_frequency(10)

    ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
    
    if Zoom:
        ax.set_xlim((Nx - boxsize)/2, (Nx + boxsize)/2)
        ax.set_ylim((Ny - boxsize)/2, (Ny + boxsize)/2)

    plt.savefig(save_directory + 'max_image_' + File[:-5] + '.jpg',bbox_inches='tight')
    plt.savefig(save_directory + 'max_image_' + File[:-5] + '.pdf',bbox_inches='tight')
    plt.close()
    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(max_image)
    HDR = HDU.header
    #Saving the correct header data. 
    i = 0
    while i < len(Header_saved_Names):
        HDR[Header_saved_Names[i]] = Header_saved_Data[i]
        i += 1
        #End while
    HDU.writeto(save_directory + 'max_image_' + File[:-5] + '.fits', clobber = True)
    #---------------------------------------------------------------------------
    #Plotting the min image.
    #---------------------------------------------------------------------------
    print 'Plotting the min image.'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = wcs)  
    cax = ax.imshow(min_image)
    ax.fill(x_beam, y_beam, alpha=1, facecolor='white', edgecolor='white', linewidth=1, zorder=2)

    ax.invert_yaxis()

    cbar = fig.colorbar(cax)
    cbar.set_label(Image_Units)

    lon = ax.coords['ra']
    lat = ax.coords['dec']

    ax.set_title('Min image of ' + Object + ' in ' + Line)
    lon.set_axislabel('J2000 Right Ascension')
    lat.set_axislabel('J2000 Declination')

    lon.set_major_formatter('hh:mm:ss')
    lat.set_major_formatter('dd:mm:ss')

    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    lat.set_minor_frequency(10)

    ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
    
    if Zoom:
        ax.set_xlim((Nx - boxsize)/2, (Nx + boxsize)/2)
        ax.set_ylim((Ny - boxsize)/2, (Ny + boxsize)/2)

    plt.savefig(save_directory + 'min_image_' + File[:-5] + '.jpg',bbox_inches='tight')
    plt.savefig(save_directory + 'min_image_' + File[:-5] + '.pdf',bbox_inches='tight')
    plt.close()
    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(min_image)
    HDR = HDU.header
    #Saving the correct header data. 
    i = 0
    while i < len(Header_saved_Names):
        HDR[Header_saved_Names[i]] = Header_saved_Data[i]
        i += 1
        #End while
    HDU.writeto(save_directory + 'min_image_' + File[:-5] + '.fits', clobber = True)
    #---------------------------------------------------------------------------
    #Plotting the standard deviation image.
    #---------------------------------------------------------------------------
    print 'Plotting the standard deviation image.'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = wcs)  
    cax = ax.imshow(standard_deviation_image)
    ax.fill(x_beam, y_beam, alpha=1, facecolor='white', edgecolor='white', linewidth=1, zorder=2)

    ax.invert_yaxis()

    cbar = fig.colorbar(cax)
    cbar.set_label(Image_Units)

    lon = ax.coords['ra']
    lat = ax.coords['dec']

    ax.set_title('Standard deviation image of ' + Object + ' in ' + Line)
    lon.set_axislabel('J2000 Right Ascension')
    lat.set_axislabel('J2000 Declination')

    lon.set_major_formatter('hh:mm:ss')
    lat.set_major_formatter('dd:mm:ss')

    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    lat.set_minor_frequency(10)

    ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
    
    if Zoom:
        ax.set_xlim((Nx - boxsize)/2, (Nx + boxsize)/2)
        ax.set_ylim((Ny - boxsize)/2, (Ny + boxsize)/2)

    plt.savefig(save_directory + 'standard_deviation_image_' + File[:-5] + '.jpg',bbox_inches='tight')
    plt.savefig(save_directory + 'standard_deviation_image_' + File[:-5] + '.pdf',bbox_inches='tight')
    plt.close()
    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(standard_deviation_image)
    HDR = HDU.header
    #Saving the correct header data. 
    i = 0
    while i < len(Header_saved_Names):
        HDR[Header_saved_Names[i]] = Header_saved_Data[i]
        i += 1
        #End while
    HDU.writeto(save_directory + 'standard_deviation_image_' + File[:-5] + '.fits', clobber = True)
    #---------------------------------------------------------------------------
    #Plotting the root mean square image.
    #---------------------------------------------------------------------------
    print 'Plotting the root mean square image.'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = wcs)  
    cax = ax.imshow(root_mean_square_image)
    ax.fill(x_beam, y_beam, alpha=1, facecolor='white', edgecolor='white', linewidth=1, zorder=2)

    ax.invert_yaxis()

    cbar = fig.colorbar(cax)
    cbar.set_label(Image_Units)

    lon = ax.coords['ra']
    lat = ax.coords['dec']

    ax.set_title('Root mean square image of ' + Object + ' in ' + Line)
    lon.set_axislabel('J2000 Right Ascension')
    lat.set_axislabel('J2000 Declination')

    lon.set_major_formatter('hh:mm:ss')
    lat.set_major_formatter('dd:mm:ss')

    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    lat.set_minor_frequency(10)

    ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
    
    if Zoom:
        ax.set_xlim((Nx - boxsize)/2, (Nx + boxsize)/2)
        ax.set_ylim((Ny - boxsize)/2, (Ny + boxsize)/2)

    plt.savefig(save_directory + 'root_mean_square_image_' + File[:-5] + '.jpg',bbox_inches='tight')
    plt.savefig(save_directory + 'root_mean_square_image_' + File[:-5] + '.pdf',bbox_inches='tight')
    plt.close()
    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(root_mean_square_image)
    HDR = HDU.header
    #Saving the correct header data. 
    i = 0
    while i < len(Header_saved_Names):
        HDR[Header_saved_Names[i]] = Header_saved_Data[i]
        i += 1
        #End while
    HDU.writeto(save_directory + 'root_mean_square_image_' + File[:-5] + '.fits', clobber = True)
    #---------------------------------------------------------------------------
    #Plotting the absolute mean deviation image.
    #---------------------------------------------------------------------------
    print 'Plotting the mean deviation image.'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = wcs)  
    cax = ax.imshow(absolute_mean_deviation_image)
    ax.fill(x_beam, y_beam, alpha=1, facecolor='white', edgecolor='white', linewidth=1, zorder=2)

    ax.invert_yaxis()

    cbar = fig.colorbar(cax)
    cbar.set_label(Image_Units)

    lon = ax.coords['ra']
    lat = ax.coords['dec']

    ax.set_title('Absolute mean deviation image of ' + Object + ' in ' + Line)
    lon.set_axislabel('J2000 Right Ascension')
    lat.set_axislabel('J2000 Declination')

    lon.set_major_formatter('hh:mm:ss')
    lat.set_major_formatter('dd:mm:ss')

    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    lat.set_minor_frequency(10)

    ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
    
    if Zoom:
        ax.set_xlim((Nx - boxsize)/2, (Nx + boxsize)/2)
        ax.set_ylim((Ny - boxsize)/2, (Ny + boxsize)/2)

    plt.savefig(save_directory + 'absolute_mean_deviation_image_' + File[:-5] + '.jpg',bbox_inches='tight')
    plt.savefig(save_directory + 'absolute_mean_deviation_image_' + File[:-5] + '.pdf',bbox_inches='tight')
    plt.close()
    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(absolute_mean_deviation_image)
    HDR = HDU.header
    #Saving the correct header data. 
    i = 0
    while i < len(Header_saved_Names):
        HDR[Header_saved_Names[i]] = Header_saved_Data[i]
        i += 1
        #End while
    HDU.writeto(save_directory + 'absolute_mean_deviation_image_' + File[:-5] + '.fits', clobber = True)
    #---------------------------------------------------------------------------
    #Plotting the sum image.
    #---------------------------------------------------------------------------
    print 'Plotting the sum image.'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = wcs)  
    cax = ax.imshow(sum_image)
    ax.contour(sum_image, Sum_Contour_Levels, colors = 'white')
    ax.fill(x_beam, y_beam, alpha=1, facecolor='white', edgecolor='white', linewidth=1, zorder=2)

    ax.invert_yaxis()

    cbar = fig.colorbar(cax)
    cbar.set_label(Image_Units)

    lon = ax.coords['ra']
    lat = ax.coords['dec']

    ax.set_title('Sum image of ' + Object + ' in ' + Line)
    lon.set_axislabel('J2000 Right Ascension')
    lat.set_axislabel('J2000 Declination')

    lon.set_major_formatter('hh:mm:ss')
    lat.set_major_formatter('dd:mm:ss')

    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    lat.set_minor_frequency(10)

    ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
    
    if Zoom:
        ax.set_xlim((Nx - boxsize)/2, (Nx + boxsize)/2)
        ax.set_ylim((Ny - boxsize)/2, (Ny + boxsize)/2)

    plt.savefig(save_directory + 'sum_image_' + File[:-5] + '.jpg',bbox_inches='tight')
    plt.savefig(save_directory + 'sum_image_' + File[:-5] + '.pdf',bbox_inches='tight')
    plt.close() 
    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(sum_image)
    HDR = HDU.header
    #Saving the correct header data. 
    i = 0
    while i < len(Header_saved_Names):
        HDR[Header_saved_Names[i]] = Header_saved_Data[i]
        i += 1
        #End while
    HDU.writeto(save_directory + 'sum_image_' + File[:-5] + '.fits', clobber = True)
    '''
    #---------------------------------------------------------------------------
    #Plotting the moment 0 map.
    #---------------------------------------------------------------------------
    print 'Plotting the moment 0 map.'
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = wcs)  
    cax = ax.imshow(integrated_image, cmap = 'magma')
    
    Cropped_Data_Contour, Contour_Levels_Contour = Create_Continuum_Contour(header)
    ax.contour(Cropped_Data_Contour, Contour_Levels_Contour, colors = 'white')
    #ax.contour(integrated_image, Integrated_Contour_Levels, colors = 'white')
    ax.fill(x_beam, y_beam, alpha=1, facecolor='white', edgecolor='white', linewidth=1, zorder=2)

    ax.invert_yaxis()

    cbar = fig.colorbar(cax)
    cbar.set_label('Jy km s$^{-1}$ beam$^{-1}$')
    
    lon = ax.coords['ra']
    lat = ax.coords['dec']

    #ax.set_title('Moment 0 map of ' + Object + ' in ' + Line)
    lon.set_axislabel('J2000 Right Ascension')
    lat.set_axislabel('J2000 Declination')

    #lon.set_major_formatter('hh:mm:ss')
    #lat.set_major_formatter('dd:mm:ss')
    
    lon.set_ticks(color='white', exclude_overlapping=True)
    lat.set_ticks(color='white', exclude_overlapping=True)

    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    lat.set_minor_frequency(10)
    
    #ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
    
    if Zoom:
        ax.set_xlim((Nx - boxsize)/2, (Nx + boxsize)/2)
        ax.set_ylim((Ny - boxsize)/2, (Ny + boxsize)/2)

    plt.savefig(save_directory + 'Moment_0_map_' + File[:-5] + '.jpg',bbox_inches='tight')
    plt.savefig(save_directory + 'Moment_0_map_' + File[:-5] + '.pdf',bbox_inches='tight')
    plt.close()
    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(integrated_image)
    HDR = HDU.header
    #Saving the correct header data. 
    i = 0
    while i < len(Header_saved_Names):
        HDR[Header_saved_Names[i]] = Header_saved_Data[i]
        i += 1
        #End while
    HDR['BUNIT'] = 'km/s'
    HDU.writeto(save_directory + 'Moment_0_map_' + File[:-5] + '.fits', clobber = True)
    #---------------------------------------------------------------------------
    #Plotting the moment 1 map.
    #---------------------------------------------------------------------------
    print 'Plotting the moment 1 map.'
    #First make the values that are not relevant white
    #velocity_field_image[velocity_field_image == 0] = np.nan

    #Then we plot.
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = wcs)  
    cmap = mpl.cm.bwr
    cmap.set_bad('k',1.)
    cax = ax.imshow(velocity_field_image, vmin = vmin_M1, vmax = vmax_M1, cmap = cmap)
    
    ax.contour(Cropped_Data_Contour, Contour_Levels_Contour, colors = 'black')
    #ax.contour(integrated_image, Integrated_Contour_Levels, colors = 'black')
    ax.fill(x_beam, y_beam, alpha=1, facecolor='white', edgecolor='white', linewidth=1, zorder=2)

    ax.invert_yaxis()

    cbar = fig.colorbar(cax)
    cbar.set_label('km s$^{-1}$')

    lon = ax.coords['ra']
    lat = ax.coords['dec']

    #ax.set_title('Moment 1 map of ' + Object + ' in ' + Line)
    lon.set_axislabel('J2000 Right Ascension')
    lat.set_axislabel('J2000 Declination')

    #lon.set_major_formatter('hh:mm:ss')
    #lat.set_major_formatter('dd:mm:ss')
    
    lon.set_ticks(color='white', exclude_overlapping=True)
    lat.set_ticks(color='white', exclude_overlapping=True)

    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    lat.set_minor_frequency(10)

    #ax.coords.grid(color='black', alpha=0.5, linestyle='solid')
    
    if Zoom:
        ax.set_xlim((Nx - boxsize)/2, (Nx + boxsize)/2)
        ax.set_ylim((Ny - boxsize)/2, (Ny + boxsize)/2)

    plt.savefig(save_directory + 'Moment_1_map_' + File[:-5] + '.jpg',bbox_inches='tight')
    plt.savefig(save_directory + 'Moment_1_map_' + File[:-5] + '.pdf',bbox_inches='tight')
    plt.close()
    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(velocity_field_image)
    HDR = HDU.header
    #Saving the correct header data. 
    i = 0
    while i < len(Header_saved_Names):
        HDR[Header_saved_Names[i]] = Header_saved_Data[i]
        i += 1
        #End while
    HDR['BUNIT'] = 'km/s'
    HDU.writeto(save_directory + 'Moment_1_map_' + File[:-5] + '.fits', clobber = True)
    #---------------------------------------------------------------------------
    #Plotting the moment 2 map.
    #--------------------------------------------------------------------------- 
    print 'Plotting the moment 2 map.'
    #First make the values that are not relevant white
    velocity_dispersion_field_image[velocity_dispersion_field_image == 0] = np.nan

    #Then plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = wcs)  
    cmap = mpl.cm.viridis
    cmap.set_bad('k',1.)
    cax = ax.imshow(velocity_dispersion_field_image[:,:], vmin = vmin_M2, vmax = vmax_M2, cmap = cmap)
    ax.fill(x_beam, y_beam, alpha=1, facecolor='white', edgecolor='white', linewidth=1, zorder=2)
    
    ax.contour(Cropped_Data_Contour, Contour_Levels_Contour, colors = 'white', linewidth = 0.25)
    
    #ax.contour(integrated_image, Integrated_Contour_Levels, colors = 'grey', linewidth = 0.25)

    ax.invert_yaxis()

    cbar = fig.colorbar(cax)
    cbar.set_label('km s$^{-1}$')

    lon = ax.coords['ra']
    lat = ax.coords['dec']

    #ax.set_title('Moment 2 map of ' + Object + ' in ' + Line)
    lon.set_axislabel('J2000 Right Ascension')
    lat.set_axislabel('J2000 Declination')
    
    #lon.set_major_formatter('hh:mm:ss')
    #lat.set_major_formatter('dd:mm:ss')
    
    lon.set_ticks(color='white', exclude_overlapping=True)
    lat.set_ticks(color='white', exclude_overlapping=True)

    lon.display_minor_ticks(True)
    lat.display_minor_ticks(True)
    lat.set_minor_frequency(10)

    #ax.coords.grid(color='white', alpha=0.5, linestyle='solid')
    
    if Zoom:
        ax.set_xlim((Nx - boxsize)/2, (Nx + boxsize)/2)
        ax.set_ylim((Ny - boxsize)/2, (Ny + boxsize)/2)

    plt.savefig(save_directory + 'Moment_2_map_' + File[:-5] + '.jpg',bbox_inches='tight')
    plt.savefig(save_directory + 'Moment_2_map_' + File[:-5] + '.pdf',bbox_inches='tight')
    plt.close()
    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(velocity_dispersion_field_image)
    HDR = HDU.header
    #Saving the correct header data. 
    i = 0
    while i < len(Header_saved_Names):
        HDR[Header_saved_Names[i]] = Header_saved_Data[i]
        i += 1
        #End while
    HDR['BUNIT'] = 'km/s'
    HDU.writeto(save_directory + 'Moment_2_map_' + File[:-5] + '.fits', clobber = True)   
    #---------------------------------------------------------------------------
    #Saving data in .txt file
    #---------------------------------------------------------------------------
    np.savetxt(save_directory + 'mean_image_' + File[:-5] + '.txt', mean_image, delimiter = ' ')

    np.savetxt(save_directory + 'median_image_' + File[:-5] + '.txt', median_image, delimiter = ' ')

    np.savetxt(save_directory + 'max_image_' + File[:-5] + '.txt', max_image, delimiter = ' ')

    np.savetxt(save_directory + 'min_image_' + File[:-5] + '.txt', min_image, delimiter = ' ')

    np.savetxt(save_directory + 'standard_deviation_image_' + File[:-5] + '.txt', standard_deviation_image, delimiter = ' ')

    np.savetxt(save_directory + 'root_mean_square_image_' + File[:-5] + '.txt', root_mean_square_image, delimiter = ' ')

    np.savetxt(save_directory + 'absolute_mean_deviation_image_' + File[:-5] + '.txt', absolute_mean_deviation_image, delimiter = ' ')

    np.savetxt(save_directory + 'sum_image_' + File[:-5] + '.txt', sum_image, delimiter = ' ')

    np.savetxt(save_directory + 'moment_0_map_' + File[:-5] + '.txt', integrated_image, delimiter = ' ')

    np.savetxt(save_directory + 'moment_1_map_' + File[:-5] + '.txt', velocity_field_image, delimiter = ' ')

    np.savetxt(save_directory + 'moment_2_map_' + File[:-5] + '.txt', velocity_dispersion_field_image, delimiter = ' ')
    
if __name__ == '__main__':
	#C180
	begin_image = 18
	end_image = 54

	Object = 'L1527'
	Line = 'C$^{18}$O'

	File = 'L1527_C18O_1.fits'
	Directory = 'C18O_1/'
	#Image_Moments(File, Directory, begin_image, end_image, Object, Line, Zoom = True, boxsize = 50)
	
	File = 'L1527_C18O_2.fits'
	Directory = 'C18O_2/'
	Image_Moments(File, Directory, begin_image, end_image, Object, Line, Zoom = True, boxsize = 160, vmin_M1 = 4, vmax_M1 = 7.8, vmin_M2 = 0, vmax_M2 = 1.8)

