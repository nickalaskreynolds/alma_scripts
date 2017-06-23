import numpy as np
import pyfits as pf
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.ndimage import interpolation as inter
from astropy import units as u
from astropy import constants as const
from astropy.coordinates import SkyCoord
from astropy.io import fits

#Make the plots nice
import Nice_Plots_2
Nice_Plots_2.set_style()

def Rotate_and_Sum(Image, Rotation, Width):
    '''
    This function rotates an given image with an angle Rotation and then
    sums the pixels in y-direction in the range 
    Ny/2-width/2:Ny/2+width/2, with Ny the length of the array in the y
    direction. It then returns an 1D array. 
    Image = 2D array 
    Rotation = Is defined positive from the x-axis to the y-axis 
    (anti-clockwise), must be given in degrees.
    Width = Total width over wich there will be summed, must be given in
    pixels.
    ''' 
    Rotated_Image = inter.rotate(Image, Rotation)

    Shape_Rotated_Image = np.shape(Rotated_Image)

    Summed_Array = np.sum(Rotated_Image[(np.round(Shape_Rotated_Image[1]/2) - np.round(Width/2)):(np.round(Shape_Rotated_Image[1]/2) + np.round(Width/2))],0)
    
    temp = np.zeros([1,Shape_Rotated_Image[1]])
    temp[0,:] = Summed_Array
    
    return temp
         
def PV_diagram(File, Directory, rotation_deg, width, Object, Line, inclination, Velocity_Curve = False, mass = 0, mass_err = 0, v_source = 0, d_source = 0, Thindisk = False, Zoom = False, v_width = 5, arcsec_width = 20, Overlay_Contour = 'None'):
    '''
    This function creates an position velocity diagram for a given fits 
    file.
    The rotation that has to be performed must be given in 
    degrees. The width is given in pixels and is assumed to be around 
    the middle. If Velocity_Curve is True then a velocity curve based on
    the given mass and mass error is produced. Both should be given in 
    solar masses.
    
    File = Name of the 3D datacube. [string]
    Directory = Directory where the file is saved.  [string]
    rotation_deg = The rotation in degrees that has to be performed to 
    align the axis over which there will be summed with the y-axis. 
    [float]
    width = The total amount of pixels over which will be summed. 
    [float]
    Object = The name of the object, only used in title figure. [string]
    Line = The molecular line of the data, only used in title figure. 
    [string] 
    Velocity_Curve = If true, a keplerian rotation curve will be drawn
    in the figure based on the parameters explained below. [boolean]
    Mass = The mass of the central object in solar masses, used for the 
    velocity curve. [float]
    Mass_Err = The error on the mass in solar masses, used for the 
    velocity curve. [float]
    v_source = The source velocity in km s^-1, used for the velocity 
    curve and the zoom function.[float]
    d_source = The distance of the source in parsec, used for the 
    velocity curve. [float]
    Thindisk = Set to true if the data is from a thindisk model, so it 
    wont draw any contour lines based on the sigma's in the image.
    [boolean]
    Zoom = If true the function will zoom in. [boolean]
    v_width = The total range in velocities (km s^-1) that will be shown
    when zoomed in. [float]
    arcsec_width = The total range in arcsec that will be shown when 
    zoomed in. [float]
    Overlay_Contour = A contour will be drawn over the figure of a 
    different dataset, which already must have a PV-diagram made by  
    this programm. If 'None', no contour will be drawn. If you want to 
    add these contours, you have to line 191 to manually add the 
    directories and filenames of the data. [string] 
    '''
    print 'Creating the PV-diagram for ' + Object
    print ' ' 
    #---------------------------------------------------------------------------
    #Getting the data from the fits file.
    #---------------------------------------------------------------------------
    
    #getting the header.
    header = pf.getheader(Directory + File)    
    
    #Getting the data.
    data = pf.getdata(Directory + File)

    #The datacubes coming from the thindisk model have 3 dimensions, 
    #while the science datacubes have 4 dimension. So we have to account
    #for that. 
    if len(np.shape(data)) == 4:
        Data = data[0,:,:,:]
    else:
        Data = data

    #Determining the shape of the data.
    Shape_Data = np.shape(Data)

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
        
    else:
        print 'I am not sure how to get the velocities.'
        raise KeyError
    
    PixelWidth_RA = header['CDELT1']
    PixelWidth_DEC = header['CDELT2']

    #The used units.
    Image_Units = header['BUNIT']

    #Creating the directory for saving the images, if it does not 
    #already exists.
    save_directory = Directory + 'PV-diagram/' 
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
    
    #---------------------------------------------------------------------------
    #Creating the array that contains the PV-diagram
    #---------------------------------------------------------------------------
    rotation_rad = np.radians(rotation_deg)
    #Creating the array of the correct size.
    y_size = np.round(np.abs(np.cos(np.abs(rotation_rad))*Shape_Data[1]) + np.abs(np.cos(np.pi/2 - np.abs(rotation_rad))*Shape_Data[2]))
    PV_Data = np.zeros([y_size, Shape_Data[0]])    

    i = 0
    while i < Shape_Data[0]:
        PV_Data[:,i] = Rotate_and_Sum(Data[i,:,:], rotation_deg, width)
        i += 1
        #End while loop
    #---------------------------------------------------------------------------
    #Contour lines
    #---------------------------------------------------------------------------
    std_PV = np.std(PV_Data[0:100,0:100]) 
    PV_Contour_Levels = np.arange(0, 50, 10)*std_PV
    PV_Contour_Levels[0] = 3*std_PV

    #---------------------------------------------------------------------------
    #Plotting
    #---------------------------------------------------------------------------
    #First we create arrays for the correct ticks
    x_values = np.arange(begin_v, begin_v + delta_v*Shape_Data[0], delta_v)
    
    #The total length in arcsec of the y axis in the new image.
    length_arcsec_new =  (np.abs(np.cos(np.abs(rotation_rad)))*np.abs(PixelWidth_RA)*Shape_Data[2]+np.abs(np.cos(np.pi/2-np.abs(rotation_rad)))*np.abs(PixelWidth_DEC)*Shape_Data[1])*3600
    y_values = np.arange(-length_arcsec_new/2, length_arcsec_new/2 + length_arcsec_new/10, length_arcsec_new/10)

    #Calculating the size of y pixel in the y direction in arcsec.
    pixel_size_y_arcsec = length_arcsec_new/y_size

    y_arcsec = np.arange(1, PV_Data.shape[0])*pixel_size_y_arcsec - length_arcsec_new/2

    #Start plotting
    print 'Starting the plotting.'
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.imshow(PV_Data, origin = 'lower', cmap = 'magma')

    #The thindisk model does not contain noise, so no contourlines are needed.
    if not Thindisk:
        ax.contour(PV_Data, PV_Contour_Levels, colors = 'white')
    
    def Make_Contour(Directory, File, Color, Line, Thindisk = False):
        '''
        This file opens the file given and plots the contours of this 
        file.
        '''
        Data = pf.getdata(Directory + File)
        
        if Thindisk:
            #If we have thindisk model data, we have no noise. So we just
            #use percentages instead.
            Max = np.max(Data)
            Data_Contour_Levels = np.arange(0.2, 1, 0.2)*Max
        else:
            std_Data = np.std(Data[0:100,0:100]) 
            Data_Contour_Levels = np.arange(0, 30, 5)*std_Data
            Data_Contour_Levels[0] = 3*std_Data
            
        ax.contour(Data, Data_Contour_Levels, colors = Color)
        ax.axvline(-10000, color = Color, label = Line)
        #ax.legend(loc = 3)
        
    if Overlay_Contour == '13CO':
        Directory = '13CO_1/PV-diagram/'
        File_Name = 'PV-Diagram_L1165-13CO_1.fits'
        Make_Contour(Directory, File_Name, 'yellow', '$^{13}$CO') 
        ax.legend(loc = 3)
    if Overlay_Contour == 'C18O':
        Directory = 'C18O_2/PV-diagram/'
        File_Name = 'PV-Diagram_L1527_C18O_2.fits'
        Make_Contour(Directory, File_Name, 'white', 'C$^{18}$O') 
        #ax.legend(loc = 3)
    if Overlay_Contour == 'Model':
        #Directory = 'Optimization_Results/Round_1/Model_Most_Common_2_52/PV-diagram/'
        #File_Name = 'PV-Diagram_Optimization_Most_Common_2_52_convolved_with_beam.fits'
        Directory =  'Optimization_Results/Round_3_d_414/Model_Best_38.3/PV-diagram/'
        File_Name = 'PV-Diagram_Optimization_Best_38.3_2_convolved_with_beam.fits'
        Make_Contour(Directory, File_Name, 'white', 'Model', Thindisk = True) 
        #ax.legend(loc = 3)
    if Overlay_Contour == 'Both':
        Directory = 'C18O_1/PV-diagram/'
        File_Name = 'PV-Diagram_L1165-C18O_1.fits'
        Make_Contour(Directory, File_Name, 'yellow', 'C$^{18}$O')
        
        Directory = '13CO_1/PV-diagram/'
        File_Name = 'PV-Diagram_L1165-13CO_1.fits'
        Make_Contour(Directory, File_Name, 'red', '$^{13}$CO')

    #ax.set_title('PV diagram of '+ Object + ' in ' + Line)
    ax.set_xlabel('Velocity [km s$^{-1}$]')
    ax.set_ylabel('Offset [$"$]')#'Position [arcsec]')$\Delta$X

    offset = 10
    #If we have correct masses we can calculate the velocity curve.
    if Velocity_Curve:
        print 'Including the velocities curves.'
        #Calculating the extreme masses within the errors, do we can also plot
        #those. 
        mass_min_err = mass - mass_err
        mass_plus_err = mass + mass_err
        
        #Creating an array containing the velocities in km s^-1. 
        velocities = (np.arange(begin_v, begin_v + delta_v*Shape_Data[0], delta_v) - v_source)

        #This function returns for a given mass (solar masses), velocity (in
        #km s^-1) and distance to the source (in pc) the radius (in arcsec)  
        #assuming Keplerian rotation.
        def Keplerian_Rotation(mass, velocity, Distance, inclination):
            radii_return =  np.sin(inclination)**2*const.G.value*mass*const.M_sun.value/(velocity*1000)/(velocity*1000)/(Distance*u.pc.to(u.m))*u.rad.to(u.arcsec) 
            #All the positive radii.
            radii_positive = radii_return[velocity < 0]
            #We also have some negative radii, so thats why we have to do this.
            radii_negative = -1*radii_return[velocity > 0]        
            return radii_positive, radii_negative
        
        #Calculate the radii.
        radii_positive, radii_negative = Keplerian_Rotation(mass, velocities, d_source, inclination)
        radii_positive_min_err, radii_negative_min_err = Keplerian_Rotation(mass_min_err, velocities, d_source, inclination)
        radii_positive_plus_err, radii_negative_plus_err = Keplerian_Rotation(mass_plus_err, velocities, d_source, inclination)

        #Changing the radii to the correct pixel coordinates for correct 
        #plotting. Plus bring the lines to the object. 
        radii_positive_pixel_coor = radii_positive/pixel_size_y_arcsec + (y_size - 1)/2  + offset
        radii_negative_pixel_coor = radii_negative/pixel_size_y_arcsec + (y_size - 1)/2  + offset

        radii_positive_min_err_pixel_coor = radii_positive_min_err/pixel_size_y_arcsec + (y_size - 1)/2 + offset 
        radii_negative_min_err_pixel_coor = radii_negative_min_err/pixel_size_y_arcsec + (y_size - 1)/2 + offset

        radii_positive_plus_err_pixel_coor = radii_positive_plus_err/pixel_size_y_arcsec + (y_size - 1)/2 + offset  
        radii_negative_plus_err_pixel_coor = radii_negative_plus_err/pixel_size_y_arcsec + (y_size - 1)/2 + offset
        
        #Plotting the velocities
        ax.plot(np.arange(0,len(radii_positive), 1), radii_positive_pixel_coor, color = 'white', label = 'Keplerian rotation', linestyle = ':')
        ax.plot(np.arange(len(radii_positive) , len(velocities), 1), radii_negative_pixel_coor, color = 'white', linestyle = ':')

        ax.plot(np.arange(0,len(radii_positive), 1), radii_positive_min_err_pixel_coor, color = 'white', linestyle = ':')
        ax.plot(np.arange(len(radii_positive) , len(velocities), 1), radii_negative_min_err_pixel_coor, color = 'white', linestyle = ':')

        ax.plot(np.arange(0,len(radii_positive), 1), radii_positive_plus_err_pixel_coor, color = 'white', linestyle = ':')
        ax.plot(np.arange(len(radii_positive) , len(velocities), 1), radii_negative_plus_err_pixel_coor, color = 'white', linestyle = ':')

        ax.axhline(np.where(y_arcsec > 0)[0][0] - 1 + offset, color = 'white', linestyle = '--')
        ax.axvline(np.where(velocities > 0)[0][0] - 0.33 , color = 'white', linestyle = '--', label = 'v$_{source}$ = ' + str(v_source) + ' km s$^{-1}$')
        #ax.legend(loc = 3)

    #If zoom is true we zoom in else we just show the entire image. 
    if Zoom:
        #First we determine at which pixel we have the source velocity.
        pix_v_source = int(np.abs((begin_v - v_source)/delta_v))
        #Then we determine what half the width of the v slice must be.
        pix_v_shift = int(v_width/delta_v/2)

        #Now we determine the central pixel for the arcsec.
        pix_arcsec_central = int(y_size/2) - 1 + offset
        pix_arcsec_shift = int(arcsec_width/pixel_size_y_arcsec/2)
        
        x = np.arange(pix_v_source - pix_v_shift, pix_v_source + pix_v_shift + 2*pix_v_shift/10., 2*pix_v_shift/10.)
        #y = np.arange(pix_arcsec_central - pix_arcsec_shift, pix_arcsec_central + pix_arcsec_shift + 2*pix_arcsec_shift/10., 2*pix_arcsec_shift/10)
        y = np.linspace(pix_arcsec_central - pix_arcsec_shift,pix_arcsec_central + pix_arcsec_shift, 11)
        
        x_label = np.round(begin_v + delta_v*x, 1)
        #y_label = np.round((y+1)*pixel_size_y_arcsec - length_arcsec_new/2., 2) + 0.01
        y_label = np.linspace(-arcsec_width/2, arcsec_width/2, 11)

        #Doing the zooming by limiting the shown x and y coordinates.
        plt.xlim([pix_v_source - pix_v_shift, pix_v_source + pix_v_shift])
        plt.ylim([pix_arcsec_central - pix_arcsec_shift, pix_arcsec_central + pix_arcsec_shift]) 

        ax.set_aspect(1.*pix_v_shift/pix_arcsec_shift)

        from matplotlib.ticker import FormatStrFormatter
        #ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

        y_label = np.array(["%.2f" % i for i in y_label])
        ax.set_xticks(x)
        ax.set_yticks(y)
        ax.set_xticklabels(x_label)
        ax.set_yticklabels(y_label)

    else:
        x = np.arange(0,Shape_Data[0] + Shape_Data[0]/10, Shape_Data[0]/10)
        y = np.arange(0,y_size + np.round(y_size/10), np.round(y_size/10))
        x_label = np.round(x_values[0::Shape_Data[0]/10], 1)
        y_label = np.round(y_values,1)
        xticks=np.linspace(-4,1,11)
        ax.set_xticks(xticks)
        ax.set_yticks(y)
        #ax.set_xticklabels(x_label)
        ax.set_yticklabels(y_label)
    
        ax.set_xlim([0, Shape_Data[0]-1])
        ax.set_ylim([0,y_size-1])
        ax.set_aspect(Shape_Data[0]/y_size)
        #ax.set_aspect((delta_v*(Shape_Data[0] - 1) )/length_arcsec_new)


    #ax.grid(color='white', alpha=0.5, linestyle='solid')
    cbar = plt.colorbar(cax)
    cbar.set_label('Jy/Beam')
    
    #ax.set_xticks(color='white', exclude_overlapping=True)
    #ax.set_yticks(color='white', exclude_overlapping=True)

    plt.savefig(save_directory + 'PV-Diagram_' + File[:-5] + '.pdf',bbox_inches='tight')
    plt.savefig(save_directory + 'PV-Diagram_' + File[:-5] + '.jpg',bbox_inches='tight')
    plt.savefig(save_directory + 'PV-Diagram_' + File[:-5] + '.ps',bbox_inches='tight')
    plt.savefig(save_directory + 'PV-Diagram_' + File[:-5] + '.eps',bbox_inches='tight')
    plt.close()
    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(PV_Data)
    HDR = HDU.header
    #Saving the correct header data. 

    #Need to create header.
    HDU.writeto(save_directory + 'PV-Diagram_' + File[:-5] + '.fits', clobber = True)

    #---------------------------------------------------------------------------
    #Saving data in .txt file
    #---------------------------------------------------------------------------    
    np.savetxt(save_directory + 'PV-Diagram_' + File[:-5] + '.txt', PV_Data, delimiter = ' ')

if __name__ == '__main__':
    rotation = -90. #degrees

    Mass = 0.3#Solar masses 
    Mass_err = 0#Solar masses
    v_source = 5.9#5#km s^-1
    distance_source = 140.#pc
    inclination = 89.
         
    #-------------------------------------------------------------------------------
    #C180
    #-------------------------------------------------------------------------------
    Object = 'L1527'
    Line = 'C$^{18}$O'

    width = 50#pixels	
    File = 'L1527_C18O_2.fits'
    Directory = 'C18O_2/'


    PV_diagram(File, Directory, rotation, width, Object, Line, inclination, Velocity_Curve = True, mass = Mass, mass_err = Mass_err, v_source = v_source, d_source = distance_source, Thindisk = False, Zoom = True, v_width = 8, arcsec_width = 10,Overlay_Contour = 'None')
    
