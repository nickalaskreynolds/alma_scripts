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
from astropy.wcs import WCS
from time import time
_TIME_=time()

#Make the plots nice
import Nice_Plots_2
Nice_Plots_2.set_style()

def Rotate_and_Sum(Filename,Image, Rotation, Width, ra, dec):
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
    #print(Filename)
    #print(Image.shape)
    imcen=Image.shape[0]/2

    w = WCS(Filename)
    #print(ra,dec)
    xcen, ycen = w.sub(['longitude','latitude']).wcs_world2pix(ra, dec,0)
    #print(xcen, ycen)
    center_shift=[imcen-ycen,imcen-xcen]
    #center_shift=[ycen-imcen,xcen-imcen]
    #print(center_shift)
    ''' 
    if center_shift[0] >0:
       center_shift[0]=center_shift[0]*(-1.0)
    if center_shift[1] >0:
       center_shift[1]=center_shift[1]*(-1.0)
    print(center_shift)
    '''
    Shifted_Image=inter.shift(Image,center_shift)
    #Shifted_image=Image[xcen-cropsize:xcen+cropsize,ycen-cropsize:ycen+cropsize]
    Rotated_Image = inter.rotate(Shifted_Image, Rotation)

    Shape_Rotated_Image = np.shape(Rotated_Image)

    Summed_Array = np.sum(Rotated_Image[(np.round(Shape_Rotated_Image[1]/2) - np.round(Width/2)):(np.round(Shape_Rotated_Image[1]/2) + np.round(Width/2))],0)
    
    temp = np.zeros([1,Shape_Rotated_Image[1]])
    temp[0,:] = Summed_Array
    #print(temp[0,:-1].shape)
    
    #return temp[0,:-1]
    return temp
        
def PV_diagram(ra,dec, File, Directory, rotation_deg, width, Object, Line, inclination, Velocity_Curve = False, mass = 0, mass_err = 0, v_source = 0, d_source = 0, Thindisk = False, Zoom = False, v_width = 5, arcsec_width = 20, Overlay_Contour = 'None',imagemin=-99.0,imagemax=-99.0,contour_interval=10.0):
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
    print('Creating the PV-diagram for ',Object)
    print(' ')
    _TIME_=time()

    font = {'family' : 'sans-serif',
        'weight' : 'bold',
        'size'   : 16}

    mpl.rc('font', **font) 
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
    print('Shape:',Shape_Data)

    #Determining the begin_v, delta_v and the amount of images N.
    print('Current 1: ',(time()-_TIME_))
    N = header['NAXIS3']
    if header['CTYPE3']  == 'VELO-LSR':
        print('True')
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
        delta_v =begin_v - begin_v_plus_one
    else:
        print('I am not sure how to get the velocities.')
        raise KeyError
    
    PixelWidth_RA = header['CDELT1']
    PixelWidth_DEC = header['CDELT2']

    #The used units.
    print('Rest freq: ',rest_freq)
    print('Vsys:',begin_v,'...Delta_v:',delta_v)
    Image_Units = header['BUNIT']

    #Creating the directory for saving the images, if it does not 
    #already exists.
    save_directory = Directory + 'PV-diagram-new/' 
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
    
    #---------------------------------------------------------------------------
    #Creating the array that contains the PV-diagram
    #---------------------------------------------------------------------------
    # NEED TO OPTIMIZE
    rotation_rad = np.radians(rotation_deg)
    #Creating the array of the correct size.
    y_size = np.round(np.abs(np.cos(np.abs(rotation_rad))*Shape_Data[1]) + np.abs(np.cos(np.pi/2 - np.abs(rotation_rad))*Shape_Data[2]))

    PV_Data = np.zeros([int(y_size), int(Shape_Data[0])])    

    i = 0
    #print('Shapes:',Shape_Data,PV_Data.shape,y_size)
    while i < Shape_Data[0]:
        #print(Directory+File)
        PV_Data[:,i] = Rotate_and_Sum(Directory+File,Data[i,:,:], rotation_deg, width,ra, dec)
        i += 1
        #End while loop
    print('Current 2: ',(time()-_TIME_))
    _TIME_=time()
    #---------------------------------------------------------------------------
    #Contour lines
    #---------------------------------------------------------------------------
    std_PV = np.std(PV_Data[0:100,0:100]) 
   
    PV_Contour_Levels = np.arange(1, 50, 1)*contour_interval*std_PV
    #print(PV_Contour_Levels)
    PV_Contour_Levels = [1000.0,2000.0]

    #---------------------------------------------------------------------------
    #Plotting
    #---------------------------------------------------------------------------
    #First we create arrays for the correct ticks
    x_values = np.arange(begin_v, begin_v + delta_v*float(Shape_Data[0]), delta_v)
    
    #The total length in arcsec of the y axis in the new image.
    length_arcsec_new =  (np.abs(np.cos(np.abs(rotation_rad)))*np.abs(PixelWidth_RA)*Shape_Data[2]+np.abs(np.cos(np.pi/2.0-np.abs(rotation_rad)))*np.abs(PixelWidth_DEC)*Shape_Data[1])*3600
    y_values = np.arange(-length_arcsec_new/2.0, length_arcsec_new/2.0 + length_arcsec_new/10.0, length_arcsec_new/10.0)

    #Calculating the size of y pixel in the y direction in arcsec.
    pixel_size_y_arcsec = length_arcsec_new/y_size

    y_arcsec = np.arange(1, PV_Data.shape[0])*pixel_size_y_arcsec - length_arcsec_new/2.0

    #Start plotting
    print('Starting the plotting.')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if imagemin !=-99.0 and imagemax !=-99.0:       
       cax = ax.imshow(PV_Data, origin = 'lower', cmap = 'magma',vmin=imagemin,vmax=imagemax,interpolation='nearest')
    else:
       cax = ax.imshow(PV_Data, origin = 'lower', cmap = 'magma',interpolation='nearest')

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
            Data_Contour_Levels = np.arange(0, 30, 1)*std_Data
            Data_Contour_Levels[0] = 3*std_Data
            
        ax.contour(Data, Data_Contour_Levels, colors = Color)
        ax.axvline(-10000, color = Color, label = Line)
        #ax.legend(loc = 3)
        
    if Overlay_Contour == '13CO':
        Directory = '13CO_1/PV-diagram-new/'
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
    ax.set_title(Object +' ' + Line)
    ax.set_xlabel('Velocity (km s$^{-1}$)')
    ax.set_ylabel('Offset ($^{\prime\prime}$)')#'Position [arcsec]')$\Delta$X

    offset = 0.0
    #If we have correct masses we can calculate the velocity curve.
    if Velocity_Curve:
        print('Including the velocities curves.')
        #Calculating the extreme masses within the errors, do we can also plot
        #those. 
        mass_min_err = mass - mass_err
        mass_plus_err = mass + mass_err
        
        #Creating an array containing the velocities in km s^-1. 
        velocities = (np.arange(begin_v, begin_v + delta_v*float(Shape_Data[0]), delta_v) - v_source)
        #print(velocities
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
        radii_positive_pixel_coor = radii_positive/pixel_size_y_arcsec + (y_size - 1.0)/2.0  + offset
        radii_negative_pixel_coor = radii_negative/pixel_size_y_arcsec + (y_size - 1.0)/2.0  + offset

        radii_positive_min_err_pixel_coor = radii_positive_min_err/pixel_size_y_arcsec + (y_size - 1.0)/2.0 + offset 
        radii_negative_min_err_pixel_coor = radii_negative_min_err/pixel_size_y_arcsec + (y_size - 1.0)/2.0 + offset

        radii_positive_plus_err_pixel_coor = radii_positive_plus_err/pixel_size_y_arcsec + (y_size - 1.0)/2.0 + offset  
        radii_negative_plus_err_pixel_coor = radii_negative_plus_err/pixel_size_y_arcsec + (y_size - 1.0)/2.0 + offset
        
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
        print('Current 4: ',(time()-_TIME_))
        _TIME_=time()
 
    #If zoom is true we zoom in else we just show the entire image. 
    if Zoom:
        #First we determine at which pixel we have the source velocity.
        pix_v_source = float(np.abs((begin_v - v_source)/delta_v))
        #Then we determine what half the width of the v slice must be.
        pix_v_shift = float(v_width/delta_v/2.0)
        #print(pix_v_source, pix_v_shift,v_width,delta_v
        #Now we determine the central pixel for the arcsec.
        pix_arcsec_central = float(y_size/2.0) - 1.0 + float(offset)
        pix_arcsec_shift = float(arcsec_width/pixel_size_y_arcsec/2.0)
        
        x = np.arange(pix_v_source - pix_v_shift, pix_v_source + pix_v_shift + 2.0*pix_v_shift/10.0, 2.0*pix_v_shift/10.0)
        #y = np.arange(pix_arcsec_central - pix_arcsec_shift, pix_arcsec_central + pix_arcsec_shift + 2*pix_arcsec_shift/10., 2*pix_arcsec_shift/10)
        y = np.linspace(pix_arcsec_central - pix_arcsec_shift,pix_arcsec_central + pix_arcsec_shift, 11.0)
       
        x_label = np.round(begin_v + delta_v*x, 1)
        #y_label = np.round((y+1)*pixel_size_y_arcsec - length_arcsec_new/2., 2) + 0.01
        y_label = np.linspace(-arcsec_width/2.0, arcsec_width/2.0, 11.0)



 
        from matplotlib.ticker import FormatStrFormatter
        #ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        y_label = np.array(["%.1f" % i for i in y_label])
        ax.set_xticks(x)
        ax.set_yticks(y)
        ax.set_xticklabels(x_label)
        ax.set_yticklabels(y_label)

        #Doing the zooming by limiting the shown x and y coordinates.
        #move after the labels, because those functions automatically adjust x and y limits.
        #print(pix_v_source - pix_v_shift, pix_v_source + pix_v_shift,pix_v_source
        #print('help!'
        ax.set_xlim([pix_v_source - pix_v_shift, pix_v_source + pix_v_shift])
        ax.set_ylim([pix_arcsec_central - pix_arcsec_shift, pix_arcsec_central + pix_arcsec_shift]) 
        #print(x
        ax.set_aspect(1.0*pix_v_shift/pix_arcsec_shift)

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
    print('Current 5: ',(time()-_TIME_))
    _TIME_=time()

    #ax.grid(color='white', alpha=0.5, linestyle='solid')

    cbar = plt.colorbar(cax)
    cbar.set_label('Jy/Beam')
    plt.draw()
    #ax.set_xticks(color='white', exclude_overlapping=True)
    #ax.set_yticks(color='white', exclude_overlapping=True)

    plt.savefig(save_directory + 'PV-Diagram_' + File[:-5] + '.pdf',bbox_inches='tight')
    #plt.savefig(save_directory + 'PV-Diagram_' + File[:-5] + '.jpg',bbox_inches='tight')
    #plt.savefig(save_directory + 'PV-Diagram_' + File[:-5] + '.ps',bbox_inches='tight')
    #plt.savefig(save_directory + 'PV-Diagram_' + File[:-5] + '.eps',bbox_inches='tight')
    plt.close()
    #---------------------------------------------------------------------------
    #Saving the image as a .fits
    #---------------------------------------------------------------------------
    HDU = fits.PrimaryHDU(PV_Data)
    HDR = HDU.header
    #Saving the correct header data. 

    #Need to create header.
    HDU.writeto(save_directory + 'PV-Diagram_' + File[:-5] + '.fits', overwrite = True)

    #---------------------------------------------------------------------------
    #Saving data in .txt file
    #---------------------------------------------------------------------------    
    np.savetxt(save_directory + 'PV-Diagram_' + File[:-5] + '.txt', PV_Data, delimiter = ' ')
    print('Finished: ',(time()-_TIME_))
    _TIME_=time()

if __name__ == '__main__':
    import numpy as np
    print('Current: ',_TIME_)
    #-------------------------------------------------------------------------------
    #C170 triplet
    #-------------------------------------------------------------------------------
    
    #rotation = 240 #degrees
    rotation = -62 #degrees
    Mass = 0.8
    Mass_err = 0.1 #Solar masses
    v_source = 4.7 #5#km s^-1
    distance_source = 230.#pc
    inclination = 45.
    aswidth = 10
    vwidth  = 10
    
    Object = 'L1448IRS3B'
    Line = 'C17O'
    width = 75 # casa wants num to be odd 21,51,<75>
    File = 'L1448IRS3B_C17O_image_taper1500k.image.fits'
    Directory = './'
    ra_orig = '03 25 36.321'
    dec_orig = '30 45 15.049'
    radecformat = 'icrs'
         
    #-------------------------------------------------------------------------------
    #C170 wide
    #-------------------------------------------------------------------------------
    '''
    rotation = -20. #degrees
    #rotation = 240. #degrees

    Mass = 10**0.187797#Solar masses 
    Mass = 0.5
    Mass_err = 0.1 #Solar masses
    v_source = 5.2 #5#km s^-1
    distance_source = 230.#pc
    inclination = -40.
    aswidth = 4
    vwidth  = 10

    Object = 'L1448IRS3A'
    Line = 'C17O'
    width = 11 # casa wants num to be odd
    File  = 'L1448IRS3B_C17O_clean_binned.fits'
    Directory = './'
    ra_orig = '03 25 36.497'
    dec_orig = '30 45 21.897'
    radecformat = 'icrs'
    '''
    #-------------------------------------------------------------------------------
    #H13CN wide
    #-------------------------------------------------------------------------------
    '''
    rotation = 350. #degrees

    Mass = 10**0.187797#Solar masses 
    Mass_err = 0.2 #Solar masses
    v_source = 5.4 #5#km s^-1
    distance_source = 230.#pc
    inclination = -40.
    inclination = 120
    aswidth = 4
    vwidth  = (11.6-5.4)*2.

    Object = 'L1448IRS3A'
    Line = 'H143CN'
    width = 21 # casa wants num to be odd
    #File  = 'L1448IRS3B_H13CN_image_taper1500k.fits'
    File = 'L1448IRS3B_H13CN_clean_image_2_binned.fits'
    Directory = './'
    ra_orig = '03 25 36.497'
    dec_orig = '30 45 21.897'
    radecformat = 'icrs'
    '''
    def icrs_hr(icrs):
        tmp0 = float(icrs[0])
        tmp01 = tmp0
        if tmp01 < 0:
            tmp0 = abs(tmp0)
        tmp1 = abs(float(icrs[1])/60.)
        tmp2 = abs(float(icrs[2])/3600.)
        icrs = tmp0 + tmp1 + tmp2
        if tmp01 < 0:
            icrs = icrs * -1.
        return icrs
    def hr_icrs(hr):
        tmp0 = int(hr)
        tmp1 = int((hr - float(tmp0))*60.)
        tmp2 = float(((hr - float(tmp0))*60. - tmp1)*60.)
        icrs = str(tmp0) + ":" + str(abs(tmp1)) + ":" + "{0:.2f}".format(round(abs(tmp2),2))
        return icrs
    while True:
        try:
            ra_fin = float(ra_orig)
            dec_fin = float(dec_orig)
        except ValueError:
            try:
                ra_fin = ra_orig.split(':')
                dec_fin = dec_orig.split(':')
                if len(ra_fin) == 1:
                    ra_fin =  ra_orig.split(' ')            
                if len(dec_fin) == 1:
                    dec_fin =  dec_orig.split(' ')
                if ((len(dec_fin) == 3) and (len(ra_fin) == 3)):
                    break
                else:
                    continue
            except ValueError:
                continue
            continue
        if radecformat in ['deg', 'hourangle', 'icrs']:
            break
        else:
            print("Please input valid format.")
            continue


    if (type(ra_fin) is list) and (radecformat == 'icrs'):
        print('converting to deg:')
        print('RA:  ' + str(icrs_hr(ra_fin)*15.))
        print('DEC: ' + str(icrs_hr(dec_fin)))
        # convert to 

    ra=icrs_hr(ra_fin)*15.
    dec=icrs_hr(dec_fin)


    rar=ra
    decr=dec

    PV_diagram(rar,decr,File, Directory, rotation, width, Object,\
     Line, inclination, Velocity_Curve = True, mass = Mass, mass_err = Mass_err,\
      v_source = v_source, d_source = distance_source, Thindisk = False, Zoom = True,\
       v_width = vwidth, arcsec_width =aswidth,Overlay_Contour = 'None')
    
    print('Mass:',Mass)
