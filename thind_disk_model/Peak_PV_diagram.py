from __future__ import division
import numpy as np
import numpy.ma
import scipy.optimize
import pyfits as pf
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.ndimage import interpolation as inter
from astropy import units as u
from astropy.coordinates import SkyCoord
from Kinematic_Centrum import *
from astropy.io import fits
from astropy import constants as const
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FormatStrFormatter

#-------------------------------------------------------------------------------
#2D gaussian fitting functions
#Source of this code:
#https://github.com/indiajoe/HandyTools4Astronomers/blob/master/GaussianFit.py
#-------------------------------------------------------------------------------
def moments2D(inpData):
    """ Returns the (amplitude, xcenter, ycenter, xsigma, ysigma, rot, bkg, e) estimated from moments in the 2d input array Data """

    bkg=np.median(np.hstack((inpData[0,:],inpData[-1,:],inpData[:,0],inpData[:,-1])))  #Taking median of the 4 edges points as background
    Data=np.ma.masked_less(inpData-bkg,0)   #Removing the background for calculating moments of pure 2D gaussian
    #We also masked any negative values before measuring moments

    amplitude=Data.max()

    total= float(Data.sum())
    Xcoords,Ycoords= np.indices(Data.shape)

    xcenter= (Xcoords*Data).sum()/total
    ycenter= (Ycoords*Data).sum()/total

    RowCut= Data[int(xcenter),:]  # Cut along the row of data near center of gaussian
    ColumnCut= Data[:,int(ycenter)]  # Cut along the column of data near center of gaussian
    xsigma= np.sqrt(np.ma.sum(ColumnCut* (np.arange(len(ColumnCut))-xcenter)**2)/ColumnCut.sum())    
    ysigma= np.sqrt(np.ma.sum(RowCut* (np.arange(len(RowCut))-ycenter)**2)/RowCut.sum())

    #Ellipcity and position angle calculation
    Mxx= np.ma.sum((Xcoords-xcenter)*(Xcoords-xcenter) * Data) /total
    Myy= np.ma.sum((Ycoords-ycenter)*(Ycoords-ycenter) * Data) /total
    Mxy= np.ma.sum((Xcoords-xcenter)*(Ycoords-ycenter) * Data) /total
    e= np.sqrt((Mxx - Myy)**2 + (2*Mxy)**2) / (Mxx + Myy)
    pa= 0.5 * np.arctan(2*Mxy / (Mxx - Myy))
    rot= np.rad2deg(pa)

    return amplitude,xcenter,ycenter,xsigma,ysigma, rot,bkg, e
    
def Gaussian2D(amplitude, xcenter, ycenter, xsigma, ysigma, rot, bkg):
    """ Returns a 2D Gaussian function with input parameters. rotation input rot should be in degress """
    rot=np.deg2rad(rot)  #Converting to radians
    Xc=xcenter*np.cos(rot) - ycenter*np.sin(rot)  #Centers in rotated coordinates
    Yc=xcenter*np.sin(rot) + ycenter*np.cos(rot)
    #Now lets define the 2D gaussian function
    def Gauss2D(x,y) :
        """ Returns the values of the defined 2d gaussian at x,y """
        xr=x * np.cos(rot) - y * np.sin(rot)  #X position in rotated coordinates
        yr=x * np.sin(rot) + y * np.cos(rot)
        return amplitude*np.exp(-(((xr-Xc)/xsigma)**2 +((yr-Yc)/ysigma)**2)/2) +bkg

    return Gauss2D 

def FitGauss2D(Data,ip=None):
    """ Fits 2D gaussian to Data with optional Initial conditions ip=(amplitude, xcenter, ycenter, xsigma, ysigma, rot, bkg)
    Example: 
    >>> X,Y=np.indices((40,40),dtype=np.float)
    >>> Data=np.exp(-(((X-25)/5)**2 +((Y-15)/10)**2)/2) + 1
    >>> FitGauss2D(Data)
    (array([  1.00000000e+00,   2.50000000e+01,   1.50000000e+01, 5.00000000e+00,   1.00000000e+01,   2.09859373e-07, 1]), 2)
     """
    if ip is None:   #Estimate the initial parameters form moments and also set rot angle to be 0
        ip=moments2D(Data)[:-1]   #Remove ellipticity from the end in parameter list

    Xcoords,Ycoords= np.indices(Data.shape)    
    def errfun(ip):
        dXcoords= Xcoords-ip[1]
        dYcoords= Ycoords-ip[2]
        Weights=np.sqrt(np.square(dXcoords)+np.square(dYcoords)) # Taking radius as the weights for least square fitting

        return np.ravel((Gaussian2D(*ip)(*np.indices(Data.shape)) - Data)/np.sqrt(Weights))  #Taking a sqrt(weight) here so that while scipy takes square of this array it will become 1/r weight.

    p, success = scipy.optimize.leastsq(errfun, ip, maxfev = 100000)
    return p,success

#-------------------------------------------------------------------------------
#My own functions
#-------------------------------------------------------------------------------
def Plot_Point_in_Channel(x_center, y_center, x_center_err, y_center_err, image_channel, channel, directory, v_obs, v_rel, boxsize = 80):
    #Creating the directory for saving the images, if it does not already exists.
    save_directory = directory + 'Channels/' 
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
            
    #Creating the figure.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    cax = ax.imshow(image_channel, cmap = 'magma', origin = 'lower')
    ax.errorbar(x_center, y_center, xerr = x_center_err, yerr= y_center_err, fmt='.', color = 'red')
    ax.set_xlabel('x [pixels]')
    ax.set_ylabel('y [pixels]')
    ax.set_title('Channel ' + str(channel) + ', v$_{obs}$ = ' + str(np.round(v_obs, 2)) + ' km s$^{-1}$, v$_{rel}$ = ' + str(np.round(v_rel, 2)) + 'km s$^{-1}$' )
    ax.set_xlim([x_center - boxsize/2, x_center + boxsize/2])
    ax.set_ylim([y_center - boxsize/2, y_center + boxsize/2])
    
    cbar = plt.colorbar(cax)
    cbar.set_label('Jy/beam')
    
    plt.savefig(save_directory + 'Channel_' + str(channel) + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Channel_' + str(channel) + '.pdf', bbox_inches='tight')
    plt.close()

def Calculate_Line(X, angle):
    #The starting y positions.
    Y = np.zeros_like(X)
    
    #Converting the angle from degrees to radians.
    angle = np.radians(angle)
    
    #rotating the line
    x = np.cos(angle)*X - np.sin(angle)*Y
    y = np.sin(angle)*X + np.cos(angle)*Y
    
    return x, y
    
def Plot_Positions(x_obj_pix, y_obj_pix, x_obj_err_pix, y_obj_err_pix, x_pos_pix, y_pos_pix, x_pos_err_pix, y_pos_err_pix, directory, delta_x_deg, delta_y_deg, boxsize_arcsec = 2):
    #Creating the directory for saving the images, if it does not already exists.
    save_directory = directory + 'Fitted_Positions/' 
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
    
    delta_x_arcsec = delta_x_deg*3600.0
    delta_y_arcsec = delta_y_deg*3600.0
    
    x_obj_arcsec = (x_obj_pix - x_obj_pix)*np.abs(delta_x_arcsec)
    y_obj_arcsec = (y_obj_pix - y_obj_pix)*np.abs(delta_y_arcsec)
    
    x_pos_arcsec = (x_pos_pix - x_obj_pix)*np.abs(delta_x_arcsec)
    y_pos_arcsec = (y_pos_pix - y_obj_pix)*np.abs(delta_y_arcsec)
    
    x = np.linspace(-boxsize_arcsec,  boxsize_arcsec, 100)
    
    PA = -90.#degrees
    
    disk_line_x, disk_line_y = Calculate_Line(x, PA)
    outflow_line_x, outflow_line_y = Calculate_Line(x, PA + 90.)
     
    #Creating the figure.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    ax.plot(disk_line_x, disk_line_y, color = 'black')
    ax.plot(outflow_line_x, outflow_line_y, color = 'black')
    
    #ax.plot(x_obj_arcsec, y_obj_arcsec, 'x', color = 'red')
    ax.plot(x_pos_arcsec, y_pos_arcsec, 'x', color = 'blue')
            
    ax.set_xlabel('$\Delta$X [arcsec]')
    ax.set_ylabel('$\Delta$Y [arcsec]')
    ax.set_title('Fitted positions of Peak PV diagram' )
    
    ax.set_aspect('equal', 'datalim')
    
    ax.set_xlim([ -boxsize_arcsec/2,  boxsize_arcsec/2])
    ax.set_ylim([ -boxsize_arcsec/2,  boxsize_arcsec/2])
    
    plt.savefig(save_directory + 'Fitted_Positions.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Fitted_Positions.pdf', bbox_inches='tight')
    plt.close()  
     
def AU_Distance(x_object, y_object, x_image, y_image, distance_source, pixel_width_x, pixel_width_y):  
    #This function returns the distance in AU of given coordinates to the center of
    #the image. The distance to the source (distance_source)must be given in parsec.
    #The pixelwidths must be given in radians. 
    
    #First the distance in radians is calculated.         
    distance_radians = np.sqrt(np.square((x_image - x_object)*pixel_width_x) + np.square((y_image - y_object)*pixel_width_y))
    
    #Now calculate the distance in AU.
    distance_AU = distance_radians*distance_source*u.pc.to(u.AU)
    return distance_AU

def AU_Distance_standard_deviation(x_object, y_object, x_image, y_image, distance_source, pixel_width_x, pixel_width_y, x_object_err, y_object_err, x_image_err, y_image_err, distance_source_err):
    #This function calculates the standard deviation / error on the distance.
    #Defining a term that will return in a lot of other terms.
    T = np.square((x_image - x_object)*pixel_width_x) + np.square((y_image - y_object)*pixel_width_y)

    #A commonly used constant.
    c = u.pc.to(u.AU)
    
    #All the different derivatives. 
    A = distance_source*c*np.square(pixel_width_x)/np.sqrt(T)*(x_image - x_object)

    B = distance_source*c*np.square(pixel_width_y)/np.sqrt(T)*(y_image - y_object)

    C = distance_source*c*np.square(pixel_width_x)/np.sqrt(T)*(x_object - x_image)

    D = distance_source*c*np.square(pixel_width_y)/np.sqrt(T)*(y_object - y_image)

    E = np.sqrt(T)*c

    #Now we add all the different terms and calculate the errors. 
    distance_AU_err = np.sqrt(np.square(A*x_image_err) + np.square(B*y_image_err) + np.square(C*x_object_err) + np.square(D*y_object_err) + np.square(E*distance_source_err))
    return distance_AU_err
    
def Gaussian_Center_Position_Standard_Deviation(var_Data, Amplitude_fit, sigma_x_fit, sigma_y_fit, beam_size_x, beam_size_y):
    #Using the information in the article 'Errors in Elliptical Gaussian Fits' 
    #written by J. J. Condon we determine the standard deviation on the position of
    #the fitted Gaussian with correlated noise. 
    theta_x = np.sqrt(8*np.log(2))*sigma_x_fit
    theta_y = np.sqrt(8*np.log(2))*sigma_y_fit

    rho_x_sq = (theta_x*theta_y)/(4*beam_size_x*beam_size_y)*(1 + np.square(beam_size_x/theta_x))**(2.5)*(1 + np.square(beam_size_y/theta_y))**(0.5)*(np.square(Amplitude_fit)/var_Data)
    rho_y_sq = (theta_x*theta_y)/(4*beam_size_x*beam_size_y)*(1 + np.square(beam_size_x/theta_x))**(0.5)*(1 + np.square(beam_size_y/theta_y))**(2.5)*(np.square(Amplitude_fit)/var_Data)

    std_x = np.sqrt(2/rho_x_sq)*sigma_x_fit
    std_y = np.sqrt(2/rho_y_sq)*sigma_y_fit

    return std_x, std_y
    
def Peak_PV_Diagram(File, Directory, Object, Line, begin_1, end_1, begin_2, end_2, v_source, distance_source, distance_source_err, Thindisk = False, Used_M0 = 'Own'):
    #This function creates an UV position velocity diagram for a given fits file. 
    #Definition of the parameters:
    #File = name of the datacube (.fits file). [string]
    #Directory = directory where the file is saved. [string]
    #Object = name of the object we are looking at. [string]
    #line = the used molecular line. [string]
    #begin_1 = [float]
    #end_1 = [float]
    #begin_2 = [float]
    #end_2 = [float]
    #v_source = [float]
    #distance_source = [float]
    #distance_source_err = [float]
    #Thindisk = [boolean] 
    #Used_M0 = the used M0 map for determining the kinematic center. [string]
    #          We have the   
    print ' '     
    print 'Creating the Peak PV-diagram for ' + File[:-5]
    #---------------------------------------------------------------------------
    #Getting the data from the fits file.
    #---------------------------------------------------------------------------
    #Open the fits file. 
    hdu = fits.open(Directory + File)[0]
    #Putting the data and header in seperate files.
    header = hdu.header
    data = hdu.data
    
    #getting the header.
    #header = pf.getheader(Directory + File)

    #Getting the data.
    #data = pf.getdata(Directory + File)

    if Thindisk:
        Data = data
    else:
        Data = data[0,:,:,:]

    #Checking for nan values and putting them on zero.
    Data[np.isnan(Data)] = 0 

    Shape_Data = np.shape(Data)

    #Determining the begin_v, delta_v and the amount of images N.
    N = header['NAXIS3']
    if header['CTYPE3']  == 'VELO-LSR':
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
        delta_v = begin_v - begin_v_plus_one
    else:
        print 'I am not sure how to get the velocities.'
        raise KeyError
    
    #Pixel dimensions in degrees.
    PixelWidth_RA = header['CDELT1']
    PixelWidth_DEC = header['CDELT2']
    
    #Converting it to radians
    PixelWidth_RA_rad = np.radians(np.abs(PixelWidth_RA))
    PixelWidth_DEC_rad = np.radians(PixelWidth_DEC)

    #Beam size
    Beam_x_deg = header['BMAJ']
    Beam_y_deg = header['BMIN']

    #Converting it to pixels
    Beam_x_pix = Beam_x_deg/np.abs(PixelWidth_RA)
    Beam_y_pix = Beam_y_deg/PixelWidth_DEC

    #Creating the directory for saving the images, if it does not already exists.
    save_directory = Directory + 'Peak_PV-diagram/' 
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)

    #---------------------------------------------------------------------------
    #Fitting the Gaussians and creating the data. 
    #--------------------------------------------------------------------------- 
    #Red shifted part of the data.    
    radius_red = np.zeros(Shape_Data[0])
    velocity_red = np.zeros(Shape_Data[0])

    #Blue shifted part of the data.
    radius_blue = np.zeros(Shape_Data[0])
    velocity_blue = np.zeros(Shape_Data[0])

    #Errors in de radii.
    radius_red_err = np.zeros(Shape_Data[0]) 
    radius_blue_err = np.zeros(Shape_Data[0])

    #Choosing the M0 map to determine the kinematic centrum.
    if Used_M0 == 'Own':
        Path_Data_Centrum = Directory + 'Moments_Images/Moment_0_map_' + File[:-5] + '.fits'
    if Used_M0 == 'Cont':
        Path_Data_Centrum = 'Cont_2/L1527_Cont_2.fits'

    x_object_deg, y_object_deg, x_object_deg_err, y_object_deg_err = Find_Centrum_M0_Map(Path_Data_Centrum, Print_Output = True)

    x_object, y_object, x_object_err, y_object_err = Degrees_to_Pixels(header, x_object_deg, y_object_deg, x_object_deg_err, y_object_deg_err)

    print 'Pixel coordinates center total object:'
    print 'x: ' + str(x_object)
    print 'y: ' + str(y_object)
    print 'The errors in the position determination:'
    print 'x error: ' + str(x_object_err)
    print 'y error: ' + str(y_object_err)
    
    x_positions = np.zeros_like(Data[:,0,0])
    x_positions_err = np.copy(x_positions)
    y_positions = np.copy(x_positions)
    y_positions_err = np.copy(x_positions)

    print 'Starting with the blueshifted images...'
    #While loop through the first set of images (red shifted).
    i = begin_image_1
    box_size = 25
    while i < end_image_1 + 1:
        print str(i) + '/' + str(end_image_1)
        #To make sure we only fit relevant data.
        Data_for_Plot = np.copy(Data[i,:,:])
        Used_Data = Data[i,:,:]
        std_Data = np.std(Used_Data)

        #y_data_used,x_data_used  = np.unravel_index(Used_Data.argmax(), Used_Data.shape)
    
        x_data_used = x_object
        y_data_used = y_object

        Used_Data[:y_data_used-box_size, :] = 0
        Used_Data[y_data_used+box_size:, :] = 0
        Used_Data[:,:y_data_used-box_size] = 0
        Used_Data[:, y_data_used+box_size:] = 0
        
        #Finding the center of the fitted Gaussian.
        temp = FitGauss2D(Used_Data)
        parameters_Gaussian = temp[0]
        parameters_err_Gaussian = temp[1]
        x_image = parameters_Gaussian[2]
        y_image = parameters_Gaussian[1]
        
        x_image_err, y_image_err = Gaussian_Center_Position_Standard_Deviation(np.var(Used_Data), parameters_Gaussian[0], parameters_Gaussian[4], parameters_Gaussian[3], Beam_x_pix, Beam_y_pix)

        #Saving the measured positions.
        x_positions[i] = x_image
        x_positions_err[i] = x_image_err
        y_positions[i] = y_image
        y_positions_err[i] = y_image_err
        
        #Calculating the radius and velocity.
        radius_blue[i] = AU_Distance(x_object, y_object, x_image, y_image, distance_source, PixelWidth_RA_rad, PixelWidth_DEC_rad)
        velocity_blue[i] = begin_v + i*delta_v - v_source
        radius_blue_err[i] = AU_Distance_standard_deviation(x_object, y_object, x_image, y_image, distance_source, PixelWidth_RA_rad, PixelWidth_DEC_rad, x_object_err, y_object_err, x_image_err, y_image_err, distance_source_err)

        #The velocities necessary of the channel plots.
        V_obs = begin_v + i*delta_v
        V_rel = V_obs - v_source

        Plot_Point_in_Channel(x_image, y_image, x_image_err, y_image_err, Data_for_Plot, i, save_directory, V_obs, V_rel)
        
        i += 1
        #End while loop  

    print 'Starting with the redshifted images...'
    #While loop through the second set of images (blue shifted).
    i = begin_image_2
    while i < end_image_2 + 1:
        print str(i) + '/' + str(end_image_2)
        #To make sure we only fit relevant data.
        Data_for_Plot = np.copy(Data[i,:,:])
        Used_Data = Data[i,:,:]
        std_Data = np.std(Used_Data)

        #Shape data.
        y_data_used, x_data_used  = np.unravel_index(Used_Data.argmax(), Used_Data.shape)
    
        Used_Data[:y_data_used-box_size, :] = 0
        Used_Data[y_data_used+box_size:, :] = 0
        Used_Data[:,:y_data_used-box_size] = 0
        Used_Data[:, y_data_used+box_size:] = 0
        
        #Finding the center of the fitted Gaussian.
        temp = FitGauss2D(Used_Data)
        parameters_Gaussian = temp[0]
        parameters_err_Gaussian = temp[1]
        x_image = parameters_Gaussian[2]
        y_image = parameters_Gaussian[1]

        x_image_err, y_image_err = Gaussian_Center_Position_Standard_Deviation(np.var(Used_Data), parameters_Gaussian[0], parameters_Gaussian[4], parameters_Gaussian[3], Beam_x_pix, Beam_y_pix)

        #Saving the measured positions.
        x_positions[i] = x_image
        x_positions_err[i] = x_image_err
        y_positions[i] = y_image
        y_positions_err[i] = y_image_err
            
        #Calculating the radius and velocity.
        radius_red[i] = AU_Distance(x_object, y_object, x_image, y_image, distance_source, PixelWidth_RA_rad, PixelWidth_DEC_rad)
        velocity_red[i] = begin_v + i*delta_v - v_source
        radius_red_err[i] = AU_Distance_standard_deviation(x_object, y_object, x_image, y_image, distance_source, PixelWidth_RA_rad, PixelWidth_DEC_rad, x_object_err, y_object_err, x_image_err, y_image_err, distance_source_err)

        #The velocities necessary of the channel plots.
        V_obs = begin_v + i*delta_v
        V_rel = V_obs - v_source

        Plot_Point_in_Channel(x_image, y_image, x_image_err, y_image_err, Data_for_Plot, i, save_directory, V_obs, V_rel)
        
        i += 1
        #End while loop 
    
    #Setting all the elements in the arrays that werent used to NaN.
    radius_red[radius_red == 0] = np.nan
    velocity_red[velocity_red == 0] = np.nan

    radius_blue[radius_blue == 0] = np.nan
    velocity_blue[velocity_blue == 0] = np.nan

    radius_red_err[radius_red_err == 0] = np.nan
    radius_blue_err[radius_blue_err == 0] = np.nan
    
    #We are only going to use the values which are not zero.
    select = (x_positions > 0)
    x_positions = x_positions[select]
    x_positions_err = x_positions_err[select]
    y_positions = y_positions[select]
    y_positions_err = y_positions_err[select]
    
    print 'Plotting the positions of the fits on sky.'
    Plot_Positions(x_object, y_object, x_object_err, y_object_err, x_positions, y_positions, x_positions_err, y_positions_err, save_directory, PixelWidth_RA, PixelWidth_DEC)

    #---------------------------------------------------------------------------
    #Determining the error bars on the data.
    #---------------------------------------------------------------------------
    velocity_err = delta_v

    #---------------------------------------------------------------------------
    #Plotting
    #---------------------------------------------------------------------------
    print 'Starting the plotting.'
    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.errorbar(radius_red, np.abs(velocity_red), xerr = radius_red_err, yerr= velocity_err, fmt='.', color = 'red')
    ax.errorbar(radius_blue, np.abs(velocity_blue), xerr = radius_blue_err, yerr= velocity_err, fmt='.', color = 'blue')

    ax.set_title('The Peak PV diagram of '+ Object + ' in ' + Line)
    ax.set_xlabel('Radius (AU)')
    ax.set_ylabel('Radial Velocity - $v_{source}$ (km s$^{-1}$)')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlim([29, 101])
    ax.set_ylim([0.9, 3.1])
    
    #Makes sure that the minor ticks have labels.
    ax.xaxis.set_minor_formatter(FormatStrFormatter("%.0f"))
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))

    #Sets the ticks to 10,100 etc instead of 10^0, 10^1, 10^2,...
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_major_formatter(ScalarFormatter())

    #Rotates the minor ticks.
    labels = ax.get_xticklabels('minor')
    for label in labels:
        label.set_rotation(45) 
    #Rotates the major ticks.
    labels = ax.get_xticklabels()
    for label in labels:
        label.set_rotation(45) 

    ax.grid(b=True, which='major', color='black', linestyle='--')
    ax.grid(b=True, which='minor', color='black', linestyle='--')

    plt.savefig(save_directory + 'Peak_PV-Diagram_' + File[:-5] + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Peak_PV-Diagram_' + File[:-5] + '.pdf', bbox_inches='tight')
    plt.close()

    #---------------------------------------------------------------------------
    #Saving data in .txt file
    #---------------------------------------------------------------------------
    save_data = np.zeros([Shape_Data[0], 8])
    save_data[:,0] = radius_red
    save_data[:,1] = radius_red_err
    save_data[:,2] = np.abs(velocity_red)
    save_data[:,3] = velocity_err  
    save_data[:,4] = radius_blue
    save_data[:,5] = radius_blue_err
    save_data[:,6] = np.abs(velocity_blue)
    save_data[:,7] = velocity_err     
    np.savetxt(save_directory + 'Peak_PV-Diagram_' + File[:-5] + '.txt', save_data, delimiter = ' ') 

if __name__ == '__main__':
    Object = 'L1527'
    v_source = 5.9#km s^-1
    distance_source = 140 #pc
    distance_source_err = 0 # pc <-is zero because it would be a systematic error.

    #C180
    Line = 'C$^{18}$O'

    begin_image_1 = 18
    end_image_1 = 30
    begin_image_2 = 44
    end_image_2 = 53

    File = 'L1527_C18O_1.fits'
    Directory = 'C18O_1/'
    #Peak_PV_Diagram(File, Directory, Object, Line, begin_image_1, end_image_1, begin_image_2, end_image_2, v_source, distance_source, distance_source_err, Used_M0 = 'Own')

    File = 'L1527_C18O_2.fits'
    Directory = 'C18O_2/'
    Peak_PV_Diagram(File, Directory, Object, Line, begin_image_1, end_image_1, begin_image_2, end_image_2, v_source, distance_source, distance_source_err, Used_M0 = 'Cont')


