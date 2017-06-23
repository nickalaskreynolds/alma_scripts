import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import constants as const
from astropy import units as u
import numpy as np
import os

def zoom(im,x,y,bb):
    '''
    Cut out a square box from image im centered on (x,y) with half-box size
     bb.
    '''
    if len(im.shape) == 3:
        return im[:,y-bb:y+bb,x-bb:x+bb]
    elif len(im.shape) == 2:
        return im[y-bb:y+bb,x-bb:x+bb]
    else:
        print 'Wrong number of dimensions! N must be 2 or 3.'

def Get_Velocities(header):
    '''
    This function determines the velocity for the first channel, the 
    velocity step and the amount of channels. 
    header = the header of the file.
    '''
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
        
    return begin_v, delta_v, N 

def Compare_Channels(Directory1, Directory2, File1, File2, Object, Line1, Line2, channel_mask, Save_Name, Contour = False, boxsize = 50):
    '''
    This function plots two datacubes, for every channel, side by side in 
    the central region. But only for specific channels. We assume that the 
    two datacubes span the same velocities and that the first file is the 
    science data and the second the model data. The residual is also 
    plotted.
    Directory1 = directory where file1 is saved. [string]
    Directory2 = directory where file2 is saved. [string]
    File1 = name of the first file. [string]
    File2 = name of the second file. [string]
    Object = Object name. [string] 
    Line1 = Molecular line of file1 [string]
    Line2 = Molecular line of file2 [string]
    Contour = if true only the contour will be plotted, if false only the data
    and no contours will be plotted [Boolean]
    channel_mask = mask for the channels that should be compared. [boolean] 
    boxsize = half the size of the box in pixels that will be shown [int]
    '''
    plt.rcParams['text.latex.preamble']=[r"\usepackage{lmodern}"]
    #Options
    params = {'text.usetex' : True,
              'font.size' : 15,
              'font.family' : 'lmodern',
              'text.latex.unicode': True,
              }
    plt.rcParams.update(params) 
    #---------------------------------------------------------------------------
    #Getting the data from the fits files.
    #---------------------------------------------------------------------------    
    #Open the fits files. 
    hdu1 = fits.open(Directory1 + File1)[0]
    hdu2 = fits.open(Directory2 + File2)[0]
    
    #Putting the data and header in seperate files.
    header1 = hdu1.header
    header2 = hdu2.header 
       
    data1 = hdu1.data
    data2 = hdu2.data
    
    #Just to check if we have to remove the first axis of the datacube.
    if len(data1.shape) == 3:
        Data1 = data1
    else:
        Data1 = data1[0,:,:,:]
    if len(data2.shape) == 3:
        Data2 = data2
    else:
        Data2 = data2[0,:,:,:]

    #Checking for nan values and putting them on zero.
    Data1[np.isnan(Data1)] = 0 
    Data2[np.isnan(Data2)] = 0 
    
    #Getting the begin velocities.
    begin_v1, delta_v1, N1 = Get_Velocities(header1)
    
    #Pixelsizes
    x_pixel_size1 = np.abs(header1['CDELT1'])#arcsec
    y_pixel_size1 = header1['CDELT2']#arcsec
    x_pixel_size2 = np.abs(header2['CDELT1'])#arcsec
    y_pixel_size2 = header2['CDELT2']#arcsec
    
    #-------------------------------------------------------------------
    #Calculating the used velocities and slicing the datacubes.
    
    Channels = np.arange(0, N1, 1)
    
    Velocities = begin_v1 + Channels[channel_mask]*delta_v1
    
    #Determining the vmin and vmax for plotting. So it is the same for 
    #all images.
    Max = np.max(Data1)#np.max(np.concatenate((Data1, Data2), axis = 0))
    Min = np.min(Data1)#np.min(np.concatenate((Data1, Data2), axis = 0))
    
    #The standard deviation of the data. Used for contours.
    std_Data1 = np.std(Data1, axis = (1,2))
    std_Data1 = std_Data1[:, np.newaxis]
    
    #The levels for the contours.
    Levels = np.arange(3, 51, 3)
    Levels = Levels[np.newaxis, :]
    
    #The final contour levels.
    Contour_Levels = Levels*std_Data1
    
    Data1 = Data1[channel_mask,:,:]  
    Data2 = Data2[channel_mask,:,:]
    
    #We only take the central part.
    xc = Data1.shape[2]/2
    yc = Data1.shape[1]/2
    Data1 = zoom(Data1,xc,yc,boxsize)
    Data2 = zoom(Data2,xc,yc,boxsize)
    Data3 = Data1 - Data2#The residuals
    #-------------------------------------------------------------------
    #Starting the plotting.
    
    #The amount of columns and rows in the figure.    
    columns = 4
    rows = int(np.round(len(Velocities)/(columns*1.)))*3

    #Size of the figure.
    ywidth = 9.7#inches    
    xwidth = 6.6#inches 
    
    #The values needed for extent.
    hmin = -boxsize*x_pixel_size1*u.deg.to(u.arcsec)
    hmax = boxsize*x_pixel_size1*u.deg.to(u.arcsec)
    vmin = -boxsize*y_pixel_size1*u.deg.to(u.arcsec)
    vmax = boxsize*y_pixel_size1*u.deg.to(u.arcsec)
    
    fig, axes = plt.subplots(nrows=rows, ncols=columns, sharex=True, sharey=True, figsize=(xwidth,ywidth))
    fig.subplots_adjust(hspace=0, wspace = 0)
    
    #Plotting all the figures at the correct positions.

    #For the data.
    j = 0
    k = 0
    m = 0
    
    #For the row loop.
    i = 0
    while i < rows:
        #For the columns.
        l = 0
        while l < columns:
            ax = axes[i,l]
            #Converting from pixel coordinates to the new coordinates.
            ytext = 0.75*boxsize*y_pixel_size1*u.deg.to(u.arcsec)
            xtext = -0.95*boxsize*x_pixel_size1*u.deg.to(u.arcsec)
            #Properties for the box we use for the text in the image.
            props = dict(facecolor= 'white')    
            if i % 3 == 0:
                if Contour:
                    ax.contour(Data1[j,:,:], Contour_Levels[j,:], origin = 'lower', extent = [hmin, hmax, vmin, vmax], colors = 'red')
                else:
                    ax.imshow(Data1[j,:,:], origin = 'lower', extent = [hmin, hmax, vmin, vmax], vmin = Min, vmax = Max, cmap = 'red')
                if j == 0:
                    ax.text(xtext, -ytext, Line1, fontsize = 8, color = 'black', bbox = props)
                    ax.text(xtext, ytext, str(np.round(Velocities[j],1)) + ' km s$^{-1}$', fontsize = 9, color = 'black', bbox = props)
                else:
                    ax.text(xtext, ytext, str(np.round(Velocities[j],1)) , fontsize = 9, color = 'black', bbox = props)
                j += 1
            elif i % 3 == 1:
                if Contour:
                    ax.contour(Data2[k,:,:], Contour_Levels[k,:], origin = 'lower', extent = [hmin, hmax, vmin, vmax], colors = 'red')                  
                else:
                    ax.imshow(Data2[k,:,:], origin = 'lower', extent = [hmin, hmax, vmin, vmax], vmin = Min, vmax = Max, cmap = 'magma')                 
                if k == 0:
                    ax.text(xtext, ytext, str(np.round(Velocities[k],1)) + ' km s$^{-1}$', fontsize = 9, color = 'black', bbox = props)
                    ax.text(xtext, -ytext, Line2, fontsize = 8, color = 'black', bbox = props)
                else:
                    ax.text(xtext, ytext, str(np.round(Velocities[k],1)), fontsize = 9, color = 'black', bbox = props)
                k += 1
            elif i % 3 == 2:
                if Contour:
                    ax.contour(Data3[m,:,:], Contour_Levels[m,:], origin = 'lower', extent = [hmin, hmax, vmin, vmax], colors = 'red')
                    ax.contour(-Data3[m,:,:], Contour_Levels[m,:], origin = 'lower', extent = [hmin, hmax, vmin, vmax], colors = 'blue')
                else:
                    Image_Channel = ax.imshow(Data3[m,:,:], origin = 'lower', extent = [hmin, hmax, vmin, vmax], vmin = Min, vmax = Max, cmap = 'magma')
                if m == len(Velocities) - columns:
                    ax.set_xlabel("$\Delta$x ['']")
                    ax.set_ylabel("$\Delta$y ['']")
                if m == len(Velocities) - 1 and Contour == False:
                    box = ax.get_position()
                    axColor = plt.axes([box.x0 + box.width, box.y0, 0.025, box.height])
                    cbar = plt.colorbar(Image_Channel, cax = axColor, orientation="vertical")
                    cbar.set_label('Jy/beam')
                if m == 0:
                    ax.text(xtext, ytext, str(np.round(Velocities[m],1)) + ' km s$^{-1}$', fontsize = 9, color = 'black', bbox = props)
                    ax.text(xtext, -ytext, 'Residual', fontsize = 8, color = 'black', bbox = props)
                else:
                    ax.text(xtext, ytext, str(np.round(Velocities[m],1)), fontsize = 9, color = 'black', bbox = props)
                m += 1
            ax.grid(color='black', alpha=1, linestyle='--')
            l += 1
        i += 1
    
    #Creating the directory for saving the images, if it does not already exists.
    save_directory = Directory2 + 'Channel_Comparison/' 
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
        
    plt.savefig(save_directory + Save_Name + '.jpg',bbox_inches='tight')
    plt.savefig(save_directory + Save_Name + '.pdf',bbox_inches='tight')
    plt.close()
    
if __name__ == '__main__':
    Object = 'L1527'
    
    Angle = '89'
    
    Directory1 = 'C18O_2/'
    Directory2 = 'Tests/Test_1/'
    
    File1 = 'L1527_C18O_2.fits'
    File2 = 'Test_1_convolved_with_beam.fits'
    
    Line1 = 'C$^{18}$O'
    Line2 = 'C$^{18}$O (Thindisk model)'
    
    begin_image = np.arange(0,96, 8)
    end_image = begin_image + 8
        
    i = 0
    while i < len(begin_image):
        print 'Channel ' + str(begin_image[i]) + ' to ' + str(end_image[i])

        mask = np.zeros(100, dtype = bool)
        mask[begin_image[i]:end_image[i]] = True
        
        save_name = 'Channel_' + str(begin_image[i]) + '-' +str(end_image[i]) + '_i_' + Angle + '_' + Line1
        
        Compare_Channels(Directory1, Directory2, File1, File2, Object, Line1, Line2, mask, save_name, Contour = False)
        
        i += 1
	
