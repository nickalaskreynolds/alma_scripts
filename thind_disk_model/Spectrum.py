import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import constants as const
import os
'''
This function calculates the spectrum of the given datacube. 
File = [string]
Directory = [string]
Boxsize = [int]
'''
def Calculate_Spectrum(File, Directory, Boxsize = 100, Thindisk = False):
    print 'Calculating the spectrum of ' + File[:-5]
    print ' ' 
    #---------------------------------------------------------------------------
    #Getting the data from the fits file.
    #---------------------------------------------------------------------------
    #Open the fits file. 
    hdu = fits.open(Directory + File)[0]
    #Putting the data and header in seperate files.
    header = hdu.header
    data = hdu.data

    if Thindisk:
        Data = data
    else:
        Data = data[0,:,:,:]

    #Checking for nan values and putting them on zero.
    Data[np.isnan(Data)] = 0 

    #Determining the shape of the data.
    Shape_Data = np.shape(Data)

    #The shapes of the axis in the data. 
    n_y = Shape_Data[1]
    n_x = Shape_Data[2]

    #Determining the begin_v, delta_v and the amount of images N.
    N = header['NAXIS3']

    if header['CTYPE3']  == 'VELO-LSR':
        begin_v = header['LSTART']
        delta_v = header['LWIDTH']
    if header['CTYPE3']  == 'FREQ':
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
    
    PixelWidth_RA = header['CDELT1']
    PixelWidth_DEC = header['CDELT2']

    #The used units.
    Image_Units = header['BUNIT']

    #First start by slicing the data.
    Data_Sliced = Data[:, (n_y - Boxsize)/2:(n_y + Boxsize)/2, (n_x - Boxsize)/2:(n_x + Boxsize)/2]

    #Then we get the spectrum.
    Spectrum = np.sum(Data_Sliced, axis = (1,2)) #Jy/beam

    #Now we need to know the velocity range.
    Velocities = np.arange(begin_v, begin_v + N*delta_v, delta_v)#km s^-1 

    #Creating the directory for saving the data, if it does not already exists.
    save_directory = Directory + 'Spectrum/' 
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)

    #Now we save the result.
    Save_File = np.zeros((len(Spectrum), 2))
    Save_File[:,0] = Velocities
    Save_File[:,1] = Spectrum

    np.savetxt(save_directory + 'Spectrum_' + Object + '.txt', Save_File, delimiter = ' ')
    #---------------------------------------------------------------------------
    #End function
    #---------------------------------------------------------------------------

'''
This function plots the spectrum calculated in the function Calculate_Spectrum.
File = [string]
Directory = [string]
Object = [string]
Line = [string]
'''
def Plot_Spectrum(File, Directory, Object, Line):
    print 'Plotting the spectrum of ' + File[:-5]
    print ' ' 

    save_directory = Directory + 'Spectrum/'
    Data = np.loadtxt(save_directory + 'Spectrum_' + Object + '.txt', delimiter = ' ')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.step(Data[:,0], Data[:,1])

    ax.set_title('The spectrum of ' + Object +  ' in ' + Line) 

    ax.set_xlabel('Velocity (km s$^{-1}$)')
    ax.set_ylabel('Jy/beam')

    ax.set_xlim([Data[0,0], Data[-1,0]])

    plt.savefig(save_directory + 'Spectrum_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Spectrum_' + Object + '.pdf', bbox_inches='tight')

    #---------------------------------------------------------------------------
    #End function
    #---------------------------------------------------------------------------

'''
This function plots multiple spectra in one figure. The spectra are all 
minimized on their largest value. 
Directories = [string]
Object = [string]
Lines = [string]
Colours = [string]
'''
def Plot_Multiple_Spectrum(Directories, Object, Lines, Colours):
    print 'Plotting the spectra of ' + Object
    print ' '

    #Creating the directory to save the result.
    save_directory = 'Spectra_' + Object + '/'
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)

    #Plotting the result.
    fig = plt.figure()
    ax = fig.add_subplot(111)

    #Adding the spectra.
    for Directory, Line, Colour in zip(Directories, Lines, Colours):
        #Loading and unpacking the data.        
        Data = np.loadtxt(Directory + 'Spectrum/Spectrum_' + Object + '.txt',
                               delimiter = ' ')
        #Unpacking the data.
        Velocities = Data[:,0]
        Spectrum = Data[:,1]

        #Normalizing the spectrum.
        Spectrum_N = Spectrum/np.max(Spectrum)

        #Plotting the line.
        ax.step(Velocities, Spectrum_N, label = Line, color = Colour)
        ax.fill_between(Velocities, 0, Spectrum_N, step = 'pre', alpha = 0.5, color = Colour)
        #End for loop

    ax.legend(loc = 'upper right', numpoints = 1)
    
    ax.set_title('The spectra of ' + Object) 

    ax.set_xlabel('Velocity (km s$^{-1}$)')
    ax.set_ylabel('Normalized intensity')

    ax.set_xlim([-6, 2])

    plt.savefig(save_directory + 'Spectra_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Spectra_' + Object + '.pdf', bbox_inches='tight')


    #---------------------------------------------------------------------------
    #End function
    #---------------------------------------------------------------------------

if __name__ == '__main__':
    Object = 'L1527'
    Line = 'C$^{18}$O'

    File = 'L1527_C18O_1.fits'
    Directory = 'C18O_1/'

    Calculate_Spectrum(File, Directory, Boxsize = 100)
    Plot_Spectrum(File, Directory, Object, Line)

    File = 'L1527_C18O_2.fits'
    Directory = 'C18O_2/'

    Calculate_Spectrum(File, Directory, Boxsize = 200)
    Plot_Spectrum(File, Directory, Object, Line)

