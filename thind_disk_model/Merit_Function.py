import numpy as np
import pyfits as pf
from thindisk_for_optimization import *

def zoom(im,x,y,bb):
    # cut out a square box from image im centered on (x,y) with half-box size bb
    if len(im.shape) == 3:
        return im[:,y-bb:y+bb,x-bb:x+bb]
    elif len(im.shape) == 2:
        return im[y-bb:y+bb,x-bb:x+bb]
    else:
        print 'Wrong number of dimensions! N must be 2 or 3.'

def Chi_Squared(Science, Model):
    '''
    This function calculates the chi squared of the Model.
    Science = [float]
    Model = [float]
    '''    
    
    #The shape of the array.
    n_z, n_y, n_x = Science.shape
    #The size of the science.
    Size = Science.size
    #Calculating the standard deviation of the science data. 
    std = np.std(Science, axis = (1,2))
    #Broadcasting it to 3 dimensions.
    std = std[:, np.newaxis, np.newaxis]
    
    #Calculating chi squared.
    temp = np.square((Science - Model)/std)
    chi_squared = np.sum(temp)/Size
    return chi_squared

def Calculate_Moment_Maps(im, std_im, delta_v, begin_v, mask):
    '''
    This function calculates the moments maps.
    '''
    #The standard deviations of the images. 
    #std_im = np.std(im, axis = (1,2)) 
    #std_im = std_im[:, np.newaxis, np.newaxis]
    
    M0 = np.sum(im*delta_v, axis = 0)
    
    numbers = np.arange(0, len(mask), 1)
    velocities = (numbers[mask]*delta_v + begin_v)
    velocities = velocities[:, np.newaxis, np.newaxis]
    
    im[im < 3*std_im] = 0
    
    M1 = np.sum(im*velocities*delta_v/M0, axis = 0)
    
    temp = im*np.square(velocities - M1)*delta_v
    M2 = np.sum(np.sqrt(temp/np.abs(M0)), axis = 0)
    
    M1[M1 == 0] = np.nan
    M2[M2 == 0] = np.nan
    
    return M0, M1, M2
    
def Merit_FunctionThindisk(Input, i, material_v):
    '''
    This is the merit function for the thindisk optimization.
    mstar = mass of the central object in solar masses.
    pa = position angle of the source in degrees.
    rc = the centrifugal radius in AU.
    size = the size of the total object in AU.
    xcen, ycen = the center of the object in pixels.
    int0 = the highest intensity of the source in JY/BEAM.
    fwhm = the fwhm of the intensity distribution in JY/BEAM.
    '''
    #Interpreting the input.
    if material_v == 'rotating_only_in_rc':
        mstar = Input[0]
        rc = Input[1]
        size = rc
        int0 = Input[2]
        fwhm = Input[3]
    elif material_v == 'fixed_size':
        mstar = Input[0]
        rc = Input[1]
        size = 800#AU
        int0 = Input[2]
        fwhm = Input[3]    
    else:
        mstar = Input[0]
        rc = Input[1]
        size = Input[2]
        int0 = Input[3]
        fwhm = Input[4]
        
    #To prevent negative parameters to be a solution.
    if mstar < 0.0 or mstar > 5:
        return 100000000
    if rc < 0:
        return 100000000
    if size < 0 or size > 2000:
        return 100000000
    if int0 < 0:
        return 100000000
    if fwhm < 0 or fwhm > 5000:
        return 100000000
    if rc > size:
        return 100000000

    #The file and directory.
    File = 'L1527_C18O_2.fits'
    Directory = 'C18O_2/'
    
    #We take the data between channel 1 and 2.
    Channel_1 = 18
    Channel_2 = 33
    
    #Channel_1 = 18
    #Channel_2 = 29
    #And between 3 and 4. 
    Channel_3 = 40
    Channel_4 = 54
    
    #Channel_3 = 44
    #Channel_4 = 53

    #---------------------------------------------------------------------------
    #The width of the zoomed box. 
    Boxsize = 152
    #First we create the thindisk datacube, given the parameters. 
    Thindisk_Datacube = Thindisk(mstar, rc, size, int0, fwhm, Boxsize, Channel_1, Channel_2, Channel_3, Channel_4, i, material_v)
    #Then we load the data 
    Science_Datacube = pf.getdata(Directory + File)[0,:,:,:]

    #Then we zoom in on central part and create a box with sides Boxsize.
    x_cen = Science_Datacube.shape[2]/2 - 1
    y_cen = Science_Datacube.shape[1]/2 - 1

    #Creating the datacube that will have the residues. 
    v_channels = Channel_2 - Channel_1 + 1 + Channel_4 - Channel_3 + 1
    Residue_Datacube = np.zeros((v_channels, Boxsize, Boxsize))

    #Cropping the science data.
    Used_Science_temp = zoom(Science_Datacube, x_cen, y_cen, Boxsize/2)
    Used_Model_temp = Thindisk_Datacube
    
    #The mask to use to get the data.
    mask = np.zeros(Used_Science_temp.shape[0], dtype = bool)
    mask[Channel_1:Channel_2] = True
    mask[Channel_3:Channel_4] = True
    
    #Slicing the datacubes to get only the channels we want. 
    Used_Science = Used_Science_temp[mask,:,:]
    Used_Model = Used_Model_temp[mask,:,:]
    
    result = Chi_Squared(Used_Science, Used_Model)
    '''
    i = 0
    while i < len(Used_Model[:,0,0]):
        import matplotlib.pyplot as plt
        plt.imshow(Used_Model[i,:,:], origin = 'lower')
        plt.colorbar()
        plt.show()
        i+=1
    '''
    #-------------------------------------------------------------------
    #Showing the results.
    #-------------------------------------------------------------------
    Show_Results = False
    
    if Show_Results == True:
        delta_v = 0.159999957279
        begin_v = 0
        std_im_science = np.std(Science_Datacube, axis = (1,2))
        std_im_science = std_im_science[:, np.newaxis, np.newaxis]
        std_im_model = np.std(Used_Model, axis = (1,2))
        std_im_model = std_im_model[:, np.newaxis, np.newaxis] 
        M0_science, M1_science, M2_science = Calculate_Moment_Maps(Used_Science, std_im_science[mask], delta_v, begin_v, mask)
        M0_model, M1_model, M2_model = Calculate_Moment_Maps(Used_Model, std_im_model, delta_v, begin_v, mask)
        
        import matplotlib.pyplot as plt
        
        fig = plt.figure()
        
        ax_1 = fig.add_subplot(321)  
        cax_1 = ax_1.imshow(M0_science, origin = 'lower')
        cbar_1 = fig.colorbar(cax_1)
        ax_1.set_title('M0 Science')
        
        ax_2 = fig.add_subplot(322)  
        cax_2 = ax_2.imshow(M0_model, origin = 'lower')
        cbar_2 = fig.colorbar(cax_2)
        ax_2.set_title('M0 Model')
        
        ax_3 = fig.add_subplot(323)  
        cax_3 = ax_3.imshow(M1_science, origin = 'lower', vmin = 1, vmax = 12.0)
        cbar_3 = fig.colorbar(cax_3)
        ax_3.set_title('M1 Science')
        
        ax_4 = fig.add_subplot(324)  
        cax_4 = ax_4.imshow(M1_model, origin = 'lower', vmin = 1, vmax = 12.0)
        cbar_4 = fig.colorbar(cax_4)
        ax_4.set_title('M1 Model')
        
        ax_5 = fig.add_subplot(325)  
        cax_5 = ax_5.imshow(M2_science, origin = 'lower', vmin = -0.5, vmax = 4)
        cbar_5 = fig.colorbar(cax_5)
        ax_5.set_title('M2 Science')
        
        ax_6 = fig.add_subplot(326)  
        cax_6 = ax_6.imshow(M2_model, origin = 'lower', vmin = -0.5, vmax = 4)
        cbar_6 = fig.colorbar(cax_6)
        ax_6.set_title('M2 Model')
        
        #Showing the images for 5 seconds and then closing them again.
        plt.draw()
        plt.pause(5)        
        plt.close()
        
    #Printing the results.
    print  '{0:9s}  {1:9s}  {2:9s}   {3:9s}  {4:9s}  {5:9s}'.format('mstar', ' rc', ' size', ' int0' , ' fwhm', ' func(x)')
    print '{0: 3.3f}   {1: 3.3f}   {2: 3.3f}   {3: 3.3f}    {4: 3.3f}   {5: 3.3f}'.format(mstar, rc, size, int0, fwhm, result)

    return result
    
if __name__ == '__main__':
    print 'testing'
    
    
'''
------------------------------------------------------------------------
Unused code
------------------------------------------------------------------------

#print 'Calculating the difference between the thindisk model and the science data.'
#Calculating the difference between the science data and the thindisk model.
i = Channel_1
while i < Channel_2 + 1:
    Used_Science = zoom(Science_Datacube[i,:,:], x_cen, y_cen, Boxsize/2)
    Used_Model = Thindisk_Datacube[i,:,:]
    Residue_Datacube[i - Channel_1,:,:] = Used_Science - Used_Model
    i += 1
    #End while

i = Channel_3
while i < Channel_4 + 1:
    Used_Science = zoom(Science_Datacube[i,:,:], x_cen, y_cen, Boxsize/2)
    Used_Model = Thindisk_Datacube[i,:,:]
    Residue_Datacube[i - Channel_3 + Channel_2 - Channel_1 + 1,:,:] = Used_Science - Used_Model
    i += 1
    #End while

result = np.sum(np.abs(Residue_Datacube))
'''
