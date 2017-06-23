import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit
mpl.rc('text', usetex=True)

#Make the plots nice
import Nice_Plots_2
Nice_Plots_2.set_style()

def list_files(path):
    # Returns a list of names (with extension, without full path) of all  
    # files in folder path.
    files = []
    for name in os.listdir(path):
        if os.path.isfile(os.path.join(path, name)):
            files.append(name)
    return files 
    
def Read_Data(List, path):
    '''
    This function opens all the files in the list List in the directory
    path. Then it put all the starting parameters and result parameters
    together in one numpy array and returns that.
    '''
    #The arrays where the data is saved.
    X0 = np.zeros((len(List),4))
    X = np.zeros((len(List),4))
    for File, Row in zip(List, np.arange(0, len(List) + 1, 1)):
        #First we load the file.
        Data = np.loadtxt(path + File, delimiter = ' ')
        #Writing to the arrays.
        X0[Row,:] = Data[:, 0]
        X[Row,:] = Data[:, 1]
            
    return X0, X
    
def Find_Best_Fit(List, path, print_output = False):
    '''
    This function gets the list with files and looks at all the file 
    names to determine the best fit (lowest function value). It then 
    prints the file name and returns the resulting parameters. 
    '''
    func_values = np.zeros(len(List))
    
    #First we need to get all the function values.
    for File, i in zip(List, np.arange(0, len(List) + 1, 1)):
        func_values[i] = File[30:-36]#float(File[41:-36])
    
    #Determing what and where the best fit is. 
    func_min = np.min(func_values)
    func_min_pos = np.where(func_values == np.min(func_values))[0][0]
    File_min = List[func_min_pos]
        
    #Information about the file.
    Angle = File[2:4]
    
    #Then we load the best file.
    Data = np.loadtxt(path + File_min, delimiter = ' ')
    bestx = Data[:, 1]
    
    if print_output:
        print '------------------------------------------------------------'
        print 'For inclination: ' + Angle + ' degrees'
        print 'The best file was: '
        print File_min
        print 'With a function value of: ' + str(func_min)
        print 'The parameters are: '    
        print  '{0:9s}  {1:9s}  {2:9s}   {3:9s}  {4:9s}'.format('mstar', ' rc', ' size', ' int0' , ' fwhm')
        print '{0: 3.3f}   {1: 3.3f}   {2: 3.3f}   {3: 3.3f}    {4: 3.3f}'.format(bestx[0], bestx[1], bestx[1], bestx[2], bestx[3])
        print '------------------------------------------------------------'
    
    #Saving the best solution.
    #Creating the directory for saving the images, if it does not already exists.
    save_directory = path + '../Best_Results/' 
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
    #Checking (if there is an file) if the best fit has been changed. 
    if os.path.isfile(save_directory + 'Best_Parameters' + '_i_' + Angle + '.txt'):
        Last = np.loadtxt(save_directory + 'Best_Parameters' + '_i_' + Angle + '.txt', delimiter = ' ')
        print ' '
        print 'For inclination: ' + Angle + ' degrees'
        
        Is_Close = np.isclose(Data, Last)*np.ones_like(Data)
        if np.sum(Is_Close) == Data.size:
            
            print 'The best fit did not change.'
        else:
            print 'The best fit changed!'

    #Saving the best fit.
    np.savetxt(save_directory + 'Best_Parameters' + '_i_' + Angle + '.txt', Data, delimiter = ' ', header = 'i = '  + Angle + '\nfunc = ' + str(func_min) + '\nFile: ' +  File_min + '\np0                       x')
    
    return bestx 

def Find_Most_Common_Fit(X, path, Angle):
    '''
    This function gets the list with files and finds the most common 
    fits and returns the parameters. 
    '''
    
    #Unpacking the data.
    M = X[:,0]
    rc = X[:,1]
    int0 = X[:,2]
    fwhm = X[:,3]    
  
    #Now we make two datasets. One below rc_select and one above.
    rc_select = 200#AU
    
    #Lower radii.
    rc1 = rc[rc < rc_select]
    M1 = M[rc < rc_select]
    int01 = int0[rc < rc_select]
    fwhm1 = fwhm[rc < rc_select]
        
    #Upper radii.
    rc2 = rc[rc > rc_select]
    M2 = M[rc > rc_select]
    int02 = int0[rc > rc_select]
    fwhm2 = fwhm[rc > rc_select]
    
    #We bin the radii to select the most common ones. 
    N_bin = 50
    
    rc_min_1 = 0
    rc_max_1 = 150
    binBoundaries_1 = np.linspace(rc_min_1,rc_max_1,N_bin + 1)
    
    #-------------------------------------------------------------------
    save_directory_fig = path + '../Result_Plots/' 
    if not os.path.exists(save_directory_fig):
        os.makedirs(save_directory_fig)
    
    
    #Lets plots these distributions.
    fig = plt.figure()
    fig.subplots_adjust(hspace=.5, wspace = .15)
    
    #The amount of fits. 
    N1 = len(rc1)
    
    #Defining the axes.
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
 
    
    ticksize = 10
       
    #Plotting the mass distribution.
    Binned_M1 = ax1.hist(M1, bins = 50)
    ax1.set_title('Mass')
    ax1.set_xlabel('M [M$_{\odot}$]')
    ax1.tick_params(axis='both', which='major', labelsize=ticksize)
    text = 'Most common results, rc $<$ ' + str(rc_select) + ' AU'
    ax1.text(0.7, 1.2,  text, transform=ax1.transAxes)
    
    #Plotting the rc distribution.
    Binned_rc1 = ax2.hist(rc1, bins = binBoundaries_1)
    ax2.set_title('Centrifugal Radius')
    ax2.set_xlabel('rc [AU]')   
    ax2.tick_params(axis='both', which='major', labelsize=ticksize)
    text = 'N = ' + str(N1)
    ax2.text(1., -0.22,  text, transform=ax2.transAxes)
    
    #Plotting the I0 distribution.
    Binned_int01 = ax3.hist(int01, bins = 50)
    ax3.set_title('Intensity')
    ax3.set_xlabel('I$_0$ [Jy/beam]')
    ax3.tick_params(axis='both', which='major', labelsize=ticksize)
    
    #Plotting the fwhm distribution.
    Binned_fwhm1 = ax4.hist(fwhm1, bins = 50)
    ax4.set_title('Fwhm')
    ax4.set_xlabel('fwhm [arcsec]') 
    ax4.tick_params(axis='both', which='major', labelsize=ticksize)
    
    plt.savefig(save_directory_fig + 'Most_Common_Results_Below_' + str(rc_select) + '_AU_i_' + Angle + '.jpg')
    plt.savefig(save_directory_fig + 'Most_Common_Results_Below_' + str(rc_select) + '_AU_i_' + Angle + '.pdf')
    plt.close()
    #-------------------------------------------------------------------
    #Lets plots these distributions.
    fig = plt.figure()
    fig.subplots_adjust(hspace=.75, wspace = .15)
    
    #The amount of fits. 
    N2 = len(rc2)
    
    #Defining the axes.
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    
    def Gaussian(x, A, x_center, sigma):
        return A*np.exp(-np.square(x - x_center)/(2*np.square(sigma)))
 
    N = 50
    
    #Plotting the mass distribution.
    Binned_M2 = ax1.hist(M2[(M2 < 1.35)*(M2 > 1.25)], bins = N, alpha=0.5, color = 'b')
    
    #Fitting a Gaussian to the bins.
    y1= Binned_M2[0]
    x1 = Binned_M2[1][0:N] + (Binned_M2[1][1] - Binned_M2[1][0])/2
    M_param = curve_fit(Gaussian, x1, y1, p0 = np.array([20, 1.28, 0.1]))[0]
    #Then create an array for plotting
    Gaussian_Fit_M = Gaussian(Binned_M2[1], M_param[0], M_param[1], M_param[2])
    #Then we plot the line.
    ax1.plot(Binned_M2[1], Gaussian_Fit_M, color = 'r')
    text = 'M $= ' + str(np.round(M_param[1],2)) + '\pm' + str(np.round(M_param[2],2))  + '$ M$_{\odot}$' 
    ax1.text(0, -0.4,  text, transform=ax1.transAxes)

    
    ax1.set_title('Mass')
    ax1.set_xlabel('M [M$_{\odot}$]')
    ax1.tick_params(axis='both', which='major', labelsize=ticksize)
    text = 'Most common results, rc $>$ ' + str(rc_select) 
    ax1.text(0.7, 1.2,  text, transform=ax1.transAxes)
    
    #Plotting the rc distribution.
    Binned_rc2 = ax2.hist(rc2[(rc2 < 600)*(rc2 > 200)], bins = N, alpha=0.5, color = 'b')#binBoundaries_2)
    
    #Fitting a Gaussian to the bins.
    y1 = Binned_rc2[0]
    x1 = Binned_rc2[1][0:N] + (Binned_rc2[1][1] - Binned_rc2[1][0])/2
    rc_param = curve_fit(Gaussian, x1, y1, p0 = np.array([50, 340, 50]))[0]
    #Then create an array for plotting
    Gaussian_Fit_rc = Gaussian(Binned_rc2[1], rc_param[0], rc_param[1], rc_param[2])
    #Then we plot the line.
    ax2.plot(Binned_rc2[1], Gaussian_Fit_rc, color = 'r')
    text = 'rc $= ' + str(np.round(rc_param[1],2)) + '\pm' + str(np.round(rc_param[2],2))  + '$ AU' 
    ax2.text(0, -0.4,  text, transform=ax2.transAxes)
    
    ax2.set_title('Centrifugal Radius')
    ax2.set_xlabel('rc [AU]')   
    ax2.tick_params(axis='both', which='major', labelsize=ticksize)
    text = 'N = ' + str(N2)
    ax2.text(.9, -0.3,  text, transform=ax2.transAxes)
        
    #Plotting the I0 distribution.
    Binned_int02 = ax3.hist(int02[(int02 < 0.5)*(int02 > 0)], bins = N, alpha=0.5, color = 'b')
    
    #Fitting a Gaussian to the bins.
    y1 = Binned_int02[0]
    x1 = Binned_int02[1][0:N] + (Binned_int02[1][1] - Binned_int02[1][0])/2
    int0_param = curve_fit(Gaussian, x1, y1, p0 = np.array([40, 0.23, 0.05]))[0]
    #Then create an array for plotting
    Gaussian_Fit_int0 = Gaussian(Binned_int02[1], int0_param[0], int0_param[1], int0_param[2])
    #Then we plot the line.
    ax3.plot(Binned_int02[1], Gaussian_Fit_int0, color = 'r')
    text = 'int0 $= ' + str(np.round(int0_param[1],2)) + '\pm' + str(np.round(int0_param[2],2))  + '$ Jy/Beam' 
    ax3.text(0, 1.2,  text, transform=ax3.transAxes)
    
    ax3.set_title('Intensity')
    ax3.set_xlabel('I$_0$ [Jy/beam]')
    ax3.tick_params(axis='both', which='major', labelsize=ticksize)
    MC_int02 = Binned_int02[1][np.where(Binned_int02[0] == np.max(Binned_int02[0]))[0][0]]
    
    #Plotting the fwhm distribution.
    Binned_fwhm2 = ax4.hist(fwhm2[(fwhm2 < 2)*(fwhm2 > 0.5)], bins = N, alpha=0.5, color = 'b')
    
    #Fitting a Gaussian to the bins.
    y1 = Binned_fwhm2[0]
    x1 = Binned_fwhm2[1][0:N] + (Binned_fwhm2[1][1] - Binned_fwhm2[1][0])/2
    fwhm_param = curve_fit(Gaussian, x1, y1, p0 = np.array([50, 1.2, 0.2]))[0]
    #Then create an array for plotting
    Gaussian_Fit_fwhm = Gaussian(Binned_fwhm2[1], fwhm_param[0], fwhm_param[1], fwhm_param[2])
    #Then we plot the line.
    ax4.plot(Binned_fwhm2[1], Gaussian_Fit_fwhm, color = 'r')
    text = 'fwhm $= ' + str(np.round(fwhm_param[1],2)) + '\pm' + str(np.round(fwhm_param[2],2))  + '"$ ' 
    ax4.text(0, 1.2,  text, transform=ax4.transAxes)
    
    ax4.set_title('Fwhm')
    ax4.set_xlabel('fwhm [arcsec]') 
    ax4.tick_params(axis='both', which='major', labelsize=ticksize)
    MC_fwhm2 = Binned_fwhm2[1][np.where(Binned_fwhm2[0] == np.max(Binned_fwhm2[0]))[0][0]]
    
    #Text with the parameters
    text_param = 'M $= ' + str(np.round(M_param[1], 3)) + '\pm' + str(np.abs(np.round(M_param[2], 3))) + '$ M$_{\odot}$ \n rc $= ' + str(np.round(rc_param[1], 2)) + '\pm' + str(np.abs(np.round(rc_param[2], 2))) + '$ AU \n I$_0$ $= ' + str(np.round(MC_int02 , 2)) + '$ Jy/beam \n fwhm $= ' + str(np.round(MC_fwhm2,2)) + '$ arcsec'
    ax3.text(0, -1.5,  text_param, transform=ax3.transAxes)
    
    plt.savefig(save_directory_fig + 'Most_Common_Results_Above_' + str(rc_select) + '_AU_i_' + Angle + '.jpg')
    plt.savefig(save_directory_fig + 'Most_Common_Results_Above_' + str(rc_select) + '_AU_i_' + Angle + '.pdf')
    plt.savefig(save_directory_fig + 'Most_Common_Results_Above_' + str(rc_select) + '_AU_i_' + Angle + '.eps')
    plt.savefig(save_directory_fig + 'Most_Common_Results_Above_' + str(rc_select) + '_AU_i_' + Angle + '.ps')
    plt.close()
    
    #-------------------------------------------------------------------
    #Lets make the plot used in the paper. 
    fig = plt.figure(figsize = (12,5))
    #fig.subplots_adjust(wspace = .15)
    
    #Which datapoints we plot, only realistic ones and the main peaks.
    select = (rc2 < 600)*(rc2 > 200)*(M2 < 1.35)*(M2 > 1.25)
    N2 = len(M2[select])
    print 'The amount of datapoints in this subset is: ' + str(N2)
    
    #Defining the axes.
    ax1 = fig.add_subplot(1,2,1)#, aspect='auto')
    ax2 = fig.add_subplot(1,2,2)#, aspect='auto')
 
    N = 30
 
    #Plotting the mass distribution.
    Binned_M2 = ax1.hist(M2[select], bins = N, alpha=0.5, color = 'b')
    ax1.set_xlabel('M [M$_{\odot}$]')
    ax1.set_ylabel('$\#$')
    ax1.tick_params(axis='both', which='major', labelsize=ticksize)
   
    
    #Plotting the rc distribution.
    Binned_rc2 = ax2.hist(rc2[select], bins = N, alpha=0.5, color = 'b')   
    ax2.set_xlabel('$R_d$ [AU]')  
    ax2.tick_params(axis='both', which='major', labelsize=ticksize)
    ax2.set_xlim([250,600])
         
    plt.savefig(save_directory_fig + 'Distribution_Plot_Paper.jpg')
    plt.savefig(save_directory_fig + 'Distribution_Plot_Paper.pdf')
    plt.savefig(save_directory_fig + 'Distribution_Plot_Paper.eps')
    plt.savefig(save_directory_fig + 'Distribution_Plot_Paper.ps')
    plt.close()
    
    fig = plt.figure(figsize = (6,5))
    ax = fig.add_subplot(1,1,1)#, aspect= 1./1200)
    ax.hist(M2[select], bins = N, alpha=0.5, color = 'b')
    ax.set_xlabel('M [M$_{\odot}$]')#, labelpad=0.6) 
    ax.set_ylabel('$\#$')
    ax.tick_params(axis='both', which='major')#, labelsize=ticksize)
    plt.savefig(save_directory_fig + 'Distribution_Plot_Mass.jpg',bbox_inches='tight')
    plt.savefig(save_directory_fig + 'Distribution_Plot_Mass.pdf',bbox_inches='tight')
    plt.close()
    
    fig = plt.figure(figsize = (6,5))
    ax = fig.add_subplot(1,1,1)#, aspect='equal')
    ax.hist(rc2[select], bins = N, alpha=0.5, color = 'b')
    ax.set_xlabel('$R_d$ [AU]')#, labelpad=20)  
    ax.tick_params(axis='both', which='major')#, labelsize=ticksize)
    ax.set_xlim([250,600])
    ax.set_ylabel('$\#$')
    plt.savefig(save_directory_fig + 'Distribution_Plot_Rd.jpg',bbox_inches='tight')
    plt.savefig(save_directory_fig + 'Distribution_Plot_Rd.pdf',bbox_inches='tight')
    plt.close()
    
    
    
    
    #-------------------------------------------------------------------
    
    MC_M1 = Binned_M1[1][np.where(Binned_M1[0] == np.max(Binned_M1[0]))[0][0]]
        
    MC_rc1 = Binned_rc1[1][np.where(Binned_rc1[0] == np.max(Binned_rc1[0]))[0][0]]
                
    MC_int01 = Binned_int01[1][np.where(Binned_int01[0] == np.max(Binned_int01[0]))[0][0]]
        
    MC_fwhm1 = Binned_fwhm1[1][np.where(Binned_fwhm1[0] == np.max(Binned_fwhm1[0]))[0][0]]
    
        
    #Combing this into one array.ticksize = 10
    MC_1 = np.array([MC_M1, MC_rc1, MC_int01, MC_fwhm1], dtype = float)
    MC_2 = np.array([M_param[1], rc_param[1], MC_int02, MC_fwhm2], dtype = float)
    
    #Saving the best solution.
    #Creating the directory for saving the images, if it does not already exists.
    save_directory = path + '../Most_Common_Results/' 
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
    
    save_name_1 = 'Most_Common_Parameters_1_' + '_i_' + Angle + '.txt'
    save_name_2 = 'Most_Common_Parameters_2_' + '_i_' + Angle + '.txt'    

    #Saving the best fit.
    np.savetxt(save_directory + save_name_1, MC_1, delimiter = ' ', header = 'i = '  + Angle + '\nrc < ' + str(rc_select) + ' AU')
    np.savetxt(save_directory + save_name_2, MC_2, delimiter = ' ', header = 'i = '  + Angle + '\nrc > ' + str(rc_select) + ' AU')
    
    #Im not sure if this part works....
    def Check_File(directory, File, Values):
        '''
        This function checks if the result has been changed. 
        '''
        #Checking (if there is an file) if the best fit has been changed. 
        if os.path.isfile(directory + File):
            Last = np.loadtxt(directory + File, delimiter = ' ')
            
            Is_Close = np.isclose(Values, Last)*np.ones_like(Values)
            if np.sum(Is_Close) == Values.size:
                print 'The most common fit did not change.'
            else:
                print 'The most common fit changed!'
                
    print ' '
    print 'For inclination: ' + Angle + ' degrees and rc < ' + str(rc_select) + ' AU'
    Check_File(save_directory, save_name_1, MC_1)  
    print ' '
    print 'For inclination: ' + Angle + ' degrees and rc > ' + str(rc_select) + ' AU'
    Check_File(save_directory, save_name_2, MC_2)    
    
    return MC_1, MC_2
    
def Update_INI(X, Path, Angle, Type):
    from ConfigParser import SafeConfigParser
    import sys
    '''
    This function modifies the .ini file with the new best parameters so
    that a new datacube can be created. And it creates the directory 
    where the datacube will be saved.
    '''
    if Type == 'Best':
        #Name = 'Optimization_Best_' + str(Angle) + '_2'
        Name = 'Optimization_' + str(Angle)
        save_directory = Path + '../Model_Best_' + str(Angle) + '/' 
    elif Type == 'MC1':
        Name = 'Optimization_Most_Common_1_' + str(Angle) + '_2' 
        save_directory = Path + '../Model_Most_Common_1_' + str(Angle) + '/' 
    elif Type == 'MC2':
        Name = 'Optimization_Most_Common_2_' + str(Angle) + '_2' 
        save_directory = Path + '../Model_Most_Common_2_' + str(Angle) + '/' 
    else:
        print 'I do not recognize the type.'
        
    #Creating the directory to save the datacube.
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
      
    #Opening the .ini file.
    parser = SafeConfigParser()
    parser.read(Name + '.ini')
    
    #Changing the parameters.
    parser.set('disk', 'mstar', str(X[0]))
    parser.set('disk', 'rc', str(X[1]))
    parser.set('disk', 'size', str(800))
    parser.set('line', 'intensity', 'gaussian, ' + str(X[2]) + ', ' + str(X[3]))
    parser.set('output', 'name', save_directory + Name)
      
    #Saving the .ini file.
    with open(Name + '.ini', 'w') as configfile:
        parser.write(configfile)

def Plot_Parameter(X0, X, Xbest, Path, Angle, N):
    '''
    This function plots the histograms of the parameters in X and X0. 
    '''

    #Unpacking the data.
    M = X[:,0]
    M0 = X0[:,0]
    Mbest = Xbest[0]
        
    rc = X[:,1]
    rc0 = X0[:,1]
    rcbest = Xbest[1]
    
    size = X[:,1]
    size0 = X0[:,1]
    sizebest = Xbest[1]
    
    int0 = X[:,2]
    int00 = X0[:,2]
    int0best = Xbest[2]
    
    fwhm = X[:,3]
    fwhm0 = X0[:,3]
    fwhmbest = Xbest[3]
    #-------------------------------------------------------------------
    
    print 'Best M = ' + str(Mbest)
    print 'Best rc = ' + str(rcbest)
    print 'Best int0 = ' + str(int0best)
    print 'Best fwhm = ' + str(fwhmbest) 
    
    xwidth = 6.6
    ywidth = 6.6
    
    fig = plt.figure(figsize=(xwidth,ywidth))
    fig.subplots_adjust(hspace=.65, wspace = .35)
        
    #Defining the axes.
    ax1 = fig.add_subplot(3,2,1)
    ax2 = fig.add_subplot(3,2,2)
    ax3 = fig.add_subplot(3,2,3)
    ax4 = fig.add_subplot(3,2,4)
    ax5 = fig.add_subplot(3,2,5)
    
    N_bin = 50
    ticksize = 10
    
    #The mass.
    M_min = 0
    M_max = 0.5
    binBoundaries = np.linspace(M_min,M_max,N_bin + 1)
    ax1.hist(M0, bins = binBoundaries, histtype='stepfilled', color='b', alpha=0.5, label='Begin Parameters')
    ax1.hist(M, bins = binBoundaries, histtype='stepfilled', color='r', alpha=0.5, label='Result Parameters')
    if Mbest < M_max and Mbest > M_min:
        ax1.axvline(Mbest, label = 'Best Parameter')
    ax1.set_title('Mass')
    ax1.set_xlabel('M [M$_{\odot}$]')
    ax1.set_ylabel('$\#$')
    ax1.tick_params(axis='both', which='major', labelsize=ticksize)
    
    #The rc.
    rc_min = 0
    rc_max = 100
    binBoundaries = np.linspace(rc_min,rc_max,N_bin + 1)
    ax2.hist(rc0, bins = binBoundaries, histtype='stepfilled', color='b', alpha=0.5, label='Begin Parameters')
    ax2.hist(rc, bins = binBoundaries, histtype='stepfilled', color='r', alpha=0.5, label='Result Parameters')
    if rcbest < rc_max and rcbest > rc_min:
        ax2.axvline(rcbest, label = 'Best Parameter')
    ax2.set_title('Centrifugal Radius')
    ax2.set_xlabel('rc [AU]')
    ax2.set_ylabel('$\#$')  
    ax2.tick_params(axis='both', which='major', labelsize=ticksize)
    
    #The size.
    size_min = 0
    size_max = 800
    binBoundaries = np.linspace(size_min,size_max,N_bin + 1)
    ax3.hist(size0, bins = binBoundaries, histtype='stepfilled', color='b', alpha=0.5, label='Begin Parameters')
    ax3.hist(size, bins = binBoundaries, histtype='stepfilled', color='r', alpha=0.5, label='Result Parameters')
    if sizebest < size_max and sizebest > size_min:
        ax3.axvline(sizebest, label = 'Best Parameter')
    ax3.set_title('Size')
    ax3.set_xlabel('size [AU]')
    ax3.set_ylabel('$\#$') 
    ax3.tick_params(axis='both', which='major', labelsize=ticksize)
    
    #The int0.
    int0_min = 0
    int0_max = 15
    binBoundaries = np.linspace(int0_min, int0_max,N_bin + 1)
    ax4.hist(int00, bins = binBoundaries, histtype='stepfilled', color='b', alpha=0.5, label='Begin Parameters')
    ax4.hist(int0, bins = binBoundaries, histtype='stepfilled', color='r', alpha=0.5, label='Result Parameters')
    if int0best < int0_max and int0best > int0_min:
        ax4.axvline(int0best, label = 'Best Parameter')
    ax4.set_title('Intensity')
    ax4.set_xlabel('I$_0$ [Jy/beam]')
    ax4.set_ylabel('$\#$')
    ax4.tick_params(axis='both', which='major', labelsize=ticksize)
    text = '$i = ' + Angle + '^{\circ}$, N = ' + str(N) + '\n \n Best parameters: \n M = ' + str(np.round(Mbest,2)) + ' M$_{\odot}$ \n rc = ' + str(np.round(rcbest,1)) + ' AU \n Size = ' + str(np.round(sizebest)) + ' AU \n I$_0$ = ' + str(np.round(int0best,2)) + ' Jy/beam \n fwhm = ' + str(np.round(fwhmbest, 2)) + ' arcsec'  
    ax4.text(1.1*int0_max, 0, text)
    #text = '$i = ' + Angle + '^{\circ}$, N = ' + str(N) + '\n Best parameters:\n M = ' + str(np.round(Mbest,2)) + ' M$_{\odot}$ \n rc = ' + str(np.round(rcbest,1)) + ' AU \n Size = ' + str(np.round(sizebest)) + ' AU \n I$_0$ = ' + str(np.round(int0best,2)) + ' Jy/beam \n fwhm = ' + str(np.round(fwhmbest, 2)) + ' arcsec'  
    #ax4.text(1.1*int0_max, 0, text)
    
    #The fwhm.
    fwhm_min = 0
    fwhm_max = 3
    binBoundaries = np.linspace(fwhm_min, fwhm_max,N_bin + 1)
    ax5.hist(fwhm0, bins = binBoundaries, histtype='stepfilled', color='b', alpha=0.5, label='Begin Parameters')
    temp = ax5.hist(fwhm, bins = binBoundaries, histtype='stepfilled', color='r', alpha=0.5, label='Result Parameters')
    if fwhmbest < fwhm_max and fwhmbest > fwhm_min:
        ax5.axvline(fwhmbest, label = 'Best Parameter')
    ax5.set_title('Fwhm')
    ax5.set_xlabel('fwhm [arcsec]')
    ax5.set_ylabel('$\#$')
    ax5.legend(loc='center right', bbox_to_anchor=(2.3, 0.35))
    ax5.tick_params(axis='both', which='major', labelsize=ticksize)
    
    #Creating the directory for saving the images, if it does not already exists.
    save_directory = Path + '../Result_Plots/' 
    if not os.path.exists(save_directory):
        os.makedirs(save_directory)
    
    plt.savefig(save_directory + 'Parameters_' + Angle + '.jpg', bbox_inches='tight')
    #plt.savefig(save_directory + 'Parameters_' + Angle + '.pdf', bbox_inches='tight')
    plt.close()
   
def Map_Results(directory):
    #Getting the list met files. 
    Files = list_files(directory)
    Len_Files = len(list_files(directory))
    
    angle = '85'
    
    #The lists where we save the files of the optimizations.
    List_52 = []
    
    #First we split all the files in different lists to have easy file 
    #management. 
    for File in Files:
        Angle = File[2:4]
        if Angle == angle:
            List_52.append(File)
        else:
            print 'I do not recognize the angle.'
    
    #The amount of fits for all the lines and inclinations.        
    N_List_52 = len(List_52)

    
    #Lets read in the data. 
    List_52_X0, List_52_X = Read_Data(List_52, directory)
 
    #Lets find the best fit.
    List_52_best = Find_Best_Fit(List_52, directory, print_output = False)
    
    #Updating the ini file.
    Update_INI(List_52_best, directory, angle, 'Best')
    
    #Lets find the most common fits. 
    #xMC_1_52, MC_2_52 = Find_Most_Common_Fit(List_52_X, directory, angle)
    
    #Update_INI(MC_2_52, directory, angle, 'MC2')
    
    #lets plot the results. 
    Plot_Parameter(List_52_X0, List_52_X, List_52_best, directory, angle, N_List_52)

def Find_Unrealistic_Fits(Path):
    '''
    This function checks all the results for masses or radii that are 
    unrealistic. So if the mass is above 5 M sun or the rc is larger 
    than 5000 AU. The files that meet these criteria are then moved to 
    another directory.
    '''
    print ' '
    print 'Checking the results.'
    #The criteria: 
    mass_limit = 2#Solar masses.
    rc_limit = 2000#AU
    fwhm_limit = 20#Arcsec
    
    #Creating the directory where the files will be moved to. 
    move_directory = Path + 'Unrealistic_Results/' 
    if not os.path.exists(move_directory):
        os.makedirs(move_directory)
        
    #Finding all the files in the directory.
    file_list = list_files(Path)
    
    #Now we check every file. 
    for File in file_list:
        Data = np.loadtxt(Path + File, delimiter = ' ')
        #Checking if the file has to be moved. 
        if Data[0,1] > mass_limit or Data[1,1] > rc_limit or Data[2,1] > fwhm_limit:
            os.rename(Path + File, move_directory + File)
        
if __name__ == '__main__':
    #The directory containing all the files. 
    Object = 'L1527'
    Line = 'C$^{18}$O (Thindisk model)'
    
    Directories = ['Optimization_Results/Round_7/Data/']
    Angle = '85'

    for directory in Directories: 
        print 'Directory: ' + directory
        Check_Results = False
        if Check_Results:
            Find_Unrealistic_Fits(directory)
     
        Map = True
        if Map:
            print 'Mapping the results.'
            Map_Results(directory)
        
        Create_Datacube = False     
        if Create_Datacube:
            from thindisk import *
            #Now using the updated .ini file we automatically create the 
            #datacube.
            print 'Creating the datacubes.'
            main('Optimization_' + Angle + '.ini')
            '''
            main('Optimization_Best_' + Angle + '_2.ini')
            main('Optimization_Most_Common_2_' + Angle + '_2.ini')
            '''
            
        Create_Moment_Maps = False
        if Create_Moment_Maps:
            from Image_Moments import *
            #Then we make the moment maps etc..
            print 'Creating the moment maps.' 
            
            begin_image = 0
            end_image = 29

            File = 'Optimization_Best_' + Angle + '_2_convolved_with_beam.fits'
            Directory = directory + '../Model_Best_' + Angle + '/'
            Image_Moments(File, Directory, begin_image, end_image, 'V883 Ori', Line, Zoom = True, boxsize = 100, vmin_M1 = 1.6, vmax_M1 = 6.4, vmin_M2 = 0.3, vmax_M2 = 1.8, Thindisk = True)
            '''
            File = 'Optimization_Most_Common_2_' + Angle + '_2_convolved_with_beam.fits'
            Directory = directory + '../Model_Most_Common_2_' + Angle + '/'
            Image_Moments(File, Directory, begin_image, end_image, 'V883 Ori', Line, Zoom = True, boxsize = 100, vmin_M1 = 1.6, vmax_M1 = 6.4, vmin_M2 = 0.3, vmax_M2 = 1.8, Thindisk = True)
            '''
            
        Create_PV_Diagrams = False           
        if Create_PV_Diagrams:
            from PV_diagram import *
            print 'Creating the PV diagrams.'
            
            Angle = '85'
            Object = 'L1527'
            Line = 'C$^{18}$O (Thindisk model)'
            
            rotation = -90 #degrees
            width = 200#pixels
            
            File = 'Optimization_85_convolved_with_beam.fits'
            Directory = directory + '../Model_Best_' + Angle + '/'
            PV_diagram(File, Directory, rotation, width, Object, Line, inclination = 85., Velocity_Curve = True, mass = 0.25, v_source = 5.9, d_source = 140, Thindisk = True, Zoom = True, v_width = 8, arcsec_width = 10, Overlay_Contour = 'C18O')
    
        Channels_Comparison = False
        if Channels_Comparison:
            from Channel_Comparison import *
            Angle = '85'
            Object = 'L1527'
            
            Line1 = 'C$^{18}$O'
            Line2 = 'C$^{18}$O (Thindisk model)'
            
            Directory1 = 'C18O_2/'
            Directory2 = directory + '../Model_Best_' + Angle + '/'
            
            File1 = 'L1527_C18O_2.fits'
            File2 = 'Optimization_' + Angle + '_convolved_with_beam.fits'
        
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
            
