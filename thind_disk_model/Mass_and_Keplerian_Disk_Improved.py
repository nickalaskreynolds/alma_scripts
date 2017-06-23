import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from astropy import units as u
from astropy import constants as const
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FormatStrFormatter
#-------------------------------------------------------------------------------
#Save directory
#-------------------------------------------------------------------------------
save_directory = 'Mass_and_Keplerian_Disk/'

#Creating the directory for saving the images, if it does not already exists.
if not os.path.exists(save_directory):
    os.makedirs(save_directory)

#Make the plots nice
import Nice_Plots_2
Nice_Plots_2.set_style()

#-------------------------------------------------------------------------------
#Fitting functions.
#-------------------------------------------------------------------------------
def f_broken_power_law(x, R, V, p_in, p_out):
    #Broken power law.
    #x are positions of the datapoints in AU and is an array.
    #R is the point where the powerlaw changes, in AU.
    #V is the velocity where the powerlaw changes, in ...
    #p_in is the power for x <= R.
    #p_out is the power for x > R.
    
    #We are working with an array, so for the if statements to work 
    #correctly we have to work element wise. 
    length = len(x)
    v = 0*x
    i = 0
    while i < length:
        if x[i] <= R:
            v[i] = V*(x[i]/R)**(-p_in)
        elif x[i] > R:
            v[i] = V*(x[i]/R)**(-p_out)
        else:
            print 'Error'
        i += 1
        #End while loop
    return v

def Chi_Squared(Science, Model):
    '''
    This function calculates the chi squared of the Model.
    Science = [float]
    Model = [float]
    '''
    #The shape of the array.
    n1, n2 = Science.shape

    Data = Science[:,0]
    Data_Err = Science[:,1]
    
    #Calculating chi squared.
    temp = np.square((Data - Model)/Data_Err)
    chi_squared = np.sum(temp)/n1
    return chi_squared

#-------------------------------------------------------------------------------
#Other functions.
#-------------------------------------------------------------------------------

def M_to_M_Sun(M):
    #The resulting constant from functions 1 and 2 is changed into solar masses.
    return np.square(M)/const.G.value*np.square(1000)*u.AU.to(u.m)/const.M_sun.value

def Random_Number(Lower_Bound, Upper_Bound):
    #This function creates an random number between the Lower and Upper bound.
    #With Lower and Upper bound Nx1 arrays. 
    length_array = len(Lower_Bound)
    return (Upper_Bound-Lower_Bound)*np.random.rand(length_array) + Lower_Bound

def Random_Starting_Parameters(N, Lower_R, Upper_R, Lower_V, Upper_V, Lower_P_in, Upper_P_in, Lower_P_out, Upper_P_out):
    #This function creates an .txt file containing N sets of random starting 
    #parameters for curve fitting of the f_broken_power_law function. The parameters
    #are random between the lower and upper bounds.
    begin_param_rand = np.random.rand(N,4)

    #Giving the beginparameters the correct units.
    begin_param_rand[:,0] = Random_Number(np.ones(N)*Lower_R, np.ones(N)*Upper_R)
    begin_param_rand[:,1] = Random_Number(np.ones(N)*Lower_V, np.ones(N)*Upper_V)
    begin_param_rand[:,2] = Random_Number(np.ones(N)*Lower_P_in, np.ones(N)*Upper_P_in)
    begin_param_rand[:,3] = Random_Number(np.ones(N)*Lower_P_out, np.ones(N)*Upper_P_out)

    np.savetxt(save_directory + 'Random_Begin_Parameters.txt', begin_param_rand, delimiter = ' ')

#-------------------------------------------------------------------------------
#Main functions
#-------------------------------------------------------------------------------
def Prepare_Data(Object, directory_1 = '/niks/', file_1 = 'niks.txt', directory_2 = '/niks/', file_2= 'niks.txt', Cutoff_Red_1 = False, cutoff_red_1 = 1, Cutoff_Blue_1 = False, cutoff_blue_1 = 1, Cutoff_Red_2 = False, cutoff_red_2 = 1, Cutoff_Blue_2 = False, cutoff_blue_2 = 1, Data1 = True, Data2 = True):
    '''
    Function that prepares the data for fitting.
    directory_1 = [string]
    file_1 = [string]
    directory_1 = [string]
    file_1 = [string]
    Object = [string]
    Cutoff_Red_1 = [Boolean]
    cutoff_red_1 = [float]
    Cutoff_Blue_1 = [Boolean]
    cutoff_blue_1 = [float]
    Cutoff_Red_2 = [Boolean]
    cutoff_red_2 = [float]
    Cutoff_Blue_2 = [Boolean]
    cutoff_blue_2 = [float]
    '''
    def Filter_Out_Nan(Data, Compare_Data):
        return Data[np.where(np.isnan(Compare_Data) == False)]  
         
    def Enforce_Cutoff(Data, Data_Velocity, Cutoff):
        return Data[np.where(Data_Velocity > Cutoff)]
    
    if Data1:
        #First read the data.    
        Data_1 = np.loadtxt(directory_1 + file_1, delimiter = ' ')  
        #The different columns are now selected.
        #Dataset 1
        radius_red_1 = Data_1[:,0]
        radius_red_err_1 = Data_1[:,1]
        velocity_red_1 = Data_1[:,2]
        velocity_red_err_1 = Data_1[:,3]  
        radius_blue_1 = Data_1[:,4]
        radius_blue_err_1 = Data_1[:,5]
        velocity_blue_1 = Data_1[:,6]
        velocity_blue_err_1 = Data_1[:,7] 
        
        #Filtering out all the nan data.
        temp_radius_red_1 = Filter_Out_Nan(radius_red_1, velocity_red_1)
        temp_radius_red_err_1 = Filter_Out_Nan(radius_red_err_1, velocity_red_1)
        temp_velocity_red_1 = Filter_Out_Nan(velocity_red_1, velocity_red_1)
        temp_velocity_red_err_1 = Filter_Out_Nan(velocity_red_err_1, velocity_red_1)

        temp_radius_blue_1 = Filter_Out_Nan(radius_blue_1, velocity_blue_1)
        temp_radius_blue_err_1 = Filter_Out_Nan(radius_blue_err_1, velocity_blue_1)
        temp_velocity_blue_1 = Filter_Out_Nan(velocity_blue_1, velocity_blue_1)
        temp_velocity_blue_err_1 = Filter_Out_Nan(velocity_blue_err_1, velocity_blue_1)
        
        #Applying the cutoff, everything below the cutoff velocities is removed.
        if Cutoff_Red_1:
            temp_radius_red_1 = Enforce_Cutoff(temp_radius_red_1, temp_velocity_red_1, cutoff_red_1)
            temp_radius_red_err_1 = Enforce_Cutoff(temp_radius_red_err_1, temp_velocity_red_1, cutoff_red_1)   
            temp_velocity_red_1 = Enforce_Cutoff(temp_velocity_red_1, temp_velocity_red_1, cutoff_red_1)
            temp_velocity_red_err_1 = Enforce_Cutoff(temp_velocity_red_err_1, temp_velocity_red_1, cutoff_red_1)
        if Cutoff_Blue_1:
            temp_radius_blue_1 = Enforce_Cutoff(temp_radius_blue_1, temp_velocity_blue_1, cutoff_blue_1)
            temp_radius_blue_err_1 = Enforce_Cutoff(temp_radius_blue_err_1, temp_velocity_blue_1, cutoff_blue_1)    
            temp_velocity_blue_1 = Enforce_Cutoff(temp_velocity_blue_1, temp_velocity_blue_1, cutoff_blue_1)
            temp_velocity_blue_err_1 = Enforce_Cutoff(temp_velocity_blue_err_1, temp_velocity_blue_1, cutoff_blue_1)
        
    if Data2:
        #First read the data.   
        Data_2 = np.loadtxt(directory_2 + file_2, delimiter = ' ')

        #The different columns are now selected.
        #Dataset 2
        radius_red_2 = Data_2[:,0]
        radius_red_err_2 = Data_2[:,1]
        velocity_red_2 = Data_2[:,2]
        velocity_red_err_2 = Data_2[:,3]  
        radius_blue_2 = Data_2[:,4]
        radius_blue_err_2 = Data_2[:,5]
        velocity_blue_2 = Data_2[:,6]
        velocity_blue_err_2 = Data_2[:,7]

        #Filtering out all the nan data.
        temp_radius_red_2 = Filter_Out_Nan(radius_red_2, velocity_red_2)
        temp_radius_red_err_2 = Filter_Out_Nan(radius_red_err_2, velocity_red_2)
        temp_velocity_red_2 = Filter_Out_Nan(velocity_red_2, velocity_red_2)
        temp_velocity_red_err_2 = Filter_Out_Nan(velocity_red_err_2, velocity_red_2)

        temp_radius_blue_2 = Filter_Out_Nan(radius_blue_2, velocity_blue_2)
        temp_radius_blue_err_2 = Filter_Out_Nan(radius_blue_err_2, velocity_blue_2)
        temp_velocity_blue_2 = Filter_Out_Nan(velocity_blue_2, velocity_blue_2)
        temp_velocity_blue_err_2 = Filter_Out_Nan(velocity_blue_err_2, velocity_blue_2)

        #Applying the cutoff, everything below the cutoff velocities is removed.
        if Cutoff_Red_2:
            temp_radius_red_2 = Enforce_Cutoff(temp_radius_red_2, temp_velocity_red_2, cutoff_red_2)
            temp_radius_red_err_2 = Enforce_Cutoff(temp_radius_red_err_2, temp_velocity_red_2, cutoff_red_2)  
            temp_velocity_red_2 = Enforce_Cutoff(temp_velocity_red_2, temp_velocity_red_2, cutoff_red_2)
            temp_velocity_red_err_2 = Enforce_Cutoff(temp_velocity_red_err_2, temp_velocity_red_2, cutoff_red_2)
        if Cutoff_Blue_2:
            temp_radius_blue_2 = Enforce_Cutoff(temp_radius_blue_2, temp_velocity_blue_2, cutoff_blue_2)
            temp_radius_blue_err_2 = Enforce_Cutoff(temp_radius_blue_err_2, temp_velocity_blue_2, cutoff_blue_2) 
            temp_velocity_blue_2 = Enforce_Cutoff(temp_velocity_blue_2, temp_velocity_blue_2, cutoff_blue_2)
            temp_velocity_blue_err_2 = Enforce_Cutoff(temp_velocity_blue_err_2, temp_velocity_blue_2, cutoff_blue_2)

    if Data1 and Data2:
        #Appending the two datasets into 1 array.
        radius_red = np.append(temp_radius_red_1, temp_radius_red_2)
        radius_red_err = np.append(temp_radius_red_err_1, temp_radius_red_err_2)
        velocity_red = np.append(temp_velocity_red_1, temp_velocity_red_2)
        velocity_red_err = np.append(temp_velocity_red_err_1, temp_velocity_red_err_2)

        radius_blue = np.append(temp_radius_blue_1, temp_radius_blue_2)
        radius_blue_err = np.append(temp_radius_blue_err_1, temp_radius_blue_err_2)
        velocity_blue = np.append(temp_velocity_blue_1, temp_velocity_blue_2)
        velocity_blue_err = np.append(temp_velocity_blue_err_1, temp_velocity_blue_err_2)
    elif Data1:
        radius_red = temp_radius_red_1
        radius_red_err = temp_radius_red_err_1
        velocity_red = temp_velocity_red_1
        velocity_red_err = temp_velocity_red_err_1

        radius_blue = temp_radius_blue_1
        radius_blue_err = temp_radius_blue_err_1
        velocity_blue = temp_velocity_blue_1
        velocity_blue_err = temp_velocity_blue_err_1   
    elif Data2:
        radius_red = temp_radius_red_2
        radius_red_err = temp_radius_red_err_2
        velocity_red = temp_velocity_red_2
        velocity_red_err = temp_velocity_red_err_2

        radius_blue = temp_radius_blue_2
        radius_blue_err = temp_radius_blue_err_2
        velocity_blue = temp_velocity_blue_2
        velocity_blue_err = temp_velocity_blue_err_2 
        
    #---------------------------------------------------------------------------
    #Also combining the blue and redshifted data for fitting. 
    #---------------------------------------------------------------------------
    radius = np.append(radius_red, radius_blue)
    radius_err = np.append(radius_red_err, radius_blue_err)
    velocity = np.append(velocity_red, velocity_blue)
    velocity_err = np.append(velocity_red_err, velocity_blue_err)

    #---------------------------------------------------------------------------
    #Saving the data.
    #---------------------------------------------------------------------------

    #First we save the redshifted data.     
    save_data_red = np.zeros([len(radius_red), 4])
    save_data_red[:,0] = radius_red
    save_data_red[:,1] = radius_red_err
    save_data_red[:,2] = velocity_red
    save_data_red[:,3] = velocity_red_err    
    np.savetxt(save_directory + 'Red_Shifted_Data_' + Object + '.txt', save_data_red, delimiter = ' ')

    #Then the blueshifted data.
    save_data_blue = np.zeros([len(radius_blue), 4])
    save_data_blue[:,0] = radius_blue
    save_data_blue[:,1] = radius_blue_err
    save_data_blue[:,2] = velocity_blue
    save_data_blue[:,3] = velocity_blue_err    
    np.savetxt(save_directory + 'Blue_Shifted_Data_' + Object + '.txt', save_data_blue, delimiter = ' ')

    #And finally all the data. 
    save_data_all = np.zeros([len(radius), 4])
    save_data_all[:,0] = radius
    save_data_all[:,1] = radius_err
    save_data_all[:,2] = velocity
    save_data_all[:,3] = velocity_err    
    np.savetxt(save_directory + 'All_Data_' + Object + '.txt', save_data_all, delimiter = ' ')   

def Fit_Data_Random_Starting_Parameters(Data, Object):
    '''
    The function that fits the data, using random starting parameters, to the broken 
    power law. 
    Data = [string]
    Object = [string]
    '''
    #---------------------------------------------------------------------------
    #Loading the data the we have to fit.
    #---------------------------------------------------------------------------
    data = np.loadtxt(save_directory + Data, delimiter = ' ')
    #Unpacking the data.
    radius = data[:,0]
    radius_err = data[:,1]
    velocity = data[:,2]
    velocity_err = data[:,3]

    #Loading the random begin parameters.
    begin_param_rand = np.loadtxt(save_directory + 'Random_Begin_Parameters.txt', delimiter = ' ') 

    #Creating the array containing the results.
    param_result = np.zeros([begin_param_rand.shape[0], begin_param_rand.shape[1]])

    #---------------------------------------------------------------------------
    #Fitting the data to the random starting parameters.
    #---------------------------------------------------------------------------
    print 'Starting with fitting the random parameters.'
    i = 0
    while i < begin_param_rand.shape[0]:
        print str(i+1) + '/' + str(begin_param_rand.shape[0])
        param_result[i,:] = curve_fit(f_broken_power_law, radius, velocity, p0 = begin_param_rand[i,:], maxfev=1000000)[0]
        i+=1

    #---------------------------------------------------------------------------
    #Saving the results. 
    #---------------------------------------------------------------------------
    print 'Saving the results.'
    np.savetxt(save_directory + 'Results_Random_Starting_Parameters_' + Object + '.txt', param_result, delimiter = ' ')

def Plot_Results_Random_Starting_Parameters(Data, Object, R_select, V_select, p_in_select, p_out_select, N):
    '''
    This function plots the results from the fits with the random starting 
    parameters.
    Data = [string]
    Object = [string]
    R_select = [float]
    V_select = [float]
    p_in_select = [float]
    p_out_select = [float]
    N = [int]
    '''
    #---------------------------------------------------------------------------    
    #Loading the data.
    #---------------------------------------------------------------------------
    param_result = np.loadtxt(save_directory + Data, delimiter = ' ')
    #Loading the random begin parameters.
    begin_param_rand = np.loadtxt(save_directory + 'Random_Begin_Parameters.txt', delimiter = ' ') 

    #---------------------------------------------------------------------------
    #Selecting data.
    #---------------------------------------------------------------------------
    begin_param_R = begin_param_rand[:,0][begin_param_rand[:,0] < R_select]
    result_param_R = param_result[:,0][param_result[:,0] < R_select] 

    begin_param_V = begin_param_rand[:,1][begin_param_rand[:,1] < V_select]     
    result_param_V = param_result[:,1][param_result[:,1] < V_select]

    begin_param_p_in = begin_param_rand[:,2][begin_param_rand[:,2] < p_in_select]
    result_param_p_in = param_result[:,2][param_result[:,2] < p_in_select]

    begin_param_p_out = begin_param_rand[:,3][begin_param_rand[:,3] < p_out_select]
    result_param_p_out = param_result[:,3][param_result[:,3] < p_out_select]

    #---------------------------------------------------------------------------
    #Histogram of the break radius.
    #---------------------------------------------------------------------------
    binBoundaries = np.linspace(0,R_select,101)
    plt.hist(begin_param_R, bins = binBoundaries, histtype='stepfilled', color='b', label='Begin Parameters')
    plt.hist(result_param_R, bins = binBoundaries, histtype='stepfilled', color='r', alpha=0.5, label='Result Parameters')
    plt.title(Object + "'s histogram of $R_{b}$, N = " + str(N))
    plt.xlabel('$R_{b}$ [$AU$]')
    plt.ylabel('#')
    plt.xlim([0,R_select])
    plt.legend()
    plt.savefig(save_directory + 'Histogram_Radius_Break_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Histogram_Radius_Break_' + Object + '.pdf', bbox_inches='tight')
    plt.close() 

    #---------------------------------------------------------------------------
    #Histogram of the break velocity.
    #---------------------------------------------------------------------------
    binBoundaries = np.linspace(0,V_select,101)
    plt.hist(begin_param_V, bins = binBoundaries, histtype='stepfilled', color='b', label='Begin Parameters')
    plt.hist(result_param_V, bins = binBoundaries, histtype='stepfilled', color='r', alpha=0.5, label='Result Parameters')
    plt.title(Object + "'s histogram of $V_{b}$, N = " + str(N))
    plt.xlabel('$V_{b}$ [$km$ $s^{-1}$]')
    plt.ylabel('#')
    plt.xlim([0, V_select])
    plt.legend()
    plt.savefig(save_directory + 'Histogram_Velocity_Break_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Histogram_Velocity_Break_' + Object + '.pdf', bbox_inches='tight')
    plt.close() 

    #---------------------------------------------------------------------------
    #Histogram of the inner power.
    #---------------------------------------------------------------------------
    binBoundaries = np.linspace(0,p_in_select,101)
    plt.hist(begin_param_p_in, bins = binBoundaries, histtype='stepfilled', color='b', label='Begin Parameters')
    plt.hist(result_param_p_in, bins = binBoundaries, histtype='stepfilled', color='r', alpha=0.5, label='Result Parameters')
    plt.title(Object + "'s histogram of $p_{in}$, N = " + str(N))
    plt.xlabel('$p_{in}$')
    plt.ylabel('#')
    plt.xlim([0,p_in_select])
    plt.legend()
    plt.savefig(save_directory + 'Histogram_P_In_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Histogram_P_In_' + Object + '.pdf', bbox_inches='tight')
    plt.close() 

    #---------------------------------------------------------------------------
    #Histogram of the outer power.
    #---------------------------------------------------------------------------
    binBoundaries = np.linspace(0,p_out_select,101)
    plt.hist(begin_param_p_out, bins = binBoundaries, histtype='stepfilled', color='b', label='Begin Parameters')
    plt.hist(result_param_p_out, bins = binBoundaries, histtype='stepfilled', color='r', alpha=0.5, label='Result Parameters')
    plt.title(Object + "'s histogram of $p_{out}$, N = " + str(N))
    plt.xlabel('$p_{out}$')
    plt.ylabel('#')
    plt.xlim([0,p_out_select])
    plt.legend()
    plt.savefig(save_directory + 'Histogram_P_Out_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Histogram_P_Out_' + Object + '.pdf', bbox_inches='tight')
    plt.close() 

def Monte_Carlo_Data(Object, Data_Random_Parameters, Data_Science, R_low, R_high, iterations):
    '''
    This function takes the result obtained above and tries to improve the results
    by doing a Monte Carlo simulation.
    Object = [string]
    Data_Random_Parameters = [string]
    Data_Science = [string]
    R_low = [float]
    R_high = [float]
    iterations = [int]
    '''
    #---------------------------------------------------------------------------
    #Loading the data.
    #---------------------------------------------------------------------------
    Data_RP = np.loadtxt(save_directory + Data_Random_Parameters, delimiter = ' ')
    Data_S = np.loadtxt(save_directory + Data_Science, delimiter = ' ')

    #Unpacking the science data.
    radius = Data_S[:,0]
    radius_err = Data_S[:,1]
    velocity = Data_S[:,2]
    velocity_err = Data_S[:,3]

    #Selecting the used data from Data_RP, which is between R_low and R_high.
    Data_RP_used = Data_RP[(Data_RP[:,0] <= R_high)*(Data_RP[:,0] >= R_low),:]

    #Shape of the array.    
    y,x = np.shape(Data_RP_used)

    #Array where the results will be saved. 
    Results = np.zeros((iterations*y,4))  

    #---------------------------------------------------------------------------
    #The Monte Carlo simulation.
    #---------------------------------------------------------------------------
    print 'Starting the Monte Carlo simulation.'
    total_iteration = 0
    i = 0
    while i < y:
        #Choosing the starting parameters.
        Param_Used = Data_RP_used[i,:]       
        j = 0
        while j < iterations:
            print 'Iteration: ' + str(total_iteration + 1) + '/' + str(iterations*y)
            #Creating random starting positions within the error bars. 
            temp_radius = Random_Number(radius - radius_err, radius + radius_err)
            temp_velocity = Random_Number(velocity - velocity_err, velocity + velocity_err)

            #Doing the fitting.
            Results[total_iteration,:] = curve_fit(f_broken_power_law, temp_radius, temp_velocity, p0 = Param_Used, maxfev=1000000)[0]
            total_iteration += 1
            j += 1
            #End while loop
        i += 1 
        #End while loop
    print 'Monte Carlo simulation is finished.'
    
    #---------------------------------------------------------------------------
    #Saving the results from the Monte Carlo simulation.
    #---------------------------------------------------------------------------
    np.savetxt(save_directory + 'Results_Monte_Carlo_' + Object + '.txt', Results, delimiter = ' ') 

def Plot_Results_Monte_Carlo(Data, Object, R_select, V_select, p_in_select, p_out_select, N):
    '''
    This function plots the results from the Monte Carlo simulation. 
    Data = [string]
    Object = [string]
    R_select = [float]
    V_select = [float]
    p_in_select = [float]
    p_out_select = [float]
    N = [int]
    '''
    #---------------------------------------------------------------------------
    #Loading the data.
    #---------------------------------------------------------------------------
    Result = np.loadtxt(save_directory + Data, delimiter = ' ')

    #---------------------------------------------------------------------------
    #Selecting data.
    #---------------------------------------------------------------------------
    result_param_R = Result[:,0][Result[:,0] < R_select] 
    
    result_param_V = Result[:,1][Result[:,1] < V_select]

    result_param_p_in = Result[:,2][Result[:,2] < p_in_select]

    result_param_p_out = Result[:,3][Result[:,3] < p_out_select]

    #---------------------------------------------------------------------------
    #Histogram of the break radius.
    #---------------------------------------------------------------------------
    binBoundaries = np.linspace(0,R_select,101)
    plt.hist(result_param_R, bins = binBoundaries, histtype='stepfilled', color='r', alpha = 0.5)
    plt.title(Object + "'s histogram of $R_{b}$, Monte Carlo N = " + str(N))
    plt.xlabel('$R_{b}$ [$AU$]')
    plt.ylabel('#')
    plt.xlim([0,R_select])
    plt.savefig(save_directory + 'Histogram_Radius_Break_Monte_Carlo_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Histogram_Radius_Break_Monte_Carlo_' + Object + '.pdf', bbox_inches='tight')
    plt.close() 

    #---------------------------------------------------------------------------
    #Histogram of the break velocity.
    #---------------------------------------------------------------------------
    binBoundaries = np.linspace(0,V_select,101)
    plt.hist(result_param_V, bins = binBoundaries, histtype='stepfilled', color='r', alpha = 0.5)
    plt.title(Object + "'s histogram of $V_{b}$, Monte Carlo N = " + str(N))
    plt.xlabel('$V_{b}$ [$km$ $s^{-1}$]')
    plt.ylabel('#')
    plt.xlim([0,V_select])
    plt.savefig(save_directory + 'Histogram_Velocity_Break_Monte_Carlo_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Histogram_Velocity_Break_Monte_Carlo_' + Object + '.pdf', bbox_inches='tight')
    plt.close() 

    #---------------------------------------------------------------------------
    #Histogram of the inner power.
    #---------------------------------------------------------------------------
    binBoundaries = np.linspace(0,p_in_select,101)
    plt.hist(result_param_p_in, bins = binBoundaries, histtype='stepfilled', color='r', alpha = 0.5)
    plt.title(Object + "'s histogram of $p_{in}$, Monte Carlo N = " + str(N))
    plt.xlabel('$p_{in}$')
    plt.ylabel('#')
    plt.xlim([0,p_in_select])
    plt.savefig(save_directory + 'Histogram_P_In_Monte_Carlo_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Histogram_P_In_Monte_Carlo_' + Object + '.pdf', bbox_inches='tight')
    plt.close() 

    #---------------------------------------------------------------------------
    #Histogram of the outer power.
    #---------------------------------------------------------------------------
    binBoundaries = np.linspace(0,p_out_select,101)
    plt.hist(result_param_p_out, bins = binBoundaries, histtype='stepfilled', color='r', alpha = 0.5)
    plt.title(Object + "'s histogram of $p_{out}$, Monte Carlo N = " + str(N))
    plt.xlabel('$p_{out}$')
    plt.ylabel('#')
    plt.xlim([0,p_out_select])
    plt.savefig(save_directory + 'Histogram_P_Out_Monte_Carlo_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Histogram_P_Out_Monte_Carlo_' + Object + '.pdf', bbox_inches='tight')
    plt.close() 
 
def Fit_Gaussian_to_Parameter_Distributions(Data, Object, R_select, V_select, p_in_select, p_out_select):
    '''
    This function fits Gaussians to the distributions of the parameters coming from
    the monte carlo simulation.
    This function is not very flexible. 
    Data = [string]
    Object = [string]
    R_select = [float]
    V_select = [float]
    p_in_select = [float]
    p_out_select = [float]
    '''
   #---------------------------------------------------------------------------
    #Loading the data.
    #---------------------------------------------------------------------------
    Result = np.loadtxt(save_directory + Data, delimiter = ' ')

    #---------------------------------------------------------------------------
    #Selecting data.
    #---------------------------------------------------------------------------
    result_param_R = Result[:,0][Result[:,0] < R_select] 
    
    result_param_V = Result[:,1][Result[:,1] < V_select]

    result_param_p_in = Result[:,2][Result[:,2] < p_in_select]

    result_param_p_out = Result[:,3][Result[:,3] < p_out_select]
    #---------------------------------------------------------------------------
    #A function defining a 1D Gaussian.
    #---------------------------------------------------------------------------
    def Gaussian(x, A, x_center, sigma):
        return A*np.exp(-np.square(x - x_center)/(2*np.square(sigma)))
        #End function

    #---------------------------------------------------------------------------
    #First we fit R_b
    #---------------------------------------------------------------------------
    #Boundaries for the bins.
    binBoundaries = np.linspace(0,R_select,101)
    #Creating the binned data, y axis and the values to plot on the x axis.
    temp = plt.hist(result_param_R, bins = binBoundaries, color = 'blue')
    y = temp[0]
    x = temp[1][0:100] + + 0.5*R_select/100

    #Fitting the values to an Gaussian.
    #We only fit the values where R < 300 AU.
    R_fit = 60
    R_fit_width = 10
    Param_Gaussian_R_b = curve_fit(Gaussian, x[abs(x - R_fit) < R_fit_width], y[abs(x - R_fit) < R_fit_width], p0 = [140000, 60, 1])[0]
    
    #The arrays for plotting.
    R_plot = np.linspace(0,R_select,1001)
    Gaussian_Fit = Gaussian(R_plot, Param_Gaussian_R_b[0], Param_Gaussian_R_b[1], Param_Gaussian_R_b[2])

    R_b = Param_Gaussian_R_b[1]
    R_b_err = np.abs(Param_Gaussian_R_b[2])    

    Label = '$r_c$ = ' + str(np.round(Param_Gaussian_R_b[1], 0)) + ' AU,\n$\sigma$ = ' + str(np.round(Param_Gaussian_R_b[2],0)) + ' AU'

    plt.plot(R_plot, Gaussian_Fit, color = 'r', linewidth = 2, label = Label)
    #plt.title('Histogram Monte Carlo $R_b$ and Gaussian fit for $|R_b - ' + str(R_fit) +'| < $ $' + str(R_fit_width) +'$ AU')
    plt.xlabel('$R_b$ [Au]')
    plt.ylabel('$\#$')
    plt.xlim([50,70])
    plt.ylim([0, 180000])
    plt.grid()
    plt.legend()
    plt.savefig(save_directory + 'Histogram_Radius_Break_Monte_Carlo_Gaussian_Fit_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Histogram_Radius_Break_Monte_Carlo_Gaussian_Fit_' + Object + '.pdf', bbox_inches='tight')
    plt.close()

    #---------------------------------------------------------------------------
    #Now we fit V_b.
    #---------------------------------------------------------------------------
    #Boundaries for the bins.
    binBoundaries = np.linspace(0,V_select,101)
    #Creating the binned data, y axis and the values to plot on the x axis.
    temp = plt.hist(result_param_V, bins = binBoundaries)
    y = temp[0]
    x = temp[1][0:100] + 0.5*V_select/100


    #Fitting the values to an Gaussian.
    #We only fit the values where V > 1.25 km s^-1.
    V_fit = 2.5
    V_fit_width = 0.5
    Param_Gaussian_V_b = curve_fit(Gaussian, x[abs(x - V_fit) < V_fit_width], y[abs(x - V_fit) < V_fit_width], p0 = [90000, 2.5, 0.25])[0]
    
    #The arrays for plotting.
    V_plot = np.linspace(0, V_select,1001)
    Gaussian_Fit = Gaussian(V_plot, Param_Gaussian_V_b[0], Param_Gaussian_V_b[1], Param_Gaussian_V_b[2])

    V_b = Param_Gaussian_V_b[1]
    V_b_err = np.abs(Param_Gaussian_V_b[2])  

    Label = '$v_c$ = ' + str(np.round(Param_Gaussian_V_b[1], 1)) + ' km s$^{-1}$,\n$\sigma$ = ' + str(np.round(Param_Gaussian_V_b[2],1)) + ' km s$^{-1}$'
    
    plt.plot(V_plot, Gaussian_Fit, color = 'r', linewidth = 2, label = Label )
    #plt.title('Histogram Monte Carlo $V_b$ and Gaussian fit for $|V_b - ' + str(V_fit) +'| < $ $' + str(V_fit_width) + '$ km s$^{-1}$')
    plt.xlabel('$V_b$ [km s$^{-1}$]')
    plt.ylabel('$\#$')
    plt.xlim([1.8, 2.8])
    plt.ylim([0, 120000])
    plt.grid()
    plt.legend()
    plt.savefig(save_directory + 'Histogram_Velocity_Break_Monte_Carlo_Gaussian_Fit_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Histogram_Velocity_Break_Monte_Carlo_Gaussian_Fit_' + Object + '.pdf', bbox_inches='tight')
    plt.close()

    #---------------------------------------------------------------------------
    #And we fit P_in.
    #---------------------------------------------------------------------------
    #Boundaries for the bins.
    binBoundaries = np.linspace(0,p_in_select,101)
    #Creating the binned data, y axis and the values to plot on the x axis.
    temp = plt.hist(result_param_p_in, bins = binBoundaries)
    y = temp[0]
    x = temp[1][0:100] + 0.5*p_in_select/100

    #Fitting the values to an Gaussian.
    #We only fit the values where P_in > 0.4 .
    P_in_fit = 0.5
    P_in_fit_width = 0.3
    Param_Gaussian_P_in = curve_fit(Gaussian, x[abs(x - P_in_fit) < P_in_fit_width], y[abs(x - P_in_fit) < P_in_fit_width], p0 = [60000, 0.5, 0.2])[0]
    
    P_in = Param_Gaussian_P_in[1]
    P_in_err = np.abs(Param_Gaussian_P_in[2])    

    #The arrays for plotting.
    P_in_plot = np.linspace(0,p_in_select,1001)
    Gaussian_Fit = Gaussian(P_in_plot, Param_Gaussian_P_in[0], Param_Gaussian_P_in[1], Param_Gaussian_P_in[2])

    Label = '$P_{in}$ = ' + str(np.round(Param_Gaussian_P_in[1], 2)) + ',\n$\sigma$ = ' + str(np.round(Param_Gaussian_P_in[2],2))
    
    plt.plot(P_in_plot, Gaussian_Fit, color = 'r', linewidth = 2, label = Label)
    #plt.title('Histogram Monte Carlo $P_{in}$ and Gaussian fit for $|P_{in} - ' + str(P_in_fit) + ' | < $ $' + str(P_in_fit_width) + '$')
    plt.xlabel('$P_{in}$')
    plt.ylabel('$\#$')
    plt.xlim([0,1])
    plt.grid()
    plt.legend(loc = 'upper right')
    plt.savefig(save_directory + 'Histogram_P_In_Monte_Carlo_Gaussian_Fit_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Histogram_P_In_Monte_Carlo_Gaussian_Fit_' + Object + '.pdf', bbox_inches='tight')
    plt.close()

    #---------------------------------------------------------------------------
    #Finally we fit P_out.
    #---------------------------------------------------------------------------
    #Boundaries for the bins.
    binBoundaries = np.linspace(0,p_out_select,101)
    #Creating the binned data, y axis and the values to plot on the x axis.
    temp = plt.hist(result_param_p_out, bins = binBoundaries)
    y = temp[0]
    x = temp[1][0:100]+ 0.5*p_out_select/100

    #Fitting the values to an Gaussian.
    #We only fit the values where P_out < 1.25 .
    P_out_fit = 1.6
    P_out_fit_width = 0.5 
    Param_Gaussian_P_out = curve_fit(Gaussian, x[abs(x - P_out_fit) < P_out_fit_width], y[abs(x - P_out_fit) < P_out_fit_width], p0 = [34000, 1.6, 0.5])[0]

    P_out = Param_Gaussian_P_out[1]
    P_out_err = np.abs(Param_Gaussian_P_out[2]) 
    
    #The arrays for plotting.
    P_out_plot = np.linspace(0,p_out_select,1001)
    Gaussian_Fit = Gaussian(P_out_plot, Param_Gaussian_P_out[0], Param_Gaussian_P_out[1], Param_Gaussian_P_out[2])
    
    Label = '$P_{out}$ = ' + str(np.round(Param_Gaussian_P_out[1], 2)) + ',\n$\sigma$ = ' + str(np.round(np.abs(Param_Gaussian_P_out[2]),2))

    plt.plot(P_out_plot, Gaussian_Fit, color = 'r', linewidth = 2, label = Label)
    #plt.title('Histogram Monte Carlo $P_{out}$ and Gaussian fit for $|P_{out} - ' + str(P_out_fit) + '  | < $ $' + str(P_out_fit_width) + '$')
    plt.xlabel('$P_{out}$')
    plt.ylabel('$\#$')
    plt.xlim([0.5,3])
    plt.grid()
    plt.legend()
    plt.savefig(save_directory + 'Histogram_P_Out_Monte_Carlo_Gaussian_Fit_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Histogram_P_Out_Monte_Carlo_Gaussian_Fit_' + Object + '.pdf', bbox_inches='tight')
    plt.close()
    #---------------------------------------------------------------------------
    #Saving the resulting parameters.
    #---------------------------------------------------------------------------
    parameters = [R_b, V_b, P_in, P_out]
    parameters_err = [R_b_err, V_b_err, P_in_err, P_out_err]    
    
    np.savetxt(save_directory + 'Parameters_' + Object + '.txt', parameters, delimiter = ' ') 
    np.savetxt(save_directory + 'Parameters_Errors_' + Object + '.txt', parameters_err, delimiter = ' ')  

def Calculate_Mass(R_b, V_b, R_b_err, V_b_err, inclination, output = False):
    '''
    Now we calculate the mass assuming Keplerian rotation.
    R_b = [float]
    V_b = [float]
    R_b_err = [float]
    V_b_err = [float]
    inclination = [float]
    '''
    #Calculating the mass.
    M = R_b*u.AU.to(u.m)*np.square(V_b*1000)/const.G.value/const.M_sun.value/(np.square(np.sin(np.radians(inclination))))#Solar masses

    #Calculating the error on the mass.
    M_err = np.sqrt(np.square(np.square(V_b*1000)/const.G.value*R_b_err*u.AU.to(u.m)) + np.square(2*R_b*u.AU.to(u.m)*(V_b*1000)/const.G.value*V_b_err*1000))/const.M_sun.value/(np.square(np.sin(np.radians(inclination))))#Solar masses
    
    #Printing the output.
    print 'Mass central object for i = ' + str(inclination) + ' deg. is: ' + str(np.round(M,2)) + '+-' + str(np.round(M_err,2)) + ' Solar masses'

    if output:
        return M, M_err

def Plot_Final_PV_Diagram(Object, data_red, data_blue, R_b, R_b_err, V_b, V_b_err, P_in, P_in_err, P_out, P_out_err, grid = False):
    '''
    Plotting the final peak PV diagram including the broken power law fit.
    Object = [string]
    data_red = [float]
    data_blue = [float]
    R_b = [float]
    R_b_err = [float]
    V_b = [float]
    V_b_err = [float]
    P_in = [float]
    P_in_err = [float]
    P_out = [float]
    P_out_err = [float]
    grid = [boolean]
    '''
    #---------------------------------------------------------------------------
    #Loading the data.
    #---------------------------------------------------------------------------
    Data_Red = np.loadtxt(save_directory + data_red, delimiter = ' ')
    Data_Blue = np.loadtxt(save_directory + data_blue, delimiter = ' ')

    #Unpacking the science data.
    radius_red = Data_Red[:,0]
    radius_red_err = Data_Red[:,1]
    velocity_red = Data_Red[:,2]
    velocity_red_err = Data_Red[:,3] 
   
    radius_blue = Data_Blue[:,0]
    radius_blue_err = Data_Blue[:,1]
    velocity_blue = Data_Blue[:,2]
    velocity_blue_err = Data_Blue[:,3]   

    #---------------------------------------------------------------------------
    #Some things necessary for plotting.
    #---------------------------------------------------------------------------
    X = np.arange(1., 150.,1.)#AU
    v_final = f_broken_power_law(X, R_b, V_b, P_in, P_out)

    #---------------------------------------------------------------------------
    #Plotting.
    #---------------------------------------------------------------------------
    fig = plt.figure()
    ax = fig.add_subplot(111)

    #Plotting the data.
    ax.errorbar(radius_red, velocity_red, xerr = radius_red_err, yerr= velocity_red_err, fmt='.', color = 'red')
    ax.errorbar(radius_blue, velocity_blue, xerr = radius_blue_err, yerr= velocity_blue_err, fmt='.', color = 'blue')

    #Plotting the broken power law.
    #Label = 'p$_{in} = $' + str(np.round(P_in, 2)) + '$\pm$' + str(np.round(P_in_err,2)) + ', p$_{out}$ = ' + str(np.round(P_out, 2)) + '$\pm$' + str(np.round(P_out_err, 2))
    Label = 'v(r) $\propto$ r$^{' + str(np.round(P_in, 2)) + '\pm' + str(np.round(P_in_err,2)) +  ' }$ if r $< ' + str(np.round(R_b))[:-2] + '\pm' + str(np.round(R_b_err))[:-2] + '$ AU \nv(r) $\propto$ r$^{' + str(np.round(P_out, 2)) + '\pm' + str(np.round(P_out_err,2)) +  ' }$ if r $> ' + str(np.round(R_b))[:-2] + '\pm' + str(np.round(R_b_err))[:-2] + '$ AU '
    ax.plot(X, v_final, label = Label, color = 'green')

    #Plotting the break point.
    Label = 'break point (' + str(np.round(R_b))[:-2] + '$\pm$' + str(np.round(R_b_err))[:-2] + ' AU , ' + str(np.round(V_b,1)) + '$\pm$' + str(np.round(V_b_err, 2)) + ' km s$^{-1}$)'
    ax.plot(R_b, V_b, 'og', label = Label)

    #Plotting the legend.
    #ax.legend(loc = 'lower left', numpoints = 1)
    
    #Set titles, labels etc.
    #ax.set_title('The Peak PV diagram of ' + Object + ' using C$^{18}$O')
    ax.set_xlabel('$r$ [AU]')
    ax.set_ylabel('$v$ [km s$^{-1}$]')
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
    '''
    #Rotates the minor ticks.
    labels = ax.get_xticklabels('minor')
    for label in labels:
        label.set_rotation(45) 
    #Rotates the major ticks.
    labels = ax.get_xticklabels()
    for label in labels:
        label.set_rotation(45) 
    '''
    #Creates an grid for minor and major ticks.
    if grid:
        ax.grid(b=True, which='major', color='black', linestyle='--')
        ax.grid(b=True, which='minor', color='black', linestyle='--')

    #ax.text(85, 2.8, 'C$^{18}$O', fontsize = 22)
    #ax.text(35, 2.2, 'v $\propto$ r$^{-' + str(np.round(P_in, 2)) + '\pm' + str(np.round(P_in_err,2)) +  '}$', fontsize = 22)
    #ax.text(70, 2.2, 'v $\propto$ r$^{-' + str(np.round(P_out, 2)) + '\pm' + str(np.round(P_out_err,2)) +  '}$', fontsize = 22)
    #ax.text(30, 1.1, 'M$_*$ = 0.34$\pm$0.03 M$_{\odot}$',fontsize = 22)
    #ax.text(30, 1.0, 'r$_c$ \ \ = ' + str(np.round(R_b))[:-2] + '$\pm$' + str(np.round(R_b_err))[:-2] + ' AU',fontsize = 22)

    #Saving.
    plt.savefig(save_directory + 'Final_Peak_PV-Diagram_' + Object + '.jpg', bbox_inches='tight')
    plt.savefig(save_directory + 'Final_Peak_PV-Diagram_' + Object + '.pdf', bbox_inches='tight')
    plt.close() 

def Calculate_Model_Datapoints_Broken_Power_Law(Data, R_b, V_b, P_in, P_out):
    radius = Data[:,0]
    
    return f_broken_power_law(radius, R_b, V_b, P_in, P_out)

#-------------------------------------------------------------------------------
#Function initializations
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    #---------------------------------------------------------------------------
    #Data information.
    #---------------------------------------------------------------------------
    Object = 'L1527'
    file_1 = 'Peak_PV-Diagram_L1527_C18O_2.txt'
    directory_1 = 'C18O_2/Peak_PV-diagram/'


    #---------------------------------------------------------------------------
    #First we prepare the data for fitting.
    #---------------------------------------------------------------------------
    Cutoff_Red_1 = False
    Cutoff_Blue_1 = False
    Cutoff_Red_2 = False
    Cutoff_Blue_2 = False
    #Prepare_Data(Object,directory_1, file_1, Data2 = False)
    
    #---------------------------------------------------------------------------
    #Then we create some randoms numbers necessary for the fitting. 
    #---------------------------------------------------------------------------
    N = 1000
    Lower_R = 0 #AU
    Upper_R = 100 #AU
    Lower_V = 0 #km s^-1
    Upper_V = 3 #km s^-1
    Lower_P_in = 0
    Upper_P_in = 2
    Lower_P_out = 0
    Upper_P_out = 2
    #Random_Starting_Parameters(N, Lower_R, Upper_R, Lower_V, Upper_V, Lower_P_in, Upper_P_in, Lower_P_out, Upper_P_out)

    #---------------------------------------------------------------------------
    #Now we fit the data to random starting parameters.
    #---------------------------------------------------------------------------
    Data_1 = 'All_Data_' + Object + '.txt'
    #Fit_Data_Random_Starting_Parameters(Data_1, Object)

    #---------------------------------------------------------------------------
    #And we plot the results.
    #---------------------------------------------------------------------------
    Data_2 = 'Results_Random_Starting_Parameters_' + Object + '.txt'
    R_select = 100 #AU
    V_select = 3 #km s^-1
    p_in_select = 2
    p_out_select = 2
    #Plot_Results_Random_Starting_Parameters(Data_2, Object, R_select, V_select, p_in_select, p_out_select, N)

    #---------------------------------------------------------------------------
    #Using these we results we also do a Monte Carlo simulation to improve the 
    #fits.
    #---------------------------------------------------------------------------
    R_low = 50 #AU
    R_high = 70 #AU
    iterations = 1000
    #Monte_Carlo_Data(Object, Data_2, Data_1, R_low, R_high, iterations)

    #---------------------------------------------------------------------------
    #We plot the results.
    #---------------------------------------------------------------------------
    Data_3 =  'Results_Monte_Carlo_' + Object + '.txt'
    R_select_2 = 100 #AU
    V_select_2 = 3 #km s^-1
    p_in_select_2 = 2
    p_out_select_2 = 2.5
    #Plot_Results_Monte_Carlo(Data_3, Object, R_select_2, V_select_2, p_in_select_2, p_out_select_2, iterations)

    #---------------------------------------------------------------------------
    #Now we fit Gaussians to the histograms.
    #---------------------------------------------------------------------------
    Fit_Gaussian_to_Parameter_Distributions(Data_3, Object, R_select_2, V_select_2, p_in_select_2, p_out_select_2)

    parameters = np.loadtxt(save_directory + 'Parameters_' + Object + '.txt', delimiter = ' ')
    parameters_err = np.loadtxt(save_directory + 'Parameters_Errors_' + Object + '.txt', delimiter = ' ')

    #Unpacking the parameters.
    R_b = parameters[0]
    R_b_err = parameters_err[0]
    V_b = parameters[1]
    V_b_err = parameters_err[1]
    P_in = parameters[2]
    P_in_err = parameters_err[2]
    P_out = parameters[3]
    P_out_err = parameters_err[3]

    #---------------------------------------------------------------------------
    #Calculating the mass.
    #---------------------------------------------------------------------------
    inclination = 89 #Degrees
    Calculate_Mass(R_b, V_b, R_b_err, V_b_err, inclination)
    
    #-------------------------------------------------------------------
    #Here we calculate the chi squared values of the fits.
    #-------------------------------------------------------------------
    All_Data = np.loadtxt(save_directory + 'All_Data_' + Object + '.txt', delimiter = ' ') 
    Broken_Power_Law_Data = Calculate_Model_Datapoints_Broken_Power_Law(All_Data, R_b, V_b, P_in, P_out)
    
    Chi_Squared_Broken_Power_Law = Chi_Squared(All_Data[:,2:4], Broken_Power_Law_Data)
    print 'The chi squared values for the broken power law fit is: ' + str(Chi_Squared_Broken_Power_Law)
 

    #---------------------------------------------------------------------------
    #Plotting the final result.
    #---------------------------------------------------------------------------
    Data_red = 'Red_Shifted_Data_' + Object + '.txt'
    Data_blue = 'Blue_Shifted_Data_' + Object + '.txt'
    Grid = False

    Plot_Final_PV_Diagram(Object, Data_red, Data_blue, R_b, R_b_err, V_b, V_b_err, P_in, P_in_err, P_out, P_out_err, Grid)
'''
#-------------------------------------------------------------------------------
Unused functions
#-------------------------------------------------------------------------------
#Keplerian rotation, power and mass (including some constantes) are fitted.
def f_1(x, A, B): 
    return A*(x**B)
    #End function

#Now we fit only the mass (including some constantes). Still Keplerian 
#rotation. 
def f_2(x, A): 
    return A*(x**(-0.5))
    #End function

#Conserved angular momentum, fit R_crit, M, C. Input needs to be in S.I. 
#units te be correct. 
def f_3(x, R_crit, M, C): 
    return C*np.sqrt(const.G.value*M*R_crit)/x
    #End function

#Conserved angular momentum, fit only the total constant A.
def f_4(x, A): 
    return A/x
    #End function


#Broken power law with fixed powers..
#x are positions of the datapoints in AU and is an array.
#R is the point where the powerlaw changes, in AU.
#V is the velocity where the powerlaw changes, in ...
def f_broken_power_law_fixed(x, R, V):
    #We are working with an array, so for the if statements to work 
    #correctly we have to work element wise. 
    length = len(x)

    #Fixing the powers.
    p_in = 0.5
    p_out = 1.0
    
    #Creating the velocity array.    
    v = 0*x
    i = 0
    while i < length:
        if x[i] <= R:
            v[i] = V*(x[i]/R)**(-p_in)
        elif x[i] > R:
            v[i] = V*(x[i]/R)**(-p_out)
        else:
            print 'Error'
        i += 1
        #End while loop
    return v
    #End function


'''
