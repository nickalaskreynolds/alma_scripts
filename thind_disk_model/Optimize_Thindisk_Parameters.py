import numpy as np
import pyfits as pf
import os
from Merit_Function import *
from scipy.optimize import minimize

def Random_Number(Lower_Bound, Upper_Bound):
    #This function creates an random number between the Lower and Upper bound.
    #With Lower and Upper bound Nx1 arrays. 
    length_array = len(Lower_Bound)
    return (Upper_Bound-Lower_Bound)*np.random.rand(length_array) + Lower_Bound

#-------------------------------------------------------------------------------
#User input
#-------------------------------------------------------------------------------
#Creating the directory for saving the images, if it does not already exists.
save_directory = 'Optimization_Results/Round_7/Data/' 
if not os.path.exists(save_directory):
    os.makedirs(save_directory)
    
#Amount of optimizations
N = 200


#The inclination
inclination = [85]#degrees

#How the velocity profile of the infaling material is.
#material_v = 'infalling'
#material_v = 'rotating'
#material_v = 'rotating_only_in_rc'
material_v = 'rotating+infalling'
#material_v = 'fixed_size'
#-------------------------------------------------------------------------------
#Optimizing the parameters.
#-------------------------------------------------------------------------------
j = 0
while j < N:
    for i in inclination:
        if material_v == 'rotating_only_in_rc' or material_v == 'fixed_size':
            #The bounds for the random starting position.
            #Mass, rc, I0, fwhm
            Lower_Bound = np.array([0.1,1,1,1])
            Upper_Bound = np.array([1,100, 10, 10])
            #Creating the random starting position.
            p0 = Random_Number(Lower_Bound, Upper_Bound)
        else:
            #The bounds for the random starting position.
            #mass, rc, size, I0, FWHM
            Lower_Bound = np.array([0.1,1,100,0,1])
            Upper_Bound = np.array([1,100, 1000, 1, 3])
            #Creating the random starting position.
            p0 = Random_Number(Lower_Bound, Upper_Bound)
        
        print 'Optimization: ' + str(j + 1)
        print 'Inclination: ' + str(i) + ' degrees'
        print 'Velocity profile: ' + material_v
        print 'p0: ' + str(p0)
        
        #Fitting the model to the data.
        result = minimize(Merit_FunctionThindisk, p0, args = (i, material_v), method='powell', options = {'disp': True})


        #The current time and data. 
        import datetime
        k = datetime.datetime.now()

        print 'Saving the result.'
        #Saving the result.
        Save_File = np.zeros((len(p0), 2))
        Save_File[:,0] = p0
        Save_File[:,1] = result.x
        np.savetxt(save_directory +  'i_' + str(i) + '_' + material_v + '_Results_Value_'+ str(result.fun) + '_Data_' +  k.isoformat() + '.txt', Save_File, delimiter = ' ', header = 'p0                       x')

    j += 1

