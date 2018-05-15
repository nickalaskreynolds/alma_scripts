import numpy as np
from math import *


print('########################')
jy = 1.e-23                  # Jansky CGS
k = 1.38064852E-16    #(cm^2)*g/(s^2)*K
kb=k
h = 6.626176E-27      #erg secs
c = 2.998E10           #cm/sec 
theta = 0.27*0.19     # arcsec^2
theta = theta*np.pi * (1./(60.*60.))**2 * (np.pi/180.)**2
nu = 337061.104*10**6  # Hz

######################
# input values
w = 945.8789 #590.146 # 945.8789 # 391.1658  # Jy* km/s
temp=25.

w = w * 1.0e5 # to Jy*cm/s
# finding temperature
t = (2.*np.log(2.)*c**2)/(k * nu**2)*jy
print('nu: ' + str(nu))
print('HPBW**2: ' + str(theta))
print('Jy to K: ' + str(t/theta))
print('Aper Flux: ', w/1E5)

w=w*t/theta

Tex=temp   # Excitation temperature of CO  in K
# finding column density
a = 2.32E-6     # s^-1
B = 56179.99E6  # Hz
z = k*Tex / (h * B)
print('Partition: ' + str(z))
j=3
g = 2*j + 1
mu = 0.11034
E = h*B*(j*(j+1))
# print(E/k)
den=w

n = (8. * np.pi * k * nu**2 * den )/(h * c**3 * a)

print('N(thin): ' + str(n))

lnN_tot = np.log(n/g) + np.log(z) + (E / (k * Tex))
N_tot=np.exp(lnN_tot)
print('log N(total): ' + str(lnN_tot))
print('N(total) (cm^-2): ' + str(N_tot))

N_tot_c17o = N_tot

#surface area of region (used 1.0" for fiducial radius, replace with exact area
#find N(CO)/N(H2) and N(C17O)/N(CO)

cell=0.02 # arcsec
area = 26856. * (cell**2)

A= area *(230.0*1.496e13)**2
#A=0.27*0.19*(distance*1.496e13)**2
print('Radius: ',(A/np.pi)**0.5/(1.496E13))
N_mol=N_tot*A
print('Area (cm^2): ' + str(A))
print('N(total) (molecules C17O): ' + str(N_mol))


abun_ratio_c17o_co = float(1./2000.)
abun_ratio_co_h2 = float(10**-4)

N_mol_h2 = N_mol/(abun_ratio_c17o_co * abun_ratio_co_h2)

print('N(total) (molecules H2): ' + str(N_mol_h2))

solmass = 1.988435E33   # grams per sol mass
mol_mass = 3.34E-24  # grams per h2 molecule
Mass_h2 = N_mol_h2 * mol_mass
mass = Mass_h2

###########

print('Mass (kg) (molecules H2): ' + str(Mass_h2))
print('Mass (Msol) (molecules H2): ' + str(Mass_h2/solmass))


# dust mass conversion

kb=1.36e-16
dgr=1.0

# Convert the two arguments from strings into numbers
kappa0= 0.899
beta = 1.75
lam0 = 1.3
lam = 0.88
flux = 1.33
distance = 230.0

distance=distance*3.09e18
temp=20.0
wave=np.array([500.0,600.0,700.0,800.0,870.0,900.0,1000.0,2000.0,3000.0,3400.0,7000.0,10000.0])
flux=flux*1.0e-23
nud=3.0e11/lam
print(nud/1.0e9, 'GHz')
#kappa0=3.5
kappa=kappa0*(lam0/lam)**beta
print(kappa, ' at ',lam,' mm')
dustmass=distance**2*flux*3e10**2/(2.0*kappa*kb*nud**2*temp)*dgr

print('Dust Mass: ', dustmass/2.0e33, 'M_sun')

gas2dust=Mass_h2/dustmass
print('Gas to Dust Mass Ratio: ',gas2dust)

###################################################################


###################################################################

print('########################')
print('Starting from the Tobin dataset')
# constants
jy = 1.e-23                  # Jansky CGS
k = 1.38064852E-16    #(cm^2)*g/(s^2)*K
kb=k
h = 6.626176E-27      #erg secs
c = 2.998E10           #cm/sec 
theta = 0.56*0.50     # arcsec^2
theta = theta*np.pi * (1./(60.*60.))**2 * (np.pi/180.)**2
nu = 219597.709*10**6  # Hz


w = 234.8819 # 945.8789 # 391.1658  # Jy* km/s
w = w * 1.0e5 # to Jy*cm/s
#JT's function for checking
#def jy_beam_to_k(bmaj=None,bmin=None,freq=None):
    #beams in degrees frequency in hertz
#    omega=(pi*bmaj*bmin)*(pi/180.0)**2
#    print('Omega: ' + str(omega))
#    wave=c/(freq)
#    jy2k=wave**2*jy/(2.0*kb*omega)*4.0*log(2.0)
#    return jy2k
#convfact=jy_beam_to_k(0.27/3600.0,0.19/3600.0,nu)
#print('Convfact: ' + str(convfact))

# finding temperature
t = (2.*np.log(2.)*c**2)/(k * nu**2)*jy
print('nu: ' + str(nu))
print('HPBW**2: ' + str(theta))
print('Jy to K: ' + str(t/theta))

Tex=temp   # Excitation temperature of CO  in K
# finding column density
a = 6.266E-8     # s^-1
B = 54891.42E6  # Hz
z = k*Tex / (h * B)
print('Partition: ' + str(z))
j=2
g = 2*j + 1
mu = 0.11049
E = h*B*(j*(j+1))

den=w

n = (8. * np.pi * k * nu**2 * den )/(h * c**3 * a)

print('N(thin): ' + str(n))

lnN_tot = np.log(n/g) + np.log(z) + (E / (k * Tex))
N_tot=np.exp(lnN_tot)
print('log N(total): ' + str(lnN_tot))
print('N(total) (cm^-2): ' + str(N_tot))

###################################################################

print('########################')
print('Starting from the Concat with range of Tobin dataset')
# range 1.25-4, 5.5,7

# constants
jy = 1.e-23                  # Jansky CGS
k = 1.38064852E-16    #(cm^2)*g/(s^2)*K
kb=k
h = 6.626176E-27      #erg secs
c = 2.998E10           #cm/sec 
theta = 0.27*0.19     # arcsec^2
theta = theta*np.pi * (1./(60.*60.))**2 * (np.pi/180.)**2
nu = 337061.104*10**6  # Hz


w = 181.9885 # 945.8789 # 391.1658  # Jy* km/s
w = w * 1.0e5 # to Jy*cm/s
#JT's function for checking
#def jy_beam_to_k(bmaj=None,bmin=None,freq=None):
    #beams in degreees frequency in hertz
#    omega=(pi*bmaj*bmin)*(pi/180.0)**2
#    print('Omega: ' + str(omega))
#    wave=c/(freq)
#    jy2k=wave**2*jy/(2.0*kb*omega)*4.0*log(2.0)
#    return jy2k
#convfact=jy_beam_to_k(0.27/3600.0,0.19/3600.0,nu)
#print('Convfact: ' + str(convfact))

# finding temperature
t = (2.*np.log(2.)*c**2)/(k * nu**2)*jy
print('nu: ' + str(nu))
print('HPBW**2: ' + str(theta))
print('Jy to K: ' + str(t/theta))

w=w*t/theta

Tex=temp  # Excitation temperature of CO  in K
# finding column density
a = 2.32E-6     # s^-1
B = 56179.99E6  # Hz
z = k*Tex / (h * B)
print('Partition: ' + str(z))
j=3
g = 2*j + 1
mu = 0.11034
E = h*B*(j*(j+1))

den=w

n = (8. * np.pi * k * nu**2 * den )/(h * c**3 * a)

print('N(thin): ' + str(n))

lnN_tot = np.log(n/g) + np.log(z) + (E / (k * Tex))
N_tot=np.exp(lnN_tot)
print('log N(total): ' + str(lnN_tot))
print('N(total) (cm^-2): ' + str(N_tot))

print('########################')

print('Ratio of C17O : C18O: ' + str(N_tot_c17o/N_tot))

print('########################')


'''
########################
Starting from the Concat dataset
nu: 3.37061104e+11
HPBW**2: 3.78806106562e-12
Jy to K: 209.701401373
Partition: 11.1265304734
N(thin): 5.61634847822e+16
log N(total): 40.1089687588
N(total) (cm^-2): 2.62484586952e+17
Area (cm^2): 1.27180397685e+32
N(total) (molecules C17O): 3.33828941547e+49
N(total) (molecules H2): 6.67657883093e+56
Mass (kg) (molecules H2): 2.22997732953e+33
Mass (Msol) (molecules H2): 1.12147358578
(340.9090909090909, 'GHz')
(1.7795733035310959, ' at ', 0.88, ' mm')
('Dust Mass: ', 0.0026868587679850847, 'M_sun')
('Gas to Dust Mass Ratio: ', 414.97851619558611)
########################
Starting from the Tobin dataset
nu: 2.19597709e+11
HPBW**2: 2.06755769664e-11
Jy to K: 90.5153707241
Partition: 11.3877245429
N(thin): 3.5130142331e+17
log N(total): 41.7504022084
N(total) (cm^-2): 1.35509360684e+18
########################
Starting from the Concat with range of Tobin dataset
nu: 3.37061104e+11
HPBW**2: 3.78806106562e-12
Jy to K: 209.701401373
Partition: 11.1265304734
N(thin): 1.73196207148e+16
log N(total): 38.9325419543
N(total) (cm^-2): 8.09446477034e+16
########################
Ratio of C17O : C18O: 3.24276643854
########################
'''
