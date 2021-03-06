
import astropy.units as u
import astropy.constants as co
import scipy as _sp
import numpy as np
import astropy.io.fits as fits
import sys
from argparse import ArgumentParser


########################################################################
def deg2rad(x):
    return x*np.pi/180.0
    
# DE-PROJECT i.e. ROTATE, INCLINE

def deproject(uv, PA, inc):
    """
    Rotate and deproject individual visibility coordinates.
    From Hughes et al. (2007) - "AN INNER HOLE IN THE DISK AROUND 
    TW HYDRAE RESOLVED IN 7 mm DUST EMISSION".
    """
    R = ( (uv**2).sum(axis=0) )**0.5
    #~ phi = _sp.arctan(uv[1]/uv[0] - deg2rad(PA))
    phi = _sp.arctan2(uv[1],uv[0]) - deg2rad(PA)
    #~ phi = _sp.arctan2( (uv[1] - deg2rad(PA) * uv[0]) , uv[0])
    newu = R * _sp.cos(phi) * _sp.cos( deg2rad(inc) )
    newv = R * _sp.sin(phi)
    newuv = _sp.array([newu, newv])
    ruv = (newuv**2).sum(axis=0)**.5
    return newuv, ruv

def rotate_field(uv, PA, U_RA_align = True):
    """
    Rotates a coordinate system (UV plane) by PA
    degrees.
    uv : 2-D array with uv[0] U and uv[1] coordinated
    PA : Position Angle, in degrees
    U_RA_align : for ALMA and PdBI the U-axis and RA are aligned
                 and thus one form of the equation must be used
                 While for SMA/CARMA (USA, meh), they are not aligned
                 and thus some sign changes have to impliemented from
                 that presented in Berger & Segransan (2007)
    
    """
    direction =  [-1, 1][int(U_RA_align)]
    u_new = uv[0] * _sp.cos( deg2rad(PA) ) + direction * uv[1] * _sp.sin( deg2rad(PA) )
    v_new = -1 * direction * uv[0] * _sp.sin( deg2rad(PA) ) + uv[1] * _sp.cos( deg2rad(PA) )
    return u_new, v_new
    
def incline(uv, inc):
    #~ ruv = ( uv[0]**2 + (uv[1] * _sp.cos(deg2rad(inc)) )**2  )**.5 
    # the PA goes from North to East in the image plane, and the 
    # Major axis is flipped 90 degrees going from
    # image to UV plane (Major in image in minor in UV)
    ruv = ( uv[0]**2 * _sp.cos(deg2rad(inc))**2 + uv[1]**2  )**.5 
    return ruv


def main(args):
    '''
    args is a dictionary of key values that are defined below
    '''
    print('Running: ',args.infile)
    outfile = '.'.join(args.infile.split('.')[:-1]) + '_deproj.uvfits'
    print('Writing to: ',outfile)
    # get the data and restfreq
    print('Reading input UV-FITS data: ' + infile)
    hdu = fits.open(infile)
    #restfreq = hdu[0].header['RESTFRQ']*u.Hz
    try:
        uu = hdu[0].data['UU']
        vv = hdu[0].data['VV']
    except IndexError:
        print('Export the MS to uv fits .')
    # deproject the coordinates
    print('Deprojecting data.')
    uuvv_new, ruv_new = deproject(np.array([uu,vv]), pa, inc)
    hdu[0].data['UU'] = uuvv_new[0]
    hdu[0].data['VV'] = uuvv_new[1]
    print('Saving output file: ' + outfile)
    hdu.writeto(outfile)


if __name__ == '__main__':

    description = 'Deprojects visibility data'\
                  'can also handle comma separate values'

    in_help = 'name of the file to deproject (FITS FORMAT) include extension'
    pa_help = 'Position angle of the source'
    d_help = 'Allow data of formats uvfits and hdf5'
    inc_help = 'Inclination of the source'

    # Initialize instance of an argument parser
    parser = ArgumentParser(description=description)

    # Add required arguments
    parser.add_argument('--input', help=in_help, type=str,dest='infile',required=True)
    parser.add_argument('-p', '--pa', type=str, help=pa_help,dest='posa',required=True)
    parser.add_argument('-i', '--inc', type=str, help=inc_help,dest='inclin',required=True)
    parser.add_argument('-dtype', type=str, help=d_help,default='uvfits')

    # Get the arguments
    args = parser.parse_args()
    assert args.dtype in ['uvfits','hdf5']
    main(args)