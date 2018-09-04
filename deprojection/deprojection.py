# this script will deproject the triplet source after the tertiary was removed
# need to do it this way since the continuum was made with mfs

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
    infile,pa,inc = args.infile,args.posa,args.inclin
    print('Running: ',infile)
    outfile = '.'.join(infile.split('.')[:-1]) + '_deproj.uvfits'
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

spws=[
['4','10','16'],
['3','9','15'],
['2','8','14'],
['0','6','12'],
['1','7','13'],
['5','11','17']]

names=[
'C17O',
'H13COp',
'H13CN',
'CO',
'SiO',
'335.5GHz']


# SUBTRACTION
sourcems = 'L1448IRS3B.cont.ms.apselfcal.concat.subclump'

class argsim(object):
    def __init__(self):
        self.infile  = ''
        self.posa   = 90-29.5  # 240,
        self.inclin = 45       # 68,
        self.dtype  = 'uvfits'

args = argsim()
dire = 'tempdepr'

allfiles = []
for x in range(len(spws)):
    for y in spws[x]:
        findir="{}/{}.{}".format(dire,names[x],y)
        if not os.path.isdir(findir):
            os.system('mkdir -p {}'.format(findir))
        else:
            os.system('rm -rf {}/*'.format(findir))
        print('Splitting')
        split(vis=sourcems,spw=y,datacolumn='data',width=999,keepflags=False,outputvis='{}/{}.{}.{}.ms'.format(findir,sourcems,names[x],y))
        print('Exporting')
        exportuvfits(vis='{}/{}.{}.{}.ms'.format(findir,sourcems,names[x],y),fitsfile='{}/{}.{}.{}.uvfits'.format(findir,sourcems,names[x],y))
        print('Modifying')
        args.infile = '{}/{}.{}.{}.uvfits'.format(findir,sourcems,names[x],y)
        print(args.infile)
        main(args)
        importuvfits(args.infile[:-7]+ '_deproj.uvfits',args.infile[:-7] + '_deproj.tmp')
        split(args.infile[:-7] + '_deproj.tmp',args.infile[:-7] + '_deproj',datacolumn='corrected')
        allfiles.append(args.infile[:-7]+ '_deproj')
print('Concat')
concat(vis=allfiles,concatvis='L1448IRS3B.cont.ms.apselfcal.concat.subclump.deproj')


# Imaging
contvis = sourcems

robusts = [-0.5,0.5]
tapers = [400,]
zipped = [[rob,tap]for rob in robusts for tap in tapers]
for x in zipped:
    robustp,taper = x
    image = contvis + '.robust'+str(robustp)+ '.taper'+str(taper)
    os.system('rm -rf '+image+'.*')

    tclean(
    vis = contvis,
    spw = '',
    field = '',
    scales = [0,5,15,20],
    phasecenter = '',
    usemask = 'auto-multithresh',
    weighting = 'briggs',
    robust = robustp,
    niter = 1000,
    threshold = '0.00008Jy',
    imagename = image,
    specmode = 'mfs',
    mask = '',
    gridder = 'standard',
    deconvolver = 'multiscale',
    imsize = 512,
    cell = '0.02arcsec',
    uvtaper = str(taper)+'klambda',
    interactive=False
    )
    print(image+'.image')
    break

