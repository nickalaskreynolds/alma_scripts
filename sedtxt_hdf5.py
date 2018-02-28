import pdspy.spectroscopy as sp
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file',default='sed.txt')
parser.add_argument('-o', '--output',default='sed.hdf5')
args = parser.parse_args()

data = np.loadtxt(args.file)
wave = data[:,0] #wavelength in microns
flux = data[:,1]/1000. #input in mJy convert to Jy
unc = 0.1*flux

sed = sp.Spectrum(wave,flux,unc)
sed.write(args.output)
