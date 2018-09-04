#!/usr/bin/env python
# Standard Modules
import argparse
import glob
from time import time
TIME=str(time)

# Custom Modules
import pdspy.interferometry as uv

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',dest='inp',required=True)
parser.add_argument('-of','--outfits',dest='of',default=False,action="store_true")
parser.add_argument('-fc', '--freqcor',dest='fc',default=False,action="store_true")
parser.add_argument('-o', '--output',dest='out',required=True)
parser.add_argument('-umi', '--uvmin',dest='uvmin',type=int,default=50000)
parser.add_argument('-uma', '--uvmax',dest='uvmax',type=int,default=500000)
args = parser.parse_args()

if not args.of:
    print('Using Directory: {}'.format(args.inp))
    if not args.fc:
        files=glob.glob('{}'.format(args.inp))
    else:
        files=glob.glob('{}/*'.format(args.inp))
    print('Found: {}'.format(files))

    data = []
    for file in files:
        if args.fc:
            data.append(uv.freqcorrect(uv.readvis(file)))
        else:
            data.append(uv.readvis(file))

    data[-1].weights[data[-1].uvdist < args.uvmin] = 0.
    data[-1].weights[data[-1].uvdist > args.uvmax] = 0.
    data[-1].weights[data[-1].weights < 0.0] = 0.
    vis = uv.concatenate(data)

    print('Writing to file: {}'.format(args.out))
    vis.write(args.out)

elif args.of:
    from fits2hdf.io.fitsio import export_fits
    from fits2hdf.io.hdfio import read_hdf
    a = read_hdf(args.inp)
    export_hdf(a, args.out)
