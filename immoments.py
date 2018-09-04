#!/bin/env python

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from sys import version_info,exit
from collections import Iterable
assert version_info[0] == 3

def typecheck(obj): 
    return not isinstance(obj, str) and isinstance(obj, Iterable)

def reading(fname):
    with fits.open(fname) as hdu:
        header = hdu[0].header
        data   = hdu[0].data
    return data,header

def saving(fname,data,header):
    fits.writeto(fname,data,header)
    pass

def immoment(data,header=None, moment=-1,axis='channel',lval=-1,uval=-1):
    if lval == -1: lval = 0
    if uval == -1: uval = data.shape[1]
    assert uval > lval
    if   moment == -1:
        momentdata = np.mean(data[0,lval:uval,:,:],axis=0)
    elif moment == 0 :
        momentdata = np.sum(data[0,lval:uval,:,:],axis=0)
    else:
        print('Other moments not supported')
        exit(0)
    return momentdata

def plotter(data,header,cmap='viridis',save=False):
    title,xlabel,ylabel,zlabel=header['OBJECT'],header['CTYPE1'],header['CTYPE2'],header['BUNIT']
    plt.figure(figsize=[20,10])
    plt.imshow(data)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.colorbar()
    plt.tight_layout()
    plt.draw()
    if save: plt.savefig(save,dpi=400)
    plt.show()

def main(args):
    data,header = reading(args.input)
    if (type(args.moments) == int) or (type(args.moments) == float):
        moment=immoment(data,header,moment=args.moments,lval=args.lval,uval=args.uval)
        plotter(moment,header,save=args.save)
    else:
        args.moments = [int(x) for x in args.moments.strip('[').strip(']').split(',')]
        for x in args.moments:
            moment=immoment(data,header,moment=x,lval=args.lval,uval=args.uval)
            plotter(moment,header,save=args.save)


if __name__ == "__main__":
    from argparse import ArgumentParser as ap
    parser = ap()
    parser.add_argument('-i', '--input',help='Input fits file name',default='image.fits',required=True)
    parser.add_argument('-s', '--save',help='Output pdf file name')
    parser.add_argument('-o', '--output',help='Output fits file name',default='cores.fits')
    parser.add_argument('-m', '--moments',help='Integer or list of moment maps to compute',default=-1)
    parser.add_argument('--lval',help='Lower channel value',default=-1)
    parser.add_argument('--uval',help='Upper channel value',default=-1)
    args = parser.parse_args()

    
    main(args)

# end of code
