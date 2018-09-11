#!/usr/bin/env python
'''
This routine creates a mask image for a fits image that excludes
vignetted regions
requires astropy

'''


import os,sys,pdb,scipy,glob
from pylab import *
from astropy.io import fits
from stsci.convolve import boxcar
verbose = True



if __name__=='__main__':
    
    filename = sys.argv[1] ## a bias-subtracted, combined flat_field is best
    hdu = fits.open(filename)
    data = hdu[0].data
    if verbose: print('Replacing data region with zeros and ones...')
    tmp = data > 25000.
    data = tmp.astype(int)
    if verbose: print('Boxcar smoothing...')
    data = boxcar(data.astype(float),(25,25))
    if verbose: print('Flipping array... ')
    tmp = data == 0.
    data = tmp.astype(int)
    if verbose: print('Boxcar smoothing again...')
    data = boxcar(data.astype(float),(51,51))
    if verbose: print('Replacing positive values to create mask...')
    tmp = data > 0.
    data = tmp.astype(int)
    if verbose: print('Writing mask image flagimg.fits...')
    out = fits.PrimaryHDU(data)
    out.writeto('flagimg.fits')
    if verbose: print('Done.')
    
    
