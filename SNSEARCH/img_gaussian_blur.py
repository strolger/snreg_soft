#!/usr/bin/env python
# 2023
# L. Strolger
# Something simple to gaussian blur HST images
"""
img_gaussian_blur.py img.fits [or @img.list]
"""

import os,sys,pdb,scipy,glob
from pylab import *
from astropy.io import fits
from scipy.ndimage import gaussian_filter

sigma = 0.5

def gaussian_blur(fitsimg, sigma=sigma):
    if not os.path.isfile(fitsimg):
        print ("%s doesn't exist... moving on..." %fitsimg)
        return()

    outimg=fitsimg.replace('.fits','_blur.fits')
    if os.path.isfile(outimg): os.remove(outimg)


    hdul = fits.open(fitsimg)
    data = hdul[0].data
    result = gaussian_filter(data, sigma=sigma)
    hdu2 = fits.PrimaryHDU(result)
    hduw = fits.HDUList([hdu2])
    hduw.writeto(outimg)
    print('%s written' %outimg)
    hdul.close()
    hduw.close()
    return()
    
    

if __name__=='__main__':
    
    img = sys.argv[1]

    imglist = []
    if img.startswith('@'):
        filename = img.strip('@')
        f = open(filename, 'r')
        lines = f.readlines()
        f.close()
        for line in lines:
            imglist.append(line.rstrip())
    else:
        imglist.append(img)

    for image in imglist:
        gaussian_blur(image)
