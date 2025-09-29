#!/usr/bin/env python
# 2022
# L. Strolger
# Something simple to copy out the science extension
"""
get_sci_img.py img.fits [or @img.list]
"""


import os,sys,pdb,scipy,glob
from pylab import *
from astropy.io import fits
from astropy import wcs


def copy_out_sci(fitsimg):
    if not os.path.isfile(fitsimg):
        print ("%s doesn't exist... moving on..." %fitsimg)
        return()

    outimg=fitsimg.replace('.fits','_sci1.fits')
    if os.path.isfile(outimg): os.remove(outimg)

    hdul = fits.open(fitsimg)
    data = hdul['SCI'].data
    wcs1 = wcs.WCS(hdul['SCI'].header)
    
    hdu2 = fits.PrimaryHDU(data, header=wcs1.to_header())
    hduw = fits.HDUList([hdu2])
    hduw.writeto(outimg)
    print('%s written' %outimg)
    hdul.close()
    hduw.close()
    return()
    
if __name__ == '__main__':
    
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
        copy_out_sci(image)


        
        
