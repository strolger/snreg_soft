#!/usr/bin/env python
# 2022
# L. Strolger
# Something simple to subtract images
"""
imadd.py template.fits image.fits
"""


import os,sys,pdb,scipy,glob
from pylab import *
from astropy.io import fits


def coadd(tmp, img):
    if not os.path.isfile(tmp):
        print ("%s doesn't exist... moving on..." %tmp)
    elif not os.path.isfile(img):
        print ("%s doesn't exist... moving on..." %img)
    else:
        outimg = img.replace('.fits','_coadd.fits')
        if os.path.isfile(outimg): os.remove(outimg)
        hdu1 = fits.open(tmp)
        try:
            dat1 = hdu1['SCI'].data
        except:
            dat1 = hdu1[0].data

        hdu2 = fits.open(img)
        try:
            dat2 = hdu2['SCI'].data
        except:
            dat2 = hdu2[0].data
            

        dato = 0.5*(dat2 + dat1)
        
        hduo = fits.PrimaryHDU(dato)
        hduw = fits.HDUList([hduo])
        hduw.writeto(outimg)
        print('%s written' %outimg)
        hdu1.close()
        hdu2.close()
        hduw.close()
    return()
    
if __name__ == '__main__':
    
    tmp = sys.argv[1]
    img = sys.argv[2]
    coadd(tmp, img)
    
    

        
        
