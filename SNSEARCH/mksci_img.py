#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from astropy.io import fits

if __name__=='__main__':

    img = sys.argv[-1]
    out = img.replace('.fits','_sci.fits')
    hduin = fits.open(img)
    data = hduin['sci'].data
    hduin.close()
    hduout = fits.PrimaryHDU(data)
    hduout.writeto(out)
    #hduout.close()
                
    print('Done')
    
