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
import subprocess

## import sextract

sigma =0.5

def gaussian_blur(fitsimg, sigma=sigma):
    if not os.path.isfile(fitsimg):
        print ("%s doesn't exist... moving on..." %fitsimg)
        return()

    origimg = None
    if fitsimg.endswith('.fz'):
        print('%s is compressed with fpack... uncompressing' %fitsimg)
        try:
            result = subprocess.run(["funpack","-v",fitsimg])
        except:
            print('problem with fpack... check if exists')
            sys.exit(1)
        origimg=fitsimg
        fitsimg=fitsimg.strip('.fz')

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
    if origimg:
        print('recompressing %s is compressed with fpack' %fitsimg)
        try:
            result = subprocess.run(["fpack","-v",fitsimg])
        except:
            print('problem with fpack... check if exists')
            sys.exit(1)
        
    return()
    
    

if __name__=='__main__':
    import getopt
    try:
        opt,arg = getopt.getopt( 
            sys.argv[1:],"v,h",
            longopts=["verbose=","quiet","help","sigma="])
    except getopt.GetoptError: 
        print ("Error : incorrect option or missing argument.")
        print (__doc__)
        sys.exit(1)
    for o, a in opt:
        if o in ["-h", "--help"]:
            print (__doc__)
            sys.exit(0)
        elif o == "-v" :
            verbose = 2
        elif o == "--quiet" :
            verbose = False
        elif o == "--verbose" :
            try: verbose = int(a)
            except ValueError : 
                print("Error : must specify verbosity level as an integer 1-10.")
                sys.exit(0)
        elif o == "--sigma" :
            sigma = float(a)

    if len(arg) < 1:
        print(__doc__)

    
    img = arg[0]
    
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
        gaussian_blur(image, sigma=sigma)
        ### the following is a kludge for remaking the catalogs for panethon
        ## paramstring='-PARAMETERS_NAME pantheon.param -c pantheon.sex -VERBOSE_TYPE Quiet '
        ## paramstring+='-NTHREADS 1 '
        ## if image.endswith('.fz'):
        ##     imfile = image.replace('.fz','')
        ## else:
        ##     imfile = image
        ## blurimg = imfile.replace('.fits','_blur.fits')
        ## sextract.runsex(blurimg, catfile=None, weightfile=None, flagfile=None,
        ##                 paramstring=paramstring, detimfile=None, verbose=1,
        ##                 clobber=True)
        
    
        
