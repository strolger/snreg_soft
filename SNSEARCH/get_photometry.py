#!/usr/bin/env python
'''
get_photometry.py [--options] image.fits

This routine generates a sextractor catalog
with the appropriate zero point to get magnitudes
and effective radii
output is image.fits.phot file
'''


import os,sys,pdb,scipy,glob
import subprocess
from subprocess import check_output
import sextract

software = os.path.dirname('/'.join(os.path.realpath(__file__).split('/')[:-1]))
sys.path.append(software)


sourcebindir = os.environ['SNPATH']
sourcedir = os.environ['SNSOURCE']

flagimg = None
zp= None
verbose = True

if __name__=='__main__':
    import getopt
    try:
        opt,arg = getopt.getopt( 
            sys.argv[1:],"v,h",
            longopts=["verbose=","quiet","help","flagimg=","zp="])
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
        elif o == "--flagimg":
            flagimg = a
        elif o == "--zp" :
            zp = float(a)

    if len(arg) < 1:
        print(__doc__)

    image = arg[0]
    if image.endswith('.fits'): image=image.strip('.fits')
    if not os.path.isfile(image+'.fits'):
        print('there is an error \n filename %s does not exist\n' %image+'.fits')
        print(__doc__)
        pdb.set_trace()
        sys.exit(1)
    if verbose > 1: print('image:\t%s.fits' %image)
       
    paramstring = ""
    paramstring +="-CATALOG_NAME %s.fits.phot" %image
    if flagimg:
        paramstring += " -PARAMETERS_NAME %s -DETECT_THRESH 2. -PHOT_APERTURES 10.0" %(os.path.join(sourcedir,"Sex/deep_find_flags.param"))
    else:
        paramstring += " -PARAMETERS_NAME %s -DETECT_THRESH 1.5 -PHOT_APERTURES 10.0" %(os.path.join(sourcedir,"Sex/deep_find.param"))
    if zp:
        if verbose: print("Assuming zeropoint is %2.1f" %zp)
        paramstring += " -MAG_ZEROPOINT %2.2f" %(zp)
    else:
        if verbose: print("Assuming no zeropoint")

    sextract.runsex(image+'.fits',verbose=verbose,
                    flagfile=flagimg, clobber=True,
                    paramstring=paramstring)
        
