#!/usr/bin/env python
'''
wcs.py [--options] image.fits

This routine just adds a WCS to the images,
using astrometry-net to find coordinate solutions
'''


import os,sys,pdb,scipy,glob
import subprocess
from subprocess import check_output
software = os.path.dirname('/'.join(os.path.realpath(__file__).split('/')[:-1]))
sys.path.append(software)


sourcebindir = os.environ['SNPATH']
sourcedir = os.environ['SNSOURCE']

image = sys.argv[1]
radius = 0.17

if __name__=='__main__':
    import getopt
    try:
        opt,arg = getopt.getopt( 
            sys.argv[1:],"v,h",
            longopts=["verbose=","quiet","help","pixscale="])
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
        elif o == "--pixscale" :
            radius = float(a)

    if len(arg) < 1:
        print(__doc__)
        

    cmd = [ os.path.join( sourcebindir, "gethead")]
    cmd.append('%s' %image)
    cmd.append('RA DEC')
    print(' '.join(cmd))
    ra, dec = check_output(' '.join(cmd),shell=True).split()
    ra = str(ra,'utf-8')
    dec = str(dec, 'utf-8')
    cmd = [ os.path.join( sourcedir, "astrometry/bin/solve-field")]
    cmd.append('%s' %image)
    cmd.append('--ra %s --dec %s --radius %f'%(ra,dec,radius))
    print(' '.join(cmd))
    subprocess.call(' '.join(cmd),shell=True)
    outimg = image.split('.fits')[0]+'.new'
    if os.path.isfile(outimg):
        print('renaming %s --> %s' %(outimg,image))
        os.rename(outimg, image)
    print("\n\n DONE !! \n\n")
