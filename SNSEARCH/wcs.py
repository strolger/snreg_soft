#!/usr/bin/env python
'''
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

cmd = [ os.path.join( sourcebindir, "gethead")]
cmd.append('%s' %image)
cmd.append('RA DEC')
print(' '.join(cmd))
ra, dec = check_output(' '.join(cmd),shell=True).split()

pdb.set_trace()
cmd = [ os.path.join( sourcedir, "astrometry/bin/solve-field")]
cmd.append('%s' %image)
cmd.append('--ra %s --dec %s --radius %f'%(ra,dec,radius))
print(' '.join(cmd))
subprocess.call(' '.join(cmd),shell=True)
if os.path.isfile(image.split('.fits')[0]+'.new'):
    print('renaming .new image')
    os.rename(image.split('.fits')[0]+'.new', image)
print("\n\n DONE !! \n\n")







