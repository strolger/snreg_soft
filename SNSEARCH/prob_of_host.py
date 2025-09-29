#!/usr/bin/env python
'''
Input:
1. FITS image
2. Coordinates of SN on image (in ra,dec u.deg)

Looking to get
1. Probability of host for any nearby sources
2. d_DLR for same nearby sources


'''
import os,sys,pdb,scipy,glob
import pandas as pd
from pylab import *
## from strolger_util import util as u
from scipy.integrate import quad

## snsearch_path = os.environ['SNSEARCH']
## sys.path.append(snsearch_path)
## import sextract

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord




def read_sexcat(filename, cmtchar='#'):
    col_names =[]
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith(cmtchar):
                col_names.append(line.split()[2])
    df = pd.read_csv(filename, comment=cmtchar, header=None, sep=r'\s+')
    df.columns=col_names
    return(df)

def make_regfile(regfile,xx, yy, df,nclose=6):
    ## select only the sources with d_DLR < 5
    df=df[df['d_DLR']<4]
    
    f=open(regfile, 'w')
    f.write('global color=green dashlist=8 3 width=1 ')
    f.write('font="helvetica 10 normal roman" select=1 ')
    f.write('highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write('image\n')
    for ii, rr in df.iloc[:nclose].iterrows():
        f.write('ellipse(%d,%d,%d,%d,%d)\n' %(rr['X_IMAGE'],
                                              rr['Y_IMAGE'],
                                              rr['A_IMAGE']*rr['KRON_RADIUS'],
                                              rr['B_IMAGE']*rr['KRON_RADIUS'],
                                              rr['THETA_IMAGE'])
                )
        f.write('# text(%.1f,%.1f) textangle=0, text={%d}\n'
                %(rr['X_IMAGE']+10, rr['Y_IMAGE']+10, rr['NUMBER']))
        f.write('point (%.1f, %.1f) # point=cross color=red\n' %(xx, yy))
    f.close()
    return()

def dv(x,*p):
    k,n=p
    I_Ie=exp(-k*(x**(1./n)-1.0))
    return(I_Ie)

def prob_of_gxy(offset,n):
    p0 = [7.67, n]
    frac_of_light = quad(dv, 0, offset, args=tuple(p0))[0]/quad(dv, 0, 1e3, args=tuple(p0))[0]
    return(frac_of_light)
    
if __name__=='__main__':

    imfile = sys.argv[1]
    ## xx = float(sys.argv[2])
    ## yy = float(sys.argv[3])

    rr = float(sys.argv[2])
    dd = float(sys.argv[3])
 
    with fits.open(imfile) as hdul:
        # Get the WCS object
        w = WCS(hdul[1].header)
        coords = SkyCoord(rr, dd, unit='deg')
        # Convert to pixel coordinates
        xx, yy = w.world_to_pixel(coords)


    clobber = True
    
    catfile = imfile.replace('.fits','.sexcat')
    regfile = imfile.replace('.fits','.reg')

    
    if not os.path.isfile(catfile) or clobber:
        sextract.runsex(imfile, verbose=True, clobber=clobber)
    
    df = read_sexcat(catfile)

    ## the easy way
    ## df['OFFSET'] = sqrt((xx-df['X_IMAGE'])**2+(yy-df['Y_IMAGE'])**2)
    ## df['RE'] = 0.5*df['KRON_RADIUS']*(df['A_IMAGE']+df['B_IMAGE'])
    ## df['d_DLR']=df['OFFSET']/df['RE']

    ## the more correct way
    to_rad = pi/180.0
    Dx = xx-df['X_IMAGE']
    Dy = yy-df['Y_IMAGE']
    df['OFFSET']=sqrt(Dx**2+Dy**2)
    phi = arctan(Dy/Dx)-to_rad*df['THETA_IMAGE'] ##in rad
    xp = df['A_IMAGE']*cos(phi)*cos(to_rad*df['THETA_IMAGE'])-df['B_IMAGE']*sin(phi)*sin(to_rad*df['THETA_IMAGE'])
    yp = df['A_IMAGE']*cos(phi)*sin(to_rad*df['THETA_IMAGE'])-df['B_IMAGE']*sin(phi)*cos(to_rad*df['THETA_IMAGE'])
    df['RE']= df['KRON_RADIUS']*sqrt(xp**2.+yp**2.)
    df['d_DLR']=df['OFFSET']/df['RE']


    df=df.sort_values(by='OFFSET', ignore_index=True)

    if not os.path.isfile(regfile) or clobber:
        make_regfile(regfile, xx, yy, df)

    ## pdb.set_trace()
    ## next to add sersic prob of host

    for nn in list([1,2,4]):
        f_hosts=[]
        for item in list(df['d_DLR']):
            tmp = 1.0-round(prob_of_gxy(item, n=nn), 5)
            if tmp < 0: pdb.set_trace()
            f_hosts.append(tmp)
        df['F_HOST_%d'%nn]=f_hosts
    ## pdb.set_trace()
    df.to_csv('probable_host.csv', index=False)
    
        
        
