#!/usr/bin/env python
'''
get_zeropoint.py [--options] image.fits
'''



import os,sys,pdb,scipy,glob
from pylab import *
import subprocess
from subprocess import check_output
import sextract, imregister


verbose = 1
clobber = True
thresh = 50
miss = 3.0
zpR = 0.0
zpB = 0.0
scale = 0.0
tol = 0.005
flagimg = None
goodmatch = 5
nstars = 50
amiss = 10
maxstars = 900
matchsets = 3
order = 3
usewcs = True

software = os.path.dirname('/'.join(os.path.realpath(__file__).split('/')[:-1]))
sys.path.append(software)

try:
    sourcebindir = os.environ['SNPATH']
except KeyError :
    sys.exit(1)

arcsec_per_radian = 206264.8062470964
def s2d(ra,dec):
    if not isinstance(ra,str):
        print ("Expecting RA in sexigesimal")
        raise
    if not isinstance(dec,str):
        print ("Expecting DEC in sexigesimal")
        raise
    try:
        ra=list(map(float,ra.split(':')))
    except:
        print ("Could not convert RA, not in the right format?")
        raise
    try:
        dec=list(map(float,dec.split(':')))
    except:
        print ("Could not convert DEC, not in the right format?")
        raise
    rd=(ra[0]+ra[1]/60.+ra[2]/3600.)*360./24.
    if dec[0] < 0:
        dd=dec[0]-dec[1]/60.-dec[2]/3600.
    else:
        dd=dec[0]+dec[1]/60.+dec[2]/3600.
    return(rd,dd)

def radecstring2rad(inra, indec):
    rd, dd = s2d(inra,indec)
    decrad = dd*3600./arcsec_per_radian
    rarad = rd *3600./arcsec_per_radian
    return(rarad,decrad)


def eq2st(plate_center_ra, plate_center_dec, object_ra, object_dec):
    ## Compute std coord offsets (arcsec) given RA and Dec and plate center in rad.
    div = (sin(object_dec)*sin(plate_center_dec)+cos(object_dec)*cos(plate_center_dec)*
           cos(object_ra - plate_center_ra))
    xi_obj = cos(object_dec)*sin(object_ra-plate_center_ra)*arcsec_per_radian/div
    eta_obj = (sin(object_dec)*cos(plate_center_dec)-cos(object_dec)*sin(plate_center_dec)*
               cos(object_ra - plate_center_ra))*arcsec_per_radian/div
    return(xi_obj,eta_obj)
    


def binmode(data,bins=None):
    from scipy import stats
    data=array(data)
    mdx = where(~isnan(data))
    data = data[mdx]
    if bins != None:
        m,mbins = histogram(data, bins=bins)
    else:
        step=1/100.
        splits = arange(0,1+step,step)
        bin_edges=stats.mstats.mquantiles(data,splits)
        bins=sort(list(set(bin_edges)))
        #warnings.simplefilter("error",RuntimeWarning)
        if (max(bins[:-1])-min(bins[1:]))*step!=0.0:
            rebins = arange(min(bins[1:]), max(bins[:-1]),
                            (max(bins[:-1])-min(bins[1:]))*step)
            m,mbins = histogram(data,bins=rebins)
        else:
            m=zeros(1); mbins=zeros(2)
    
    mdx = where(m == max(m))
    mbin =0.5*(mbins[mdx[0][0]+1]+mbins[mdx[0][0]])
    return (mbin,
            array(zip(*vstack([m,mbins[:-1]])[::-1])))

if __name__=='__main__':
    import getopt
    try:
        opt,arg = getopt.getopt( 
            sys.argv[1:],"v,h",
            longopts=["verbose=","quiet","help","flagimg=","thresh=","scale=",
                      "force","useiraf","order="])
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
        elif o == "--thresh" :
            thresh = float(a)
        elif o == "--scale" :
            scale = float(a)
        elif o == "--force":
            clobber = True
        elif o == "--flagimg":
            flagimg = a
        elif o == "--useiraf":
            useiraf=True
        elif o == "--order":
            order = int(a)

    if len(arg) < 1:
        print(__doc__)

    
    # first arg is the image file
    image = arg[0]
    if image.endswith('.fits'): image=image.strip('.fits')
    if not os.path.isfile(image+'.fits'):
        print('there is an error \n filename %s does not exist\n' %image+'.fits')
        print(__doc__)
        pdb.set_trace()
        sys.exit(1)
    if verbose > 1: print('image:\t%s.fits' %image)
    if not os.path.isfile(image+'.sexcat') or not os.path.isfile(image+'.fits.stars') or clobber:
        sextract.runsex(image+'.fits',verbose=verbose,
                        flagfile=flagimg, clobber=True,
                        paramstring="-DETECT_THRESH %.2f -PHOT_APERTURES 10.0" %thresh)
        f = open(image+'.sexcat')
        lines = f.readlines()
        f.close()
        w = open(image+'.fits.stars','w')
        for line in lines:
            if line.startswith('#'): continue
            if flagimg:
                if int(line.split()[-2])==0:
                    w.write(line)
                else:
                    continue
            else:
                w.write(line)
        w.close()

    #Get header info
    cmd = [ os.path.join( sourcebindir, "gethead -u")]
    cmd.append('%s.fits' %image)
    cmd.append('RA DEC EQUINOX NAXIS1 NAXIS2')# CTYPE1 CTYPE2')
    print(' '.join(cmd))
    ra, dec, epoch, xaxis, yaxis = check_output(' '.join(cmd),shell=True).split()
    ra = str(ra,'utf-8')
    dec = str(dec, 'utf-8')
    epoch = str(epoch, 'utf-8')
    xaxis = int(xaxis)
    yaxis = int(yaxis)
    print(ra, dec, epoch, xaxis, yaxis)
    rar, decr = radecstring2rad(ra,dec)
    print('Read from Header: %.5f, %.5f' %(rar,decr))

    if usewcs: #use WCS to get center RA and DEC
        epoch = 200.0
        midx = int(xaxis)/2.
        midy = int(yaxis)/2.
        cmd = [ os.path.join( sourcebindir, "xy2sky")]
        cmd.append('%s.fits' %image)
        cmd.append('%f %f' %(midx,midy))
        print(' '.join(cmd))
        ra,dec,j1,j2,j3 = check_output(' '.join(cmd),shell=True).split()
        ra = str(ra,'utf-8')
        dec = str(dec, 'utf-8')
        print('%s %s' %(ra,dec))
        rar2, decr2 = radecstring2rad(ra, dec)
        midx+=100 # offset by 100 pixels to get scale
        midy+=100
        cmd = [ os.path.join( sourcebindir, "xy2sky")]
        cmd.append('%s.fits' %image)
        cmd.append('%f %f' %(midx,midy))
        print(' '.join(cmd))
        ra1,dec1,j1,j2,j3 = check_output(' '.join(cmd),shell=True).split()
        ra1 = str(ra1,'utf-8')
        dec1 = str(dec1, 'utf-8')
        print('%s %s' %(ra1,dec1))
        rar1, decr1 = radecstring2rad(ra1, dec1)
        dtheta,dphi = eq2st(rar2,decr2,rar1,decr1)
        if not scale:
            scale = sqrt(dtheta**2.+dphi**2.)/sqrt(20000.)
        print('Read from WCS: %.2f arcsec/pixel' %scale)
    size1 = xaxis*scale/sqrt(2.)
    size2 = yaxis*scale/sqrt(2.)
    if size1 > size2:
        gscsize = size1
    else:
        gscsize = size2



    OKmatch = sqrt(scale**2+miss**2.)

    if scale == 0.0:
        gscsize = max(xaxis,yaxis)*0.5*1.2 ## just a guess
        OKmatch = 50.

    if os.path.isfile('%s.fits.catalog' %image): os.remove('%s.fits.catalog' %image)
    cmd = [ os.path.join( sourcebindir, "scat")]
    cmd.append('-c ub1 -r %f -n 1000 -j %s %s %s' %(gscsize,ra,dec,epoch))
    cmd.append('| sort -k 5 >> %s.fits.catalog' %image)
    print(' '.join(cmd))
    subprocess.call(' '.join(cmd), shell=True)
        

    ## make a magnitude sorted list on tangent plane
    ## I'll fix the tangent plane stuff later
    if os.path.isfile('%s.fits.xieta' %image):os.remove('%s.fits.xieta' %image)
    f = open('%s.fits.xieta' %image,'w')
    g = open('%s.fits.catalog' %image,'r')
    lines = g.readlines()
    g.close()
    cnt =0
    ii = 0
    for line in lines:
        if line.startswith('#'): continue
        c = line.rstrip().split()
        ra1, dec1, magB, magR = c[1], c[2], float(c[4]), float(c[6])

        cmd = [ os.path.join( sourcebindir, "sky2xy")]
        cmd.append('%s.fits' %image)
        cmd.append('%s %s' %(ra1,dec1))
        ## print(' '.join(cmd))
        c = check_output(' '.join(cmd),shell=True).split()
        x, y = float(c[4]),float(c[5])

        ## rar1, decr1 = radecstring2rad(ra1,dec1)
        ## x,y = eq2st(rar, decr, rar1, decr1)
        ## x=x+xaxis/2.
        ## y=y+yaxis/2.

        ## rar1 *=180./pi
        ## decr1 *=180./pi
        cnt += 1
        if (magB  < 20. and x > 0 and y > 0 and x < float(xaxis) and y < float(yaxis)):
            f.write('%f %f %d %f %f %f %f\n' %(x,y,cnt, magB, magR, rar1, decr1))
    f.close()
    
    imregister.runmatch(image+'.fits.xieta',image+'.fits.stars',image+'.match',
                        order=order, scale = 0.0, tol=tol, maxstars=maxstars,
                        amiss=10, cbmiss = 0.01)
    data1 = loadtxt(image+'.fits.stars')
    data2 = loadtxt(image+'.match')
    diffs = []
    for item in data2:
        offsets = sqrt((item[0]-data1[:,0])**2.+(item[1]-data1[:,1])**2.)
        idx = where(offsets==min(offsets))
        if item[7] <30.:
            diffs.append(item[7]-data1[idx][:,3])
    diffs = array (diffs)
    print ('ave diffs = %2.1f\t mode = %2.1f\t stdev = %2.1f' %(average(diffs), binmode(diffs,bins=10)[0], std(diffs)))
    
    ax = subplot(111)
    ax.hist(diffs)#,bins=20)
    ax.axvline(average(diffs), color='red',
               label = r'$\overline{\delta}= %2.1f\,, \sigma = %2.1f$' %(average(diffs), std(diffs)))
    ax.axvline(binmode(diffs,bins=10)[0], color='red', linestyle='--',
               label = r'$Mo= %2.1f$' %(binmode(diffs,bins=10)[0]))
    ax.set_title('Histogram of zeropoints')
    ax.set_xlabel('Instrumental magnitude offset from USNO B1')
    ax.legend(loc=2,frameon=False)
    savefig(image+'_zmag.png')
    
    
                        
          
