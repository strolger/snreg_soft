#!/usr/bin/env python
# 2018
# L. Strolger
# Adapted from older perl version of imregister for NGSS/GOODS

#
"""
imregister.py [--options] template.fits image.fits
"""
import os,sys,pdb,scipy,glob
from pylab import *
import sextract


verbose = 1
clobber = False
thresh = 10
scale = 0.0
flagimg = None
goodmatch = 5
nstars = 50
tol = 0.005
maxstars = 900
matchsets = 3
order = 3
useiraf = False

software = os.path.dirname('/'.join(os.path.realpath(__file__).split('/')[:-1]))
sys.path.append(software)



def runmatch(starlist1, starlist2, matchlist, goodmatch=goodmatch, nstars=nstars, scale=scale, tol=tol,
           maxstars=maxstars, matchsets=matchsets, order=order, amiss=0.0, cbmiss=0.0,verbose=verbose):
    '''
    script to make the match from a pair of star lists
    '''
    import subprocess
    from subprocess import check_output
    
    try:
        sourcebindir = os.environ['SNPATH']
    except KeyError :
        sexbindir = ''  

    matchcmd = [ os.path.join( sourcebindir, "starmatch")]
    matchcmd.append('%s' %starlist2)
    matchcmd.append('%s' %starlist1)
    matchcmd.append('%f' %scale)
    matchcmd.append('%f 0' %tol)
    matchcmd.append('%d' %matchsets)
    matchcmd.append('%d' %nstars)
    if amiss and cbmiss: matchcmd.append('%f %f 0.05' %(amiss,cbmiss))
    if verbose: print(' '.join(matchcmd))
    solution = check_output(' '.join(matchcmd),shell=True)
    solution = reshape(array(list(map(float,solution.rstrip().split()))),(2,3))
    if not solution.any():
        print('starmatch fails.')
        return()

    if verbose:
        il = solution[0,:]
        print("\nx = %4.3f*x0 + %4.3f*y0 + %4.3f" %(il[0], il[1], il[2]))
        il = solution[1,:]
        print("y = %4.3f*x0 + %4.3f*y0 + %4.3f\n" %(il[0], il[1], il[2]))


    matchcmd = [ os.path.join( sourcebindir, "finalmatch")]
    matchcmd.append('%s' %starlist2)
    matchcmd.append('%s' %starlist1)
    matchcmd.append('%s' %matchlist)
    for item in solution.flatten():
        matchcmd.append('%f' %item)
    matchcmd.append('%d' %goodmatch)
    matchcmd.append('%d' %maxstars)
    matchcmd.append('0 0')
    if verbose: print(' '.join(matchcmd))
    subprocess.call(' '.join(matchcmd), shell=True)

    ## if amiss and cbmiss:
    ##     print('prelim matching done')
    ##     return()
    

    ## linear pass
    if verbose: print('Linear Pass')
    matchcmd = [ os.path.join( sourcebindir, "jtxform")]
    matchcmd.append('%s' %matchlist)
    matchcmd.append('tmp.map 3 fitinfo maxres=2.8')
    if verbose: print(' '.join(matchcmd))
    subprocess.call(' '.join(matchcmd), shell=True)
        
    matchcmd = [ os.path.join( sourcebindir, "rematch")]
    matchcmd.append('%s' %starlist1)
    matchcmd.append('%s' %starlist2)
    matchcmd.append('%s' %matchlist)
    matchcmd.append('tmp.map')
    matchcmd.append('%d' %goodmatch)
    matchcmd.append('%d' %maxstars)
    matchcmd.append('0 0')
    if verbose: print(' '.join(matchcmd))
    subprocess.call(' '.join(matchcmd), shell=True)

    ## quadratic pass
    wc = sum(1 for line in open(matchlist))
    if wc > 8 and order > 3:
        if verbose: print('Quadratic Pass')
        matchcmd = [ os.path.join( sourcebindir, "jtxform")]
        matchcmd.append('%s' %matchlist)
        matchcmd.append('tmp.map 6 fitinfo maxres=1.8')
        if verbose: print(' '.join(matchcmd))
        subprocess.call(' '.join(matchcmd), shell=True)
        
        matchcmd = [ os.path.join( sourcebindir, "rematch")]
        matchcmd.append('%s' %starlist1)
        matchcmd.append('%s' %starlist2)
        matchcmd.append('%s' %matchlist)
        matchcmd.append('tmp.map')
        matchcmd.append('%d' %goodmatch)
        matchcmd.append('%d' %maxstars)
        matchcmd.append('0 0')
        if verbose: print(' '.join(matchcmd))
        subprocess.call(' '.join(matchcmd), shell=True)
    ## elif os.path.isfile('jtxform.out'):
    ##     os.rename('jtxform.out', matchlist)

    ## cubic pass
    wc = sum(1 for line in open(matchlist))
    if wc > 14 and order > 3:
        if verbose : print ('Cubic pass')
        matchcmd = [ os.path.join( sourcebindir, "jtxform")]
        matchcmd.append('%s' %matchlist)
        matchcmd.append('tmp.map 10 fitinfo maxres=1.5')
        if verbose: print(' '.join(matchcmd))
        subprocess.call(' '.join(matchcmd), shell=True)
        
        matchcmd = [ os.path.join( sourcebindir, "rematch")]
        matchcmd.append('%s' %starlist1)
        matchcmd.append('%s' %starlist2)
        matchcmd.append('%s' %matchlist)
        matchcmd.append('tmp.map')
        matchcmd.append('%d' %goodmatch)
        matchcmd.append('%d' %maxstars)
        matchcmd.append('0 0')
        if verbose: print(' '.join(matchcmd))
        subprocess.call(' '.join(matchcmd), shell=True)
    ## elif os.path.isfile('jtxform.out'):
    ##     os.rename('jtxform.out', matchlist)

    if verbose: print('Final Matching of star lists done')
    return()

def runtransform(template, image, outimg, useiraf=False, matchlist=None, order=order, clobber=clobber):
    import subprocess
    from subprocess import check_output
    from shutil import copyfile
    if useiraf:
        import pyraf
        from pyraf import iraf
        from iraf import geomap, geotran, wcscopy
    
    try:
        sourcebindir = os.environ['SNPATH']
    except KeyError :
        sexbindir = ''

    if os.path.isfile(outimage) and not clobber:
        if verbose: print('%s exists... not overwritting' %outimg)
        return()
    
    transcmd = [ os.path.join( sourcebindir, "gethead")]
    transcmd.append('%s' %template)
    transcmd.append('naxis1 naxis2')
    if verbose: print(' '.join(transcmd))
    nx, ny = map(int,check_output(' '.join(transcmd),shell=True).split())

    if not useiraf:
        transcmd = [ os.path.join( sourcebindir, "jtxform")]
        transcmd.append('%s' %matchlist)
        transcmd.append('tmp.map')
        transcmd.append('%d' %order)
        transcmd.append('fitinfo maxres=5')
        if verbose: print(' '.join(transcmd))
        subprocess.call(' '.join(transcmd), shell=True)


        transcmd = [ os.path.join( sourcebindir, "jtmap")]
        transcmd.append('-src %s' %image)
        transcmd.append('-im %s' %outimg)
        transcmd.append('-nx %d -ny %d' %(nx,ny))
        transcmd.append('-map %s' %(image.replace('.fits','.map')))
        transcmd.append('-sub 2 -jacobian -bitpix -32')
        if verbose: print(' '.join(transcmd))
        subprocess.call(' '.join(transcmd), shell=True)
        if os.path.isfile(outimg):
            print('Successfully wrote %s' %outimg)
        else:
            print('Sorry... %s not written' %outimg)
    else:
        f = open(matchlist)
        lines = f.readlines()#[1:]
        f.close()
        if os.path.isfile('iraf.match'): os.remove('iraf.match')
        w = open('iraf.match','w')
        for line in lines:
            c = line.split()
            w.write('%s %s %s %s\n' %(c[2], c[3], c[0], c[1]))
        w.close()
        
        if order == 3: iorder = 2
        if order == 6: iorder = 3
        if order == 10: iorder = 4
        
        wc = sum(1 for line in open(matchlist))
        if (wc < 8):
            iorder =2
        elif (wc < 25 & order ==10):
            iorder = 3
        elif (wc < 50 & order ==10):
            iorder = 4
        elif (wc > 60 & order ==10):
            iorder = 5
        print ("\n\n %d matches using IRAF order: %d\n" %(wc,iorder))
        if os.path.isfile('%s.geomap'%image):os.remove('%s.geomap'%image)

        geomap('iraf.match','%s.geomap' %image, 1, nx, 1, ny,
               fitgeom='general', function='chebyshev',
               xxorder=iorder, xyorder=iorder, xxterms='none',
               yxorder=iorder, yyorder=iorder, yxterms='none',
               reject=2., calctype='real', verbose=True, interactive=False)
        
        geotran(image, outimg, '%s.geomap' %image, 'iraf.match', geometry='geometric',
                interpolant='linear', boundary='constant',constant=0.,
                fluxconserve=True, xmin=1., ymin=1.,xmax=nx,ymax=ny,verbose=True)

        wcscopy(outimg, template)
    return()

  

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
        
    # first arg is the template file
    template = arg[0]
    if template.endswith('.fits'): template=template.strip('.fits')
    if not os.path.isfile(template+'.fits'):
        print('there is an error \n filename %s does not exist\n' %template+'.fits')
        print(__doc__)
        pdb.set_trace()
        sys.exit(1)
    if verbose > 1: print('template:\t%s.fits' %template)
    if not os.path.isfile(template+'.sexcat') or not os.path.isfile(template+'.fits.stars') or clobber:
        sextract.runsex(template+'.fits',verbose=verbose, flagfile=flagimg, clobber=clobber, paramstring="-DETECT_THRESH %.2f" %thresh)
        f = open(template+'.sexcat')
        lines = f.readlines()
        f.close()
        w = open(template+'.fits.stars','w')
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
        

    img = arg[1]

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
        if image.endswith('.fits'): image=image.strip('.fits')
        if not os.path.isfile(image+'.fits'):
            print('there is an error \n filename %s does not exist\n' %image+'.fits')
            print(__doc__)
            sys.exit(1)
        if verbose > 1 : print ('image:\t%s.fits' %image)
        outimage = image+'.shift.fits'

        if not os.path.isfile(image+'.sexcat') or not os.path.isfile(image+'.fits.stars') or clobber:
            sextract.runsex(image+'.fits',verbose=verbose, flagfile=flagimg, clobber=clobber, paramstring="-DETECT_THRESH %.2f" %thresh)
            if verbose > 1: print ('writing %s.fits.stars' %image)
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

        runmatch(template+'.fits.stars',image+'.fits.stars',image+'.match',order=order, scale=scale)
        if os.path.isfile('tmp.map'): os.rename('tmp.map','%s.map'%image)
        runtransform(template+'.fits',image+'.fits',image+'.shift.fits', matchlist=image+'.match',useiraf=useiraf,order=order,clobber=clobber)

    if verbose: print('Done.')
