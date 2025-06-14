#! /usr/bin/env python
# 2018
# L. Strolger
# Adapted from code written for FFSN
# by S. Rodney


# 2010
# S.Rodney
#
# Defines a python class for Sextractor catalogs
# and functions to run sextractor on an image 
# with an input bad pixel mask, and to clean the 
# resulting catalog by rejecting obvious bad sources
# 2012.07.17 master zeropoint dicts kept in tools.hstsnphot
"""
sextract.py [--options] image.fits [sex-options]

Command-line access to run sextractor with HSTSNpipeline defaults
on a single image.  If the image has suffix _sci.fits or _masked.fits
and the corresponding _weight.fits and _badpix.fits files are present,
then they are automatically used.

--verbose X   :  run with verbose level X, 0-10
--checkaper   :  make a checkimage showing apertures 
--clobber     :  overwrite existing images

Any options given after the image name are passed directly to sextractor.
e.g.:   sextract.py --checkaper some_image_sci.fits -DETECT_MINAREA 3 -DETECT_THRESH 2
"""

class Catalog(object):
    """ A catalog of astronomical sources """
    def __init__( self, catfile=None, imfile=None ):
        """ Initiate a catalog of sources with
        (x,y) pixel coordinates and/or (ra,dec)
        sky coordinates, as well as measured 
        magnitudes and/or fluxes.
        If provided, the 'imfile' is an image file 
        from which to extract WCS information for
        connecting (x,y) to (ra,dec) and the 
        'catfile' is the name of an ASCII catalog 
        file for reading/writing the list of sources.
        """
        self.coldata = {} # dictionary of data, keyed by column name
        self.colindex = {}
        if catfile : 
            import os
            catfile = os.path.abspath( catfile )
            self.readcatfile( catfile ) 
            self.filename = catfile 

    @property
    def colnames( self ):
        """ return a list of column names """
        return( [ colname for colname in self.coldata.keys() ] )

    @property
    def nrows( self ):
        if len( self.coldata.values() )>0  : 
            return( len( self.coldata.values()[0] ) )
        else : 
            return( 0 )

    @property
    def ncols( self ):
        return( len( self.coldata.keys() ) )

    def getrow( self, rownum ) :
        """ extract the given row from the catalog data"""
        row = []
        for colindex in range(self.ncols) : 
            row.append( self.data[colindex][rownum] )
        return( row )

    def getrowdict( self, rownum ) :
        """ extract the given row from the catalog data
        as a dictionary keyed by column names"""
        rowdict = {}
        for colname in self.colnames : 
            rowdict[colname] = self.data[colname][rownum]
        return( rowdict )

    @property
    def ra( self ):
        for k in ['X_WORLD','RA','ra','R.A.','ALPHA_J2000']: 
            if k in self.colnames : 
                return( self.coldata[k] )
        return( [] )

    @property
    def dec( self ):
        for k in ['Y_WORLD','DEC','dec','Dec','Decl','DELTA_J2000']: 
            if k in self.colnames : 
                return( self.coldata[k] )
        return( [] )

    @property
    def x( self ):
        if 'X_IMAGE' in self.colnames : 
            return( self.coldata['X_IMAGE'] )
        else : 
            return( [] )

    @property
    def y( self ):
        if 'Y_IMAGE' in self.colnames : 
            return( self.coldata['Y_IMAGE'] )
        else : 
            return( [] )



    # # TODO : UPDATE THIS
    # def appendrow( self, row ):
    #     """ append a single row of data,
    #     The row of data may be given as a 
    #     python list or a dict, keyed by 
    #     column names.
    #     """
    #     self.data.insert( 1, -1 )
    #     self.updaterow( self.nrows-1, row )

    # # TODO : UPDATE THIS
    # def updaterow( self, index, row ):
    #     """ update a single row of data at 
    #     the given row index. 
    #     The row of data may be given as a 
    #     python list or a dict, keyed by 
    #     column names.
    #     """
    #     if type(row) == list : 
    #         for colindex in range(len(row)) : 
    #             self.data[colindex][index] = row[colindex]
    #     elif type(row) == dict : 
    #         for column in row.keys()  : 
    #             self.data[column][index] = row[column]

    def delrow( self, index ):
        self.data.delete( index )


    def readcatfile( self, catfile ) :
        """ read a sextractor catalog """
        # TODO : handle different table types
        from numpy import loadtxt, array, genfromtxt
        # read the sex-cat header
        fin = open( catfile, 'r' )
        header = fin.readlines()
        fin.close()
        # TODO : handle catalogs with a single row of data !!
        #catdat = loadtxt( catfile )
        catdat = genfromtxt( catfile )
        self.coldata = {}
        self.colindex = {}
        # guard against single-line catalogs: 
        if len(catdat.shape)==1 : catdat = array([ catdat ])
        idfirst=True
        for line in header :
            if not line.startswith('#'): break
            try : 
                icol = int( line.split()[1] ) - 1
            except ValueError : 
                try: 
                    # is this a reverse-formatted header? 
                    idfirst = False
                    icol = int( line.split()[-1].strip('()') ) - 1
                except e: 
                    continue
            if idfirst : colname =  line.split()[2]
            elif '(' in line : 
                colname = line.lstrip('#').split('(')[0].strip().replace(' ','_')
            self.colindex[icol+1] = colname

            if colname in ['NUMBER', 'FLAGS'] : 
                coldata = catdat[:,icol].astype( int )
            else :
                coldata = catdat[:,icol].astype( float )
            self.coldata[colname] = coldata 
        #self.data = asciidata.open( catfile )


    def writecatfile( self, catfile, colInfo=None, headComment=None, selectrows=[] ) :
        """ write out catalog as a sextractor .cat file.
        selectrows : a sequence of row numbers to select for output. 
            If empty, all rows are printed.
        """
        fout = open( catfile, 'w' )
        iclist = self.colindex.keys()
        iclist.sort()
        for icol in iclist : 
            print >> fout,"# %3i %15s "%(icol, self.colindex[icol] )
        for irow in range(self.nrows) :
            if len(selectrows) and irow not in selectrows : continue
            outstring = '  '
            for icol in iclist :
                colname = self.colindex[icol]
                datum = self.coldata[colname][irow] 
                if datum % 1 :
                    if abs(datum)>1e-4 : 
                        outstring += ' %14.7f'%datum
                    else : 
                        outstring += ' %14.7e'%datum
                else : 
                    outstring += ' %-14i'%datum
            print >> fout, outstring
        fout.close()
        return( catfile )

    def wds9reg( self, regfile=None, color='red', circlesize=10, width=1,
                 xcol='X_IMAGE', ycol='Y_IMAGE', textcol=None, criteria=None):
        """ write sources to a ds9 region file.
        circlesize may be a scalar or an array with one size 
        for each object in the catalog. Criteria can now be used to only write certian regions 
        to the file from the catelog. criteria is a tuple with values (colname,comparison_operator,value).
        Eg.     criteria = ('EXT_IMAGE','==',1)
        """
        import numpy as np
        import os
        #import exceptions

        if( xcol not in self.colnames or ycol not in self.colnames ):
            print("ERROR: cannot make a ds9 region file without %s,%s"%(xcol,ycol))
            return( -1 )

        xcolname = xcol.lower() 
        convertRADEC=False
        if ('alpha' in xcolname) or ('ra' in xcolname) or ('r.a.' in xcolname) or ('world' in xcolname) : 
            if isinstance( self.coldata[xcol][0], float) :  
                try : 
                    from astLib import astCoords
                    print( "Converting RA,DEC from decimal degrees to hh:mm:ss,dd:mm:ss")
                    convertRADEC=True
                except : 
                    print( "Can't load astLib for RA,DEC conversion!" )
                    return(-1)

        if not regfile : 
            regfile = os.path.splitext( self.filename )[0] +'.reg'
 
        fout = open(regfile,'w')
        print>>fout,"# Region file format: DS9 version 5.4"
        print>>fout,"""global color=%s font="times 16 bold" select=1 edit=1 move=0 delete=1 include=1 fixed=0 width=%i"""%(color,width)
        for irow in range(self.nrows):
            if np.iterable(circlesize): cs = circlesize[irow]
            else : cs = circlesize

            textstr = ''
            if textcol : 
                # textcol may be a comma-sep'd list of col names.
                textvals = [] 
                for tc in textcol.split(',') : 
                    tval = self.coldata[tc][irow]
                    if type(tval)in [int,np.int,np.int16,np.int32,np.int64,np.int8,np.int_] : 
                        textvals.append('%i'%tval)
                    elif type(tval) in [float,np.float,np.float32,np.float64,np.float128,np.float_]  : 
                        textvals.append('%.2f'%tval)
                    else : textvals.append(str(tval))
                textstr = '# text={%s}'%':'.join( textvals )
           
            #Only write the line to the file if it passes the criteria passed to this function. 
            #This allows us to separate extensions for makeing flt region files in the fakes
            if criteria is None: write_line = True
            else:write_line = eval("self.coldata['%s'][irow] %s %s"%criteria)
            
            if write_line: 
                cstr = '%.1f'%cs
                x, y = self.coldata[xcol][irow], self.coldata[ycol][irow]
                if isinstance( x, float) :
                    xstr, ystr = '%.8f'%x, '%.8f'%y
                    if convertRADEC : 
                        xstr = astCoords.decimal2hms( x, delimiter=':' )
                        ystr = astCoords.decimal2dms( y, delimiter=':' )
                        if cs > 1 : cs = min( 0.75, cs*0.06 )
                        cstr = '%.3f"'%( cs )
                print>>fout,'circle(%s,%s,%s)%s'%( xstr, ystr, cstr, textstr )

        fout.close() 
        return( regfile )

    def wxycat( self, catfile ):
        """ write sources to a simple x,y catalog file.
        (useful for use with the Vista marktv command)

        NEEDS UPDATING !!!
        """
        from numpy import iterable

        XCOL,YCOL = 'X_IMAGE','Y_IMAGE'
        if XCOL in self.colnames : 
            xcol = XCOL
        else : 
            for colname in self.colnames: 
                if colname.find(XCOL)>=0 :
                    xcol = colname
                    break
        if YCOL in self.colnames : 
            ycol = YCOL
        else : 
            for colname in self.colnames: 
                if colname.find(YCOL)>=0 :
                    ycol = colname
                    break
        fout = open(catfile,'w')
        for i in range(self.nrows):
            print>>fout,"%12.1f %12.1f"%(
                self.coldata[xcol][i],
                self.coldata[ycol][i] )


def wht2var( infile, outfile=None, platescale='default',
             verbose=False ):
    """
    convert a _wht.fits weight map from multi-drizzle 
    into a variance map suitable for use as a 
    sextractor weight map.
    following Casertano et al 2006.
    """
    #import astropy.io.fits as pyfits
    from astropy.io import fits as pyfits
    from numpy import nan_to_num

    # size of the input pixels in arcsec
    # default is to read it from the config file
    if platescale=='default': pscalein = PLATESCALE

    # get size of the output pixel in arcsec
    # and pixfrac in units of input pixels from header
    try : 
        pscaleout = 3600*( 
            abs(pyfits.getval(infile,'CD1_1')) + 
            abs(pyfits.getval(infile,'CD2_2')))/2.  
    except KeyError :
        pscaleout = 3600*( 
            abs(pyfits.getval(infile,'CD1_1',ext=1)) + 
            abs(pyfits.getval(infile,'CD2_2',ext=1)))/2.  
    try : 
        pixfrac = pyfits.getval(infile,'PIXFRAC')
    except KeyError: 
        try: pixfrac = pyfits.getval(infile,'D001PIXF')
        except KeyError: 
            # TODO : FIX THIS WHEN MAKING WHT IMAGES!!!!
            pixfrac=1

    l = pscaleout / pscalein
    p = pixfrac  
    if l<=p : Fa = (l/p)**2 * (1-l/(3*p))**2
    else : Fa = (1 - p/(3*l))**2

    # open file, check for multiple extensions
    inim = pyfits.open( infile )
    extnames = [ext.name for ext in inim]
    whtext=0
    if len(inim)>1:
        if 'WHT' in extnames : whtext = extnames.index('WHT')
        else : whtext=2
    
    # make the variance map
    vardata = nan_to_num(  Fa / inim[whtext].data )
    varhdu = pyfits.ImageHDU( 
        data=vardata, header=inim[whtext].header, name='VAR' )
    varhdu.header['EXTNAME']='VAR'

    # write it to a new extension or a new file
    if outfile==None or outfile==infile : 
        if 'VAR'  in extnames : inim['VAR'] = varhdu
        else : inim.insert( whtext+1, varhdu )
        if verbose : print("adding VAR ext to %s"%infile)
        inim.writeto( infile, clobber=True )
    else : 
        if verbose : print("writing var map to %s"%outfile)
        pyfits.writeto( outfile, vardata, clobber=True )

    
def runsex(imfile, catfile=None, weightfile=None, flagfile=None,
           paramstring='', detimfile=None, returncat=False, verbose=False, clobber=False):
    '''
    Run sextractor on a single image.
    if returncat==True : Returns a Catalog object.
    else : returns a string with the catalog filename
    specify detimfile to run in double-image mode. 
    ''' 
    import shutil
    import os
    import subprocess
    #import exceptions
    #import astropy.io.fits as pyfits
    from astropy.io import fits as pyfits
    ## from tools import hstsnphot

    # locate the sextractor binary and the config file dir
    try : sexbindir = os.environ['SNPATH']
    except KeyError : sexbindir = ''  
    try : 
        # if we got here via startpipe.py :
        sexconfigdir = os.path.join(os.environ['SNSOURCE'],'Sex')
    except KeyError :
        # if we're working from a command-line call using the pipeline installation: 
        sexconfigdir = os.path.join( os.path.dirname(os.path.abspath(os.path.realpath(__file__))),'sexconfig')
    if not os.path.isdir( sexconfigdir ) : 
        
        # otherwise, fingers crossed...
        sexconfigdir = ''

    if verbose>2 : 
        print( 'using %s with config files from %s'%(   
                os.path.join(sexbindir,'sex'), sexconfigdir ) )

    imfile = os.path.abspath( imfile )
    if catfile==None : 
        imroot = os.path.splitext( imfile )[0]
        catfile = imroot + ".sexcat"

    if os.path.isfile( catfile ) and not clobber : 
        print("%s exists. Not clobbering."%catfile )
        if returncat : 
            sexcat = Catalog( imfile=imfile, catfile=catfile )
            return( sexcat )
        else : 
            return( catfile )

    # Using some generic cofiguration files, good for ground work.
    configfile = os.path.join( sexconfigdir,"default.sex")
    paramfile = os.path.join( sexconfigdir,"daofind.param")
    convfile = os.path.join( sexconfigdir,"default.conv")
    nnwfile = os.path.join( sexconfigdir,"default.nnw")


    # run sextractor 
    sexcommand = [ os.path.join( sexbindir, "sex" ) ]
    if detimfile : 
        sexcommand.append( "%s,%s"%(detimfile,imfile) )
    else : 
        sexcommand.append("%s"%imfile)
    sexcommand.append("-c %s"%configfile) 
    sexcommand.append("-FILTER_NAME %s"%convfile)
    sexcommand.append("-STARNNW_NAME %s"%nnwfile)
    sexcommand.append("-CATALOG_NAME %s"%catfile)
    sexcommand.append("-CATALOG_TYPE ASCII_HEAD")

    if weightfile : 
        sexcommand.append("-WEIGHT_IMAGE %s"%weightfile)
        sexcommand.append("-WEIGHT_TYPE MAP_WEIGHT")
    else : 
        sexcommand.append("-WEIGHT_TYPE NONE")

    if flagfile : 
        paramfile = os.path.join( sexconfigdir,"daofind_flags.param")
        sexcommand.append("-FLAG_IMAGE %s"%flagfile)
    else : 
        sexcommand.append("-FLAG_IMAGE NONE")
    sexcommand.append("-PARAMETERS_NAME %s"%paramfile)

    if verbose>3: 
        sexcommand.append("-VERBOSE_TYPE FULL")
    elif verbose>1: 
        sexcommand.append("-VERBOSE_TYPE NORMAL")
    else :
        sexcommand.append("-VERBOSE_TYPE QUIET")
    sexcommand.append(paramstring)

    if verbose : print( ' '.join(sexcommand) )
    subprocess.call( ' '.join(sexcommand), shell=True )

    if not os.path.isfile( catfile ) :
        raise RuntimeError( 
            "sextractor failed to produce expected catalog file %s!!"%catfile)

    # # read in the catalog
    if returncat : 
        sexcat = Catalog( imfile=imfile, catfile=catfile )
        return( sexcat )
    else : 
        return( catfile )


def cleancat( catfile, newcatfile=None, 
              fwhm_min=1.8, flag_max=1, Bmin=1,
              isogood_min=0.3,
              verbose=False, clobber=False ): 
    """ 
    read in the given catfile, 
    clean it of obvious cosmic rays and bad sources
    write it out to newcatfile and return the 
    name of the newcatfile. 
    """
    # TODO : this is painfully slow
    # should be far faster with numpy 
    import os
    cat = Catalog( catfile=catfile )
    newcat = Catalog( )

    if not newcatfile :
        newcatfile = os.path.splitext( catfile )[0] +"_clean.cat"

    if os.path.isfile( newcatfile ) and not clobber: 
        print("%s already exists. Not clobbering."%newcatfile)
        return( newcatfile )

    newcat.filename = newcatfile 
    for colname in cat.colnames : 
        newcat.coldata[colname] = []
    for icol in cat.colindex : 
        newcat.colindex[icol] = cat.colindex[icol]
    for irow in range(cat.nrows) : 
        if fwhm_min :
            if cat.coldata['FWHM_IMAGE'][irow] < fwhm_min : continue
        if flag_max : 
            if cat.coldata['FLAGS'][irow] > flag_max : continue
        if Bmin : 
            if cat.coldata['B_IMAGE'][irow] < Bmin : continue

        if 'NIMAFLAGS_ISO' in cat.colnames and 'ISOAREA_IMAGE' in cat.colnames : 
            isogood=float(cat.coldata['NIMAFLAGS_ISO'][irow]) / cat.coldata['ISOAREA_IMAGE'][irow]
            if isogood < isogood_min :  continue

        for colname in newcat.colnames : 
            newcat.coldata[colname].append( cat.coldata[colname][irow] )
    newcat.writecatfile( newcatfile )
    return( newcatfile )
        


def astromatch(cat1, cat2, radius, coords='sky'):
    '''
    Simple astrometric catalog matching between two catalog objects.  
    Returns an Nx4 array with the N sources in cat1 that have matches
    within 'radius' arcesconds in cat2. 
    The coords keyword specifies the coordinate system used for 
    the return array : 
       coords='sky'  :   [ra1,dec1,ra2,dec2]
       coords='xy'   :   [x1, y1, x2, y2] 
    '''
    import numpy as np
    if type( cat1 ) == str : 
        cat1 = Catalog( cat1 ) 
    if type( cat2 ) == str : 
        cat2 = Catalog( cat2 ) 
    ra1, dec1, x1, y1  = cat1.ra,cat1.dec,cat1.x,cat1.y
    ra2, dec2, x2, y2  = cat2.ra,cat2.dec,cat2.x,cat2.y

    outlist = []
    for i in range(0,cat1.nrows):
        ra_scale = np.cos( dec1[i]*np.pi/180.)
        delta_ra = 3600.*(ra1[i] - ra2)*ra_scale
        delta_dec = 3600.*(dec1[i] - dec2) 
        dist = np.sqrt(delta_ra**2+delta_dec**2)
        dmin = np.amin(dist)
        jmin = np.where(dist==dmin)
        if dmin <= radius: 
            if coords=='sky' : 
                outlist.append( [ra1[i],dec1[i],ra2[jmin],dec2[jmin]] )
            elif coords=='xy' :
                outlist.append( [x1[i],y1[i],x2[jmin],y2[jmin]] )
    return np.array( outlist )



if __name__ == '__main__':
    import sys
    import getopt
    import os

    verbose=1
    checkaper=False
    clobber=False

    # read in arguments and options
    try:
        opt,arg = getopt.getopt( 
            sys.argv[1:],"v,h",
            longopts=["verbose=","quiet","help","checkaper","clobber" ] )
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
        elif o == '--checkaper' :
            checkaper=True
        elif o == '--clobber' :
            clobber=True

    # first arg is the image filename    
    imfile = arg[0]

    ## # Check for _weight.fits and _badpix.fits images
    ## sfx = '.fits'
    ## if imfile.endswith('_sci.fits') : sfx = '_sci.fits'
    ## elif imfile.endswith('_masked.fits') : sfx = '_masked.fits'
    ## whtfile = imfile.replace( sfx, '_weight.fits') 
    ## if not os.path.isfile( whtfile ) : 
    ##     whtfile = imfile.replace( sfx, '_wht.fits') 
    ##     if not os.path.isfile( whtfile ) : 
    ##         whtfile = None
    ## bpxfile = imfile.replace( sfx, '_badpix.fits') 
    ## if not os.path.isfile( bpxfile ) : 
    ##     bpxfile = imfile.replace( sfx, '_bpx.fits') 
    ##     if not os.path.isfile( bpxfile ) : 
    ##         bpxfile = None

    whtfile = None
    bpxfile = None
    
    # any additional args get passed straight to sextractor
    paramstring = ''
    if len(arg) > 1 : 
        paramstring = ' '.join(arg[1:])

    if checkaper : 
        paramstring+=' -CHECKIMAGE_TYPE APERTURES -CHECKIMAGE_NAME %s '%imfile.replace('.fits','_aper.fits')

    runsex(imfile, catfile=None, weightfile=whtfile, flagfile=bpxfile,
           paramstring=paramstring, detimfile=None, verbose=verbose, clobber=clobber)
    

