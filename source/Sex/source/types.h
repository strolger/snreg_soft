 /*
 				types.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP, Leiden observatory & ESO)
*
*	Contents:	global type definitions.
*
*	Last modify:	11/08/98
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include <stdio.h>

/*-------------------------------- flags ------------------------------------*/

#define		OBJ_CROWDED	0x0001
#define		OBJ_MERGED	0x0002
#define		OBJ_SATUR	0x0004
#define		OBJ_TRUNC	0x0008
#define		OBJ_APERT_PB	0x0010
#define		OBJ_ISO_PB	0x0020
#define		OBJ_DOVERFLOW	0x0040
#define		OBJ_OVERFLOW	0x0080

/*---------------------------- preanalyse flags -----------------------------*/

#define		ANALYSE_FAST		0
#define		ANALYSE_FULL		1
#define		ANALYSE_ROBUST		2

/*--------------------------------- typedefs --------------------------------*/
typedef	unsigned char	BYTE;			/* a byte */
typedef	unsigned short	USHORT;			/* 0 to 65535 integers */
typedef	float		PIXTYPE;		/* pixel type */
typedef unsigned int	FLAGTYPE;		/* flag type */
typedef	char		pliststruct;		/* Dummy type for plist */

#ifdef INT64
typedef	int		LONG;
typedef	unsigned int	ULONG;
#else
typedef	long		LONG;
typedef	unsigned long	ULONG;
#endif

typedef  enum {CHECK_NONE, CHECK_IDENTICAL, CHECK_BACKGROUND,
        CHECK_BACKRMS, CHECK_MINIBACKGROUND, CHECK_MINIBACKRMS,
        CHECK_SUBTRACTED, CHECK_FILTERED, CHECK_OBJECTS, CHECK_APERTURES,
	CHECK_SEGMENTATION, CHECK_ASSOC, CHECK_SUBOBJECTS,
	CHECK_SUBPROTOS, CHECK_PROTOS, CHECK_MAPSOM}	checkenum;
	/* CHECK_IMAGE type */

typedef  enum {WEIGHT_NONE, WEIGHT_FROMBACK, WEIGHT_FROMRMSMAP,
		WEIGHT_FROMVARMAP, WEIGHT_FROMWEIGHTMAP, WEIGHT_FROMINTERP}
		weightenum;				/* WEIGHT_IMAGE type */

/*--------------------------------- objects ---------------------------------*/
/* I: "PIXEL" parameters */

typedef struct
  {
/* ---- basic parameters */
  int		number;				/* ID */
  int		fdnpix;				/* nb of extracted pix */
  int		dnpix;				/* nb of pix above thresh  */
  int		npix;				/* "" in measured frame */
  float		fdflux;				/* integrated ext. flux */
  float		dflux;				/* integrated det. flux */
  float		flux;				/* integrated mes. flux */
  float		fluxerr;			/* integrated variance */
  PIXTYPE	fdpeak;				/* peak intensity (ADU) */
  PIXTYPE	dpeak;				/* peak intensity (ADU) */
  PIXTYPE	peak;				/* peak intensity (ADU) */
/* ---- astrometric data */
  int		peakx,peaky;			/* pos of brightest pix */
  double       	mx, my;				/* barycenter */
  double	poserr_mx2, poserr_my2,
		poserr_mxy;			/* Error ellips moments */
/* ---- morphological data */			
  int		xmin,xmax,ymin,ymax,ycmin,ycmax;/* x,y limits */
  PIXTYPE	*blank, *dblank; 	       	/* BLANKing sub-images  */
  int		*submap;			/* Pixel-index sub-map */
  int		subx,suby, subw,subh;		/* sub-image pos. and size */
  short		flag;				/* extraction flags */
  FLAGTYPE	imaflag[MAXFLAG];		/* flags from FLAG-images */
  BYTE		singuflag;			/* flags for singularities */
  int		imanflag[MAXFLAG];     		/* number of MOST flags */
  double	mx2,my2,mxy;			/* variances and covariance */
  float		a, b, theta, abcor;		/* moments and angle */
  float		cxx,cyy,cxy;			/* ellipse parameters */
  int		firstpix;			/* ptr to first pixel */
  int		lastpix;			/* ptr to last pixel */
  float		bkg, dbkg, sigbkg;		/* Background stats (ADU) */
  float		thresh;				/* measur. threshold (ADU) */
  float		dthresh;		       	/* detect. threshold (ADU) */
  float		mthresh;		       	/* max. threshold (ADU) */
  int		iso[NISO];			/* isophotal areas */
  float		fwhm;				/* IMAGE FWHM */
  }	objstruct;

/* II: "BLIND" parameters */
typedef struct
  {
/* ---- photometric data */
  float		flux_iso;			/* ISO integrated flux */
  float		fluxerr_iso;			/* RMS error on ISO flux */
  float		mag_iso;			/* ISO mag */
  float		magerr_iso;			/* ISO mag uncertainty */
  float		flux_isocor;			/* ISOCOR integrated flux */
  float		fluxerr_isocor;			/* RMS error on ISOCOR flux */
  float		mag_isocor;			/* ISOCOR mag */
  float		magerr_isocor;			/* ISOCOR mag uncertainty */
  float		kronfactor;			/* kron parameter */
  float		flux_auto;			/* AUTO integrated flux */
  float		fluxerr_auto;			/* RMS error on AUTO flux */
  float		mag_auto;			/* AUTO mag */
  float		magerr_auto;			/* AUTO mag uncertainty */
  float		flux_best;			/* BEST integrated flux */
  float		fluxerr_best;			/* RMS error on BEST flux */
  float		mag_best;			/* BEST mag */
  float		magerr_best;			/* BEST mag uncertainty */
  float		*flux_aper;			/* APER flux vector */
  float		*fluxerr_aper;			/* APER flux error vector  */
  float		*mag_aper;			/* APER magnitude vector */
  float		*magerr_aper;			/* APER mag error vector */
/* ---- astrometric data */
  double	posx,posy;			/* "FITS" pos. in pixels */
  double	mamaposx,mamaposy;		/* "MAMA" pos. in pixels */
  float		sposx,sposy;			/* single precision pos. */
  float		poserr_a, poserr_b,
		poserr_theta;			/* Error ellips parameters */
  float		poserr_cxx, poserr_cyy,
		poserr_cxy;			/* pos. error ellipse */
  double	poserr_mx2w, poserr_my2w,
		poserr_mxyw;			/* WORLD error moments */
  float		poserr_aw, poserr_bw,
		poserr_thetaw;			/* WORLD error parameters */
  float		poserr_thetas;			/* native error pos. angle */
  float		poserr_theta2000;		/* J2000 error pos. angle */
  float		poserr_theta1950;		/* B1950 error pos. angle */
  float		poserr_cxxw, poserr_cyyw,
		poserr_cxyw;			/* WORLD error ellipse */
  double	mx2w,my2w,mxyw;			/* WORLD var. and covar. */
  double	peakxw, peakyw;			/* WORLD of brightest pix */
  double	mxw, myw;			/* WORLD barycenters */
  double	alphas, deltas;			/* native alpha, delta */
  float		thetas;				/* native position angle E/N*/
  double	peakalphas, peakdeltas;		/* native for brightest pix */
  double	peakalpha2000, peakdelta2000;	/* J2000 for brightest pix */
  double	peakalpha1950, peakdelta1950;	/* B1950 for brightest pix */
  double	alpha2000, delta2000;		/* J2000 alpha, delta */
  float		theta2000;			/* J2000 position angle E/N */
  double	alpha1950, delta1950;		/* B1950 alpha, delta */
  float		theta1950;			/* B1950 position angle E/N */
  float		aw, bw;				/* WORLD ellipse size */
  float		thetaw;				/* WORLD position angle */
  float		cxxw,cyyw,cxyw;			/* WORLD ellipse parameters */
  float		npixw, fdnpixw;			/* WORLD isophotal areas */
  float		threshmu;			/* det. surface brightnees */
  float		maxmu;				/* max. surface brightnees */
  float		elong;				/* elongation */
  float		ellip;				/* ellipticity */
  float		sprob;				/* Stellarity index */
  float		fwhmw;				/* WORLD FWHM */
  float		*assoc;				/* ASSOCiated data */
  int		assoc_number;			/* nb of ASSOCiated objects */
  float		*vignet;			/* Pixel data */
  float		*vigshift;			/* (Shifted) pixel data */
/* ---- SOM fitting */
  float		flux_somfit;			/* Fitted amplitude */
  float		fluxerr_somfit;			/* RMS error on SOM flux */
  float		mag_somfit;			/* Magnitude from SOM fit */
  float		magerr_somfit;			/* Mag. err. from SOM fit */
  float		stderr_somfit;			/* Fitting reduced error */
  float		*vector_somfit;			/* SOM fit vector position */
/* ---- Growth curves and stuff */
  float		*flux_growth;			/* Cumulated growth_curve */
  float		flux_growthstep;		/* Growth-curve step */
  float		*mag_growth;			/* Cumulated growth_curve */
  float		mag_growthstep;			/* Growth-curve step */
  float		*flux_radius;			/* f-light-radii */
/* ---- PSF-fitting */
  float		flux_psf;			/* Flux from PSF-fitting */
  float		fluxerr_psf;			/* RMS error on PSF flux */
  float		mag_psf;			/* Mag from PSF-fitting */
  float		magerr_psf;			/* RMS mag from PSF-fitting */
  float		x_psf,y_psf;			/* Coords from PSF-fitting */
  short		niter_psf;			/* # of PSF-fitting iterat. */
  float		chi2_psf;			/* Red. chi2 of PSF-fitting */
  double	xw_psf, yw_psf;			/* WORLD coords */
  double	alphas_psf, deltas_psf;		/* native alpha, delta */
  double	alpha2000_psf, delta2000_psf;	/* J2000 alpha, delta */
  double	alpha1950_psf, delta1950_psf;	/* B1950 alpha, delta */
  double	poserrmx2_psf, poserrmy2_psf,
		poserrmxy_psf;			/* Error ellips moments */
  float		poserra_psf, poserrb_psf,
		poserrtheta_psf;		/* Error ellips parameters */
  float		poserrcxx_psf, poserrcyy_psf,
		poserrcxy_psf;			/* pos. error ellipse */
  double	poserrmx2w_psf, poserrmy2w_psf,
		poserrmxyw_psf;			/* WORLD error moments */
  float		poserraw_psf, poserrbw_psf,
		poserrthetaw_psf;		/* WORLD error parameters */
  float		poserrthetas_psf;		/* native error pos. angle */
  float		poserrtheta2000_psf;		/* J2000 error pos. angle */
  float		poserrtheta1950_psf;		/* B1950 error pos. angle */
  float		poserrcxxw_psf, poserrcyyw_psf,
		poserrcxyw_psf;			/* WORLD error ellipse */
  }	obj2struct;

/*----------------------------- lists of objects ----------------------------*/
typedef struct
  {
  int		nobj;			/* number of objects in list */
  objstruct	*obj;			/* pointer to the object array */
  int		npix;			/* number of pixels in pixel-list */
  pliststruct	*plist;			/* pointer to the pixel-list */
  PIXTYPE	dthresh;		/* detection threshold */
  PIXTYPE	thresh;			/* analysis threshold */
  }	objliststruct;


/*----------------------------- image parameters ----------------------------*/
typedef struct
  {
  char		filename[MAXCHAR];	/* pointer to the image filename */
  char		*rfilename;		/* pointer to the reduced image name */
  char		ident[MAXCHAR];		/* field identifier (read from FITS)*/
  char		rident[MAXCHAR];	/* field identifier (relative) */
  FILE		*file;			/* pointer the image file structure */
  char		*fitshead;		/* pointer to the FITS header */
  int		fitsheadsize;		/* FITS header size */
/* ---- main image parameters */
  int		bitpix, bytepix;	/* nb of bits and bytes per pixel */
  int		bitsgn;			/* non-zero if signed integer data */
  int		width, height;		/* x,y size of the field */
  size_t	npix;			/* total number of pixels */
  double	bscale, bzero;		/* FITS scale and offset */
  double	ngamma;			/* normalized photo gamma */
  int		nlevels;		/* nb of quantification levels */
  float		pixmin, pixmax;		/* min and max values in frame */
  int		y;			/* y current position in field */
  int		ymin;			/* y limit (lowest accessible) */
  int		ymax;			/* y limit (highest accessible+1) */
  int		yblank;			/* y blanking limit (highest+1) */
  PIXTYPE	*strip;			/* pointer to the image buffer */
  PIXTYPE	*copystrip;	       	/* pointer to another image buffer */
  FLAGTYPE	*fstrip;		/* pointer to the FLAG buffer */
  int		stripheight;		/* height  of a strip (in lines) */
  int		stripmargin;		/* number of lines in margin */
  int		stripstep;		/* number of lines at each read */
  int		stripy;			/* y position in buffer */
  int		stripylim;		/* y limit in buffer */
  int		stripysclim;		/* y scroll limit in buffer */
/* ---- image (de-)compression */
  enum {COMPRESS_NONE, COMPRESS_BASEBYTE, COMPRESS_PREVPIX}
		compress_type;		/* image compression type */
  char		*compress_buf;		/* de-compression buffer */
  char		*compress_bufptr;	/* present pixel in buffer */
  int		compress_curval;	/* current pixel or checksum value */
  int		compress_checkval;	/* foreseen pixel or checksum value */
  int		compress_npix;		/* remaining pixels in buffer */
/* ---- basic astrometric parameters */
   double	pixscale;		/* pixel size in arcsec.pix-1 */
   double	epoch;			/* epoch of coordinates */
/* ---- background parameters */
  float		*back;			/* ptr to the background map in mem */
  float		*dback;			/* ptr to the background deriv. map */
  float		*sigma;			/* ptr to the sigma map */
  float		*dsigma;		/* Ptr to the sigma deriv. map */
  int		backw, backh;		/* x,y size of a bkgnd mesh */
  int		nbackp;			/* total nb of pixels per bkgnd mesh */
  int		nbackx, nbacky;		/* x,y number of bkgnd meshes */
  int		nback;			/* total number of bkgnd meshes */
  int		nbackfx, nbackfy;	/* x,y size of bkgnd filtering mask */
  float		backmean;		/* median bkgnd value in image */
  float		backsig;		/* median bkgnd rms in image */
  float		sigfac;			/* scaling RMS factor (for WEIGHTs) */
  PIXTYPE	*backline;		/* current interpolated bkgnd line */
  PIXTYPE	dthresh;		/* detection threshold */
  PIXTYPE	thresh;			/* analysis threshold */
/* ---- astrometric parameters */
  struct structastrom	*astrom;	/* astrometric data */
  struct structassoc	*assoc;		/* ptr to the assoc-list */
  int		flags;			/* flags defining the field type */
/* ---- image interpolation */
  int		interp_flag;		/* interpolation for this field? */
  PIXTYPE	*interp_backup;		/* backup line for interpolation */
  PIXTYPE	weight_thresh;		/* interpolation threshold */
  int		*interp_ytimeoutbuf;	/* interpolation timeout line buffer */
  int		interp_xtimeout;	/* interpolation timeout value in x */
  int		interp_ytimeout;	/* interpolation timeout value in y */
  }	picstruct;


/*------------------------------- preferences -------------------------------*/
typedef struct
  {
  char		prefs_name[MAXCHAR];			/* prefs filename*/
  char		*(image_name[2]);			/* image filenames */
  int		nimage_name;				/* nb of params */
  char		cat_name[MAXCHAR];			/* catalog filename*/
/*----- thresholding */
  enum	{T_SIGMA, T_MAGNITUDE}		threshold_type;	/* threshold unit */
  double	dthresh[2];				/* detect. threshold */
  int		ndthresh;				/* (1 or 2 entries) */
  double	thresh[2];				/* analysis thresh. */
  int		nthresh;				/* (1 or 2 entries) */
/*----- extraction */
  int		dimage_flag;				/* detect. image ? */
  int		ext_minarea;				/* min area in pix. */
  int		deb_maxarea;				/* max deblend. area */
  int		filter_flag;				/* smoothing on/off */
  char		filter_name[MAXCHAR];			/* mask filename */
  double	filter_thresh[2];			/* Filter thresholds */
  int		nfilter_thresh;				/* nb of params */
  int		deblend_nthresh;			/* threshold number */
  double	deblend_mincont;			/* minimum contrast */
  double	satur_level;				/* saturation level */
  enum	{CCD, PHOTO}			detect_type;	/* detection type */
/*----- Flagging */
  char		*(fimage_name[MAXFLAG]);		/* flagmap filenames */
  int		nfimage_name;				/* nb of params */
  enum	{FLAG_OR, FLAG_AND, FLAG_MIN, FLAG_MAX, FLAG_MOST}
				flag_type[MAXFLAG];	/* flag combination */
  int		imaflag_size;				/* requested iso nb1 */
  int		imanflag_size;				/* requested iso nb2 */
  int		nimaisoflag;				/* effective iso nb */
  int		nimaflag;				/* effective ima nb */
/*----- cleaning */
  int		clean_flag;				/* allow cleaning ? */
  double	clean_param;				/* cleaning effic. */
/*----- Weighting */
  char		*(wimage_name[2]);       		/* weight filenames */
  int		nwimage_name;				/* nb of params */
  weightenum	weight_type[2];				/* weighting scheme */
  int		nweight_type;				/* nb of params */
  int		weight_flag;				/* do we weight ? */
  int		dweight_flag;				/* detection weight? */

/*----- photometry */
  enum	{CAT_NONE, ASCII, ASCII_HEAD, ASCII_SKYCAT,
	FITS_LDAC, FITS_BINIMHEAD, FITS_NOIMHEAD, FITS_10}
		cat_type;				/* type of catalog */
  enum	{PNONE, FIXED, AUTO}		apert_type;	/* type of aperture */
  double	apert[MAXNAPER];			/* apert size (pix) */
  int		naper;					/* effective apert. */
  int		flux_apersize, fluxerr_apersize;	/* requested apert. */
  int		mag_apersize, magerr_apersize;		/* requested apert. */
  double	autoparam[2];				/* Kron parameters */
  int		nautoparam;				/* nb of Kron params */
  double	autoaper[2];				/* minimum apertures */
  int		nautoaper;				/* nb of min. aperts */
  double	mag_zeropoint;				/* magnitude offsets */
  double	mag_gamma;				/* for emulsions */
  double	gain;					/* only for CCD */
/*----- S/G separation */
  double	pixel_scale;				/* in arcsec */
  double	seeing_fwhm;				/* in arcsec */
  char		nnw_name[MAXCHAR];			/* nnw filename */
/*----- background */
  enum	{IMAGE, AFILE}			back_origin;	/* origin of bkgnd */
  char		back_name[MAXCHAR];			/* bkgnd filename */
  int		backsize[2];				/* bkgnd mesh size */
  int		nbacksize;				/* nb of params */
  int		backfsize[2];				/* bkgnd filt. size */
  int		nbackfsize;				/* nb of params */
  enum	{GLOBAL, LOCAL}			pback_type;	/* phot. bkgnd type */
  int		pback_size;				/* rect. ann. width */
/*----- memory */
  int		clean_stacksize;			/* size of buffer */
  int		mem_pixstack;				/* pixel stack size */
  int		mem_bufsize;				/* strip height */
/*----- catalog output */
  char		param_name[MAXCHAR];			/* param. filename */
/*----- miscellaneous */
  int		pipe_flag;				/* allow piping ? */
  enum	{QUIET, NORM, WARN, FULL}      	verbose_type;	/* display type */
/*----- CHECK-images */
  int		check_flag;				/* CHECK-image flag */
  checkenum    	check_type[MAXCHECK];		       	/* check-image types */
  int		ncheck_type;				/* nb of params */
  char		*(check_name[MAXCHECK]);       		/* check-image names */
  int		ncheck_name;				/* nb of params */
  struct structcheck	*(check[MAXCHECK]);		/* check-image ptrs */
/*----- ASSOCiation */
  int		assoc_flag;				/* ASSOCiation flag */
  char		assoc_name[MAXCHAR];			/* ASSOC-list name */
  int		assoc_param[3];				/* ASSOC param cols */
  int		nassoc_param;				/* nb of params */
  int		assoc_data[MAXNASSOC];		       	/* ASSOC data cols */
  int		nassoc_data;				/* nb of params */
  enum	{ASSOC_FIRST, ASSOC_MEAN, ASSOC_MAGMEAN,
	 ASSOC_SUM, ASSOC_MAGSUM, ASSOC_MIN, ASSOC_MAX}
		assoc_type;				/* type of assoc. */
  enum	{ASSOCSELEC_ALL, ASSOCSELEC_MATCHED, ASSOCSELEC_NOMATCHED}
		assocselec_type;		       	/* type of assoc. */
  double	assoc_radius;				/* ASSOC range */
  int		assoc_size;				/* nb of parameters */
  char		retina_name[MAXCHAR];			/* retina filename */
  int		vignetsize[2];				/* vignet size */
  int		vigshiftsize[2];			/* shift-vignet size */
  int		cleanmargin;				/* CLEANing margin */
  char		som_name[MAXCHAR];			/* SOM filename */
  int		somfit_flag;				/* ASSOCiation flag */
  int		somfit_vectorsize;			/* SOMfit vec. size */
/*----- masking */
  enum {MASK_NONE, MASK_BLANK, MASK_CORRECT}
		mask_type;				/* type of masking */
  int		blank_flag;				/* BLANKing flag */
  double	weight_thresh[2];      			/* weight threshlds */
  int		nweight_thresh;				/* nb of params */
/*----- interpolation */
  enum	{INTERP_NONE, INTERP_VARONLY, INTERP_ALL}
		interp_type[2];				/* interpolat. type */
  int		ninterp_type;				/* nb of params */
  int		interp_xtimeout[2];   			/* interp. x timeout */
  int		ninterp_xtimeout;       	        /* nb of params */
  int		interp_ytimeout[2];   			/* interp. y timeout */
  int		ninterp_ytimeout;       		/* nb of params */
/*----- astrometry */
  int		world_flag;				/* WORLD required */
/*----- growth curve */
  int		growth_flag;				/* gr. curve needed */
  int		flux_growthsize;       			/* number of elem. */
  int		mag_growthsize;       			/* number of elem. */
  int		flux_radiussize;       			/* number of elem. */
  double	growth_step;				/* step size (pix) */
  double	flux_frac[MAXNAPER];			/* for FLUX_RADIUS */
  int		nflux_frac;       			/* number of elem. */
/*----- PSF-fitting */
  char		psf_name[MAXCHAR];			/* PSF filename */
/*----- customize */
  double	mama_corflex;
  int		fitsunsigned_flag;			/* Force unsign FITS */
  }	prefstruct;


/*-------------------------------- catalog  ---------------------------------*/

typedef struct
  {
  int		ndetect;				/* nb of detections */
  int		ntotal;					/* Total object nb */
  int		nparam;					/* Nb of parameters */
/*----- Misc. strings defining the extraction */
  char		prefs_name[MAXCHAR];			/* Prefs filename*/
  char		image_name[MAXCHAR];			/* image filename*/
  char		nnw_name[MAXCHAR];			/* NNW name */
  char		filter_name[MAXCHAR];			/* Filter name */
  char		soft_name[MAXCHAR];			/* Sextractor version */
/*----- time */
  char		ext_date[12],ext_time[12];		/* date and time */
  double	ext_elapsed;				/* processing time */
  }		sexcatstruct;

