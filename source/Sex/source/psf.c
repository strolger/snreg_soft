 /*
 				psf.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN (IAP, Leiden observatory & ESO)
*
*	Contents:	Fit the PSF to a detection.
*
*	Last modify:	20/07/98
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	"define.h"
#include	"globals.h"
#include	"fitscat.h"
#include	"check.h"
#include	"poly.h"
#include	"psf.h"

/*------------------------------- variables ---------------------------------*/

static double	*data, *weight, *mat, *vmat, *wmat;
static PIXTYPE	*datah, *weighth, *checkmask;
static double	*psfmasks, *psfmaskx,*psfmasky;
extern keystruct	objkey[];
extern objstruct	outobj;
extern obj2struct	outobj2;

/********************************* psf_init **********************************/
/*
Allocate memory and stuff for the PSF-fitting.
*/
void	psf_init(psfstruct *psf)
  {
   double	*w, *px,*py, hw,hh, yb;
   int		width, height, npix, p, x,y;

  width = PSF_WIDTH;
  height = PSF_HEIGHT;
  npix = width*height;
  QMALLOC(psfmasks, double, npix);
/* Initialize x-moments, y-moments and weights */
  QMALLOC(weighth, PIXTYPE, npix);
  QMALLOC(weight, double, npix);
  QMALLOC(psfmaskx, double, npix);
  QMALLOC(psfmasky, double, npix);
  hw = (width-1)/2.0;
  hh = (height-1)/2.0;
  w = weight;
  px = psfmaskx;
  py = psfmasky;
  for (y=0; y<height; y++)
    {
    yb = y-hh;
    for (x=0; x<width; x++)
      {
      *(w++) = 1.0;
      *(px++) = x-hw;
      *(py++) = yb;
      }
    }

/* Allocate working space */
  QMALLOC(datah, PIXTYPE, npix);
  QMALLOC(data, double, npix);
  QMALLOC(mat, double, npix*PSF_NA);
  if (prefs.check[CHECK_SUBPROTOS] || prefs.check[CHECK_PROTOS])
    QMALLOC(checkmask, PIXTYPE, npix);

  return;
  }  


/********************************* psf_end ***********************************/
/*
Free memory occupied by the PSF-fitting stuff.
*/
void	psf_end(psfstruct *psf)
  {
   int	d, ndim;

  QFREE(psfmasks);
  QFREE(psfmaskx);
  QFREE(psfmasky);
  QFREE(datah);
  QFREE(data);
  QFREE(weighth);
  QFREE(weight);
  QFREE(data);
  QFREE(mat);	
  QFREE(vmat);	
  QFREE(wmat);	

  if (prefs.check[CHECK_SUBPROTOS] || prefs.check[CHECK_PROTOS])
    QFREE(checkmask);

  ndim = psf->poly->ndim;
  for (d=0; d<ndim; d++)
    free(psf->contextname[d]);
  free(psf->contextname);
  free(psf->context);
  free(psf->contexttyp);
  poly_end(psf->poly);
  free(psf->maskcomp);
  free(psf->maskloc);
  free(psf->masksize);
  free(psf);

  return;
  }


/********************************* psf_load *********************************/
/*
Read the PSF data from a FITS file.
*/
psfstruct	*psf_load(char *filename)
  {
   static objstruct	saveobj;
   static obj2struct	saveobj2;
   psfstruct		*psf;
   catstruct		*cat;
   tabstruct		*tab;
   keystruct		*key;
   char			*head, *ci,*co, str[MAXCHAR];
   int			ndim, dim[POLY_MAXDIM],
			i,k;

/* Open the cat (well it is not a "cat", but simply a FITS file */
  if (!(cat = read_cat(filename)))
    error(EXIT_FAILURE, "*Error*: PSF file not found: ", filename);

/* OK, we now allocate memory for the PSF structure itself */
  QCALLOC(psf, psfstruct, 1);

/* Store a short copy of the PSF filename */
  if (ci=strrchr(filename, '/'))
    strcpy(psf->name, ci+1);
  else
    strcpy(psf->name, filename);

/* Look if a polynomial description of the PSF is available */
  if (tab = name_to_tab(cat, "PSF_POLY", 0))
    {
/*-- So we have a variable PSF! */

    head = tab->headbuf;

/*-- Dimension of the polynomial */
    if (fitsread(head, "POLNAXIS", &ndim, H_INT,T_LONG) != RETURN_OK)
      goto headerror;
    if (ndim > POLY_MAXDIM)
        {
        sprintf(str, "*Error*: The POLNAXIS parameter in %s exceeds %d",
		psf->name, POLY_MAXDIM);
        error(EXIT_FAILURE, str, "");
        }

    QMALLOC(psf->contextname, char *, ndim);
    QMALLOC(psf->context, double *, ndim);
    QMALLOC(psf->contexttyp, t_type, ndim);
    QMALLOC(psf->contextoffset, double, ndim);
    QMALLOC(psf->contextscale, double, ndim);

/*-- We will have to use the outobj structs, so we first save their content */
    saveobj = outobj;
    saveobj2 = outobj2;
/*-- outobj's are used as FLAG arrays, so we initialize them to 0 */
    memset(&outobj, 0, sizeof(outobj));
    memset(&outobj2, 0, sizeof(outobj2));
    for (i=0; i<ndim; i++)
      {
/*---- Polynomial degrees */
      sprintf(str, "POLAXIS%1d", i+1);
      if (fitsread(head, str, &dim[i], H_INT,T_LONG) != RETURN_OK)
        goto headerror;

/*---- Contexts */
      QMALLOC(psf->contextname[i], char, 80);
      sprintf(str, "POLNAME%1d", i+1);
      if (fitsread(head,str,psf->contextname[i],H_STRING,T_STRING)!=RETURN_OK)
        goto headerror;
      if ((k = findkey(psf->contextname[i], (char *)objkey,
		sizeof(keystruct)))==RETURN_ERROR)
        {
        sprintf(str, "*Error*: %s CONTEXT parameter in %s unknown",
		psf->contextname[i], psf->name);
        error(EXIT_FAILURE, str, "");
        }
      key = objkey+k;
      psf->context[i] = key->ptr;
      psf->contexttyp[i] = key->ttype;
/*---- Declare the parameter "active", to trigger computation by SExtractor */
      *((char *)key->ptr) = (char)'\1';

/*---- Scaling of the context parameter */
      sprintf(str, "POLZERO%1d", i+1);
      if (fitsread(head, str, &psf->contextoffset[i], H_EXPO, T_DOUBLE)
		!=RETURN_OK)
        goto headerror;
      sprintf(str, "POLSCAL%1d", i+1);
      if (fitsread(head, str, &psf->contextscale[i], H_EXPO, T_DOUBLE)
		!=RETURN_OK)
        goto headerror;
      }

    psf->poly = poly_init(dim, ndim);
/*-- Load the PSF polynomial data */
    key = read_key(tab, "POLY_COEFFS");
    psf->poly->coeff = key->ptr;

/*-- Don't touch my arrays!! */
    blank_keys(tab);

/*-- Update the permanent FLAG arrays (that is, perform an "OR" on them) */
    for (ci=(char *)&outobj,co=(char *)&flagobj,i=sizeof(objstruct); i--;)
      *(co++) |= *(ci++);
    for (ci=(char *)&outobj2,co=(char *)&flagobj2,i=sizeof(obj2struct); i--;)
      *(co++) |= *(ci++);
/*-- Restore previous outobj contents */
    outobj = saveobj;
    outobj2 = saveobj2;
    }
  else
    {
/*-- This is a simple, constant PSF */
    *dim = 1;
    psf->poly = poly_init(dim, 1);
    psf->context = NULL;
    }

/* Now the PSF pixels themselves */
  if (!(tab = name_to_tab(cat, "PSF_MASK", 0)))
    error(EXIT_FAILURE, "*Error*: PSF_MASK table not found in catalog ",
	filename);

/* Load important scalars (which are stored as FITS keywords) */
  head = tab->headbuf;

/* Dimensionality of the PSF mask */
  if (fitsread(head, "PSFNAXIS", &psf->maskdim, H_INT, T_LONG) != RETURN_OK)
    goto headerror;
  if (psf->maskdim<2 || psf->maskdim>3)
    error(EXIT_FAILURE, "*Error*: wrong dimensionality for the PSF "
	"mask in ", filename);
  QMALLOC(psf->masksize, int, psf->maskdim);
  for (i=0; i<psf->maskdim; i++)
    psf->masksize[i] = 1;
  psf->masknpix = 1;
  for (i=0; i<psf->maskdim; i++)
    {
    sprintf(str, "PSFAXIS%1d", i+1);
    if (fitsread(head, str, &psf->masksize[i], H_INT,T_LONG) != RETURN_OK)
      goto headerror;
    psf->masknpix *= psf->masksize[i];
    }

/* Load the PSF mask data */
  key = read_key(tab, "PSF_MASK");
  psf->maskcomp = key->ptr;

  QMALLOC(psf->maskloc, double, psf->masksize[0]*psf->masksize[1]);

/* But don't touch my arrays!! */
  blank_keys(tab);

  free_cat(cat, 1);

  return psf;

headerror:
  error(EXIT_FAILURE, "*Error*: Incorrect or obsolete PSF data in ", filename);
  }


/******************************** psf_fit ***********************************/
void	psf_fit(psfstruct *psf, picstruct *field, picstruct *wfield,
		objstruct *obj)
  {
   checkstruct		*check;
   static double	x2[PSF_NPSF],y2[PSF_NPSF],xy[PSF_NPSF],sum[PSF_NPSF],
			deltax[PSF_NPSF],deltay[PSF_NPSF],
			sol[PSF_NA], covmat[PSF_NA*PSF_NA];
   double		*d, *m, *w, *ps, *px, *py, *vmat,*wmat,
			dx,dy, x1,y1, mx2,my2,mxy,m0,pix, wsum,wthresh,val,
			backnoise2, gain, radmin2,radmax2,satlevel, chi2;
   float		*dh, *wh;
   PIXTYPE		*cpix;
   int			i,j,p, npix, ix,iy,niter;

  gain = prefs.gain;
  backnoise2 = field->backsig*field->backsig;
  satlevel = prefs.satur_level - obj->bkg;
  wthresh = wfield?wfield->weight_thresh:BIG;
  radmin2 = PSF_MINSHIFT*PSF_MINSHIFT;
  radmax2 = PSF_MAXSHIFT*PSF_MAXSHIFT;
  npix = PSF_NPIX;
  ix = (int)(obj->mx+0.4999);
  iy = (int)(obj->my+0.4999);
  copyimage(field, datah, PSF_WIDTH, PSF_HEIGHT, ix, iy);

/* Compute weights */
  wsum = 0.0;
  if (wfield)
    {
    copyimage(wfield, weighth, PSF_WIDTH, PSF_HEIGHT, ix, iy);
    for (wh=weighth ,w=weight, p=npix; p--;)
      wsum += (*(w++) = (pix=*(wh++))<wthresh? sqrt(pix): 0.0);
    }
  else
    for (w=weight, dh=datah, p=npix; p--;)
      wsum += (*(w++) = ((pix = *(dh++))>-BIG && pix<satlevel)?
		1.0/sqrt(backnoise2+(pix>0.0?pix/gain:0.0))
		:0.0);

/* Special action if all the weights are zero!! */
  if (wsum == 0.0)
    {
    obj2->flux_psf = 0.0;
    obj2->x_psf = 999999.0;
    obj2->y_psf = 999999.0;
    obj2->niter_psf = 0;
    return;
    }

/* Weight the data */
  dh = datah;
  d = data;
  w = weight;
  for (p=npix; p--;)
    *(d++) = *(dh++)*(*(w++));

/* Get the local PSF */
  psf_build(psf);

/* Initialize PSF shifts */
  for (j=0; j<PSF_NPSF; j++)
    deltax[j] =deltay[j] = 0.0;

  niter = 0;
  for (i=0; i<PSF_NITER; i++)
    {
    niter++;
    m = mat;
    for (j=0; j<PSF_NPSF; j++)
      {
      dx = deltax[j];
      dy = deltay[j];
      psf_shift(psf, psfmasks, PSF_WIDTH, PSF_HEIGHT, dx, dy);
/*---- 0th order (weighted) */
      ps = psfmasks;
      w = weight;
      for (p=npix; p--;)
        *(m++) = *(ps++)**(w++);
/*---- 1st order (weighted) moment in x */
      ps = psfmasks;
      px = psfmaskx;
      w = weight;
      for (p=npix; p--;)
        *(m++) = *(ps++)*(*(px++)-dx)**(w++);
/*---- 1st order (weighted) moment in y */
      ps = psfmasks;
      py = psfmasky;
      w = weight;
      for (p=npix; p--;)
        *(m++) = *(ps++)*(*(py++)-dy)**(w++);
/*---- 2nd order moments */
      ps = psfmasks;
      px = psfmaskx;
      py = psfmasky;
      m0 = mx2 = my2 = mxy = 0.0;
      for (p=npix; p--;)
        {
        mxy += *ps*(x1=*(px++)-dx)*(y1=*(py++)-dy);
        mx2 += *ps*x1*x1;
        my2 += *ps*y1*y1;
        m0 += *(ps++);
        }
      x2[j] = mx2;
      y2[j] = my2;
      xy[j] = mxy;
      sum[j] = m0;
      }
    svdfit(mat, data, npix, PSF_NA, sol, &vmat, &wmat);
/*-- Update the PSF shifts */
    for (j=0; j<PSF_NPSF; j++)
      {
      dx = (sol[j*PSF_NA+1]*x2[j]+sol[j*PSF_NA+2]*xy[j])/sol[j*PSF_NA];
      dy = (sol[j*PSF_NA+2]*y2[j]+sol[j*PSF_NA+1]*xy[j])/sol[j*PSF_NA];
      deltax[j] += dx;
      deltay[j] += dy;
/*---- Exit if too much decentering */
      if ((dx*dx+dy*dy) < radmin2
	|| (deltax[j]*deltax[j] + deltay[j]*deltay[j]) > radmax2)
        i = PSF_NITER;
      }
    }

  dx = deltax[0];
  dy = deltay[0];
  psf_shift(psf, psfmasks, PSF_WIDTH, PSF_HEIGHT, dx, dy);
  obj2->flux_psf = sol[0];
  obj2->x_psf = ix+dx+1.0;
  obj2->y_psf = iy+dy+1.0;
  obj2->niter_psf = niter;

/* Compute chi2 if asked to */
  if (FLAG(obj2.chi2_psf))
    {
    chi2 = 0.0;
    val = obj2->flux_psf;
    for (d=data,w=weight,ps=psfmasks,p=npix; p--;)
      pix = *(d++) - *(ps++)**(w++)*val, chi2 += pix*pix;
    obj2->chi2_psf = chi2/(npix - 3);
    }

/* Error estimates on parameters */
  if (FLAG(obj2.fluxerr_psf))
    {
     double	*var, vara,varb,covab;

/*-- Compute variances and covariances if asked to */
    svdvar(vmat, wmat, PSF_NA, covmat);
    var = covmat;
    for (j=0; j<PSF_NPSF; j++)
      {
/*---- First, the error on the flux estimate */
      obj2->fluxerr_psf = sqrt(*var);
      if (FLAG(obj2.poserrmx2_psf))
        {
/*------ Variances and covariance along x and y */
        vara = *(var += PSF_NA+1);
        covab = *(++var);
        varb = *(var += PSF_NA);
        var += PSF_NA+1;
        obj2->poserrmx2_psf = (vara*x2[j]*x2[j]+varb*xy[j]*xy[j]
				+2*covab*x2[j]*xy[j])/(sol[0]*sol[0]);
        obj2->poserrmy2_psf = (varb*y2[j]*y2[j]+vara*xy[j]*xy[j]
				+2*covab*y2[j]*xy[j])/(sol[0]*sol[0]);
        obj2->poserrmxy_psf = (vara*x2[j]*xy[j]+varb*y2[j]*xy[j]
				+covab*(x2[j]*y2[j]+xy[j]*xy[j]))
				/(sol[0]*sol[0]);

/*------ If requested, translate variances to major and minor error axes... */
        if (FLAG(obj2.poserra_psf))
          {
           double	pmx2,pmy2,temp,theta;

          if (fabs(temp=obj2->poserrmx2_psf-obj2->poserrmy2_psf) > 0.0)
            theta = atan2(2.0 * obj2->poserrmxy_psf,temp) / 2.0;
          else
            theta = PI/4.0;

          temp = sqrt(0.25*temp*temp+obj2->poserrmxy_psf*obj2->poserrmxy_psf);
          pmy2 = pmx2 = 0.5*(obj2->poserrmx2_psf+obj2->poserrmy2_psf);
          pmx2+=temp;
          pmy2-=temp;

          obj2->poserra_psf = (float)sqrt(pmx2);
          obj2->poserrb_psf = (float)sqrt(pmy2);
          obj2->poserrtheta_psf = theta*180.0/PI;
          }

/*------ ...Or ellipse parameters */
        if (FLAG(obj2.poserr_cxx))
          {
           double	xm2,ym2, xym, temp;

          xm2 = obj2->poserrmx2_psf;
          ym2 = obj2->poserrmy2_psf;
          xym = obj2->poserrmxy_psf;
          obj2->poserrcxx_psf = (float)(ym2/(temp=xm2*ym2-xym*xym));
          obj2->poserrcyy_psf = (float)(xm2/temp);
          obj2->poserrcxy_psf = (float)(-2*xym/temp);
          }
        }
      else
        var += 3*PSF_NA+3;
      }
    }

  if (prefs.check[CHECK_SUBPROTOS] || prefs.check[CHECK_PROTOS])
    for (j=0; j<PSF_NPSF; j++)
      {
      cpix = checkmask;
      ps = psfmasks;
      for (p=npix; p--;)
        *(cpix++) = (PIXTYPE)*(ps++);
      if (check = prefs.check[CHECK_SUBPROTOS])
        addcheck(check, checkmask,PSF_WIDTH,PSF_HEIGHT,ix,iy,-obj2->flux_psf);
      if (check = prefs.check[CHECK_PROTOS])
        addcheck(check, checkmask,PSF_WIDTH,PSF_HEIGHT,ix,iy,obj2->flux_psf);
      }

  return;
  }


/******************************* psf_build **********************************/
/*
Build the local PSF (function of "context").
*/
void	psf_build(psfstruct *psf)
  {
   static double	pos[POLY_MAXDIM];
   double	*pl, *basis, fac;
   float	*ppc;
   int		i,n,p, ndim, npix;

  npix = psf->masksize[0]*psf->masksize[1];
/* Reset the Local PSF mask */
  memset(psf->maskloc, 0.0, npix*sizeof(double));

/* Grab the context vector */
  ndim = psf->poly->ndim;
  for (i=0; i<ndim; i++)
    {
    pos[i] = (*(double *)ttypeconv(psf->context[i],psf->contexttyp[i],T_DOUBLE)
	- psf->contextoffset[i]) / psf->contextscale[i];
    }
  poly_func(psf->poly, pos);

  basis = psf->poly->basis;

  ppc = psf->maskcomp;
/* Sum each component */
  for (n = (psf->maskdim>2?psf->masksize[2]:1); n--;)
    {
    pl = psf->maskloc;
    fac = *(basis++);
    for (p=npix; p--;)
      *(pl++) +=  fac**(ppc++);
    }

  return;
  }


/******************************* psf_shift **********************************/
/*
Shift the PSF from an arbitrary amount.
*/
int	psf_shift(psfstruct *psf, double *pix2, int w2, int h2,
		double dx, double dy)
  {
   double	*mask,*maskt, *pix1, *pix12, *pixin,*pixin0, *pixout,*pixout0,
		xc1,xc2,yc1,yc2, xs1,ys1, x1,y1, x,y, dxm,dym,
		val,step2;
   int		i,j,k,n,t, *start,*startt, *nmask,*nmaskt,
		ixs2,iys2, ix2,iy2, dix2,diy2, nx2,ny2, iys1a, ny1, hmw,hmh,
		ix,iy, ix1,iy1,w1,h1;

  step2 = 1.0;
  pix1 = psf->maskloc;
  w1 = psf->masksize[0];
  h1 = psf->masksize[1];

  xc1 = (double)((w1-1)/2);	/* Im1 center x-coord*/
  xc2 = (double)((w2-1)/2);	/* Im2 center x-coord*/
  xs1 = xc1 - dx - xc2*step2;	/* Im1 start x-coord */
  if ((int)xs1 >= w1)
    return RETURN_ERROR;
  ixs2 = 0;			/* Int part of Im2 start x-coord */
  if (xs1<0.0)
    {
    dix2 = (int)(1-xs1/step2);
/*-- Simply leave here if the images do not overlap in x */
    if (dix2 >= w2)
      return RETURN_ERROR;
    ixs2 += dix2;
    xs1 += dix2*step2;
    }
  nx2 = (int)((w1-xs1)/step2);/* nb of interpolated Im2 pixels along x */
  if (nx2>(ix2=w2-ixs2))
    nx2 = ix2;
  yc1 = (double)((h1-1)/2);	/* Im1 center y-coord */
  yc2 = (double)((h2-1)/2);	/* Im2 center y-coord */
  ys1 = yc1 - dy - yc2*step2;	/* Im1 start y-coord */
  if ((int)ys1 >= h1)
    return RETURN_ERROR;
  iys2 = 0;			/* Int part of Im2 start y-coord */
  if (ys1<0.0)
    {
    diy2 = (int)(1-ys1/step2);
/*-- Simply leave here if the images do not overlap in y */
    if (diy2 >= h2)
      return RETURN_ERROR;
    iys2 += diy2;
    ys1 += diy2*step2;
    }
  ny2 = (int)((h1-ys1)/step2);/* nb of interpolated Im2 pixels along y */
  if (ny2>(iy2=h2-iys2))
    ny2 = iy2;

/* Set the yrange for the x-resampling with some margin for interpolation */
  iys1a = (int)ys1;		/* Int part of Im1 start y-coord with margin */
  hmh = INTERPH/2 - 1;		/* Interpolant start */
  if (iys1a<0 || ((iys1a -= hmh)< 0))
    iys1a = 0;
  ny1 = (int)(ys1+ny2*step2)+INTERPW-hmh;	/* Interpolated Im1 y size */
  if (ny1>h1)					/* with margin */
    ny1 = h1;
/* Express everything relative to the effective Im1 start (with margin) */
  ny1 -= iys1a;
  ys1 -= (double)iys1a;

/* Allocate interpolant stuff for the x direction */
  QMALLOC(mask, double, nx2*INTERPW);	/* Interpolation masks */
  QMALLOC(nmask, int, nx2);		/* Interpolation mask sizes */
  QMALLOC(start, int, nx2);		/* Int part of Im1 conv starts */
/* Compute the local interpolant and data starting points in x */
  hmw = INTERPW/2 - 1;
  x1 = xs1;
  maskt = mask;
  nmaskt = nmask;
  startt = start;
  for (j=nx2; j--; x1+=step2)
    {
    ix = (ix1=(int)x1) - hmw;
    dxm = ix1 - x1 - hmw;	/* starting point in the interpolation func */
    if (ix < 0)
      {
      n = INTERPW+ix;
      dxm -= (double)ix;
      ix = 0;
      }
    else
      n = INTERPW;
    if (n>(t=w1-ix))
      n=t;
    *(startt++) = ix;
    *(nmaskt++) = n;
    for (x=dxm, i=n; i--; x+=1.0)
      *(maskt++) = INTERPF(x);
    }

  QCALLOC(pix12, double, nx2*ny1);	/* Intermediary frame-buffer */

/* Make the interpolation in x (this includes transposition) */
  pixin0 = pix1+iys1a*w1;
  pixout0 = pix12;
  for (k=ny1; k--; pixin0+=w1, pixout0++)
    {
    maskt = mask;
    nmaskt = nmask;
    startt = start;
    pixout = pixout0;
    for (j=nx2; j--; pixout+=ny1)
      {
      pixin = pixin0+*(startt++);
      val = 0.0; 
      for (i=*(nmaskt++); i--;)
        val += *(maskt++)**(pixin++);
      *pixout = val;
      }
    }

/* Reallocate interpolant stuff for the y direction */
  QREALLOC(mask, double, ny2*INTERPH);	/* Interpolation masks */
  QREALLOC(nmask, int, ny2);		/* Interpolation mask sizes */
  QREALLOC(start, int, ny2);		/* Int part of Im1 conv starts */

/* Compute the local interpolant and data starting points in y */
  hmh = INTERPH/2 - 1;
  y1 = ys1;
  maskt = mask;
  nmaskt = nmask;
  startt = start;
  for (j=ny2; j--; y1+=step2)
    {
    iy = (iy1=(int)y1) - hmh;
    dym = iy1 - y1 - hmh;	/* starting point in the interpolation func */
    if (iy < 0)
      {
      n = INTERPH+iy;
      dym -= (double)iy;
      iy = 0;
      }
    else
      n = INTERPH;
    if (n>(t=ny1-iy))
      n=t;
    *(startt++) = iy;
    *(nmaskt++) = n;
    for (y=dym, i=n; i--; y+=1.0)
      *(maskt++) = INTERPF(y);
    }

/* Initialize destination buffer to zero */
  memset(pix2, 0, w2*h2*sizeof(double));

/* Make the interpolation in y  and transpose once again */
  pixin0 = pix12;
  pixout0 = pix2+ixs2+iys2*w2;
  for (k=nx2; k--; pixin0+=ny1, pixout0++)
    {
    maskt = mask;
    nmaskt = nmask;
    startt = start;
    pixout = pixout0;
    for (j=ny2; j--; pixout+=w2)
      {
      pixin = pixin0+*(startt++);
      val = 0.0; 
      for (i=*(nmaskt++); i--;)
        val += *(maskt++)**(pixin++);
      *pixout = val;
      }
    }

/* Free memory */
  free(pix12);
  free(mask);
  free(nmask);
  free(start);

  return RETURN_OK;
  }


/******************************** svdfit ************************************/
/*
General least-square fit A.x = b, based on Singular Value Decomposition (SVD).
Loosely adapted from Numerical Recipes in C, 2nd Ed. (p. 671).
Note: the a and v matrices are transposed with respect to the N.R. convention.
*/
void svdfit(double *a, double *b, int m, int n, double *sol,
	double **vout, double **wout)
  {
#define MAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))
#define	PYTHAG(a,b)	((at=fabs(a)) > (bt=fabs(b)) ? \
				  (ct=bt/at,at*sqrt(1.0+ct*ct)) \
				: (bt ? (ct=at/bt,bt*sqrt(1.0+ct*ct)): 0.0))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define	TOL		1.0e-11

   int			flag,i,its,j,jj,k,l,nm,mmi,nml;
   static double	rv1[PSF_NA], tmp[PSF_NA], vmat[PSF_NA*PSF_NA],
			wmat[PSF_NA];
   double		c,f,h,s,x,y,z,
			anorm, g, scale,
			at,bt,ct,maxarg1,maxarg2,
			thresh, wmax,
			*w,*ap,*ap0,*ap1,*ap10,*rv1p,*vp,*vp0,*vp1,*vp10,
			*bp,*tmpp;

  anorm = g = scale = 0.0;
  if (m < n)
    error(EXIT_FAILURE, "*Error*: Not enough rows for solving the system ",
	"in svdfit()");
  
  for (i=0;i<n;i++)
    {
    l = i+1;
    nml = n-l;
    rv1[i] = scale*g;
    g = s = scale = 0.0;
    if ((mmi = m - i) > 0)
      {
      ap = ap0 = a+i*(m+1);
      for (k=mmi;k--;)
        scale += fabs(*(ap++));
      if (scale)
        {
        for (ap=ap0,k=mmi; k--; ap++)
          {
          *ap /= scale;
          s += *ap**ap;
          }
        f = *ap0;
        g = -SIGN(sqrt(s),f);
        h = f*g-s;
        *ap0 = f-g;
        ap10 = a+l*m+i;
        for (j=nml;j--; ap10+=m)
          {
          for (s=0.0,ap=ap0,ap1=ap10,k=mmi; k--;)
            s += *(ap1++)**(ap++);
          f = s/h;
          for (ap=ap0,ap1=ap10,k=mmi; k--;)
            *(ap1++) += f**(ap++);
          }
        for (ap=ap0,k=mmi; k--;)
          *(ap++) *= scale;
        }
      }
    wmat[i] = scale*g;
    g = s = scale = 0.0;
    if (i < m && i+1 != n)
      {
      ap = ap0 = a+i+m*l;
      for (k=nml;k--; ap+=m)
        scale += fabs(*ap);
      if (scale)
        {
        for (ap=ap0,k=nml;k--; ap+=m)
          {
          *ap /= scale;
          s += *ap**ap;
          }
        f=*ap0;
        g = -SIGN(sqrt(s),f);
        h=f*g-s;
        *ap0=f-g;
        rv1p = rv1+l;
        for (ap=ap0,k=nml;k--; ap+=m)
          *(rv1p++) = *ap/h;
        ap10 = a+l+m*l;
        for (j=m-l; j--; ap10++)
          {
          for (s=0.0,ap=ap0,ap1=ap10,k=nml; k--; ap+=m,ap1+=m)
            s += *ap1**ap;
          rv1p = rv1+l;
          for (ap1=ap10,k=nml;k--; ap1+=m)
            *ap1 += s**(rv1p++);
          }
        for (ap=ap0,k=nml;k--; ap+=m)
          *ap *= scale;
        }
      }
    anorm=MAX(anorm,(fabs(wmat[i])+fabs(rv1[i])));
    }

  for (i=n-1;i>=0;i--)
    {
    if (i < n-1)
      {
      if (g)
        {
        ap0 = a+l*m+i;
        vp0 = vmat+i*n+l;
        vp10 = vmat+l*n+l;
        g *= *ap0;
        for (ap=ap0,vp=vp0,j=nml; j--; ap+=m)
          *(vp++) = *ap/g;
        for (j=nml; j--; vp10+=n)
          {
          for (s=0.0,ap=ap0,vp1=vp10,k=nml; k--; ap+=m)
            s += *ap**(vp1++);
          for (vp=vp0,vp1=vp10,k=nml; k--;)
            *(vp1++) += s**(vp++);
          }
        }
      vp = vmat+l*n+i;
      vp1 = vmat+i*n+l;
      for (j=nml; j--; vp+=n)
        *vp = *(vp1++) = 0.0;
      }
    vmat[i*n+i]=1.0;
    g=rv1[i];
    l=i;
    nml = n-l;
    }

  for (i=(m<n?m:n); --i>=0;)
    {
    l=i+1;
    nml = n-l;
    mmi=m-i;
    g=wmat[i];
    ap0 = a+i*m+i;
    ap10 = ap0 + m;
    for (ap=ap10,j=nml;j--;ap+=m)
      *ap=0.0;
    if (g)
      {
      g=1.0/g;
      for (j=nml;j--; ap10+=m)
        {
        for (s=0.0,ap=ap0,ap1=ap10,k=mmi; --k;)
              s += *(++ap)**(++ap1);
        f = (s/(*ap0))*g;
        for (ap=ap0,ap1=ap10,k=mmi;k--;)
          *(ap1++) += f**(ap++);
        }
      for (ap=ap0,j=mmi;j--;)
        *(ap++) *= g;
      }
    else
      for (ap=ap0,j=mmi;j--;)
        *(ap++)=0.0;
    ++(*ap0);
    }

  for (k=n; --k>=0;)
      {
      for (its=0;its<100;its++)
        {
        flag=1;
        for (l=k;l>=0;l--)
          {
          nm=l-1;
          if (fabs(rv1[l])+anorm == anorm)
            {
            flag=0;
            break;
            }
          if (fabs(wmat[nm])+anorm == anorm)
            break;
          }
        if (flag)
          {
          c=0.0;
          s=1.0;
          ap0 = a+nm*m;
          ap10 = a+l*m;
          for (i=l; i<=k; i++,ap10+=m)
            {
            f=s*rv1[i];
            if (fabs(f)+anorm == anorm)
              break;
            g=wmat[i];
            h=PYTHAG(f,g);
            wmat[i]=h;
            h=1.0/h;
            c=g*h;
            s=(-f*h);
            for (ap=ap0,ap1=ap10,j=m; j--;)
              {
              z = *ap1;
              y = *ap;
              *(ap++) = y*c+z*s;
              *(ap1++) = z*c-y*s;
              }
            }
          }
        z=wmat[k];
        if (l == k)
          {
          if (z < 0.0)
            {
            wmat[k] = -z;
            vp = vmat+k*n;
            for (j=n; j--; vp++)
              *vp = (-*vp);
            }
          break;
          }
        if (its == 99)
          error(EXIT_FAILURE, "*Error*: No convergence in 100 SVD iterations ",
		"in svdfit()");
        x=wmat[l];
        nm=k-1;
        y=wmat[nm];
        g=rv1[nm];
        h=rv1[k];
        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
        g=PYTHAG(f,1.0);
        f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
        c=s=1.0;
        ap10 = a+l*m;
        vp10 = vmat+l*n;
        for (j=l;j<=nm;j++,ap10+=m,vp10+=n)
          {
          i=j+1;
          g=rv1[i];
          y=wmat[i];
          h=s*g;
          g=c*g;
          z=PYTHAG(f,h);
          rv1[j]=z;
          c=f/z;
          s=h/z;
          f=x*c+g*s;
          g=g*c-x*s;
          h=y*s;
          y=y*c;
          for (vp=(vp1=vp10)+n,jj=n; jj--;)
            {
            z = *vp;
            x = *vp1;
            *(vp1++) = x*c+z*s;
            *(vp++) = z*c-x*s;
            }
          z=PYTHAG(f,h);
          wmat[j]=z;
          if (z)
            {
            z=1.0/z;
            c=f*z;
            s=h*z;
            }
          f=c*g+s*y;
          x=c*y-s*g;
          for (ap=(ap1=ap10)+m,jj=m; jj--;)
            {
            z = *ap;
            y = *ap1;
            *(ap1++) = y*c+z*s;
            *(ap++) = z*c-y*s;
            }
          }
        rv1[l]=0.0;
        rv1[k]=f;
        wmat[k]=x;
        }
      }

  wmax=0.0;
  w = wmat;
  for (j=n;j--; w++)
    if (*w > wmax)
      wmax=*w;
  thresh=TOL*wmax;
  w = wmat;
  for (j=n;j--; w++)
    if (*w < thresh)
      *w = 0.0;

  w = wmat;
  ap = a;
  tmpp = tmp;
  for (j=n; j--; w++)
    {
    s=0.0;
    if (*w)
      {
      bp = b;
      for (i=m; i--;)
        s += *(ap++)**(bp++);
      s /= *w;
      }
    else
      ap += m;
    *(tmpp++) = s;
    }

  vp0 = vmat;
  for (j=0; j<n; j++,vp0++)
    {
    s=0.0;
    tmpp = tmp;
    for (vp=vp0,jj=n; jj--; vp+=n)
      s += *vp**(tmpp++);
    sol[j]=s;
    }

/* Return pointers to the vmat and wmat matrices if requested (for covar.) */
  if (vout)
    *vout = vmat;
  if (wout)
    *wout = wmat;

  return;
  }

#undef SIGN
#undef MAX
#undef PYTHAG
#undef TOL

/******************************** svdvar ************************************/
/*
Computation of the covariance matrix from the SVD vmat and wmat matrices.A
dapted from Numerical Recipes in C, 2nd Ed. (p. 679).
*/
void svdvar(double *v, double *w, int n, double *cov)
  {
   static double	wti[PSF_NA];
   double		sum;
   int			i,j,k;

  for (i=0; i<n; i++)
    wti[i] = w[i]? 1.0/(w[i]*w[i]) : 0.0;

  for (i=0; i<n; i++)
    for (j=0; j<=i; j++)
      {
      for (sum=0.0,k=0; k<n; k++)
        sum += v[k*n+i]*v[k*n+j]*wti[k];
      cov[j*n+i] = cov[i*n+j] = sum;
      }

  return;
  }

