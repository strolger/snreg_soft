/*
 				makeit.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN, IAP & Leiden observatory
*
*	Contents:	main program.
*
*	Last modify:	13/07/98
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<time.h>

#include	"define.h"
#include	"globals.h"
#include	"fitscat.h"
#include	"check.h"
#include	"field.h"
#include	"filter.h"
#include	"growth.h"
#include	"interpolate.h"
#include	"psf.h"
#include	"retina.h"
#include	"som.h"
#include	"weight.h"

/******************************** makeit *************************************/
/*
Manage the whole stuff.
*/
void	makeit()

  {
   checkstruct		*check;
   picstruct		*dfield, *field,*pffield[MAXFLAG], *wfield,*dwfield;
   static time_t        thetime1, thetime2;
   int			i;

/*Initial time measurement*/
  time(&thetime1);

  dfield = field = wfield = dwfield = NULL;

/* Initialize globals variables */
  initglob();

  NFPRINTF(OUTPUT, "Setting catalog parameters");
  readcatparams(prefs.param_name);

  if (FLAG(obj2.flux_psf))
    {
    NFPRINTF(OUTPUT, "Reading PSF information");
    thepsf = psf_load(prefs.psf_name); 
    }

  useprefs();			/* update things accor. to prefs parameters */

  if (prefs.filter_flag)
    {
    NFPRINTF(OUTPUT, "Reading detection filter");
    getfilter(prefs.filter_name);	/* get the detection filter */
    }

  if (FLAG(obj2.sprob))
    {
    NFPRINTF(OUTPUT, "Initializing Neural Network");
    neurinit();
    NFPRINTF(OUTPUT, "Reading Neural Network Weights");
    getnnw(); 
    }

  if (prefs.somfit_flag)
    {
     int	margin;

    thesom = som_load(prefs.som_name);
    if ((margin=(thesom->inputsize[1]+1)/2) > prefs.cleanmargin)
      prefs.cleanmargin = margin;
    if (prefs.somfit_vectorsize>thesom->neurdim)
      {
      prefs.somfit_vectorsize = thesom->neurdim;
      sprintf(gstr,"%d", prefs.somfit_vectorsize);
      warning("Dimensionality of the SOM-fit vector limited to ", gstr);
      }
    }

/* Prepare growth-curve buffer */
  if (prefs.growth_flag)
    initgrowth();

  if (prefs.dimage_flag)
    {
    dfield = newfield(prefs.image_name[0], DETECT_FIELD);
    field = newfield(prefs.image_name[1], MEASURE_FIELD);
    if ((field->width!=dfield->width) || (field->height!=dfield->height))
      error(EXIT_FAILURE, "*Error*: Frames have different sizes","");
/*-- Prepare interpolation */
    if (prefs.weight_flag && prefs.interp_type[0] == INTERP_ALL)
      init_interpolate(dfield, -1, -1);
    if (prefs.interp_type[1] == INTERP_ALL)
      init_interpolate(field, -1, -1);
/*-- Copy field structures to static ones (for catalog info) */
    thefield1 = *field;
    thefield2 = *dfield;
    }
  else
    {
    field = newfield(prefs.image_name[0], DETECT_FIELD | MEASURE_FIELD);
/*-- Prepare interpolation */
    if (prefs.weight_flag && prefs.interp_type[0] == INTERP_ALL)
      init_interpolate(field, -1, -1);       /* 0.0 or anything else */
    thefield1 = thefield2 = *field;
    }

  for (i=0; i<prefs.nimaflag; i++)
    {
    pffield[i] = newfield(prefs.fimage_name[i], FLAG_FIELD);
    if ((pffield[i]->width!=field->width)
	|| (pffield[i]->height!=field->height))
      error(EXIT_FAILURE,
	"*Error*: Incompatible FLAG-map size in ", prefs.fimage_name[i]);
    }

/* Init the WEIGHT-images */
  if (prefs.weight_flag) 
    {
     weightenum	wtype;
     PIXTYPE	interpthresh;

    if (prefs.nweight_type>1)
      {
/*---- Double-weight-map mode */
      if (prefs.weight_type[1] != WEIGHT_NONE)
        {
/*------ First: the "measurement" weights */
        wfield = newweight(prefs.wimage_name[1],field,prefs.weight_type[1]);
        wtype = prefs.weight_type[1];
        interpthresh = prefs.weight_thresh[1];
/*------ Convert the interpolation threshold to variance units */
        if (wtype == WEIGHT_FROMWEIGHTMAP)
          weight_to_var(wfield, &interpthresh, 1);
        else if (wtype == WEIGHT_FROMRMSMAP || wtype == WEIGHT_FROMBACK)
          rms_to_var(wfield, &interpthresh, 1);
        wfield->weight_thresh = interpthresh;
        if (prefs.interp_type[1] != INTERP_NONE)
          init_interpolate(wfield,
		prefs.interp_xtimeout[1], prefs.interp_ytimeout[1]);
        }
/*---- The "detection" weights */
      if (prefs.weight_type[0] != WEIGHT_NONE)
        {
        if (prefs.weight_type[0] == WEIGHT_FROMINTERP)
          dwfield=newweight(prefs.wimage_name[0],wfield,prefs.weight_type[0]);
        else
          dwfield = newweight(prefs.wimage_name[0], dfield?dfield:field,
		prefs.weight_type[0]);
/*------ Interpolated weight-fields lack info about the orig. weight-type */
        wtype = (prefs.weight_type[0] == WEIGHT_FROMINTERP)?
		prefs.weight_type[1]: prefs.weight_type[0];
        interpthresh = prefs.weight_thresh[0];
/*------ Convert the interpolation threshold to variance units */
        if (wtype == WEIGHT_FROMWEIGHTMAP)
          weight_to_var(dwfield, &interpthresh, 1);
        else if (wtype == WEIGHT_FROMRMSMAP || wtype == WEIGHT_FROMBACK)
          rms_to_var(dwfield, &interpthresh, 1);
        dwfield->weight_thresh = interpthresh;
        if (prefs.interp_type[0] != INTERP_NONE)
          init_interpolate(dwfield,
		prefs.interp_xtimeout[0], prefs.interp_ytimeout[0]);
        }
      }
    else
      {
/*---- Single-weight-map mode */
      wfield = newweight(prefs.wimage_name[0], dfield?dfield:field,
			prefs.weight_type[0]);
      wtype = prefs.weight_type[0];
      interpthresh = prefs.weight_thresh[0];
/*---- Convert the interpolation threshold to variance units */
      if (wtype == WEIGHT_FROMWEIGHTMAP)
        weight_to_var(wfield, &interpthresh, 1);
      else if (wtype == WEIGHT_FROMRMSMAP || wtype == WEIGHT_FROMBACK)
        rms_to_var(wfield, &interpthresh, 1);
      wfield->weight_thresh = interpthresh;
      if (prefs.interp_type[0] != INTERP_NONE)
        init_interpolate(wfield,
		prefs.interp_xtimeout[0], prefs.interp_ytimeout[0]);
      }
    }

/* Init the CHECK-images */
  if (prefs.check_flag)
   {
    checkenum	c;

   NFPRINTF(OUTPUT, "Initializing check-image(s)");
    for (i=0; i<prefs.ncheck_type; i++)
      if ((c=prefs.check_type[i]) != CHECK_NONE)
        {
        if (prefs.check[c])
          error(EXIT_FAILURE, "*Error*: 2 CHECK_IMAGEs cannot have the same ",
			" CHECK_IMAGE_TYPE");
        prefs.check[c] = initcheck(field,prefs.check_name[i],
				prefs.check_type[i]);
        free(prefs.check_name[i]);
        }
   }

/* Initialize PSF workspace */
  if (FLAG(obj2.flux_psf))
    psf_init(thepsf);

  NFPRINTF(OUTPUT, "Initializing catalog");
  initcat(field);

/* Start the extraction pipeline */
  NFPRINTF(OUTPUT, "Scanning image");
  scanimage(field, dfield, pffield, prefs.nimaflag, wfield, dwfield);

  NFPRINTF(OUTPUT, "Closing files");

/* End the CHECK-image processing */
  if (prefs.check_flag)
    for (i=0; i<MAXCHECK; i++)
      if (check=prefs.check[i])
        endcheck(field, check);

/*Final time measurements*/
  if (time(&thetime2)!=-1)
    {
    if (!strftime(cat.ext_date, 10, "%d/%m/%y", localtime(&thetime2)))
      error(EXIT_FAILURE, "*Internal Error*: Date string too long ","");
    if (!strftime(cat.ext_time, 10, "%H:%M:%S", localtime(&thetime2)))
      error(EXIT_FAILURE, "*Internal Error*: Time/date string too long ","");
#ifdef SUN_OS
    cat.ext_elapsed = (double)thetime2 - (double)thetime1;
#else
    cat.ext_elapsed = difftime(thetime2, thetime1);
#endif
    }

  endcat();

  for (i=0; i<prefs.nimaflag; i++)
    endfield(pffield[i]);
  endfield(field);
  if (dfield)
    endfield(dfield);
  if (wfield)
    endfield(wfield);
  if (dwfield)
    endfield(dwfield);

  if (prefs.filter_flag)
    endfilter();

  if (prefs.somfit_flag)
    som_end(thesom);

  if (prefs.growth_flag)
    endgrowth();

  if (FLAG(obj2.flux_psf))
    psf_end(thepsf);

  if (FLAG(obj2.sprob))
    neurclose();

  if (prefs.verbose_type != QUIET)
    fprintf(OUTPUT, "Objects: detected %-8d / sextracted %-8d\n",
	cat.ndetect, cat.ntotal);
/*
  if (FLAG(obj.retinout))
    endretina(theretina);
*/
  return;
  }


/******************************** initglob ***********************************/
/*
Initialize a few global variables
*/
void	initglob()
  {
   int	i;

  for (i=0; i<37; i++)
    {
    ctg[i] = cos(i*PI/18);
    stg[i] = sin(i*PI/18);
    }


  return;
  }

/*
int matherr(struct exception *x)
{
printf("***MATH ERROR***: %d %s %f %f\n",
x->type, x->name, x->arg1, x->retval);
return (0);
}

*/