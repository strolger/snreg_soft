/* Various utility routines for superflat.c
 * John Tonry Oct 1999
 */

#include <stdio.h>
#include <math.h>
#include "superflat.h"

#define EXTRAHEAD (32)	    /* Extra header lines allocated for insertion */

#define NFITS (2880)
#define MAXBUF (65536)

static int buf[MAXBUF];

static float *srtbuf=NULL;
static int nsrtbuf=0;

/* Return a mode over a range of counts */
/* NOTE: assumes medcount has just been run previously */
modecount(ctr, sig)
int ctr;			/* Central data value */
int sig;			/* 1-sigma deviation */
{
  int i, m, p, c, mode, isig;
  double mdble;
  isig = MAX(1, (int)(0.5*sig+0.5));
  
  for(i=ctr-3*isig-1, m=0; i<=ctr-isig-1;   i++) m += buf[i];
  for(i=ctr-isig,     c=0; i<=ctr+isig;     i++) c += buf[i];
  for(i=ctr+isig+1,   p=0; i<=ctr+3*isig+1; i++) p += buf[i];
  mdble = ((double)(m-p) / (2*(m+p)-4*c)) * (2*isig+1) + ctr;
  if(mdble < 0) mode = (int)(mdble-0.5);
  else 		mode = (int)(mdble+0.5);
  return(mode);
}

/* Return a median over a range of counts */
/* NOTE: assumes medcount has just been run previously */
mclipcount(n1, n2)
int n1;			/* Lower data value */
int n2;			/* Upper data value */
{
  int i, j, median, ncut;
  for(i=n1, ncut=0; i<=n2; i++) ncut += buf[i];
  for(median=n1, i=0; median<=n2; median++) {
    i += buf[median];
    if(i > ncut/2) break;
  }
  return(median);
}

/* Return a median derived from counting the specified subarray */
medcount(sx, sy, nx, ny, frame, sig)
int sx, sy, nx, ny;
struct ushortfits *frame;	/* image data structure */
int *sig;			/* 1-sigma of distribution */
{
  int i, j, ntot, median, nsx;

  for(i=0; i<65536; i++) buf[i] = 0;
  nsx = frame->nx;
  for(j=sy, ntot=0; j<sy+ny; j++) {
    for(i=sx; i<sx+nx; i++) {
      if((frame->data)[i+nsx*j] > 0) {
	buf[(frame->data)[i+nsx*j]]++;
	ntot++;
      }
    }
  }
  for(median=1, i=0; median<65536; median++) {
    i += buf[median];
    if(i < ntot/6) j = median;
    if(i > ntot/2) break;
  }
  *sig = MAX(median-j, 1);

  return(median);
}

/* Return a median which is derived from sorting the specified subarray */
medsort(sx, sy, nx, ny, frame)
int sx, sy, nx, ny;
struct ushortfits *frame;	/* image data structure */
{
  int i, j, nsrt;
  if(nx*ny > nsrtbuf) {
    if(srtbuf != NULL) free(srtbuf);
    srtbuf = (float *)calloc(nx*ny, sizeof(float));
  }
  for(j=sy, nsrt=0; j<sy+ny; j++) {
    for(i=sx; i<sx+nx; i++) {
      if((frame->data)[i+frame->nx*j] > 0) {
	srtbuf[nsrt++] = (frame->data)[i+frame->nx*j];
      }
    }
  }
  qsort4(nsrt, srtbuf);
  return((int)(0.5*(srtbuf[(nsrt-1)/2] + srtbuf[nsrt/2])));
}

/*
 * copy all but the data, and allocate space for the data
 */
copyframe(src, dest, npix)
     struct ushortfits *src, **dest;	/* image data structure */
     int npix;				/* Number of data pixels to alloc */
{
  int i, j, k, l, sum, nsx, idx;

/* Create a destination frame */
  if(*dest == NULL) newframe(dest);
  else eraseframe(dest);

/* Copy everything but the data */
/* Note that the header is *NOT* kept consistent */
  (*dest)->fname = (char *)malloc(strlen(src->fname)+1);
  strcpy((*dest)->fname, src->fname);

  (*dest)->header = (char *)malloc(80*(src->nhead+EXTRAHEAD)+1);
  (*dest)->headstore = 80*(src->nhead+EXTRAHEAD) + 1;
  strncpy((*dest)->header, src->header, 80*(src->nhead+EXTRAHEAD)+1);

  (*dest)->data = (unsigned short *)calloc(npix, sizeof(unsigned short));
  (*dest)->datastore = 2*npix;

  (*dest)->nhead = src->nhead;
  (*dest)->nx = src->nx;
  (*dest)->ny = src->ny;
  (*dest)->bitpix = src->bitpix;
  (*dest)->bias = src->bias;
  (*dest)->sky = src->sky;
  (*dest)->scrunch = src->scrunch;
}

/*
 * scrunch a fits file down by a factor
 */
scrunchit(src, dest, sx, sy, nx, ny, nborder, factor)
     struct ushortfits *src, **dest;	/* image data structure */
     int sx, sy, nx, ny;		/* location to be extracted */
     int nborder;			/* boundary to be ignored */
     int factor;			/* scrunch factor */
{
  int i, j, k, l, sum, nsx, idx, npix;

  npix = (nx*ny+(factor*factor)/2) / (factor*factor) + 1;
  copyframe(src, dest, npix);
  (*dest)->nx = nx / factor;
  (*dest)->ny = ny / factor;
  (*dest)->scrunch = src->scrunch * factor;

  src->ix = nx;
  src->iy = ny;
  src->sx = sx;
  src->sy = sy;

/* Copy the data, scrunched down */
  nsx = src->nx;
  for(l=0; l<ny/factor; l++) {
    for(k=0; k<nx/factor; k++) {
      sum = npix = 0;
      for(j=(l==0?nborder:0); j<(l==ny/factor-1?factor-nborder:factor); j++) {
	idx = sx + (j+sy)*nsx + factor*(k + l*nsx);
	for(i=(k==0?nborder:0); i<(k==nx/factor-1?factor-nborder:factor); i++) {
	  if((src->data)[idx] > 0) {
	    npix++;
	    sum += (src->data)[idx];
	  }
	  idx++;
	}
      }
      ((*dest)->data)[k+l*nx/factor] = (sum + npix/2) / MAX(1,npix);
    }
  }

/* Sanitize the header's NAXIS1, NAXIS2, and BITPIX */
  chfitshead(&i, (*dest)->header, "NAXIS1  ", NULL, nx/factor, 0.0);
  chfitshead(&i, (*dest)->header, "NAXIS2  ", NULL, ny/factor, 0.0);

  return(0);
}

/* Allocate memory for a new data frame */
newframe(frame)
struct ushortfits **frame;	/* image data structure */
{
  if(*frame == NULL) {
    *frame = (struct ushortfits *)malloc(sizeof(struct ushortfits));
    (*frame)->fname = NULL;
    (*frame)->headstore = (*frame)->datastore = 0;
  }
}

/* Wipe out memory allocated for an existing data frame */
eraseframe(frame)
struct ushortfits **frame;	/* image data structure */
{
  if((*frame)->fname != NULL) free((*frame)->fname);
  if( (*frame)->headstore > 0) free((*frame)->header);
  if((*frame)->datastore > 0) free((*frame)->data);
  (*frame)->headstore = (*frame)->datastore = 0;
  (*frame)->fname = NULL;
}








#define ERR_CANT_OPEN_FILE		1
#define ERR_CANT_READ_HEADER		2
#define ERR_CANT_READ_DATA		3
#define ERR_CANT_WRITE_HEADER		4
#define ERR_CANT_WRITE_DATA		5
#define ERR_INSUFFICIENT_HEADER		6
#define ERR_NO_NAXIS			7
#define ERR_NO_END			8
#define ERR_INCONSISTENT_NAXIS		9
#define ERR_BUFFER_OVERFLOW	       10
#define ERR_NOT_A_FITS_FILE	       11
#define ERR_BAD_BITPIX		       12

/*
 * rfits16() reads FITS image from disk, returning ushort data
 */
rfits16(frame, file)
char *file;			/* file name */
struct ushortfits **frame;	/* image data structure */
{
  int fd, mode=0, err, nbyte, nhead, npix, nx, ny;
  double scale, zero;	/* BSCALE and BZERO parameters from header */
  int bitpix;		/* BITPIX for disk and resultant array */
  int mefits, naxis, dims[10], i;
  int bp_dest, fitswab;

/* First see if we need to allocate memory for this frame */
  if(*frame == NULL) newframe(frame);

/* Save the file name */
  if((*frame)->fname != NULL) free((*frame)->fname);
  (*frame)->fname = (char *)malloc(strlen(file)+1);
  strcpy((*frame)->fname, file);

  /* See if this is a fits file */
  testfitslen_(file, &mefits, &bitpix, &naxis, dims, &nhead, strlen(file));
  if(mefits != 1) {
    fprintf(stderr, "File '%s' does not appear to be a FITS file!\n", file);
    return(ERR_NOT_A_FITS_FILE);
  }

/* Create an array for the FITS header if necessary */
  if(nhead*80 > (*frame)->headstore) {
    if((*frame)->headstore > 0) free((*frame)->header);
    (*frame)->header = (char *)malloc(80*(nhead+EXTRAHEAD)+1);
    (*frame)->headstore = 80*(nhead+EXTRAHEAD) + 1;
  }
  (*frame)->nhead = nhead;
  if( (err=openc_(&fd, file, &mode, 80)) != 0) return(err);
  if( (err=rhead_(&fd, &nhead, (*frame)->header, 80)) != 0) {
    close(fd); 
    return(err);
  }

  if( (err=parsehead_(&nhead, (*frame)->header, 
		      &bitpix, &nx,&ny, &scale, &zero,80)) != 0) {
    close(fd);   return(err);  }

  (*frame)->nx = nx;
  (*frame)->ny = ny;
  (*frame)->bitpix = -16;
  (*frame)->bias = (*frame)->sky = 0;
  (*frame)->scrunch = 1;
  chfitshead(&i, (*frame)->header, "BITPIX  ", NULL, -16, 0.0);
  chfitshead(&i, (*frame)->header, "BSCALE  ", "FLOAT", 0, 1.0);
  chfitshead(&i, (*frame)->header, "BZERO   ", "FLOAT", 0, 0.0);

  npix = nx * ny;
  nbyte = 2*((ABS(bitpix)*npix+15)/16);

/* Create an array for the data if necessary */
  if(nbyte > (*frame)->datastore) {
    if((*frame)->datastore > 0) free((*frame)->data);
    (*frame)->data = (unsigned short *)calloc(nbyte/2+1, 
					      sizeof(unsigned short));
    (*frame)->datastore = nbyte+2;
  }

  if( (err=read(fd, (*frame)->data, nbyte)) != nbyte) 
    {close(fd); return(ERR_CANT_READ_DATA);}
  close(fd);

  fitswab = 1;
  bp_dest = -16;
  fitsconvert_(&fitswab, &npix, &scale, &zero, 
	       &bitpix, (*frame)->data, &bp_dest, (*frame)->data);
  return(0);
}

wfits16(frame, file)
char *file;			/* file name */
struct ushortfits **frame;	/* image data structure */
{
  int fd, mode=1, err, nx, ny, bitpix;
  int nwrite, nbyte, imagebyte, nwrit;
  int nhead, nconv, fitswab, bp_src;
  double scale, zero;

  if( (err=openc_(&fd, file, &mode, strlen(file)+1)) != 0) return(err);

  nhead = (*frame)->nhead;

  if( (err=parsehead_(&nhead,(*frame)->header,
		      &bitpix,&nx,&ny,&scale,&zero,80)) != 0) {
    close(fd);    return(err);  }

  /*
  fprintf(stderr,"bitpix,nx,ny,scale,zero = %d %d %d %f %f\n", 
	  bitpix,nx,ny,scale,zero);
  */

  if( (err=whead_(&fd, &nhead, (*frame)->header, nhead*80+1)) != 0) {
    close(fd);  return(err); }

  imagebyte = 2*((ABS(bitpix)*nx*ny+15)/16);
  nbyte = imagebyte;		       	/* Bytes to write */

/* Set scale and zero for output conversion routines */
  if(scale != 1.0) scale = 1/scale;
  if(zero != 0.0) zero = -zero;
  fitswab = 2;
  bp_src = -16;

  nwrit = 0;
  while(nbyte > 0) {
    nwrite = MIN(MAXBUF, nbyte);

/* Convert from machine specific floating point to 16 bit format */
    nconv = nwrite/2;

    fitsconvert_(&fitswab, &nconv, &scale, &zero, 
	       &bp_src, ((*frame)->data)+nwrit/2, &bitpix, buf);

    if( (err=write(fd, buf, nwrite)) != nwrite) {
      close(fd); return(ERR_CANT_WRITE_DATA); }
    nbyte -= nwrite;
    nwrit += nwrite;
  }

  nwrite = NFITS - (imagebyte%NFITS);
  if(nwrite < NFITS) {
    bzero(buf,nwrite);
    if( (err=write(fd, buf, nwrite)) != nwrite) {
      close(fd); return(ERR_CANT_WRITE_DATA); }
  }

  close(fd);
  return(0);
}
