/* Program to get scale from image1 to image2*/
/* Brian Schmidt - 18 August 1999 */
/*
 * Syntax: remap [options]
 *		-src fname	input image to be warped
 *		-im fname	output file name
 *		-map fname	parameter file with map coefficients
 *		-nx N		Size of output image in x direction
 *		-ny N		Size of output image in y direction
 *		-sub N		subsample illumination calculation by N
 */

#include <stdio.h>
#include <math.h>

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))

#define MAXSIZE 1024
#define MAXPAR 10
#define CASES_FOR_SPEED	/* Usually slightly faster rebinning */
#define POISON 		/* Enable poison pixels */

main(argc,argv)
int argc;
char **argv;
{
  float *source, *warp;
  float *xmap, *ymap;
  char *header, headin[1000*80];
  float x, y, xc, yc, xp, yp, poison=0.0;

  int npix, err, nsx, nsy, nhead;
  int i, j;
  double atof();

  double xpar[MAXPAR], ypar[MAXPAR], xmid, ymid;
  char *fimage=NULL, *ftemplate=NULL, *fstarfile=NULL;
  FILE *input;
  int nx=MAXSIZE, ny=MAXSIZE;
  int subsample=1;
  char line[256], *success;

  FILE *fin;

  fsrc = "n4434k.nog";
  fwarp = "warp.fits";

/* Parse the arguments */
  for(i=1; i<argc; i++) {
    if(strncmp(argv[i], "-ny", 3) == 0) {
      ny = atoi(argv[++i]);
    } else if(strncmp(argv[i], "-sub", 4) == 0) {
      subsample = atoi(argv[++i]);
    } else if(strncmp(argv[i], "-t", 2) == 0) {
      filet = argv[++i];
    } else if(strncmp(argv[i], "-i", 2) == 0) {
      filei = argv[++i];
    } else {
      syntax(*argv);
      exit(0);
    }   
  }

  /* Read in template */
  rfitsreal(headin, &nhead, &nsxt, &nsyt, &template, ftemplate);
  printf("Read %s, size = (%d, %d)...\n", fsrc, nsx, nsy);

  /*read in image*/
  rfitsreal(headin, &nhead, &nsxi, &nsyi, &image, fimage);
  printf("Read %s, size = (%d, %d)...\n", fsrc, nsx, nsy);
  

  if((input = fopen(starfile,"r")) == NULL) {
    fprintf(stderr,"Cannot open starfile file %s\n", starfile);
    exit(1);
  }
  while (feof(input)==0) {
    if( (success = fgets(line, 256, input)) == line) {
      if (sscanf(line, "%lf %lf %lf", &dummy,&xmid, &ymid)==3) {
	if (xmid > 10 && xmid <nsxt-10 && ymid >10 && ymid < nsyt-10) {
	  i++;
	  getbackratio(template,nsxt,nsyt,image,nsxi,nsyi,radius,&back[i],&ratio[i]);
	}
      }
    }
  }
  

  fclose(input);
  }

/* Allocate some memory for the maps */
  xmap = (float *) calloc(nsx*nsy, sizeof(float));
  ymap = (float *) calloc(nsx*nsy, sizeof(float));

/* Calculate mapping of the source image to destination coordinates */
  for(j=0; j<nsy; j++) {
    for(i=0; i<nsx; i++) {
      x = (i+0.5) / xmid - 2;
      y = (j+0.5) / ymid - 2;
      xmap[i+j*nsx] =   xpar[0] +   xpar[1]*x + xpar[2]*y +
        xpar[3]*x*x +  xpar[4]*x*y + xpar[5]*y*y +
        xpar[6]*x*x*x +  xpar[7]*x*x*y + xpar[8]*x*y*y + xpar[9]*y*y*y;
      ymap[i+j*nsx] =   ypar[0] +   ypar[1]*x + ypar[2]*y +
        ypar[3]*x*x +  ypar[4]*x*y + ypar[5]*y*y +
        ypar[6]*x*x*x +  ypar[7]*x*x*y + ypar[8]*x*y*y + ypar[9]*y*y*y;
    }
  }
  printf("Map arrays created...\n");

/* Rebin uniform illumination of lens plane to observer's plane */
  printf("Rebinning...\n");
  rebin(nsx,nsy,source,nx,ny,warp,subsample,xmap,ymap,poison);

  printf("Writing output file %s...\n", fwarp);
  /* Create a FITS header to write image */
  updatefitshead(&headin, -32, nx, ny);
  wfitsreal(headin, warp, fwarp);

}

/* Rebin image src to image dest using maps xmap and ymap */
rebin(mx,my,src,nx,ny,dest,subsample,xmap,ymap,poison)
float *src, *dest;	/* source and destination images */
int mx, my, nx, ny;	/* size of source and destination images */
int subsample;		/* subsample factor (larger = better and slower) */
float *xmap;		/* Pixel [i,j] of source goes to x = xmap[i,j] */
float *ymap;		/* Pixel [i,j] of source goes to y = ymap[i,j] */
float poison;		/* Nuke destination pixels which have this source */
{
  float rtmp, sx, sy, dp, *xvary, *yvary;
  float x0, x1, y0, y1, area, dy0, dyp, dym, dx0, dxp, dxm;
  int i, j, ii, jj, rotated;
  int ix0, ix1, iy0, iy1;
  register float *dptr, flux;
  register int nl, l, m, dstride;

  for(i=0; i<nx*ny; i++) dest[i] = 0;

  dxp = xmap[mx/2+1 +  (my/2)*mx] - xmap[mx/2-1 +  (my/2)*mx];
  dyp = ymap[mx/2+1 +  (my/2)*mx] - ymap[mx/2-1 +  (my/2)*mx];
  if(ABS(dyp) < ABS(dxp)) {
    rotated = 0;
    xvary = xmap;
    yvary = ymap;
  } else {
    rotated = 1;
    xvary = ymap;
    yvary = xmap;
  }

  /*   printf("Rotated = %d...\n", rotated); */

  dp = 1.0 / subsample;
  for(jj=subsample; jj<subsample*(my-2); jj++) {
    j = jj / subsample;
    sy = (jj%subsample) * dp;
    for(ii=subsample; ii<subsample*(mx-2); ii++) {
      i = ii / subsample;
      sx = (ii%subsample) * dp;

/* Bilinear interpolation of x, x-1, and x+1 */
      dx0 =
	xvary[i   +  j*mx]    * (1-sx) * (1-sy) +
	xvary[i+1 +  j*mx]    *    sx  * (1-sy) +
	xvary[i+1 + (j+1)*mx] *    sx  *    sy  +
	xvary[i   + (j+1)*mx] * (1-sx) *    sy;

      dxm =
	xvary[i-1 +  j*mx]    * (1-sx) * (1-sy) +
	xvary[i   +  j*mx]    *    sx  * (1-sy) +
	xvary[i   + (j+1)*mx] *    sx  *    sy  +
	xvary[i-1 + (j+1)*mx] * (1-sx) *    sy;

      dxp =
	xvary[i+1 +  j*mx]    * (1-sx) * (1-sy) +
	xvary[i+2 +  j*mx]    *    sx  * (1-sy) +
	xvary[i+2 + (j+1)*mx] *    sx  *    sy  +
	xvary[i+1 + (j+1)*mx] * (1-sx) *    sy;

/* Bilinear interpolation of y, y-1, and y+1 */
      dy0 =
	yvary[i   +  j*mx]    * (1-sx) * (1-sy) +
	yvary[i+1 +  j*mx]    *    sx  * (1-sy) +
	yvary[i+1 + (j+1)*mx] *    sx  *    sy  +
	yvary[i   + (j+1)*mx] * (1-sx) *    sy;

      dym =
	yvary[i   + (j-1)*mx] * (1-sx) * (1-sy) +
	yvary[i+1 + (j-1)*mx] *    sx  * (1-sy) +
	yvary[i+1 + j*mx]     *    sx  *    sy  +
	yvary[i   + j*mx]     * (1-sx) *    sy;

      dyp =
	yvary[i   + (j+1)*mx] * (1-sx) * (1-sy) +
	yvary[i+1 + (j+1)*mx] *    sx  * (1-sy) +
	yvary[i+1 + (j+2)*mx] *    sx  *    sy +
	yvary[i   + (j+2)*mx] * (1-sx) *    sy;

/* But are we rotated?  If so swap x's and y's */
      if(rotated) {
	{rtmp = dx0;   dx0 = dy0;   dy0 = rtmp;}
	{rtmp = dxm;   dxm = dym;   dym = rtmp;}
	{rtmp = dxp;   dxp = dyp;   dyp = rtmp;}
      }

/* Determine junctures between adjacent destination pixels */
      x0 = 0.5 * (dx0+dxm);
      x1 = 0.5 * (dx0+dxp);
      y0 = 0.5 * (dy0+dym);
      y1 = 0.5 * (dy0+dyp);

/* Fudge things if this is a fold in the map */
      if( (x0-dx0) * (x1-dx0) > 0) {x0 = 0.5 * (x0+x1);   x1 = 2*dx0 - x0;}
      if( (y0-dy0) * (y1-dy0) > 0) {y0 = 0.5 * (y0+y1);   y1 = 2*dy0 - y0;}

/* Make sure that x0 < x1 and y0 < y1 */
      if(x0 > x1) {rtmp = x0;   x0 = x1;   x1 = rtmp;}
      if(y0 > y1) {rtmp = y0;   y0 = y1;   y1 = rtmp;}

      if(x1 < 0 || y1 < 0) continue;

      if(x0 < 0) ix0 = 0;
      else ix0 = (int)x0;
      if(ix0 >= nx) continue;

      ix1 = (int)x1;
      if(ix1 > nx-1) ix1 = nx-1;

      if(y0 < 0) iy0 = 0;
      else iy0 = (int)y0;
      if(iy0 >= ny) continue;

      iy1 = (int)y1;
      if(iy1 > ny-1) iy1 = ny-1;

/* How much flux gets distributed from the source image? */
      if(src != NULL) {
	flux = src[i+mx*j] / (subsample*subsample*(x1-x0)*(y1-y0));
#ifdef POISON
	if(src[i+mx*j] == poison) flux = poison;
#endif
      }	else {
	flux = 1.0 / (subsample*subsample*(x1-x0)*(y1-y0));
      }

/* This is what we are doing... But break it into cases for speed... */
#ifndef CASES_FOR_SPEED
      dptr = &dest[ix0+iy0*nx];
      for(m=iy0; m<=iy1; m++) {
	nl = ix1 - ix0 + 1;
	for(l=ix0; l<=ix1; l++) {
	  area = (MIN(x1,l+1) - MAX(x0,l)) * (MIN(y1,m+1) - MAX(y0,m));
	  *dptr++ += area * flux;
#ifdef POISON
	  if(flux == poison) *(dptr-1) = poison;
#endif
	}
	dptr += nx - nl;
      }

#else
/* But do it faster */
      if(ix1 == ix0) {
	if(iy1 == iy0) {
/* Entirely within a pixel */
	  dest[ix0+iy0*nx] += flux * (x1-x0) * (y1-y0);
	} else {
/* Strip in y */
	  flux *= (x1-x0);
	  dptr = &dest[ix0+iy0*nx];
	  for(l=iy0; l<=iy1; l++) {
	    area = MIN(y1,l+1) - MAX(y0,l);
	    *dptr += area * flux;
	    dptr += nx;
	  }
	}

      } else {

	if(iy1 == iy0) {
/* Strip in x */
	  flux *= (y1-y0);
	  dptr = &dest[ix0+iy0*nx];
	  for(l=ix0; l<=ix1; l++) {
	    area = MIN(x1,l+1) - MAX(x0,l);
	    *dptr++ += area * flux;
	  }
	} else {
/* Corners */
	  dest[ix0+iy0*nx] += flux * (ix0+1-x0) * (iy0+1-y0);
	  dest[ix1+iy0*nx] += flux * (x1-ix1)   * (iy0+1-y0);
	  dest[ix1+iy1*nx] += flux * (x1-ix1)   * (y1-iy1);
	  dest[ix0+iy1*nx] += flux * (ix0+1-x0) * (y1-iy1);

/* Edges */
	  if(ix1-ix0 < 2 && iy1-iy0 < 2) continue;
	  for(l=iy0+1; l<=iy1-1; l++) {
	    dest[ix0+l*nx] += flux * (ix0+1-x0);
	    dest[ix1+l*nx] += flux * (x1-ix1);
	  }

	  for(l=ix0+1; l<=ix1-1; l++) {
	    dest[l+iy0*nx] += flux * (iy0+1-y0);
	    dest[l+iy1*nx] += flux * (y1-iy1);
	  }

/* Center */
	  dptr = &dest[ix0+1+(iy0+1)*nx];
	  nl = ix1 - ix0 - 1;
	  dstride = nx - nl;
	  m = iy1 - iy0 - 1;
	  while(m--) {
	    l = nl;
	    while(l--) *dptr++ += flux;
	    dptr += dstride;
	  }
	}
      }
#endif
    }
  }
}

#define NFITS (2880)

/*
 * newfitshead() will create a minimal FITS header
 */
newfitshead(header,bitpix,nx,ny,object)
char **header;		/* header */
int bitpix, nx, ny;	/* BITPIX, NAXIS1, NAXIS2 header entries */
char *object;		/* object name */
{
  int i;
  *header = (char *) malloc(NFITS);
  for(i=0; i<NFITS; i++) (*header)[i] = ' ';
  wfitem(0, *header, "SIMPLE  ", "T", 0);
  wfitem(1, *header, "BITPIX  ", NULL, bitpix);
  wfitem(2, *header, "NAXIS   ", NULL, 2);
  wfitem(3, *header, "NAXIS1  ", NULL, nx);
  wfitem(4, *header, "NAXIS2  ", NULL, ny);
  wfitem(5, *header, "OBJECT  ", object, 0);
  wfitem(6, *header, "EXPTIME ", NULL,1);
  wfitem(7, *header, "END     ", NULL, 1);

  return(0);
}

/*
 * updatefitshead() will update FITS header with nx,by, and bitpixinfo
 */
updatefitshead(header,bitpix,nx,ny)
char *header;		/* header */
int bitpix, nx, ny;	/* BITPIX, NAXIS1, NAXIS2 header entries */
{
  int i;
  /*  *header = (char *) malloc(NFITS);
  for(i=0; i<NFITS; i++) (*header)[i] = ' '; */
  wfitem(0, header, "SIMPLE  ", "T", 0);
  wfitem(1, header, "BITPIX  ", NULL, bitpix);
  wfitem(2, header, "NAXIS   ", NULL, 2);
  wfitem(3, header, "NAXIS1  ", NULL, nx);
  wfitem(4, header, "NAXIS2  ", NULL, ny);
  return(0);
}



/*
 * wfitem() writes a specific line to a FITS header
 */
wfitem(n, header, keyword, cvalue, ivalue)
char *header, *cvalue, *keyword;
int n, ivalue;
{
  int i;
  if(strncmp(keyword, "END     ", 8) == 0) {
    i = sprintf(&header[80*n], "%-8s", keyword);
  } else if(strncmp(keyword, "SIMPLE  ", 8) == 0) {
    i = sprintf(&header[80*n], "%-8s= %20s", keyword,cvalue);
  } else {
    if(cvalue == NULL) {
      i = sprintf(&header[80*n], "%-8s= %20d", keyword, ivalue);
    } else {
      i = sprintf(&header[80*n], "%-8s= '%s'", keyword, cvalue);
    }
  }
  header[80*n+i] = ' ';
  return(0);
}

/*
 * counthead() returns the number of lines in a FITS header
 */
counthead(header)
char *header;
{
  int i;
  for(i=0; i<80*1000; i+=80) {
    if(strncmp(&header[i], "END     ", 8) == 0) break;
  }
  return(i/80+1);
}

/*
 * wfitsreal() will write a real FITS image onto disk.
 */

#define ERR_CANT_READ_DATA		3
#define ERR_CANT_WRITE_DATA		5
#define ERR_INCONSISTENT_NAXIS		9
#define ERR_BUFFER_OVERFLOW	       10
#define MAXBUF (65536)

static char buf[MAXBUF];

wfitsreal(head,data,file)
char *head;		/* header */
char *file;		/* file name */
char *data;		/* image data */
{
  int fd, mode=1, err, nx, ny, bitpix;
  int nwrite, nbyte, imagebyte, nwrit;
  int nhead;
  double scale, zero;

  if( (err=openc_(&fd, file, &mode, strlen(file)+1)) != 0) return(err);

  nhead = counthead(head);

  if( (err=parsehead_(&nhead,head,&bitpix,&nx,&ny,&scale,&zero,80)) != 0) {
    close(fd);    return(err);  }

  if( (err=whead_(&fd, &nhead, head, strlen(head)+1)) != 0) {
    close(fd);  return(err); }

  imagebyte = 2*((ABS(bitpix)*nx*ny+15)/16);
  nbyte = imagebyte;		       	/* Bytes to write */

  nwrit = 0;
  while(nbyte > 0) {
    nwrite = MIN(MAXBUF,nbyte);
/* Convert from machine specific floating point to FITS IEEE format */
    nx = nwrite/4;
    fpieee_(&nx, &data[nwrit], buf);

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

rfitsreal(head,nhead,nx,ny,data,file)
char *head;		/* header */
char *file;		/* file name */
float **data;		/* image data */
int *nx, *ny;		/* NAXIS1 and NAXIS2 from header */
int *nhead;             /* number of header lines */
{
  int fd, mode=0, err, nbyte,npix;
  double scale, zero;	/* BSCALE and BZERO parameters from header */
  int bitpix;		/* BITPIX for disk and resultant array */

  if( (err=openc_(&fd, file, &mode, 80)) != 0) return(err);

  *nhead = counthead(head);
  if( (err=rhead_(&fd, nhead, head, 80)) != 0) {close(fd); return(err);}
  if( (err=parsehead_(nhead,head,&bitpix,nx,ny,&scale,&zero,80)) != 0) {
    close(fd);    return(err);  }

  npix = (*nx) * (*ny);
  nbyte = 2*((ABS(bitpix)*npix+15)/16);

  *data = (float *) calloc(npix+17, sizeof(float));

  if( (err=read(fd, *data, nbyte)) != nbyte) 
    {close(fd); return(ERR_CANT_READ_DATA);}
  close(fd);

  if(bitpix == -32) {
    ieeefp_(&npix, *data, *data);

  } else if(bitpix == 16) {
    shortfp_(&npix, &scale, &zero, *data, *data);

  } else if(bitpix == -16) {
    unsignedfp_(&npix, &scale, &zero, *data, *data);

  } else {
    fprintf(stderr,"I cannot deal with BITPIX = %d\n", bitpix);
    return(1);
  }

  return(0);
}

syntax(s)
char *s;
{
  printf("Syntax: remap [options]\n");
  printf("       -src fname      Source image to be lensed\n");
  printf("       -im fname       Image of source seen through lens\n");
  printf("       -map fname	 Parameter file with map coefficients\n"); 
  printf("       -nx N           Size of output image in x direction\n");
  printf("       -ny N           Size of output image in y direction\n");
  printf("       -sub N          subsample illumination calculation by N\n");
}





