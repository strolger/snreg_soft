/* Program to remap one image into another */
/*    05 Oct 2000 (JT): Removed 0.5 pixel bug when sub > 1 */
/*    11 Jul 2000 (JT): added -bloat option */
/*    26 Jun 2000 (JT): added -region option */
/* John Tonry - 18 June 1999 */
/*
 * Syntax: remap [options]
 *		-src fname	input image to be warped
 *		-im fname	output file name
 *   -or-     -dest fname     output file name
 *		-map fname	parameter file with map coefficients
 *		-nx N		Size of output image in x direction
 *		-ny N		Size of output image in y direction
 *		-sx N		Start of output image in x direction
 *		-sy N		Start of output image in y direction
 *		-cx N		Center of output image in x direction
 *		-cy N		Center of output image in y direction
 *   -or-       -region [x1:x2,y1:y2] specify output image location
 *              -bloat x        Expand output image by factor x
 *              -bitpix N       Set output file format, N = 16/-16/-32
 *              -jacobian       Apply a Jacobian (conserve surface brightness)
 *              -poison x       Set poison data value to x (default 0.0)
 *		-sub N		subsample illumination calculation by N
 *		-halfpixel	if map is based on integer centers for pixels
13,14d9
 */

#include <stdio.h>
#include <math.h>

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))

#define MAXSIZE 1024
#define MAXPAR 10
/* #define CASES_FOR_SPEED	/* Usually slightly faster rebinning */
#define POISON 		/* Enable poison pixels */

main(argc,argv)
int argc;
char **argv;
{
  float *source, *warp;
  float *xmap, *ymap, *jacobian=NULL;
  char *header, *headin;
  float x, y, poison=0.0;

  int npix, err, nsx, nsy, nhead, jacob;
  int i, j, n;
  int mx0=0, my0=0, mx1=-1, my1=-1, cx=-1, cy=-1;
  double atof(), scale, zero;

  double xpar[MAXPAR], ypar[MAXPAR], xm, xp, ym, yp, xrng, yrng;
  char *fsrc=NULL, *fwarp=NULL, *fmap=NULL;
  double dudx, dudy, dvdx, dvdy;
  double bloat=1.0;
  FILE *input;
  int nx=MAXSIZE, ny=MAXSIZE, bitpix=0;
  int subsample=1, halfpixel=0;
  char line[256], *success;

  FILE *fin;

/* Parse the arguments */
  for(i=1; i<argc; i++) {

    if(strncmp(argv[i], "-nx", 3) == 0) {
      nx = atoi(argv[++i]);

    } else if(strncmp(argv[i], "-ny", 3) == 0) {
      ny = atoi(argv[++i]);

    } else if(strncmp(argv[i], "-sx", 3) == 0) {
      mx0 = atoi(argv[++i]);

    } else if(strncmp(argv[i], "-sy", 3) == 0) {
      my0 = atoi(argv[++i]);

    } else if(strncmp(argv[i], "-cx", 3) == 0) {
      cx = atoi(argv[++i]);

    } else if(strncmp(argv[i], "-cy", 3) == 0) {
      cy = atoi(argv[++i]);

    } else if(strncmp(argv[i], "-region", 7) == 0) {
      if( sscanf(argv[++i], "[%d:%d,%d:%d]", &mx0, &mx1, &my0, &my1) != 4) {
        syntax(*argv);
        exit(0);
      }

    } else if(strncmp(argv[i], "-sub", 4) == 0) {
      subsample = atoi(argv[++i]);

    } else if(strncmp(argv[i], "-bitpix", 7) == 0) {
      bitpix = atoi(argv[++i]);

    } else if(strncmp(argv[i], "-src", 4) == 0) {
      fsrc = argv[++i];

    } else if(strncmp(argv[i], "-map", 4) == 0) {
      fmap = argv[++i];

    } else if(strncmp(argv[i], "-jacobian", 9) == 0) {
      jacob = 1;

    } else if(strncmp(argv[i], "-poison", 7) == 0) {
      poison = atof(argv[++i]);

    } else if(strncmp(argv[i], "-bloat", 6) == 0) {
      bloat = atof(argv[++i]);

    } else if(strncmp(argv[i], "-half", 5) == 0) {
      halfpixel = -1;		/* I.e. center of LL pixel at (0.0,0.0) */

    } else if(strncmp(argv[i], "+half", 5) == 0) {
      halfpixel = +1;		/* I.e. center of LL pixel at (1.0,1.0) */

    } else if(strncmp(argv[i], "-im", 3) == 0) {
      fwarp = argv[++i];

    } else if(strncmp(argv[i], "-dest", 5) == 0) {
      fwarp = argv[++i];
    } else {
      syntax(*argv);
      exit(0);
    }   
  }
  if (fsrc == NULL || fwarp == NULL) {
    syntax(*argv);
    exit(0);
  }

/* Sanitize conflicting size arguments */

/* -region has priority */
  if(mx1 >= 0 && my1 >= 0) {
    nx = mx1 - mx0 + 1;
    ny = my1 - my0 + 1;

/* -cx -cy has next priority (nx, ny assumed) */
  } else if(cx >= 0 && cy >= 0) {
    mx0 = cx - nx/2;
    my0 = cy - ny/2;
  }

/* Allocate some memory for the destination image */
  npix = nx * ny;
  warp = (float *) calloc(npix, sizeof(float));

/* Read in source image */
  rfitsreal(&headin, &nsx, &nsy, &source, fsrc);
  printf("Read %s, size = (%d, %d)...\n", fsrc, nsx, nsy);

/* Read in coordinate mapping parameters */
  if(fmap == NULL) {
    xm = 0;
    ym = 0;
    xp = nsx;
    yp = nsy;
    for(i=0; i<MAXPAR; i++) {
      xpar[i] = ypar[i] = 0.0;
    }
    xpar[1] = 1.0;
    ypar[2] = 1.0;
  } else {
    if((input = fopen(fmap,"r")) == NULL) {
      fprintf(stderr,"Cannot open map file %s\n", fmap);
      exit(1);
    }
    if( (success = fgets(line, 256, input)) == line)
      sscanf(line, "%lf %lf %lf %lf", &xm, &xp, &ym, &yp);
    for(i=0; i<MAXPAR; i++) {
      if( (success = fgets(line, 256, input)) == line)
	sscanf(line, "%lf %lf", &xpar[i], &ypar[i]);
    }
    fclose(input);
  }

/* Allocate some memory for the maps */
  xmap = (float *) calloc(nsx*nsy, sizeof(float));
  ymap = (float *) calloc(nsx*nsy, sizeof(float));
  if(jacob) jacobian = (float *) calloc(nsx*nsy, sizeof(float));

/* Calculate mapping of the source image to destination coordinates */
  xrng = (xp-xm) / 2;
  yrng = (yp-ym) / 2;
  for(j=0; j<nsy; j++) {
    for(i=0; i<nsx; i++) {
      if(halfpixel == 0) {
	x = ((i    )-0.5*(xp+xm)) / xrng;
	y = ((j    )-0.5*(yp+ym)) / yrng;
      } else if(halfpixel == -1) {
	x = ((i-0.5)-0.5*(xp+xm)) / xrng;
	y = ((j-0.5)-0.5*(yp+ym)) / yrng;
      } else if(halfpixel == +1) {
	x = ((i+0.5)-0.5*(xp+xm)) / xrng;
	y = ((j+0.5)-0.5*(yp+ym)) / yrng;
      }
      xmap[i+j*nsx] =   xpar[0] +   xpar[1]*x + xpar[2]*y +
        xpar[3]*x*x +  xpar[4]*x*y + xpar[5]*y*y +
        xpar[6]*x*x*x +  xpar[7]*x*x*y + xpar[8]*x*y*y + xpar[9]*y*y*y;
      ymap[i+j*nsx] =   ypar[0] +   ypar[1]*x + ypar[2]*y +
        ypar[3]*x*x +  ypar[4]*x*y + ypar[5]*y*y +
        ypar[6]*x*x*x +  ypar[7]*x*x*y + ypar[8]*x*y*y + ypar[9]*y*y*y;
      xmap[i+j*nsx] *= bloat;
      ymap[i+j*nsx] *= bloat;
      if(jacob) {
	dudx =   xpar[1] + 2*xpar[3]*x + xpar[4]*y + 3*xpar[6]*x*x +
	  2*xpar[7]*x*y + xpar[8]*y*y;
	dudy = xpar[2] + xpar[4]*x + 2*xpar[5]*y + xpar[7]*x*x + 
	  2*xpar[8]*x*y + 3*xpar[9]*y*y;
	dvdx =   ypar[1] + 2*ypar[3]*x + ypar[4]*y + 3*ypar[6]*x*x +
	  2*ypar[7]*x*y + ypar[8]*y*y;
	dvdy = ypar[2] + ypar[4]*x + 2*ypar[5]*y + ypar[7]*x*x + 
	  2*ypar[8]*x*y + 3*ypar[9]*y*y;
	jacobian[i+j*nsx] = (dudx*dvdy-dudy*dvdx) / (xrng*yrng);
	if(jacobian[i+j*nsx] < 0.0) jacobian[i+j*nsx] *= -1;
	if(jacobian[i+j*nsx] == 0.0) {
	  fprintf(stderr,"Jacobian is zero at %d,%d; setting to 1\n",i,j);
	  jacobian[i+j*nsx] = 1.0;
	}
        jacobian[i+j*nsx] *= bloat*bloat;
	if(halfpixel == -1) {
	  xmap[i+j*nsx] += 0.5;
	  ymap[i+j*nsx] += 0.5;
	} else if(halfpixel == +1) {
	  xmap[i+j*nsx] -= 0.5;
	  ymap[i+j*nsx] -= 0.5;
	}
      }
    }
  }
  printf("Map arrays created...\n");

/* Rebin uniform illumination of lens plane to observer's plane */
  printf("Rebinning...\n");
  rebin(nsx,nsy,source,mx0,my0,nx,ny,warp,subsample,xmap,ymap,jacobian,poison);

  printf("Writing output file %s\n", fwarp);
/* Create a FITS header to write image */
  /*
  newfitshead(&header, -32, nx, ny, "Warped image");
  */
/* Alter the old FITS header for the warped image */
  chfitshead(&n, headin, "NAXIS1  ", NULL, nx, 0.0);
  chfitshead(&n, headin, "NAXIS2  ", NULL, ny, 0.0);
  if(bitpix != 0) {
    chfitshead(&n, headin, "BITPIX  ", NULL, bitpix, 0.0);
  } else {
    if(ifitshead(headin, "BITPIX  ", &bitpix) != 0) bitpix = -32;
  }
  chfitshead(&n, headin, "CNPIX1  ", NULL, 0, 0.0);
  chfitshead(&n, headin, "CNPIX2  ", NULL, 0, 0.0);

  /*
  for(i=0; i<16; i++) {
    fprintf(stderr,"%.30s\n", headin+80*i);
  }
  */
  if(bitpix == -32) {
    chfitshead(&n, headin, "BSCALE  ", "FLOAT", 0, 1.0);
    chfitshead(&n, headin, "BZERO   ", "FLOAT", 0, 0.0);
    wfitsreal(headin, warp, fwarp);
  } else if( ABS(bitpix) == 16 ) {
    wfitsshort(headin, warp, fwarp);
  }
}

/* Rebin image src to image dest using maps xmap and ymap */
rebin(mx,my,src,fx,fy,nx,ny,dest,subsample,xmap,ymap,jacobian,poison)
float *src, *dest;	/* source and destination images */
int mx, my;		/* size of source image */
int fx, fy, nx, ny;	/* start and size of destination image */
int subsample;		/* subsample factor (larger = better and slower) */
float *xmap;		/* Pixel [i,j] of source goes to x = xmap[i,j] */
float *ymap;		/* Pixel [i,j] of source goes to y = ymap[i,j] */
float *jacobian;	/* Jacobian at pixel [i,j] (or NULL if not used) */
float poison;		/* Nuke destination pixels which have this source */
{
  float rtmp, sx0, sy0, sxm, sym, sxp, syp, dp;
  float x0, x1, y0, y1, area, dy0, dyp, dym, dx0, dxp, dxm;
  int jm, jp, im, ip;
  int i, j, ii, jj, rotated;
  int ix0, ix1, iy0, iy1;
  register float *dptr, flux;
  register int nl, l, m, dstride;

#ifdef POISON
/* If poison = 0 have to initialize to something other than zero! */
  if(poison == 0.0) {
    for(i=0; i<nx*ny; i++) dest[i] = 0.125;
  } else {
    for(i=0; i<nx*ny; i++) dest[i] = 0.0;
  }
#else
  for(i=0; i<nx*ny; i++) dest[i] = 0.0;
#endif

  /*   printf("Rotated = %d...\n", rotated); */

  dp = 1.0 / subsample;
  for(jj=subsample; jj<subsample*(my-2); jj++) {
    jm = (jj-1) / subsample;
    j  = (jj  ) / subsample;
    jp = (jj+1) / subsample;
    sym = (((jj-1)%subsample) + 0.5) * dp;
    sy0 = (((jj  )%subsample) + 0.5) * dp;
    syp = (((jj+1)%subsample) + 0.5) * dp;
    for(ii=subsample; ii<subsample*(mx-2); ii++) {
      im = (ii-1) / subsample;
      i  = (ii  ) / subsample;
      ip = (ii+1) / subsample;
      sxm = (((ii-1)%subsample) + 0.5) * dp;
      sx0 = (((ii  )%subsample) + 0.5) * dp;
      sxp = (((ii+1)%subsample) + 0.5) * dp;

/* Bilinear interpolation of x, x-1, and x+1 */
      dxm =
	xmap[im   +  jm*mx]    * (1-sxm) * (1-sym) +
	xmap[im+1 +  jm*mx]    *    sxm  * (1-sym) +
	xmap[im+1 + (jm+1)*mx] *    sxm  *    sym  +
	xmap[im   + (jm+1)*mx] * (1-sxm) *    sym;

      dx0 =
	xmap[i    +  j*mx]     * (1-sx0) * (1-sy0) +
	xmap[i+1  +  j*mx]     *    sx0  * (1-sy0) +
	xmap[i+1  + (j+1)*mx]  *    sx0  *    sy0  +
	xmap[i    + (j+1)*mx]  * (1-sx0) *    sy0;

      dxp =
	xmap[ip   +  jp*mx]    * (1-sxp) * (1-syp) +
	xmap[ip+1 +  jp*mx]    *    sxp  * (1-syp) +
	xmap[ip+1 + (jp+1)*mx] *    sxp  *    syp  +
	xmap[ip   + (jp+1)*mx] * (1-sxp) *    syp;

/* Bilinear interpolation of y, y-1, and y+1 */
      dym =
	ymap[im   +  jm*mx]    * (1-sxm) * (1-sym) +
	ymap[im+1 +  jm*mx]    *    sxm  * (1-sym) +
	ymap[im+1 + (jm+1)*mx] *    sxm  *    sym  +
	ymap[im   + (jm+1)*mx] * (1-sxm) *    sym;

      dy0 =
	ymap[i    +  j*mx]     * (1-sx0) * (1-sy0) +
	ymap[i+1  +  j*mx]     *    sx0  * (1-sy0) +
	ymap[i+1  + (j+1)*mx]  *    sx0  *    sy0  +
	ymap[i    + (j+1)*mx]  * (1-sx0) *    sy0;

      dyp =
	ymap[ip   +  jp*mx]    * (1-sxp) * (1-syp) +
	ymap[ip+1 +  jp*mx]    *    sxp  * (1-syp) +
	ymap[ip+1 + (jp+1)*mx] *    sxp  *    syp  +
	ymap[ip   + (jp+1)*mx] * (1-sxp) *    syp;

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

      if(x1 < fx || y1 < fy) continue;

      if(x0 < fx) ix0 = fx;
      else ix0 = (int)x0;
      if(ix0 >= nx+fx) continue;

      ix1 = (int)x1;
      if(ix1 > nx+fx-1) ix1 = nx+fx-1;

      if(y0 < fy) iy0 = fy;
      else iy0 = (int)y0;
      if(iy0 >= ny+fy) continue;

      iy1 = (int)y1;
      if(iy1 > ny+fy-1) iy1 = ny+fy-1;

/* How much flux gets distributed from the source image? */
      if(src != NULL) {
	flux = src[i+mx*j] / (subsample*subsample*(x1-x0)*(y1-y0));
#ifdef POISON
	if(src[i+mx*j] == poison) flux = poison;
#endif
      }	else {
	flux = 1.0 / (subsample*subsample*(x1-x0)*(y1-y0));
      }
/* Multiply in the Jacobian (if non-NULL) */
      if(jacobian != NULL) flux *= jacobian[i+mx*j];

/* This is what we are doing... But break it into cases for speed... */
#ifndef CASES_FOR_SPEED
      dptr = &dest[ix0-fx+(iy0-fy)*nx];
      for(m=iy0; m<=iy1; m++) {
	nl = ix1 - ix0 + 1;
	for(l=ix0; l<=ix1; l++) {
	  area = (MIN(x1,l+1) - MAX(x0,l)) * (MIN(y1,m+1) - MAX(y0,m));
#ifdef POISON
	  if(flux == poison || *dptr == poison) *dptr++ = poison;
	  else *dptr++ += area * flux;
#else
	  *dptr++ += area * flux;
#endif
	}
	dptr += nx - nl;
      }

#else
/* But do it faster */
      if(ix1 == ix0) {
	if(iy1 == iy0) {
/* Entirely within a pixel */
	  dest[ix0-fx+(iy0-fy)*nx] += flux * (x1-x0) * (y1-y0);
	} else {
/* Strip in y */
	  flux *= (x1-x0);
	  dptr = &dest[ix0-fx+(iy0-fy)*nx];
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
          dptr = &dest[ix0-fx+(iy0-fy)*nx];
	  for(l=ix0; l<=ix1; l++) {
	    area = MIN(x1,l+1) - MAX(x0,l);
	    *dptr++ += area * flux;
	  }
	} else {
/* Corners */
          dest[ix0-fx+(iy0-fy)*nx] += flux * (ix0+1-x0) * (iy0+1-y0);
          dest[ix1-fx+(iy0-fy)*nx] += flux * (x1-ix1)   * (iy0+1-y0);
          dest[ix1-fx+(iy1-fy)*nx] += flux * (x1-ix1)   * (y1-iy1);
          dest[ix0-fx+(iy1-fy)*nx] += flux * (ix0+1-x0) * (y1-iy1);

/* Edges */
	  if(ix1-ix0 < 2 && iy1-iy0 < 2) continue;
	  for(l=iy0+1; l<=iy1-1; l++) {
            dest[ix0-fx+(l-fy)*nx] += flux * (ix0+1-x0);
            dest[ix1-fx+(l-fy)*nx] += flux * (x1-ix1);
	  }

	  for(l=ix0+1; l<=ix1-1; l++) {
            dest[l-fx+(iy0-fy)*nx] += flux * (iy0+1-y0);
            dest[l-fx+(iy1-fy)*nx] += flux * (y1-iy1);
	  }

/* Center */
          dptr = &dest[ix0-fx+1+(iy0-fy+1)*nx];
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
#ifdef POISON
  if(poison == 0.0) {
    for(i=0; i<nx*ny; i++) {
      if(dest[i] != poison) dest[i] -= 0.125;
    }
  }
#endif
}

syntax(s)
char *s;
{
  printf("Syntax: remap [options]\n");
  printf("       -src fname   Source image to be lensed\n");
  printf("       -im fname    Image of source seen through lens\n");
  printf(" -or-  -dest fname  Image of source seen through lens\n");
  printf("       -map fname   Parameter file with map coefficients\n"); 
  printf("       -nx N        Size of output image in x direction\n");
  printf("       -ny N        Size of output image in y direction\n");
  printf("       -sx N        Start of output image in x direction\n");
  printf("       -sy N        Start of output image in y direction\n");
  printf("       -cx N        Center of output image in x direction\n");
  printf("       -cy N        Center of output image in y direction\n");
  printf("       -region [x1:x2,y1:y2] specify output image location\n");
  printf("       -jacobian    Apply a Jacobian (conserve surface brightness, not flux)\n");
  printf("       -bloat x     Expand output image by factor x\n");
  printf("       -poison x    Set poison data value to x (default 0.0)\n");
  printf("       -bitpix N    Set output file format, N = 16/-16/-32\n");
  printf("       -sub N       subsample illumination calculation by N\n");
  printf("       -halfpixel   if map uses integer centers for pixels\n");
}
