/* getfluxbk                              Thu Aug 19 12:23:13 EST 1999 sj  *
 *                                                                         *
 *   gets flux and background from FITS image at specified positions       *
 *   requires median.c and rwfits.c                                        *
 *                                                                         *
 ***************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include "median.h"

#define INDEX(x,y,nx) ((y*nx)+x)
#define BUF 256

/* background annulus from aperture+start to aperture+start+width */


float BKG_ANNULUS_START=20.0;
float BKG_ANNULUS_WIDTH=8.0;

void usage(char *name);

int 
main(int argc, char *argv[]) {

  char *headin;
  int nhead, nx, ny, snx, sny;
  float *image, x, y;
  float junkf, inx, iny;
  float sub, *nzb;
  int i,j,ix,iy;
  int subsz;
  float ap2, bi2, bo2, r2,bad=0.;
  int suma, sumb;
  float maxa, flux, background, bkginap,mina;
  float bkg_inner_rad, bkg_outer_rad;  
  int warn = 0;
  float aperture = 5.0, saturation = 30000.0;
  FILE *in;
  char s[BUF+1],junk[BUF+1];
    
  float xoff=0.0, yoff=0.0;	/* JT 00102 */


  if (argc < 3) {
    usage(argv[0]);
    exit(1);
  }

  if (argc > 3) { /* parse optional command line arguments */
    for (i=3;i<argc;i++) {

      if ( strcasecmp(argv[i],"-warn") == 0 ) {
	warn = 1;
	continue;    /* goto next argument */
      }
      
      if ( strcasecmp(argv[i],"-aper") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-aper option requires argument: radius in pixels\n");
	  usage(argv[0]);
	  exit(1);
	}
	aperture = atof(argv[i]);
	if (aperture <= 0 || aperture > 5000 ) {
	  fprintf(stderr,"invalid aperture %f pixels\n",aperture);
	  usage(argv[0]);
	  exit(1);
	}
	continue;
      }

      if ( strcasecmp(argv[i],"-an") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-an option requires argument: inner radius in pixels\n");
	  usage(argv[0]);
	  exit(1);
	}
	BKG_ANNULUS_START=atof(argv[i]);
	if (BKG_ANNULUS_START <= 0 || BKG_ANNULUS_START > 5000 ) {
	  fprintf(stderr,"invalid annulus %f pixels\n",BKG_ANNULUS_START);
	  usage(argv[0]);
	  exit(1);
	}
	continue;
      }

      if ( strcasecmp(argv[i],"-xoff") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-an option requires argument: x offset added to coords\n");
	  usage(argv[0]);
	  exit(1);
	}
	xoff = atof(argv[i]);
	continue;
      }

      if ( strcasecmp(argv[i],"-yoff") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-an option requires argument: y offset added to coords\n");
	  usage(argv[0]);
	  exit(1);
	}
	yoff = atof(argv[i]);
	continue;
      }

      if ( strcasecmp(argv[i],"-width") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-an option requires argument: width of annulus in pixels\n");
	  usage(argv[0]);
	  exit(1);
	}
	BKG_ANNULUS_WIDTH=atof(argv[i]);
	if (BKG_ANNULUS_WIDTH <= 1 || BKG_ANNULUS_WIDTH > 5000 ) {
	  fprintf(stderr,"invalid annulus %f pixels\n",BKG_ANNULUS_WIDTH);
	  usage(argv[0]);
	  exit(1);
	}
	continue;
      }
      
      if ( strcasecmp(argv[i],"-bad") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-bad option requires argument: bad pixel value\n");
	  usage(argv[0]);
	  exit(1);
	}
	bad = atof(argv[i]);
	continue;
      }


      if ( strcasecmp(argv[i],"-sat") == 0 ) {
	if ( ++i >= argc ) {
	  fprintf(stderr,"-sat option requires argument: level in counts\n");
	  usage(argv[0]);
	  exit(1);
	}
	saturation = atof(argv[i]);
	if (saturation <= 0 ) {
	  fprintf(stderr,"invalid saturation level %f counts\n",saturation);
	  usage(argv[0]);
	  exit(1);
	}
	continue;
      }

      /* unrecognized command line option */
      fprintf(stderr,"unrecognized command line option: %s\n", argv[i]);
      usage(argv[0]);
      exit(1);
    }     

  }

  bkg_inner_rad =  BKG_ANNULUS_START;
  bkg_outer_rad = bkg_inner_rad + BKG_ANNULUS_WIDTH; 
  subsz = (int) (bkg_outer_rad + 2.0);   /* size of subraster */

  if (! (in = fopen(argv[2],"r"))) {
    fprintf(stderr,"Unable to open input file %s.\n",argv[2]);
    exit(1);
  }

  
  if (rfitsreal(&headin, &nx, &ny, &image, argv[1])) {
    fprintf(stderr,"Unable to open image %s\n",argv[1]);
    exit(1);
  }
  /*  printf("Read %s, size = (%d, %d)...\n", argv[1], nx, ny); */

  ap2 = aperture*aperture;
  bi2 = bkg_inner_rad*bkg_inner_rad;
  bo2 = bkg_outer_rad*bkg_outer_rad;
  
  /* Subraster size: this is the size of the grid we actually examine */
  snx = sny = 2*subsz + 1;

  /* create array to hold nonzero background values */
  /* it's bigger than it has to be, but that's okay */
  nzb = (float *) malloc( (size_t) (snx*sny*sizeof(float)));

  while (fgets(s,BUF,in)) {

    if (s[0] == '#') { /* ignore lines starting with # */
      continue;    /* get next input line */
    } 

    if ( sscanf(s,"%f%f%f%s",&junkf,&inx,&iny,junk) < 3 ) {
      if (warn) fprintf(stderr,"Unable to parse line in input file:\n%s",s);
      continue;    /* get next input line */
    }
  
    x = inx - 1.0 + xoff;   
    y = iny - 1.0 + yoff;

    if ( (x-subsz < 1)    || (y-subsz < 1) || 
	 (x+subsz > nx-1) || (y+subsz > ny-1) ) {
      if (warn) fprintf(stderr,"%9.3f %9.3f is too close to the edge!\n",inx,iny);
      printf("%9.3f %9.3f %12.3f %12.3f\n",inx,iny,0.0,0.0);
      continue;      /* get next input line */
    }

    suma = sumb = 0;
    maxa = flux = background = 0.0;
    mina=1e32;
    for (i=-subsz;i<=subsz;i++)
      for (j=-subsz;j<=subsz;j++) {
	ix = (int) (x + i);
	iy = (int) (y + j);
	sub = image[INDEX(ix,iy,nx)];
	r2 = ((float)ix-x)*((float)ix-x)+((float)iy-y)*((float)iy-y);
	
	if (r2 <= ap2) {    /* in the aperture */
	  suma++;
	  flux += sub;
	  if (sub > maxa) maxa = sub;
	  if (sub < mina) mina = sub;
	}
	
	
	if (bi2 <= r2 && r2 <= bo2) {   /* in the background */
	  nzb[sumb++] = sub;
	}
      }
  
    if (maxa > saturation || mina <= bad) {
      if (warn) fprintf(stderr,"%9.3f %9.3f has saturated/bad pixel in aperture!\n",inx,iny);
      printf("%9.3f %9.3f %12.3f %12.3f\n",inx,iny,0.0,0.0);
      continue;
    }
    
    median(nzb, sumb, &background);  /* median background per pixel */
    
    /* normalize to exactly pi*aperture^2 pixels */
    /*  this is necessary because the exact number of summed pixels */
    /*  varies based on the sub-pixel position of the center */
    flux *= (M_PI*ap2/(float)suma);

    bkginap = background * (M_PI*ap2);  /* how much background in aperture */
    
    printf("%9.3f %9.3f %12.3f %12.3f\n",inx,iny,flux-bkginap,background);
    
  }

  /* clean up */
  free(nzb);
  fclose(in);

  return 0;

}


void
usage(char *name) {
  
    fprintf(stderr,"Usage: %s image coordfile %s",
	    name, "[-aper pix] [-sat level] [-warn]\n\n");
    fprintf(stderr,"-aper pix:  aperture radius in pixels  (default 5)\n");
    fprintf(stderr,"-sat level: saturation level in counts (default 30000)\n");
    fprintf(stderr,"-bad level: bad pixel level in counts  (default 0)\n");
    fprintf(stderr,"-warn:      print warnings on stderr\n");
    fprintf(stderr,"-an:        background annulus radius  (default 20)\n");
    fprintf(stderr,"-width:     background annulus width   (default 8)\n");
    fprintf(stderr,"\n");
    fprintf(stderr,"coordfile format:  number x y [ignored junk...]\n");
    fprintf(stderr,"output on stdout:  x  y  flux-background  background/pixel\n");
    fprintf(stderr,"\n");
    return;
}












