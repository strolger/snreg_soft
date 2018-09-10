/* Program to create a superflat from a pile of FITS images */
/* NEED:
 *    Need to search among accumulators for peak
 *    "Curious, nothing written" needs mean and expanded sigma, some how
 *    
 */
/* John Tonry - 26 Oct 99, 29 May 00 */
/*
 * Syntax: superflat [options] source_files...
 *		-out fname	output file name
 *		-imx N		Size of input images in x direction
 *		-imy N		Size of input images in y direction
 *		-sx N		Starting point of image in x direction
 *		-sy N		Starting point of image in y direction
 *		-border N	Number of bad border pixels
 *		-bias N		Set bias value
 *		-nmem N		Number of frames which fit in memory
 *		-scrunch N	Use scrunch factor N
 *		-test		Write test output files
 *		-scramble	Scramble order of use of input files
 */

#include <stdio.h>
#include <math.h>
#include "superflat.h"

int TEST=0;

/* Smaller RMSCALE means more accuracy but can overflow if nmem is small */
/* RMSZERO is provided in case RMSCALE*rms comes out less than 0.5, oh dear */
#define RMSCALE	2.0	/* Size of accumulator bin = rms*RMSCALE+RMSZERO */
#define RMSZERO	1.0	/* Size of accumulator bin = rms*RMSCALE+RMSZERO */

struct ushortfits mask[MAXFILE];
struct ushortfits frame[MAXFILE];

main(argc,argv)
int argc;
char **argv;
{
  int i, j, n, arg;
  double atof();
  char *fin=NULL, *fout=NULL;
  int imx=0, imy=0;	/* Dimensions of actual image */
  int sx=0, sy=0;		/* Starting points of actual image */
  int nborder=0;		/* Number of bad border pixels */
  int scrunch=8;		/* Scrunch factor */
  int nmem=5;			/* Number of frames which fit in memory */
  int setbias=-1;		/* Force bias to this level */
  int median, mode, medc, sig, nfile, rms, sky, bias;
  char fname[80];
  double gain=0.0;		/* CCD gain */
  double atof();
  int err;
  int SCRAMBLE=0, nextarg[MAXFILE];

  struct ushortfits *fptr=NULL, *mptr;
  struct ushortfits maskflat;

/* Parse the arguments */
  for(arg=1; arg<argc; arg++) {

    if(strncmp(argv[arg], "-imx", 4) == 0) {
      imx = atoi(argv[++arg]);

    } else if(strncmp(argv[arg], "-imy", 4) == 0) {
      imy = atoi(argv[++arg]);

    } else if(strncmp(argv[arg], "-sx", 3) == 0) {
      sx = atoi(argv[++arg]);

    } else if(strncmp(argv[arg], "-sy", 3) == 0) {
      sy = atoi(argv[++arg]);

    } else if(strncmp(argv[arg], "-nmem", 5) == 0) {
      nmem = atoi(argv[++arg]);

    } else if(strncmp(argv[arg], "-border", 7) == 0) {
      nborder = atoi(argv[++arg]);

    } else if(strncmp(argv[arg], "-scrunch", 8) == 0) {
      scrunch = atoi(argv[++arg]);

    } else if(strncmp(argv[arg], "-bias", 5) == 0) {
      setbias = atoi(argv[++arg]);

    } else if(strncmp(argv[arg], "-out", 4) == 0) {
      fout = argv[++arg];

    } else if(strncmp(argv[arg], "-gain", 5) == 0) {
      gain = atof(argv[++arg]);

    } else if(strncmp(argv[arg], "-test", 5) == 0) {
      TEST = 1;

    } else if(strncmp(argv[arg], "-scramble", 9) == 0) {
      SCRAMBLE = 1;

    } else {
      if(argv[arg][0] == '-') {
	syntax(*argv);
	exit(0);
      }
      break;
    }   
  }

/* We have to have at least 3 frames to accumulate mode and 1 for reads */
  nmem = MAX(nmem, 4);

  nfile = argc - arg;

  if(nfile > MAXFILE) {
    fprintf(stderr,"Maximum number of images is %d; but %d requested\n",
	    MAXFILE, nfile);
    exit(1);
  }

  if(nfile <= nmem) {
    printf("Setting nmem to NFILE-1 = %d\n", nfile-1);
    nmem = nfile - 1;
  }

  if(SCRAMBLE == 0) {
    for(i=0; i<nfile; i++) nextarg[i] = arg + i;
  } else {
    n = pow(2.0, (int)(log((double)nfile)/log(2.0) - 0.001)) + 0.5;
    j = 0;
    nextarg[j++] = 0;
    while(n > 0) {
      for(i=1; i<nfile; i+=2) {
	if(i*n >= nfile) break;
	nextarg[j++] = i*n;
      }
      n /= 2;
    }
  }
  if(TEST) {
    printf("Scramble sequence:\n");
    for(i=0; i<nfile; i++) printf("%3d%c", nextarg[i], i==nfile-1?'\n':' ');
  }
/* 
 * Loop over all source files.  After nmem have been accumulated, pause
 * to scrunch them up, get a scrunched estimate of the flatfield, create
 * masks, and initialize the mode accumulation.
 * Subsequent images can be scrunched and masked directly, and contributions
 * added to the mode.
 */
  for(i=0; i<nfile; i++) {

    fin = argv[arg + nextarg[i]];

/* 
 * Zero out frame unless we've got the mode established
 * Once that happens, frame[0] is our accumulator for subsequent reads,
 * frame[1] is our estimate of the mode
 * sig becomes our estimate of the...
 * and frame[2:nmem-1] have their data areas junked to become mode counters
 */
    if(i <= nmem-1) {
      fptr = &frame[i];
      frame[i].fname = NULL;
      frame[i].headstore = frame[i].datastore = 0;
    } else {
      fptr = &frame[0];
    }

/* Read in source image */
    err = rfits16(&fptr, fin);
    printf("Read %s, size = (%d,%d) bitpix = %d %s%d)\n", 
	   fin, fptr->nx, fptr->ny, fptr->bitpix, 
	   err?"(error = ":"(", err);

    if(i == 0 && (fptr->nx%scrunch) != 0) {
      fprintf(stderr,"Well, folks, there's currently a bug whereby NX has to be an even\n");
      fprintf(stderr,"multiple of SCRUNCH.  So it's probably best to quit right now...\n");
      fprintf(stderr,"Sorry 'bout that.  Offer JT inducements to fix it.\n");
      exit(1);
    }

/*
    for(j=0; j<20; j++) printf("%6d",fptr->data[j]);   printf("\n");
*/
    if(imx == 0) imx = fptr-> nx;
    if(imy == 0) imy = fptr-> ny;

/* Get gain from the header (if available) */
    if(gain == 0.0) {
      if(rfitshead(frame[0].header, "GAIN    ", &gain) != 0) gain = 1.0;
      printf("Gain = %.2f\n", gain);
      if(gain == 0.0) {
	fprintf(stderr,"Something's wrong, gain still is 0, setting to 1\n");
	gain = 1.0;
      }
    }

/*
  for(i=0; i<fptr->nhead; i++) printf("%.80s\n", fptr->header+80*i);
*/

/*
  wfits16(&fptr, "test.big");
*/

/* antique -- we probably want to count, not sort */
    /*
    fptr->bias = medsort(fptr->nx-1-4, sy, 4, imy, fptr);
    printf("bias = %d\n", fptr->bias);
    */

/* Count the median in the bias strip to get the bias level */
    if(setbias == -1) {
      bias = medcount(fptr->nx-1-4, sy, 4, imy, fptr, &sig);
      printf("Bias median = %d, sigma = %d\n", bias, sig);
    } else {
      bias = setbias;
    }

/* antique -- median is enough */
    /*
    medc = mclipcount(bias-3*sig, bias+3*sig);
    mode = modecount(medc, sig);
    printf("bias median = %d, sigma = %d, clip_med = %d, mode = %d\n", 
	   bias, sig, medc, mode);
    */

/*
    median = medcount(sx, sy, imx, imy, fptr, &sig);
    medc = mclipcount(median-3*sig, median+3*sig);
    mode = modecount(medc, sig);

    printf("data median = %d, sigma = %d, clip_med = %d, mode = %d\n", 
	   median, sig, medc, mode);
*/

/* Initialize the mask array */
    mptr = &mask[i];
    mask[i].fname = NULL;
    mask[i].headstore = mask[i].datastore = 0;

/* Fill the mask array with a scrunched down image */
    scrunchit(fptr, &mptr, sx,sy,imx,imy,nborder, scrunch);

/* Get the mode of the mask array as a sky estimate */
    median = medcount(nborder/scrunch, nborder/scrunch, 
		      mask[i].nx-2*(nborder/scrunch), 
		      mask[i].ny-2*(nborder/scrunch), 	 mptr, &sig); 
    medc = mclipcount(median-3*sig, median+3*sig);
    mode = modecount(medc, sig);
    mask[i].bias = bias;
    if(ABS(mode-median) < sig) {
      mask[i].sky = mode - bias;
    } else {
      printf("Funny mode from tilted sky?  Using median.\n");
      mask[i].sky = median - bias;
    }

    printf("Scrunched median = %d, sigma = %d, mode = %d, sky = mode-bias = %d\n", 
	   median, sig, mode, mask[i].sky);

    /*
    printf("%s: size = (%d,%d) bitpix = %d bias = %d sky = %d\n",
	   mask[i].fname, mask[i].nx, mask[i].ny, mask[i].bitpix, 
	   mask[i].bias, mask[i].sky);
    */

    /*
    strcpy(fname, mask[i].fname);
    strcat(fname, ".scr");
    wfits16(&mptr, fname);
    */


/* When enough have accumulated (N = nmem)
 *   create the mask flatfield
 *   create the first N masks
 *   get the mode and sig arrays
 *   assemble the first N contributions to the mode
 */
    if(i == nmem-1) {
      printf("Computing median flatfield\n");
      mptr = &maskflat;
      maskflat.fname = NULL;
      maskflat.headstore = maskflat.datastore = 0;

      medianflat(nmem, mask, &mptr);

      if(TEST) wfits16(&mptr, "test.flat");

      printf("flatfield written\n");
      printf("nmem = %d\n", nmem);

      sky = 0;
      for(j=0; j<nmem; j++) {
	printf("creating mask %d\n", j);
	createmask(&mask[j], &maskflat, gain, 5.0, 3.0);
	/*
	mptr = &mask[j];
	wfits16(&mptr, "test.mask");
	*/
	printf("mask %d created\n", j);
	sky += mask[j].sky;
      }

      sky = sky / nmem;
      rms = MAX(1, (int)(((float)MEANFLAT)/sky *
			 sqrt( 0.5 * gain * sky) / gain + 0.5));

      printf("Mean sky = %d, when scaled to %d rms = %d\n", sky, MEANFLAT, rms);

      createmode(nmem, mask, frame, rms);
      fptr = &frame[1];
      if(TEST) wfits16(&fptr, "test.mode");
      printf("mode %d created\n", j);


/* When i >= nmem, just keep up the mask creation and mode contributions */
    } else if(i >= nmem-1) {
      createmask(&mask[i], &maskflat, gain, 5.0, 3.0);
      printf("mask %d created\n", i);
      contribmode(nmem, &mask[i], frame, rms);
    }
  }

/* 
 * Now modify the mode frame according to the counts which have been
 * accumulated.
 */
  printf("Now refining the mode...\n");
  refinemode(nmem, frame, rms);
  fptr = &frame[1];
  if(fout != NULL) {
    wfits16(&fptr, fout);
    printf("Output file %s written\n", fout);
  } else {
    wfits16(&fptr, "super.flat");
    printf("Output file super.flat written\n");
  }

  if(TEST) {
    fptr = &frame[0];
    wfits16(&fptr, "test.resid");
    if(nmem > 4) {
      fptr = &frame[2];
      wfits16(&fptr, "test.minus");
      fptr = &frame[3];
      wfits16(&fptr, "test.ctr");
      fptr = &frame[4];
      wfits16(&fptr, "test.plus");
    }
    if(nmem > 5) {
      fptr = &frame[5];
      wfits16(&fptr, "test.sum");
    }
    if(nmem > 6) {
      fptr = &frame[6];
      wfits16(&fptr, "test.choice");
    }
  }
}

medianflat(nmem, mask, flat)
int nmem;
struct ushortfits *mask, **flat;
{
  float buf[MAXFILE];
  float invsky[MAXFILE];
  int i, j, k, nmed, scaled, nx, ny;

  nx = mask[0].nx;
  ny = mask[0].ny;
  copyframe(mask, flat, nx*ny+1);
  printf("frame copied %d %d\n", nx, ny);

  for(k=0; k<nmem; k++) {
    invsky[k] = 1.0 / (float)mask[k].sky;
    /*
    printf("%d, bias = %d, sky = %d, data = %d @ %d\n", 
	   k, mask[k].bias, mask[k].sky, (mask[k].data)[0], mask[k].data);
    */
  }

  for(j=0; j<ny; j++) {
    for(i=0; i<nx; i++) {
      for(k=0, nmed = 0, buf[0]=0.0; k<nmem; k++) {
	scaled = (mask[k].data)[i+j*nx];
	if(scaled > 0) {
	  buf[nmed++] = (scaled-mask[k].bias) * invsky[k];
	  /*
	  printf("%d %6d %9.4f\n",k, scaled, buf[nmed-1]);
	  */
	}
      }
      qsort4(nmed, buf);
      (*flat)->data[i+j*nx] = MEANFLAT * 0.5*(buf[nmed/2]+buf[(nmed-1)/2]);
      /*
      printf("%4d %4d %4d %8.4f %6d\n",i, j, nmed, 
	     0.5*(buf[nmed/2]+buf[(nmed-1)/2]), (*flat)->data[i+j*nx]);
      */
    }
  }
}

static unsigned char *aptr[2*MAXFILE];

createmode(nmem, mask, frame, rms)
int nmem, rms;
struct ushortfits *mask, *frame;
{
  float buf[MAXFILE], invsky[MAXFILE];
  
  int nsig[MAXFILE];
  int i, j, k, nmed, scaled, nx, ny, scrunch, median;
  int sx, sy, ix, iy, idx, imx;
  double rmscinv;
  int ncurious;

  nx = frame[0].nx;
  ny = frame[0].ny;
  ix = frame[0].ix;
  iy = frame[0].iy;
  sx = frame[0].sx;
  sy = frame[0].sy;
  scrunch = mask[0].scrunch;
  rmscinv = 1.0 / (RMSCALE*rms+RMSZERO);

/* Create pointers to the accumulation arrays */
/* The aptr are uchar pointers to the (now no longer needed) storage which
 * was allocated for frames 2, 3, ... nmem-1
 * Note the clever interleaving of the arrays...
 */
  for(i=0; i<2*(nmem-2); i++) {
    aptr[i] = (unsigned char *)(frame[2+i/2]).data;
    if(i & 1) aptr[i] += 1;
  }

  for(k=0; k<nmem; k++) {
    invsky[k] = 1.0 / (float)mask[k].sky;
    /*
    printf("%d, bias = %d, sky = %d, data = %d @ %d\n", 
	   k, mask[k].bias, mask[k].sky, (mask[k].data)[0], mask[k].data);
    */
  }

  for(j=0; j<iy; j++) {
    for(i=0; i<ix; i++) {
      idx = i+sx + (j+sy)*nx;
      imx = i/scrunch + (j/scrunch)*(ix/scrunch);
      for(k=0, nmed = 0; k<nmem; k++) {
	scaled = ((frame[k].data)[idx]-mask[k].bias) * (mask[k].data)[imx];
	if(scaled > 0) {
	  buf[nmed++] = MEANFLAT * scaled * invsky[k];
	  /*
	  printf("%d %6d %9.4f\n",k, scaled, buf[nmed-1]);
	  */
	}
      }
      qsort4(nmed, buf);
      median = 0.5*(buf[nmed/2]+buf[(nmed-1)/2]) + 0.5;
/* Note the bias offset of nmem-2 to make + for NINT and to prep for save */
      for(k=0; k<nmem; k++) {
	nsig[k] = rmscinv * (buf[k] - median) + nmem-2 + 0.5;
	(frame[k].data)[idx] = 0;
      }
/* Now destroy frame[1] with the median */
      (frame[1].data)[idx] = median;
/* Add in counts to the mode accumulators */
      median = 0;
      for(k=0; k<nmed; k++) {
	if(nsig[k] >= 0 && nsig[k] < 2*(nmem-2)) {
	  *(aptr[nsig[k]]+2*idx) += 1;
	  median++;
	}
      }
      if(median==0) {
	ncurious++;
	if((ncurious%100) == 0 && nmed>0) {
	  printf("Curious... nothing written at %d,%d out of %d\n", i,j, nmed);
	  for(k=0; k<nmed; k++) printf(" %8.1f", buf[k]); printf("\n");
	  for(k=0; k<nmed; k++) printf(" %8d", nsig[k]); printf("\n");
	}
      }	

/* For testing, save the entire sorted arrays */
/*
      for(k=0; k<nmem; k++) {
	(frame[k].data)[idx] = k < nmed ? MEANFLAT * buf[k] : 0;
      }
*/
      /*
      printf("%4d %4d %4d %8.4f %6d\n",i, j, nmed, 
	     0.5*(buf[nmed/2]+buf[(nmed-1)/2]), (*flat)->data[i+j*nx]);
      */
    }
  }
}

contribmode(nmem, mask, frame, rms)
int nmem, rms;
struct ushortfits *mask, *frame;
{
  unsigned char *ptr[MAXFILE];
  float invsky;
  int i, j, k, nmed, scaled, nx, ny, scrunch, median, bias;
  int sx, sy, ix, iy, idx, imx, nsig;
  double rmscinv;

  nx = frame[0].nx;
  ny = frame[0].ny;
  ix = frame[0].ix;
  iy = frame[0].iy;
  sx = frame[0].sx;
  sy = frame[0].sy;
  bias = mask->bias;
  scrunch = mask->scrunch;
  rmscinv = 1.0 / (RMSCALE*rms+RMSZERO);

  invsky = 1.0 / (float)mask->sky;
  /*
  printf("bias = %d, sky = %d, data = %d @ %d\n", 
	   mask->bias, mask->sky, (mask->data)[0], mask->data);
  */
/* Create pointers to the accumulation arrays */
  for(i=0; i<2*(nmem-2); i++) {
    ptr[i] = (unsigned char *)(frame[2+i/2]).data;
    if(i & 1) ptr[i]++;
  }

  for(j=0; j<iy; j++) {
    for(i=0; i<ix; i++) {
      idx = i+sx + (j+sy)*nx;
      imx = i/scrunch + (j/scrunch)*(ix/scrunch);
      scaled = MEANFLAT * invsky * 
	((frame->data)[idx] - bias) * (mask->data)[imx] + 0.5;
      median = (frame[1].data)[idx];

/* Note the bias offset of nmem-2 to make + for NINT and to ready for save */
      nsig = rmscinv * (scaled - median) + nmem-2 + 0.5;

/* Add in counts to the mode accumulators */
      if(nsig >= 0 && nsig < 2*(nmem-2)) *(aptr[nsig]+2*idx) += 1;
    }
  }
}

refinemode(nmem, frame, rms)
int rms;
struct ushortfits *frame;
{
  float invsky;
  int i, j, k, l, nmed, scaled, nx, ny, mode, sum, big;
  int sx, sy, ix, iy, idx, imx, m, c, p;
  unsigned char *mptr, *cptr, *pptr;
  double rmscale;

  nx = frame[0].nx;
  ny = frame[0].ny;
  ix = frame[0].ix;
  iy = frame[0].iy;
  sx = frame[0].sx;
  sy = frame[0].sy;
  rmscale = RMSCALE*rms + RMSZERO;

  /*
  mptr = aptr[nmem-3];
  cptr = aptr[nmem-2];
  pptr = aptr[nmem-1];
  */

  for(j=0; j<iy; j++) {
    for(i=0; i<ix; i++) {
      idx = i+sx + (j+sy)*nx;

/* Find the accumulator with the largest count */
      for(k=0, big=0; k<2*(nmem-2); k++) {
	if(*(aptr[k]+2*idx) > big) {
	  big = *(aptr[k]+2*idx);
	  l = k;
	}
      }

/* Get the counts */
      if(l == 0) {
	mode = -rmscale * (nmem-2);
	if(TEST) {
	  m = -1;
	  c = *(aptr[l]   + 2*idx);
	  p = *(aptr[l+1] + 2*idx);
	}
      } else if(l == 2*(nmem-2)-1) {
	mode = +rmscale * (nmem-1);
	if(TEST) {
	  m = *(aptr[l-1] + 2*idx);
	  c = *(aptr[l]   + 2*idx);
	  p = -1;
	}

      } else {
	m = *(aptr[l-1] + 2*idx);
	c = *(aptr[l]   + 2*idx);
	p = *(aptr[l+1] + 2*idx);

	if( (c >= m && c > p) || (c > m && c >= p) ) {
	  mode = ((double)(m-p) / (2*(m+p)-4*c)) * rmscale + 0.5;
	} else {
	  if(m > p) mode = -rmscale;
	  else      mode =  rmscale;
	}
	mode += rmscale * (l - (nmem-2));

      }

/* Add in counts to the mode accumulators */
      /*
      m = *(mptr + 2*idx);
      c = *(cptr + 2*idx);
      p = *(pptr + 2*idx);

      if( (c >= m && c > p) || (c > m && c >= p) ) {
	mode = ((double)(m-p) / (2*(m+p)-4*c)) * rmscale + 0.5;
      } else {
	if(m > p) mode = -rmscale;
	else      mode =  rmscale;
      }
      */

      (frame[1].data)[idx] += mode;
      if(TEST) {
/* For honks, save the mode correction */
	(frame[0].data)[idx] = mode + MEANFLAT;
/* For more honks, save the location of the mode */
	if(nmem > 6) {
	  (frame[6].data)[idx] = l;
	}
/* For yet more honks, save the total counts */
	if(nmem > 5) {
	  for(k=0, sum=0; k<2*(nmem-2); k++) sum += *(aptr[k]+2*idx);
	  (frame[5].data)[idx] = sum;
	}
/* For ultimate honks, save the counts from the m, c, and p accumulators */
	if(nmem > 4) {
	  (frame[2].data)[idx] = m;
	  (frame[3].data)[idx] = c;
	  (frame[4].data)[idx] = p;
	}
      }
    }
  }
}

createmask(mask, flat, gain, badsig, poorsig)
struct ushortfits *mask, *flat;
double gain, badsig, poorsig;
{
  int i, j, k, nx, ny, bias, sky;
  double nsig, rms;

  nx = mask->nx;
  ny = mask->ny;
  k = mask->scrunch*mask->scrunch;
  /*   fprintf(stderr,"createmask: k = %d, gain = %.3f", k, gain); */
  rms = sqrt(gain * mask->sky * k) / (k*gain);
  bias = mask->bias;
  sky = mask->sky;

  printf("rms = %.2f bias = %d sky = %d...", rms, bias, sky);

  for(j=0; j<ny; j++) {
    for(i=0; i<nx; i++) {
      if((flat->data)[i+j*nx] - sky > 0) {
	nsig = ( (((int)((mask->data)[i+j*nx]-bias)) * MEANFLAT) / 
		 (flat->data)[i+j*nx] - sky) / rms;
      } else {
	nsig = 0.0;
      }
      /*
      (mask->data)[i+j*nx] = 10*nsig + 10000;
      */
      if(nsig > badsig) {
	(mask->data)[i+j*nx] = 0;
      } else if(nsig > poorsig && i > 0 && i < nx-1 && j > 0 && j < ny-1) {
	(mask->data)[i+j*nx] = 2;
      } else {
	(mask->data)[i+j*nx] = 1;
      }
    }
  }

  printf("finding neighbors...", rms, bias, sky);
  for(j=1; j<ny-1; j++) {
    for(i=1; i<nx-1; i++) {
      if( (mask->data)[i+j*nx] == 2) {
	if( (mask->data)[i+1+ j   *nx] == 0 ||
	    (mask->data)[i+1+(j+1)*nx] == 0 ||
	    (mask->data)[i  +(j+1)*nx] == 0 ||
	    (mask->data)[i-1+(j+1)*nx] == 0 ||
	    (mask->data)[i-1+ j   *nx] == 0 ||
	    (mask->data)[i-1+(j-1)*nx] == 0 ||
	    (mask->data)[i  +(j-1)*nx] == 0 ||
	    (mask->data)[i+1+(j-1)*nx] == 0) {
	  (mask->data)[i+j*nx] = 0;
	} else {
	  (mask->data)[i+j*nx] = 1;
	}
      }
    }
  }
}

syntax(s)
char *s;
{
  printf("Syntax: superflat [options] source_files...\n");
  printf("	-out fname	output file name\n");
  printf("	-imx N		Size of input images in x direction\n");
  printf("	-imy N		Size of input images in y direction\n");
  printf("	-sx N		Starting point of image in x direction\n");
  printf("	-sy N		Starting point of image in y direction\n");
  printf("	-border N	Number of bad border pixels\n");
  printf("	-nmem N		Number of frames which fit in memory\n");
  printf("	-scrunch N	Use scrunch factor N\n");
  printf("	-bias X		Use X for bias value\n");
  printf("	-gain X		Use X e/ADU\n");
  printf("	-scramble	Scramble order of use of input files\n");
  printf("	-test		Write test output files\n");
}
