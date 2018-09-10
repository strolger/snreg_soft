#include<stdio.h>
#include<math.h>
#include<string.h>
#define HEADLINES 1000
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))

main (argc,argv)
     int argc;
     char *argv[];
{  

  FILE *stream1,*stream2,*stream3,*stream4;
  char maskim[100],superim[100],imageout[100],indexx;
  char *header,*header1,*header2;
  int n,xaxis,yaxis,xaxis1,yaxis1,xaxis2,yaxis2,sf=0;
  int bitpix,nxorig,nyorig,min,max;
  double scale,zero;

  long i,N;
  float low,hi,mid;
  double X,Y,XX,XY,alpha,beta; 
  float *mask,*outim,*superflat,percentimage();

  
  if (argc < 4) {
    syntax(*argv);
    exit(0);
  }   

  sprintf(maskim,argv[1]);
  sprintf(superim,argv[2]);
  sprintf(imageout,argv[3]);
  scale = atoi(argv[4]);

  
  rfitsreal(&header, &xaxis, &yaxis, &mask, maskim);
  printf("Read %s, size = (%d, %d)...\n", maskim, xaxis, yaxis);
  rfitsreal(&header1, &xaxis1, &yaxis1, &superflat, superim);
  printf("Read %s, size = (%d, %d)...\n", superim, xaxis1, yaxis1);
  
  if (xaxis != xaxis1 || yaxis != yaxis1) {
    printf("Images do not have same dimensions...Aborting\n");
    exit(-1);
  }
  
  for (i=0;i<xaxis*yaxis;i++) {
    if (superflat[i] > 0) {    superflat[i]=scale/superflat[i]*mask[i]; }  
    else { superflat[i]=0;}
  }
  printf("Writing output file %s...\n", imageout);
  bitpix=-32;
  /* Create a FITS header to write image */

/* Alter the old FITS header for the defringed image */

  /* if bitpix defined above, change the header */
  if(bitpix != 0) {
    chfitshead(&n, header, "BITPIX  ", NULL, bitpix, 0.0);
  } else { /* make it what it was*/
    if(ifitshead(header, "BITPIX  ", &bitpix) != 0) bitpix = -32;
  }
  /*bscale and bzero are 1,0 for real*/
  if(bitpix == -32) {
    chfitshead(&n, header, "BSCALE  ", "FLOAT", 0, 1.0);
    chfitshead(&n, header, "BZERO   ", "FLOAT", 0, 0.0);
    wfitsreal(header, superflat, imageout);
  } else if( ABS(bitpix) == 16 ) {
    wfitsshort(header, superflat, imageout);
  }
}
 
syntax(s)
char *s;
{
  printf("Syntax: mksuperflat maskfile superflat outfile scale\n\n");
  printf("Make flat for multiply = scale/superflat*mask\n");
}


