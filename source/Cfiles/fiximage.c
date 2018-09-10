#include<stdio.h>
#include<stdlib.h>
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
  char imagein[100],fringeim[100],imageout[100],sfim[100],indexx;
  char *header,*header1,*header2;
  int n,nt,xaxis,yaxis,xaxis1,yaxis1,xaxis2,yaxis2,sf=0;
  int bitpix,nxorig,nyorig;
  float min,max;
  float scale,offset;

  long i,N,ctlow,cthi,ctmed;
  float low,hi,mid;
  double X,Y,XX,XY,alpha,beta; 
  float *image,*fringe,*superflat,percentimage();

  
  if (argc != 8) {
    syntax(*argv);
    exit(0);
  }   

  ctlow = cthi = ctmed = 0;

  sprintf(imagein,argv[1]);
  sprintf(imageout,argv[2]);
  min = atof(argv[3]);
  max = atof(argv[4]);
  offset = atof(argv[5]);
  scale = atof(argv[6]);

  if (!strcmp(argv[7],"real")) {
    bitpix = -32;
  } else if (!strcmp(argv[7],"ushort")) {
    bitpix = 16;
  } else if (!strcmp(argv[7],"autoclean")) {
    bitpix = -16;
  } else {
    printf("\n     pixtype must be real, ushort or autoclean!\n\n");
    syntax(*argv);
    exit(1);
  }
   
  printf("reading %s...\n",imagein);
  if (rfitsreal(&header, &xaxis, &yaxis, &image, imagein)) {
    printf("\n     error reading input file!\n");
    exit(1);
  }
  printf("read %s, size = (%d, %d)...\n", imagein, xaxis, yaxis);
   
  printf("adding %.1f to all pixels\n",offset);
  printf("then multiplying all pixels by %.5f...\n",scale);
  printf("then fixing all pixels outside range of %.1f and %.1f\n",min,max);
  for (i=0;i<xaxis*yaxis;i++) {


    image[i] += offset;
    image[i] *= scale;
    
    if (image[i] < min) {
      image[i] = min;
      ctlow++;
    } else if (image[i] > max) {
      image[i] = max;
      cthi++;
    } else {
      ctmed++;
    }

  }
  printf("pixels low: %d  good: %d  high: %d\n",ctlow,ctmed,cthi);
  printf("writing output file %s...\n", imageout);
  /* Create a FITS header to write image */

/* Alter the old FITS header for the new image */
  chfitshead(&n, header, "BITPIX  ", NULL, bitpix, 0.0);

  chfitshead(&n, header, "BSCALE  ", "FLOAT", 0, 1.0);
  if (n == -1) {        /* no BSCALE in header */
    n = counthead(header);
    addfitshead(&n, header, "BSCALE  ", "FLOAT", 0, 1.0);
  }

  /*bscale and bzero are 1,0 for real*/
  if(bitpix == -32) {
 
    chfitshead(&n, header, "BZERO   ", "FLOAT", 0, 0.0);
    if (n == -1) {        /* no BZERO in header */
      n = counthead(header);
      addfitshead(&n, header, "BZERO   ", "FLOAT", 0, 0.0);
    }

    wfitsreal(header, image, imageout);

  } else if( bitpix == 16 ) {

    chfitshead(&n, header, "BZERO   ", "FLOAT", 0, 32768.0);
    if (n == -1) {        /* no BZERO in header */
      n = counthead(header);
      addfitshead(&n, header, "BZERO   ", "FLOAT", 0, 32768.0);
     }

    wfitsshort(header, image, imageout);
    
  } else if ( bitpix == -16) {

    chfitshead(&n, header, "BZERO   ", "FLOAT", 0, 0.0);
    if (n == -1) {        /* no BZERO in header */
      n = counthead(header);
      addfitshead(&n, header, "BZERO   ", "FLOAT", 0, 0.0);
    }

    wfitsshort(header, image, imageout);
    
  } else {  /* this should never happen */
    printf ("Invalid bitpix: %d, aborting!\n",bitpix);
    exit (1);
  }
  
}
 
syntax(s)
char *s;
{
  printf("Syntax: fiximage infile outfile min max offset scale pixtype\n\n");
  printf(" infile, outfile: fits image names (they can be the same)\n");
  printf(" min, max, offset, scale: floating point numbers\n");
  printf(" pixtype:\n");
  printf("\t\treal       floating point data output\n");
  printf("\t\tushort     IRAF fake ushort (signed 16-bit with bscale=32k)\n");
  printf("\t\tautoclean  true ushort (unsigned 16-bit, 0-65535)\n\n");
  printf(" first reads infile, \n");
  printf(" then adds offset to all pixels\n");
  printf(" then multiplies all pixels by scale\n");
  printf(" then sets all pixels < min (> max) to min (max) [truncation]\n");
  printf(" and finally writes out a pixtype image as outfile\n\n");
}

