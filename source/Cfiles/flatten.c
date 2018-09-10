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
  char imagein[100],flatin[100],headstr[80],indexx;
  char *header,*fheader;
  int n,nt,xaxis,yaxis,fxaxis,fyaxis,sf=0;
  int bitpix,nxorig,nyorig;
  float min=0.0,max=65535.0;
  float scale,offset;

  long i,N,j;
  float low,hi,mid;
  double X,Y,XX,XY,alpha,beta; 
  float *image,*flat;

  
  if (argc < 4) {
    syntax(*argv);
    exit(0);
  }   

  sprintf(flatin,argv[1]);
  scale = (float) atof(argv[2]);

  bitpix = 16;
   
  printf("reading %s...\n",flatin);
  if (rfitsreal(&fheader, &fxaxis, &fyaxis, &flat, flatin)) {
    printf("\n     error reading flatfield file!\n");
    exit(1);
  }
  printf("read %s, size = (%d, %d)...\n", flatin, fxaxis, fyaxis);


  for (j=3;j<argc;j++) {

    sprintf(imagein,argv[j]);
    printf("reading %s...\n",imagein);
    if (rfitsreal(&header, &xaxis, &yaxis, &image, imagein)) {
      printf("\n     error reading input file!\n");
      exit(1);
    }
    printf("read %s, size = (%d, %d)...\n", imagein, xaxis, yaxis);
  
    if ( (xaxis != fxaxis) || (yaxis != fyaxis) ) {
      printf("\n %s has different size than flatfield, skipping it!\n",imagein);
      continue;
    }
 
    printf("dividing %s by %s and multiplying by %g\n",imagein,flatin,scale);
    printf("then truncating to ushort (0-65535)\n");

    for (i=0;i<xaxis*yaxis;i++) {

      image[i] *= (scale/flat[i]);
    
      if (image[i] < min) {
	image[i] = min;
      } else if (image[i] > max) {
	image[i] = max;
      }

    }
  
    printf("writing output file %s...\n", imagein);
    /* Create a FITS header to write image */

    /* Alter the old FITS header for the new image */
    chfitshead(&n, header, "BITPIX  ", NULL, bitpix, 0.0);

    chfitshead(&n, header, "BSCALE  ", "FLOAT", 0, 1.0);
    if (n == -1) {        /* no BSCALE in header */
      n = counthead(header);
      addfitshead(&n, header, "BSCALE  ", "FLOAT", 0, 1.0);
    }

    chfitshead(&n, header, "BZERO   ", "FLOAT", 0, 32768.0);
    if (n == -1) {        /* no BZERO in header */
      n = counthead(header);
      addfitshead(&n, header, "BZERO   ", "FLOAT", 0, 32768.0);
    }

    sprintf(headstr,"flattened by %s, scaled by %g",flatin,scale);
    n = counthead(header);
    addfitshead(&n, header, "HISTORY ", headstr, 0, 0); 
    wfitsshort(header, image, imagein);
    free(image);
    free(header);
  }

}
 
syntax(s)
char *s;
{
  printf("Syntax: flatten flatfield scale image1 image2 ...\n\n");
  printf(" flattens images by dividing by flatfield\n");
  printf(" then multiplies by scale\n");
  printf(" then truncates to ushort (0-65535)\n\n");
}

