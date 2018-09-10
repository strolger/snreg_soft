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
  char imagein[100],fringeim[100],imageout[100],sfim[100],indexx;
  char *header,*header1,*header2;
  int n,xaxis,yaxis,xaxis1,yaxis1,xaxis2,yaxis2,sf=0;
  int bitpix,nxorig,nyorig;
  double scale,zero;

  long i,N;
  float low,hi,mid;
  double X,Y,XX,XY,alpha,beta; 
  float *image,*fringe,*superflat,percentimage();

  
  if (argc < 4) {
    syntax(*argv);
    exit(0);
  }   

  sprintf(imagein,argv[1]);
  sprintf(fringeim,argv[2]);
  sprintf(imageout,argv[3]);

  for(i=4; i<argc; i++) {
    if(strncmp(argv[i], "-bitpix", 7) == 0) {
      bitpix = atoi(argv[++i]);
    }
    else {
      if(strncmp(argv[i], "-superflat", 10) == 0) {
	sf=1;
	strcpy(sfim,argv[++i]);
      }
      else {
	syntax(*argv);
	exit(-1);
      }
    }
  } 
  
  rfitsreal(&header, &xaxis, &yaxis, &image, imagein);
  printf("Read %s, size = (%d, %d)...\n", imagein, xaxis, yaxis);
  
  rfitsreal(&header1,&xaxis1, &yaxis1, &fringe, fringeim);
  printf("Read %s, size = (%d, %d)...\n", fringeim, xaxis1, yaxis1);
  
  if (sf ==1) {
    rfitsreal(&header2,&xaxis2, &yaxis2, &superflat, sfim);
    printf("Read %s, size = (%d, %d)...\n", sfim, xaxis2, yaxis2);
    if (xaxis != xaxis2 || yaxis != yaxis2) {
      printf("superflat and image are not the same size, abort!\n"); 
      exit (-1);
    }
  }
  if (xaxis != xaxis1 || yaxis != yaxis1) {
    printf("fringe and image are not the same size, abort!\n"); 
    exit (-1);
  }
  /*First Flatfield if given*/
  if (sf ==1) {
    printf("Applying Flatfield...");
    for (i=0;i<xaxis*yaxis;i++) {
      image[i] *= superflat[i];
    }
    printf("Done\n");
  }  


  low=percentimage(image,xaxis/2.,yaxis/2.,xaxis,yaxis,(int)xaxis/4,0.0027);
  mid=percentimage(image,xaxis/2.,yaxis/2.,xaxis,yaxis,(int)xaxis/4,0.5);
  hi=percentimage(image,xaxis/2.,yaxis/2.,xaxis,yaxis,(int)xaxis/4,0.9973);

  N=X=XX=Y=XY=0.;
  for (i=0;i<xaxis*yaxis;i++) {
    if (image[i] > low && image[i] < hi) { /* only use those near median*/
      X +=fringe[i];
      XX+=fringe[i]*fringe[i];
      Y +=image[i];
      XY+=fringe[i]*image[i];
      N++;
    }
  }
  /*  Least squares to solve
      Image = beta + alpha*fringe
      i.e. (Image - background) - alpha*fringe=0; 
  */
  beta = (Y*XX-X*XY)/(N*XX-X*X);
  alpha =  (N*XY-X*Y)/(N*XX-X*X);
  if (alpha > 3) alpha=2.99;
  if (alpha < 0) alpha=0.;
  printf ("Image background: %7.3f (%7.3f from median) Fringe Scale: %7.3f\n",beta,mid,alpha);   
/*correct for above fit*/
  for (i=0;i<xaxis*yaxis;i++) {
    if (image[i] !=0) image[i] -= alpha*fringe[i];
  }
  printf("Writing output file %s...\n", imageout);
  
  /* Create a FITS header to write image */

/* Alter the old FITS header for the defringed image */

/* if bitpix defined above, change the header */
  if(bitpix != 0) {
    chfitshead(&n, header, "BITPIX  ", NULL, bitpix, 0.0);
  } else { /* make it what it was*/
    if(ifitshead(header, "BITPIX  ", &bitpix) != 0) bitpix = -32;
  }
  if(bitpix == 16) { /*want ushorts if at all possible*/
    chfitshead(&n, header, "BSCALE  ", "FLOAT", 0, 1.0);
    chfitshead(&n, header, "BZERO   ", "FLOAT", 0, 32768.);
  }
  /*bscale and bzero are 1,0 for real*/
  if(bitpix == -32) {
    chfitshead(&n, header, "BSCALE  ", "FLOAT", 0, 1.0);
    chfitshead(&n, header, "BZERO   ", "FLOAT", 0, 0.0);
    wfitsreal(header, image, imageout);
  } else if( ABS(bitpix) == 16 ) {
    wfitsshort(header, image, imageout);
  }
}
   

float percentimage(data,x,y,nx,ny,s,frac)
     int nx,ny,s;  
     float x,y,frac;
     float *data;
{
  int j,k,jlo,jhi,klo,khi,index[15000],N,i;
  float background,medpix[15000];
  
  jlo = y-s;
  jhi = y+s;		
  klo = x-s;
  khi = x+s;
  
  if (jlo < 0) jlo =0;
  if (klo < 0) klo = 0;
  if (jhi >= ny) jhi = ny;
  if (khi >= nx) khi = nx;
  
  
  for (N=1,i=klo; i<=khi ;i+=(s/50+1)) { /* go over x*/
    for (j=jlo; j<=jhi;j+=(s/50+1)) { /* go over y*/
      if (data[j*nx+i] != 0) { /* only use good pixels */
	medpix[N] = data[j*nx+i];
	N++;
      }
    }
  }  
  N--;
  if (N > 10) {
    indexx(N,medpix,index); /* sort to get median,etc.*/
    background = medpix[(int)(N*frac)+1];
    return(background);
  }
  else return(-9999);
}
 
syntax(s)
char *s;
{
  printf("Syntax: defringe infile fringefile outfile -superflat file -bitpix N  (16/-16/-32)\n");
  printf("fringe frame will be scaled to image and subtracted, output into outfile in same format (real/short) as input if not specified\n. Superflat will be multiplied through if exists");
}


