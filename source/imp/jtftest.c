/* Q+D to test jtfits.c */

#include <stdio.h>

#define ABS(a) (((a) > 0) ? (a) : -(a))

main(argc,argv)
int argc;
char **argv;
{
  int i, j, n, err, arg;
  int nhead, bitpix, bp_data, naxis, dims[16];
  int size[16], offset[16];
  int npix, maxhead, maxdata, bp_disk;
  double scale, zero, scout, zout, atof();
  char *head, *data;
  char *fin, *fout;

  bp_data = bp_disk = -32;
  scout = 1.0;
  zout = 0.0;

  size[0] = offset[0] = size[1] = offset[1] = 0;

/* Parse the arguments */
  for(arg=1; arg<argc; arg++) {

    if(strncmp(argv[arg], "-data", 5) == 0) {
      bp_data = atoi(argv[++arg]);

    } else if(strncmp(argv[arg], "-disk", 5) == 0) {
      bp_disk = atoi(argv[++arg]);

    } else if(strncmp(argv[arg], "-in", 3) == 0) {
      fin = argv[++arg];

    } else if(strncmp(argv[arg], "-out", 3) == 0) {
      fout = argv[++arg];

    } else if(strncmp(argv[arg], "-scale", 6) == 0) {
      scout = atof(argv[++arg]);

    } else if(strncmp(argv[arg], "-zero", 5) == 0) {
      zout = atof(argv[++arg]);

    } else if(strncmp(argv[arg], "-x0", 3) == 0) {
      offset[0] = atoi(argv[++arg]);

    } else if(strncmp(argv[arg], "-nx", 3) == 0) {
      size[0] = atoi(argv[++arg]);

    } else if(strncmp(argv[arg], "-y0", 3) == 0) {
      offset[1] = atoi(argv[++arg]);

    } else if(strncmp(argv[arg], "-ny", 3) == 0) {
      size[1] = atoi(argv[++arg]);

    } else {
      fprintf(stderr, "Unknown arg: %s\n", argv[arg]);
      break;
    }   
  }

/* Have we got a FITS file, and how much space do we need? */
/* What is the scale and zero in the header?  (Not that I care...) */
  if((err=testfits(fin, &nhead, &bitpix, &naxis, dims, &scale, &zero)) != 0) {
    fprintf(stderr, "testfits failed with error %d\n", err);
    exit(1);
  }
  printf("bitpix = %d  naxis = %d  nhead = %d\n",
	 bitpix, naxis, nhead);
  printf("dims =");
  for(i=0; i<naxis; i++) printf("%5d x", dims[i]); printf("\n");

/* Allocate some space for header and data */
  allocfits(nhead+30, &head, naxis, dims, bp_data, &data);
  maxhead = 80 * (nhead+30);
  maxdata = (ABS(bp_data)/sizeof(char));  
  for(i=0; i<naxis; i++) maxdata *= dims[i];
  printf("space allocated for header (%d) and data (%d)\n", maxhead, maxdata);

/* Now read in the data and header */
  if(size[0] == 0 || size[1] == 0) {
    rfits(maxhead, &nhead, head, maxdata, bp_data, &naxis, dims, data, fin);
  } else {
    rfsubarray(maxhead, head, maxdata, bp_data, offset, size, data, fin);
  }

  printf("Data read:  nhead = %d  naxis = %d\n", nhead, naxis);
  printf("dims =");
  for(i=0; i<naxis; i++) printf("%5d x", dims[i]); printf("\n");

/* Show me the first 8 pixels, please */
  switch (bp_data) {
  case -32:
    for(i=0; i<8; i++) printf(" %8g", *((float *)(data+i*sizeof(float))));
    break;
  case 32:
    for(i=0; i<8; i++) printf(" %8d", *((int *)(data+i*sizeof(int))));
    break;
  case -16:
    for(i=0; i<8; i++) printf(" %8d", *((unsigned short int *)(data+i*sizeof(unsigned short int))));
    break;
  case 16:
    for(i=0; i<8; i++) printf(" %8d", *((short *)(data+i*sizeof(short))));
    break;
  default:
    fprintf(stderr,"Cannot handle bp_data = %d\n", bp_data);
  }
  printf("\n");

/* Set the desired output SCALE and ZERO */
  if(chfitshead(&n, head, "BSCALE  ", "FLOAT", 0, scout) != 0)
    addfitshead(&n, head, "BSCALE  ", "FLOAT", 0, scout);
  if(chfitshead(&n, head, "BZERO   ", "FLOAT", 0, zout) != 0)
    addfitshead(&n, head, "BZERO   ", "FLOAT", 0, zout);
  printf("Output scale/zero changed to %g, %g\n", scout, zout);

/* Set the desired output BITPIX */
  chfitshead(&n, head, "BITPIX  ", NULL, bp_disk, 0.0);
  printf("Output bitpix changed to %d\n", bp_disk);

/* Write the output file */
  wfits(maxhead, head, bp_data, data, fout);
}

rfitsreal(head, nx,ny, data, file)
char **head;		/* header */
char *file;		/* file name */
float **data;		/* image data */
int *nx, *ny;		/* NAXIS1 and NAXIS2 from header */
{
  int nhead, bitpix, bp_data, naxis, dims[16];
  int err, maxhead, maxdata;
  double scale, zero;

  bp_data = -32;

/* Have we got a FITS file, and how much space do we need? */
  if((err=testfits(file, &nhead, &bitpix, &naxis, dims, &scale, &zero)) != 0) {
    return(err);
  }
  *nx = dims[0];
  *ny = dims[1];

/* Allocate some space for header and data */
  allocfits(nhead+30, head, naxis, dims, bp_data, data);
  maxhead = 80 * (nhead+30);
  maxdata = (*nx) * (*ny) * (ABS(bp_data)/sizeof(char));  

/* Now read in the data and header */
  err = rfits(maxhead, &nhead, *head, maxdata, bp_data, &naxis, dims, *data, file);
  return(err);
}

/*
 * wfitsreal() will write a real FITS image onto disk.
 */
wfitsreal(head, data, file)
char *head;		/* header */
char *file;		/* file name */
char *data;		/* image data */
{
  int maxhead, n;

/* Set the desired output BITPIX */
  chfitshead(&n, head, "BITPIX  ", NULL, -32, 0.0);

/* Write the output file */
  maxhead = 80*counthead(head);
  wfits(maxhead, head, -32, data, file);
}
