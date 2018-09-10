/* Various utility routines for jtfits.c */

#include <stdio.h>

#define ABS(a) (((a) > 0) ? (a) : -(a))

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

wfitsreal_(head, data, file, headlen, flen)
char *head;		/* header */
char *data;		/* image data */
char *file;		/* file name */
int headlen, flen;
{
  char fname[256];
  if(flen >= 256) {
    fprintf(stderr, "Cannot write such a long file name %s\n", file);
    return(1);
  }
  while(--flen > 0 && (file[flen] == ' ' || file[flen] == '\0'));
  strncpy(fname, file, flen+1);
  fname[flen+1] = '\0';

  wfitsreal(head, data, fname);
}

newfhead_(n, header, nx, ny, headlen)
int *n, *nx, *ny, headlen;
char *header;
{
  int i;
  for(i=0; i<6*80; i++) header[i] = ' ';
  sprintf(header+0*80, "SIMPLE  =                    T");
  sprintf(header+1*80, "BITPIX  =                   16");
  sprintf(header+2*80, "NAXIS   =                    2");
  sprintf(header+3*80, "NAXIS1  = %20d", *nx);
  sprintf(header+4*80, "NAXIS2  = %20d", *ny);
  sprintf(header+5*80, "END");
  for(i=0; i<5; i++) header[i*80+30] = ' ';
  header[5*80+3] = ' ';
  *n = 6;
  return(0);
}
