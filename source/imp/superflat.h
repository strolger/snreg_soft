/* Program to create a superflat from a pile of FITS images */
/* John Tonry - 26 October 1999 */

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))

#define MAXFILE 256
#define MEANFLAT 10000

struct ushortfits {
  char *fname;			/* Associated file name */
  int nhead;			/* Number of header records */
  char *header;			/* Header records */
  int nx;			/* Size in x direction */
  int ny;			/* Size in y direction */
  int ix;			/* Image size in x direction */
  int iy;			/* Image size in y direction */
  int sx;			/* Image start in x direction */
  int sy;			/* Image start in y direction */
  int bitpix;			/* Bitpix parameter */
  int scrunch;			/* Scrunch factor */
  unsigned short int *data;	/* Image data */
  int headstore;		/* Storage allocated for header (bytes) */
  int datastore;		/* Storage allocated for data (bytes) */
  int bias;			/* Bias level */
  int sky;			/* Mean sky level */
};

