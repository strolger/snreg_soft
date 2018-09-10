/* jtfits.c: routines to read and write FITS files
 * 
 * jtfits 1.01: 00-09-12   working out the bugs...
 * jtfits 1.00: 99-12-05   initial rev
 *
 * branched off from rwfits 2.1 99-12-02
 *
 * rwfits 2.1 - 11/02/99 JT Merge in bug fixes from Saurabh Jha
 * rwfits 2.0 - 9/21/99 JT Add header and convenience routines like rfitsreal()
 * rwfits 1.1 - 9/07/99 Add code for overlapping swab()
 * rwfits 1.0 - 7/08/92  initial rev
 *
 * Routines to read and write disk FITS format data.
 * These files have a header which consists of 36*n 80 byte header records
 * which are a duplicate of the FITS tape format, with a mandatory END record.
 * This is followed by N*2880 bytes of data.
 *
 * The conversion routines look to see whether LOWENDIAN and VAXFP are defined
 * Set them if necessary!
 * E.g.    VAX:        -DLOWENDIAN -DVAXFP
 *         DEC MIPS:   -DLOWENDIAN
 *         Intel:      -DLOWENDIAN
 *         Sun:
 *
 * NOTE: problems have arisen with DEC Alpha's running OSF1 because OSF1
 * has a swab() which doesn't work in place!  You can define OSF_ALPHA to
 * cause a swab() to be compiled here which works right:
 *         OSF1 ALPHA:  -DLOWENDIAN -DOSF_ALPHA
 *
 *
 * rwfits.c: John Tonry 7/8/92
 *
 *
 */

#include <stdio.h>
#include <math.h>

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))

#define NFITS (2880)
#define MAXBUF (4*65536)

static char buf[MAXBUF];

#define ERR_CANT_OPEN_FILE		1
#define ERR_CANT_READ_HEADER		2
#define ERR_CANT_READ_DATA		3
#define ERR_CANT_WRITE_HEADER		4
#define ERR_CANT_WRITE_DATA		5
#define ERR_INSUFFICIENT_HEADER		6
#define ERR_NO_NAXIS			7
#define ERR_NO_END			8
#define ERR_INCONSISTENT_NAXIS		9
#define ERR_BUFFER_OVERFLOW	       10
#define ERR_NOT_A_FITS_FILE	       11
#define ERR_ILLEGAL_BITPIX	       12
#define ERR_CANT_SEEK		       13

/*
 * rfits() reads FITS image from disk into existing header and data array
 */
rfits(maxhead, nhead, head, maxdata, bp_data, naxis, dim, data, file)
int maxhead;		/* header storage available (bytes) */
int *nhead;		/* header count */
char *head;		/* header */
int maxdata;		/* data storage available (bytes) */
int bp_data;		/* BITPIX desired for data */
int *naxis;		/* number of dimensions of FITS image */
int *dim;		/* sizes of each dimension */
char *data;		/* image data */
char *file;		/* file name */
{
  int err, npix, i;
  double scale, zero;	/* BSCALE and BZERO parameters from header */
  int bitpix;		/* BITPIX of data in disk file */
  int fd;

  if( (fd=open(file, 0)) == -1) return(ERR_CANT_OPEN_FILE);

  if( (err=rhead(fd, maxhead, head)) != 0) {
    close(fd); 
    return(err);
  }

  if((err=parsehead(maxhead,nhead,head,&bitpix,naxis,dim,&scale,&zero))!= 0) {
    close(fd);    
    return(err);
  }
  
  for(i=0, npix=1; i<*naxis; i++) npix *= dim[i];

  if(2*((ABS(bp_data)*npix+15)/16) > maxdata) return(ERR_BUFFER_OVERFLOW);

  if( (err=diskread(fd, npix, bitpix, bp_data, scale, zero, data)) != 0) {
    close(fd);
    return(err);
  }

  close(fd);
  return(0);
}

/*
 * rfsubarray reads a portion of a FITS image from disk into a data array
 */
rfsubarray(maxhead, head, maxdata, bp_data, offset, size, data, file)
int maxhead;		/* header storage available (bytes) */
char *head;		/* header */
int maxdata;		/* data storage available (bytes) */
int bp_data;		/* BITPIX desired for data */
int *offset;		/* offset in each dimension of desired read */
int *size;		/* sizes in each dimension of desired read */
char *data;		/* image data */
char *file;		/* file name */
{
  int err, npix, i, j, nhead, naxis, dim[16], n;
  double scale, zero;	/* BSCALE and BZERO parameters from header */
  int bitpix;		/* BITPIX of data in disk file */
  int fd;
  char keyword[9];
  int seekloc, curloc, dtot, stot, stride;

  if( (fd=open(file, 0)) == -1) return(ERR_CANT_OPEN_FILE);
  if( (err=rhead(fd, maxhead, head)) != 0) {
    close(fd); 
    return(err);
  }
  if((err=parsehead(maxhead,&nhead,head,&bitpix,&naxis,dim,&scale,&zero))!= 0){
    close(fd);
    return(err);
  }
  
  for(i=0, npix=1; i<naxis; i++) npix *= size[i];
  stride = size[0];
  if(size[0] == dim[0]) stride *= size[1];

  if(2*((ABS(bp_data)*npix+15)/16) > maxdata) return(ERR_BUFFER_OVERFLOW);

  for(i=0, curloc=0; i<npix; i+=stride) {
    for(j=0, seekloc=0, dtot=1, stot=1; j<naxis; j++) {
      seekloc += dtot*(offset[j] + (i/stot)%size[j]);
      dtot *= dim[j];
      stot *= size[j];
    }
    seekloc *= ABS(bitpix) / 8;

    if(lseek(fd, seekloc-curloc, 1) < seekloc-curloc) {
      close(fd);
      return(ERR_CANT_SEEK);
    }
    if( (err=diskread(fd, stride, bitpix, bp_data, scale, zero, 
		      &data[(i*ABS(bp_data))/8])) != 0) 
      return(err);
    curloc = seekloc + (ABS(bitpix)*stride)/8;
  }

  close(fd);
/* OF COURSE IF WE WERE REALLY CAREFUL WE WOULD CHANGE NAXIS IF 0 SIZE */
/* Fix up the header to conform with the new reality */
  for(j=0; j<naxis; j++) {
/* Set the NAXISn dimensions */
    sprintf(keyword,"NAXIS%d  ", j+1);
    chfitshead(&i, head, keyword, NULL, size[j], 0.0);
/* and use CNPIXn as the offset */
    sprintf(keyword,"CNPIX%d  ", j+1);
    if(ifitshead(head, keyword, &n) != 0) n = 0;
    n += offset[j];
    if(chfitshead(&i, head, keyword, NULL, n, 0.0) != 0)
      addfitshead(&i, head, keyword, NULL, n, 0.0);
  }
  return(0);
}

/*
 * wfits() writes a FITS image onto disk.
 */
wfits(maxhead, head, bp_data, data, file)
int maxhead;		/* storage allocated for header (bytes) */
char *head;		/* header */
int bp_data;		/* BITPIX of data array */
char *data;		/* image data */
char *file;		/* file name */
{
  int nhead;		/* header count */
  int fd, mode=1, err, bitpix, naxis, dim[16];
  int nwrite, nbyte, imagebyte, nwritten, npix, i;
  double scale, zero;

  if( (fd=creat(file, 0644)) == -1) return(ERR_CANT_OPEN_FILE);

/* Find out desired output BITPIX, BSCALE, and BZERO */
  if( (err=parsehead(maxhead,&nhead,head,&bitpix,&naxis,dim,&scale,&zero)) != 0) {
    close(fd);    return(err);  }

  if( (err=whead(fd, nhead, head)) != 0) {
    close(fd);  return(err); }

  for(i=0, npix=1; i<naxis; i++) npix *= dim[i];

/* Set scale and zero for output conversion routines */
  if(scale != 1.0) scale = 1/scale;
  if(zero != 0.0) zero = -zero;

  imagebyte = 2*((ABS(bitpix)*npix+15)/16);
  nbyte = imagebyte;		       	/* Bytes still to write */
  nwritten = 0;				/* Bytes written */
  while(nbyte > 0) {
    nwrite = MIN(MAXBUF,nbyte);
/* Convert from machine specific type to desired output FITS format */
    npix = (8*nwrite)/ABS(bitpix);
    fitsconvert(2, npix, scale, zero, bp_data, 
		 &data[(ABS(bp_data)*nwritten)/ABS(bitpix)], bitpix, buf);

    if( (err=write(fd, buf, nwrite)) != nwrite) {
      close(fd); 
      return(ERR_CANT_WRITE_DATA);
    }
    nbyte -= nwrite;
    nwritten += nwrite;
  }

  nwrite = NFITS - (imagebyte%NFITS);
  if(nwrite < NFITS) {
    bzero(buf,nwrite);
    if( (err=write(fd, buf, nwrite)) != nwrite) {
      close(fd); return(ERR_CANT_WRITE_DATA); }
  }

  close(fd);
  return(0);
}

/* 
 * Read npix pixels from the current pointer in fd; format into data array
 */
diskread(fd, npix, bitpix, bp_data, scale, zero, data)
int fd, npix, bitpix, bp_data;
double scale, zero;
char *data;
{
  int nbyte;					/* Bytes yet to read */
  int totpix = 0;				/* Pixels read */
  int thisread, thispix;

  nbyte = 2*((ABS(bitpix)*npix+15)/16);

  while(nbyte > 0) {
    thisread = MIN(MAXBUF, nbyte);
    if( read(fd, buf, thisread) != thisread) {
      close(fd); 
      return(ERR_CANT_READ_DATA);
    }
/* Convert from machine specific type to desired output FITS format */
    thispix = (8*thisread)/ABS(bitpix);
    fitsconvert(1, thispix, scale, zero, bitpix, buf, 
		bp_data, &data[(ABS(bp_data)/8)*totpix]);
    nbyte -= thisread;
    totpix += thispix;
  }
  return(0);
}

fitsconvert(fitswab, npix, scale, zero, bp_src, src, bp_dest, dest)
int fitswab;		/* bit 0,1 = 0/1 for fitsorder on src/dest */
int npix;		/* Number of pixels */
int bp_src, bp_dest;	/* BITPIX of source and destination */
double scale, zero;	/* BSCALE and BZERO parameters from header */
char *src, *dest;
{
/* Check for disallowed BITPIX */
  if(ABS(bp_dest) != 32 || ABS(bp_src) != 32) {
    if(bp_dest != -32 && bp_dest != -16 && 
       bp_dest != 16 &&  bp_dest != 1) {
      fprintf(stderr,"Cannot convert destination BITPIX = %d\n", bp_dest);
      return(ERR_ILLEGAL_BITPIX);
    }
    if(bp_src != -32 && bp_src != -16 && 
       bp_src != 16 &&  bp_src != 1) {
      fprintf(stderr,"Cannot convert source BITPIX = %d\n", bp_src);
      return(ERR_ILLEGAL_BITPIX);
    }
  }

  /*
  fprintf(stderr, "fitsconvert: fitswab = %d  npix = %d  bp_src = %d  bp_dest = %d  s,z = %g %g\n", fitswab, npix, bp_src, bp_dest, scale, zero);
  */

/* If input is FITSendian, do the necessary swabbing */
  if(fitswab & 1) {
/* Separate out the case involving IEEE -> fp because of VAX possibility */
    if(bp_src == -32) ieeefp(npix, src, src);
#ifdef LOWENDIAN
    else FITSorder(bp_src, npix, src);
#endif
  }

  switch(bp_src) {
  case -32:
    if(bp_dest == -32) {
      if(src != dest) bcopy(src, dest, sizeof(float)*npix);
    } else if(bp_dest == 32) {
      FIconvert(npix, scale, zero, src, dest);
    } else if(bp_dest == -16) {
      FUconvert(npix, scale, zero, src, dest);
    } else if(bp_dest == 16) {
      FSconvert(npix, scale, zero, src, dest);
    } else if(bp_dest == 1) {
      FBconvert(npix, scale, zero, src, dest);
    }
    break;

  case 32:
    if(bp_dest == -32) {
      IFconvert(npix, scale, zero, src, dest);
    } else if(bp_dest == 32) {
      if(src != dest) bcopy(src, dest, sizeof(int)*npix);
    }
    break;

  case -16:
    if(bp_dest == -32) {
      UFconvert(npix, scale, zero, src, dest);
    } else if(bp_dest == -16) {
      UUconvert(npix, scale, zero, src, dest);
    } else if(bp_dest == 16) {
      USconvert(npix, scale, zero, src, dest);
    } else if(bp_dest == 1) {
      SBconvert(npix, scale, zero, src, dest);
    }
    break;

  case 16:
    if(bp_dest == -32) {
      SFconvert(npix, scale, zero, src, dest);
    } else if(bp_dest == -16) {
      SUconvert(npix, scale, zero, src, dest);
    } else if(bp_dest == 16) {
      SSconvert(npix, scale, zero, src, dest);
    } else if(bp_dest == 1) {
      SBconvert(npix, scale, zero, src, dest);
    }
    break;

  case 1:
    if(bp_dest == -32) {
      BFconvert(npix, scale, zero, src, dest);
    } else if(bp_dest == -16) {
      BSconvert(npix, scale, zero, src, dest);
    } else if(bp_dest == 16) {
      BSconvert(npix, scale, zero, src, dest);
    } else if(bp_dest == 1) {
      if(src != dest) bcopy(src, dest, 2*((npix+15)/16) );
    }
    break;

  default:
    fprintf(stderr,"I cannot deal with source BITPIX = %d\n", bp_src);
    return(ERR_ILLEGAL_BITPIX);
  }

/* If output is FITSendian, do the necessary swabbing */
  if(fitswab & 2) {
/* Separate out the case involving fp -> IEEE because of VAX possibility */
    if(bp_dest == -32) fpieee(npix, dest, dest);
#ifdef LOWENDIAN
    else FITSorder(bp_dest, npix, dest);
#endif
  }
  return(0);
}

/**************************************/
/* Below are header altering routines */
/**************************************/

/*
 * newfitshead() will create a minimal FITS header, don't try to add to it!
 */
newfitshead(header,bitpix,nx,ny,object)
char **header;		/* header */
int bitpix, nx, ny;	/* BITPIX, NAXIS1, NAXIS2 header entries */
char *object;		/* object name */
{
  int i;
  *header = (char *) malloc(NFITS+1);
  for(i=0; i<NFITS; i++) (*header)[i] = ' ';
  (*header)[NFITS] = '\0';	/* Just to keep a strlen() happy */
  wfitem(0, *header, "SIMPLE  ", "T", 0, 0.0);
  wfitem(1, *header, "BITPIX  ", NULL, bitpix, 0.0);
  wfitem(2, *header, "NAXIS   ", NULL, 2, 0.0);
  wfitem(3, *header, "NAXIS1  ", NULL, nx, 0.0);
  wfitem(4, *header, "NAXIS2  ", NULL, ny, 0.0);
  wfitem(5, *header, "OBJECT  ", object, 0, 0.0);
  wfitem(6, *header, "END     ", NULL, 0, 0.0);
  return(0);
}

/*
 * addfitshead() adds a new line to a FITS header
 */
addfitshead(nhead, header, keyword, cvalue, ivalue, rvalue)
int *nhead;	/* Number of header entries */
char *header;	/* Header */
char *keyword;  /* Keyword to write */
char *cvalue;	/* (possible) character value to set (if non NULL) */
int ivalue;	/* (possible) integer value to set (if cvalue = "INTEGER") */
double rvalue;	/* (possible) real value to set (if cvalue = "FLOAT") */
{
  int i;
  if(strncmp(&header[80*(*nhead-1)], "END     ", 8) != 0) {
    /*
    fprintf(stderr,"addfitshead: header count is wrong; adjusting count\n");
    */
    for(i=0; i<strlen(header)/80; i++) {
      if(strncmp(&header[80*i], "END     ", 8) == 0) break;
    }
    *nhead = i + 1;
  }
  wfitem(*nhead-1, header, keyword, cvalue, ivalue, rvalue);
  wfitem(*nhead,   header, "END     ", NULL, 0, 0.0);
  *nhead += 1;
  return(0);
}

/*
 * chfitshead() changes a line in a FITS header
 */
chfitshead(n, header, keyword, cvalue, ivalue, rvalue)
int *n;		/* Line where keyword was found in header (nhead if not found) */
char *header;	/* Header */
char *keyword;  /* Keyword to write */
char *cvalue;	/* (possible) character value to set (if non NULL) */
int ivalue;	/* (possible) integer value to set (if cvalue = "INTEGER") */
double rvalue;	/* (possible) real value to set (if cvalue = "FLOAT") */
{
  int i;
  for(i=0; i<strlen(header)/80; i++) {
    if(strncmp(&header[80*i], keyword, 8) == 0) break;
/* Keyword never found in header? */
    if(strncmp(&header[80*i], "END     ", 8) == 0) {
      *n = i + 1;
      return(1);
    }
  }

  wfitem(i, header, keyword, cvalue, ivalue, rvalue);
  *n = i;
  return(0);
}

/*
 * ifitshead() reads an integer component of a FITS header
 */
ifitshead(header, keyword, ivalue)
char *header;	/* Header */
char *keyword;  /* Keyword to write */
int *ivalue;	/* integer value retrieved */
{
  int i;
  for(i=0; i<strlen(header)/80; i++) {
    if(strncmp(&header[80*i], keyword, 8) == 0) break;
/* Keyword never found in header? */
    if(strncmp(&header[80*i], "END     ", 8) == 0) return(1);
  }
  *ivalue = atoi(&header[80*i+9]);
  return(0);
}

/*
 * rfitshead() reads a real component of a FITS header
 */
rfitshead(header, keyword, value)
char *header;	/* Header */
char *keyword;  /* Keyword to write */
double *value;	/* integer value retrieved */
{
  int i;
  double atof();
  for(i=0; i<strlen(header)/80; i++) {
    if(strncmp(&header[80*i], keyword, 8) == 0) break;
/* Keyword never found in header? */
    if(strncmp(&header[80*i], "END     ", 8) == 0) return(1);
  }
  *value = atof(&header[80*i+9]);
  return(0);
}

/*
 * wfitem() writes a specific line to a FITS header
 */
wfitem(n, header, keyword, cvalue, ivalue, rvalue)
int n;		/* Index of header entry to overwrite */
char *header;	/* Header */
char *keyword;  /* Keyword to write */
char *cvalue;	/* (possible) character value to set (if non NULL) */
int ivalue;	/* (possible) integer value to set (if cvalue = "INTEGER") */
double rvalue;	/* (possible) real value to set (if cvalue = "FLOAT") */
{
  int i;
  if(strncmp(keyword, "END     ", 8) == 0) {
    i = sprintf(&header[80*n], "%-8s", keyword);
  } else if(strncmp(keyword, "SIMPLE  ", 8) == 0) {
    i = sprintf(&header[80*n], "%-8s= %20s", keyword,cvalue);
  } else {
    if(cvalue == NULL || (strncmp(cvalue, "INTEGER", 7)==0)) {
      i = sprintf(&header[80*n], "%-8s= %20d", keyword, ivalue);
    } else if(strncmp(cvalue, "FLOAT", 5) == 0) {
      i = sprintf(&header[80*n], "%-8s= %20.8e", keyword, rvalue);
    } else {
      i = sprintf(&header[80*n], "%-8s= %20s", keyword, cvalue);
    }
  }
  if (i < 0) {
    fprintf(stderr,"wfitem: Error writing header string '%s'!\n", keyword);
    return(1);
  }
  for(  ; i<80; i++) {    /* overwrite any lurking \0's */
    if (header[80*n+i] == '\0') header[80*n+i] = ' ';
  }
  return(0);
}

/*
 * counthead() returns the number of lines in a FITS header
 */
counthead(header)
char *header;
{
  int i;
  for(i=0; i<strlen(header); i+=80) {
    if(strncmp(&header[i], "END     ", 8) == 0) break;
  }
  return(i/80+1);
}




/*******************************************/
/* Below are the basic read/write routines */
/*******************************************/

/* Read just a little bit to see if the essential keywords are present */
testfits(file, nhead, bitpix, naxis, dim, scale, zero)
char *file;
int *naxis, *dim, *nhead, *bitpix;
double *scale, *zero;
{
  int fd, i, nread, mode=0, n, nmax=0, setscale=0, setzero=0;
  double atof();
  *bitpix = *naxis = 0;

  if( (fd = open(file, 0)) == -1) return(ERR_CANT_OPEN_FILE);
  if((nread=read(fd,buf,80)) != 80) {
/*    fprintf(stderr,"testfits: error sniffing header\n"); */
    close(fd);
    return(ERR_CANT_READ_HEADER);
  }
  if(strncmp(buf,"SIMPLE  ",8) != 0) return(0);

  *nhead = 0;

/* Now look for NAXIS, NAXISn, and END */
  for(i=1; i<10000; i++) {
    if((nread=read(fd, buf, 80)) != 80) {
      close(fd);
      return(ERR_CANT_READ_HEADER);
    }
/*    printf("%3d %.20s\n", i, buf); */
    if(strncmp(buf, "END     ", 8) == 0) {
      *nhead = i + 1;
      break;
    } else if(strncmp(buf, "BITPIX  ", 8) == 0) {
      *bitpix = atoi(buf+10);
    } else if(strncmp(buf, "NAXIS   ", 8) == 0) {
      *naxis = atoi(buf+10);
    } else if( (strncmp(buf, "NAXIS", 5) == 0) &&
	     (strncmp(buf+6, "  ", 2) == 0) &&
	     (sscanf(buf+5, "%1d", &n) == 1) && n <= *naxis && n > 0) {
      dim[n-1] = atoi(buf+10);
      nmax = MAX(n, nmax);
    } else if( strncmp(buf, "BSCALE  ", 8) == 0) {
      *scale = atof(buf+10);
      setscale += 1;
    } else if( strncmp(buf, "BZERO   ", 8) == 0)  {
      *zero = atof(buf+10);
      setzero += 1;
    }
  }
  close(fd);
  if(*nhead == 0) return(ERR_NO_END);
  if(*naxis == 0 || nmax != *naxis) return(ERR_NO_NAXIS);
  if(*bitpix == 0) return(ERR_NOT_A_FITS_FILE);
  if( setscale != 1 ) *scale = 1.0;
  if( setzero != 1) *zero = 0.0;
  return(0);
}

/* Allocate space for a header and data */
allocfits(nhead, head, naxis, dim, bitpix, data)
int nhead;		/* Number of header entries to allocate */
char **head;		/* Pointer to header pointer */
int naxis;		/* Number of axes */
int *dim;		/* Array of axis dimensions */
int bitpix;		/* BITPIX of eventual desired data format */
char **data;		/* Pointer to data array */
{
  int nalloc;
  *head = (char *)malloc(80*nhead + 1);
  (*head)[80*nhead] = '\0';		/* Terminate in case of strlen() */
  nalloc = 1;  
  while(--naxis >= 0) nalloc *= dim[naxis];
  nalloc = (2*((ABS(bitpix)*nalloc+15)/16)+sizeof(float)-1) / sizeof(float);
/* calloc floats for alignment, and toss in one extra for paranoia */
  *data = (char *)calloc(nalloc+1, sizeof(float));
  return(0);
}

/* Read just the file's header */
rhead(fd, maxhead, head)
int fd, maxhead;
char *head;
{
  int i, nread, extra, end;

  end = 0;
  for(i=0; i<maxhead; i+=80) {
    if((nread=read(fd,head+i,80)) != 80) {
      fprintf(stderr,"rfhead: error reading header\n"); 
      return(ERR_CANT_READ_HEADER); }
    if(strncmp(&head[i],"END     ",8) == 0) {
      end = 1; break;
    }
  }
  if(end == 0) {
    fprintf(stderr,"rfhead: insufficient buffer for header\n"); 
    return(ERR_INSUFFICIENT_HEADER); }

  i+= 80;
  extra = NFITS - (i%NFITS);
  if((i%NFITS) !=0) {
    if((nread = read(fd,buf,extra)) != extra) {
      fprintf(stderr,"rfhead: error reading header\n"); 
      return(ERR_CANT_READ_HEADER); }
  }
/* Paranoia in case there are nulls in the header */
  while(--i >= 0) {
    if(head[i] == '\0') head[i] = ' ';
  }
  return(0);
}

/* Write just the file's header */
whead(fd,nhead,head)
int fd, nhead;
char *head;
{
  int i, nwrite, extra, end;

  end = 0;
  for(i=0; i<nhead*80; i+=80) {

    if((nwrite=write(fd,head+i,80)) != 80) {
      fprintf(stderr,"wfhead: error writing header\n"); 
      return(ERR_CANT_WRITE_HEADER); }
    if(strncmp(&head[i],"END     ",8) == 0) {
      end = 1; break;
    }
  }

  if(end == 0) {
    fprintf(stderr,"wfhead: END missing from header\n"); 
    return(ERR_NO_END); }
  i+= 80;

  extra = NFITS - (i%NFITS);
  if((i%NFITS) !=0) {
    for(i=0;i<extra;i++) buf[i] = ' ';
    if((nwrite = write(fd,buf,extra)) != extra) { 
      fprintf(stderr,"wfhead: error writing header\n"); 
      return(ERR_CANT_WRITE_HEADER); }
  }
  return(0);
}

parsehead(maxhead,nhead,head,bitpix,naxis,dim,scale,zero)
int maxhead;
int *nhead, *bitpix, *naxis, *dim;
char *head;
double *scale, *zero;
{
  int i, setscale=0, setzero=0, n, nmax=0;
  double atof();
  *bitpix = *naxis = *nhead = 0;

/* Now look for NAXIS, NAXISn, and END */
  for(i=0; i<maxhead; i+=80) {
    if(strncmp(head+i, "END     ", 8) == 0) {
      *nhead = i / 80 + 1;
      break;
    }
    if(strncmp(head+i, "BITPIX  ", 8) == 0) *bitpix = atoi(head+i+10);
    if(strncmp(head+i, "NAXIS   ", 8) == 0) *naxis = atoi(head+i+10);
    if( (strncmp(head+i, "NAXIS", 5) == 0) &&
	     (strncmp(head+i+6, "  ", 2) == 0) &&
	     (sscanf(head+i+5, "%1d", &n) == 1) && n <= *naxis && n > 0) {
      dim[n-1] = atoi(head+i+10);
      nmax = MAX(n, nmax);
    }
    if( strncmp(head+i, "BSCALE  ", 8) == 0) {
      *scale = atof(head+i+10);
      setscale += 1; }
    if( strncmp(head+i, "BZERO   ", 8) == 0)  {
      *zero = atof(head+i+10);
      setzero += 1; }
  }
  if(*nhead == 0) return(ERR_NO_END);
  if(*naxis == 0 || nmax != *naxis) return(ERR_NO_NAXIS);
  if( setscale != 1 ) *scale = 1.0;
  if( setzero != 1) *zero = 0.0;
  return(0);
}

/**********************************/
/* Basic type conversion routines */
/**********************************/

#define NINT(x) ( ((x) < 0.0) ? (int)((x)-0.5) : (int)((x)+0.5) )
#define TRUNCATE(low,high,x) ( ((x) > (low)) ? (((x) < (high)) ? (x) : (high)) : (low) )

/*
 * Various format conversion routines
 * Routines converting to FP (or bit to short) work from the end of 
 * the array so the input and output arrays may be identical.
 */

/* Convert (unsigned) short integers to floating point */
UFconvert(npix, scale, zero, src, dest)
int npix;
unsigned short int *src;
float *dest;
double scale, zero;
{
  int i;
  float tol=0.00001*scale;
  if(scale == 1.0) {
    if(zero == 0.0) {
      for(i=(npix-1);i>=0;i--) dest[i] = src[i];
    } else {
      for(i=(npix-1);i>=0;i--) dest[i] = (src[i]+zero);
    }
  } else {
    for(i=(npix-1);i>=0;i--) {
      dest[i] = zero + scale*src[i];
/* If the value is REALLY close to zero, it probably is supposed to be zero */
      if(dest[i] > -tol && dest[i] < tol) dest[i] = 0.0;
    }
  }
}

/* Convert short integers to floating point */
SFconvert(npix, scale, zero, src, dest)
int npix;
short int *src;
float *dest;
double scale, zero;
{
  int i;
  float tol=0.00001*scale;
  if(scale == 1.0) {
    if(zero == 0.0) {
      for(i=(npix-1);i>=0;i--) dest[i] = src[i];
    } else {
      for(i=(npix-1);i>=0;i--) dest[i] = src[i] + zero;
    }
  } else {
    for(i=(npix-1);i>=0;i--) {
      dest[i] = zero + scale*src[i];
/* If the value is REALLY close to zero, it probably is supposed to be zero */
      if(dest[i] > -tol && dest[i] < tol) dest[i] = 0.0;
    }
  }
}

/* Convert floating point to short integers */
FSconvert(npix,scale,zero,src,dest)
int npix;
short int *dest;
float *src;
double scale, zero;
{
  int i;
  double temp;
  if(scale == 1.0 && zero == 0.0) {
    for(i=0;i<npix;i++) {
      temp = *src++;
      temp = TRUNCATE(-32768,32767,temp);
      *dest++ = NINT(temp);
    }
  } else {
    for(i=0;i<npix;i++) {
      temp = zero + scale*(*src++);
      temp = TRUNCATE(-32768,32767,temp);
      *dest++ = NINT(temp);
    }
  }
}

/* Convert floating point to unsigned short integers */
FUconvert(npix,scale,zero,src,dest)
int npix;
float *src;
unsigned short int *dest;
double scale, zero;
{
  int i;
  double temp;
  if(scale == 1.0 && zero == 0.0) {
    for(i=0; i<npix; i++) {
      temp = *src++;
      temp = TRUNCATE(0,65535,temp);
      *dest++ = NINT(temp);
      /*
      printf("%10d %12g %12g %12d\n", i, *(src-1), temp, *(dest-1));
      */
    }
  } else {
    for(i=0; i<npix; i++) {
      temp = zero + scale*(*src++);
      temp = TRUNCATE(0,65535,temp);
      *dest++ = NINT(temp);
    }
  }
}

/* Convert short ints to unsigned short integers */
SUconvert(npix,scale,zero,src,dest)
int npix;
short *src;
unsigned short int *dest;
double scale, zero;
{
  int i, itemp, izero;
  double temp;
  if(scale == 1.0) {
    if(zero == 0.0) {
      for(i=0; i<npix; i++) {
	itemp = *src++;
	*dest++ = MAX(0, itemp);
      }
    } else {
      izero = NINT(zero);
      for(i=0; i<npix; i++) {
	itemp = izero + *src++;
	*dest++ = TRUNCATE(0, 65535, itemp);
      }
    }
  } else {
    for(i=0; i<npix; i++) {
      temp = zero + scale*(*src++);
      temp = TRUNCATE(0, 65535, temp);
      *dest++ = NINT(temp);
    }
  }
}

/* Convert unsigned short ints to short integers */
USconvert(npix,scale,zero,src,dest)
int npix;
unsigned short *src;
short int *dest;
double scale, zero;
{
  int i, itemp, izero;
  double temp;
  if(scale == 1.0) {
    if(zero == 0.0) {
      for(i=0; i<npix; i++) {
	itemp = *src++;
	*dest++ = MIN(32767, itemp);
      }
    } else {
      izero = NINT(zero);
      for(i=0; i<npix; i++) {
	itemp = izero + *src++;
	*dest++ = TRUNCATE(-32768, 32767, itemp);
      }
    }
  } else {
    for(i=0; i<npix; i++) {
      temp = zero + scale*(*src++);
      temp = TRUNCATE(-32768, 32767, temp);
      *dest++ = NINT(temp);
    }
  }
}

/* Convert unsigned short ints to unsigned short integers */
UUconvert(npix,scale,zero,src,dest)
int npix;
unsigned short int *dest;
unsigned short *src;
double scale, zero;
{
  int i, itemp, izero;
  double temp;
  if(scale == 1.0) {
    if(zero == 0.0) {
      bcopy(src, dest, sizeof(unsigned short int)*npix);
    } else {
      izero = NINT(zero);
      for(i=0; i<npix; i++) {
	itemp = izero + *src++;
	*dest++ = TRUNCATE(0, 65535, itemp);
      }
    }
  } else {
    for(i=0; i<npix; i++) {
      temp = zero + scale*(*src++);
      temp = TRUNCATE(0, 65535, temp);
      *dest++ = NINT(temp);
    }
  }
}

/* Convert short ints to short integers */
SSconvert(npix,scale,zero,src,dest)
int npix;
short int *dest;
short *src;
double scale, zero;
{
  int i, itemp, izero;
  double temp;
  if(scale == 1.0) {
    if(zero == 0.0) {
      bcopy(src, dest, sizeof(short int)*npix);
    } else {
      izero = NINT(zero);
      for(i=0; i<npix; i++) {
	itemp = izero + *src++;
	*dest++ = TRUNCATE(-32768, 32767, itemp);
      }
    }
  } else {
    for(i=0; i<npix; i++) {
      temp = zero + scale*(*src++);
      temp = TRUNCATE(-32768, 32767, temp);
      *dest++ = NINT(temp);
    }
  }
}

/* Convert long integers to floating point */
IFconvert(npix, scale, zero, src, dest)
int npix;
int *src;
float *dest;
double scale, zero;
{
  int i;
  float tol=1e-10*scale;
  if(scale == 1.0) {
    if(zero == 0.0) {
      for(i=0; i<npix; i++) dest[i] = src[i];
    } else {
      for(i=0; i<npix; i++) dest[i] = src[i] + zero;
    }
  } else {
    for(i=0; i<npix; i++) {
      dest[i] = zero + scale*src[i];
/* If the value is REALLY close to zero, it probably is supposed to be zero */
      if(dest[i] > -tol && dest[i] < tol) dest[i] = 0.0;
    }
  }
}

/* Convert floating point to long integers */
FIconvert(npix,scale,zero,src,dest)
int npix;
int *dest;
float *src;
double scale, zero;
{
  int i;
  double temp;
  if(scale == 1.0 && zero == 0.0) {
    for(i=0;i<npix;i++) {
      temp = *src++;
      temp = TRUNCATE(-2147483648.0, 2147483647.0, temp);
      *dest++ = NINT(temp);
    }
  } else {
    for(i=0;i<npix;i++) {
      temp = zero + scale*(*src++);
      temp = TRUNCATE(-2147483648.0, 2147483647.0, temp);
      *dest++ = NINT(temp);
    }
  }
}

/* Convert bitmap to short integers */
BSconvert(npix,scale,zero,src,dest)
int npix;
unsigned char *src;
short *dest;
double scale, zero;
{
  int i;
  unsigned char b;
  if( ((npix-1)%16) != 15 ) b = src[(npix-1)/16];
  for(i=(npix-1); i>=0; i--) {
    if( (i%16) == 15 ) b = src[i/16];
    dest[i] = (b & 0x8000) ? 1 : 0;
    b = b << 1;
  }
}

/* Convert bitmap to floats */
BFconvert(npix,scale,zero,src,dest)
int npix;
unsigned short *src;
float *dest;
double scale, zero;
{
  int i;
  unsigned short b;
  if( ((npix-1)%16) != 15 ) b = src[(npix-1)/16];
  for(i=(npix-1); i>=0; i--) {
    if( (i%16) == 15 ) b = src[i/16];
    dest[i] = (b & 0x8000) ? 1.0 : 0.0;
    b = b << 1;
  }
}

/* Convert short integers to a bitmap */
SBconvert(npix,scale,zero,src,dest)
int npix;
short *src;
unsigned short *dest;
double scale, zero;
{
  unsigned short *dp=dest;
  int i;
  *dp = 0;
  for(i=0; i<npix; i++) {
    *dp = (*dp >> 1) | ((*src++ == 0) ? 0x0000 : 0x8000);
    if( (i%16) == 15 ) dp++;
  }
}

/* Convert floats to a bitmap */
FBconvert(npix,scale,zero,src,dest)
int npix;
float *src;
unsigned short *dest;
double scale, zero;
{
  unsigned short *dp=dest;
  int i;
  *dp = 0;
  for(i=0; i<npix; i++) {
    *dp = (*dp >> 1) | ((*src++ == 0.0) ? 0x0000 : 0x8000);
    if( (i%16) == 15 ) dp++;
  }
}

/* Convert IEEE format floating point to floating point */
ieeefp(npix,data,fp)
int npix;
char *data;
char *fp;		/* Really FP, but may need to access bytes */
{
  int i;
#ifdef VAXFP
  register char *dp=data, temp;
  for(i=0;i<npix;i++) {
    temp  = (*dp == 0) ? 0 : (*dp+1) ;
    *fp++ = *(dp+1);
    *fp++ = temp;
    temp  = *(dp+2);
    *fp++ = *(dp+3);
    *fp++ = temp;
    dp += 4;
  }
#else
#ifdef LOWENDIAN
  FITSorder(32,npix,data);
#endif
  if(fp != data) bcopy(data, fp, sizeof(float)*npix);
#endif
}

/* Convert floating point to IEEE format floating point */
fpieee(npix,fp,data)
int npix;
char *data;
char *fp;		/* Really FP, but may need to access bytes */
{
  int i;
#ifdef VAXFP
  register char *dp=data, temp;
  for(i=0;i<npix;i++) {
    temp  = *(fp);
    *dp++ = (*(fp+1) == 0) ? 0 : (*(fp+1)-1);
    *dp++ = temp;
    temp  = *(fp+2);
    *dp++ = *(fp+3);
    *dp++ = temp;
    fp += 4;
  }
#else
  if(fp != data) bcopy(fp, data, sizeof(float)*npix);
#ifdef LOWENDIAN
  FITSorder(32,npix,data);
#endif
#endif
}

#ifdef OSF_ALPHA
swab(src, dest, n)
char *src, *dest;
int n;
{
  register char temp, *s=src, *d=dest;
  while(n > 1) {
    temp = *s++;
    *d++ = *s++;
    *d++ = temp;
    n -= 2;
  }
}
#endif

FITSorder(bitpix,npix,data)
int bitpix, npix;
short *data;
{
  register short temp, *dp=data;
  register int i=npix;

/* Swap bytes first; work on an even number of bytes */
  swab(data, data, 2*((npix*ABS(bitpix)+15)/16));
  if(ABS(bitpix) <= 16) return(0);

/* Swap words to complete the 32 bit byte swap */
  while(i--) {
    temp = *dp;
    *dp = *(dp+1);    dp++;
    *dp++ = temp;
  }
  return(0);
}
