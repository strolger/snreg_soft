/*
* Routines to read and write disk FITS format data.
* These files have a header which consists of 36*n 80 byte header records
* which are a duplicate of the FITS tape format, with a mandatory END record.
* This is followed by 1440*n short integers of data.
*/

#include <stdio.h>

#define NFITS (2880)
#define ERR_CANT_OPEN_FILE		1 
#define ERR_CANT_READ_HEADER		2 
#define ERR_CANT_READ_DATA		3 
#define ERR_CANT_WRITE_HEADER		4 
#define ERR_CANT_WRITE_DATA		5 
#define ERR_INSUFFICIENT_HEADER		6 
#define ERR_NO_NAXIS			7 
#define ERR_NO_END			8 
#define ERR_INCONSISTENT_NAXIS		9 

#define SWAB 0		/* Swap bytes? */


static char buf[NFITS];

/* Read a short integer file's header and data */

rffile_(n,head,nx,ny,data,file,headlen,filelen)
     int *n, *nx, *ny, headlen, filelen;
     short *data;
     char *head, *file;
{
  int err;
  err = rfitsns(16,n,head,nx,ny,data,file,headlen,filelen);
  if(SWAB) swab(data,data,2*(*nx)*(*ny));
  return(err);
}

/* Read a 32 bit file's header and data */

rffint_(n,head,nx,ny,data,file,headlen,filelen)
     int *n, *nx, *ny, headlen, filelen;
     unsigned short *data;
     char *head, *file;
{
  int err;
  register int i;
  unsigned short temp;
  err = rfitsns(32,n,head,nx,ny,data,file,headlen,filelen);
  if(SWAB) {
    swab(data,data,4*(*nx)*(*ny));
    for(i=0;i<2*(*nx)*(*ny);i+=2) {
      temp = data[i];
      data[i] = data[i+1];
      data[i+1] = temp;
    }
  }
  return(err);
}

/* Read a bitmap file's header and data */

rfbitfile_(n,head,nx,ny,data,file,headlen,filelen)
     int *n, *nx, *ny, headlen, filelen;
     short *data;
     char *head, *file;
{
  int err;
  err = rfitsns(1,n,head,nx,ny,data,file,headlen,filelen);
  if(SWAB) swab(data,data,((*nx)*(*ny)+7)/8);
  return(err);
}

/* Read a file's header and data ... no byte swap version */

rffilens_(n,head,nx,ny,data,file,headlen,filelen)
     int *n, *nx, *ny, headlen, filelen;
     short *data;
     char *head, *file;
{
  int err;
  err = rfitsns(16,n,head,nx,ny,data,file,headlen,filelen);
  return(err);
}


/* Read a file of int's don't swap bytes */

rffintns_(n,head,nx,ny,data,file,headlen,filelen)
     int *n, *nx, *ny, headlen, filelen;
     int *data;
     char *head, *file;
{
  int err;
  err = rfitsns(32,n,head,nx,ny,data,file,headlen,filelen);
  return(err);
}

/* Just read the bits, with bitsper bits per data element */
rfitsns(bitsper,n,head,nx,ny,data,file,headlen,filelen)
     int bitsper;
     int *n, *nx, *ny, headlen, filelen;
     unsigned char *data;
     char *head, *file;
{
  int fd, err, need, i, nread;

  while(--filelen > 0 && (file[filelen]==' ' || file[filelen]==0));
  for(i=0;i<=filelen;i++) buf[i] = file[i] ;
  buf[filelen+1] = '\0';

  if((fd=open(buf,0)) == -1) {
    fprintf(stderr,"rffile: cannot open file %s\n",buf);
    (*nx)=(*ny)=0;
    return(ERR_CANT_OPEN_FILE);
  }

/* Read the header, looking for NAXIS1, NAXIS2 and END */
  if((err = getheader(fd,n,head,nx,ny)) != 0) {
    close(fd);
    (*nx)=(*ny)=0;
    return(err);
  }

/* Read the data according to the values for NX, NY */
  need = (bitsper*(*nx)*(*ny)+7) / 8;
  if((nread=read(fd,data,need)) != need) {
    fprintf(stderr,"rffile: end of file before data read\n");
    (*nx)=(*ny)=0;
    close(fd);
    return(ERR_CANT_READ_DATA);
  }
  close(fd);
  return(0);
}


/* Read just the file's header */

rfhead_(n,head,file,headlen,filelen)
     int *n, headlen, filelen;
     char *head, *file;
{
  int fd, err, nread, nx, ny, i;

  while(--filelen > 0 && (file[filelen]==' ' || file[filelen]==0));
  for(i=0;i<=filelen;i++) buf[i] = file[i] ;
  buf[filelen+1] = '\0';

  if((fd=open(buf,0)) == -1) {
    fprintf(stderr,"rffile: cannot open file %s\n",buf);
    return(ERR_CANT_OPEN_FILE);
  }

/* Read the header, looking for NAXIS1, NAXIS2 and END */
  if((err = getheader(fd,n,head,&nx,&ny)) != 0) {
    close(fd);
    return(err);
  }
  close(fd);
  return(0);
}

/* Write a short integer fits format file */

wffile_(n,head,nx,ny,data,file,headlen,filelen)
     int *n, *nx, *ny, headlen, filelen;
     short *data;
     char *head, *file;
{
  int err;
  if(SWAB) swab(data,data,2*(*nx)*(*ny));
  err = wfitsns(16,n,head,nx,ny,data,file,headlen,filelen);
  if(SWAB) swab(data,data,2*(*nx)*(*ny));
  return(err);
}

/* Write a 32 bit file's header and data */

wffint_(n,head,nx,ny,data,file,headlen,filelen)
     int *n, *nx, *ny, headlen, filelen;
     unsigned short *data;
     char *head, *file;
{
  int err;
  register int i;
  unsigned short temp;
  if(SWAB) {
    swab(data,data,4*(*nx)*(*ny));
    for(i=0;i<2*(*nx)*(*ny);i+=2) {
      temp = data[i];
      data[i] = data[i+1];
      data[i+1] = temp;
    }
  }
  err = wfitsns(32,n,head,nx,ny,data,file,headlen,filelen);
  if(SWAB) {
    for(i=0;i<2*(*nx)*(*ny);i+=2) {
      temp = data[i];
      data[i] = data[i+1];
      data[i+1] = temp;
    }
    swab(data,data,4*(*nx)*(*ny));
  }
  return(err);
}

/* Write a bitmap fits format file */

wfbitfile_(n,head,nx,ny,data,file,headlen,filelen)
     int *n, *nx, *ny, headlen, filelen;
     short *data;
     char *head, *file;
{
  int err;
  if(SWAB) swab(data,data,2*((*nx)*(*ny)+15)/16);
  err = wfitsns(1,n,head,nx,ny,data,file,headlen,filelen);
  if(SWAB) swab(data,data,2*((*nx)*(*ny)+15)/16);
  return(err);
}


/* Write a fits format file ... no byte swap version */

wffilens_(n,head,nx,ny,data,file,headlen,filelen)
     int *n, *nx, *ny, headlen, filelen;
     short *data;
     char *head, *file;
{
  int err;
  err = wfitsns(16,n,head,nx,ny,data,file,headlen,filelen);
  return(err);
}

/* Write a fits format file ... no byte swap version */

wffintns_(n,head,nx,ny,data,file,headlen,filelen)
     int *n, *nx, *ny, headlen, filelen;
     short *data;
     char *head, *file;
{
  int err;
  err = wfitsns(32,n,head,nx,ny,data,file,headlen,filelen);
  return(err);
}

wfitsns(bitsper,n,head,nx,ny,data,file,headlen,filelen)
     int bitsper;
     int *n, *nx, *ny, headlen, filelen;
     short *data;
     char *head, *file;
{
  int fd, end, nwrite, need;
  int naxis1, naxis2;
  int i, extra;

/* Consistency check, nx = NAXIS1 and ny = NAXIS2 or else! */
  naxis1 = 0;
  naxis2 = 0;
  for(i=0;i<(*n) && (naxis1==0 || naxis2==0);i++) {
    if(strncmp(head+80*i,"NAXIS1",6) == 0) naxis1 = atoi(head+80*i+10);
    if(strncmp(head+80*i,"NAXIS2",6) == 0) naxis2 = atoi(head+80*i+10);
  }
  if(naxis1==0 || naxis2==0) {
    fprintf(stderr,"wffile: cannot find NAXIS1, NAXIS2\n");
    return(ERR_NO_NAXIS);
  }   
  if(naxis1 != *nx || naxis2 != *ny) {
    fprintf(stderr,"wffile: nx, NAXIS1, ny, NAXIS2 = %d %d %d %d\n",
	    *nx, naxis1, *ny, naxis2);
    return(ERR_INCONSISTENT_NAXIS);
  }

/* Open the file */
  while(--filelen > 0 && (file[filelen]==' ' || file[filelen]==0));
  for(i=0;i<=filelen;i++) buf[i] = file[i] ;
  buf[filelen+1] = '\0';
/*
  fprintf(stderr,"file %s filelen %d\n",buf,filelen);
*/
  if((fd=creat(buf,0644)) == -1) {
    fprintf(stderr,"wffile: cannot open file %s\n",buf);
    return(ERR_CANT_OPEN_FILE);
  }

/* Write the header, quitting as soon as END is found */
  end = 0;
  for(i=0;i<(*n)*80;i+=80) {
    if((nwrite=write(fd,head+i,80)) != 80) {
      fprintf(stderr,"wffile: error writing header\n");
      close(fd);
      return(ERR_CANT_WRITE_HEADER);
    }
    if(head[i]=='E' && head[i+1]=='N' && head[i+2]=='D') {
      end = 1;
      break;
    }
  }
  if(end == 0) {
    fprintf(stderr,"wffile: END missing from header\n");
    close(fd);
    return(ERR_NO_END);
  }

  i += 80;
  extra = NFITS - (i%NFITS);
  if(extra < NFITS) {
    for(i=0;i<extra;i++) buf[i] = ' ';
    if((nwrite = write(fd,buf,extra)) != extra) {
      fprintf(stderr,"wffile: error writing header\n");
      close(fd);
      return(ERR_CANT_WRITE_HEADER);
    }
  }

/* Write the data to the file */

  need = (bitsper*(*nx)*(*ny)+7) / 8;
  if((nwrite = write(fd,data,need)) != need) {
    fprintf(stderr,"wffile: error writing data\n");
    close(fd);
    return(ERR_CANT_WRITE_DATA);
  }

  extra = NFITS - (need%NFITS);
  if( extra < NFITS ) {
    for(i=0;i<extra;i++) buf[i] = 0;
    if((nwrite = write(fd,buf,extra)) != extra) {
      fprintf(stderr,"wffile: error writing data tail\n");
      close(fd);
      return(ERR_CANT_WRITE_DATA);
    }
  }
  close(fd);
  return(0);
}

getheader(fd,n,head,nx,ny)
     int fd, *n, *nx, *ny;
     char *head;
{
  int i, nread, extra, end;

  end = 0;
  for(i=0;i<(*n)*80;i+=80) {
    if((nread=read(fd,head+i,80)) != 80) {
      fprintf(stderr,"rffile: error reading header; END not found\n");
      return(1);
    }
    if(head[i]=='E' && head[i+1]=='N' && head[i+2]=='D') {
      end = 1;
      break;
    }
  }
  if(end == 0) {
    fprintf(stderr,"rffile: insufficient buffer for header\n");
    return(ERR_INSUFFICIENT_HEADER);
  }

  i+= 80;
  extra = NFITS - (i%NFITS);
  if((i%NFITS) !=0) {
    if((nread = read(fd,buf,extra)) != extra) {
      fprintf(stderr,"rffile: error reading header, file too short\n");
      return(ERR_CANT_READ_HEADER);
    }
  }
  *n = i / 80;
  *nx = 0;
  *ny = 0;
  for(i=0;i<(*n) && (*nx==0 || *ny==0);i++) {
    if(strncmp(head+80*i,"NAXIS1",6) == 0) *nx = atoi(head+80*i+10);
    if(strncmp(head+80*i,"NAXIS2",6) == 0) *ny = atoi(head+80*i+10);
  }
  if(*nx==0 || *ny==0) {
    fprintf(stderr,"rffile: cannot find NX, NY\n");
    return(ERR_NO_NAXIS);
  }   
  return(0);
}

/* Open a file and read a number of bytes from it */

rawreadfile_(nbyte,data,swapbyte,file,filelen)
int *nbyte, filelen, *swapbyte;
char *data, *file;
{
  int fd, nread, i;

  while(--filelen > 0 && (file[filelen]==' ' || file[filelen]==0));
  for(i=0;i<=filelen;i++) buf[i] = file[i] ;
  buf[filelen+1] = '\0';

  if((fd=open(buf,0)) == -1) {
    fprintf(stderr,"rffile: cannot open file %s\n",buf);
    *nbyte = 0;
    return(ERR_CANT_OPEN_FILE);
  }

/* Read the data according to the requested number of bytes */
  if((nread=read(fd,data,*nbyte)) != *nbyte) {
    fprintf(stderr,"rawreadfile: end of file before data read\n");
    *nbyte = nread;
    close(fd);
    return(ERR_CANT_READ_DATA);
  }
  close(fd);
  if(*swapbyte == 1) swab(data,data,*nbyte);
  return(0);
}

swapbytes_(n,s)
     int *n;
     short *s;
{
  swab(s,s,*n);
}

/* Convert short to bitmap */
shortbit_(npix,sh,data)
int *npix;
unsigned char *data;
short *sh;
{
  unsigned char *dp=data;
  int i;
  *dp = 0;
  for(i=0;i<*npix;i++) {
    *dp = (*dp >> 1) | ((*sh++ == 0) ? 0x00 : 0x80);
    if( (i%8) == 7 ) dp++;
  }
/* It seems silly always to SWAB data, even on BIGENDIAN machines, but my
 * current definition of a bitmap is in terms of 16-bit words, with LOW bit
 * first.  It's wrong, but it's the law.
 */
/*  if(SWAB)  */
  swab(data,data,2*((*npix)+15)/16);
}
