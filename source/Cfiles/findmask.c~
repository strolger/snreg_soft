#include<stdio.h>
#include<math.h>
#define NFITS (2880)
#define SWAB 0		/* Swap bytes?  Set to 1 for DecStations*/
#define HEADLINES 800

static char buf[NFITS];

main (argc,argv)
     int argc;
     char *argv[];
{
  FILE *stream1,*stream2;
  short *svector(),*image,*imagetwo,*image2,*image3; 
  int nhead=HEADLINES,nx,ny,i,j,oblitsize;
  int nhead1=HEADLINES,nx1,ny1,oblitvalue,saturated;
  int nhead2=HEADLINES,X,Y;
  float s,s1,sigma,m1,b1,tmin,tmax,mean;
  double x,y,xx,xy,m,b,temp,temp1,temp2,temp3,n;
  char file3[40],file1[40],file2[40],header[HEADLINES][80],header1[HEADLINES][80],header2[HEADLINES][80],file4[40];  
  void wffile(), rffile(),oblit(),rmback();
  
  
  if (argc == 6) {
    sprintf(file1,argv[1]);
    sprintf(file2,argv[2]);
    sscanf(argv[3],"%d",&saturated);
    sscanf(argv[4],"%d",&oblitvalue);
    sscanf(argv[5],"%d",&oblitsize);
    
    image=svector(0,5000*5000);
      
    strcat(file2,".mask");

    stream2 = fopen(file2,"w");
    rffile(&nhead,header,&nx,&ny,image,file1,HEADLINES,40); /*read in images*/
    for (i=0; i<nx*ny; i++) { /*obliterate of saturated regions */
      if ((image[i] > saturated) && image[i] !=32000 ) {
   	Y = (int)(i/nx);
    	X = (i-Y*nx);
	fprintf(stream2,"%d  %d\n",X+1,Y+1);
      }
    }
  }
  else {
    printf("\n findmask: in1 fieldname saturated_value oblit_value oblit_size\n Writes Mask file masking out pixels greater than \n saturated_value and replaces them with oblit_value.\n Images must be datatype SHORT (16)\n");
    exit();
  }
}



/* Read a short integer file's header and data */

void rffile(n,head,npix,nline,data,file,headlen,filelen)
     int *n, *npix, *nline, headlen, filelen;
     short *data;
     char *head, *file;
{
  rfitsns(16,n,head,npix,nline,data,file,headlen,filelen);
  if(SWAB) swab(data,data,2*(*npix)*(*nline));
}

/* Read a 32 bit file's header and data */

void rffint(n,head,npix,nline,data,file,headlen,filelen)
     int *n, *npix, *nline, headlen, filelen;
     unsigned short *data;
     char *head, *file;
{
  register int i;
  unsigned short temp;
  rfitsns(32,n,head,npix,nline,data,file,headlen,filelen);
  if(SWAB) {
    swab(data,data,4*(*npix)*(*nline));
    for(i=0;i<2*(*npix)*(*nline);i+=2) {
      temp = data[i];
      data[i] = data[i+1];
      data[i+1] = temp;
    }
  }
}

/* Read a bitmap file's header and data */

void rfbitfile(n,head,npix,nline,data,file,headlen,filelen)
     int *n, *npix, *nline, headlen, filelen;
     short *data;
     char *head, *file;
{
  rfitsns(1,n,head,npix,nline,data,file,headlen,filelen);
  if(SWAB) swab(data,data,2*(*npix)*(*nline));
}

/* Read a file's header and data ... no byte swap version */

void rffilens(n,head,npix,nline,data,file,headlen,filelen)
     int *n, *npix, *nline, headlen, filelen;
     short *data;
     char *head, *file;
{
  rfitsns(16,n,head,npix,nline,data,file,headlen,filelen);
}


/* Read a file of int's don't swap bytes */

void rffintns(n,head,npix,nline,data,file,headlen,filelen)
     int *n, *npix, *nline, headlen, filelen;
     int *data;
     char *head, *file;
{
  rfitsns(32,n,head,npix,nline,data,file,headlen,filelen);
}

/* Just read the bits, with bitsper bits per data element */
rfitsns(bitsper,n,head,npix,nline,data,file,headlen,filelen)
     int bitsper;
     int *n, *npix, *nline, headlen, filelen;
     unsigned char *data;
     char *head, *file;
{
  int fd, err, need, i, nread;

  while(--filelen > 0 && (file[filelen]==' ' || file[filelen]==0));
  for(i=0;i<=filelen;i++) buf[i] = file[i] ;
  buf[filelen+1] = '\0';

  if((fd=open(buf,0)) == -1) {
    fprintf(stderr,"rffile: cannot open file %s\n",buf);
    (*npix)=(*nline)=0;
    return;
  }

/* Read the header, looking for NAXIS1, NAXIS2 and END */
  if((err = getheader(fd,n,head,npix,nline)) != 0) {
    close(fd);
    (*npix)=(*nline)=0;
    return;
  }

/* Read the data according to the values for NPIX, NLINE */
  need = (bitsper*(*npix)*(*nline)+7) / 8;
  if((nread=read(fd,data,need)) != need) {
    fprintf(stderr,"rffile: end of file before data read\n");
    (*npix)=(*nline)=0;
  }
  close(fd);
}


/* Read just the file's header */

void rfhead(n,head,file,headlen,filelen)
     int *n, headlen, filelen;
     char *head, *file;
{
  int fd, err, nread, npix, nline, i;

  while(--filelen > 0 && (file[filelen]==' ' || file[filelen]==0));
  for(i=0;i<=filelen;i++) buf[i] = file[i] ;
  buf[filelen+1] = '\0';

  if((fd=open(buf,0)) == -1) {
    fprintf(stderr,"rffile: cannot open file %s\n",buf);
    return;
  }

/* Read the header, looking for NAXIS1, NAXIS2 and END */
  if((err = getheader(fd,n,head,&npix,&nline)) != 0) {
    close(fd);
    return;
  }
  close(fd);
}

/* Write a short integer fits format file */

void wffile(n,head,npix,nline,data,file,headlen,filelen)
     int *n, *npix, *nline, headlen, filelen;
     short *data;
     char *head, *file;
{
  if(SWAB) swab(data,data,2*(*npix)*(*nline));
  wfitsns(16,n,head,npix,nline,data,file,headlen,filelen);
  if(SWAB) swab(data,data,2*(*npix)*(*nline));
}

/* Write a 32 bit file's header and data */

void wffint(n,head,npix,nline,data,file,headlen,filelen)
     int *n, *npix, *nline, headlen, filelen;
     unsigned short *data;
     char *head, *file;
{
  register int i;
  unsigned short temp;
  if(SWAB) {
    swab(data,data,4*(*npix)*(*nline));
    for(i=0;i<2*(*npix)*(*nline);i+=2) {
      temp = data[i];
      data[i] = data[i+1];
      data[i+1] = temp;
    }
  }
  wfitsns(32,n,head,npix,nline,data,file,headlen,filelen);
  if(SWAB) {
    for(i=0;i<2*(*npix)*(*nline);i+=2) {
      temp = data[i];
      data[i] = data[i+1];
      data[i+1] = temp;
    }
    swab(data,data,4*(*npix)*(*nline));
  }
}

/* Write a bitmap fits format file */

void wfbitfile(n,head,npix,nline,data,file,headlen,filelen)
     int *n, *npix, *nline, headlen, filelen;
     short *data;
     char *head, *file;
{
  if(SWAB) swab(data,data,2*((*npix)*(*nline)+15)/16);
  wfitsns(1,n,head,npix,nline,data,file,headlen,filelen);
  if(SWAB) swab(data,data,2*((*npix)*(*nline)+15)/16);
}


/* Write a fits format file ... no byte swap version */

void wffilens(n,head,npix,nline,data,file,headlen,filelen)
     int *n, *npix, *nline, headlen, filelen;
     short *data;
     char *head, *file;
{
  wfitsns(16,n,head,npix,nline,data,file,headlen,filelen);
}

/* Write a fits format file ... no byte swap version */

void wffintns(n,head,npix,nline,data,file,headlen,filelen)
     int *n, *npix, *nline, headlen, filelen;
     short *data;
     char *head, *file;
{
  wfitsns(32,n,head,npix,nline,data,file,headlen,filelen);
}

wfitsns(bitsper,n,head,npix,nline,data,file,headlen,filelen)
     int bitsper;
     int *n, *npix, *nline, headlen, filelen;
     short *data;
     char *head, *file;
{
  int fd, end, nwrite, need;
  int naxis1, naxis2;
  int i, extra;

/* Consistency check, npix = NAXIS1 and nline = NAXIS2 or else! */
  naxis1 = 0;
  naxis2 = 0;
  for(i=0;i<(*n) && (naxis1==0 || naxis2==0);i++) {
    if(strncmp(head+80*i,"NAXIS1",6) == 0) naxis1 = atoi(head+80*i+10);
    if(strncmp(head+80*i,"NAXIS2",6) == 0) naxis2 = atoi(head+80*i+10);
  }
  if(naxis1==0 || naxis2==0) {
    fprintf(stderr,"wffile: cannot find NAXIS1, NAXIS2\n");
    return;
  }   
  if(naxis1 != *npix || naxis2 != *nline) {
    fprintf(stderr,"wffile: npix, NAXIS1, nline, NAXIS2 = %d %d %d %d\n",
	    *npix, naxis1, *nline, naxis2);
    return;
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
    return;
  }

/* Write the header, quitting as soon as END is found */
  end = 0;
  for(i=0;i<(*n)*80;i+=80) {
    if((nwrite=write(fd,head+i,80)) != 80) {
      fprintf(stderr,"wffile: error writing header\n");
      close(fd);
      return;
    }
    if(head[i]=='E' && head[i+1]=='N' && head[i+2]=='D') {
      end = 1;
      break;
    }
  }
  if(end == 0) {
    fprintf(stderr,"wffile: END missing from header\n");
    close(fd);
    return;
  }

  i += 80;
  extra = NFITS - (i%NFITS);
  if(extra < NFITS) {
    for(i=0;i<extra;i++) buf[i] = ' ';
    if((nwrite = write(fd,buf,extra)) != extra) {
      fprintf(stderr,"wffile: error writing header\n");
      close(fd);
      return;
    }
  }

/* Write the data to the file */

  need = (bitsper*(*npix)*(*nline)+7) / 8;
  if((nwrite = write(fd,data,need)) != need) {
    fprintf(stderr,"wffile: error writing data\n");
    close(fd);
    return;
  }

  extra = NFITS - (need%NFITS);
  if( extra < NFITS ) {
    for(i=0;i<extra;i++) buf[i] = 0;
    if((nwrite = write(fd,buf,extra)) != extra) {
      fprintf(stderr,"wffile: error writing data tail\n");
    }
  }
  close(fd);
  return(0);
}

getheader(fd,n,head,npix,nline)
     int fd, *n, *npix, *nline;
     char *head;
{
  int i, nread, extra, end,j;

  end = 0;
  for(i=0;i<(*n)*80;i+=80) {
    if((nread=read(fd,head+i,80)) != 80) {
      fprintf(stderr,"rffile: error reading header; END not found\n");
      return(1);
    }
    if(head[i]=='B' && head[i+1]=='I' && head[i+2]=='T' && head[i+3]=='P' && head[i+4]=='I' && head[i+5] == 'X') {
       
      for (j=6; j<80; j++) {
	if (head[i+j] == '-') {
	  head[i+j] = ' ';
	  head[i+j+1] = '1';
	  head[i+j+2] = '6';
	  j=80;
	}
      }
    }
    
    if(head[i]=='E' && head[i+1]=='N' && head[i+2]=='D') {
      end = 1;
      break;
    }
  }
  if(end == 0) {
    fprintf(stderr,"rffile: insufficient buffer for header\n");
    return(2);
  }

  i+= 80;
  extra = NFITS - (i%NFITS);
  if((i%NFITS) !=0) {
    if((nread = read(fd,buf,extra)) != extra) {
      fprintf(stderr,"rffile: error reading header, file too short\n");
      return(1);
    }
  }
  *n = i / 80;
  *npix = 0;
  *nline = 0;
  for(i=0;i<(*n) && (*npix==0 || *nline==0);i++) {
    if(strncmp(head+80*i,"NAXIS1",6) == 0) *npix = atoi(head+80*i+10);
    if(strncmp(head+80*i,"NAXIS2",6) == 0) *nline = atoi(head+80*i+10);
  }
  if(*npix==0 || *nline==0) {
    fprintf(stderr,"rffile: cannot find NPIX, NLINE\n");
    return(3);
  }   
  return(0);
}

/* Open a file and read a number of bytes from it */

void rawreadfile(nbyte,data,swapbyte,file,filelen)
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
    return;
  }

/* Read the data according to the requested number of bytes */
  if((nread=read(fd,data,*nbyte)) != *nbyte) {
    fprintf(stderr,"rawreadfile: end of file before data read\n");
    *nbyte = nread;
  }
  close(fd);
  if(*swapbyte == 1) swab(data,data,*nbyte);
}

void swapbytes(n,s)
     int *n;
     short *s;
{
  swab(s,s,*n);
}


short *svector(nl,nh)
int nl,nh;
{
	 short *v;

	v=( short *)malloc((unsigned) (nh-nl+1)*sizeof( short));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}
