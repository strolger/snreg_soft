#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"


void phot(DATA_TYPE *image,DATA_TYPE *image2,DATA_TYPE *Bpsf,int xi,int yi,char *fname)
{
  int    i,j,ip,jp;
  double a,b,c,ps,im,weight,a2,b2,x,y,r,err1,err2;
  FILE   *outfile;


 
 a=b=c=a2=b2=err1=err2=0.0;
 for(i=xi-conv_size/2;i<=xi+conv_size/2;i++)
  for(j=yi-conv_size/2;j<=yi+conv_size/2;j++)
  {
   ip = i-xi+conv_size/2;
   jp = j-yi+conv_size/2;
   ps = Bpsf[ip+jp*conv_size];
   im = image[i+image_width*j];
   weight = image2[i+image_width*j];
   x=(double)(i-xi);
   y=(double)(j-yi);
   r=sqrt(x*x+y*y);

   if(r <= radphot)
   {
    a += ps*im;
    b += ps*ps;
    c += im*im;
    a2 += ps*im/weight;
    b2 += ps*ps/weight;
    err1 += ps*ps*weight;
    err2 += ps*ps/weight;
   }

  }
  err1=sqrt(err1/(double)adu_el)/b;
  err2=sqrt(err2/(double)adu_el)/b2;

 printf("file: %s\n", fname);

 outfile=fopen(fname,"a+");
 fprintf(outfile,"%15.6lf %12.4lf %12.4lf %12.4lf %12.4lf %9.4lf\n",exp_date, a2/b2,err2,a/b,err1,fabs(a)/sqrt(b*c));
 fclose(outfile);
 return;
}
