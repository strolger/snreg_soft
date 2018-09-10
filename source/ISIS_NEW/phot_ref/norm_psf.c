#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void norm_psf()
{
 int       i,j;
 DATA_TYPE r,rx;
 double    sum,x,y;
 FILE      *file;
 char      s[256];


 rx=(DATA_TYPE)stack_width/2-1.0;
 sum=0.0;

 for(i=0;i<stack_width;i++)
 {
  for(j=0;j<stack_width;j++)
  {
   x=(double)(i-stack_width/2);
   y=(double)(j-stack_width/2);
   r=sqrt(x*x+y*y);
   if(r<=rx) sum += psf[i+stack_width*j];
  }
 }

 sum = 1.0/sum; 

 for(i=0;i<stack_width*stack_width;i++) psf[i] *= sum*10000.0;
 sprintf(s,"toto");
 file=fopen(s,"wb");
 fwrite(psf,stack_width*stack_width*sizeof(DATA_TYPE),1,file);
 fclose(file);
 return;
}
