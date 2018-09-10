#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "short.h"


void write_image2(DATA_TYPE *image, FILE *file,int cx, int cy)
{
 int   i,j,offset;
 float pf,*pfloat;
 char  cp[8],*cp2;

 offset=(cx+cy*width0)*4+2880;
 printf("Offset: %i %i %i %i %i\n", offset, cx, cy, bitpix, (int)fabs((double)bitpix)/8);

 fseek(file,offset,SEEK_SET);

 offset=(width0-width)*4;

 for (j=0; j<height; j++)
 {
  if(!swap_flag)
  {
   for(i=0;i<width;i++)
   {
    pf = (float)image[i+width*j];
    fwrite(&pf,4,1,file);
   }
   fseek(file,offset,SEEK_CUR);  
  }
  else
  {
   for(i=0;i<width;i++)
   {
    pf = (float)image[i+width*j];
    cp2 = (char *)&pf;
    cp[0]=cp2[3]; cp[1]=cp2[2]; cp[2]=cp2[1]; cp[3]=cp2[0];
    pfloat=(float *)&cp[0];
    fwrite(pfloat,4,1,file);
   }
   fseek(file,offset,SEEK_CUR);  
  }
 }

 return;
}
