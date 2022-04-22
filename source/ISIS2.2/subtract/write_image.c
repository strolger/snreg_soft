#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "short.h"


void write_image(DATA_TYPE *image, FILE *file,int cx, int cy)
{
 int   i,j,offset;
 float pf,*pfloat;
 char  cp[8],*cp2;

 offset=(cx+(cy+mesh_size/2)*width0+mesh_size/2)*4+2880;
 printf("offset: %i %i %i %i %i\n", offset, cx, cy, bitpix, (int)fabs((double)bitpix)/8);

 fseek(file,offset,SEEK_SET);

 offset=(width0-width+mesh_size-1)*4;

 for (j=mesh_size/2; j<height-mesh_size/2; j++)
 {
  if(!swap_flag)
  {
   for(i=mesh_size/2;i<width-mesh_size/2;i++)
   {
    pf = (float)image[i+width*j];
    fwrite(&pf,4,1,file);
   }
   fseek(file,offset,SEEK_CUR);  
  }
  else
  {
   for(i=mesh_size/2;i<width-mesh_size/2;i++)
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
