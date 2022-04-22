#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"

char max(int x, int y, int mesh, DATA_TYPE *image, DATA_TYPE *sum)
{
  int           i,j;
  DATA_TYPE     pix,pix2;
 
 
 
  pix=image[x+width*y];
  *sum=0.0;
  for(i=x-mesh;i<=x+mesh;i++)
   for(j=y-mesh;j<=y+mesh;j++)
   {
    pix2=image[i+width*j];
    *sum += pix2;
    if(pix2>pix || pix2>SATURATION || pix2<PIX_MIN) return 0;
   }
 
   

 return 1;
}
