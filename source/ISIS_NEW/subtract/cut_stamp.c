#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

void cut_stamp(int xi, int yi, DATA_TYPE *image,  DATA_TYPE *imc, double *sum)
{
  int       i,j,xij;
 double     q;


   
  *sum=0.0;
 

  for(i=xi-stamp_size/2;i<=xi+stamp_size/2;i++)
  {
    for(j=yi-stamp_size/2;j<=yi+stamp_size/2;j++)
    {
      xij=i-xi+stamp_size/2+stamp_size*(j-yi+stamp_size/2);
      q=image[i+width*j];
      imc[xij] = q;
      *sum += q;
    }
  }

   


 return;

}
