#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"



void smooth(DATA_TYPE *image, int n, int mesh)
{
 int i,j,ii,jj;
 DATA_TYPE norm;

 for(i=mesh;i<width-mesh;i++)
  for(j=mesh;j<n-mesh;j++) temp_smooth[i+width*j]=0.0;
 

 for(i=mesh;i<width-mesh;i++)
  for(j=mesh;j<n-mesh;j++)
  {
   for(ii=-mesh;ii<=mesh;ii++)
    for(jj=-mesh;jj<=mesh;jj++)
     temp_smooth[i+width*j] += image[i+ii+width*(j+jj)];
  }


 norm=1.0/((double)((2*mesh+1)*(2*mesh+1)));

 for(i=mesh;i<width-mesh;i++)
  for(j=mesh;j<n-mesh;j++) image[i+width*j]=temp_smooth[i+width*j]*norm;

 

 return;
}
