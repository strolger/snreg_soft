#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>


#include "types.h"

void smooth(DATA_TYPE *image)
{
  int    i,j,k,l,m2;
  double moy; 


  
  m2=1;
 
  for(i=m2;i<width-m2;i++)
   for(j=m2;j<height-m2;j++)
    {
      moy=0.0;
      for(k=-m2;k<=m2;k++)
       for(l=-m2;l<=m2;l++)
        moy += image[i+k+width*(j+l)];
      
        temp[i+width*j] = moy/(DATA_TYPE)((2*m2+1)*(2*m2+1));
 
    }


   for(i=m2;i<width-m2;i++)
    for(j=m2;j<height-m2;j++)
      image[i+width*j]=temp[i+width*j];
  

 return;
}

