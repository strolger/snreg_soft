#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

void xy_conv_stamp(DATA_TYPE *image,int n,int istamp,char test)
{
  int       i,j,xc,xij,sub_width,xi,yi;
  char      s[256];
  FILE      *file;
  DATA_TYPE *v0,*imc;


 
  xi  = stamps[istamp].x;
  yi  = stamps[istamp].y;
  imc = stamps[istamp].vectors[n];


  sub_width=stamp_size+mesh_size-1;

  

  for(i=xi-stamp_size/2-mesh_size/2;i<=xi+stamp_size/2+mesh_size/2;i++)
  {
   
    for(j=yi-stamp_size/2;j<=yi+stamp_size/2;j++)
    {
      xij=i-xi+sub_width/2+sub_width*(j-yi+stamp_size/2);
      temp[xij]=0.0;
      for(xc=-mesh_size/2;xc<=mesh_size/2;xc++)
      temp[xij] += image[i+width*(j+xc)]*filter_y[mesh_size/2-xc];
    }
  }



  for(j=-stamp_size/2;j<=stamp_size/2;j++)
  {
    for(i=-stamp_size/2;i<=stamp_size/2;i++)
    {  
      xij=i+stamp_size/2+stamp_size*(j+stamp_size/2);
      imc[xij]=0.0;
      for(xc=-mesh_size/2;xc<=mesh_size/2;xc++)
      imc[xij] += temp[i+xc+sub_width/2+sub_width*(j+stamp_size/2)]*filter_x[mesh_size/2-xc];
    }
  }

 if(test) 
 {
  v0 = stamps[istamp].vectors[0];
  for(i=0;i<stamp_size*stamp_size;i++) imc[i] -= v0[i];
 }
 

 return;
}
