#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"

char *defect_map(DATA_TYPE *imref,DATA_TYPE *image)
{
 char  *map,*mask;
 int   i,j,k,l,mesh1,mesh2,ii,jj;
 char  s[256];
 FILE  *file;

 mesh2=(stamp_size+mesh_size)/2;
 mesh1=stamp_size/2+1;

 map=def_map;
 mask=def_mask;

 for(i=0;i<width*height;i++) map[i]=mask[i]=1;


 j=0;
 for(i=0;i<width;i++)
 { 

  for(j=0;j<height;j++)
  {  
    if(imref[i+width*j]>=SATURATION || image[i+width*j]>=SATURATION)
    {
     for(k=-mesh2;k<=mesh2;k++)
     {
      ii=i+k; 
      for(l=-mesh2;l<=mesh2;l++)
      {
        jj=j+l;

        if(ii>-1 && jj>-1 && ii<width && jj<height) map[ii+width*jj]=0;
      }
     }
    }

    if(imref[i+width*j]<PIX_MIN || image[i+width*j]<PIX_MIN)
    {

     for(k=-mesh1;k<=mesh1;k++)
     {
      ii=i+k; 
      for(l=-mesh1;l<=mesh1;l++)
      {
        jj=j+l;
        if(ii>-1 && jj>-1 && ii<width && jj<height) map[ii+width*jj]=0;
      }
     }
    }
  }
 }


 

 return map;
}
