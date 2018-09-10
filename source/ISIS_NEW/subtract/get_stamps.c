#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"

void get_stamps(DATA_TYPE *image, char *map)
{
 int       i,j,ix,jx,nb,nsx,nsy,area_size_x,area_size_y,
           px,py,imax,jmax,mesh2;
 DATA_TYPE pixmax,pix,sum;
 
  nb=0;
 

 
  nsx=nstamps_x;
  nsy=nstamps_y;

  area_size_x=ceil((double)(width-stamp_size-mesh_size)/(double)nsx);
  area_size_y=ceil((double)(height-stamp_size-mesh_size)/(double)nsy);
  


  mesh2=stamp_size/2+mesh_size/2;
  ++nsx; ++nsy;


  for(ix=0;ix<nsx;ix++)
   for(jx=0;jx<nsy;jx++)
   {
     pixmax=max_stamp_thresh;
     px=stamp_size+ix*area_size_x;
     py=stamp_size+jx*area_size_y;
     for(i=px;i<px+area_size_x;i++)
     {
      if(i>=width-mesh2) break;

        for(j=py;j<py+area_size_y;j++)
        {
         if(j>=height-mesh2) break;

         if(max(i,j,2,image,&sum) && map[i+width*j])
         {
              pix = image[i+width*j];
              if(sum>pixmax && pix<SATURATION) {imax=i; jmax=j; pixmax=sum;}
         }
        }
       }
      if(pixmax>max_stamp_thresh)
      {
        stamps[nb].x=imax;
        stamps[nb++].y=jmax;     
      }
     
   }
   stamp_number=nb;
 printf("nb: %i\n", nb);

 return;
}
