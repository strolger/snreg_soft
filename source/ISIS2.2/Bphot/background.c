#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"



void init_bg(double rad1,double rad2)
{
 int i,j,k,rad,r1,r2;
  

 bg_edge=(int)ceil(rad2);
 
  r1=rad1*rad1;
  r2=rad2*rad2;
  k=0;
  if(!(bg_list=(int *)malloc((2*bg_edge+1)*(2*bg_edge+1)*sizeof(int))))
          printf("Malloc error at init_bg\n");
  for(i=-bg_edge;i<=bg_edge;i++)
   for(j=-bg_edge;j<=bg_edge;j++)
   {
     rad=i*i+j*j;
     if(rad>r1 && rad<r2) bg_list[k++]=i+width*j;
   }
  if(!(bg_temp=(DATA_TYPE *)malloc(k*sizeof(DATA_TYPE))))
          printf("Malloc error at init_bg\n");
  if(!(bg_index=(int *)malloc(k*sizeof(int))))
          printf("Malloc error at init_bg\n");
  nbg=k;  
 printf("nbg: %i %i\n", nbg, bg_edge);

 return;
}


double background(int x,int y,DATA_TYPE *image,double *sigma)
{
 double bg,moy,sig,moy2,sigma2;
 int    i,j,n,offset,ix,iy,off,count,count0;


 
 offset=x+y*width;
if(x>=bg_edge && x<width-bg_edge && y>=bg_edge && y<height-bg_edge) 
{

 for(i=0;i<nbg;i++) 
    bg_temp[i]=image[bg_list[i]+offset];
   
 n=nbg;

 quick_sort(bg_temp, bg_index, n);
 bg=bg_temp[bg_index[n/2]];


  
  *sigma=0.0;
  for(i=0;i<nbg;i++) *sigma += bg_temp[i]/ADU_EL;
   *sigma=sqrt(*sigma/(double)nbg);
 
}
else
{
  n=0;
  for(i=0;i<nbg;i++)
  {

   off=bg_list[i]+offset;
   iy=off/width;
   ix=off-iy*width;

   if(ix>=bg_edge && ix<width-bg_edge && iy>=bg_edge && iy<height-bg_edge)
   bg_temp[n++]=image[off];
  }
if(n>1)
{
 quick_sort(bg_temp, bg_index, n);
 bg=bg_temp[bg_index[n/2]];
}
else bg=0.0;
*sigma=0.0;

}


 return bg;
}
