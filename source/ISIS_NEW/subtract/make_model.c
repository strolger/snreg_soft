#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

void make_model(int stamp_num,double *kernel_sol)
{
 
 int       i1,k,ix,iy,i,kernel_out,xi,yi;
 double    ax,ay,coeff;
 DATA_TYPE *vector;
 char      s[256];
 FILE      *file;


 
 xi=stamps[stamp_num].x;
 yi=stamps[stamp_num].y;
 for(i=0;i<stamp_size*stamp_size;i++) temp[i]=0.0;

 
  vector=stamps[stamp_num].vectors[0];
   coeff=kernel_sol[1];

 for(i=0;i<stamp_size*stamp_size;i++)  temp[i] += coeff*vector[i];

  
 k=2;
 for(i1=1;i1<ncomp_kernel;i1++)
 {
   vector=stamps[stamp_num].vectors[i1];
   coeff=0.0; 
   ax=1.0;
   for(ix=0;ix<=deg_spatial;ix++)
   {
    ay=1.0;
    for(iy=0;iy<=deg_spatial-ix;iy++)
    {
      coeff += kernel_sol[k++]*ax*ay;
      ay *= (double)yi;
    }
     ax *= (double)xi;
   }

    
  for(i=0;i<stamp_size*stamp_size;i++)
  {
    temp[i] += coeff*vector[i];
  }
  
}

 
 return;
}
 
