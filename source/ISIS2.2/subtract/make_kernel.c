#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

void make_kernel(int xi, int yi,double *kernel_sol)
{
 
 int    i1,k,ix,iy,i,kernel_out;
 double ax,ay,yd,xd;
 char      s[256];
 FILE      *file;


 xd=(double)xi;
 yd=(double)yi;

 k=2;

 for(i1=1;i1<ncomp_kernel;i1++)
 {
   kernel_coeffs[i1]=0.0;
   ax=1.0;
   for(ix=0;ix<=deg_spatial;ix++)
   {
    ay=1.0;
    for(iy=0;iy<=deg_spatial-ix;iy++)
    {
      kernel_coeffs[i1] += kernel_sol[k++]*ax*ay;
      ay *= (double)yd;
    }
     ax *= (double)xd;
   }
 }

 kernel_coeffs[0]=kernel_sol[1]; 
 
 
 for(i=0;i<mesh_size*mesh_size;i++) kernel[i] = 0.0;
  sum_kernel=0.0;
  for(i=0;i<mesh_size*mesh_size;i++)
  {
     for(i1=0;i1<ncomp_kernel;i1++)
     {
          kernel[i] += kernel_coeffs[i1]*kernel_vec[i1][i];
     }
    sum_kernel += kernel[i];    
  }


 
 return;
}
 
