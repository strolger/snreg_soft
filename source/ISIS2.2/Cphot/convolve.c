#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void convolve(DATA_TYPE *psf,DATA_TYPE *kernel,DATA_TYPE *new_psf)
{
 int      i,j,ic,jc,ik,jk,ip,jp;
 double   q,ps,sum;


     sum=0.0;
     for(i=stack_width/2-conv_size/2;i<=stack_width/2+conv_size/2;i++)
      {
       for(j=stack_width/2-conv_size/2;j<=stack_width/2+conv_size/2;j++)
       {
        q=0.0;
        for(ic=i-half_mesh_size;ic<=i+half_mesh_size;ic++)
         for(jc=j-half_mesh_size;jc<=j+half_mesh_size;jc++)
	 {
          ik=i-ic+half_mesh_size;
          jk=j-jc+half_mesh_size;
          if(ic<0 || ic>= stack_width || jc<0 || jc>= stack_width) ps=0.0;
          else ps=psf[ic+stack_width*jc];

          q  += ps*kernel[ik+jk*mesh_size];
         }
         ip=i-stack_width/2+conv_size/2;
         jp=j-stack_width/2+conv_size/2;
         new_psf[ip+conv_size*jp] = q;
         sum += q;
       }
      }
     
 

 return;
}
