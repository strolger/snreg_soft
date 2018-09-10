#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

void  kernel_convolve(DATA_TYPE *image, double *kernel_sol)
{
 int       i1,j1,mesh2,i2,j2,conv_step,conv_step2,
           nsteps_x,nsteps_y,i,j,i0,j0,ic,jc,ik,jk;
 char      s[256];
 FILE      *file;
 double    background,q,nq,nq0,qim,qk;

 conv_step=mesh_size;
 
 conv_step2=conv_step/2;
 mesh2=mesh_size/2;
 
 nsteps_x=ceil((double)(width-mesh_size*0)/(double)conv_step);
 nsteps_y=ceil((double)(height-mesh_size*0)/(double)conv_step);

 printf("nxy: %i %i %i %i\n",nsteps_x,nsteps_y,conv_step,height-mesh2);


 for(i1=0;i1<nsteps_x;i1++)
  for(j1=0;j1<nsteps_y;j1++)
  {
    i0=i1*conv_step+mesh2;
    j0=j1*conv_step+mesh2;
    make_kernel(i0+mesh2,j0+mesh2,kernel_sol);

    for(i2=0;i2<conv_step;i2++)
    {
     i=i0+i2;
     if(i>=width-mesh2) break;
     for(j2=0;j2<conv_step;j2++)
     {
       j=j0+j2;
       if(j>=height-mesh2) break;
       q=0.0;
     if(image[i+width*j]>PIX_MIN && image[i+width*j]<SATURATION)
     {
      nq=nq0=0.0;
      for(ic=i-mesh2;ic<=i+mesh2;ic++)
      {
       for(jc=j-mesh2;jc<=j+mesh2;jc++)
       {
        ik=i-ic+mesh2;
        jk=j-jc+mesh2;
        qim = image[ic+width*jc];
        qk = kernel[ik+jk*mesh_size];
        nq0 += qk;
        if(qim>PIX_MIN && qim<SATURATION)
        {
         q  += qim*qk;
         nq += qk;
        }
       }
      }

      conv_image[i+width*j] = q*nq0/nq;
}
else conv_image[i+width*j]=image[i+width*j];


      
     }
    }
  }



  

 return;
}
