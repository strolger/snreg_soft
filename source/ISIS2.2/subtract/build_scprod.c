#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

void build_scprod(DATA_TYPE *image)
{

 int       istamp,xc,yc,xi,yi,i1,i2,k,ibg,ncomp1,ncomp2,
           ncomp,nbg_vec,i,ii;
 double    p0,q;
 DATA_TYPE **vec;

 

 ncomp1=ncomp_kernel-1;

 ncomp2 = ((deg_spatial+1)*(deg_spatial+2))/2;
 ncomp=ncomp1*ncomp2;

 nbg_vec = ((deg_bg+1)*(deg_bg+2))/2;

 for(i=0;i<=ncomp+nbg_vec+1;i++) scprod[i]=0.0;


 for(istamp=0;istamp<stamp_number;istamp++)
 {

  while(!stamps[istamp].keep)
  {
    ++istamp; 
   if(istamp >= stamp_number) break;
  }

 if(istamp >= stamp_number) break;

  vec=stamps[istamp].vectors;
  xi=stamps[istamp].x;
  yi=stamps[istamp].y;
   
    p0=stamps[istamp].scprod[1];
    scprod[1] += p0;


   for(i1=1;i1<ncomp1+1;i1++)
   {
      p0=stamps[istamp].scprod[i1+1];
      for(i2=0;i2<ncomp2;i2++)
      {
        ii=(i1-1)*ncomp2+i2+1;
        scprod[ii+1] += p0*wxy[istamp][i2];
      }
   }

   for(ibg=0;ibg<nbg_vec;ibg++)
   {
    q=0.0;
    for(xc=-stamp_size/2;xc<=stamp_size/2;xc++)
     for(yc=-stamp_size/2;yc<=stamp_size/2;yc++)
     {
      k=xc+stamp_size/2+stamp_size*(yc+stamp_size/2);
      q += vec[ncomp1+ibg+1][k]*image[xc+xi+width*(yc+yi)];
     }
     scprod[ncomp+ibg+2] += q;
   }

 }


 /*for(i=0;i<ncomp+nbg_vec+1;i++)
   printf("scprod[%i]: %lf\n",i,scprod[i+1]);*/


 return;
}
