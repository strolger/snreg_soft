#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

void build_scprod0(DATA_TYPE *image)
{

 int       istamp,xc,yc,xi,yi,i1,i2,k,ncomp1,ncomp2,ncomp,nbg_vec,i;
 double    p0,q;
 DATA_TYPE **vec;

 

 ncomp1=ncomp_kernel;

 ncomp2 = ((deg_spatial+1)*(deg_spatial+2))/2;
 ncomp=ncomp1*ncomp2;

 nbg_vec = ((deg_bg+1)*(deg_bg+2))/2;



 for(istamp=0;istamp<stamp_number;istamp++)
 {
  
  vec=stamps[istamp].vectors;
  xi=stamps[istamp].x;
  yi=stamps[istamp].y;
   

   for(i1=0;i1<ncomp1;i1++)
   {
    p0=0.0;
    for(xc=-stamp_size/2;xc<=stamp_size/2;xc++)
     for(yc=-stamp_size/2;yc<=stamp_size/2;yc++)
     {
      k=xc+stamp_size/2+stamp_size*(yc+stamp_size/2);
      p0 += vec[i1][k]*image[xc+xi+width*(yc+yi)];
     }
      stamps[istamp].scprod[i1+1]=p0;
   }
      
  
    q=0.0;
    for(xc=-stamp_size/2;xc<=stamp_size/2;xc++)
     for(yc=-stamp_size/2;yc<=stamp_size/2;yc++)
     {
      k=xc+stamp_size/2+stamp_size*(yc+stamp_size/2);
      q += vec[ncomp1][k]*image[xc+xi+width*(yc+yi)];
     }
     stamps[istamp].scprod[ncomp1+1]=q;
   

  
 }


 
 return;
}
