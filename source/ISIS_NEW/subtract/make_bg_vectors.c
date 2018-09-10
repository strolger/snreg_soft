#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

void   make_bg_vector(int ns, int nvec)
{
 int    i,j,xi,yi,idegx,idegy,di,dj,nbg_vectors,nv,nk;
 double pix,ax,ay;
 DATA_TYPE *im;

 nbg_vectors=((deg_bg+1)*(deg_bg+2))/2;

 

  xi=stamps[ns].x;
  yi=stamps[ns].y;


  di=xi-stamp_size/2;
  dj=yi-stamp_size/2;
  for(i=xi-stamp_size/2;i<=xi+stamp_size/2;i++)
   for(j=yi-stamp_size/2;j<=yi+stamp_size/2;j++)
   {
    ax=1.0;
    nv=nvec;
    for(idegx=0;idegx<=deg_bg;idegx++)
    {
     ay=1.0;
     for(idegy=0;idegy<=deg_bg-idegx;idegy++)
     {
       im=stamps[ns].vectors[nv];
       im[i-di+stamp_size*(j-dj)] = ax*ay;
      ay *= (double)j;
      ++nv;
     }
     ax *= (double)i;

    }
   }


 return;
}
