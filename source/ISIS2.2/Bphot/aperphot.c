#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void aperphot(object *objlist,DATA_TYPE *image)
{
 int       i,xi,yi,ii,jj,mesh,*index,nb,nk;
 DATA_TYPE x,y,r,dx,dy,phot,phot2,bg,pixmax,*temp,ratio,sig,sigma;
 object    obj;
 

 mesh=ceil(radphot)+1;
 nb=0;


  
  xi     = objlist->xi;
  yi     = objlist->yi;
  x      = objlist->x;
  y      = objlist->y;
  bg     = objlist->bg;
  pixmax = objlist->pixmax;
  sigma  = objlist->sigma;

  phot=sig=phot2=0.0;
  for(ii=xi-mesh;ii<=xi+mesh;ii++)
   for(jj=yi-mesh;jj<=yi+mesh;jj++)
    {
      phot2 += image[ii+width*jj]-bg;
      dx=x-(DATA_TYPE)ii;
      dy=y-(DATA_TYPE)jj;
      r=sqrt(dx*dx+dy*dy);
      if(r<=radphot) 
      {
       phot += image[ii+width*jj]-bg;
       sig += sigma*sigma;
      }
    }



 
  objlist->aperphot=phot;
 
 
  
 


 return;
}
