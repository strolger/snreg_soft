#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void aperphot(object *objlist,DATA_TYPE *image,DATA_TYPE rad,int *nn)
{
 int       i,xi,yi,ii,jj,mesh,*index,nb,nk;
 DATA_TYPE x,y,r,dx,dy,phot,phot2,bg,pixmax,*temp,ratio,sig,sigma;
 object    *nobjlist,obj;
 

 mesh=stack_width/2;
 nobjlist=(object *)malloc(NSTARS*sizeof(object));
 nb=0;

 for(i=0;i<NSTARS;i++)
 {
  
  xi     = objlist[i].xi;
  yi     = objlist[i].yi;
  x      = objlist[i].x;
  y      = objlist[i].y;
  bg     = objlist[i].bg;
  pixmax = objlist[i].pixmax;
  sigma  = objlist[i].sigma;

  phot=sig=phot2=0.0;
  for(ii=xi-mesh;ii<=xi+mesh;ii++)
   for(jj=yi-mesh;jj<=yi+mesh;jj++)
    {
      phot2 += image[ii+width*jj]-bg;
      dx=x-(DATA_TYPE)ii;
      dy=y-(DATA_TYPE)jj;
      r=sqrt(dx*dx+dy*dy);
      if(r<=rad) 
      {
       phot += image[ii+width*jj]-bg;
       sig += sigma*sigma;
      }
    }

 ratio = pixmax/phot2;
 sig   = sqrt(sig);

 if(ratio>0.0 && ratio<ratio_max && phot>100.0*sig)
 {
  nobjlist[nb]=objlist[i];
  nobjlist[nb].aperphot=phot;
  nobjlist[nb].aperphot2=phot2;
  nobjlist[nb++].ratio=ratio;
 }
  
 }
  temp = (DATA_TYPE *)malloc(nb*sizeof(DATA_TYPE));
  index = (int *)malloc(nb*sizeof(int));

  for(i=0;i<nb;i++) temp[i]=nobjlist[i].ratio;
  quick_sort(temp,index,nb);

 nk=0;
 for(i=Nfirst;i<Nfirst+Nkeep;i++)
 {
   if(i>=nb) break;
   obj=nobjlist[index[nb-i-1]];
   objlist[nk++] = obj;
 } 

  free(temp);
  free(index);

 *nn=nk; 

 return;
}
