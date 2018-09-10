#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void get_gaussian(double *axx, double *ayy, double *nteta)
{
 int    i,j,iter;
 double x,y,r,rx,xr,yr,teta,ax,ay,m0,
        mxx,myy,mxy,diff,difft,l1,l2,
        ax_old,ay_old,delta;

 rx=(DATA_TYPE)stack_width/2-1.0;

 mxx=myy=mxy=0.0;
 ax=ax0;
 ay=ay0;
 teta=0.0;
 m0  = 0.0;
 ax_old=ay_old=0.0;
 iter=0;

 while(iter<10)
 { 
 for(i=0;i<stack_width;i++)
  for(j=0;j<stack_width;j++)
  {
   x=(double)(i-stack_width/2);
   y=(double)(j-stack_width/2);
   r=sqrt(x*x+y*y);
  
   if(r<=rx)
   {
    diff  = psf[i+stack_width*j];
    difft = exp(-x*x*ax-y*y*ay);

    mxx  += diff*difft*x*x;
    myy  += diff*difft*y*y;
    mxy  += diff*difft*y*x;
    m0   += diff*difft;
   }
  }  

  ax_old=ax;
  ay_old=ay;
  delta=fabs((mxx+myy)*(mxx+myy)-4.0*(mxx*myy-mxy*mxy));
  l1=0.5*(mxx+myy+sqrt(delta))/m0;
  l2=0.5*(mxx+myy-sqrt(delta))/m0;
  ax=1.0/l1/2.0-ax_old;
  ay=1.0/l2/2.0-ay_old;
  delta=l1-myy;
  if(delta != 0.0) teta=atan(mxy/delta);
  else teta = asin(1.0);
  ++iter;
}

 *axx=ax;
 *ayy=ay;
 *nteta=teta; 
 return;
}
