#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void make_phot(object *objlist,DATA_TYPE *im,double ax,double ay, double teta,int nn,double rad,char *name)
{

 int       xi,yi,i,xi2,yi2,mesh;
 DATA_TYPE x,y;
 double    x1,y1,xr,yr,norm,s1,s2,pho,bg,
           phot,g,r,s3,qx,correl;
 FILE      *outfile;




 outfile=fopen(name,"w"); 
 norm=asin(1.0)*2.0/sqrt(ax*ay);
 mesh=ceil(rad)+1;

 for(i=0;i<nn;i++)
 {
   x=objlist[i].x;
   y=objlist[i].y;
   xi=objlist[i].xi;
   yi=objlist[i].yi;
   bg=objlist[i].bg;
    s1=s2=s3=0.0;

   for(xi2=xi-mesh;xi2<=xi+mesh;xi2++)
    for(yi2=yi-mesh;yi2<=yi+mesh;yi2++)
    {
     x1=x-(double)xi2;
     y1=y-(double)yi2;
     r = sqrt(x1*x1+y1*y1);
     xr = cos(teta)*x1+sin(teta)*y1;
     yr = cos(teta)*y1-sin(teta)*x1;
     g  = exp(-xr*xr*ax-yr*yr*ay)/norm;
     if(r<=rad)
     {
      qx=(im[xi2+yi2*width]-bg);
      s1 += qx*g;
      s2 += g*g;
      s3 += qx*qx;
     }
    }

    phot=s1/s2;
    correl=s1/sqrt(s2*s3);
 fprintf(outfile,"%9.3f %9.3f %5i %5i %12.3f %12.3f %9.3f\n",x,y,xi,yi,objlist[i].aperphot,phot,bg);
 
 }
 fclose(outfile);

 return;
}
