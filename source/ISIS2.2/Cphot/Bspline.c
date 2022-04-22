#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void Bspline(DATA_TYPE *psf,DATA_TYPE *Bpsf,double dx, double dy)
{

  int xi1,yi1;


   splie2(x1a,x2a,psf,stack_width,stack_width,y2a);



   for(xi1=0;xi1<stack_width;xi1++)
   {
     for(yi1=0;yi1<stack_width;yi1++)
     { 
      xxp[yi1]=(DATA_TYPE)xi1+dx;
      yyp[yi1]=(DATA_TYPE)yi1+dy;
     }

     splin2(x1a,x2a,psf,y2a,stack_width,stack_width,xxp,yyp,yy,uv,ytmp,yytmp);
 
     for(yi1=0;yi1<stack_width;yi1++) 
     {
      Bpsf[xi1+yi1*stack_width]=yy[yi1];
     }
   }


 return;
}
