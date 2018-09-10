#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>


#include "types.h"

double get_background(int i,int j,bg_struct bgnorm)
{
 int    ii,jj,n;
 double ax,ay,bg;


 
  n=0;
  ax=1.0;
  bg=0.0;



   for(ii=0;ii<=deg_bg;ii++)
   {
     ay=1.0;
     for(jj=0;jj<=deg_bg-ii;jj++)
     {
       bg += bgnorm.coeff[n++]*ax*ay;
       ay *= (double)j;
     }
      ax *= (double)i;
   }



 return bg;
}
