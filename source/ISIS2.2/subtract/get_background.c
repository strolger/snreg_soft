#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

double get_background(int xi, int yi, double *scprod)
{
 double  background,ax,ay,xd,yd;
 int     i,j,k;


 background=0.0;
 k=1;
 xd=(double)xi;
 yd=(double)yi;


 ax=1.0;
 for(i=0;i<=deg_bg;i++)
 {
  ay=1.0;
  for(j=0;j<=deg_bg-i;j++)
  {
    background += scprod[ncomp_background+k++]*ax*ay;
    ay *= yd;
  }
  ax *= xd;
 }



 return background;

}
