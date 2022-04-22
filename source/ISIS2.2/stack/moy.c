#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>


#include "types.h"

void moy(DATA_TYPE **stack, int n, int line_buffer)
{
 int    i,j,k,nn;
 double moy,pix; 
 char   test;



 for(i=0;i<width;i++)
 {
  for(j=0;j<line_buffer;j++)
   {
    moy=0.0;
    nn=0;
    test=1;
    for(k=0;k<n;k++)
    {
     pix=stack[k][i+width*j];
     if(pix>pixmin)
     {
       moy += stack[k][i+width*j];
       ++nn;
     }
     else test=0;
    }
    /*printf("pix: %f %i %i %i\n", stack[k][i+width*j],k,i,j);}*/
    if(test) ref[i+(old_buffer+j)*width] = moy/(double)nn;
    else ref[i+(old_buffer+j)*width] = pixmin;
    /*printf("moy: %lf\n",  moy/(double)n);*/
   }
 }


 old_buffer += line_buffer;

 return;
}
