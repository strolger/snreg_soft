#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>


#include "types.h"

void moy(DATA_TYPE **stack, DATA_TYPE *ref, int n, int line_buffer, char flag)
{
 int    i,j,k,nn;
 double moy,pix; 
 char   test;



 for(i=0;i<width;i++)
 {
  for(j=0;j<line_buffer && j<height-old_buffer;j++)
   {
    moy=0.0;
    nn=0;
    for(k=0;k<n;k++)
    {
       pix=stack[k][i+width*j];
       moy += stack[k][i+width*j];
       ++nn;
    }
    /*printf("pix: %f %i %i %i\n", stack[k][i+width*j],k,i,j);}*/
    ref[i+(old_buffer+j)*width] = moy/(double)nn;
    /*printf("moy: %lf\n",  moy/(double)n);*/
   }
 }


 if(flag) old_buffer += line_buffer;

 return;
}
