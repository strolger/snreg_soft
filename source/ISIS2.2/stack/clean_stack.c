#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>


#include "types.h"

void clean_stack(DATA_TYPE **stack, int n, int line_buffer)
{
  int       i,j,k,nc,nc2;
  DATA_TYPE median,moy,nmoy,sigs;
 
 printf("Making median 0\n");


 for(i=0;i<width;i++)
 {
   for(j=0;j<line_buffer;j++)
   {
     for(k=0;k<n;k++) buffer[k] = stack[k][i+width*j];
     quick_sort (buffer,ndex,n);
     median = buffer[ndex[n/2]];
     moy=median;
     nc2=n;
     nc=n+1;
     
     while(nc>nc2 && nc2>2)
     {
      nc=nc2;
      nmoy=sigs=0.0;
      for(k=0;k<n;k++) sigs += fabs(buffer[k]-moy);
      sigs /= (double)n;

      for(k=0,nc2=0;k<n;k++)
      {
       if(fabs(buffer[k]-moy) < 3.0*sigs)
       {
         ++nc2;
         nmoy += buffer[k];
       }
      }
      if(nc2>0) moy = nmoy/(double)nc2;
      else {moy=0.0; break;}
     }
     /*printf("median: %lf moy: %lf n: %i %i\n",median,moy,nc,n);*/
     ref[i+(old_buffer+j)*width] = moy;
    
   }
 }

 old_buffer += line_buffer;

 return;
}

