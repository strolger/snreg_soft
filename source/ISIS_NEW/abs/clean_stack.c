#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>


#include "types.h"

void clean_stack(DATA_TYPE **stack, DATA_TYPE *ref,int n, int line_buffer,char flag)
{
  int    i,j,k;
  double mean,sum;
 
 printf("Making median: %i frames buff: %i\n", n, old_buffer);

 
 for(i=0;i<width;i++)
 {
   nlines_read=0;
   for(j=0;j<line_buffer && j<height-old_buffer;j++)
   {
     sum=0.0;
     for(k=0;k<n;k++) {buffer[k] = stack[k][i+width*j]; sum += buffer[k];}
   
     quick_sort (buffer,ndex,n);
     /*for(k=0,mean=0.0;k<n_reject;k++) mean += buffer[ndex[n-1-k]];
       mean /= (double)n_reject;*/
      mean=buffer[ndex[n-1]];
     if( buffer[ndex[n-1-n_reject]]>mean*0.5)
     ref[i+(old_buffer+j)*width] = sum/(double)n;
     else ref[i+(old_buffer+j)*width] = (sum-mean*(double)n_reject)/(double)n;
     ++nlines_read;
   }
 }

 if(flag) old_buffer += line_buffer;

 return;
}

