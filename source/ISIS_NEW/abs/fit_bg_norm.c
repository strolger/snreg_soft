#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>


#include "types.h"

void fit_bg_norm(DATA_TYPE **stack,int n,int line_buffer, bg_struct *bgnorm)
{

 int i;


 for(i=1;i<n;i++)
 {
  bgnorm[i] = fit(stack[0],stack[i],line_buffer);
  printf("norm: %lf bg: %lf\n",bgnorm[i].norm,bgnorm[i].coeff[0]);
 }

 
  j0b += line_buffer;

 return;
}
