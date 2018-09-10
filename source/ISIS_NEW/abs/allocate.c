#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>


#include "types.h"

void allocate()
{
 int i,nc;


 i=deg_bg+2;

 nc=((i+1)*(i+2))/2;


 printf("nc: %i\n", nc);

 mat = (double **)malloc(nc*sizeof(double *));
 for(i=0;i<nc;i++)
  mat[i] = (double *)malloc(nc*sizeof(double));

 vec = (double *)malloc(nc*sizeof(double));
 lsvec = (double *)malloc(nc*sizeof(double));
 indx  = (int *)malloc(nc*sizeof(int));

}
