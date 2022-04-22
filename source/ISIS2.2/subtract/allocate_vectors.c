#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

void allocate_vectors()
{
 int i,j,nc,nbg_vectors;

 nc=ncomp_kernel+2;
 nbg_vectors=((deg_bg+1)*(deg_bg+2))/2;

 for(i=0;i<stamp_number;i++)
 {

   if(!(stamps[i].vectors=(DATA_TYPE **)malloc((ncomp_kernel+nbg_vectors)*sizeof(DATA_TYPE *))))
    {
     printf("Cannot Allocate Stamp Vector\n"); 
     exit(0);
    }

   for(j=0;j<ncomp_kernel+nbg_vectors;j++)
   {
    if(!(stamps[i].vectors[j]=(DATA_TYPE *)malloc(stamp_size*stamp_size*sizeof(DATA_TYPE))))
    {
     printf("Cannot Allocate Stamp Vector\n"); 
     exit(0);
    }
   }
    if(!(stamps[i].area=(DATA_TYPE *)malloc(stamp_size*stamp_size*sizeof(DATA_TYPE))))
    {
     printf("Cannot Allocate Stamp Vector\n"); 
     exit(0);
    } 


   stamps[i].mat = (double **)malloc(nc*sizeof(double *));
   for(j=0;j<nc;j++) stamps[i].mat[j] = (double *)malloc(nc*sizeof(double));
   stamps[i].scprod =  (double *)malloc(nc*sizeof(double));
 }
  

 def_map=(char *)malloc(width*height*sizeof(char));
 def_mask=(char *)malloc(width*height*sizeof(char));

  check_mat=(double **)malloc(nc*sizeof(double *));
  for(j=0;j<nc;j++) 
  check_mat[j] = (double *)malloc(nc*sizeof(double));
  check_vec=(double *)malloc(nc*sizeof(double));

 return;
}
