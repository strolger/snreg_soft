#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

DATA_TYPE *kernel_vector(int n, int deg_x, int deg_y,int ig,char *test)
{
 DATA_TYPE *vector,*kernel0;
 int       i,j,dx,dy,ix;
 double    sum,sum_x,sum_y,x,qe;

 vector=(DATA_TYPE *)malloc(mesh_size*mesh_size*sizeof(DATA_TYPE));
 dx=(deg_x/2)*2-deg_x;
 dy=(deg_y/2)*2-deg_y;
 sum_x=sum_y=sum=0.0;
 *test=0;



 for(ix=0;ix<mesh_size;ix++)
 {
   x            = (double)(ix-mesh_size/2);
   qe           = exp(-x*x*sigma_gauss[ig]);
   filter_x[ix] = qe*power(x,deg_x);
   filter_y[ix] = qe*power(x,deg_y);
   sum_x       += filter_x[ix];
   sum_y       += filter_y[ix];
 }

 if(n>0) kernel0=kernel_vec[0];
    sum_vectors[0]=1.0;


 if(dx == 0 && dy == 0)
 {
  
  for(ix=0;ix<mesh_size;ix++) filter_x[ix] /= sum_x;
  for(ix=0;ix<mesh_size;ix++) filter_y[ix] /= sum_y;


  for(i=0;i<mesh_size;i++)
   for(j=0;j<mesh_size;j++)
   {
    vector[i+mesh_size*j]=filter_x[i]*filter_y[j];
   }
 
  if(n>0) 
  {
   for(i=0;i<mesh_size*mesh_size;i++)
   {
    vector[i] -= kernel0[i];
    sum += vector[i];
   }
   sum_vectors[n]=0.0;
   *test=1;
  }
  


 }
 else
 {
  for(i=0;i<mesh_size;i++)
   for(j=0;j<mesh_size;j++)
   {
    vector[i+mesh_size*j]=filter_x[i]*filter_y[j];
    sum += vector[i+mesh_size*j];
   }

   sum_vectors[n]=0.0;
 }

 return vector;
}
