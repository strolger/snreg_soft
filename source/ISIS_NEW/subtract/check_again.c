#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

char check_again(double *kernel_sol,DATA_TYPE *image)
{
 int       i,istamp,number;
 double    diff,sig,bg,mean,scatter;
 char      check;
 DATA_TYPE *im;
 
 check=0;
 mean=scatter=0.0;
 number=0;
 for(istamp=0;istamp<stamp_number;istamp++)
 {
  if(stamps[istamp].keep)
  {
   im=stamps[istamp].area;
   bg=get_background(stamps[istamp].x, stamps[istamp].y, kernel_sol);
   make_model(istamp, kernel_sol);
   sig=0.0;
   for(i=0;i<stamp_size*stamp_size;i++) 
   {
    diff = temp[i]-im[i]+bg; 
    sig += diff*diff/im[i];
   }
   sig /= (double)(stamp_size*stamp_size);
   stamps[istamp].chi2=sig;
   mean += sig;
   scatter += sig*sig;
   ++number;
   printf("x: %i y: %i sig: %lf keep: %i\n", stamps[istamp].x, stamps[istamp].y,sig,stamps[istamp].keep);
  }
 }

  mean /= (double)number;
 scatter = scatter-(double)number*mean*mean; 
 scatter = sqrt(scatter/(double)number);

 printf("mean: %lf scatter: %lf\n", mean, scatter); 

 for(istamp=0;istamp<stamp_number;istamp++)
 {
  if(stamps[istamp].keep)
  {
    if((stamps[istamp].chi2-mean)>1.0*scatter)
      {
      printf("xr: %i y: %i sig: %lf keep: %i\n", stamps[istamp].x, stamps[istamp].y,stamps[istamp].chi2,stamps[istamp].keep);
       stamps[istamp].keep=0;
       check=1;
      }
  }
 }
 

 return check;
}
