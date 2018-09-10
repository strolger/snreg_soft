#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

void  check_stamps()
{
 int    ncomps,i,j,im,jm,n,n_old;
 double d,sum,mean,sigma,sigma0;

  ncomps=ncomp_kernel+1;
  mean=sigma=0.0;
 
 
for(i=0;i<stamp_number;i++)
{
 

  for(im=1;im<=ncomps;im++)
  {
    check_vec[im]=stamps[i].scprod[im];
    for(jm=1;jm<=im;jm++)
    {

     check_mat[im][jm]=stamps[i].mat[im][jm];
     check_mat[jm][im]=check_mat[im][jm];

    }
  }
  
  ludcmp(check_mat,ncomps,indx,&d);
  lubksb(check_mat,ncomps,indx,check_vec);

  for(j=0,sum=0.0;j<ncomp_kernel;j++) 
  {
   sum += check_vec[j+1]*sum_vectors[j];
  }
  check_stack[i] = sum;
  stamps[i].norm=sum;
  mean += sum;
  sigma += sum*sum;
   printf("x: %i y: %i %lf\n", stamps[i].x,stamps[i].y,sum);
}
  mean /= (double)stamp_number;


  quick_sort (check_stack, indx, stamp_number);
   mean = check_stack[indx[stamp_number/2]];
    printf("median1: %lf\n",mean);

  sigma=0.0;
  for(i=0;i<stamp_number;i++)
  {
   stamps[i].diff=fabs((stamps[i].norm-mean)*sqrt(stamps[i].sum));
   check_stack[i] = stamps[i].diff;
   sigma +=stamps[i].diff*stamps[i].diff;;
 
  }

  sigma0 = sqrt(sigma/(double)stamp_number);
  n_old    = stamp_number+1;
  n           = stamp_number;
  printf("sigma00: %lf n: %i\n",sigma0,n);

  while(n<n_old)
  {
    n_old  = n;
    sigma = 0.0; 
    n        = 0;
   for(i=0;i<stamp_number;i++)
   {
    if(stamps[i].diff<2.5*sigma0) 
    {
     sigma +=  stamps[i].diff*stamps[i].diff;
     ++n;
    }
   }
   sigma0=sqrt(sigma/(double)n);
  printf("sigma0: %lf n: %i\n",sigma0,n);
  }

 

 for(i=0;i<stamp_number;i++)
 {
   stamps[i].keep=1;
   if(stamps[i].diff>2.5*sigma0)
     {
      stamps[i].keep=0;
      printf("Reject: %lf %i %i\n",stamps[i].diff,stamps[i].x,stamps[i].y);
     }
 }

 

 return;
}
 
