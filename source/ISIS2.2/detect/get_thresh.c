#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

DATA_TYPE get_thresh(DATA_TYPE *image, DATA_TYPE *sigma)
{
 int       bxi,byi,bxs,bys,ntemp,*index,i,j,k;
 DATA_TYPE *temp_im,threshold;
 double    sg,sgg;


 bxi = width/2-32;
 byi = height/2-32;
 bxs = width/2+32;
 bys = height/2+32;
 
 if(bxi<0) bxi=0;
 if(byi<0) byi=0;
 if(bxs<0) bxs=width;
 if(bys<0) bys=height;

 ntemp=(bxs-bxi+1)*(bys-byi+1);
 if(!(temp_im = (DATA_TYPE *)malloc(ntemp*sizeof(DATA_TYPE))))
   {printf("Cannot allocate memory\n"); exit(0);}
 if(!(index = (int *)malloc(ntemp*sizeof(int))))
   {printf("Cannot allocate memory\n"); exit(0);}
 
 k=0;
 for(i=bxi;i<bxs;i++)
  for(j=byi;j<bys;j++)
   {
    temp_im[k++] = image[i+width*j];
   }


  quick_sort (temp_im,index,k);
  threshold=temp_im[index[k/2]];



  for(i=0;i<k;i++)
  {
   
    temp_im[i] = fabs(temp_im[i]-threshold);
  }

  quick_sort (temp_im,index,k);
  *sigma=temp_im[index[k/2]];


 free(temp_im);
 free(index);

 return threshold;
}
