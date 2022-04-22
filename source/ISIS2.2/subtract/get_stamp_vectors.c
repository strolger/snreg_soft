#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

void get_stamp_vectors(DATA_TYPE *imref,DATA_TYPE *image)
{
  int       i,idegx,idegy,ix,ig,nvec; 
  double    x,qe;
  char      s[256],test;
  FILE      *file;

 sum_vectors=(DATA_TYPE *)malloc(ncomp_kernel*sizeof(DATA_TYPE));


/***************************************************************************************/
/********************* Bloc principal: Vecteurs Gaussienne*polynome ********************/
/***************************************************************************************/

  nvec=0;
  sum_vectors[0]=1.0;

  for(i=0;i<stamp_number;i++) 
    cut_stamp(stamps[i].x, stamps[i].y, image,stamps[i].area,&stamps[i].sum);

  for(ig=0;ig<ngauss;ig++)
  {
   
   for(idegx=0;idegx<=deg_fixe[ig];idegx++)
    for(idegy=0;idegy<=deg_fixe[ig]-idegx;idegy++)
    {
      
      kernel_vec[nvec]=kernel_vector(nvec,idegx,idegy,ig,&test);
      
       for(i=0;i<stamp_number;i++)
       {
        xy_conv_stamp(imref,nvec,i,test);
       }
      ++nvec;
    }
  }

/***************************************************************************************/
/******************** Vecteurs correspondants au fond de ciel **************************/
/***************************************************************************************/
  
   for(i=0;i<stamp_number;i++)
  {
    make_bg_vector(i, nvec);
  }

 
  

 return;
}
