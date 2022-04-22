#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"


void make_vectors()
{
 int          ig,idegx,idegy,nvec,i;
 char         test;
 
  



 
 mesh_size=2*half_mesh_size+1;
 kernel_vec = (DATA_TYPE **)malloc(ncomp_kernel*sizeof(DATA_TYPE *));
 filter_x=(DATA_TYPE *)malloc(mesh_size*sizeof(DATA_TYPE));
 filter_y=(DATA_TYPE *)malloc(mesh_size*sizeof(DATA_TYPE));
 sum_vectors=(DATA_TYPE *)malloc(ncomp_kernel*sizeof(DATA_TYPE));
 
 



 nvec=0;
 for(ig=0;ig<ngauss;ig++)
  {
   for(idegx=0;idegx<=deg_fixe[ig];idegx++)
    for(idegy=0;idegy<=deg_fixe[ig]-idegx;idegy++)
    {

      kernel_vec[nvec]=kernel_vector(nvec,idegx,idegy,ig,&test);
      ++nvec;
    }
  }


 free(filter_x);
 free(filter_y);
 
 return;
}
