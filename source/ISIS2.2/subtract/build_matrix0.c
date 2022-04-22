#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"


void build_matrix0()
{
 int       ncomp1,ncomp2,nbg_vec,mat_size,i,j,pixnb,
           istamp,k,ideg1,ideg2,ncomp,i1,i2,j1,j2,
           ivecbg;
 double    **matrix,**matrix0,a1,a2,p0,qw,q;
 DATA_TYPE **vec,xstamp,ystamp;


 
 ncomp1  = ncomp_kernel;
 ncomp2  = ((deg_spatial+1)*(deg_spatial+2))/2;
 ncomp   = ncomp1*ncomp2;
 nbg_vec = ((deg_bg+1)*(deg_bg+2))/2;



 
/****************************************************************************/
/****************************** Initilisations ******************************/
/****************************************************************************/

 
 

  pixnb=stamp_size*stamp_size;


 for(istamp=0;istamp<stamp_number;istamp++)
 {

  vec=stamps[istamp].vectors;
  xstamp=(double)stamps[istamp].x;
  ystamp=(double)stamps[istamp].y;



  for(i=0;i<ncomp1;i++)
   for(j=0;j<=i;j++)
   {
     q=0.0;
     for(k=0;k<pixnb;k++)
     q += vec[i][k]*vec[j][k];
     stamps[istamp].mat[i+1][j+1]=q;
   }

  
 

   for(i1=0;i1<ncomp1;i1++)
    { 
      i=ncomp;
      ivecbg=ncomp1;

      p0=0.0;
      for(k=0;k<pixnb;k++)
       p0 += vec[i1][k]*vec[ivecbg][k];
        stamps[istamp].mat[ncomp1+1][i1+1] = p0;
    }  


     for(k=0,q=0.0;k<pixnb;k++)
      q += vec[ivecbg][k]*vec[ncomp1][k];
      stamps[istamp].mat[ncomp1+1][ncomp1+1]=q;
    

 }



 return;
}
