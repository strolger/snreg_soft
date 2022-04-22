#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"


double **build_matrix()
{
 int       ncomp1,ncomp2,nbg_vec,mat_size,i,j,pixnb,
           istamp,k,ideg1,ideg2,ncomp,i1,i2,j1,j2,
           ibg,jbg,ivecbg,jj;
 double    **matrix,**matrix0,a1,a2,p0,qw,q;
 DATA_TYPE **vec,xstamp,ystamp;


 
 ncomp1=ncomp_kernel-1;

 ncomp2 = ((deg_spatial+1)*(deg_spatial+2))/2;
 ncomp=ncomp1*ncomp2;

 nbg_vec = ((deg_bg+1)*(deg_bg+2))/2;


/****************************************************************************/
/************ Espace pour la matrice et la matrice temporaire ***************/
/****************************************************************************/

 mat_size = ncomp1*ncomp2+nbg_vec+1;
      printf("Mat_size: %i ncomp2: %i  ncomp1: %i nbg_vec: %i \n",mat_size,ncomp2,ncomp1,nbg_vec);


 matrix = (double **)malloc((mat_size+1)*sizeof(double *));
 for(i=0;i<=mat_size;i++)
 matrix[i] = (double *)malloc((mat_size+1)*sizeof(double));



 wxy = (double **)malloc(stamp_number*sizeof(double *));
 for(i=0;i<stamp_number;i++)
 wxy[i] = (double *)malloc(ncomp2*sizeof(double));

 
/****************************************************************************/
/****************************** Initilisations ******************************/
/****************************************************************************/

 
 for(i=0;i<=mat_size;i++)
   for(j=0;j<=mat_size;j++)
     matrix[i][j] = 0.0;



  pixnb=stamp_size*stamp_size;


 for(istamp=0;istamp<stamp_number;istamp++)
 {


  while(!stamps[istamp].keep)
  {
    ++istamp; 
    if(istamp >= stamp_number) break;
  }

  if(istamp >= stamp_number) break;


  vec=stamps[istamp].vectors;
  xstamp=(double)stamps[istamp].x;
  ystamp=(double)stamps[istamp].y;



   k=0;
   a1=1.0;
   for(ideg1=0;ideg1<=deg_spatial;ideg1++)
   {
      a2=1.0;
      for(ideg2=0;ideg2<=deg_spatial-ideg1;ideg2++)
      {
	 wxy[istamp][k++] = a1*a2;
         a2 *= ystamp;
      }
      a1 *= xstamp;
   }


     matrix0=stamps[istamp].mat;
     for(i=0;i<ncomp;i++)
     {
       i1=i/ncomp2;
       i2=i-i1*ncomp2;

        for(j=0;j<=i;j++)
        {
         j1=j/ncomp2;
          j2=j-j1*ncomp2;
          matrix[i+2][j+2] += wxy[istamp][i2]*wxy[istamp][j2]*matrix0[i1+2][j1+2 ];
        }
      }

   matrix[1][1] += matrix0[1][1];

   for(i=0;i<ncomp;i++)
   {
       i1=i/ncomp2;
       i2=i-i1*ncomp2;
      matrix[i+2][1] += wxy[istamp][i2]*matrix0[i1+2][1];
   }
  


   for(ibg=0;ibg<nbg_vec;ibg++)
   {
    i=ncomp+ibg+1;
    ivecbg=ncomp1+ibg+1;
    for(i1=1;i1<ncomp1+1;i1++)
    { 
      p0=0.0;
      for(k=0;k<pixnb;k++)
       p0 += vec[i1][k]*vec[ivecbg][k];
       
      for(i2=0;i2<ncomp2;i2++)
      {
       jj=(i1-1)*ncomp2+i2+1;
       matrix[i+1][jj+1] +=  p0*wxy[istamp][i2];
      }
    }

      p0=0.0;
      for(k=0;k<pixnb;k++)
       p0 += vec[0][k]*vec[ivecbg][k];
      matrix[i+1][1] +=  p0;

    for(jbg=0;jbg<=ibg;jbg++)
    {
     for(k=0,q=0.0;k<pixnb;k++)
     q += vec[ivecbg][k]*vec[ncomp1+jbg+1][k];
      matrix[i+1][ncomp+jbg+2] += q;
    }
    
    }



 }

  for(i=0;i<mat_size;i++)
    for(j=0;j<=i;j++){
    matrix[j+1][i+1]=matrix[i+1][j+1];
    /*printf("matrix[%i][%i]: %lf\n", i,j,matrix[i+1][j+1]);*/}
  
  



 return matrix;
}
