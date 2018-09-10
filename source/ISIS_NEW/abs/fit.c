#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

bg_struct fit(DATA_TYPE *image1,DATA_TYPE *image2,int line_buffer)
{

 int       i,j,ii,jj,n,nc;
 double    ax,ay,d;
 bg_struct bg;


 
 i=deg_bg+2;
 nc=((i+1)*(i+2))/2;


  for(ii=0;ii<nc;ii++)
   {
    for(jj=0;jj<nc;jj++)
    {
     mat[ii][jj] = 0.0;
    }

     vec[ii] = 0.0;
   }



 for(i=0;i<width;i++)
 {
  for(j=j0b;j<j0b+line_buffer;j++)
  {  
   lsvec[0]=image1[i+width*(j-j0b)];
   n=1;
   ax=1.0;
   for(ii=0;ii<=deg_bg;ii++)
   {
     ay=1.0;
     for(jj=0;jj<=deg_bg-ii;jj++)
     {
       lsvec[n++] = ax*ay;
       ay *= (double)j;
     }
      ax *= (double)i;
   }

   if(image1[i+width*(j-j0b)]<saturation && image2[i+width*(j-j0b)]<saturation  && image1[i+width*(j-j0b)]>pixmin && image2[i+width*(j-j0b)]>pixmin)
   {

    for(ii=0;ii<n;ii++)
    {
     for(jj=0;jj<=n-ii;jj++)
     {
      mat[ii+1][jj+1] += lsvec[ii]*lsvec[jj];
     }
     vec[ii+1] += lsvec[ii]*image2[i+width*(j-j0b)];
    }
   }
  }
 }


 
   for(ii=0;ii<n;ii++)
     for(jj=0;jj<=n-ii;jj++)
       mat[jj+1][ii+1]=mat[ii+1][jj+1];
      
    
   

  ludcmp(mat,n,indx,&d);
  lubksb(mat,n,indx,vec);

  bg.coeff = (double *)malloc(n*sizeof(double));
  for(i=0;i<n-1;i++) bg.coeff[i]=vec[i+2];
  bg.norm =vec [1];
 


 return bg;
}

