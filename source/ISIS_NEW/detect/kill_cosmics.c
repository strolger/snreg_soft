#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"


void kill_cosmics(DATA_TYPE *image,DATA_TYPE thresh)
{
  int        i,j,k,l,n,mesh,test,*index,mesh2,n2,mesh_moy,mesh_smooth,
            kk,ll,mesh_m,ncosmics; 
  DATA_TYPE  pix0,slope[16],*pix,mslope,med,rbg,bg,q,qd,sigma,moy,cosmic_tresh;
 

  printf("adu_el %f\n", adu_el);
  cosmic_tresh=100.0;
  mesh=1;
  rbg=6.0;
  mesh2=ceil(rbg+3.0);
  mesh_moy=1;
  mesh_smooth=1;
  ncosmics=0;
  pix=(DATA_TYPE *)malloc((2*mesh2+1)*(2*mesh2+1)*sizeof(DATA_TYPE));
  index=(int *)malloc((2*mesh2+1)*(2*mesh2+1)*sizeof(int));

  for(i=mesh2;i<width-mesh2;i++)
   for(j=mesh2;j<height-mesh2;j++)
   {
      pix0=image[i+width*j];

      if(pix0>thresh && pix0<saturation)
      {
          test=1;
          for(k=-mesh;k<=mesh;k++)
           for(l=-mesh;l<=mesh;l++)
	   {
             if(image[i+k+width*(j+l)]>pix0) {test=0; break;}
           }
           if(test)
	   {
            n=0;
            for(k=-mesh;k<=mesh;k++)
             for(l=-mesh;l<=mesh;l++)
	     {
	       /*if(image[i+k+width*(j+l)]<pix0)
		 {*/
                   pix[n] = image[i+k+width*(j+l)];
                   slope[n++] = pix0-image[i+k+width*(j+l)];
	       
                 
	     }
	    if(n<1){printf("pouet\n"); exit(0);} 
             quick_sort (slope,index,n);
             mslope=slope[index[n/2]];
             quick_sort (pix,index,n);
             med=pix[index[n/2]];

	   if(mslope > 5.0*sqrt(fabs(med)/adu_el))
	   {
            n2=0;
            for(k=-mesh2;k<=mesh2;k++)
             for(l=-mesh2;l<=mesh2;l++)
	     {
               if(sqrt((double)(k*k)+(double)(l*l)) > rbg)
	       {
                   pix[n2++] = image[i+k+width*(j+l)];
	       }    
	     }

             quick_sort (pix,index,n2);
             bg=pix[index[n2/2]];
             qd=med-bg;
             sigma=sqrt(fabs(med)/adu_el+fabs(bg)/adu_el);
             qd += 5.0*sigma;
            
	     if(mslope>cosmic_tresh*qd) 
	     {
              ++ncosmics;
              for(k=-mesh_smooth;k<=mesh_smooth;k++)
               for(l=-mesh_smooth;l<=mesh_smooth;l++)
	       {
                n=0;
                mesh_m=mesh_moy;
               if(k==0 && l==0) mesh_m=mesh_moy+1; 
               for(kk=-mesh_m;kk<=mesh_m;kk++)
                for(ll=-mesh_m;ll<=mesh_m;ll++)
	         {
                        pix[n++] = image[i+k+kk+width*(j+l+ll)];
		 }
                  
		  quick_sort (pix,index,n);
                  n=index[n/2];
                 
		  image[i+k+width*(j+l)] = pix[n];
	       }

	      }
	     }

	  }  
      }
   }


 free(pix);
 free(index);

 printf("%i cosmics have been removed\n",ncosmics);

 return;
}
