#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void max(DATA_TYPE *image,DATA_TYPE bg,DATA_TYPE sigma,int mesh,char *fileout)
{

 int        i,j,m2,k,l,test,buff_obj0,buff_obj,nb_obj,*index,
            ind,nbw,mesh_phot,ix,iy,sub_x,sub_y,nsub_x,nsub_y,
            *area_count,narea_max;
 DATA_TYPE  pix0,*pixx,*pixm,*xxp,*yyp,thresh,flux,med,sig,subwidth_x,
            subwidth_y;
 CENT_TYPE  x,y;
 FILE       *file;
 

 nsub_x=3;
 nsub_y=3;
 subwidth_x=(double)width/(double)nsub_x+1.0;
 subwidth_y=(double)height/(double)nsub_y+1.0;

 area_count = (int *)malloc(nsub_x*nsub_y*sizeof(int));
 for(i=0;i<nsub_x*nsub_y;i++) area_count[i]=0;

 thresh = bg + 5.0*sigma;
 file=fopen(fileout,"w");
 buff_obj0=buff_obj=1000;
 m2=15;
 mesh_phot=5;
 xxp=(DATA_TYPE *)malloc(buff_obj0*sizeof(DATA_TYPE));
 yyp=(DATA_TYPE *)malloc(buff_obj0*sizeof(DATA_TYPE));
 pixx=(DATA_TYPE *)malloc(buff_obj0*sizeof(DATA_TYPE));
 pixm=(DATA_TYPE *)malloc(buff_obj0*sizeof(DATA_TYPE));
 nb_obj=0;


  for(i=m2;i<width-m2;i++)
   for(j=m2;j<height-m2;j++)
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
           flux=0.0;
           for(ix=-mesh_phot;ix<=mesh_phot;ix++) 
              for(iy=-mesh_phot;iy<=mesh_phot;iy++)
               flux += image[i-ix+width*(j-iy)];
            flux /= (DATA_TYPE)((2*mesh_phot+1)*(2*mesh_phot+1));
            flux = pix0-flux;
           if(flux > 0.0)
           {
            cent(i,j,&x,&y,image); 
            if(nb_obj>buff_obj-1)
            {
              buff_obj += buff_obj0;
              if(!(xxp=(DATA_TYPE *)realloc(xxp,buff_obj*sizeof(DATA_TYPE))))
              printf("Realloc error at local_max\n");
              if(!(yyp=(DATA_TYPE *)realloc(yyp,buff_obj*sizeof(DATA_TYPE))))
              printf("Realloc error at local_max\n");
              if(!(pixx=(DATA_TYPE *)realloc(pixx,buff_obj*sizeof(DATA_TYPE))))
              printf("Realloc error at local_max\n");
	      if(!(pixm=(DATA_TYPE *)realloc(pixm,buff_obj*sizeof(DATA_TYPE))))
              printf("Realloc error at local_max\n");
            }
              xxp[nb_obj]=x;
              yyp[nb_obj]=y;
              pixx[nb_obj] = flux;
            
              pixm[nb_obj] = (pix0-bg)/(fabs(flux)+1.0);
              pixm[nb_obj] = image[i+width*j];
              ++nb_obj;
          }
         }
      }
    }


printf("Alloc: %i\n",nb_obj);
 index=(int *)malloc(nb_obj*sizeof(int));
printf("Sorting\n");


 quick_sort (pixx,index,nb_obj);

printf("Sorting Done\n");




  narea_max=(int)floor((double)nb_obj/(double)(nsub_x*nsub_y));
  narea_max=100.0;

  printf("Writing\n");

  for(i=0,nbw=0;i<nb_obj;i++)
  {           
   ind = index[nb_obj-i-1];
   sub_x = (int)floor(xxp[ind]/subwidth_x);
   sub_y = (int)floor(yyp[ind]/subwidth_y);
   ++area_count[sub_x+nsub_x*sub_y];

   if(area_count[sub_x+nsub_x*sub_y]<narea_max)
   {
    fprintf(file,"%lf %lf %lf %lf %lf\n",xxp[ind],yyp[ind],pixx[ind],pixm[ind]);
    ++nbw;
   }
  }

  /*for(i=0;i<nsub_x*nsub_y;i++) printf("N: %i\n", area_count[i]);*/

 printf("%i objects detected %i objects written\n", nb_obj, nbw);

 fclose(file);
 free(index);
 free(xxp);
 free(yyp);
 free(pixx);
 free(area_count);

 return;
}
