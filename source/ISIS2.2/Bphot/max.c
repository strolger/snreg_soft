#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

object *max(DATA_TYPE *image,DATA_TYPE thresh,int mesh,DATA_TYPE rad11,DATA_TYPE rad22)
{

 int       i,j,m2,k,l,test,nb_obj,buff_obj,buff_obj0;
 int       *index;
 DATA_TYPE *temp,pix0;
 CENT_TYPE x,y;
 FILE      *file;
 object    *obj,*objlist;
 double    bg,sigma,rad1,rad2;
 
 rad1=rad11;
 rad2=rad22;
 ADU_EL=1.0;





 m2=stack_width+2;
 buff_obj0=1000;
 buff_obj=buff_obj0;
 if(!(obj=(object *)malloc(buff_obj0*sizeof(object)))) 
   printf("Malloc error at local_max\n");
 
  nb_obj=0;



  for(i=m2;i<width-m2;i++)
   for(j=m2;j<height-m2;j++)
    {

     pix0=image[i+width*j];

      if(pix0>thresh && pix0<saturation && pix0>pix_min)
      {
          test=1;
          for(k=-mesh;k<=mesh;k++)
           for(l=-mesh;l<=mesh;l++)
	   {
	     if(image[i+k+width*(j+l)]>pix0) {test=0; break;}
           }


	  if(test) 
	  {
           for(k=-mesh2;k<=mesh2;k++)
            for(l=-mesh2;l<=mesh2;l++)
	    {
	     if(image[i+k+width*(j+l)]>pix0 || image[i+k+width*(j+l)]<pix_min) {test=0;}
            }
	  }
 
           if(test) 
	   { 
             if(nb_obj>buff_obj-1)
             {
              buff_obj += buff_obj0;
              if(!(obj=(object *)realloc(obj,buff_obj*sizeof(object))))
              printf("Realloc error at local_max\n");
             }
	     bg=background(i,j,image,&sigma);
            if(sigma>0.0001)
	    {
             cent(i,j,&x,&y,image);

             obj[nb_obj].x=x;
             obj[nb_obj].y=y;
             obj[nb_obj].xi=i;
             obj[nb_obj].yi=j;
             obj[nb_obj].bg=bg;

             obj[nb_obj].sigma=sigma;
             obj[nb_obj].pixmax=pix0-bg;
             ++nb_obj;
	    }

          }
      }
    }

printf("nb_obj: %i %i\n",nb_obj,NSTARS);

  temp = (DATA_TYPE *)malloc(nb_obj*sizeof(DATA_TYPE));

  index = (int *)malloc(nb_obj*sizeof(int));

  for(i=0;i<nb_obj;i++) temp[i]=obj[i].pixmax;
  quick_sort(temp,index,nb_obj);

   nnstars=nb_obj;
   if(nnstars>NSTARS_MAX) nnstars=NSTARS_MAX;
  objlist = (object *)malloc(nnstars*sizeof(object));




  for(i=0;i<nnstars;i++)
  {
    objlist[i]=obj[index[nb_obj-i-1]];
    aperphot(&objlist[i],image);
  }


  free(obj);    
  free(temp);
  free(index);

 return objlist;
}
