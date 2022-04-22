#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void max(DATA_TYPE *image,DATA_TYPE *image2,DATA_TYPE thresh,int mesh,char *fileout)
{

 int        i,j,m2,k,l,test,buff_obj0,buff_obj,nb_obj,*index,ind,nbw;
 DATA_TYPE  pix0,*pixx,*xxp,*yyp;
 CENT_TYPE  x,y;
 FILE       *file;
 char       var_file[256];
 

 file=fopen(fileout,"w");
 buff_obj0=buff_obj=1000;
 m2=15;

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
            cent(i,j,&x,&y,image2); 
            sprintf(var_file,"lc%i.data",+nb_obj);
            
            fprintf(file,"%lf %lf %i %i %s %lf %lf\n",x,y,i,j,var_file,pix0,image2[i+width*j]);
             ++nb_obj;
          }
      }
    }

 

 

 

 printf("%i objects detected\n", nb_obj);

 fclose(file);
 
 
 return;
}
