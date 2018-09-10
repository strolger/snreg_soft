#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void max(DATA_TYPE *image,DATA_TYPE thresh,int mesh,char *fileout)
{

 int       i,j,m2,k,l,test;
 DATA_TYPE pix0;
 CENT_TYPE x,y;
 FILE      *file;
 

 file=fopen(fileout,"w");

 m2=15;

  for(i=m2;i<width-m2;i++)
   for(j=m2;j<height-m2;j++)
    {
      if(image[i+width*j]>thresh)
      {
          pix0=image[i+width*j];
          test=1;
          for(k=-mesh;k<=mesh;k++)
           for(l=-mesh;l<=mesh;l++)
	   {
	     if(image[i+k+width*(j+l)]>pix0) {test=0; break;}
           }
	  if(test) 
	  {
            cent(i,j,&x,&y,image); 
            fprintf(file,"%lf %lf %lf\n",x,y,pix0);
          }
      }
    }

   fclose(file);

 return;
}
