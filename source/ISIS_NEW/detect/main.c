#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

int main(int argc,char *argv[])
{
  int       i;          
  DATA_TYPE *im,*im2,thresh,thresh_n,sigma;
  char      filename[256],filename2[256],fileout[256],config_file[256];
  FILE      *file;


  thresh=0.0;
     
     sprintf(fileout,"bright.data");
     sprintf(config_file, "default_config");
      for (i=1; i<argc; i++)
        switch(tolower(argv[i][1]))
        {
           case 'c': sprintf(config_file, "%s", argv[++i]);
                        break; 
           case 'i': sprintf(filename, "%s", argv[++i]);
                        break; 
           case 'o': sprintf(fileout, "%s", argv[++i]);
                        break; 
           case 'r': sprintf(filename2, "%s", argv[++i]);
                        break; 
           case 't': thresh=atof(argv[++i]);
                        break; 
        }


 read_config(config_file);


 im = readfits(filename);
 im2 = readfits(filename2);
 if(!(temp = (DATA_TYPE *)malloc(width*height*sizeof(DATA_TYPE))))
   {printf("Cannot make malloc\n"); exit(0);}
 printf("Image read\n");
 thresh_n=get_thresh(im,&sigma);
 printf("thresh: %f sigma: %f\n",thresh,sigma);
 /*smooth(im);
   smooth(im2);*/

 printf("thresh: %f sigma: %f\n",thresh,sigma);
 max(im,im2,thresh,2,fileout);

 sprintf(filename,"toto");
 file=fopen(filename,"wb");
 fwrite(im,width*height*4,1,file);
 fclose(file);
}
