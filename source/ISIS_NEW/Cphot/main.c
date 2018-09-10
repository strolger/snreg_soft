#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

int main(int argc,char *argv[])
{
  int       i,width,height;          
  DATA_TYPE *im,*im2,image;
  char      filename[256],fileout[256],file_phot[256],config_file[256],
            config_file_phot[256],file2[256];

   
      for (i=1; i<argc; i++)
        {
          if(!(strncmp(argv[i],"-c",2))) sprintf(config_file,"%s",argv[i+1]);
          if(!(strncmp(argv[i],"-d",2))) sprintf(config_file_phot,"%s",argv[i+1]);
          if(!(strncmp(argv[i],"-e",2))) exp_date=atof(argv[i+1]);
          if(!(strncmp(argv[i],"-i",2))) sprintf(filename,"%s",argv[i+1]);
          if(!(strncmp(argv[i],"-j",2))) sprintf(file2,"%s",argv[i+1]);
          if(!(strncmp(argv[i],"-o",2))) sprintf(file_phot,"%s",argv[i+1]);
        }

      read_config(config_file);
      read_phot_conf(config_file_phot);
      im = readfits(&image,filename,&width,&height,1);
      im2 = readfits(&image,file2,&width,&height,1);
      image_width=width;

      make_phot(im,im2,file_phot,filename);

}


