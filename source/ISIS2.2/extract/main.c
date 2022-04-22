#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

int main(int argc,char *argv[])
{
  int       i;          
  DATA_TYPE *im,thresh,sigma,cthresh;
  char      filename[256],fileout[256],config_file[256];




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
           case 's': cthresh=atof(argv[++i]);
	     break;
        }

 read_config(config_file);




 im = readfits(filename);
 printf("Image read\n");
 thresh=get_thresh(im,&sigma);
 printf("thresh: %f sigma: %f\n",thresh,sigma);

 kill_cosmics(im,thresh,cthresh);
 printf("thresh: %f sigma: %f\n",thresh,sigma);
 max(im,thresh,sigma,2,fileout);
 sprintf(filename,"temp.fits");
 writefits(im,filename);
}
