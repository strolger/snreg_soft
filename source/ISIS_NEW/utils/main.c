#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

int main(int argc,char *argv[])
{
  int       i;          
  DATA_TYPE *im;
  char      filename[256],fileout[256];

   
     sprintf(fileout,"image.fits");
      for (i=1; i<argc; i++)
        switch(tolower(argv[i][1]))
        {
           case 'i': sprintf(filename, "%s", argv[++i]);
                        break; 
           case 'o': sprintf(fileout, "%s", argv[++i]);
        }


 im = readfits(filename);
 writefits(im,fileout);
}
