#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"
#include "short.h"

void writefits(DATA_TYPE *image, int w1, int w2, char *filename)
{
 int   i,j;
 float pf,*pfloat;
 char  cp[8],*cp2;
 FILE  *file;
 char  head[36][80]=
        {"SIMPLE  =                    T /",
         "BITPIX  =                  -32 /",
         "NAXIS   =                    2 /",
         "NAXIS1  =                      /",
         "NAXIS2  =                      /",
         "BSCALE  =                  1.0 /",
         "BZERO   =                  0.0 /",
         "BUNIT   = '        '           /",
         "OBJECT  = '                    /",
         "CDELT1  =                  1.0 /",
         "CDELT2  =                  1.0 /",
         "CRVAL1  =                  0.0 /",
         "CRVAL2  =                  0.0 /",     
         "END                             "};


 printf("size: %i %i %s\n",w1,w2,filename);
 file=fopen(filename,"wb");
 sprintf(&head[3][26],"%i /",w1);
 sprintf(&head[4][26],"%i /",w2);

   for (i=0; i<36; i++)
    for (j=0; j<80; j++) head[i][j] = (head[i][j])?head[i][j]:' ';
 
 fwrite (head, 2880, 1, file);        

 for (i=0; i<w1*w2; i++)
 {
  if(!swap_flag)
  {
   pf = (float)image[i];
   fwrite(&pf,4,1,file);
  }
  else
  {
   pf = (float)image[i];
   cp2 = (char *)&pf;
   cp[0]=cp2[3]; cp[1]=cp2[2]; cp[2]=cp2[1]; cp[3]=cp2[0];
   pfloat=(float *)&cp[0];
   fwrite(pfloat,4,1,file);
  }
 }

 fclose(file);
 return;
}
