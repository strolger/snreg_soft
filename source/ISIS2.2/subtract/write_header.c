#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"

void write_header(FILE *file)
{
 int   i,j;
 float fl;
 char  head[36][80]=
        {"SIMPLE  =                    T /",
         "BITPIX  =                  -32 /",
         "NAXIS   =                    2 /",
         "NAXIS1  =                      /",
         "NAXIS2  =                      /",
         "BSCALE  =                  1.0 /",
         "BZERO   =                  0.0 /",
         "BUNIT   = '        '           /",
         "OBJECT  = '        '           /",
         "CDELT1  =                  1.0 /",
         "CDELT2  =                  1.0 /",
         "CRVAL1  =                  0.0 /",
         "CRVAL2  =                  0.0 /",     
         "END                             "};



 sprintf(&head[3][26],"%i /",width0);
 sprintf(&head[4][26],"%i /",height0);

   for (i=0; i<36; i++)
    for (j=0; j<80; j++) head[i][j] = (head[i][j])?head[i][j]:' ';
 
 fl=0.0;
 fwrite(head, 2880, 1, file);  
 for(i=0;i<width0*height0;i++) fwrite(&fl,4,1,file); 
 fclose(file);     

 return;
}
