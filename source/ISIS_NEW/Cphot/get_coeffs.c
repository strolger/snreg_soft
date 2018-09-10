#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"


void get_coeffs(double *coeffs,table_struct *kernel_table,int xi,int yi,int nkernels)
{
 int i;
 FILE *file;
 char *filename;


 for(i=0;i<nkernels;i++)
 {
  if(xi>=kernel_table[i].xmin && xi<kernel_table[i].xmax && yi>=kernel_table[i].ymin && yi<kernel_table[i].ymax)
   {
    filename=kernel_table[i].filename;
   }
 }

  printf("kernel_file: %s\n",filename);

  file=fopen(filename,"rb");
  fread(coeffs,(ncomp_total+1)*sizeof(double),1,file);

  fclose(file);
  
 return;
}
