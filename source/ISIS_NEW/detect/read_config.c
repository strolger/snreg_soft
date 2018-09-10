#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void read_config(char *filename)
{
 FILE *file;
 double df;
 char   s[256],s2[256];


if(!(file=fopen(filename,"r"))) 
     {printf("Oops...cannot find any configuration file...Goodbye !\n"); exit(0)
;}

while(feof(file) == 0)
 {
   fgets(s,256,file);
   sscanf(s,"%s",s2);
   
   if(!(strncmp(s2,"saturation",10)))
   {
    sscanf(s,"%s %lf",s2,&df);
    saturation=df; 
   }

   if(!(strncmp(s2,"nb_adu_el",9)))
   {
    sscanf(s,"%s %lf",s2,&df);
    adu_el=df; 
   }
 }

 fclose(file);

 return;
}
