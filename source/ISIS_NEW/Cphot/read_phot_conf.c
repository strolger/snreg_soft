#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void read_phot_conf(char *filename)
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
   
   if(!(strncmp(s2,"radphot",7)))
   {
    sscanf(s,"%s %lf",s2,&df);
    radphot=df; 
   }
   if(!(strncmp(s2,"psf_width",9)))
   {
    sscanf(s,"%s %i",s2,&stack_width);
   }
   
    if(!(strncmp(s2,"rad_aper",8)))
   {
    sscanf(s,"%s %lf",s2,&df);
    rad_aper=df; 
   }

   if(!(strncmp(s2,"nb_adu_el",9)))
   {
    sscanf(s,"%s %i",s2,&adu_el);
   }
 }

 fclose(file);

 return;
}
