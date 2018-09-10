#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"


FILE  *read_header(char *name)
{
 FILE      *ficin;
 int       test,nb,i;
 char      *header,s2[256],*cp,*cp2;
 
 
  if(!(ficin=fopen(name,"r"))) 
     {printf("Cannot Find File: %s\n", name); exit(0);}

   cp=(char *)malloc(16*sizeof(char));
   cp2=(char *)malloc(16*sizeof(char));

   
   
   if(!(header=(char *)malloc(80)))
   {printf("Not enough memory to load Header\n"); exit(0);} 


 



 test=0; nb=1; 

 while(!test)
 {

   fread(header,80,1,ficin);

    for(i=0;i<1;i++)
    {
     if(!(strncmp(&header[i*80],"NAXIS1",6)))  
     {
      sscanf(&header[i*80+10],"%s",s2);    
      width=atoi(s2);
     }
    

     if(!(strncmp(&header[i*80],"NAXIS2",6)))  
     {
      sscanf(&header[i*80+10],"%s",s2);    
      height=atoi(s2);
     }


     if(!(strncmp(&header[i*80],"BITPIX",6)))  
     {
      sscanf(&header[i*80+10],"%s",s2);    
      bitpix=atoi(s2);
     }

    
   

     if(!(strncmp(&header[i*80],"END ",4)))  
     {
      test=1;
     } 


   if(!test) ++nb;
  
   }

  }


  nb=nb-((int)floor((double)(nb-1)/(double)36))*36;
  nb=36-nb;

   for(i=0;i<nb;i++)  fread(header,80,1,ficin);


  
 return ficin;
}
