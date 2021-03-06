#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"
#include "short.h"


DATA_TYPE *readfits(char *name)
{
 FILE      *ficin;
 int       test,nb,i,bitpix,pi,*pint;
 char      *header,s2[256],*cp,*cp2;
 short_type pi2,*pint2;
 DATA_TYPE *image;
 float     *pfloat,pf;
 double    pd,*pdouble,bscale,bzero;


   cp=(char *)malloc(16*sizeof(char));
   cp2=(char *)malloc(16*sizeof(char));

  if(!(ficin=fopen(name,"r"))) 
     {printf("Cannot Find File: %s\n", name); exit(0);}

   
   
   if(!(header=(char *)malloc(80)))
   {printf("Not enough memory to load Header\n"); exit(0);} 


 
 printf("name: %s\n", name);


 bscale=1.0;
 bzero=0.0;

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

     if(!(strncmp(&header[i*80],"BSCALE",6)))  
     {
      sscanf(&header[i*80+10],"%s",s2);    
      bscale=atof(s2);
     }
     if(!(strncmp(&header[i*80],"BZERO",5)))  
     {
      sscanf(&header[i*80+10],"%s",s2);    
      bzero=atof(s2);
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

  printf("width: %i height: %i bitpix: %i\n", width, height, bitpix);
  image=(DATA_TYPE *)malloc(width*height*sizeof(DATA_TYPE));

  if(!swap_flag)
  {
   if(bitpix == -32)
   {
      for(i=0;i<width*height;i++)
      {
       fread(&pf,sizeof(float),1,ficin);
       image[i]=pf*bscale+bzero;
      }
   }

   if(bitpix == 32)
   {
      for(i=0;i<width*height;i++)
      {
       fread(&pi,sizeof(int),1,ficin);
       image[i]=(DATA_TYPE)pi*bscale+bzero;
      }
   } 

   if(bitpix == 16)
   {
      for(i=0;i<width*height;i++)
      {
       fread(&pi2,sizeof(short_type),1,ficin);
       image[i]=(DATA_TYPE)pi2*bscale+bzero;
      }
   }
  
   if(bitpix == -64)
   {
      for(i=0;i<width*height;i++)
      {
       fread(&pd,sizeof(double),1,ficin);
       image[i]=(DATA_TYPE)pd*bscale+bzero;
      }
   }
   
  }
  else
  {
   if(bitpix == -32)
   {
    for(i=0;i<width*height;i++)
    {
     fread(cp2,4,1,ficin);  
     cp[0]=cp2[3]; cp[1]=cp2[2]; cp[2]=cp2[1]; cp[3]=cp2[0]; 
     pfloat=(float *)&cp[0]; image[i]=(DATA_TYPE) (*pfloat)*bscale+bzero;
    }
   }
   if(bitpix == 32)
   {
    for(i=0;i<width*height;i++)
    {
     fread(cp2,4,1,ficin);  
     cp[0]=cp2[3]; cp[1]=cp2[2]; cp[2]=cp2[1]; cp[3]=cp2[0]; 
     pint=(int *)&cp[0]; image[i]=(DATA_TYPE) (*pint)*bscale+bzero;
    }
   }
  if(bitpix == 16)
   {
    for(i=0;i<width*height;i++)
    {
     fread(cp2,2,1,ficin);  
     cp[0]=cp2[1]; cp[1]=cp2[0]; 
     pint2=(short_type *)&cp[0]; image[i]=(DATA_TYPE) (*pint2)*bscale+bzero;
    }
   }if(bitpix == -64)
   {
    for(i=0;i<width*height;i++)
    {
     fread(cp2,8,1,ficin);  
     cp[0]=cp2[7]; cp[7]=cp2[0]; cp[1]=cp2[6]; cp[6]=cp2[1];
     cp[2]=cp2[5]; cp[5]=cp2[2]; cp[3]=cp2[4]; cp[4]=cp2[3];
     
     pdouble=(double *)&cp[0]; image[i]=(DATA_TYPE) (*pdouble)*bscale+bzero;
    }
   }
  

  
  }
  fclose(ficin);

  sprintf(name,"toto");
  ficin=fopen(name,"w");
  fwrite(image,sizeof(DATA_TYPE)*width*height,1,ficin);
  fclose(ficin);


 return image;
}
