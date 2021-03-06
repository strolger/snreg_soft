#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"

void read_image(char *name)
{
 FILE      *ficin;
 int       test,nb,i,bitpix,pi,*pint;
 char      *header,s2[256],*cp,*cp2;
 short_type pi2,*pint2;
 DATA_TYPE *image;
 float     *pfloat,pf;
 double    pd,*pdouble,bscale,bzero;

  if(!(ficin=fopen(name,"r"))) 
     {printf("Cannot Find File: %s\n", name); exit(0);}

   cp=(char *)malloc(16*sizeof(char));
   cp2=(char *)malloc(16*sizeof(char));

   

 
 printf("name: %s\n", name);



 test=0; nb=1; 

 

  if(!swap_flag)
  {
   if(bitpix == -32)
   {
      for(i=0;i<width*height;i++)
      {
       fread(&pf,sizeof(float),1,ficin);
       sub_image[i]=pf*bscale+bzero;
      }
   }

   if(bitpix == 32)
   {
      for(i=0;i<width*height;i++)
      {
       fread(&pi,sizeof(int),1,ficin);
       sub_image[i]=(DATA_TYPE)pi*bscale+bzero;
      }
   } 

   if(bitpix == 16)
   {
      for(i=0;i<width*height;i++)
      {
       fread(&pi2,sizeof(short_type),1,ficin);
       sub_image[i]=(DATA_TYPE)pi2*bscale+bzero;
      }
   }
  
   if(bitpix == -64)
   {
      for(i=0;i<width*height;i++)
      {
       fread(&pd,sizeof(double),1,ficin);
       sub_image[i]=pd*bscale+bzero;
      }
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
  }
  else
  {
   if(bitpix == -32)
   {
    for(i=0;i<width*height;i++)
    {
     fread(cp2,4,1,ficin);  
     cp[0]=cp2[3]; cp[1]=cp2[2]; cp[2]=cp2[1]; cp[3]=cp2[0]; 
     pfloat=(float *)&cp[0]; 
     sub_image[i]=(DATA_TYPE) (*pfloat)*bscale+bzero;
    }
   }
   if(bitpix == 32)
   {
    for(i=0;i<width*height;i++)
    {
     fread(cp2,4,1,ficin);  
     cp[0]=cp2[3]; cp[1]=cp2[2]; cp[2]=cp2[1]; cp[3]=cp2[0]; 
     pint=(int *)&cp[0]; 
     sub_image[i]=(DATA_TYPE) (*pint)*bscale+bzero;
    }
   }
  if(bitpix == 16)
   {
    for(i=0;i<width*height;i++)
    {
     fread(cp2,2,1,ficin);  
     cp[0]=cp2[1]; cp[1]=cp2[0]; 
     pint2=(short_type *)&cp[0]; 
     sub_image[i]=(DATA_TYPE) (*pint2)*bscale+bzero;
    }
   }if(bitpix == -64)
   {
    for(i=0;i<width*height;i++)
    {
     fread(cp2,8,1,ficin);  
     cp[0]=cp2[7]; cp[7]=cp2[0]; cp[1]=cp2[6]; cp[6]=cp2[1];
     cp[2]=cp2[5]; cp[5]=cp2[2]; cp[3]=cp2[4]; cp[4]=cp2[3];
     
     pdouble=(double *)&cp[0]; 
     sub_image[i]=(DATA_TYPE) (*pdouble)*bscale+bzero;
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
