#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"
#include "short.h"

void readfits(FILE *ficin,int line_buffer,DATA_TYPE *image,int mesh_smooth)
{
 int        test,nb,i,pi,*pint;
 char       *header,s2[256],*cp,*cp2;
 short_type pi2,*pint2;
 float      *pfloat,pf;
 double     pd,*pdouble;


   cp=(char *)malloc(16*sizeof(char));
   cp2=(char *)malloc(16*sizeof(char));

 if(!swap_flag)
  {
   if(bitpix == -32)
   {
      for(i=0;i<width*line_buffer;i++)
      {
       fread(&pf,sizeof(float),1,ficin);
       image[i]=pf;
      }
   }

   if(bitpix == 32)
   {
      for(i=0;i<width*line_buffer;i++)
      {
       fread(&pi,sizeof(int),1,ficin);
       image[i]=pi;
      }
   } 

   if(bitpix == 16)
   {
      for(i=0;i<width*line_buffer;i++)
      {
       fread(&pi2,sizeof(short_type),1,ficin);
       image[i]=pi2;
      }
   }
  
   if(bitpix == -64)
   {
      for(i=0;i<width*line_buffer;i++)
      {
       fread(&pd,sizeof(double),1,ficin);
       image[i]=pd;
      }
   }
   
  }
  else
  {
   if(bitpix == -32)
   {
    for(i=0;i<width*line_buffer;i++)
    {
     fread(cp2,4,1,ficin);  
     cp[0]=cp2[3]; cp[1]=cp2[2]; cp[2]=cp2[1]; cp[3]=cp2[0]; 
     pfloat=(float *)&cp[0]; image[i]=(DATA_TYPE) (*pfloat);
    }
   }
   if(bitpix == 32)
   {
    for(i=0;i<width*line_buffer;i++)
    {
     fread(cp2,4,1,ficin);  
     cp[0]=cp2[3]; cp[1]=cp2[2]; cp[2]=cp2[1]; cp[3]=cp2[0]; 
     pint=(int *)&cp[0]; image[i]=(DATA_TYPE) (*pint);
    }
   }
  if(bitpix == 16)
   {
    for(i=0;i<width*line_buffer;i++)
    {
     fread(cp2,2,1,ficin);  
     cp[0]=cp2[1]; cp[1]=cp2[0]; 
     pint2=(short_type *)&cp[0]; image[i]=(DATA_TYPE) (*pint2);
    }
   }if(bitpix == -64)
   {
    for(i=0;i<width*line_buffer;i++)
    {
     fread(cp2,8,1,ficin);  
     cp[0]=cp2[7]; cp[7]=cp2[0]; cp[1]=cp2[6]; cp[6]=cp2[1];
     cp[2]=cp2[5]; cp[5]=cp2[2]; cp[3]=cp2[4]; cp[4]=cp2[3];
     
     pdouble=(double *)&cp[0]; image[i]=(DATA_TYPE) (*pdouble);
    }
   }
  }
 i=(long)(-mesh_smooth*width*fabs((double)(bitpix/8)))*2;
 fseek(ficin,i,SEEK_CUR);


 return;
}
