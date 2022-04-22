#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"
#include "short.h"


void read_image(FILE *ficin,DATA_TYPE *image,int cx, int cy, int offset_header,int bitpix1,double bscale,double bzero)
{
 int       test,nb,i,j,k,pi,*pint,offset;
 char      *header,s2[256],*cp,*cp2;
 short_type pi2,*pint2;
 float     *pfloat,pf;
 double    pd,*pdouble;

  
   cp=(char *)malloc(16*sizeof(char));
   cp2=(char *)malloc(16*sizeof(char));

   
   bitpix=bitpix1;

 
 printf("bitpix: %i %i \n", bitpix,(int)fabs((double)bitpix)/8);
 printf("width: %i height: %i %i %i\n",width,height,width0,height0);
 offset=(cx+cy*width0)*(int)fabs((double)bitpix)/8;

 offset=(cx+cy*width0)*(int)fabs((double)bitpix)/8+offset_header;

 fseek(ficin,offset,SEEK_SET);

 test=0; nb=1; 

     offset=(width0-width)*(int)fabs((double)bitpix)/8;


  if(!swap_flag)
  {
   if(bitpix == -32)
   {
     for(k=0;k<height;k++) 
     {
      for(j=0;j<width;j++)
      {
       fread(&pf,sizeof(float),1,ficin);
       image[j+k*width]=pf;
      } 
      fseek(ficin,offset,SEEK_CUR);  
     }

   }

   if(bitpix == 32)
   {
     for(k=0;k<height;k++) 
     {
      for(j=0;j<width;j++)
      {
       fread(&pi,sizeof(int),1,ficin);
       image[j+k*width]=pi;
      }
      fseek(ficin,offset,SEEK_CUR);  
     }
   } 

   if(bitpix == 16)
   {
     for(k=0;k<height;k++) 
     {
      for(j=0;j<width;j++)
      {
       fread(&pi2,sizeof(short_type),1,ficin);
       image[j+k*width]=pi2;
      }
       fseek(ficin,offset,SEEK_CUR);  
     }
   }
  
   if(bitpix == -64)
   {
     for(k=0;k<height;k++) 
     {
      for(j=0;j<width;j++)
      {
       fread(&pd,sizeof(double),1,ficin);
       image[j+k*width]=pd;
      }
      fseek(ficin,offset,SEEK_CUR);  
     }
   }
   
  }
  else
  {
   if(bitpix == -32)
   {
     for(k=0;k<height;k++) 
     {
      for(j=0;j<width;j++)
      {
        fread(cp2,4,1,ficin);  
        cp[0]=cp2[3]; cp[1]=cp2[2]; cp[2]=cp2[1]; cp[3]=cp2[0]; 
        pfloat=(float *)&cp[0]; image[j+k*width]=(DATA_TYPE) (*pfloat);
      }
        fseek(ficin,offset,SEEK_CUR);  
     }
   }

   if(bitpix == 32)
   {
    for(k=0;k<height;k++) 
     {
      for(j=0;j<width;j++)
      {
       fread(cp2,4,1,ficin);  
       cp[0]=cp2[3]; cp[1]=cp2[2]; cp[2]=cp2[1]; cp[3]=cp2[0]; 
        pint=(int *)&cp[0]; image[j+k*width]=(DATA_TYPE) (*pint);
      } 
       fseek(ficin,offset,SEEK_CUR);  

    }       
   }

  if(bitpix == 16)
   {
    for(k=0;k<height;k++) 
    {
     for(j=0;j<width;j++)
     {
      fread(cp2,2,1,ficin);  
      cp[0]=cp2[1]; cp[1]=cp2[0]; 
      pint2=(short_type *)&cp[0]; image[j+k*width]=(DATA_TYPE) (*pint2);
     }
      fseek(ficin,offset,SEEK_CUR);
    }
   }

   if(bitpix == -64)
   {
    for(k=0;k<height;k++) 
    {
     for(j=0;j<width;j++)
     {
      fread(cp2,8,1,ficin);  
      cp[0]=cp2[7]; cp[7]=cp2[0]; cp[1]=cp2[6]; cp[6]=cp2[1];
      cp[2]=cp2[5]; cp[5]=cp2[2]; cp[3]=cp2[4]; cp[4]=cp2[3];
      pdouble=(double *)&cp[0]; image[j+k*width]=(DATA_TYPE) (*pdouble);
     }
      fseek(ficin,offset,SEEK_CUR);
    }
   }
  

  
  }
  printf("offset: %i %i %lf\n", offset, bitpix, image[width*85+135]);



 return;
}
