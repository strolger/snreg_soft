#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

int get_psf(DATA_TYPE *psf,table_struct *psf_table,int xi,int yi,int npsf)
{
 DATA_TYPE *ps;
 int        i,psf_size,test;
 char       *filename,name[256];
 


 test=0;
 for(i=0;i<npsf;i++)
 {
  if(xi>=psf_table[i].xmin && xi<psf_table[i].xmax && yi>=psf_table[i].ymin && yi<psf_table[i].ymax)
   {
    filename=psf_table[i].filename;
    test=1;
   }
 }


 
 if(test)
 {
  sprintf(name,"psf.fits");
  ps=readfits(psf,filename,&psf_size,&psf_size,0);
  
  writefits(psf,psf_size,psf_size,name);
 }

 return test;
}
