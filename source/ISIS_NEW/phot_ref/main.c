#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

int main(int argc,char *argv[])
{
  int       i,nn;          
  DATA_TYPE *im,rad,rad1,rad2,thresh,sigma;
  char      filename[256],fileout[256],fileconf[256];
  object    *objlist;
  double    ax,ay,teta;
  FILE      *file;


    
     sprintf(fileout,"psf.data");
     sprintf(fileconf,"config_default");
      for (i=1; i<argc; i++)
        switch(tolower(argv[i][1]))
        { 
           case 'c': sprintf(fileconf, "%s", argv[++i]);
                        break;
           case 'i': sprintf(filename, "%s", argv[++i]);
                        break; 
           case 'o': sprintf(fileout, "%s", argv[++i]);
                        break; 
        }

 
 ax0=ay0= 0.5;
 read_config(fileconf);
 rad=radphot;
 rad1=rad1bg;
 rad2=rad2bg;
 
 
 stack = (DATA_TYPE **)malloc(NSTARS*sizeof(DATA_TYPE *));
 for(i=0;i<NSTARS;i++) 
 stack[i]=(DATA_TYPE *)malloc(stack_width*stack_width*sizeof(DATA_TYPE));
 psf = (DATA_TYPE *)malloc(stack_width*stack_width*sizeof(DATA_TYPE));
 
 im = readfits(filename);
 thresh=get_thresh(im,&sigma);

 objlist=max(im,thresh+10.0*sigma,2,rad1,rad2);
 aperphot(objlist,im,rad,&nn);
 make_psf(objlist,im,nn);
 norm_psf();
  

 get_gaussian(&ax,&ay,&teta);
 
 make_phot(objlist,im,ax,ay,teta,nn,rad,fileout);

 
 printf("mean: %lf ax: %lf ay: %lf teta: %lf\n",2.0*sqrt(log(2.0)*(0.5/ax+0.5/ay)),2.0*sqrt(log(2.0)/ax),2.0*sqrt(log(2.0)/ay),teta/sin(1.0)*2.0*180.0);

 file=fopen(fileout,"wb");
 fwrite(psf,stack_width*stack_width*4,1,file);
 fclose(file);

}
