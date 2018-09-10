#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

int main(int argc,char *argv[])
{
 FILE      *file1,*file,*outf,*file_psf,*psf_table;
 char      filename[256],config_file[256],outfile[256],fout[256],psf_file[256],           psf_table_name[256];
 int       i,j,iload_x,iload_y,offset,cx,cy,mesh,nj,smesh,frame_count;
 DATA_TYPE *image,*psf_map,thresh;
 object    *obj;
 double    bscale,bzero;


 /*swap_flag=0;*/
       sprintf(config_file, "default_config");
       sprintf(outfile, "psf.fits");
       sprintf(fout, "ref.data");
       sprintf(psf_table_name, "psf_table");
       outf=fopen(fout,"w");
       psf_table= fopen(psf_table_name,"w");      

       for (i=1; i<argc; i++)
        {
          if(!(strncmp(argv[i],"-c",2))) sprintf(config_file,"%s",argv[i+1]);
          if(!(strncmp(argv[i],"-i",2))) sprintf(filename,"%s",argv[i+1]);
          if(!(strncmp(argv[i],"-o",2))) sprintf(outfile,"%s",argv[i+1]);
          /*if(!(strncmp(argv[i],"-s",2))) swap_flag=1;*/
        }


 

 read_config(config_file);
 stack = (DATA_TYPE **)malloc(NSTARS*sizeof(DATA_TYPE *));
 for(i=0;i<NSTARS;i++) 
 stack[i]=(DATA_TYPE *)malloc(stack_width*stack_width*sizeof(DATA_TYPE));
 psf = (DATA_TYPE *)malloc(stack_width*stack_width*sizeof(DATA_TYPE));
 psf_map = (DATA_TYPE *)malloc(stack_width*stack_width*nload_x*nload_y*sizeof(DATA_TYPE));
 file1=read_header(filename,&offset,&bscale,&bzero); 

 thresh=0.0;
 frame_count=0;
 smesh=stack_width/2;
 nstars_buffer=NSTARS;
 width0=width;
 height0=height;
 width = floor((double)(width0+2*(nload_x-1)*smesh)/(double)nload_x);
 height = floor((double)(height0+2*(nload_y-1)*smesh)/(double)nload_y);

printf("width: %i height: %i %i %i %i\n",width0,height0,width,height,smesh);

 image=(DATA_TYPE *)malloc(width*height*sizeof(DATA_TYPE));
 init_bg(rad1bg,rad2bg);
 mesh=2;

 nj=0;
 for(iload_y=0;iload_y<nload_y;iload_y++)
  for(iload_x=0;iload_x<nload_x;iload_x++)
  {
   cx=(iload_x*(width-2*smesh));
   cy=(iload_y*(height-2*smesh));
   read_image(file1,image,cx,cy,offset,bitpix,bscale,bzero);
   obj=max(image,thresh,mesh,rad1bg,rad2bg);
   sprintf(psf_file,"psf_file%i.fits",frame_count);
 
  fprintf(psf_table,"%s %5i %5i %5i %5i %5i %5i\n",psf_file,cx+smesh,cy+smesh,width-smesh+cx,height-smesh+cy,cx,cy); 

   make_psf(obj,image,nnstars,cx,cy,frame_count++,outf,psf_file);
  
  

   for(i=0;i<stack_width;i++)
    for(j=0;j<stack_width;j++)
     psf_map[i+iload_x*stack_width+nload_x*stack_width*(j+iload_y*stack_width) ] = psf[i+stack_width*j];
  }
 
   writefits(psf_map,stack_width*nload_x,stack_width*nload_y,outfile);
   fclose(outf);
   fclose(psf_table);

}
