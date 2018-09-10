#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>


#include "types.h"

int main(int argc,char *argv[])
{
  int       i,j,nargc,infile,xi,yi,i_max,j_max,mesh,var_number,iv,
            line_buffer,nbuffers,mline_buffer,reste;          
  DATA_TYPE **stack,pixmax,pix;
  char      filename[256],fileout[256],filein[256],s[1000],
            file_phot[256],fileconf[256],med_flag;
  FILE      **file_list,*file2,*outfile;
  double    x,y;
  bg_struct *bgnorm;

     infile=0;
     mesh=2;
     nargc=argc;
     line_buffer=256;
     old_buffer=0;
     deg_bg=1;
     j0b=0;
     med_flag=1;

     sprintf(fileout,"image.fits");
     sprintf(file_phot,"phot.data");
     sprintf(fileconf,"phot_config");

      for (i=1; i<argc; i++)
        {
           if(!(strncmp(argv[i],"-i",2))) 
	   {
             sprintf(filein,"%s",argv[i+1]);
             nargc -= 2;
             infile=1;
           }
           if(!(strncmp(argv[i],"-o",2))) 
	   {
             sprintf(fileout,"%s",argv[i+1]);
             nargc -= 2;
           }
           if(!(strncmp(argv[i],"-p",2))) 
	   {
             sprintf(file_phot,"%s",argv[i+1]);
             nargc -= 2;
           }

           if(!(strncmp(argv[i],"-c",2))) 
	   {
             sprintf(fileconf,"%s",argv[i+1]);
             nargc -= 2;
           }
           if(!(strncmp(argv[i],"-m",2))) 
	   {
             med_flag=0;
             nargc -= 1;
           }
         }


 
 file_list = (FILE **)malloc(nargc*sizeof(FILE *));
 ndex = (int *)malloc(nargc*sizeof(int));
 buffer = (DATA_TYPE *)malloc(nargc*sizeof(DATA_TYPE));

 printf("Nargc: %i\n", nargc);
 read_config(fileconf);
 printf("sat: %f\n", saturation);
 allocate();
 

 

 for(i=1;i<nargc;i++)
 {
  sprintf(filename,"%s",argv[i]);
  printf("i: %i filename: %s\n", i,filename);
  file_list[i-1]=read_header(filename);
 }

 
 nbuffers=floor((double)height/(double)line_buffer);
 line_buffer=(int)floor((double)height/(double)nbuffers);
 reste = height-line_buffer*nbuffers;
 mline_buffer = line_buffer;

 
 printf("nbuffers: %i %i %i %i %i\n", nbuffers,nargc-1,line_buffer,mline_buffer,reste); 
 
 
 stack = (DATA_TYPE **)malloc(nargc*sizeof(DATA_TYPE *));
 for(i=0;i<nargc;i++)  
 stack[i] = (DATA_TYPE *)malloc(width*mline_buffer*sizeof(DATA_TYPE));
 ref=(DATA_TYPE *)malloc(width*height*sizeof(DATA_TYPE));
 bgnorm=(bg_struct *)malloc(nargc*sizeof(bg_struct));
 
 

 for(j=0;j<nbuffers;j++)
 {
  for(i=0;i<nargc-1;i++)
  {
   readfits(file_list[i],line_buffer,stack[i]);
  
  }
   if(med_flag)  clean_stack(stack,nargc-1,line_buffer);
   else          moy(stack,nargc-1,line_buffer);
 }

 
  for(i=0;i<nargc-1;i++)
  {
   readfits(file_list[i],mline_buffer,stack[i]);
  }


if(reste>0)
{
 if(med_flag)  clean_stack(stack,nargc-1,reste);
 else          moy(stack,nargc-1,reste); 
}

/*if(infile)
{
 var_number=0; 
 file2=fopen(filein,"r");
 outfile=fopen(file_phot,"w");
 while(!feof(file2))
 {
     fgets(s,1000,file2);
     ++var_number;
 }

 rewind(file2);

 for(iv=0;iv<var_number-1;iv++)
 {
   fgets(s,1000,file2);
   sscanf(s,"%i %i",&xi,&yi);
   pixmax=fabs(stack[xi+width*yi]);
   i_max=xi;
   j_max=yi;
   for(i=xi-mesh;i<=xi+mesh;i++)
    for(j=yi-mesh;j<=yi+mesh;j++)
    {
      pix=fabs(image[i+width*j]);
      if(pix>pixmax) {i_max=i; j_max=j; pixmax=pix;}
    }
   cent(i_max,j_max,&x,&y,stack);
 
 sprintf(s,"lc%i.data",iv);
 fprintf(outfile,"%f %f %i %i %s\n",x,y,i_max,j_max,s);
 }

  fclose(file2);
  fclose(outfile);
}*/


 writefits(ref,fileout);
 

 for(i=0;i<nargc;i++) free(stack[i]);
 free(stack);
 free(ref);
 free(buffer); 
 free(ndex);
 free(bgnorm);
}
