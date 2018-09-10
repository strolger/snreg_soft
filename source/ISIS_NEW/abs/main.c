#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>


#include "types.h"

int main(int argc,char *argv[])
{
  int       i,j,nargc,infile,xi,yi,i_max,j_max,mesh,var_number,iv,
            line_buffer,nbuffers,mline_buffer,reste,ist,line_buffer_smooth,
            mesh_smooth,kl,delta;          
  DATA_TYPE **stack,**stack2,pixmax,pix,*imtemp,*ref,*ref2;
  char      filename[256],fileout[256],filein[256],s[1000],
            file_phot[256],fileconf[256],med_flag,pfile[256];
  FILE      **file_list,**file_list2,*file2,*outfile;
  double    x,y;
  bg_struct *bgnorm;

     infile=0;
     mesh=2;
     nargc=argc;
     line_buffer=64;
     old_buffer=0;
     deg_bg=1;
     j0b=0;
     med_flag=1;
     mesh_smooth=3;
     n_reject=2;
     sprintf(fileout,"image.fits");
     sprintf(file_phot,"phot.data");
     sprintf(fileconf,"phot_config");

       for (i=1; i<argc; i++)
         {
          
           if(!(strncmp(argv[i],"-o",2))) 
	   {
             sprintf(fileout,"%s",argv[i+1]);
             nargc -= 2;
           }
           if(!(strncmp(argv[i],"-t",2))) 
	   {
             n_reject=atoi(argv[i+1]);
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
           if(!(strncmp(argv[i],"-s",2))) 
	   {
             mesh_smooth=atoi(argv[i+1]);
             nargc -= 2;
           }

         }

 
 file_list  = (FILE **)malloc(nargc*sizeof(FILE *));
 file_list2 = (FILE **)malloc(nargc*sizeof(FILE *));
 ndex = (int *)malloc(nargc*sizeof(int));
 buffer = (DATA_TYPE *)malloc(nargc*sizeof(DATA_TYPE));

 printf("Nargc: %i\n", nargc);
 read_config(fileconf);
 printf("sat %f mesh: %i rej: %i\n", saturation,mesh_smooth,n_reject);
 allocate();
 

 printf("medflag: %i \n",med_flag);

 for(i=1;i<nargc;i++)
 {
  sprintf(filename,"%s",argv[i]);
  printf("i: %i filename: %s\n", i,filename);
  file_list[i-1]=read_header(filename);
  sprintf(pfile,"interp_%s",argv[i]+5);
  file_list2[i-1]=read_header(pfile);
 }

 
 nbuffers=floor((double)height/(double)line_buffer);
 line_buffer=(int)floor((double)height/(double)nbuffers);
 reste = height-line_buffer*nbuffers;
 mline_buffer = line_buffer-mesh_smooth;
 line_buffer_smooth=line_buffer+2*mesh_smooth;
 
 printf("nbuffers: %i %i %i %i %i\n", nbuffers,nargc-1,line_buffer,mline_buffer,reste); 
 
 imtemp=(DATA_TYPE *)malloc(height*line_buffer_smooth*sizeof(DATA_TYPE));
 temp_smooth=(DATA_TYPE *)malloc(height*line_buffer_smooth*sizeof(DATA_TYPE));

 stack = (DATA_TYPE **)malloc(nargc*sizeof(DATA_TYPE *));
 stack2 = (DATA_TYPE **)malloc(nargc*sizeof(DATA_TYPE *));

 for(i=0;i<nargc;i++){ 
 stack[i] = (DATA_TYPE *)malloc(width*line_buffer_smooth*sizeof(DATA_TYPE));
 stack2[i] = (DATA_TYPE *)malloc(width*line_buffer_smooth*sizeof(DATA_TYPE));}

 ref=(DATA_TYPE *)malloc(height*width*sizeof(DATA_TYPE));
 ref2=(DATA_TYPE *)malloc(height*width*sizeof(DATA_TYPE));


 bgnorm=(bg_struct *)malloc(nargc*sizeof(bg_struct));
 
 old_buffer=0;

 for(j=0;j<nbuffers;j++)
 {
  for(i=0;i<nargc-1;i++)
  {
   kl=line_buffer_smooth;
   if(j<1) kl = line_buffer_smooth-mesh_smooth;
   if(kl+old_buffer>height) kl=height-old_buffer; 

   readfits(file_list[i],kl,stack[i],mesh_smooth);
   readfits(file_list2[i],kl,imtemp,mesh_smooth);
   smooth(stack[i],kl,mesh_smooth);
     smooth(imtemp,kl,mesh_smooth);
   delta=0;
   if(j>0) delta=width*mesh_smooth;
   for(ist=0;ist<width*line_buffer;ist++)
   {   
     stack2[i][ist]=fabs(stack[i][ist+delta]);
    stack[i][ist]=fabs(stack[i][ist+delta])/sqrt(fabs(imtemp[ist+delta])+0.001);    
   }
  }

   if(med_flag)
   {
     clean_stack(stack,ref,nargc-1,line_buffer,0);
     clean_stack(stack2,ref2,nargc-1,line_buffer,1);
   }

   else
   {
      moy(stack,ref,nargc-1,line_buffer,0);  
      moy(stack2,ref2,nargc-1,line_buffer,1);
   }


 }

 

 

for(i=0;i<nargc-1;i++)
{
if(reste>0)
{
 if(med_flag)
   {
     readfits(file_list[i],reste+mesh_smooth,stack[i],0); 
     readfits(file_list2[i],reste+mesh_smooth,imtemp,0);
     smooth(stack[i],kl,mesh_smooth);
     smooth(imtemp,kl,mesh_smooth);
   delta=height*mesh_smooth;
   for(ist=0;ist<height*reste;ist++)
   {   
    stack2[i][ist]=fabs(stack[i][ist+delta]);
    stack[i][ist]=fabs(stack[i][ist+delta])/sqrt(fabs(imtemp[ist+delta])+0.001);    
   }
     clean_stack(stack,ref,nargc-1,reste,0);
     clean_stack(stack2,ref2,nargc-1,reste,0);
   }

   else
   {
      moy(stack,ref,nargc-1,reste,0);  
      moy(stack2,ref2,nargc-1,reste,0);
   } 
}
}

printf("writing\n");
 writefits(ref,fileout);
 sprintf(fileout,"abs.fits");
 writefits(ref2,fileout);

printf("writing done\n");


for(i=0;i<nargc;i++){free(stack[i]); free(stack2[i]);}
 free(stack);
 free(stack2);
 free(ref);
 free(ref2);
 free(buffer); 
 free(ndex);
 free(bgnorm);
}
