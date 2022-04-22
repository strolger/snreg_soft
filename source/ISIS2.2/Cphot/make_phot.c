#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void make_phot(DATA_TYPE *image,DATA_TYPE *image2,char *filename,char *conv_file)
{
  FILE         *phot_file;
  char         s[256],name[256],fname[256];
  double       x,y,*coeffs,dx,dy,dphot;
  int          xi,yi,ncomp2,i,ncomp_bg,nkernels,npsf,nto_phot;
  DATA_TYPE    *kernel,*psf,*new_psf,*Bpsf;
  table_struct *kernel_table,*psf_table;
  FILE         *table_file;

  if(!(phot_file=fopen(filename,"r"))) 
    {printf("Cannot find file:%s\n",filename); exit(0);}

  ncomp_kernel=0;
  mesh_size=2*half_mesh_size+1;
  conv_size=((int)radphot)*2+5;

  for(i=0;i<ngauss;i++) ncomp_kernel += ((deg_fixe[i]+1)*(deg_fixe[i]+2))/2;
  ncomp2=((deg_spatial+1)*(deg_spatial+2))/2;
  ncomp_bg = ((deg_bg+1)*(deg_bg+2))/2;
  ncomp_total = ncomp_kernel*ncomp2+ncomp_bg;
  printf("ncomp_total: %i\n",ncomp_total);

  coeffs=(double *)malloc((ncomp_total+1)*sizeof(double));
  kernel=(DATA_TYPE *)malloc(mesh_size*mesh_size*sizeof(DATA_TYPE));
  kernel_buff = (double *)malloc(ncomp_kernel*sizeof(double));
  psf=(DATA_TYPE *)malloc(stack_width*stack_width*sizeof(DATA_TYPE));
  new_psf=(DATA_TYPE *)malloc(conv_size*conv_size*sizeof(DATA_TYPE));
  Bpsf=(DATA_TYPE *)malloc(stack_width*stack_width*sizeof(DATA_TYPE));
  spline_allocate();


 sprintf(name,"kt_%s",&conv_file[5]);
 if(!(table_file=fopen(name,"r"))) 
    {printf("Cannot find file: %s\n",name); exit(0);}
 nkernels=0;
 while(!(feof(table_file)))
 {
  fgets(s,256,table_file);
  ++nkernels;
 }
 --nkernels;
 rewind(table_file);



 kernel_table=(table_struct *)malloc(nkernels*sizeof(table_struct));
 for(i=0;i<nkernels;i++)
 {
  fgets(s,256,table_file);
  sscanf(s,"%s %i %i %i %i %i %i",kernel_table[i].filename,&(kernel_table[i].xmin),&(kernel_table[i].ymin),&(kernel_table[i].xmax),&(kernel_table[i].ymax),&(kernel_table[i].offset_x),&(kernel_table[i].offset_y));
 }
 fclose(table_file);


 sprintf(name,"psf_table");
 if(!(table_file=fopen(name,"r"))) 
    {printf("Cannot find file: %s\n",name); exit(0);}
 npsf=0;
 while(!(feof(table_file)))
 {
  fgets(s,256,table_file);
  ++npsf;
 }
 --npsf;
 printf("npsf: %i\n",npsf);
 rewind(table_file);


 psf_table=(table_struct *)malloc(npsf*sizeof(table_struct));
 for(i=0;i<npsf;i++)
 {
  fgets(s,256,table_file);
  sscanf(s,"%s %i %i %i %i %i %i",psf_table[i].filename,&(psf_table[i].xmin),&(psf_table[i].ymin),&(psf_table[i].xmax),&(psf_table[i].ymax),&(psf_table[i].offset_x),&(psf_table[i].offset_y));
 }
 
 fclose(table_file);



 make_vectors();

  nto_phot=0;
  while(!(feof(phot_file)))
  {
    fgets(s,256,phot_file);
    ++nto_phot;
  }

  rewind(phot_file);

  for(i=0;i<nto_phot-1;i++)
  {
    fgets(s,256,phot_file);
    sscanf(s,"%lf %lf %i %i %s",&x,&y,&xi,&yi,fname);
    get_coeffs(coeffs,kernel_table,xi,yi,nkernels);
    make_kernel(kernel,xi,yi,coeffs);
    printf("x: %lf y: %lf\n",x,y);
    if(get_psf(psf,psf_table,x,y,npsf))
    {
     printf("psf done\n");
     dx=(double)xi-x;
     dy=(double)yi-y;
     Bspline(psf,Bpsf,dx,dy);
     printf("spline done\n");
     convolve(Bpsf,kernel,new_psf);
     printf("convolve done\n");
     phot(image,image2,new_psf,xi,yi,fname);
    }
  }

  sprintf(s,"new_psf.fits");
  writefits(new_psf,conv_size,conv_size,s);
  sprintf(s,"kernel.fits");
  writefits(kernel,mesh_size,mesh_size,s);
  

  fclose(phot_file);
  free(kernel);
  free(coeffs);
  free(kernel_buff);
  free(psf);
  free(new_psf);

 return;
}
