#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include <string.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

int main(int argc,char *argv[])
{
  int         i,j,k,iload_x,iload_y,nload_x,nload_y,cx,cy,frame_number,
              offset1,offset2,i1,i2,ii,ncomp1,ncomp2,bitpix1,bitpix2
              ,ncomp_total0,pos,cx2,cy2;
  double      **matrix,d,background,bscale1,bzero1,bscale2,bzero2;
  char        s[256],s2[256],config_file[256],*map,check,kernel_name[256],
              table_name[256],*spt;
  FILE        *file,*file0,*file1,*file2,*kernel_table,*coeff_file;


 
 

 
  sprintf(config_file, "default_config");
  pos=get_index(argv[2],"/");
  spt=&argv[2][pos+1];
  if(strlen(spt) > 8) sprintf(kernel_name, "kt_%s",&spt[7]);
  else sprintf(kernel_name, "kt_%s",spt);
  kernel_table=fopen(kernel_name,"w");
      
        for (i=1; i<argc; i++)
        {
           if(!(strncmp(argv[i],"-c",2))) sprintf(config_file,"%s",argv[i+1]);
        }

 printf("Reading Image1\n");
 file1=read_header(argv[1],&offset1,&bscale1,&bzero1); bitpix1=bitpix;
 file2=read_header(argv[2],&offset2,&bscale2,&bzero2); bitpix2=bitpix;
 printf("offsets: %i %i\n",offset1,offset2);

 frame_number=0;
 width0  = width;
 height0 = height;
 sprintf(s,"conv.fits");
 file=fopen(s,"wb");
 write_header(file);
 sprintf(s2,"conv0.fits");
 file0=fopen(s2,"wb");
 write_header(file0);

 file=fopen(s,"r+");
 file0=fopen(s2,"r+");


 
 read_config(config_file); 


 


 mesh_size=half_mesh_size*2+1;
 stamp_size=half_stamp_size*2+1;
 nstamps=nstamps_x*nstamps_y;
 printf("ss: %i %i %i %i\n",sub_width,sub_height,width,height);

 sub_width = floor((double)(width-mesh_size+1)/(double)sub_width+(double)mesh_size-1.0);
 sub_height = floor((double)(height-mesh_size+1)/(double)sub_height+(double)mesh_size-1.0);

 printf("ss: %i %i %i %i\n",sub_width,sub_height,width,height);
 


 printf("nsx: %i nsy: %i\n", nstamps_x, nstamps_y);
 printf("mesh_size: %i stam_size: %i\n", mesh_size, stamp_size);
 printf("ng: %i %i\n", ngauss, deg_fixe[0]);
 printf("sg: %lf\n",sigma_gauss[0]);
 printf("sat: %lf pix_min: %lf max_stamp %i\n",SATURATION,PIX_MIN,max_stamp_thresh);
 printf("deg_spatial: %i deg_bg: %i \n", deg_spatial, deg_bg);
 printf("subs: %i %i\n", sub_height,sub_width);
 
 

 ncomp_kernel=0;
 for(i=0;i<ngauss;i++) ncomp_kernel += ((deg_fixe[i]+1)*(deg_fixe[i]+2))/2;
 ncomp1=ncomp_kernel;
 ncomp2=((deg_spatial+1)*(deg_spatial+2))/2;
 ncomp_kernel_xy = ncomp_kernel*ncomp2;
 ncomp_bg = ((deg_bg+1)*(deg_bg+2))/2;
 ncomp_total0=ncomp_total = ncomp_kernel_xy+ncomp_bg;
 ncomp_background =  (ncomp_kernel-1)*ncomp2+1;




 width   = sub_width;
 height  = sub_height;
 nload_x = width0/(sub_width-mesh_size+1); 
 nload_y = height0/(sub_height-mesh_size+1);

   if(!(stamps=(stamp_struct *)malloc(nstamps*sizeof(stamp_struct))))
    {
     printf("Cannot Allocate Stamp List\n"); 
     exit(0);
    }

 printf("Allocate: %i\n",nstamps);
 stamp_number=nstamps;
 allocate_vectors();
 printf("Getting Stamps\n");

 filter_x=(DATA_TYPE *)malloc(mesh_size*sizeof(DATA_TYPE));
 filter_y=(DATA_TYPE *)malloc(mesh_size*sizeof(DATA_TYPE));
 temp=(DATA_TYPE *)malloc((stamp_size+mesh_size)*stamp_size*sizeof(DATA_TYPE));
 scprod = (double *)malloc((ncomp_total+1)*sizeof(double));
 indx  = (int *)malloc((ncomp_total+1+100)*sizeof(int));
 kernel_vec = (DATA_TYPE **)malloc(ncomp_kernel*sizeof(DATA_TYPE *));
 kernel_coeffs = (double *)malloc(ncomp_kernel*sizeof(double));
 kernel = (double *)malloc(mesh_size*mesh_size*sizeof(double));
 if(!(sub_ref = (DATA_TYPE *)malloc(sub_width*sub_height*sizeof(DATA_TYPE))))
    {
     printf("Cannot Allocate  Ref Image\n"); 
     exit(0);
    };

 if(!(sub_image = (DATA_TYPE *)malloc(sub_width*sub_height*sizeof(DATA_TYPE))))
    {
     printf("Cannot Allocate  Image\n"); 
     exit(0);
    };
 if(!(conv_image=(DATA_TYPE *)malloc(width*height*sizeof(DATA_TYPE))))
    {
     printf("Cannot Allocate  Conv Image\n"); 
     exit(0);
    };
 check_stack=(DATA_TYPE *)malloc(nstamps*sizeof(DATA_TYPE));



 
for(iload_y=0;iload_y<nload_y;iload_y++)
 for(iload_x=0;iload_x<nload_x;iload_x++)
 {
  cx=(iload_x*(sub_width-mesh_size+1));
  cy=(iload_y*(sub_height-mesh_size+1));
  cx2=(iload_x*sub_width);
  cy2=(iload_y*sub_height);
 

 read_image(file1,sub_ref,cx,cy,offset1,bitpix1,bscale1,bzero1);
 read_image(file2,sub_image,cx,cy,offset2,bitpix2,bscale2,bzero2);

 printf("sub_ref: %lf %lf\n",sub_ref[1000],sub_image[1000]);

 for(i=0;i<width*height;i++) conv_image[i]=0.0;
 printf("Making defect map\n");
 map=defect_map(sub_ref,sub_image);
printf("Looking for Stamps\n");
 get_stamps(sub_ref,map); 
 get_stamp_vectors(sub_ref,sub_image);


printf("Building Zero Order Matrix\n");


 build_matrix0();
 build_scprod0(sub_image);

printf("Checking Stamps\n");
 check_stamps();
  
printf("Expanding Matrix\n");

 matrix=build_matrix();
printf("Expanding Matrix Done\n");
 build_scprod(sub_image);

  ncomp_total = (ncomp1-1)*ncomp2+1+ncomp_bg;
  ludcmp(matrix,ncomp_total,indx,&d);
  lubksb(matrix,ncomp_total,indx,scprod);

     

 


 printf("Checking again\n");
 check=check_again(scprod,sub_image);


 printf("Checking Done\n");
  if(check)
   {
    printf("Expanding Matrix\n");
    matrix=build_matrix();
    printf("Expanding Matrix Done\n");

    build_scprod(sub_image);
    printf("Expanding Matrix Done\n");
    ludcmp(matrix,ncomp_total,indx,&d);
    lubksb(matrix,ncomp_total,indx,scprod);
   }

  printf("Convolving\n");
  kernel_convolve(sub_ref,scprod);

    make_kernel(width/2,height/2,scprod);

 
  
  
 for(i=mesh_size/2;i<width-mesh_size/2;i++) 
  for(j=mesh_size/2;j<height-mesh_size/2;j++) 
   {
    background=get_background(i, j, scprod);
    conv_image[i+width*j] = (conv_image[i+width*j]+background);
    if(sub_ref[i] < PIX_MIN || sub_ref[i+width*j] < PIX_MIN )
    conv_image[i] = 0.0;
   }

 
    write_image(conv_image,file0,cx,cy); 

 

  for(i=mesh_size/2;i<width-mesh_size/2;i++) 
   for(j=mesh_size/2;j<height-mesh_size/2;j++) 
   {
    background=get_background(i,j,scprod);
    conv_image[i+width*j] = (conv_image[i+width*j]-sub_image[i+width*j])/sum_kernel;
     if(sub_image[i+width*j] < PIX_MIN || sub_ref[i+width*j] < PIX_MIN)
      conv_image[i+width*j]  = 0.0;
   }

 printf("Writing Convolution\n");

 make_kernel(0,0,scprod);
 printf("sum_kernel 0 0: %lf\n", sum_kernel);
 make_kernel(width,height,scprod);
 printf("sum_kernel width height: %lf\n", sum_kernel);

 write_image(conv_image,file,cx,cy); 

 /**************************************************************************/
 /**************************************************************************/

  pos=get_index(argv[2],"/");
  sprintf(kernel_name,"kc_%i%s",frame_number,&argv[2][pos+1]);
printf("kernel: %s\n",kernel_name);
  fprintf(kernel_table,"%s %5i %5i %5i %5i %5i %5i\n",kernel_name,mesh_size/2+
  cx ,mesh_size/2+cy,-mesh_size/2+cx+sub_width,-mesh_size/2+cy+sub_height,cx,
  cy);
printf("kernel\n");

 

  ++frame_number;
  coeff_file=fopen(kernel_name,"wb");
  fwrite(scprod,(ncomp_total0+1)*sizeof(double),1,coeff_file);
  fclose(coeff_file);
printf("kernel\n");

 /**************************************************************************/
 /**************************************************************************/



 }

 fclose(file);
 fclose(file0);
 fclose(file1);
 fclose(file2);
 fclose(kernel_table);
 free(stamps);
 free(filter_x);
 free(filter_y);
 free(temp);

}
