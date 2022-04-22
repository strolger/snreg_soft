#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"
#include "vectors.h"

char read_config(char *filename)
{
 char   check,s[1024],s2[256];
 FILE   *file;
 double  nd,qx;
 int     ni;

 check=1;

  if(!(file=fopen(filename,"r"))) 
     {printf("Oops...cannot find any configuration file...Goodbye !\n"); exit(0);}


 while(feof(file) == 0)
 {
   fgets(s,256,file);
   sscanf(s,"%s",s2);
   if(!(strncmp(s2,"nstamps_x",9)))
   {
    sscanf(s,"%s %i",s2,&nstamps_x);
   }
   if(!(strncmp(s2,"nstamps_y",9)))
   {
    sscanf(s,"%s %i",s2,&nstamps_y);
   }
   if(!(strncmp(s2,"half_mesh_size",14)))
   {
    sscanf(s,"%s %i",s2,&half_mesh_size);
   }
   if(!(strncmp(s2,"half_stamp_size",15)))
   {
    sscanf(s,"%s %i",s2,&half_stamp_size);
   }
   if(!(strncmp(s2,"ngauss",6)))
   {
    sscanf(s,"%s %i",s2,&ngauss);
   }
   if(!(strncmp(s2,"deg_gauss1",10)))
   {
    sscanf(s,"%s %i",s2,&deg_fixe[0]);
   }
   if(!(strncmp(s2,"deg_gauss2",10)))
   {
    sscanf(s,"%s %i",s2,&deg_fixe[1]);
   }
   if(!(strncmp(s2,"deg_gauss3",10)))
   {
    sscanf(s,"%s %i",s2,&deg_fixe[2]);
   }
   if(!(strncmp(s2,"deg_spatial",11)))
   {
    sscanf(s,"%s %i",s2,&deg_spatial);
   } 
   if(!(strncmp(s2,"sigma_gauss1",12)))
   {
    sscanf(s,"%s %lf",s2,&nd);
    sigma_gauss[0]=1.0/(2.0*nd*nd);
   }
   if(!(strncmp(s2,"sigma_gauss2",12)))
   {
    sscanf(s,"%s %lf",s2,&nd);
    sigma_gauss[1]=1.0/(2.0*nd*nd);
   }
   if(!(strncmp(s2,"sigma_gauss3",12)))
   {
    sscanf(s,"%s %lf",s2,&nd);
    sigma_gauss[2]=1.0/(2.0*nd*nd);
   }
   if(!(strncmp(s2,"saturation",10)))
   {
    sscanf(s,"%s %lf",s2,&qx);
    SATURATION=qx;
   }
   if(!(strncmp(s2,"pix_min",7)))
   {
    sscanf(s,"%s %lf",s2,&qx);
    PIX_MIN=qx;
   }
   if(!(strncmp(s2,"min_stamp_center",16)))
   {
    sscanf(s,"%s %i",s2,&max_stamp_thresh);
   }
   if(!(strncmp(s2,"deg_bg",6)))
   {
    sscanf(s,"%s %i",s2,&deg_bg);
   }
   if(!(strncmp(s2,"sub_x",5)))
   {
    sscanf(s,"%s %i",s2,&sub_width);
   }
   if(!(strncmp(s2,"sub_y",5)))
   {
    sscanf(s,"%s %i",s2,&sub_height);
   }
   
 }



 fclose(file);

 return check;
}
