#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void make_psf(object *objlist,DATA_TYPE *im,int nn,int cx,int cy,int frame_count,FILE *outfile,char *psf_file)
{
 int       i,ii,j,xi,yi,xi2,yi2,mesh,xi1,yi1,*index,
           *index2,nn2,ind,jx,jy;
 DATA_TYPE x,y,dx,dy,y1,y2,y3,y4,bg,phot,*temp,*data;
 DATA_TYPE *x1a,*x2a,*ytmp,*y2a,*xxp,*yyp,*yy,*yytmp,*pix,*uv;
 double    c0,c1,c2,correl,ph,rx,c1p,c2p,norm,c0p;
 

 


 mesh=stack_width/2;
 temp = (DATA_TYPE *)malloc(nn*sizeof(DATA_TYPE));
 index = (int *)malloc(nn*sizeof(int));
 index2 = (int *)malloc(nn*sizeof(int));

 if(!(y2a=(DATA_TYPE  *)malloc(stack_width*stack_width*sizeof(DATA_TYPE))))
 {printf("Cannot allocate\n"); exit(0);}
 if(!(pix=(DATA_TYPE  *)malloc(stack_width*stack_width*sizeof(DATA_TYPE))))
 {printf("Cannot allocate\n"); exit(0);}
 x1a=(DATA_TYPE  *)malloc(stack_width*sizeof(DATA_TYPE));
 x2a=(DATA_TYPE  *)malloc(stack_width*sizeof(DATA_TYPE));
 uv=(DATA_TYPE  *)malloc(stack_width*sizeof(DATA_TYPE));
 ytmp=(DATA_TYPE  *)malloc(stack_width*sizeof(DATA_TYPE));
 if(!(yytmp=(DATA_TYPE  *)malloc(stack_width*sizeof(DATA_TYPE))))
 {printf("Cannot allocate\n"); exit(0);}
 xxp=(DATA_TYPE  *)malloc(stack_width*sizeof(DATA_TYPE));
 yyp=(DATA_TYPE  *)malloc(stack_width*sizeof(DATA_TYPE));
 yy=(DATA_TYPE  *)malloc(stack_width*sizeof(DATA_TYPE));
 for(i=0;i<stack_width;i++) {x1a[i]=(DATA_TYPE)i; x2a[i]=(DATA_TYPE)i;}

 if(nn>nstars_buffer) 
 {
  stack = (DATA_TYPE **)realloc(stack,nn*sizeof(DATA_TYPE *));
  for(i=nstars_buffer;i<nn;i++) 
  stack[i]=(DATA_TYPE *)malloc(stack_width*stack_width*sizeof(DATA_TYPE));
  nstars_buffer=nn;
 }

 for(i=0;i<nn;i++)
 {
   x=objlist[i].x;
   y=objlist[i].y;
   xi=objlist[i].xi;
   yi=objlist[i].yi;
   bg=objlist[i].bg;
   phot=objlist[i].aperphot;
   dx=x-(DATA_TYPE)xi;
   dy=y-(DATA_TYPE)yi;
 for(xi2=xi-mesh;xi2<=xi+mesh;xi2++)
  for(yi2=yi-mesh;yi2<=yi+mesh;yi2++)
   {
     xi1=xi2-xi+mesh;
     yi1=yi2-yi+mesh;
     pix[xi1+yi1*stack_width] = im[xi2+width*yi2];
   }
   for(ii=0;ii<stack_width*stack_width;ii++) y2a[ii]=0.0;

        splie2(x1a,x2a,pix,stack_width,stack_width,y2a);

      for(xi1=0;xi1<stack_width;xi1++)
      {
         for(yi1=0;yi1<stack_width;yi1++)
         {
            xxp[yi1]=(DATA_TYPE)xi1+dx;
            yyp[yi1]=(DATA_TYPE)yi1+dy;
         }

         splin2(x1a,x2a,pix,y2a,stack_width,stack_width,xxp,yyp,yy,uv,ytmp,yytmp); 
         for(yi1=0;yi1<stack_width;yi1++) 
         {
           stack[i][xi1+yi1*stack_width]=yy[yi1]-bg;
         }
      }


 }

 nn2=NSTARS;
 if(nn<NSTARS) nn2=nn;
 c0=c0p=0.0;

  for(jx=0;jx<stack_width;jx++)
  {
   for(jy=0;jy<stack_width;jy++)
   {
    j =  jx+jy*stack_width;
    for(i=0;i<nn2;i++) temp[i]=stack[i][j]/objlist[i].aperphot;
    quick_sort(temp,index,nn2);

    psf[j]=temp[index[nn2/2]];
    c0 += psf[j]*psf[j];
    rx=(jx-stack_width/2)*(jx-stack_width/2)+(jy-stack_width/2)*(jy-stack_width/2);
    rx=sqrt(rx);
    if(rx<=radphot)
    {
         c0p += psf[j]*psf[j];
    }
   }
  }
 

 for(i=0;i<nn;i++)
 {
  data=stack[i];
  c1=c2=c1p=c2p=0.0;
  for(jx=0;jx<stack_width;jx++)
  {
   for(jy=0;jy<stack_width;jy++)
   {
    j =  jx+jy*stack_width;
    c1 += data[j]*psf[j];
    c2 += data[j]*data[j];
    rx=(jx-stack_width/2)*(jx-stack_width/2)+(jy-stack_width/2)*(jy-stack_width/2);
    rx=sqrt(rx);
    if(rx<=radphot)
    {
     c1p += data[j]*psf[j];
     c2p += data[j]*data[j];
    }
   }
  }
  correl=c1/sqrt(c2*c0);
  objlist[i].phot = c1p/c0p;
  printf("correl: %f phot: %f %f\n",correl,c1p/c0p,objlist[i].aperphot);
  temp[i]=correl;
 }



 quick_sort(temp,index2,nn);
 c0=norm=0.0;
 
 for(jx=0;jx<stack_width;jx++)
  {
   for(jy=0;jy<stack_width;jy++)
   {
    j =  jx+jy*stack_width;
    
    for(i=0;i<nn2;i++){ind=index2[nn-i-1]; temp[i]=stack[ind][j]/objlist[ind].phot;}
    quick_sort(temp,index,nn2);
    psf[j]=temp[index[nn2/2]];
    rx=(jx-stack_width/2)*(jx-stack_width/2)+(jy-stack_width/2)*(jy-stack_width/2);
    rx=sqrt(rx);
    if(rx<=radphot)
    {
     c0 += psf[j]*psf[j];
    }
    if(rx<=rad_aper) norm += psf[j];
   }
  }

 for(j=0;j<stack_width*stack_width;j++) psf[j] /= norm;
  c0 /= norm*norm;

  printf("Writing Fits\n");
  writefits(psf,stack_width,stack_width,psf_file);

 for(i=0;i<nn;i++)
 {
  ind=index2[nn-i-1];
  data=stack[ind];
  c1=c2=c1p=c2p=0.0;
  for(jx=0;jx<stack_width;jx++)
  {
   for(jy=0;jy<stack_width;jy++)
   {
    j =  jx+jy*stack_width;
   
    rx=(jx-stack_width/2)*(jx-stack_width/2)+(jy-stack_width/2)*(jy-stack_width/2);
    rx=sqrt(rx);
    if(rx<=radphot)
    {
     c1 += data[j]*psf[j];
     c2 += data[j]*data[j];
    }
   }
  }

 


  correl=c1/sqrt(c2*c0);
  printf("x: %f y: %fq correl: %lf %lf %f\n",objlist[ind].x+(DATA_TYPE)cx,objlist[ind].y+(DATA_TYPE)cy,correl,c1/c0,objlist[ind].aperphot);
  fprintf(outfile,"%9.3f %9.3f %7i %7i %12.3f %12.3f %9.3f %5i\n",objlist[ind].x+(DATA_TYPE)cx,objlist[ind].y+(DATA_TYPE)cy,
objlist[ind].xi,objlist[ind].yi,objlist[ind].aperphot,c1/c0,objlist[ind].bg,frame_count);
 
 }


 
 

 free(temp);
 free(index);
 free(index2);
 free(y2a);
 free(pix);
 free(x1a);
 free(x2a);
 free(uv);
 free(ytmp);
 free(yytmp);
 free(xxp);
 free(yyp);
 free(yy);



 return;
}
