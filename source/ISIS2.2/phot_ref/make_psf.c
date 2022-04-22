#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void make_psf(object *objlist,DATA_TYPE *im,int nn)
{
 int       i,ii,j,xi,yi,xi2,yi2,mesh,xi1,yi1,*index;
 DATA_TYPE x,y,dx,dy,y1,y2,y3,y4,bg,phot,*temp;
 DATA_TYPE *x1a,*x2a,*ytmp,*y2a,*xxp,*yyp,*yy,*yytmp,*pix,*uv;

 mesh=stack_width/2;
 temp = (DATA_TYPE *)malloc(nn*sizeof(DATA_TYPE));
 index = (int *)malloc(nn*sizeof(int));
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
           stack[i][xi1+yi1*stack_width]=(yy[yi1]-bg)/phot;
         }
      }

 }

 

 for(j=0;j<stack_width*stack_width;j++)
 {
  for(i=0;i<nn;i++) temp[i]=stack[i][j];
  quick_sort(temp,index,nn);
  psf[j]=temp[index[nn/2]];
 }




 free(temp);
 free(index);
 return;
}
