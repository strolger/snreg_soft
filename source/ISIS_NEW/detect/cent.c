#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"


 void cent(int i,int j,CENT_TYPE *x,CENT_TYPE *y,DATA_TYPE *v)
 {
    CENT_TYPE marg[3],dx,dy,num;
    int       m,k,l;


    m=1;

       for(k=-m;k<=m;k++)
       {
        marg[k+1]=0.0;
        for(l=-m;l<=m;l++)
	{
          marg[k+1] += v[i+k+width*(j+l)];

        }
       }

       num=(marg[2]+marg[0]-2.0*marg[1]);
       if(num!=0.0)
       dx=0.5*(marg[0]-marg[2])/
                 (marg[2]+marg[0]-2.0*marg[1]);
       else
       dx=0.0;

       for(l=-m;l<=m;l++)
       {
        marg[l+1]=0.0;
        for(k=-m;k<=m;k++)
	{
          marg[l+1] += v[i+k+width*(j+l)];
        }
       }


       num=(marg[2]+marg[0]-2.0*marg[1]);
       if(num!=0.0)
       dy=0.5*(marg[0]-marg[2])/
                 (marg[2]+marg[0]-2.0*marg[1]);

       else
       dy=0.0;



       if(fabs(dx)>1.0) dx=0.0; if(fabs(dy)>1.0) dy=0.0;


        *x=dx+(CENT_TYPE)i;
        *y=dy+(CENT_TYPE)j;


       

       return;
 }
