#include<stdio.h>
#include<math.h>

#include "types.h"


void splie2(DATA_TYPE *x1a,DATA_TYPE *x2a,DATA_TYPE *ya,int m,int n,DATA_TYPE *y2a)
{
	int j;
        DATA_TYPE *u,y,qm;
        qm=1.0;

        if(!(u=(DATA_TYPE *)malloc(m*sizeof(DATA_TYPE))))
           {printf("Cannot allocate\n"); exit(0);}

	for (j=0;j<n;j++)
        {
		spline(x2a,&ya[j*m],m,qm,qm,&y2a[j*m],u);
        }

       free(u);
return;
}
