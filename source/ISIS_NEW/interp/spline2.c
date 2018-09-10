#include "types.h"


void splin2(DATA_TYPE *x1a,DATA_TYPE *x2a,DATA_TYPE *ya,DATA_TYPE *y2a,int m,int n,DATA_TYPE *x2,DATA_TYPE *x1,DATA_TYPE *y,DATA_TYPE *u,DATA_TYPE *ytmp,DATA_TYPE *yytmp)
{
	int j;
	DATA_TYPE *vector();



	for (j=0;j<n;j++)
		splint(x2a,&ya[j*m],&y2a[j*m],m,x2[j],&yytmp[j]);


	spline(x1a,yytmp,n,1.0e30,1.0e30,ytmp,u);


	for (j=0;j<n;j++){
          
          if(x1[j]>-1.0 && x1[j]<(double)n && x2[j]>-1.0 && x2[j]<(double)m) splint(x1a,yytmp,ytmp,n,x1[j],&y[j]);
          else y[j]=10.0; /*printf("j: %i y: %lf %lf %lf %i %i\n",j,y[j],x1[j],x2[j],n,m);*/}

	return;

}
