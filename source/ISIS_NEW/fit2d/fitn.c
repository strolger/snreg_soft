#include<math.h>
#include<stdio.h>
#include "fit.h"

extern double sigma();

void fit(double *x,double *y,double *z,int ndeg,int ndata,double *v,int flag)
{

	int     i,j,k,njk,nterm,indx[1000];
	double  **a,m[1000],f[1000],d;


	nterm = (ndeg+1)*(ndeg+2)/2;

	a = (double **)malloc((nterm+1)*sizeof(double *));
        for(i=0;i<=nterm;i++)
        a[i] = (double *)malloc((nterm+1)*sizeof(double));


	


	for(j=1;j<=nterm;j++)
	{
	  v[j] = 0.0;
	 for(k=1;k<=nterm;k++)
	 {
	   a[j][k] = 0.0;	 
	 }
	}	



	

        

	for(i=0;i<ndata;i++)
	{

/************ Build Moment Table **************/


		f[0] = 1.0;
	        njk  = 1;

	     for(j=1;j<=ndeg;j++)
	     {
		for(k=0;k<j;k++)
	        {
		   f[njk] = f[njk-j]*y[i];
		   ++njk;
	        }

		 f[njk] = f[njk-j-1]*x[i];
		 ++njk;
	     }

/*********************************************/
	
		


	for(j=0;j<nterm;j++)
	{
	  v[j+1] += z[i]*f[j];

	 for(k=0;k<=j;k++)
	 {
           	a[j+1][k+1] += f[j]*f[k];	 
	 }
	}

     }

      



	for(j=1;j<=nterm;j++)
	{
	 for(k=1;k<j;k++)
	 {
           		a[k][j] = a[j][k];	 
	 }
	}

  

/*********************************************/
	
		


	printf("flag: %i\n",flag);

       if(nterm>1)
       {
	ludcmp(a,nterm,indx,&d);
	lubksb(a,nterm,indx,v);
       }
       else {       v[1]=0.0;
       if(flag) for(i=0;i<ndata;i++) v[1] += z[i]-x[i];
       else for(i=0;i<ndata;i++) v[1] += z[i]-y[i];
       v[1] /= (double)ndata;}

    	for(j=1;j<=nterm;j++)
	{
	      printf("v[%i]: %lf\n",j,v[j]);
	}


        for(i=0;i<=nterm;i++)
	free(a[i]);
	free(a);

	return;
}
