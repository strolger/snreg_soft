#include<math.h>
#include<stdio.h>


double   poly(double x, double y,int ndeg,double *coeff)
{
	double  pol,f[1000];
	int     j,k,njk,nterm;
	

	pol   = coeff[1];
	nterm = (ndeg+1)*(ndeg+2)/2; 


		f[0] = 1.0;
	        njk  = 1;

	     for(j=1;j<=ndeg;j++)
	     {
		for(k=0;k<j;k++)
	        {
		   f[njk] = f[njk-j]*y;
		   pol   += coeff[njk+1]*f[njk];
		   ++njk;
	        }

		 f[njk] = f[njk-j-1]*x;
		  pol   += coeff[njk+1]*f[njk];
		 ++njk;
	     }



	return pol;


}