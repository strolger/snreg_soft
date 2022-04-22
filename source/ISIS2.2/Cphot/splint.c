#include<stdio.h>
#include<math.h>

#include "types.h"


void splint(DATA_TYPE *xa,DATA_TYPE *ya,DATA_TYPE *y2a,int n,DATA_TYPE x,DATA_TYPE *y)
{
	int klo,khi,k;
	DATA_TYPE h,b,a;

	klo=0;
	khi=n-1;


	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;}
	


	h=xa[khi]-xa[klo];
	if (h == 0.0) printf("Bad XA input to routine SPLINT\n");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;


	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
  return;
}
