#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "image.h"
#include "functions.h"

double power(double x, int n)
{
 int    i;
 double p;

 p=1.0; 
 if(n>0) for(i=0;i<n;i++)   p *= x;
 
 return p;
}
