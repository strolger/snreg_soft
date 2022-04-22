#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

void spline_allocate()
{
 int i;



 if(!(y2a=(DATA_TYPE  *)malloc(stack_width*stack_width*sizeof(DATA_TYPE))))
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


 return;
}
