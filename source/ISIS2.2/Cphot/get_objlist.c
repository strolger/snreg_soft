#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>

#include "types.h"

object *get_objlist(char *name,int *nn)
{
 object    *obj;
 FILE      *file;
 DATA_TYPE x,y,aperphot,phot,bg;
 int       xi,yi,n;
 char      s[1000];

 file=fopen(name,"r");
 obj=(object *)malloc((NSTARS_MAX+2)*sizeof(object));

 n=0;
 while(!feof(file))
 {
  fgets(s,1000,file);
  sscanf(s,"%f %f %i %i %f %f %f",&x,&y,&xi,&yi,&aperphot,&phot,&bg);
  obj[n].x=x;
  obj[n].y=y;
  obj[n].xi=xi;
  obj[n].yi=yi;
  obj[n].aperphot=aperphot;
  obj[n].phot=phot;
  obj[n++].bg=bg;
 }
 fclose(file);

 printf("n: %i NSTARS_MAX: %i\n",n,NSTARS_MAX);

 *nn=n;

 return obj;
}
