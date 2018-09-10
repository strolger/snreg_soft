#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>


int get_index(char *s1,char *s2)
{

 int  i,n,n2,pos;


 n=strlen(s1);
 n2=strlen(s2);
 pos=-1;

 for(i=n-n2;i>-1;i--) 
 {
  if(!(strncmp(&s1[i],s2,n2))) 
  {
   pos=i;
   break;
  }
	   
 }


 return pos;
}
