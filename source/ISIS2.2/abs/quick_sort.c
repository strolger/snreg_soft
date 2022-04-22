#include<stdio.h>
#include "types.h"


typedef DATA_TYPE LIST_TYPE;

void quick_sort (LIST_TYPE *list, int *nindex, int n)
{
   int i;
   void quick_sort_1();

   for(i=0;i<n;i++) nindex [i]=i;
   quick_sort_1 (list, nindex, 0, n-1);

   return;
}



void quick_sort_1(LIST_TYPE *list, int *nindex, int left_end, int right_end)
{
 int i,j,k,temp;
 LIST_TYPE chosen;
 
 
 chosen = list[nindex[(left_end + right_end)/2]];
 i = left_end-1;
 j = right_end+1;

 for(;;)
 {

   while(list [nindex[++i]] < chosen);
   while(list [nindex[--j]] > chosen);
   


   if (i < j){
   temp=nindex [j];
    nindex [j] = nindex [i];
    nindex [i] = temp;}
         else if (i == j) {
	   ++i; break;}
         else break;

   
 } 

 if (left_end < j)  quick_sort_1 (list, nindex, left_end, j);
 if (i < right_end) quick_sort_1 (list, nindex, i, right_end);

    return;


}
