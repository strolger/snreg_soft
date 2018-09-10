#include<stdio.h>
#include "types.h"


typedef DATA_TYPE LIST_TYPE;

void quick_sort (LIST_TYPE *list, int *index, int n)
{
   int i;
   void quick_sort_1();

   for(i=0;i<n;i++) index [i]=i;
   quick_sort_1 (list, index, 0, n-1);

   return;
}



void quick_sort_1(LIST_TYPE *list, int *index, int left_end, int right_end)
{
 int i,j,k,temp;
 LIST_TYPE chosen;
 
 
 chosen = list[index[(left_end + right_end)/2]];
 i = left_end-1;
 j = right_end+1;

 for(;;)
 {

   while(list [index[++i]] < chosen);
   while(list [index[--j]] > chosen);
   


   if (i < j){
   temp=index [j];
    index [j] = index [i];
    index [i] = temp;}
         else if (i == j) {
	   ++i; break;}
         else break;

   
 } 

 if (left_end < j)  quick_sort_1 (list, index, left_end, j);
 if (i < right_end) quick_sort_1 (list, index, i, right_end);

    return;


}
