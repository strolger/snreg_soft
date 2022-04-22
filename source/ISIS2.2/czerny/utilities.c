#include <math.h>
#include "complex.h"

 void add(complex_mat *a, complex_mat *b, complex_mat *result)
 {
   int i;
  

   for(i=0;i<a->nb;i++){
   result->mat[i].real=a->mat[i].real+b->mat[i].real;
   result->mat[i].imaginary=a->mat[i].imaginary+b->mat[i].imaginary;}
   return;
 }


void sub(complex_mat *a, complex_mat *b, complex_mat *result)
 {
   int i;
  

   for(i=0;i<a->nb;i++){
   result->mat[i].real=a->mat[i].real - b->mat[i].real;
   result->mat[i].imaginary=a->mat[i].imaginary - b->mat[i].imaginary;}
   return;
 }

void mult(complex_mat *a, complex_mat *b, complex_mat *result)
 {
   int i;

   for(i=0;i<a->nb;i++){
   result->mat[i].real=a->mat[i].real*b->mat[i].real - a->mat[i].imaginary*b->mat[i].imaginary;
   result->mat[i].imaginary=a->mat[i].imaginary*b->mat[i].real + b->mat[i].imaginary*a->mat[i].real;}
   return;
 }

void mult_conj(complex_mat *a, complex_mat *b, complex_mat *result)
 {
   int i;

   for(i=0;i<a->nb;i++){
   result->mat[i].real=a->mat[i].real*b->mat[i].real + a->mat[i].imaginary*b->mat[i].imaginary;
   result->mat[i].imaginary=a->mat[i].imaginary*b->mat[i].real - b->mat[i].imaginary*a->mat[i].real;}
   return;
 }


void cmult(complex_mat *a, complex *b, complex_mat *result)
 {
   int i;

   for(i=0;i<a->nb;i++){
   result->mat[i].real=a->mat[i].real*b->real - a->mat[i].imaginary*b->imaginary;
   result->mat[i].imaginary=a->mat[i].imaginary*b->real + b->imaginary*a->mat[i].real;}
   return;
 }

complex scal_prod(complex_mat *a, complex_mat *b, double *gk)
{
 complex result;
 int i;

  result.real=result.imaginary=0.0;

 for(i=0;i<a->nb;i++){
  result.real      += (a->mat[i].real*b->mat[i].real + a->mat[i].imaginary*b->mat[i].imaginary)*gk[i];
  result.imaginary += (a->mat[i].imaginary*b->mat[i].real - b->mat[i].imaginary*a->mat[i].real)*gk[i];}

 return result;

}

complex comp_div(complex a, complex b)
{
 complex result;
 double den;

  den=b.real*b.real + b.imaginary*b.imaginary;
  
  result.real      = (a.real*b.real + a.imaginary*b.imaginary)/den;
  result.imaginary = (a.imaginary*b.real - b.imaginary*a.real)/den;

 return result;
}

 void copy(complex_mat *a, complex_mat *result)
 {
   int i;
  

   for(i=0;i<a->nb;i++){
   result->mat[i].real=a->mat[i].real;
   result->mat[i].imaginary=a->mat[i].imaginary;}
   return;
 }


void allocate(complex_mat *a)
{
 a->mat=(complex *)malloc(a->nb*sizeof(complex));
 return;
}

void init(complex_mat *a,double x, double y)
{
   int i;
     
    
   for(i=0;i<a->nb;i++){
    a->mat[i].real=x;
    a->mat[i].imaginary=y;
   }
   return;
}
