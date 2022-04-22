#include <math.h>
#include "complex.h"
#include "funcs.h"
#include "allocate.h"
 
double F(double omega,double *X,double *t,double *gk, int N, int nbharm, double Var_xx)
{
 complex     alpha,c,c1,c2;
 double      den,sum;
 int         i,j,nb;



 /****************************************************************/
 /***** Phi=z^N*X ************************************************/
 /****************************************************************/

  nb=(int)floor((double)nbharm/2.0-0.4);
  for(i=0;i<N;i++) {Phi.mat[i].real=cos(omega*t[i]*(double)nb)*X[i];
   Phi.mat[i].imaginary=sin(omega*t[i]*(double)nb)*X[i];}

 /****************************************************************/
 
 for(i=0;i<N;i++){ z.mat[i].real=cos(omega*t[i]);
 z.mat[i].imaginary=sin(omega*t[i]);}
 sum=0.0;

 for(i=0;i<nbharm;i++){
 if(i<1) init(&ph,1.0,0.0); 
  else{
   mult_conj(&q,&phm1,&temp1);     /***** z^(n-1)*ph_(n-1) ******/
   cmult(&temp1,&alpha,&temp3);    /***** alpha*z^(n-1)*ph_(n-1) ****/
   mult(&z,&phm1,&temp1);          /***** z*ph_(n-1) ******/

   sub(&temp1,&temp3,&ph);         /**** z*ph_(n-1)-alpha*z^(n-1)*ph_(n-1) ***/
  }

 

  if(i<1)  init(&q,1.0,0.0);
   else{
   mult(&q,&z,&temp1);             /*** z*q **********/
   copy(&temp1,&q);}               /**** q=q*z *******/

  mult(&z,&ph,&temp1);      
  mult_conj(&q,&ph,&temp2);
  c1=scal_prod(&temp1,&ph,gk);
  c2=scal_prod(&temp2,&ph,gk);
  alpha=comp_div(c1,c2);
  den=0.0;
  for(j=0;j<N;j++) den += (ph.mat[j].real*ph.mat[j].real+
  ph.mat[j].imaginary*ph.mat[j].imaginary)*gk[j];
  den=sqrt(den);
  c=scal_prod(&Phi,&ph,gk);
  c.real /= den;
  c.imaginary /= den;
  sum += c.real*c.real + c.imaginary*c.imaginary;

  copy(&ph,&phm1);
 }

 /*sum=(double)(N-nbharm-2)*sum/((double)(nbharm-1)*(Var_xx-sum));*/
 
 sum /= Var_xx;

 return sum;
}

