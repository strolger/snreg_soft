#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "complex.h"
#include "funcs.h"
#include "allocate.h"

extern double F();

double fit(double freq1,double freq2,double *ti,double *Xi,double *gk,int N,int nbharm,double *anova,double accuracy)
{
 int    i,nbren;
 double omega,moy,magmin,magmax,jd_min,jd_max,Var_xx,
        freq,pas,pi,Czerny,Czerny_max,period,*X,*t,freq_out; 
 FILE   *outfile;
 char   fileout[256];


 sprintf(fileout,"spectra.data");

  if(!(outfile=fopen(fileout,"w")))
  {
   printf("\n ! ERROR IN OPENING : %s\n\n", fileout);
   exit(0);
  }


         pi=asin(1.0)*2.0;
         X=(double *)malloc(N*sizeof(double));
         t=(double *)malloc(N*sizeof(double));
         for(i=0;i<N;i++){ X[i]=Xi[i]; t[i]=ti[i];}

        magmin=jd_min=1.0E22;
        magmax=jd_max=-1.0E22;
        moy=0.0;

        for(i=0;i<N;i++)
        {     
          if(magmin>X[i]) magmin=X[i];
          if(magmax<X[i]) magmax=X[i];
          if(jd_min>t[i]) jd_min=t[i];
          if(jd_max<t[i]) jd_max=t[i];
          moy += X[i];
        }
        moy /= (double)N;



        for(i=0;i<N;i++) { t[i] -= jd_min; X[i] -= moy;}
        for(i=0,Var_xx=0.0;i<N;i++) Var_xx += X[i]*X[i]*gk[i];
       
        pas   = accuracy/(jd_max-jd_min);
        nbren = (freq1-freq2)/pas;

 omega=freq2*2.0*pi;
 pas *= 2.0*pi;

 Phi.nb=ph.nb=phm1.nb=z.nb=temp2.nb=temp1.nb=temp3.nb=q.nb=N;
 allocate(&Phi); allocate(&ph); allocate(&phm1); allocate(&z);
 allocate(&temp1); allocate(&temp2); allocate(&temp3), allocate(&q); 


 Czerny_max=0.0;
 for(i=0;i<nbren;i++)
 {
     Czerny=F(omega,X,t,gk,N,nbharm,Var_xx);
     if(Czerny_max<Czerny) {Czerny_max=Czerny; freq_out=omega/(2.0*pi);}
     fprintf(outfile,"%f %f \n", omega/(2.0*pi), Czerny);
     omega += pas;
 }
 fclose(outfile);

 free(Phi.mat); free(ph.mat); free(phm1.mat); free(z.mat); free(temp1.mat); 
 free(temp2.mat); free(temp3.mat); free(q.mat); free(X); free(t);
 *anova=Czerny_max;

 return freq_out;
}


