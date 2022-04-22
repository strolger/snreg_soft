#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "complex.h"
#include "funcs.h"
#include "allocate.h"


extern double fit();
 
int main(int argc,char *argv[])
{
 FILE        *infile,*outfile;
 int         nbharm,i,j,harmonics,N,qflag,a,post,nbren,
             nx,ny;
 double      *gk,*t,*X,omega,pi,period1,period2,freq1,
             freq2,jd_min,jd_max,pas,magmin,magmax,mmag,
             jjulday,error,moy,Czerny,Czerny_max,period,
             freq,m1,m2,anova,accuracy;
 char        filein[256],fileout[256],s[1024],cbuf[1024];
 double       *dismag,*disphase;

         qflag=post=0;
         period1=100.0;
         period2=1000.0;
         harmonics=2;
         sprintf(fileout, "lc.data");
	 for (a=1; a<argc; a++)
           switch(tolower(argv[a][1]))
              {
                case 'a': period1=atof(argv[++a]);
                        break;
                case 'b': period2=atof(argv[++a]);
                        break;
                case 'i': sprintf(filein, "%s", argv[++a]);
                        break;
                case 'n': harmonics=atoi(argv[++a]);
                        break;
 		case 'o': sprintf(fileout, "%s", argv[++a]);
                        break;
                case 'p': post=1;
                        break;       
                case 'q': qflag=1;
                        break;
                             
              }

        if(!(infile=fopen(filein,"r")))
        {
         printf("\n ! ERROR IN OPENING : %s\n\n", filein);
         exit(0);
        }
        if(!(outfile=fopen(fileout,"w")))
        {
         printf("\n ! ERROR IN OPENING : %s\n\n", fileout);
         exit(0);
        }

        fgets(s,256,infile);
        i=0;
        while(feof(infile) == 0)
        {
           fgets(s,256,infile); ++i;
        }
        N=i-1;
        rewind(infile);

         pi=asin(1.0)*2.0;
         t=(double *)malloc(N*sizeof(double));
         gk=(double *)malloc(N*sizeof(double));
         X=(double *)malloc(N*sizeof(double));
         dismag=(double *)malloc(N*2*sizeof(double));
         disphase=(double *)malloc(N*2*sizeof(double));


        magmin=jd_min=1.0E22;
        magmax=jd_max=-1.0E22;
        fgets(s,256,infile);
        moy=0.0;
        for(i=0;i<N;i++)
        {
          fgets(s,256,infile);
          sscanf(s," %lf %lf %lf",&jjulday,&mmag ,&error);
          X[i]    = mmag;
          /*jjulday /= 3600.0*24;*/
          t[i] = jjulday;
          if(magmin>mmag) magmin=mmag;
          if(magmax<mmag) magmax=mmag;
          if(jjulday<jd_min) jd_min=jjulday;
          if(jjulday>jd_max) jd_max=jjulday;
          gk[i]=1.0/error/error;
          moy += mmag;
         
          if(!qflag) printf("%lf %lf %lf\n",t[i],X[i],gk[i]);
        }       
         moy /= (double)N;




 freq1=1.0/period1;
 freq2=1.0/period2;
 
 nbharm=harmonics*2+1;
 accuracy=0.01;
 freq=fit(freq1,freq2,t,X,gk,N,nbharm,&anova,accuracy);
 pas   = accuracy/(jd_max-jd_min);
 period=1.0/freq;

 printf("freq: %lf\n", freq);

 

 freq1=freq+5.0*pas;
 freq2=freq-5.0*pas;
 accuracy=0.001;
  pas   = accuracy/(jd_max-jd_min);
 freq=fit(freq1,freq2,t,X,gk,N,nbharm,&anova,accuracy);

 printf("freq: %lf\n", freq);

printf("pas: %lf %lf %lf\n",accuracy,pas*(jd_max-jd_min),jd_max-jd_min);
 period=1.0/freq;


 printf("Period: %15.9lf %8.5lf %s\n", period,anova,filein);

 
        for(i=0;i<N;i++)
        {
          disphase[i]=t[i]*freq-floor(t[i]*freq);
          dismag[i]=X[i];
          dismag[i+N]=dismag[i];
          disphase[i+N]=disphase[i]+1.0;
          fprintf(outfile,"%f %f %lf\n",disphase[i],dismag[i],t[i]);
        }


  m1=magmin-(magmax-magmin)/2.25;
  m2=magmax+(magmax-magmin)/3.5;

  /*if(post) metafl("post");
  else metafl("xwin");
  setpag("da4l");
  disini();
  pagera();
  complx();
  color("blue");
  graf(0.,2.,0.,0.25,m1,m2,m1,1.0);
  title();
  incmrk(-1);
  marker(21);
  hsymbl(15);
  color("blue");
  curve(disphase,dismag,2*N);
  hsymbl(1);
  legini(cbuf,1,8);
  nx = nxposn(0.1);
  ny = nyposn(m1+(magmax-magmin)/4.0);
  legpos(nx,ny);
  sprintf(s,"%5.1lf",1.0/freq);
  leglin(cbuf,s,1);
  legtit(filein);
  legend(cbuf,2);
  disfin();*/
 

  
  fclose(infile);
  fclose(outfile);


 free(X);
 free(t);


 
}


