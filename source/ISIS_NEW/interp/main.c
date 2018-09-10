#include <math.h>
#include <stdio.h>
#include "types.h"

main(int argc, char *argv[])
{
   FILE   *ficin,*ast,*ficout;
   int    i,j,degree,nbcoeff,w1,w2,test,nb,bytenb,
          jnr,knr,k,lsup,jo,nkbl,ij,ik,n1,n2,qt;
   double coeffx[64],coeffy[64],xi,yi,t,u,
          y1,y2,y3,y4,pixmin,pixmax;
   DATA_TYPE  *pix,*pix_int,*y2a,*x2a,*x1a,*uv,*ytmp,*yytmp,y,
          xp,yp,*xxp,*yy,*yyp,*flt,*pix_bl,*err;
   char   astro[256],infic[256],header[2880*15],s2[256],cp[8],
          s2y[4],*cpp,cp2[8];

 pixmin=1.0;
 pixmax=1.0e022;
 
            
         for (i=1; i<argc; i++)
	   switch(tolower(argv[i][1]))
	      {
	        case 'i': sprintf(infic,"%s",argv[++i]);
			break;  
               
              }

  

          /**************** Get Offsets ****************/

             sprintf(astro,"coeff");
             ast=fopen(astro,"rb");


                i=0;
                fread(&j,4,1,ast);
                printf("nbcoeff %i\n",j);
                nbcoeff=j;
                degree=(int)floor(-1.5+sqrt(0.25+2.0*(double)j));
printf("degre: %i\n",degree);
                for(i=0;i<j;i++)
                {
                 fread(&coeffx[i+1],8,1,ast);
                 printf("coeffx[%i] %lf\n",i,coeffx[i+1]);
                }

                for(i=0;i<j;i++)
                {
                 fread(&coeffy[i+1],8,1,ast);
                 printf("coeffy[%i] %lf\n",i,coeffy[i+1]);
                }

 fclose(ast);

        /**********************************************/

 
 pix=readfits(infic);
 


 w1=width;
 w2=height;
 printf("w1: %i w2: %i\n", w1, w2);
 lsup=w1; if(w2>w1) lsup=w2;

 bytenb=w1*w2*sizeof(DATA_TYPE);


 if(!(y2a=(DATA_TYPE *)malloc(bytenb))){printf("Cannot allocate\n"); exit(0);}
 x1a=(DATA_TYPE *)malloc(lsup*sizeof(DATA_TYPE));
 x2a=(DATA_TYPE *)malloc(lsup*sizeof(DATA_TYPE));
 uv=(DATA_TYPE *)malloc(lsup*sizeof(DATA_TYPE));
 ytmp=(DATA_TYPE *)malloc(lsup*sizeof(DATA_TYPE));
 if(!(yytmp=(DATA_TYPE *)malloc(lsup*sizeof(DATA_TYPE)))){printf("Cannot allocate\n"); exit(0);}
 xxp=(DATA_TYPE *)malloc(lsup*sizeof(DATA_TYPE));
 yyp=(DATA_TYPE *)malloc(lsup*sizeof(DATA_TYPE));
 yy=(DATA_TYPE *)malloc(lsup*sizeof(DATA_TYPE));
 pix_int=(DATA_TYPE *)malloc(bytenb);
 pix_bl=(DATA_TYPE *)malloc(bytenb);
 err=(DATA_TYPE *)malloc(bytenb);

 for(i=0;i<w1*w2;i++) {y2a[i]=0.01; pix_int[i]=0.01;}

 

 for(i=0;i<lsup;i++) {x1a[i]=(DATA_TYPE)i; x2a[i]=(DATA_TYPE)i;}

 printf("Spline Tabulation \n");
        splie2(x1a,x2a,pix,w1,w2,y2a);
 printf("Spline Tabulation Done: %i %i %i\n",w1,w2,lsup);



 for(i=0;i<w1;i++)
 {

   for(j=0;j<lsup;j++)
   {
      xi=(double)i; yi=(double)j;
     if(degree>0)
     {
      xp=poly(xi,yi,degree,coeffx);
      yp=poly(xi,yi,degree,coeffy);
     }
     else
     {
       xp=xi+coeffx[1];
       yp=yi+coeffy[1];
     }
      xxp[j]=xp;
      yyp[j]=yp;
      jnr=(int)floor((double)xp); 
      knr=(int)floor((double)yp);     




/*********************************************************************************/

     if(jnr>-1 && (jnr+1)<w1 && knr>-1 && (knr+1)<w2)
     {
      
      t=xp-floor(xp);
      u=yp-floor(yp);

      y1=pix[jnr+knr*w1];
      y2=pix[jnr+1+knr*w1];
      y3=pix[jnr+1+(knr+1)*w1];
      y4=pix[jnr+(knr+1)*w1];
      
      pix_bl[j]=(1.0-t)*(1.0-u)*y1+t*(1.0-u)*y2+t*u*y3+(1.0-t)*u*y4;
      err[j]=pix[jnr+knr*w1];
     }
      else { pix_bl[j]=0.01; err[j]=1.0;}

/*********************************************************************************/
   }
 nkbl=0;


        splin2(x1a,x2a,pix,y2a,w1,w2,xxp,yyp,yy,uv,ytmp,yytmp); 


        for (k=0;k<w2;k++)
        {
         if(pix_bl[k]>0.02)
	 {
	  if(fabs(yy[k]- pix_bl[k]) < 2.0e12*sqrt(fabs(err[k])+1.0)) 
	  pix_int[i+w1*k]=yy[k];
	  else pix_int[i+w1*k]=pix_bl[k];
	 }
	}
	 
 }



 for (i=1;i<w1-1;i++) 
  for (k=1;k<w2-1;k++)
  {
      xi=(double)i; yi=(double)k;
      xp=poly(xi,yi,degree,coeffx);
      yp=poly(xi,yi,degree,coeffy);
      jnr=(int)xp; knr=(int)yp;
   if(jnr>-1 && jnr<w1 && knr>-1 &&  knr<w2)    
   if(pix[jnr+w1*knr]<=pixmin || pix[jnr+w1*knr]>pixmax)
    {
       for(ij=-1;ij<=1;ij++) 
        for(ik=-1;ik<=1;ik++)
        { 
         n1=i+ij; n2=k+ik; 
         if(n1<0) n1=0; if(n1>w1-1) n1=w1-1; 
         if(n2<0) n2=0; if(n2>w2-1) n2=w2-1;
         pix_int[n1+w1*n2]=1.0;
	}
    }
   }




  for(i=0;i<w1*w2;i++) {if(pix_int[i]<0.01) pix_int[i]=0.01;
  /*printf("pix: %lf\n", pix_int[i]);*/}



  sprintf(infic,"interp.fits");
  writefits(pix_int, infic);

 free(pix_int);
 free(pix);


}



