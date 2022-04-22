#include<math.h>
#include<stdio.h>
#include "fit.h"



main(argc,argv)
int  argc;
char *argv[];
{	
	FILE     *fic,*ficout;
	int     i,j,n,nb,ndeg,a,nterm,qflag,ik,ik0,deltik,flag;
	double  x[30000],y[30000],zx[30000],zy[30000],
		nx[30000],ny[30000],nzx[30000],nzy[30000],
		coeffx[1000],coeffy[1000],deltax,deltay,delta,clip,
		sx,sy;
	char    s[1000],filein[256],fileout[256];
	
	
	sprintf(fileout,"coeff");
	ndeg  = 2; 
	
	    qflag=0;
	   for (a=1; a<argc; a++)
           switch(tolower(argv[a][1]))
              {
		case 'd': ndeg = atoi(argv[++a]);
                          break;

		case 'q': qflag=1;			 
                          break;

                case 'i': sprintf(filein,"%s",argv[++a]);			 
                          break;

		case 'o': sprintf(fileout,"%s",argv[++a]);
                          break;
		
	      }

		if(!(fic = fopen(filein,"r")))
		nnerror(filein);
		if(!(ficout = fopen(fileout,"w")))
		nnerror(fileout);	

	n  = 0;
	
	nterm = (ndeg+1)*(ndeg+2)/2;

	while(feof(fic) == 0)
	{
	 fgets(s,1000,fic);
	 sscanf(s,"%i  %lf %lf %lf %lf",&nb,&zx[n],&zy[n],&x[n],&y[n]);
	 ++n;	  
	}
         

	fit(x,y,zx,ndeg,n-1,coeffx,1);
	fit(x,y,zy,ndeg,n-1,coeffy,0);

	 sx = sigma(x,y,zx,ndeg,n-1,coeffx) ;
	 sy = sigma(x,y,zy,ndeg,n-1,coeffy) ;
	 clip = 9.0*( sx*sx + sy *sy)/2.0;

                 if(clip<1.0e-10) clip=1.0e-10;


	ik=n-1;
	deltik=1;
	while(deltik>0)
	{
	  ik0=ik;
	  ik=0;
	  for(i=0;i<n-1;i++)
	  {
	   deltax = poly(x[i],y[i],ndeg,coeffx)-zx[i];
	   deltay = poly(x[i],y[i],ndeg,coeffy)-zy[i];
	   delta = deltax*deltax + deltay*deltay;

	   nx[ik] = x[i];
	   ny[ik] = y[i];
	   nzx[ik] = zx[i];
	   nzy[ik] = zy[i];
	   
	   if(delta < clip) ++ik;
	  } 
		
		deltik = ik0-ik;

          if(ik<1) break;
	  fit(nx,ny,nzx,ndeg,ik,coeffx,1);
	  fit(nx,ny,nzy,ndeg,ik,coeffy,0);
          if(ndeg>0) 
	  {
	   sx = sigma(nx,ny,nzx,ndeg,ik-1,coeffx);
           sy = sigma(nx,ny,nzy,ndeg,ik-1,coeffy);
          }
          else
	  {
            sx=sy=0.0;
            for(i=0;i<ik;i++) sx += (nzx[i]-x[i]-coeffx[1])*(nzx[i]-x[i]-coeffx[1]);
            for(i=0;i<ik;i++) sy += (nzy[i]-y[i]-coeffy[1])*(nzy[i]-y[i]-coeffy[1]);
           sx=sqrt(sx/(double)ik);
           sy=sqrt(sy/(double)ik);

          }
	  clip = 9.0*( sx*sx + sy *sy)/2.0;
	  printf("sigmax: %lf sigmay: %lf ndata: %i nrest: %i\n",sx,sy,n-1,ik-1);
	}

	fwrite(&nterm,4,1,ficout);
	for(i=1;i<=nterm;i++)
	{
	   fwrite(&coeffx[i],8,1,ficout);
	}
	for(i=1;i<=nterm;i++)
	{
	   fwrite(&coeffy[i],8,1,ficout);
	}

	if(!qflag)
	for(i=0;i<n-1;i++)
	{
	   printf("zx[i] %lf poly %lf\n",zx[i],poly(x[i],y[i],ndeg,coeffx));
	   printf("zy[i] %lf poly %lf\n",zy[i],poly(x[i],y[i],ndeg,coeffy));
	}
        

	fclose(fic);
}

void nnerror(s)
char *s;
{
	printf("Error in opening: %s\n",s);
	exit(0);
	return;
}
