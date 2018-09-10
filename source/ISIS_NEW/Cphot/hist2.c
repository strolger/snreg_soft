#include<stdio.h>
#include<math.h>



typedef struct {
                float x[10000];
		float y[10000];
		float   pixnb[10000];
		} st;


main(argc,argv)
int argc;
char *argv[];
{
	FILE    *ex1,*ex2,*ex3;

	int     i,j,k,l,ij,nhistx,nhisty,i1,i2,nb,a,imax,no,
		qflag,lflag,nbf,width,NBIN,BIN,NMAX;

	char    exin1[256],exin2[256],exout[256],s[1000];
	float   xx,yy,xin,yin,x0,m0,m1,m2,m3,thrs,chi,cor,pixnb;
	double  offx,offy,max,diffy,diffx,*correl;
	long    place;
	st      st1,st2;


	place = (long)0;
	qflag = 0;
	lflag = 0;
        BIN   = 2;
        width = 500;
        NMAX  = 500;

	for (a=1; a<argc; a++)
           switch(tolower(argv[a][1]))
              {
                case 'i': sprintf(exin1,argv[++a]);
                        break;
                case 'o': sprintf(exin2,argv[++a]);
                        break;
		case 'q': qflag=1;
			 break;
		case 'b': BIN=atoi(argv[++a]);
			 break;
		default : {
			   printf("SYNTAX: rbin -i File -l first line to read \n");
			   exit(0);
			  }
	      }

	sprintf(exout,"coeff");

	ex1=fopen(exin1,"r");
	ex2=fopen(exin2,"r");
	ex3=fopen(exout,"w");
	
	
	

	
/************** Read Files*****************/
	
	i1=0;
	while(!feof(ex1) && i1<NMAX)
	{
		fgets(s,1000,ex1);
                sscanf(s,"%f %f %f",&xx,&yy,&pixnb);
              if(xx<width-1 && yy<width-1)
	      {
	      	st1.x[i1]    = xx;
		st1.y[i1]    = yy;	
		st1.pixnb[i1] = pixnb;
		++i1;
              }
	}


	i2=0;
	while(!feof(ex2) && i2<NMAX)
	{

		fgets(s,1000,ex2);
                sscanf(s,"%f %f %f",&xx,&yy,&pixnb);
               if(xx<width-1 && yy<width-1)
	       {
		st2.x[i2]    = xx;
		st2.y[i2]    = yy;	
		st2.pixnb[i2] = pixnb;
		++i2;
               }
	}



/*************** Build X Histogram ***********************/

     max = 0;
     NBIN=2*width/BIN+1;   
     correl=(double *)malloc(NBIN*NBIN*sizeof(double));
   

	for(i=0;i<NBIN*NBIN;i++) correl[i] = 0.0; 


       
             for(i=0;i<i1;i++)
               for(j=0;j<i2;j++){
                 diffx=st1.x[i]-st2.x[j]+(double)width;
                 diffy=st1.y[i]-st2.y[j]+(double)width;
                   nhistx=(int)floor(diffx/(double)BIN);
                   nhisty=(int)floor(diffy/(double)BIN); 

                   correl[nhistx+nhisty*NBIN] += 1.;
         }       
          


	
	for(i=0;i<NBIN*NBIN;i++)
	{

	  if( correl[i] > max ) 
	  { 
		max  = correl[i];
	        imax = i;

	  }
	  
	}


        k=imax/NBIN;
        l=imax-k*NBIN;
	offx =l*BIN-width+BIN/2;
        offy=k*BIN-width+BIN/2;
	printf("offx: %f offy: %f\n",offx,offy);

        i=1;
	
	fwrite(&i,4,1,ex3);
	fwrite(&offx,8,1,ex3);
	fwrite(&offy,8,1,ex3);
	
	

	
	printf("number of sources read: %i %i\n", i1,i2);
	fclose(ex1);
	fclose(ex2);
	fclose(ex3);
}
