
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>




#define   NBSTARMAX   500
#define   MAXLINE     500


#define	        CLIP    2.5
#define		SYNTAX \
"box -i Picture1 -j outfile name -t image size \n\
   options are:                  -s box size\n\
 			         -n maximun number of stars in a box\n"

extern void   gcross(),plcirc();
extern double poly();


struct mag {
	double x;
	double y;       
        double mag;
	};

struct fd {
        struct mag *object;
        int        nbref;

        };



typedef struct mag    ptst;
typedef struct fd     findstruct;



int
main(int argc,char *argv[])
{
	FILE      *fic1,*fic2,*fic3,*ast,*ficout,*ys;

	ptst      *dataref,*sm,*sg,sx;



	int       i,j,n,*nbtry1,widbox,nxy,n1,nx,ny,mno,a,icross,ij,ma,mfit,lista[10],
		  n2,imsiz,offset,ifit,iclip,ivar,nbcof,qflag,**boxref,find,*writeit,
		  nadd,iclip0,nnref,inref,ntab[30000],ichisq3,numchi[100],
		  numplate,pass,degree,fitflag,nbcoeff,ncent,ic,imsiz_y;

	long      nb,nbb,nbelt;

	size_t    place;

	double    boxsiz,*sig,chixy,tchif,radius,test,sigmax,
		  chisq,chisq2,chisqx,chisqy,**covar,coeffx[50],
		  coeffy[50],X,Y,f[20],offpavx,pxx,pyy,mag,max_fit,
		  offpavy,num,magtest,offpavx0,offpavy0,test1,test2
                  ,xcent,ycent,ad,b
                  ,diffxcent,diffycent,rmagcent,magcent,mx,my,maa
                  ,mdd,min_x,max_x,min_y,max_y,pixsize;

        double    *magc,*magref,*color,*yref,*xref,*diffx,*diffy;

	float     mea0,mea1,mea2,mea3,thrs2,chisq3[100]
		  ,mcor,mchi,thrs,flag1,flag2,Xm,Ym,tmag,dcent;

	char      s[1000],ficnamin1[256],ficnamin2[256],outfic[256],
		  astro[256],ofic[256],chiflag,corflag,chiwrite;
	findstruct findit;

	short     intmag;

 printf("Memory Allocation\n");	
        magc   = malloc(MAXLINE*sizeof(double));
        magref = malloc(MAXLINE*sizeof(double));
        color  = malloc(MAXLINE*sizeof(double));
        yref   = malloc(MAXLINE*sizeof(double));
        xref   = malloc(MAXLINE*sizeof(double));
        diffx  = malloc(MAXLINE*sizeof(double));
        diffy  = malloc(MAXLINE*sizeof(double));
 printf("Memory Allocation Done\n");	      

        pixsize=3600.0;
        boxsiz = 8.0;
	nb     = (long)15;
        sigmax = 1.0;
	max_fit = 3.0;	

/*******************************************************************/

	fitflag=0;
	sprintf(outfic, "dayfile");
	sprintf(astro, "coeff");

 	if(argc < 5) printf("%s\n", SYNTAX);
	qflag=0;
	
        radius=2.0;
	for (a=1; a<argc; a++)
	   switch(tolower(argv[a][1]))
	      {
	      case 'c':	sprintf(astro, "%s", argv[++a]);
			break;
              case 'd':	degree=atoi(argv[++a]);
			break;
	      case 'f': fitflag=1;
			break;
	      case 'i':	sprintf(ficnamin1, "%s", argv[++a]);
			break;
	      case 'j': sprintf(ficnamin2, "%s", argv[++a]);
			break;
              case 'l': max_fit=atof(argv[++a]);
                        break;               
	      case 'o': sprintf(outfic, "%s", argv[++a]);
			break;
	      case 'm':	sigmax=atof(argv[++a]);
			break;
	      case 'n':	nb=(long)atoi(argv[++a]);
			break;
	      case 'p':	num=atof(argv[++a]);
			break;     
	      case 'q':	qflag=1;
			break;	
              case 'r':	radius=atof(argv[++a]);
			break;   
	      case 's':	boxsiz=atof(argv[++a]);
			break;
	      case 't':	imsiz=atoi(argv[++a]);
			break;
	      }
	
 printf("Reading Files\n");	
		
		if(!(fic1   = fopen(ficnamin1,"r")))
		exit(0);
		if(!(fic2   = fopen(ficnamin2,"r")))
		exit(0);
		if(!(fic3   = fopen(astro,"r")))
		exit(0);
		ficout = fopen(outfic,"wb");


	
/***************************************************************/	

         fgets(s,1000,fic1);

	 sscanf(s,"%lf %lf %lf",&mx,&my,&mag);
         min_x=max_x=mx;
         min_y=max_y=my;


	while(feof(fic1) == 0)
	{
	 fgets(s,1000,fic1);
	 sscanf(s,"%lf %lf %lf",&mx,&my,&mag);
          if(min_x>mx) min_x=mx;
          if(min_y>my) min_y=my;
          if(max_x<mx) max_x=mx;
          if(max_y<my) max_y=my;
        }


  imsiz=ceil((max_x-min_x)*3600.0/pixsize);
  imsiz_y=ceil((max_y-min_y)*3600.0/pixsize);

 printf("%lf %lf %lf %lf imsiz: %i imsiz_y: %i\n",min_x,min_y,max_x,max_y,imsiz,imsiz_y);
 

	
/***************************************************/
/*    Memory allocation using pointers arrays      */
/***************************************************/
	
	  /*number of stars in a box*/

	widbox =  (int)floor((double)imsiz/boxsiz);
	++widbox;
	nbb    = (long)(widbox*(imsiz_y+1)); /*number of boxes*/
	nbelt  = (long)(nbb*nb);
       

if(!qflag)
printf("imsiz: %i number of box: %i widbox %i nbelt %i\n",imsiz,nbb,widbox,nbelt);

	sm   = (ptst *)malloc((size_t)1*(size_t)sizeof(ptst));

/*******for nbtry1 ***********/
	nbtry1 = (int *)malloc((size_t)(nbb*sizeof(int)));

/*******for dataref **********/

	dataref = (ptst *)malloc((size_t)(NBSTARMAX*sizeof(ptst)));

/*******for boxref **********/	

	boxref    = (int **)malloc((size_t)(nbb*sizeof(int *)));
	for (i=0; i< nbb; i++)
	boxref[i] = (int* )malloc((size_t)(nb*sizeof(int)));




	covar = (double **)malloc(100*sizeof(double *));	
	for (i=0; i< 100; i++)
	covar[i] = (double *)malloc((10) * sizeof(double));

	chixy = 0.0;



  
	/**************** Get Offsets ****************/

                i=0;
	        fread(&j,4,1,fic3);
		printf("nbcoeff %i\n",j);
                nbcoeff=j;
                degree=(int)floor(-1.5+sqrt(0.25+2.0*(double)j));
printf("degre: %i\n",degree);
                for(i=0;i<j;i++)
                {
                 fread(&coeffx[i+1],8,1,fic3);
                 printf("coeffx[%i] %lf\n",i,coeffx[i+1]);
                }

                for(i=0;i<j;i++)
                {
                 fread(&coeffy[i+1],8,1,fic3);
                 printf("coeffy[%i] %lf\n",i,coeffy[i+1]);
                }

	

    	/**********************************************/

/******************************************/
/******Box the data for initial File ******/
/******************************************/


/***Initialise number of stars in a box***/

if(!qflag)
printf("initialising\n");
	for(i=0;i<nbb;i++)
        nbtry1[i]=0;

if(!qflag)
printf("boxing\n");

 min_x=min_y=0.0;

	n1   = 0;
        rewind(fic1);
	while(feof(fic1) == 0 && n1<MAXLINE)
	{
		 
	 
	 fgets(s,1000,fic1);
	 sscanf(s,"%lf %lf %lf",&mx,&my,&mag);

         mx=(mx-min_x)*3600.0/pixsize;
         my=(my-min_y)*3600.0/pixsize;   

	 nx = floor(mx/boxsiz);
	 ny = floor(my/boxsiz);
	 nxy= nx+widbox*ny;

if(nxy <0) nxy=0;
if(nxy >(nbb-1)) nxy=nbb-1;


         dataref[n1].x      = mx;
	 dataref[n1].y      = my;
         dataref[n1].mag    = mag;


if(nbtry1[nxy]>(nb-1)) nbtry1[nxy]=nb-1;
         boxref[nxy][nbtry1[nxy]] = n1;

if(nbtry1[nxy]>nb-2) {nbtry1[nxy]=0; }
	 ++nbtry1[nxy];
	 ++n1;

	} 

/******************************************/
/***************Cross files****************/
/******************************************/

if(!qflag)
printf("crossing\n");

	n2=0;
	icross=ic=0;
	inref =0;
        test=0.0;
        ncent=0;
        dcent=10000.0;

	while(feof(fic2) == 0 && n2 < MAXLINE)
	{
	
	 n2++;


		 fgets(s,1000,fic2);

	 sscanf(s,"%lf %lf %lf",&mx,&my,&mag);


         mx=(mx-min_x)*3600.0/pixsize;
         my=(my-min_y)*3600.0/pixsize;   

		


		X = (double)mx;
                Y = (double)my;

            if(nbcoeff==1)
            {
                sm->x=X+coeffx[1];
                sm->y=Y+coeffy[1];
            }

            else
            {
		sm->x =  poly(X,Y,degree,coeffx);
		sm->y =  poly(X,Y,degree,coeffy);
            }
            



	  gcross(dataref,boxref,sm,&findit,boxsiz,widbox,nbtry1,nbb,nb,&find,radius);


	  sg = findit.object;
        
	    if(find)
            {

               if(sg->mag>0.01)   sg->mag  =  2.5*log10(sg->mag);
                 else  sg->mag = 0.01;

               if(mag>0.01)       mag      =  2.5*log10(mag);
                 else  mag = 0.01;

                       
 
                  magc[icross]    =  mag-2.5*log10(4.0);
                  magref[icross]  =  sg->mag-2.5*log10(4.0);

                  color[icross]   =  (magc[icross]-magref[icross]);
                 
                 diffx[icross]   =  (sg->x-sm->x);
                  diffy[icross]   =  (sg->y-sm->y);
                  xref[icross]    =  sg->x;
                  yref[icross]    =  sg->y;
                      if(fitflag)    
                fprintf(ficout,"2 %lf %lf %lf %lf\n",sg->x,sg->y,mx,my);  
                     else
             fprintf(ficout,"2 %lf %lf %lf %lf\n",sg->x,sg->y,sg->mag,mag);  
                ++icross; 

             
	 }
        }

	fclose(fic1);
	fclose(fic2);
	fclose(ficout);


       
 printf("infile1: %i infile2: %i stars found in both files: %i sigma: \n"
    ,--n1,--n2,--icross);

	

	for(i=0;i<(int)nbb;i++)
	{
	 free(boxref[i]);
	}
	for(i=0;i<100;i++)
	{
	 free(covar[i]);
	}


	 free(nbtry1);
	 free(dataref);
	 free(covar);
	 free(boxref);   


}



void    gcross(ptst *dataref,int **boxref,ptst *sm,findstruct *pt,double boxsiz,int widbox,int *nbtry1,long nbb,long nb,int *find,double tol)
{
	int     nx,ny,nxy,test,i,cx,cy;
	double  difx,dify,x,y,xref,yref,ray,ray0;


	
	test = 1;
	ray0=tol;
	*find = 1;
	

	nx  = floor(sm->x/boxsiz);
	ny  = floor(sm->y/boxsiz);
	nxy = nx+widbox*ny;
	if((sm->x - (double)nx*boxsiz-boxsiz/2.0) > 0.0) cx = 1.0;
	else cx = -1.0;
	if((sm->y - (double)ny*boxsiz-boxsiz/2.0) > 0.0) cy = 1.0;
	else cy = -1.0;


	if(nxy<(int)nbb && nxy >-1)
	 {
	i=0;
	while( i< nbtry1[nxy] )
	{

	   difx = sm->x - dataref[boxref[nxy][i]].x;
	   dify = sm->y - dataref[boxref[nxy][i]].y;
	   ray  = sqrt(difx*difx+dify*dify);



	 if(ray<ray0)
         {
	  test=0;
	  ray0=ray;
	  pt->object = &dataref[boxref[nxy][i]];
	  pt->nbref  = boxref[nxy][i];

	 }
	 i += 1;
	}
	}


	

	
	if(test)
	{
	
	 nxy = nx+cx+widbox*ny;

	 if(nxy<(int)nbb && nxy >-1)
	 {
	   i=0;
	   while( i< nbtry1[nxy] )
	  {
	   
	   
	   difx = sm->x - dataref[boxref[nxy][i]].x;
	   dify = sm->y - dataref[boxref[nxy][i]].y;
	   ray  = sqrt(difx*difx+dify*dify);

	   if(ray<ray0)
	   {
	    test=0;
	    ray0=ray;
	    pt->object = &dataref[boxref[nxy][i]];
	    pt->nbref  = boxref[nxy][i];


	    
	   }
	   i += 1;
	  }
	 }
	}



	i=0;
	if(test)
	{
	
	 nxy = nx+widbox*(ny+cy);
	 if(nxy<(int)nbb  && nxy >-1)
	 {
	   while(i< nbtry1[nxy] )
	  {
	   
	   
	   difx = sm->x - dataref[boxref[nxy][i]].x;
	   dify = sm->y - dataref[boxref[nxy][i]].y;
	   ray  = sqrt(difx*difx+dify*dify);
	   if(ray<tol)
	   {
	    test=0;
	    ray0=ray;
	    
	    /*pt = &dataref[boxref[nxy][i]];*/

	    pt->object = &dataref[boxref[nxy][i]];
	    pt->nbref  = boxref[nxy][i];

	   }
	   i += 1;
	  }
         }
	}



	i=0;
	if(test)
	{

	 nxy = nx+cx+widbox*(ny+cy);
	 if(nxy<(int)nbb &&  nxy >-1)
	 { 
	   while(i< nbtry1[nxy] && test>0)
	  {
	   
	   difx = sm->x - dataref[boxref[nxy][i]].x;
	   dify = sm->y - dataref[boxref[nxy][i]].y;
	   ray  = sqrt(difx*difx+dify*dify);
	   if(ray<ray0)
	   {
	    test=0;
	    ray0=ray;
	    pt->object = &dataref[boxref[nxy][i]];
	    pt->nbref  = boxref[nxy][i];

	   }
	   i += 1;
	  }
	 }
	}



	if(test) 
	{	 
         *find = 0;
	}

	

	return ;

}


