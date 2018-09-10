#include <math.h>
#include <malloc.h>
#include <stdio.h>

#define MOBJ 5001
#define MTRI 25005
#define NPAR 7
#define MERRMAX 0.1
#define ERRMAX 1.0
#define ROTERR 0.02


main (argc,argv)
     int argc;
     char *argv[];
{
  
  FILE *stream1,*stream2,*stream3;

  float *x1, *y1, *amag1, *rsum1, *x2,*y2, *amag2,*rsum2;
  float *xm1, *ym1,*xm2,*ym2,*xf1,*yf1,*xf2,*yf2,*xh1,*yh1,*xh2,*yh2,*amagdiff;
  float *ratiolist,*ratsky,*sky1,*sky2,*diffx,*diffy,*xpos1,*ypos1,*mag1,*emag1;
  float **r2,**r1,**rabs1,**rabs2;
  int *index1, *indexr1,  *index2,*indexr2,*index, **ijk1,**ijk2;
  
  float *vector(),**matrix();
  int   *ivector(),**imatrix(), matchstars();
  void indexx();


  int nobj,nstart1,nstart2,i,j,k,l,m,ii,jj,kk,ll,mm,itype,nstar1,nstar2,q,r,rmax,firsttime;
  int jbest,nm,nf,irepeat,anf,igotavg,iend,id,imtype,nsum1,nsum2,jlow,jhigh;
  float errmax,sepmin,tolerance,rtolerance,sigcut,rsummin; 
  float err,r12,r13,r23,xhunt,tol,delta,temp;
  float sum,sumrs,xtest,ytest,sumr2,sumr1,avgdiff,gamma;
  float ax,bx,cx,ay,by,cy,sumdx,sumdy,sumdx2,sumdy2;
  float sdx,sdy,xcut,ycut,XC,YC,xc,yc,TILT,AMINOR,AMAJOR,AREA,FMAG;
  float xx,yy,bmag,scaleratio,srtolerance,ratio,x,y,dummy;
  float ratiomap[10],starpar[NPAR],staravg[NPAR];
  float AMINOR1,AMAJOR1,TILT1;
  char ians,line[150],file1[80],file2[80],file3[80],file4[80];
  int ns1,ns2,nstar,nbr,nbl,nbb,nbt,gotit=0;
  float mag2,emag2,FERR,xa1,ya1,xa2,ya2,xc1,xc2,yc1,yc2;


  void elliptint(),transform();

  if ( argc != 9 ) {
     printf( "Usage: %s templist objlist ax bx cx ay by cy\n", argv[0] );
     exit( 1 );
  }

  
  x1=vector(0,MOBJ);
  y1=vector(0,MOBJ);
  rsum1=vector(0,MTRI);
  x2=vector(0,MOBJ);
  y2=vector(0,MOBJ);
  rsum2=vector(0,MTRI);
  xm1=vector(0,MTRI);
  ym1=vector(0,MTRI);
  xm2=vector(0,MTRI);
  ym2=vector(0,MTRI);
  xf1=vector(0,MTRI);
  yf1=vector(0,MTRI);
  xf2=vector(0,MTRI);
  yf2=vector(0,MTRI);
  xh1=vector(0,MTRI);
  yh1=vector(0,MTRI);
  xh2=vector(0,MTRI);
  yh2=vector(0,MTRI);
  xpos1=vector(0,MTRI);
  ypos1=vector(0,MTRI);
  ratiolist=vector(0,MTRI);
  ratsky=vector(0,MOBJ);
  sky1=vector(0,MOBJ);
  sky2=vector(0,MOBJ);
  diffx=vector(0,MOBJ);
  diffy=vector(0,MOBJ);
  index1=ivector(0,MOBJ);
  indexr1=ivector(0,MTRI);
  index2=ivector(0,MOBJ);
  indexr2=ivector(0,MTRI);
  index=ivector(0,MOBJ);
  r2=matrix(0,3,0,MTRI);
  r1=matrix(0,3,0,MTRI);
  rabs1=matrix(0,3,0,MTRI);
  rabs2=matrix(0,3,0,MTRI);
  ijk1=imatrix(0,3,0,MTRI);
  ijk2=imatrix(0,3,0,MTRI);  



  nobj =  50;         /* Max # Objects used for matching */
  errmax = 0.25;      /* Max phot error used for stars in match */
  sepmin = 20.0;      /* Min Separation in pixels for Triangle */
  tolerance = 0.001;  /* Tolerance in Similarity of Triangle */
  rtolerance = 0.001; /* Tolerance in Similarity of Triangle */
  sigcut = 3.0;       /* Get Rid of outliers which are sigcut out*/
  rsummin = 1.10;     /* Min obliqueness of triangle */
  

  ns2=10000;
  

  sprintf(file1,argv[1]);
  sprintf(file2,argv[2]);
  bx = atof( argv[3] );
  cx = atof( argv[4] );
  ax = atof( argv[5] );
  by = atof( argv[6] );
  cy = atof( argv[7] );
  ay = atof( argv[8] );

  if ((stream2 = fopen(file2,"r"))== NULL) {
    printf("\nFile %s Not found...exiting before I dump the core\n", file2);
    exit(0);
  }
  
  if ((stream1 = fopen(file1,"r"))== NULL) {
    printf("File %s Not found - exiting...\n",file1);
    exit(0);
  } 
  
  
  /* load up first file */ 
  for (i=1; feof(stream1) == 0 && i < 5000; i++) {
    fgets(line,150, stream1); 
    if (sscanf(line,"%f %f",&x1[i],&y1[i]) <1) i--;
  }
  nstar1 = i-1;
  if (nstar1 < 3) {
     exit( 1 );
  }
  
  /* load up 2nd file */
  for (i=1; feof(stream2) == 0 && i < 5000; i++) {
     fgets(line,150, stream2); 
     if (sscanf(line,"%f %f",&x2[i],&y2[i]) < 1 ) i--;
  }
  nstar2 = i-1;
  
  fclose(stream1);
  fclose(stream2);

  if (ax!=0 || bx != 0) {
  nf = matchstars( x1, y1, nstar1, x2, y2, nstar2, xf1, yf1, xf2, yf2,
	      ax, bx, cx, ay, by, cy, ERRMAX );
  }
  else { /* no transform, assume stars are pre-aligned */
    if (nstar1 != nstar2) {
      printf ("nstar1 != nstar2, and no transformation supplied\n\n");
      exit (1);
    }
    for (i=1;i<=nstar1;i++) {
      xf1[i]=x1[i];
      yf1[i]=y1[i];
      xf2[i]=x2[i];
      yf2[i]=y2[i];
    }
    nf=nstar1;
  }
  
  /* Do Multvariate regression to find tranformation coefficients */
  irepeat =1;
  while (irepeat == 1 && nf > 2) {
     transform(nf,xf1,yf1,xf2,yf2,&ax,&bx,&cx,&ay,&by,&cy);

     sumdx = 0.0;
     sumdy = 0.0;
     sumdx2 = 0.0;
     sumdy2 = 0.0;
	
     for (i=1; i<=nf;i++) {       /* Calc sigma */
	x = ax +bx*xf2[i] + cx*yf2[i];
	y = ay +by*xf2[i] + cy*yf2[i];
	
	diffx[i] = x - xf1[i];
	diffy[i] = y - yf1[i];
	  
	sumdx += diffx[i];
	sumdx2 += diffx[i]*diffx[i];
	sumdy += diffy[i];
	sumdy2 += diffy[i]*diffy[i];
	
	xh1[i] = xf1[i];
	yh1[i] = yf1[i];
	xh2[i] = xf2[i];
	yh2[i] = yf2[i];
     }  
     anf = nf;
     sdx = sqrt((1.0/(anf-1.0))*(sumdx2-(sumdx*sumdx/anf)));
     sdy = sqrt((1.0/(anf-1.0))*(sumdy2-(sumdy*sumdy/anf)));
	
     xcut = sigcut*sdx;
     ycut = sigcut*sdy;
	
     kk=1;
	
     irepeat = 0;
	
     for (i=1; i<=nf;i++) { /* cull those who are too far out */
	if (fabs(diffx[i]) > xcut || fabs(diffy[i]) > ycut) irepeat = 1;
	else {
	   xf1[kk] = xh1[i];
	   yf1[kk] = yh1[i];
	   xf2[kk] = xh2[i];
	   yf2[kk] = yh2[i];
	   kk ++;
	}
     }	
     if(irepeat == 1) nf = kk-1; /* got rid of a star, iterate */
  } /*End of WHile */
      
  printf("%f %f %f %f %f %f\n",bx,cx,ax,by,cy,ay);

}


void transform(nstar,x2,y2,x1,y1,ax,bx,cx,ay,by,cy)   
     int nstar;
     float x2[],y2[],x1[],y1[],*ax,*bx,*cx,*ay,*by,*cy;
{
  int i;
  double sum=0.,sumx=0.,sumy=0.,sumxy=0.,sumu=0.,sumv=0.,sumx2=0.,sumy2=0.;
  double sumux=0.,sumuy=0.,sumvx=0.,sumvy=0.;
  double delta,sumdx,sumdy,sumabsdx,sumabsdy;
  double x1p,y1p,dx,dy;
  for (i=1; i<=nstar;i++) {
    sum ++;
    sumx += x1[i];
    sumy += y1[i];
    sumu += x2[i];
    sumv += y2[i];
    sumx2 += x1[i]*x1[i];
    sumxy += x1[i]*y1[i];
    sumy2 += y1[i]*y1[i];
    sumux += x1[i]*x2[i];
    sumvx += x1[i]*y2[i];
    sumuy += y1[i]*x2[i];
    sumvy += y1[i]*y2[i];
  }
  delta = sum*(sumx2*sumy2 - sumxy*sumxy) + sumx*(sumxy*sumy - sumx*sumy2) + sumy*(sumx*sumxy-sumx2*sumy);
  *ax = (sumu*(sumx2*sumy2 - sumxy*sumxy) + sumx*(sumxy*sumuy - sumux*sumy2) + sumy*(sumux*sumxy - sumx2*sumuy))/delta;
  *bx = (sum*(sumux*sumy2 - sumxy*sumuy) + sumu*(sumxy*sumy - sumx*sumy2) + sumy*(sumuy*sumx - sumy*sumux))/delta;
  *cx = (sum*(sumx2*sumuy - sumxy*sumux) + sumx*(sumux*sumy - sumx*sumuy) + sumu*(sumx*sumxy - sumx2*sumy))/delta;
  *ay = (sumv*(sumx2*sumy2 - sumxy*sumxy) + sumx*(sumxy*sumvy - sumvx*sumy2) + sumy*(sumvx*sumxy - sumx2*sumvy))/delta;
  *by = (sum*(sumvx*sumy2 - sumxy*sumvy) + sumv*(sumxy*sumy - sumx*sumy2) + sumy*(sumvy*sumx - sumy*sumvx))/delta;
  *cy = (sum*(sumx2*sumvy - sumxy*sumvx) + sumx*(sumvx*sumy - sumx*sumvy) + sumv*(sumx*sumxy - sumx2*sumy))/delta;
  
  sumdx = 0.0;
  sumdy = 0.0;
  sumabsdx = 0.0;
  sumabsdy = 0.0;
  
  for (i=1;i<=nstar;i++) {
    x1p = *ax + *bx * x1[i] + *cx * y1[i];
    y1p = *ay + *by * x1[i] + *cy * y1[i];
    dx = x2[i] - x1p;
    dy = y2[i] - y1p;
    sumdx = sumdx + dx;
    sumdy = sumdy + dy;
    sumabsdx = sumabsdx + abs(dx);
    sumabsdy = sumabsdy + abs(dy);
  }
  sumdx /= nstar;
  sumdy /= nstar;
  sumabsdx /= nstar;
  sumabsdy /= nstar;
}

void elliptint(AMAJOR,AMINOR,TILT,AREA,starpar)
     float starpar[],*AMAJOR,*AMINOR,*TILT,*AREA;
{
  float ROOT,ROOT1,ROOT2,A7MNA5,A7PLA5,A7,A5;
  
  if (*AMAJOR < 0.01) *AMAJOR =0.01;
  if (*AMINOR < 0.01) *AMINOR = 0.01;
  
  ROOT1 = (2.35482 / (*AMINOR))*(2.35482 / (*AMINOR));
  ROOT2 = (2.35482 / (*AMAJOR))*(2.35482 / (*AMAJOR));
  
  *AREA = 6.2832 / sqrt(ROOT1*ROOT2);
  
  starpar[6] = 0.5*(ROOT2 - ROOT1) * sin(2*(*TILT/57.29578));
  A7MNA5 = (ROOT1-ROOT2) * cos(2*(*TILT/57.29578));
  A7PLA5 = ROOT+ROOT2;
  
  A7 = 0.5*(A7PLA5 + A7MNA5);
  A5 = 0.5*(A7PLA5 - A7MNA5);
  
  starpar[5] = 1./A5;
  starpar[7] = 1./A7;
}


int matchstars( x1, y1, n1, x2, y2, n2, xf1, yf1, xf2, yf2, ax, bx, cx,
		 ay, by, cy, limit )
     float *x1, *y1, *x2, *y2, *xf1, *yf1, *xf2, *yf2, ax, bx, cx, ay, by, cy, limit;
     int n1, n2;
{
   int i, j, nout, bestj;
   float x, y, dist, bestx, besty, bestdist;

   nout = 0;
   for ( i = 1; i <= n1; i++ ) {
      bestdist = 999.99;
      for ( j = 1; j <= n2; j++ ) {
	 x = x2[j] * bx + y2[j] * cx + ax;
	 y = x2[j] * by + y2[j] * cy + ay;
	 dist = sqrt( ( x1[i] - x ) * ( x1[i] - x ) + ( y1[i] - y ) * ( y1[i] - y ) );
	 if ( dist < bestdist ) {
	    bestdist = dist;
	    bestj = j;
	 }
      }
      if ( bestdist < limit ) {
	 nout ++;
	 xf1[nout] = x1[i];
	 yf1[nout] = y1[i];
	 xf2[nout] = x2[bestj];
	 yf2[nout] = y2[bestj];
      }
   }
   return nout;
}


void nrerror(error_text)
char error_text[];
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}



float *vector(nl,nh)
int nl,nh;
{
	float *v;

	v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

int *ivector(nl,nh)
int nl,nh;
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

double *dvector(nl,nh)
int nl,nh;
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}



float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i;
	float **m;

	m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
	int i,**m;

	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}



float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
float **a;
int oldrl,oldrh,oldcl,oldch,newrl,newcl;
{
	int i,j;
	float **m;

	m=(float **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(float*));
	if (!m) nrerror("allocation failure in submatrix()");
	m -= newrl;

	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	return m;
}



void free_vector(v,nl,nh)
float *v;
int nl,nh;
{
	free((char*) (v+nl));
}

void free_ivector(v,nl,nh)
int *v,nl,nh;
{
	free((char*) (v+nl));
}

void free_dvector(v,nl,nh)
double *v;
int nl,nh;
{
	free((char*) (v+nl));
}



void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}



void free_submatrix(b,nrl,nrh,ncl,nch)
float **b;
int nrl,nrh,ncl,nch;
{
	free((char*) (b+nrl));
}



float **convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
int nrl,nrh,ncl,nch;
{
	int i,j,nrow,ncol;
	float **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	m = (float **) malloc((unsigned) (nrow)*sizeof(float*));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m -= nrl;
	for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
	return m;
}



void free_convert_matrix(b,nrl,nrh,ncl,nch)
float **b;
int nrl,nrh,ncl,nch;
{
	free((char*) (b+nrl));
}
void indexx(n,arrin,indx)
int n,indx[];
float arrin[];
{
	int l,j,ir,indxt,i;
	float *temp,q;
	void free_vector();
	float *vector();

	temp = vector(1,n);
	

	for (j=1;j<=n;j++) indx[j]=j;
	l=(n >> 1) + 1;
	ir=n;
	for (;;) {
		if (l > 1)
			q=arrin[(indxt=indx[--l])];
		else {
			q=arrin[(indxt=indx[ir])];
			indx[ir]=indx[1];
			if (--ir == 1) {
				indx[1]=indxt;
				for (j=1;j<=n;j++) temp[j]=arrin[j];
				for (j=1;j<=n;j++) arrin[j] = temp[indx[j]];
				free_vector(temp,1,n);
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && arrin[indx[j]] < arrin[indx[j+1]]) j++;
			if (q < arrin[indx[j]]) {
				indx[i]=indx[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		indx[i]=indxt;
	}
}

