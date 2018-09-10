#include<stdio.h>
#include<math.h>
#define MOBJ 5001
#define MTRI 100000
#define NPAR 7
#define MERRMAX 1
#define ERRMAX 2
#define SQUARE 0.04

main (argc,argv)
     int argc;
     char *argv[];
{
  
  FILE *stream1,*stream2,*stream3,*stream4;

  float *x1, *y1, *amag1, *rsum1, *x2,*y2, *amag2,*rsum2;
  float *xm1, *ym1,*xm2,*ym2,*xf1,*yf1,*xf2,*yf2,*xh1,*yh1,*xh2,*yh2,*amagdiff;
  float *ratiolist,*ratsky,*sky1,*sky2,*diffx,*diffy,*xpos1,*ypos1,*mag1,*emag1;
  float *alist,*blist,*clist,*dlist,*elist,*flist;
  float **r2,**r1,**rabs1,**rabs2;
  int *index1, *indexr1,  *index2,*indexr2,*index, **ijk1,**ijk2;
  float *vector(),**matrix();
  int   *ivector(),**imatrix();



  int nobj,nstart1,nstart2,i,j,k,l,m,ii,jj,kk,ll,mm,itype,nstar1,nstar2,q,r=1,rmax,firsttime;
  int iii,jjj,kkk,rnew;
  int jbest,nm,nf,irepeat,anf,igotavg,iend,id,imtype,nsum1,nsum2,jlow,jhigh;
  int STARMATCHES,skiprow;
  float abest=0.,bbest=0.,cbest=0.,dbest=0.,ebest=0.,fbest=0.,amiss,cbmiss,a;
  float errmax,sepmin,tolerance,rtolerance,sigcut,rsummin; 
  float err,r12,r13,r23,xhunt,tol,delta,temp;
  float sum,sumrs,xtest,ytest,sumr2,sumr1,avgdiff,gamma;
  float ax,bx,cx,ay,by,cy,sumdx,sumdy,sumdx2,sumdy2;
  float sdx,sdy,xcut,ycut,XC,YC,xc,yc,TILT,AMINOR,AMAJOR,AREA,FMAG;
  float xx,yy,bmag,scaleratio,srtolerance,ratio,x,y,dummy;
  float ratiomap[10],starpar[NPAR],staravg[NPAR];
  float AMINOR1,AMAJOR1,TILT1,ROTERR;
  char ians,line[150],file1[80],file2[80],file3[80],file4[80],dummyc[1000];
  int ns1,ns2,nstar,nbr,nbl,nbb,nbt,gotit=0,diagnose;
  float mag2,emag2,FERR,xa1,ya1,xa2,ya2,xc1,xc2,yc1,yc2;


  void elliptint(),transform();

  
  x1=vector(0,MOBJ);
  y1=vector(0,MOBJ);
  amag1=vector(0,MOBJ);
  rsum1=vector(0,MTRI);
  x2=vector(0,MOBJ);
  y2=vector(0,MOBJ);
  amag2=vector(0,MOBJ);
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
  mag1=vector(0,MTRI);
  emag1=vector(0,MTRI);
  amagdiff=vector(0,MTRI);
  ratiolist=vector(0,MTRI);
  alist=vector(0,MTRI);
  blist=vector(0,MTRI);
  clist=vector(0,MTRI);
  dlist=vector(0,MTRI);
  elist=vector(0,MTRI);
  flist=vector(0,MTRI);
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


  if (argc <2) {
    printf ("starmatch: objfiletemplate objfileimage scaleratio scaletol skipcol(0 1) STARMATCHES nob amiss cbmiss debug\n\nobjfiles have\n    x y mag\non each line\nif scaleratio is 0, finds scale\n");
    exit(0);
  }

  stream4 = fopen("param_dophot","w");          

  ns2=10000;
  

  nobj =  50;         /* Max # Objects used for matching */

  sprintf(file1,argv[1]);
  sprintf(file2,argv[2]);
  sscanf(argv[3],"%f",&scaleratio);
  sscanf(argv[4],"%f",&srtolerance);
  if (argc>5)     sscanf(argv[5],"%d",&skiprow);
  else skiprow=0;
  if (argc>6)     sscanf(argv[6],"%d",&STARMATCHES);
  else STARMATCHES=2;
  if (argc>7)     sscanf(argv[7],"%d",&nobj);
  if (argc>8)     sscanf(argv[8],"%f",&amiss);
  else amiss=5;
  if (argc>9)     sscanf(argv[9],"%f",&cbmiss);
  else cbmiss=0.05;
  if (argc>10) sscanf(argv[10],"%f",&ROTERR);
  else ROTERR=9999;
  if (argc>11) diagnose=1;
  else diagnose=0;





  errmax = 1.5;      /* Max phot error used for stars in match */
  sepmin = 8.0;      /* Min Separation in pixels for Triangle */
  tolerance = 0.005;  /* Tolerance in Similarity of Triangle */
  rtolerance = 0.005;   /* Tolerance in Similarity of Triangle */
  sigcut = 3.0;       /* Get Rid of outliers which are sigcut out*/
  rsummin = 1.10;     /* Min obliqueness of triangle */


   
  if ((stream2 = fopen(file2,"r"))== NULL) {
    printf("\nTemplate File Not found...exiting before I dump the core\n");
    exit(0);
  }
  
  if ((stream1 = fopen(file1,"r"))== NULL) {
    printf("File %s Not found - exiting...\n",file1);
    exit(0);
  } 
  
  
  /* load up first file */ 
  for (i=1; feof(stream1) == 0 && i < 5000; i++) {
    fgets(line,150, stream1); 
    if(strncmp(line, "#", 1) != 0) {
      if (skiprow == 0)  sscanf(line,"%f %f %f",&x1[i],&y1[i],&amag1[i]);
      else sscanf(line,"%s %f %f %f",dummyc,&x1[i],&y1[i],&amag1[i]);
      if (amag1[i] ==0) amag1[i]=99.99;
      if (x1[i]==0 && y1[i]==0) i--;
    }
    else i--;
  }
  nstar1 = i-1;
  if (nstar1 < 3)   printf("\n%s only has %d Star/s...skipping it\n",file1,nstar1);
  else {

    indexx(nstar1,amag1,index1); /* sort by magnitudes  - give sort into index1 */
    
    /* load up 2nd file */
    for (i=1; feof(stream2) == 0 && i < 5000; i++) {
      fgets(line,150, stream2); 
      if(strncmp(line, "#", 1) != 0) {
	if (skiprow == 0)        sscanf(line,"%f %f %f",&x2[i],&y2[i],&amag2[i]);
	else  sscanf(line,"%s %f %f %f",dummyc,&x2[i],&y2[i],&amag2[i]);
      if (amag2[i] ==0) amag2[i]=99.99;
      if (x2[i]==0 && y2[i]==0) i--;
      }
      else i--;
    }
    nstar2 = i-1;
    
    indexx(nstar2,amag2,index2); /* sort by magnitudes  - give sort into index1 */
    
    fclose(stream1);
    fclose(stream2);
    
    
    /* Calculate Triangle Ratios */
    
    for (nstart2=1;nstart2 < nstar2 && gotit !=1;nstart2+=nobj) {          
      for (nstart1=1;nstart1 < nstar1 && gotit !=1 ;nstart1 +=nobj) {
	/* 1st file */ 
	for (kk=1,i=nstart1; i<nstart1+nobj-2-1 && i <= nstar1; i++) {
	  for (j=i+1; j < nstart1+nobj-1-1 && j <= nstar1; j++) {
	    r12 = sqrt(((x1[index1[i]]-x1[index1[j]])*(x1[index1[i]]-x1[index1[j]])) +
		       ((y1[index1[i]]-y1[index1[j]])*(y1[index1[i]]-y1[index1[j]])));
	    if (r12 > sepmin && i != j) {
	      for(k=j+1;k<nstart1+nobj-1 && k <= nstar1;k++) {
		r13 = sqrt(((x1[index1[i]]-x1[index1[k]])*(x1[index1[i]]-x1[index1[k]])) +
			   ((y1[index1[i]]-y1[index1[k]])*(y1[index1[i]]-y1[index1[k]])));
		if (r13 > sepmin && i != k ) {
		  r23 = sqrt(((x1[index1[j]]-x1[index1[k]])*(x1[index1[j]]-x1[index1[k]])) 
			     + ((y1[index1[j]]-y1[index1[k]])*(y1[index1[j]]-y1[index1[k]])));
		  if (r23 > sepmin && j != k ) {
		    if (r12 > r13) { /*sort from to small large*/
		      temp=r12;
		      r12=r13;
		      r13=temp;
		    }
		    
		    if (r12 > r23) {		  
		      temp=r12 ;
		      r12=r23;
		      r23=temp;
		    }
		    
		    if (r13 > r23) {		  
		      temp=r13;
		      r13=r23;
		      r23=temp;
		    }
		    
		    
		    r1[1][kk] = r12/r13;
		    r1[2][kk] = r12/r23;
		    r1[3][kk] = r13/r23;
		    
		    rsum1[kk] = r1[1][kk]+r1[2][kk]+r1[3][kk];
		    
		    if (rsum1[kk] > rsummin) {
		      rabs1[1][kk] = r12;
		      rabs1[2][kk] = r13;
		      rabs1[3][kk] = r23;
		      ijk1[1][kk] = i;
		      ijk1[2][kk] = j;
		      ijk1[3][kk] = k;
		      kk ++;
		    }
		  }
		}
	      }
	    }
	  }
	}
	nsum1 = kk - 1;
	if (nsum1 <2) nsum1 =2; 
	indexx(nsum1,rsum1,indexr1);
	
	/* 2nd file */
	
	for (kk=1,i=nstart2; i<nstart2+nobj-2-1 && i <= nstar2; i++) {
	  for (j=i+1; j < nstart2+nobj-1-1 && j <=nstar2; j++) {
	    r12 = sqrt(((x2[index2[i]]-x2[index2[j]])*(x2[index2[i]]-x2[index2[j]])) +
		       ((y2[index2[i]]-y2[index2[j]])*(y2[index2[i]]-y2[index2[j]])));
	    if (r12 > sepmin && i != j) {
	      for(k=j+1;k<nstart2+nobj-1 && k <= nstar2;k++) {
		r13 = sqrt(((x2[index2[i]]-x2[index2[k]])*(x2[index2[i]]-x2[index2[k]])) +
			   ((y2[index2[i]]-y2[index2[k]])*(y2[index2[i]]-y2[index2[k]])));
		if (r13 > sepmin && i !=k) {
		  r23 = sqrt(((x2[index2[j]]-x2[index2[k]])*(x2[index2[j]]-x2[index2[k]])) +
			     ((y2[index2[j]]-y2[index2[k]])*(y2[index2[j]]-y2[index2[k]])));
		  if (r23 > sepmin && j !=k) {
		    if (r12 > r13) { /*sort from to small large*/
		      temp=r12;
		      r12=r13;
		      r13=temp;
		    }
		    
		    if (r12 > r23) {		  
		      temp=r12 ;
		      r12=r23;
		      r23=temp;
		    }
		    
		    if (r13 > r23) {		  
		      temp=r13;
		      r13=r23;
		      r23=temp;
		    }
		    
		    
		    r2[1][kk] = r12/r13;
		    r2[2][kk] = r12/r23;
		    r2[3][kk] = r13/r23;
		    
		    rsum2[kk] = r2[1][kk]+r2[2][kk]+r2[3][kk];
		    if (rsum2[kk] > rsummin) {
		      rabs2[1][kk] = r12;
		      rabs2[2][kk] = r13;
		      rabs2[3][kk] = r23;
		      ijk2[1][kk] = i;
		      ijk2[2][kk] = j;
		      ijk2[3][kk] = k;
		      kk ++;
		    }
		  }
		}
	      }
	    }
	  }
	}
	nsum2 = kk -1;
	if (nsum2 < 2) nsum2 =2; 
	indexx(nsum2,rsum2,indexr2);
	
	firsttime = 0;
	while ((firsttime==0) || (firsttime== 1 && abest!=0 && bbest !=0 && cbest !=0)) {
	  rnew=0;
	  firsttime ++;
	  jlow = nsum2/2;
	  jhigh = jlow;
	  kk =1;
	  for (i=1;i <= nsum2 && rnew<100; i++) {
	    
	    xhunt = rsum2[i] - tolerance; /*Search list in template file for closest */
	    hunt(rsum1,nsum1,xhunt,&jlow); /* Match to the sums in active file */
	    
	    if (jlow > 1) jlow--;
	    xhunt = rsum2[i] + tolerance;
	    hunt(rsum1,nsum1,xhunt,&jhigh);
	    if (jhigh < nsum1) jhigh++;
	    
	    tol = tolerance;
	    jbest = 0;
	    
	    for (jbest=jlow; jbest <= jhigh+1; jbest++) {  /* find best in the range give by Hunt */
	      
	      if (jbest > 0) {
		ll=0;
		mm=0;
		for (l=1;l<=3;l++) {   /* Do triangles match */
		  delta = fabs(r1[l][indexr1[jbest]]-r2[l][indexr2[i]]);
		  if (delta <= rtolerance) ll ++;
		}
		if (ll >= 3 ) {  /* Do triangles have 0 or 90degree Rotationo? */
		  ll=3;
		  
		  xf1[1] = x1[index1[ijk1[1][indexr1[jbest]]]];
		  yf1[1] = y1[index1[ijk1[1][indexr1[jbest]]]];
	      
		  xf1[2] = x1[index1[ijk1[2][indexr1[jbest]]]];
		  yf1[2] = y1[index1[ijk1[2][indexr1[jbest]]]];
		      
		  xf1[3] = x1[index1[ijk1[3][indexr1[jbest]]]];
		  yf1[3] = y1[index1[ijk1[3][indexr1[jbest]]]];

		  for (iii=1;iii<=3;iii++) {
		    for (jjj=1;jjj<=3;jjj++) {
		      for (kkk=1;kkk<=3;kkk++) {
			if (iii !=jjj && iii!=kkk && jjj !=kkk) {
			
			  xf2[iii] = x2[index2[ijk2[1][indexr2[i]]]];
			  yf2[iii] = y2[index2[ijk2[1][indexr2[i]]]];
			  xf2[jjj] = x2[index2[ijk2[2][indexr2[i]]]];
			  yf2[jjj] = y2[index2[ijk2[2][indexr2[i]]]];
			  xf2[kkk] = x2[index2[ijk2[3][indexr2[i]]]];
			  yf2[kkk] = y2[index2[ijk2[3][indexr2[i]]]];
			  
			  transform(3,xf1,yf1,xf2,yf2,&ax,&bx,&cx,&ay,&by,&cy);
			  
			  if ( fabs((bx*bx+cx*cx)/(by*by+cy*cy) - 1.) < SQUARE && fabs(fabs(bx)-fabs(cy)) < SQUARE && fabs(fabs(cx)-fabs(by)) < SQUARE && (bx*bx/(cx*cx+0.000001) < ROTERR || cx*cx/(bx*bx+0.000001) < ROTERR) ) {
			    iii=999;
			    if (scaleratio != 0) {  /* Scale ratio defined, check!*/
			      ratio = (sqrt(bx*bx+cx*cx)+sqrt(by*by+cy*cy))/2;
			      if (fabs(ratio - scaleratio) <= srtolerance){
				if (diagnose ==1)  printf("%7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n",ax,bx,cx,ay,by,cy,xf1[1],xf2[1]);
				if (abest ==0.0 && bbest==0.0 && cbest==0.0 && r < MTRI/5) { /*find best*/
				  alist[r]=ax;
				  blist[r]=bx;
				  clist[r]=cx;
				  dlist[r]=ay;
				  elist[r]=by;
				  flist[r]=cy;
				  r++;
				  rnew++;
				  ll=0;
				}
				else { /*2nd iteration see if this fits*/
				  if (fabs(ax-abest) < amiss && fabs(bx-bbest)<cbmiss &&
				      fabs(cx-cbest)<cbmiss && fabs(ay-dbest) < amiss && fabs(by-ebest)<cbmiss &&
				      fabs(cy-fbest)<cbmiss) {
				    ll=3; 
				  }
				  else ll=0; 
				}
			      }
			      else ll=0;
			    }
			    else { /*no prior on scale, do the best we can without this info*/
			      if (abest ==0.0 && bbest==0.0 && cbest==0.0) { /*find best*/
				alist[r]=ax;
				blist[r]=bx;
				clist[r]=cx;
				rnew++;
				r++;
				ll=0;
			      }
			      else {/*2nd iteration*/
				if (fabs(ax-abest) < amiss && fabs(bx-bbest)<cbmiss &&
				    fabs(cx-cbest)<cbmiss) {
				  ll=3;
				}
				else ll=0;
			      }
			    }
			    if (ll == 3 && kk < 1000) {  /* We have a matching triangle, Make a list */
			
			      xm1[kk] = xf1[1];
			      ym1[kk] = yf1[1];
			      xm2[kk] = xf2[1];
			      ym2[kk] = yf2[1];
			      kk++;
			      xm1[kk] = xf1[2];
			      ym1[kk] = yf1[2];
			      xm2[kk] = xf2[2];
			      ym2[kk] = yf2[2];
			      kk++;
			      xm1[kk] = xf1[3];
			      ym1[kk] = yf1[3];
			      xm2[kk] = xf2[3];
			      ym2[kk] = yf2[3];
			      kk++;
			    }
			  }
			}
		      }
		    }
		  }		     
		}	      
	      }
	    }
	  }
	  if (abest ==0 && bbest ==0 && cbest==0 && r !=0) {
	    for (rmax=STARMATCHES-1,l=1; l < r; l++) {
	      for (q=0,m=1; m < r && rmax<20; m++) {
		if ((fabs(alist[l]-alist[m]) <= amiss) && (fabs(blist[l]-blist[m]) <= cbmiss) && (fabs(clist[l]-clist[m]) <= cbmiss) && (fabs(dlist[l]-dlist[m]) <= amiss) && (fabs(elist[l]-elist[m]) <= cbmiss) && (fabs(flist[l]-flist[m]) <= cbmiss) && (l!=m)) {
		  if (diagnose ==1)  printf("match[%d]->%7.3f %7.3f %7.3f =",q+1,alist[l],blist[l],clist[l]);
		  if (diagnose ==1)  printf("       %7.3f %7.3f %7.3f \n",alist[m],blist[m],clist[m]);
		  q++;
		}
	      }
	      if (q >= rmax) {
		rmax = q;
		abest = alist[l];
		bbest = blist[l];
		cbest = clist[l];
		dbest = dlist[l];
		ebest = elist[l];
		fbest = flist[l];
		
	      }
	    }
	  }	      
	  
	  if (abest != 0  && cbest !=0 && kk > STARMATCHES*3-1) gotit=1;
	}
      }
    }
    
    nm = kk -1;
    if (nm >= 3) {
      indexx(nm,ym1,index); /* sort by y */
      
      xf1[1]=xtest = xm1[index[1]];
      yf1[1]=ytest = ym1[1];
      xf2[1]=xm2[index[1]];         /* Get rid of stars in list more than once */
      yf2[1]=ym2[index[1]];
      kk = 2;
      
      for (i=2;i <=nm; i++) {
	if (xm1[index[i]] != xtest && ym1[i] != ytest) {
	  xf1[kk] = xm1[index[i]];
	  yf1[kk] = ym1[i];
	  xf2[kk] = xm2[index[i]];
	  yf2[kk] = ym2[index[i]];
	  
	  xtest = xf1[kk];
	  ytest = yf1[kk];
	  if (diagnose ==1)  printf("star[%d]->%7.3f %7.3f %7.3f %7.3f\n",kk,xf1[kk],yf1[kk],xf2[kk],yf2[kk]);
	  kk ++;
	}
      }
      nf = kk - 1;
      
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
      
      stream1 = fopen(file1,"r");
      stream2 = fopen(file2,"r");
      stream3 = fopen("matchstar.out","w");          
      /*      fprintf(stream3,"# Scaleratio=%f\n",scaleratio); */
      kk = 1;
      igotavg = 0;
      
      iend = 0;
      i=1; 
      
      /* read in file 1 again*/
      while (feof(stream1) ==0) {
	fgets(line, 150, stream1); 
	sscanf(line,"%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f"  
	       ,&id,&itype,&XC,&YC,&FMAG,&FERR,&starpar[1],&AMAJOR,
	       &AMINOR,&TILT,&dummy,&dummy,&dummy,&dummy,&dummy);
	xpos1[i] = XC ;
	ypos1[i] = YC ;
	mag1[i] = FMAG;
	emag1[i] = FERR;
	if (itype==1 || itype == 11 && FERR < MERRMAX) /* good star, keep */
	  i++ ;
      }
      ns1 = i;
      
      for (i=1;feof(stream2) == 0; i++ ) {  /*read in second file again*/
	fgets(line,150, stream2); 
	sscanf(line,"%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f"
	       ,&id,&itype,&xc,&yc,&FMAG,&FERR,&starpar[1],&AMAJOR,
	       &AMINOR,&TILT,&dummy,&dummy,&dummy,&dummy,&dummy);
	
	x = xc;
	y = yc;
	mag2 = FMAG;
	emag2 = FERR;
	XC = ax + bx*x +cx*y;  
	YC = ay + by*x + cy*y; /* determine this star's position in other file */
	
	for (j=1; j < ns1; j++) { /* does it match a star in file 1?, if yes, write it out*/
	  if ( sqrt((XC-xpos1[j])*(XC-xpos1[j]) + (YC-ypos1[j])*(YC-ypos1[j])) < ERRMAX && itype == 1 ) {
	    fprintf(stream3,"\n%9.2f%9.2f%9.2f%9.2f",xc,yc,xpos1[j],ypos1[j]);
	    j = 99999;
	  }
	} 
      }  
      /*      fprintf(stream3,"\n");*/
      /*      printf("%30s Scale Ratio=%7.3f Nstars=%d \n X = %f + x*%f y*%f\n Y = %f + x*%f + y*%f\n",file1,scaleratio,nf,ax,bx,cx,ay,by,cy); */
      printf("%f %f %f %f %f %f\n",bx,cx,ax,by,cy,ay); 
    }
    else     {
      printf("\nMATCH STAR FAILED\n"); 
      exit(0);
    }
    fclose(stream1);
    fclose(stream2);
    fclose(stream3);
  }
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
