#include<stdio.h>
#include<math.h>
#define SPACING 15.
#define MAGSPACE 0.375
#define MAG 17.
#define HEADLINES 800
main (argc,argv)
     int argc;
     char *argv[];
{  
  FILE *stream1,*stream2,*stream3,*stream4;
  char input[100],field[100],info[100],REF[100],CONV[100],fieldname[100];
  char section[100],line[150],word1[100],word2[100],stars[100],line1[150];
  int i,j,k,nstar1,nstar2,ii,itype;
  float dummy,FWHM1,FWHM2,FWHM1A,FWHM2A,FWHM1B,FWHM2B,C1x,C1y;
  float S1x,S1y,S1z,S2x,S2y,S2z,S3x,S3y,S3z,x,y,chi,width;
  int x1,x2,y1,y2,flag,gx1,gx2,gy1,gy2;
  void findstar();	

  sprintf(field,argv[1]);

  strcpy(info,field);
  strcpy(stars,field);
  strcat(stars,".stars");
  strcat(info,".info");

  stream1 = fopen(info,"w");   
  
  stream2 = fopen("obj.edited","r");

  while( feof(stream2) == 0 && itype !=1 ) {
    fscanf(stream2,"%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f"
	   ,&ii,&itype,&dummy,&dummy,&dummy,&dummy,&dummy,&FWHM1A,
	   &FWHM2A,&dummy,&dummy,&dummy,&dummy,&dummy,&dummy);
  }
  fclose(stream2);

  stream2 = fopen("obj.edited","r");
  stream3 = fopen(stars,"w");
  while( feof(stream2) == 0 ) {
    fgets(line,150, stream2); 
    sscanf(line,"%d %d %f %f %f %f",&ii,&itype,&dummy,&dummy,&dummy,&width);
    if (itype==1 || itype == 3 || (itype==7 && dummy < -9.0) || (itype==6 && dummy < -9.0) || (itype==2 && width < 2.) || (itype==9 && dummy < -9)) {
      fprintf(stream3,"%s",line);
    }
  }

  fclose(stream2);
  fclose(stream3);

  findstar(field,&S1x,&S1y,&S1z,&S2x,&S2y,&S2z,&S3x,&S3y,&S3z,&gx1,&gx2,&gy1,&gy2);

  fprintf(stream1,"FIELD NAME: %s\n",field);
  fprintf(stream1,"FWHM: %f %f\n",FWHM1A,FWHM2A);
  fprintf(stream1,"STAR1: %f %f %f\n",S1x,S1y,S1z);     
  fprintf(stream1,"STAR2: %f %f %f\n",S2x,S2y,S2z);     
  fprintf(stream1,"STAR3: %f %f %f\n",S3x,S3y,S3z);     
  fprintf(stream1,"SECTION:[%d:%d,%d:%d]\n",(int)S1x-10,(int)S1x+10,(int)S1y-10,(int)S1y+10);     
  fprintf(stream1,"CONSTANT:%f %f %f \n",0.,0.,0.);     
  fprintf(stream1,"SENSPOS:%f %f %f \n",(gx1+gx2)/2.,(gy1+gy2)/2.,SPACING);     

  fclose(stream1);


  stream4=fopen("edited.coo.1","w");
  fprintf(stream4,"%f %f\n%f %f\n%f %f\n",S1x,S1y,S2x,S2y,S3x,S3y);
  fclose(stream4);

  stream4=fopen("edited.pst.1","w");
  fprintf(stream4,"%f %f\n%f %f\n%f %f\n",S1x,S1y,S2x,S2y,S3x,S3y);
  fclose(stream4);
 
  stream4=fopen("edited.pk.1","w");

fprintf(stream4,"#N ID    XCENTER   YCENTER   MAG         MERR          MSKY           NITER    \\\n");
fprintf(stream4,"#U ##    pixels    pixels    magnitudes  magnitudes    counts         ##       \\\n");
fprintf(stream4,"#F %%-9d  %%-10.3f   %%-10.3f   %%-12.3f     %%-14.3f       %%-15.7g        %%-6d\n");

fprintf(stream4,"\n#N         SHARPNESS   CHI         PIER  PERROR                                \\\n");
fprintf(stream4,"#U         ##          ##          ##    perrors                               \\\n");
fprintf(stream4,"#F         %%-23.3f     %%-12.3f     %%-6d  %%-13s                                  \n"); 
      
  for (j=1;j<4;j++) {/*write out positions of sensitivity check stars*/
    for (k=1;k<4;k++) {
      fprintf(stream4,"%-9d%-10.3f%-10.3f%-12.3f%-14.3f%-15.7f%-6d   \\\n",1,(gx1+gx2)/2.+SPACING*(k-2),(gy1+gy2)/2.+SPACING*(j-2),MAG+((j-1)*3+k-1)*MAGSPACE,0.01,0.,3);
      fprintf(stream4,"%-23.3f%-12.3f%-6d%-13s\n",.01,1.,0,"No_error");
    }
  }
fclose(stream4);
}
  
  
 

void findstar (field,S1x,S1y,S1z,S2x,S2y,S2z,S3x,S3y,S3z,gx1,gx2,gy1,gy2)
     char field[];
     float *S1x,*S1y,*S1z,*S2x,*S2y,*S2z,*S3x,*S3y,*S3z;
     int *gx1,*gx2,*gy1,*gy2;		

 {
  
   FILE *stream1,*stream2,*stream3,*stream4,*stream5,*stream6;
   short *svector(),*image;
   char searchfile[25],reffile[25],info[100],fieldname[100],subname[100];
   char gal[100],line[150],remotefile[50],name[50],num[5],batchname[50],stars[50];
   char header[HEADLINES][80],header1[HEADLINES][80];
   int i,j,k,nstar1,nstar2,ii,itype[10000],itype2[10000],l,RN=0;
   int testSN(),nhead=HEADLINES,nhead1=HEADLINES,index[10000],flagobj[10000];
   float err[10000], elong,x1[10000],y1[10000],amag1[10000],x2[10000],y2[10000],amag2[10000],score[10000],sum[10000];
   float dummy,C,Cx,Cy,Cm,min=-1000,tmin,round[10000],sharp[10000];
   float FWHM1,FWHM2,FWHM1A,FWHM2A,FWHM1B,FWHM2B,FWHM,max[10000];
   float temp,POSERR;
   float summax,minmax,maxmin,THRESH,dist,seeing,background,bkgrnd();
   int X1,X2,Y1,Y2,flag=0,xaxis,yaxis,galN,n3,remote=0,flaggal=0;
   int bx1[500],bx2[500],by1[500],by2[500],star[50],maxRN,x,y,xx,yy;
   int mask1[100],mask2[100],mask3[100],mask4[100],nmask,xbox[5000],ybox[5000];
   void itoa(),reverse();
   float teststar(),testbox();
   image=svector(0,5000*5000);
   

   rffile(&nhead,header,&xaxis,&yaxis,image,"edited.fits",HEADLINES,40); /*read in images*/
   
   strcpy(info,field);
   strcat(info,".info");
  
   strcpy(stars,field);
   strcat(stars,".stars");
  
   
   stream3 = fopen(stars,"r");
   
   
   for (i=1; feof(stream3) == 0 && flag != 1; ) { /*Read in objects from obj.image_stage1*/
     fgets(line,150, stream3); 
     j=sscanf(line,"%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f",&ii,&itype[i],&x1[i],&y1[i],&amag1[i],&err[i],&dummy,&FWHM,&dummy,&dummy,&elong,&dummy,&dummy,&dummy,&dummy);
     if (itype[i] == 1 || itype[i] == 11 && j > 3) {
       if (seeing==0)  seeing=FWHM;
       i++;
     }
   } 
   nstar1 = i-1;
   fclose(stream3);

   background=bkgrnd(image,(float)xaxis/2,(float)yaxis/2,xaxis,yaxis,xaxis/2);
   for (i=1;i<=nstar1;i++) {
     score[i]=teststar(image,x1[i],y1[i],xaxis,yaxis,seeing,&max[i],background);
     score[i]+=(err[i]/0.01);
   }
   indexx(nstar1,score,index); /* sort to get best first*/   
   *S1x=x1[index[1]];
   *S1y=y1[index[1]];
   *S1z=max[index[1]];
   *S2x=x1[index[2]];
   *S2y=y1[index[2]];
   *S2z=max[index[2]];
   *S3x=x1[index[3]];
   *S3y=y1[index[3]];
   *S3z=max[index[3]];  
   
   for (k=0,i=xaxis/4;i<xaxis*3/4;i+=30) {
     for (j=yaxis/4;j<yaxis*3/4;j+=30) {
	k++;
	score[k]=testbox(image,i,j,xaxis,yaxis,background);
	xbox[k]=i;
	ybox[k]=j;
     }
   }
   indexx(k,score,index); /* sort to get best first*/   
 
   *gx1=xbox[index[1]]-30;
   *gx2=xbox[index[1]]+30;
   *gy1=ybox[index[1]]-30;
   *gy2=ybox[index[1]]+30;
   
   

   
 }
   

void itoa(n,s)
     char s[];
     int n;
{
  int i,o, sign;
  void reverse();
  
  if ((sign = n) < 0)
    n = -n;
  i = 0;
  do {
    s[i++] = n % 10 + '0';
  } while ((n /= 10) > 0);
  if (sign < 0)
    s[i++] = '-';
  s[i] = '\0';
  reverse(s);
}



void reverse(s)
     char s[];
{
  int c,i,j;
  
  for (i=0, j = strlen(s) -1; i < j; i++, j--) {
    c = s[i];
    s[i] =s [j];
    s[j] =c;
  }
}






float teststar(data,x,y,nx,ny,FWHM,max,back)   /*Check for bad pixels, etc..*/
     int nx,ny;
     float FWHM,x,y,*max,back;
     short *data;
{
  int j,jlo,jhi,klo,khi,index[42*42],N,i,backlo,backhi;
  float medpix[42*42],background;
  float sum=0,score=0,bkgrnd();

  
  background=bkgrnd(data,x,y,nx,ny,30);
  score = fabs((back-background)/back)*300;
  if (fabs((back-background)/back) > .025) 
    score *= score;

  jlo = y-2;
  jhi = y+2;		
  klo = x-2;
  khi = x+2;
    
  if (jlo < 0) jlo =0;
  if (klo < 0) klo = 0;
  if (jhi >= ny) jhi = ny;
  if (khi >= nx) khi = nx;
  
    
  for (N=1,i=klo; i<=khi ;i++) { /* go over x*/
    for (j=jlo; j<=jhi;j++) { /* go over y*/
      medpix[N] = data[j*nx+i]-background;
      N++;
    }
  }  
  N--;
  indexx(N,medpix,index); /* sort to get median,etc.*/
  *max=medpix[N];
  score += 100*(((1.8/FWHM)*(1.8/FWHM)*20000)-medpix[N])/((1.8/FWHM)*(1.8/FWHM)*20000)*(((1.8/FWHM)*(1.8/FWHM)*20000)-medpix[N])/((1.8/FWHM)*(1.8/FWHM)*20000);
  if ((medpix[N] > 22000) || ((((1.8/FWHM)*(1.8/FWHM)*20000)-medpix[N])/((1.8/FWHM)*(1.8/FWHM)*20000)) < -.15)
    score+=100;
       
  jlo = y-20;
  jhi = y+20;		
  klo = x-20;
  khi = x+20;
  
  if (jlo < 0) jlo =0;
  if (klo < 0) klo = 0;
  if (jhi >= ny) jhi = ny;
  if (khi >= nx) khi = nx;
       
      
  for (N=1,i=klo; i<=khi ;i++) { /* go over x*/
    for (j=jlo; j<=jhi;j++) { /* go over y*/
      if (sqrt((float)(i-x)*(i-x)+(float)(j-y)*(j-y)) >1.5*FWHM) {
	medpix[N] = data[j*nx+i]-background;
	sum+=medpix[N];
	N++;
      }
    }
  }  
  N--;
  if (N > 5) {
    indexx(N,medpix,index); /* sort to get median,etc.*/
  }    
  score += (sum/N/2*sum/N/2);
  score += (medpix[N]/75)*(medpix[N]/75);
  score += (x-nx/2)/(nx/6)*(x-nx/2)/(nx/6);
  score += (y-ny/2)/(ny/6)*(y-ny/2)/(ny/6);
  if (x < 200 || y < 200 || nx-x < 200 || ny-y < 200)
    score+=99;
  return(score);
}
  

float bkgrnd(data,x,y,nx,ny,s)
     int nx,ny,s;  
     float x,y;
     short *data;
{
  int j,k,jlo,jhi,klo,khi,index[500],N,i;
  float background,medpix[500];
  
  jlo = y-s;
  jhi = y+s;		
  klo = x-s;
  khi = x+s;
  
  if (jlo < 0) jlo =0;
  if (klo < 0) klo = 0;
  if (jhi >= ny) jhi = ny;
  if (khi >= nx) khi = nx;
  
  
  for (N=1,i=klo; i<=khi ;i+=s/10) { /* go over x*/
    for (j=jlo; j<=jhi;j+=s/10) { /* go over y*/
      if (data[j*nx+i] > -100) { /* only use good pixels */
	medpix[N] = data[j*nx+i];
	N++;
      }
    }
  }  
  N--;
  if (N > 10) {
    indexx(N,medpix,index); /* sort to get median,etc.*/
    background = medpix[N/2];
    return(background);
  }
  else return(-9999);
}
 
  
float testbox(data,X,Y,nx,ny,back)   /*Check for bad pixels, etc..*/
     int nx,ny,X,Y;
     float back;
     short *data;
{
  int j,jlo,jhi,klo,khi,index[60*60],N,i,backlo,backhi;
  float medpix[60*60],background,x,y;
  float sum=0,score=0;
  
  x=X; 
  y=Y;		

  background=bkgrnd(data,x,y,nx,ny,60);
  score = fabs((back-background)/back)*300;
  if (fabs((back-background)/back) > .025) 
    score *= score;

  jlo = y-29;
  jhi = y+29;		
  klo = x-29;
  khi = x+29;
    
  if (jlo < 0) jlo =0;
  if (klo < 0) klo = 0;
  if (jhi >= ny) jhi = ny;
  if (khi >= nx) khi = nx;
  
    
  for (N=1,i=klo; i<=khi ;i++) { /* go over x*/
    for (j=jlo; j<=jhi;j++) { /* go over y*/
      medpix[N] = data[j*nx+i]-background;
      N++;
    }
  }  
  N--;
  indexx(N,medpix,index); /* sort to get median,etc.*/
  for (sum=0,i=1;i<=N;i++) {
    sum+=fabs((float)medpix[i]);
  }
  score += (sum/N)/sqrt(background)*10;
  score += medpix[3]/-100;
  score += medpix[N-2]/100;

  return(score);
}
  
  
  


