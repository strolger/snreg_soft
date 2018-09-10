#include<stdio.h>
#include<math.h>
#include<string.h>
#define HEADLINES 800
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))

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
  float MAG,MAGSPACE,SPACING;
  int x1,x2,y1,y2,flag,gx1,gx2,gy1,gy2,n2stars;
  void findstar();	

  if (argc < 4) {
    syntax(*argv);
    exit(0);
  }   

  sprintf(field,argv[1]);
  n2stars = atoi(argv[2]);
  MAG = atof(argv[3]);
  MAGSPACE=atof(argv[4]);
  SPACING=atof(argv[5]);


  
  strcpy(info,field);
  strcat(info,".info");


  stream1 = fopen(info,"w");   

  findstar(field,&gx1,&gy1,SPACING,n2stars);

  fprintf(stream1,"FIELD NAME: %s\n",field);
  fprintf(stream1,"SENSPOS:%d %d %f %d %f %f\n",gx1,gy1,SPACING,n2stars,MAG,MAGSPACE);     

  fclose(stream1);


  stream4=fopen("tosubtract.coo","w");

  fprintf(stream4,"#N ID    XCENTER   YCENTER   MAG         MERR          MSKY           NITER    \\\n");
  fprintf(stream4,"#U ##    pixels    pixels    magnitudes  magnitudes    counts         ##       \\\n");
  fprintf(stream4,"#F %%-9d  %%-10.3f   %%-10.3f   %%-12.3f     %%-14.3f       %%-15.7g        %%-6d\n");
  
  fprintf(stream4,"\n#N         SHARPNESS   CHI         PIER  PERROR                                \\\n");
  fprintf(stream4,"#U         ##          ##          ##    perrors                               \\\n");
  fprintf(stream4,"#F         %%-23.3f     %%-12.3f     %%-6d  %%-13s                                  \n"); 
      
  

  for (j=1;j<=n2stars;j++) {/*write out positions of sensitivity check stars*/
    for (k=1;k<=n2stars;k++) {
      
        fprintf(stream4,"%-9d%-10.3f%-10.3f%-12.3f%-14.3f%-15.7f%-6d   \\\n",1,(gx1)+SPACING*(k-1),(gy1)+SPACING*(j-1),MAG+MAGSPACE*((k-1)*n2stars+j-1),0.01,0.,3);
	fprintf(stream4,"%-23.3f%-12.3f%-6d%-13s\n",.01,1.,0,"No_error");
      }
  }
  fclose(stream4);
}
  
  
 

void findstar (field,gx1,gy1,SPACING,n2stars)
     char field[];
     float SPACING;
     int *gx1,*gy1,n2stars;		

 {
  
   FILE *stream1,*stream2,*stream3,*stream4,*stream5,*stream6;
   short *data;
   char searchfile[25],reffile[25],info[100],fieldname[100],subname[100];
   char gal[100],line[150],remotefile[50],name[50],num[5],batchname[50],stars[50];
   char *header,infile[100];
   int i,j,k,nstar1,nstar2,ii,itype[10000],itype2[10000],l,RN=0;
   int testSN(),nhead=HEADLINES,nhead1=HEADLINES,index[10000],flagobj[10000];
   float err[10000], elong,x1[10000],y1[10000],amag1[10000],x2[10000],y2[10000],amag2[10000],score[10000],sum[10000];
   float dummy,C,Cx,Cy,Cm,min=-1000,tmin,round[10000],sharp[10000];
   float FWHM1,FWHM2,FWHM1A,FWHM2A,FWHM1B,FWHM2B,FWHM,max[10000];
   float temp,POSERR;
   float summax,minmax,maxmin,THRESH,dist,seeing,background,bkgrnd(),testbox();
   int X1,X2,Y1,Y2,flag=0,xaxis,yaxis,galN,n3,remote=0,flaggal=0;
   int bx1[500],bx2[500],by1[500],by2[500],star[50],maxRN,x,y,xx,yy;
   int mask1[100],mask2[100],mask3[100],mask4[100],nmask,xbox[5000],ybox[5000];
   void itoa(),reverse();
   float *image;
   
   strcpy(infile,field);
   strcat(infile,".fits");

   rfitsreal(&header, &xaxis, &yaxis, &image, infile); 
   printf("Read %s, size = (%d, %d)...\n", field, xaxis, yaxis);

   /*finds background for image */
   background=bkgrnd(image,(float)xaxis/2,(float)yaxis/2,xaxis,yaxis,xaxis/2);
   
   /*divide the image into pieces and score them for paucity of activity*/
   for (k=0,i=xaxis/4;i<xaxis*3/4;i+=n2stars*SPACING) {
     for (j=yaxis/4;j<yaxis*3/4;j+=n2stars*SPACING) {
       k++;
       score[k]=testbox(image,(float)i,(float)j,xaxis,yaxis,background,n2stars,SPACING);
       xbox[k]=i;
       ybox[k]=j;
     }
   }
   indexx(k,score,index); /* sort to get best first*/   
 
   *gx1=xbox[index[1]]-(n2stars-1)/2*SPACING; /*return best x and y position*/
   *gy1=ybox[index[1]]-(n2stars-1)/2*SPACING;
   
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

float testbox(data,X,Y,nx,ny,back,n2,spacing)  
     int nx,ny,n2;
     float X,Y,back,spacing,*data;
{
  int j,jlo,jhi,klo,khi,index[50000],N,i,backlo,backhi;
  float medpix[50000],background,x,y,bkgrnd();
  float sum=0,score=0;
  
  x=X; 
  y=Y;		

  background=bkgrnd(data,x,y,nx,ny,(int)spacing*n2); /*get background of this particular region*/
  if (background<=0) return(1e5);
  score = fabs((back-background)/back)*300; /*add to score if it not close*/
  if (fabs((back-background)/back) > .025) /*as you go away from the mean background, make it goe as the square of the score*/
    score *= score;


/* load in the pixels in this region*/
  jlo = y-spacing*(float)n2/2;
  jhi = y+spacing*(float)n2/2;		
  klo = x-spacing*(float)n2/2;;
  khi = x+spacing*(float)n2/2;;
    
  if (jlo < 0) jlo =0;
  if (klo < 0) klo = 0;
  if (jhi >= ny) jhi = ny;
  if (khi >= nx) khi = nx;
  
    
  for (N=1,i=klo; i<=khi ;i++) { /* go over x*/
    for (j=jlo; j<=jhi;j++) { /* go over y*/
      medpix[N] = data[j*nx+i]-background; /*what is this minus background*/
      N++;
    }
  }  
  N--;
  indexx(N,medpix,index); /* sort to get median,etc.*/
  for (sum=0,i=1;i<=N;i++) {
    sum+=fabs((float)medpix[i]);
  }
   
  score += (sum/N)/sqrt(background)*10; /*score goes as average pixel/some estimate of the noise*/
  score += medpix[3]/-100; /*stay clear if > 3 bad pixels low*/
  score += medpix[N-2]/100; /*stay clear if > 3 pixel bad high*/

  return(score);
}



float bkgrnd(data,x,y,nx,ny,s)
     int nx,ny,s;  
     float x,y;
     float *data;
{
  int j,k,jlo,jhi,klo,khi,index[5000],N,i;
  float background,medpix[5000];
  
  jlo = y-s;
  jhi = y+s;		
  klo = x-s;
  khi = x+s;
  
  if (jlo < 0) jlo =0;
  if (klo < 0) klo = 0;
  if (jhi >= ny) jhi = ny;
  if (khi >= nx) khi = nx;
  
  
  for (N=1,i=klo; i<=khi ;i+=(s/10+1)) { /* go over x*/
    for (j=jlo; j<=jhi;j+=(s/10+1)) { /* go over y*/
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
 
 /*scores a box for paucity of activity..*/             

  

syntax(s)
char *s;
{
  printf("Syntax: premkt \n");
  printf("       field filename  .fits will be added to this to open image\n");
  printf("       n2stars 	 NxN grid of stars to subtract n"); 
  printf("       MAG             Magnitude of brightest star to subtract\n");
  printf("       MAGSPACE        Magnitude Spacing between stars\n");
  printf("       SPACING         Spacing in pixels between subtracted stars\n");
}



