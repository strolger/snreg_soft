/* Program to apply the inverse of a jtxform map to two coords */
/* John Tonry - 13 Oct 2000 */
/*
 * Syntax: revmap mapfile x y
 */

#include <stdio.h>
#include <math.h>

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))

#define MAXPAR 10

main(argc,argv)
int argc;
char **argv;
{
  double x, y, xdst, ydst, xsrc, ysrc, xmap, ymap;
  double xgx, xgy, ygx, ygy, xprv, yprv;
  char line[256], *success;
  FILE *input;
  double xpar[MAXPAR], ypar[MAXPAR], xm, xp, ym, yp, xrng, yrng;
  int i, secret=0;
  double atof();

/* Parse the arguments */
  if(argc < 4) {
    syntax(*argv);
    exit(1);
  }

  if(argc >= 5) secret = 1;

/* What are the desired destination coords for which we want the source coords? */
  xdst = atof(argv[2]);
  ydst = atof(argv[3]);

/* Get the map */
  if((input = fopen(argv[1],"r")) == NULL) {
    fprintf(stderr,"Cannot open map file %s\n", argv[1]);
    exit(1);
  }
  if( (success = fgets(line, 256, input)) == line)
    sscanf(line, "%lf %lf %lf %lf", &xm, &xp, &ym, &yp);
  for(i=0; i<MAXPAR; i++) {
    if( (success = fgets(line, 256, input)) == line)
      sscanf(line, "%lf %lf", &xpar[i], &ypar[i]);
  }
  fclose(input);
  xrng = (xp-xm) / 2;
  yrng = (yp-ym) / 2;

/* Show us the frontwards map just to make us feel safe */
  if(secret) {
    x = (xdst-0.5*(xp+xm)) / xrng;
    y = (ydst-0.5*(yp+ym)) / yrng;
    xmap =   xpar[0] +   xpar[1]*x + xpar[2]*y +
      xpar[3]*x*x +  xpar[4]*x*y + xpar[5]*y*y +
      xpar[6]*x*x*x +  xpar[7]*x*x*y + xpar[8]*x*y*y + xpar[9]*y*y*y;
    ymap =   ypar[0] +   ypar[1]*x + ypar[2]*y +
      ypar[3]*x*x +  ypar[4]*x*y + ypar[5]*y*y +
      ypar[6]*x*x*x +  ypar[7]*x*x*y + ypar[8]*x*y*y + ypar[9]*y*y*y;
    printf("Direct map: %10.2f %10.2f -> %10.2f %10.2f\n", xdst, ydst, xmap, ymap);
    exit(0);
  }

/* Calculate mapping of the destination image back to source coordinates */
  xsrc = 0.5*(xm+xp);
  ysrc = 0.5*(ym+yp);
  x = (xsrc-0.5*(xp+xm)) / xrng;
  y = (ysrc-0.5*(yp+ym)) / yrng;

  for(i=0; i<10; i++) {
    xprv = xsrc;
    yprv = ysrc;
    xmap =   xpar[0] +   xpar[1]*x + xpar[2]*y +
      xpar[3]*x*x +  xpar[4]*x*y + xpar[5]*y*y +
      xpar[6]*x*x*x +  xpar[7]*x*x*y + xpar[8]*x*y*y + xpar[9]*y*y*y;
    ymap =   ypar[0] +   ypar[1]*x + ypar[2]*y +
      ypar[3]*x*x +  ypar[4]*x*y + ypar[5]*y*y +
      ypar[6]*x*x*x +  ypar[7]*x*x*y + ypar[8]*x*y*y + ypar[9]*y*y*y;

    xgx = xpar[1] + 2*xpar[3]*x +  xpar[4]*y +
      3*xpar[6]*x*x +  2*xpar[7]*x*y + xpar[8]*y*y;
    xgy = xpar[2] + xpar[4]*x + 2*xpar[5]*y +
      xpar[7]*x*x + 2*xpar[8]*x*y + 3*xpar[9]*y*y;
    ygx = ypar[1] + 2*ypar[3]*x +  ypar[4]*y +
      3*ypar[6]*x*x +  2*ypar[7]*x*y + ypar[8]*y*y;
    ygy = ypar[2] + ypar[4]*x + 2*ypar[5]*y +
      ypar[7]*x*x + 2*ypar[8]*x*y + 3*ypar[9]*y*y;

    x += ((xdst-xmap)*xgx-(ydst-ymap)*ygx)/(xgx*xgx-ygx*ygx);
    y += ((xdst-xmap)*xgy-(ydst-ymap)*ygy)/(xgy*xgy-ygy*ygy);

    xsrc = x*xrng + 0.5*(xp+xm);
    ysrc = y*yrng + 0.5*(yp+ym);
    if(ABS(xprv-xsrc) < 0.0001 && ABS(yprv-ysrc) < 0.0001) break;
  }
  printf("%10.2f %10.2f\n", xsrc, ysrc);

}

syntax(s)
char *s;
{
  printf("Syntax: revmap mapfile x y\n");
}
