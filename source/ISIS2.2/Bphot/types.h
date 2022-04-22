#include "short.h"

#define   CENT_TYPE double

int       width,height,width0,height0,bitpix,nload_x,nload_y;
double    ADU_EL;
int       *bg_index,*bg_list,nbg,bg_edge,nstars_buffer;
DATA_TYPE *bg_temp,**stack,*psf;
DATA_TYPE saturation,radphot,rad1bg,rad2bg,pix_min,rad_aper;
int       stack_width,NSTARS,mesh2,nnstars,NSTARS_MAX;

typedef struct
{
  int       xi,yi;
  DATA_TYPE x,y;
  DATA_TYPE aperphot,aperphot2,bg,sigma;
  DATA_TYPE pixmax,ratio,phot;
  DATA_TYPE *stack;
} object;



void      read_image(FILE *,DATA_TYPE *,int,int,int,int,double,double);
void      read_config(char *);
FILE      *read_header(char *,int *,double *,double *);
object    *max(DATA_TYPE *,DATA_TYPE,int,DATA_TYPE,DATA_TYPE);
void      cent(int ,int ,CENT_TYPE *,CENT_TYPE *,DATA_TYPE *);
void      init_bg(double,double);
double    background(int,int,DATA_TYPE *,double *);
void      aperphot(object *,DATA_TYPE *);
void      make_psf(object *,DATA_TYPE *,int,int,int,int,FILE *,char *);
void      splie2(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int ,int ,DATA_TYPE *);
void      splin2(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int,int ,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *);
void      spline(DATA_TYPE *,DATA_TYPE *,int,DATA_TYPE,DATA_TYPE,DATA_TYPE *,DATA_TYPE *);
void      splint(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int,DATA_TYPE ,DATA_TYPE *);
void      writefits(DATA_TYPE *, int, int, char *);
