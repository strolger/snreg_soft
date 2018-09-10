#include "short.h"

#define   CENT_TYPE double
int       width,height,bitpix;
int       NSTARS;
DATA_TYPE saturation,radphot,rad1bg,rad2bg;
double    ADU_EL;
int       *bg_index,*bg_list,nbg,bg_edge;
DATA_TYPE *bg_temp;
DATA_TYPE ratio_max;
int       Nfirst,Nkeep;
DATA_TYPE **stack;
int       stack_width;
DATA_TYPE *psf;
int       mesh2;
double    ax0,ay0;

typedef struct
{
  int       xi,yi;
  DATA_TYPE x,y;
  DATA_TYPE aperphot,aperphot2,bg,sigma;
  DATA_TYPE pixmax,ratio;
  DATA_TYPE *stack;
} object;




DATA_TYPE *readfits(char *name);
object    *max(DATA_TYPE *,DATA_TYPE,int,DATA_TYPE,DATA_TYPE);
void      cent(int ,int ,CENT_TYPE *,CENT_TYPE *,DATA_TYPE *);
void      init_bg(double,double);
double    background(int,int,DATA_TYPE *,double *);
void      aperphot(object *,DATA_TYPE *,DATA_TYPE,int *);
void      make_psf(object *,DATA_TYPE *,int);
void      norm_psf();
void      get_gaussian(double *,double *,double *);
void      make_phot(object *,DATA_TYPE *,double,double,double,int,double,char *);
void      read_config(char *);
void    splie2(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int ,int ,DATA_TYPE *);
void    splin2(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int,int ,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *);
void    spline(DATA_TYPE *,DATA_TYPE *,int,DATA_TYPE,DATA_TYPE,DATA_TYPE *,DATA_TYPE *);
void    splint(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int,DATA_TYPE ,DATA_TYPE *);
void    splin2(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int m,int n,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *);
DATA_TYPE get_thresh(DATA_TYPE *,DATA_TYPE *);
