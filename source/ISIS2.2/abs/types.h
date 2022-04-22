#include<stdio.h>
#include "short.h"

#define   CENT_TYPE  double

int        width,height,bitpix,old_buffer,deg_bg,j0b,nlines_read,
           n_reject;
int       *ndex,*indx;
DATA_TYPE *buffer,saturation,pixmin,*temp_smooth;

double    *lsvec,**mat,*vec;


typedef struct
{
 double norm;
 double *coeff;
} bg_struct;


FILE      *read_header(char *);
void      readfits(FILE *,int,DATA_TYPE *,int);
void      writefits(DATA_TYPE *,char *);
void      cent(int,int,CENT_TYPE *,CENT_TYPE *,DATA_TYPE *);
void      clean_stack(DATA_TYPE **,DATA_TYPE *, int, int, char);
void      moy(DATA_TYPE **,DATA_TYPE *, int, int,char);
void      quick_sort(DATA_TYPE *,int *,int);
bg_struct fit(DATA_TYPE *,DATA_TYPE *,int);
void      fit_bg_norm(DATA_TYPE **,int,int, bg_struct *);
void      allocate();
void      ludcmp(double  **,int ,int *,double  *);
void      lubksb(double  **,int,int *,double  *);
void      read_config(char *);
double    get_background(int ,int ,bg_struct);
void      smooth(DATA_TYPE *, int, int);
