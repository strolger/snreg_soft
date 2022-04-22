#include  "short.h"
#define    CENT_TYPE double

int        bitpix;
int        ncomp_kernel,ngauss,deg_fixe[16];
DATA_TYPE  **kernel_vec,sigma_gauss[16],saturation,*filter_x,*filter_y,
           *sum_vectors,sum_kernel,rad_aper,radphot;
int        nstamps_x,nstamps_y,half_mesh_size,ngauss,deg_spatial,deg_bg,
           stack_width,sub_width,sub_height,mesh_size,ncomp_total,conv_size,
           image_width,adu_el;
double     *kernel_buff,exp_date;
DATA_TYPE  *x1a,*x2a,*ytmp,*y2a,*xxp,*yyp,*yy,*yytmp,*uv;



typedef struct
{
 char filename[256];
 int  xmin;
 int  ymin;
 int  xmax;
 int  ymax;
 int  offset_x;
 int  offset_y;
} table_struct;



DATA_TYPE   *readfits(DATA_TYPE *,char *name,int *,int *,int);
void        writefits(DATA_TYPE * ,int ,int ,char *);
void        make_phot(DATA_TYPE *,DATA_TYPE *,char *,char *);
DATA_TYPE   *kernel_vector(int,int,int,int,char *);
void        make_vectors();
double      power(double, int);
void        get_coeffs(double *coeffs,table_struct *,int,int,int);
void        make_kernel(DATA_TYPE *,int,int,double *);
int         get_psf(DATA_TYPE *,table_struct *,int ,int,int);
void        read_config(char *);
void        read_config_phot(char *);
void        convolve(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *);
void        spline_allocate();
void        Bspline(DATA_TYPE *,DATA_TYPE *,double, double);
void        splie2(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int ,int ,DATA_TYPE *);
void        splin2(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int,int ,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *);
void        spline(DATA_TYPE *,DATA_TYPE *,int,DATA_TYPE,DATA_TYPE,DATA_TYPE *,DATA_TYPE *);
void        splint(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int,DATA_TYPE ,DATA_TYPE *);
void        phot(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int,int,char *fname);
