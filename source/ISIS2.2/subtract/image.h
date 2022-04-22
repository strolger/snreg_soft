#include "short.h"

int       width,height;
char      *header;
int       nstamps,stamp_size,mesh_size,half_mesh_size,half_stamp_size;
DATA_TYPE SATURATION,PIX_MIN;
int       max_stamp_thresh;
int       nstamps_x,nstamps_y;
int       sub_width,sub_height;
int       width0,height0,bitpix;
DATA_TYPE *sub_image,*sub_ref;
DATA_TYPE *check_stack;
char      *def_map,*def_mask;
double    **check_mat,*check_vec;

typedef struct
{
  int       x,y;
  DATA_TYPE **vectors;
  DATA_TYPE *area;
  double    **mat;
  double    *scprod;
  double    sum;
  double    chi2;
  double    norm;
  double    diff;
  char      keep;
} stamp_struct;


stamp_struct *stamps;
int          stamp_number;
