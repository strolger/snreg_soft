#include "short.h"


#define   CENT_TYPE double
int       width,height,bitpix;
DATA_TYPE saturation,adu_el,*temp;
DATA_TYPE *readfits(char *name);
void      max(DATA_TYPE *,DATA_TYPE *,DATA_TYPE,int,char *);
void      cent(int ,int ,CENT_TYPE *,CENT_TYPE *,DATA_TYPE *);
DATA_TYPE get_thresh(DATA_TYPE *,DATA_TYPE *);
void      kill_cosmics(DATA_TYPE *,DATA_TYPE);
void      smooth(DATA_TYPE *);
