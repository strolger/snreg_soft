#define   DATA_TYPE float
#define   CENT_TYPE double
int       width,height,bitpix;
DATA_TYPE *readfits(char *name);
void      writefits(DATA_TYPE *,char *name);
