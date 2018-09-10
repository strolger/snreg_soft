#include "short.h"
int       width,height,bitpix;
DATA_TYPE *readfits(char *name);
extern double  poly(double , double ,int ,double *);
extern void    splie2(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int ,int ,DATA_TYPE *);
extern void    splin2(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int,int ,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *);
extern void    writefits(DATA_TYPE *, char *);
extern void    spline(DATA_TYPE *,DATA_TYPE *,int,DATA_TYPE,DATA_TYPE,DATA_TYPE *,DATA_TYPE *);
extern void    splint(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int,DATA_TYPE ,DATA_TYPE *);
extern void    splin2(DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,int m,int n,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *,DATA_TYPE *);
