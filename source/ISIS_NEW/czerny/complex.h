#include<math.h>



struct comp{
 double real;
 double imaginary;
};

typedef struct comp complex;

struct comp_mat{
 complex *mat;
 int      nb;
};

typedef struct comp_mat complex_mat;
