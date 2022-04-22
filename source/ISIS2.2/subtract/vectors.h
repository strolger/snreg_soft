int       deg_fixe[16],ngauss,ncomposantes,deg_bg,deg_spatial,
          ncomp_kernel,ncomp_kernel_xy,ncomp_bg,ncomp_total,*indx,
          ncomp_background;
DATA_TYPE *filter_x,*filter_y,sigma_gauss[16],*temp,**kernel_vec;
double    **wxy,*scprod,*kernel_coeffs,*kernel,sum_kernel;
DATA_TYPE *conv_image,*sum_vectors;
