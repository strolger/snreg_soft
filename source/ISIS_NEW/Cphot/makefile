ORIG  = .			# where the sources lie
DEST  = ..			# where to put the executable
#COPTS =  -fast           	# options for the C compiler
DEST  = ../bin
STDH  = types.h
LIBS  =  -lm 


 ALL   =  main.o writefits.o readfits.o make_phot.o make_vectors.o \
          kernel_vector.o read_config.o power.o get_coeffs.o make_kernel.o \
	  get_psf.o read_phot_conf.o convolve.o spline_allocate.o Bspline.o \
	  splie2.o spline.o splint.o spline2.o phot.o


 Cphot: $(ALL)
	$(CC) $(ALL) -o Cphot $(LIBS) $(COPTS)
	cp Cphot $(DEST)

 

 main.o: $(STDH) main.c
	$(CC) $(COPTS)  -c main.c

 writefits.o: $(STDH) writefits.c
	$(CC) $(COPTS)  -c writefits.c

 readfits.o: $(STDH) readfits.c
	$(CC) $(COPTS)  -c readfits.c

 make_phot.o: $(STDH) make_phot.c
	$(CC) $(COPTS)  -c make_phot.c

 make_vectors.o:  $(STDH) make_vectors.c
	$(CC) $(COPTS) -c make_vectors.c

 kernel_vector.o:  $(STDH) kernel_vector.c
	$(CC) $(COPTS) -c kernel_vector.c

 read_config.o:  $(STDH) read_config.c
	$(CC) $(COPTS) -c read_config.c

 power.o:  $(STDH) power.c
	$(CC) $(COPTS) -c power.c

 get_coeffs.o:  $(STDH) get_coeffs.c
	$(CC) $(COPTS) -c get_coeffs.c

 make_kernel.o:  $(STDH) make_kernel.c
	$(CC) $(COPTS) -c make_kernel.c

 get_psf.o:  $(STDH)  get_psf.c
	$(CC) $(COPTS) -c  get_psf.c

 read_phot_conf.o:  $(STDH) read_phot_conf.c
	$(CC) $(COPTS) -c read_phot_conf.c

 convolve.o:  $(STDH) convolve.c
	$(CC) $(COPTS) -c convolve.c

 spline_allocate.o:  $(STDH) spline_allocate.c
	$(CC) $(COPTS) -c spline_allocate.c

 Bspline.o:  $(STDH) Bspline.c
	$(CC) $(COPTS) -c Bspline.c

 splie2.o: $(STDH) splie2.c
	$(CC) $(COPTS)  -c splie2.c

 spline2.o: $(STDH) spline2.c
	$(CC) $(COPTS)  -c spline2.c

 spline.o: $(STDH) spline.c
	$(CC) $(COPTS)  -c spline.c

 splint.o: $(STDH) splint.c
	$(CC) $(COPTS)  -c splint.c

 phot.o:  $(STDH) phot.c
	$(CC) $(COPTS) -c phot.c


# .o:  $(STDH) .c
#	$(CC) $(COPTS) -c .c
#	$(CC) $(COPTS) -c .c
#	$(CC) $(COPTS) -c .c
#	$(CC) $(COPTS) -c .c
