ORIG  = .			# where the sources lie
DEST  = ..			# where to put the executable
#CC    = cc 			# the C compiler
#COPTS =  -fast            	# options for the C compiler
DEST  = ../bin
STDH  = types.h
LIBS  =  -lm 


 ALL   =  cent.o main.o max.o readfits.o quick_sort.o background.o \
          aperphot.o make_psf.o norm_psf.o get_gaussian.o make_phot.o\
	  read_config.o spline.o splint2.o splie2.o spline2.o get_thresh.o


 phot_ref: $(ALL)
	$(CC) $(ALL) -o phot_ref $(LIBS) $(COPTS)
	cp phot_ref $(DEST)


 cent.o: $(STDH)  cent.c
	$(CC) $(COPTS)  -c  cent.c

 main.o: $(STDH) main.c
	$(CC) $(COPTS)  -c main.c

 max.o: $(STDH) max.c
	$(CC) $(COPTS)  -c max.c

 readfits.o: $(STDH) readfits.c
	$(CC) $(COPTS)  -c readfits.c

 quick_sort.o: $(STDH) quick_sort.c
	$(CC) $(COPTS)  -c quick_sort.c

 background.o: $(STDH) background.c
	$(CC) $(COPTS)  -c background.c

 aperphot.o: $(STDH) aperphot.c
	$(CC) $(COPTS)  -c aperphot.c

 make_psf.o: $(STDH) make_psf.c
	$(CC) $(COPTS)  -c make_psf.c

 norm_psf.o: $(STDH) norm_psf.c
	$(CC) $(COPTS)  -c norm_psf.c

 get_gaussian.o: $(STDH) get_gaussian.c
	$(CC) $(COPTS)  -c get_gaussian.c

 make_phot.o: $(STDH) make_phot.c
	$(CC) $(COPTS)  -c make_phot.c

 read_config.o: $(STDH) read_config.c
	$(CC) $(COPTS)  -c read_config.c

 spline.o: $(STDH) spline.c
	$(CC) $(COPTS)  -c spline.c

 splint2.o: $(STDH) splint2.c
	$(CC) $(COPTS)  -c splint2.c

 splie2.o: $(STDH) splie2.c
	$(CC) $(COPTS)  -c splie2.c
 
 spline2.o: $(STDH) spline2.c
	$(CC) $(COPTS)  -c spline2.c

 get_thresh.o: $(STDH) get_thresh.c
	$(CC) $(COPTS)  -c get_thresh.c

 #.o: $(STDH) .c
 #	$(CC) $(COPTS)  -c .c
 
 
 
 
