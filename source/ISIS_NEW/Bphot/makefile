ORIG  = .                       # where the sources lie
DEST  = ..                      # where to put the executable
#CC    = cc                     # the C compiler
##COPTS =  -fast                    # options for the C compiler
DEST  = ../bin
STDH  = types.h
LIBS  =  -lm 

 ALL   =  main.o read_image.o read_header.o read_config.o cent.o \
	  background.o make_psf.o splie2.o spline.o splint.o \
	  spline2.o quick_sort.o max.o aperphot.o writefits.o

 Bphot: $(ALL)
	$(CC) $(ALL) -o Bphot $(LIBS) $(COPTS)
	cp Bphot $(DEST)

 main.o: $(STDH) main.c
	$(CC) $(COPTS)  -c main.c

 read_image.o: $(STDH) read_image.c
	$(CC) $(COPTS)  -c read_image.c

 read_header.o: $(STDH) read_header.c
	$(CC) $(COPTS)  -c read_header.c

 read_config.o: $(STDH) read_config.c
	$(CC) $(COPTS)  -c read_config.c

 cent.o: $(STDH) cent.c
	$(CC) $(COPTS)  -c cent.c

 background.o: $(STDH) background.c
	$(CC) $(COPTS)  -c background.c

 make_psf.o: $(STDH) make_psf.c
	$(CC) $(COPTS)  -c make_psf.c

 splie2.o: $(STDH) splie2.c
	$(CC) $(COPTS)  -c splie2.c

 spline2.o: $(STDH) spline2.c
	$(CC) $(COPTS)  -c spline2.c

 spline.o: $(STDH) spline.c
	$(CC) $(COPTS)  -c spline.c

 splint.o: $(STDH) splint.c
	$(CC) $(COPTS)  -c splint.c

 quick_sort.o: $(STDH) quick_sort.c
	$(CC) $(COPTS)  -c quick_sort.c

 max.o: $(STDH) max.c
	$(CC) $(COPTS)  -c max.c

 aperphot.o: $(STDH) aperphot.c
	$(CC) $(COPTS)  -c aperphot.c

 writefits.o: $(STDH) writefits.c
	$(CC) $(COPTS)  -c writefits.c

#.o: $(STDH) .c
#	$(CC) $(COPTS)  -c .c
#	$(CC) $(COPTS)  -c .c
#	$(CC) $(COPTS)  -c .c
#	$(CC) $(COPTS)  -c .c
#	$(CC) $(COPTS)  -c .c
