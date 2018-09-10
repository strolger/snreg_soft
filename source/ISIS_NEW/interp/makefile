ORIG  = .			# where the sources lie
DEST  = ../bin			# where to put the executable
#CC    = cc 			# the C compiler
#COPTS =  -fast            	# options for the C compiler
LIBS  =  -lm 
STDH  = types.h short.h

ALL   =  main.o poly.o spline.o splint2.o splie2.o spline2.o readfits.o \
         writefits.o


 interp: $(ALL)
	$(CC) $(ALL) -o interp $(LIBS) $(COPTS)
	cp interp $(DEST)


 main.o: $(STDH) main.c
	$(CC) $(COPTS)  -c main.c

 poly.o: $(STDH) poly.c
	$(CC) $(COPTS)  -c poly.c
 
 spline.o: $(STDH) spline.c
	$(CC) $(COPTS)  -c spline.c

 splint2.o: $(STDH) splint2.c
	$(CC) $(COPTS)  -c splint2.c

 splie2.o: $(STDH) splie2.c
	$(CC) $(COPTS)  -c splie2.c
 
 spline2.o: $(STDH) spline2.c
	$(CC) $(COPTS)  -c spline2.c

 readfits.o: $(STDH) readfits.c
	$(CC) $(COPTS)  -c readfits.c

 writefits.o: $(STDH) writefits.c
	$(CC) $(COPTS)  -c writefits.c
	$(CC) $(COPTS)  -c writefits.c
	$(CC) $(COPTS)  -c writefits.c
	$(CC) $(COPTS)  -c writefits.c
