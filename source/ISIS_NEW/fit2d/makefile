 #CC    = cc 			# the C compiler
 #COPTS =  -fast            	# options for the C compiler
 LIBS  = -lm
 DEST  = ../bin

 ALL   =    main.o fitn.o lubksb.o ludcmp.o poly.o \
	   sigma.o


 fitn:  $(ALL)
	$(CC) $(ALL) -o fitn $(LIBS) $(COPTS)
	cp fitn $(DEST)

 main.o: $(STDH) main.c 
	$(CC) $(COPTS)  -c main.c

 fitn.o: $(STDH) fitn.c 
	$(CC) $(COPTS)  -c fitn.c

 lubksb.o: $(STDH) lubksb.c 
	$(CC) $(COPTS)  -c lubksb.c

 ludcmp.o: $(STDH)   ludcmp.c 
	$(CC) $(COPTS)  -c ludcmp.c

 poly.o: $(STDH)    poly.c 
	$(CC) $(COPTS)  -c  poly.c

 sigma.o: $(STDH)    sigma.c 
	$(CC) $(COPTS)  -c  sigma.c
	$(CC) $(COPTS)  -c  sigma.c
	$(CC) $(COPTS)  -c  sigma.c
	$(CC) $(COPTS)  -c  sigma.c
