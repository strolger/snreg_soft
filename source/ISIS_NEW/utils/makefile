ORIG  = .			# where the sources lie
DEST  = ..			# where to put the executable
#COPTS =  -fast            	# options for the C compiler
DEST  = ../bin

LIBS  =  -lm 


 ALL   =  main.o writefits.o readfits.o


 convert: $(ALL)
	$(CC) $(ALL) -o convert $(LIBS) $(COPTS)
	cp convert $(DEST)


 

 main.o: $(STDH) main.c
	$(CC) $(COPTS)  -c main.c

 writefits.o: $(STDH) writefits.c
	$(CC) $(COPTS)  -c writefits.c

 readfits.o: $(STDH) readfits.c
	$(CC) $(COPTS)  -c readfits.c









