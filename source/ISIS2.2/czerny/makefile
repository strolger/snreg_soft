ORIG  = .			# where the sources lie
DEST  = ../bin			# where to put the executable
#COPTS =  -fast          	# options for the C compiler
DEST  = ../bin
#CC    = gcc
LIBS  =  -lm 


 ALL   =  main.o series.o utilities.o fit.o


 czerny: $(ALL)
	$(CC) $(ALL) -o czerny $(LIBS) $(COPTS)
	cp czerny $(DEST)



 main.o: $(STDH) main.c
	$(CC) $(COPTS)  -c main.c

 fit.o: $(STDH)  fit.c
	$(CC) $(COPTS)  -c fit.c

 series.o: $(STDH) series.c
	$(CC) $(COPTS)  -c series.c

 utilities.o: $(STDH) utilities.c
	$(CC) $(COPTS)  -c utilities.c




