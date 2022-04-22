ORIG  = .			# where the sources lie
DEST  = ..			# where to put the executable
#CC    = cc 			# the C compiler
#COPTS =  -fast           	# options for the C compiler
DEST  = ../bin

LIBS  =  -lm 


 ALL   =  cent.o  main.o max.o readfits.o read_config.o quick_sort.o \
          get_thresh.o kill_cosmics.o writefits.o

 extract: $(ALL)
	$(CC) $(ALL) -o extract $(LIBS) $(COPTS)
	cp extract $(DEST)


 cent.o: $(STDH)  cent.c
	$(CC) $(COPTS)  -c  cent.c

 main.o: $(STDH) main.c
	$(CC) $(COPTS)  -c main.c

 max.o: $(STDH) max.c
	$(CC) $(COPTS)  -c max.c

 readfits.o: $(STDH) readfits.c
	$(CC) $(COPTS)  -c readfits.c

 read_config.o: $(STDH) read_config.c
	$(CC) $(COPTS)  -c read_config.c

 quick_sort.o: $(STDH) quick_sort.c
	$(CC) $(COPTS)  -c quick_sort.c

 get_thresh.o: $(STDH) get_thresh.c
	$(CC) $(COPTS)  -c get_thresh.c

 kill_cosmics.o: $(STDH) kill_cosmics.c
	$(CC) $(COPTS)  -c kill_cosmics.c

 writefits.o: $(STDH) writefits.c
	$(CC) $(COPTS)  -c writefits.c