ORIG  = .			# where the sources lie
DEST  = ..			# where to put the executable
#CC    = gcc 			# the C compiler
#COPTS =  -fast            	# options for the C compiler
DEST  = ../bin
STDH  = types.h
LIBS  =  -lm 


 ALL   =  main.o writefits.o readfits.o read_header.o cent.o quick_sort.o \
	  clean_stack.o  allocate.o moy.o read_config.o get_background.o


 stack: $(ALL)
	$(CC) $(ALL) -o stack $(LIBS) $(COPTS)
	cp stack $(DEST)


 

 main.o: $(STDH) main.c
	$(CC) $(COPTS)  -c main.c

 writefits.o: $(STDH) writefits.c
	$(CC) $(COPTS)  -c writefits.c

 readfits.o: $(STDH) readfits.c
	$(CC) $(COPTS)  -c readfits.c

 read_header.o: $(STDH) read_header.c
	$(CC) $(COPTS)  -c read_header.c

 cent.o: $(STDH) cent.c
	$(CC) $(COPTS)  -c cent.c

 clean_stack.o: $(STDH)  clean_stack.c
	$(CC) $(COPTS)  -c  clean_stack.c

 quick_sort.o: $(STDH) quick_sort.c
	$(CC) $(COPTS)  -c quick_sort.c

 moy.o: $(STDH) moy.c
	$(CC) $(COPTS)  -c moy.c

 allocate.o: $(STDH) allocate.c
	$(CC) $(COPTS)  -c allocate.c

 read_config.o: $(STDH) read_config.c
	$(CC) $(COPTS)  -c read_config.c

 get_background.o: $(STDH) get_background.c
	$(CC) $(COPTS)  -c get_background.c
	$(CC) $(COPTS)  -c get_background.c
	$(CC) $(COPTS)  -c get_background.c
	$(CC) $(COPTS)  -c get_background.c
