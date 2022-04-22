ORIG  = .			# where the sources lie
DEST  = ../bin			# where to put the executable
#CC    = cc 			# the C compiler
#COPTS =  -fast             	# options for the C compiler
LIBS  =  -lm 


ALL   =   newbincross2.o  poly.o



cross: $(ALL)
	$(CC) $(ALL) -o cross $(LIBS) $(COPTS)
	cp cross $(DEST)


newbincross2.o: $(STDH) newbincross2.c
	$(CC) $(COPTS)  -c newbincross2.c


poly.o: $(STDH) poly.c
	$(CC) $(COPTS)  -c poly.c







