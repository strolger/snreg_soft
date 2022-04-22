#CC    = cc                  # the C compiler
#COPTS =  -fast             # options for the C compiler
DEST  = ../bin
LIBS  =  -lm 
STDH  = functions.h image.h vectors.h

 ALL   = main.o allocate_vectors.o get_stamp_vectors.o get_stamps.o \
 make_bg_vectors.o max.o power.o read_image.o xy_conv_stamp.o build_matrix.o \
 build_scprod.o ludcmp.o lubksb.o kernel_vector.o make_kernel.o \
 kernel_convolve.o get_background.o cut_stamp.o defect_map.o \
 check_stamps.o build_matrix0.o build_scprod0.o check_again.o \
 make_model.o read_config.o read_header.o write_header.o \
 write_image.o quick_sort.o get_index.o write_image2.o

 
 mrj_phot: $(ALL)
	$(CC) $(ALL) -o mrj_phot $(LIBS) $(COPTS)
	cp mrj_phot $(DEST)


 main.o: $(STDH) main.c
	$(CC) $(COPTS)  -c main.c

 allocate_vectors.o: $(STDH) allocate_vectors.c
	$(CC) $(COPTS)  -c allocate_vectors.c

 get_stamp_vectors.o: $(STDH) get_stamp_vectors.c
	$(CC) $(COPTS)  -c get_stamp_vectors.c

 get_stamps.o: $(STDH) get_stamps.c
	$(CC) $(COPTS)  -c get_stamps.c

 make_bg_vectors.o: $(STDH) make_bg_vectors.c
	$(CC) $(COPTS)  -c make_bg_vectors.c
 
 max.o: $(STDH) max.c
	$(CC) $(COPTS)  -c max.c

 power.o: $(STDH) power.c
	$(CC) $(COPTS)  -c power.c

 read_image.o: $(STDH) read_image.c
	$(CC) $(COPTS)  -c read_image.c

 xy_conv_stamp.o: $(STDH) xy_conv_stamp.c
	$(CC) $(COPTS)  -c xy_conv_stamp.c

 build_matrix.o: $(STDH) build_matrix.c
	$(CC) $(COPTS)  -c build_matrix.c

 build_scprod.o: $(STDH) build_scprod.c
	$(CC) $(COPTS)  -c build_scprod.c

 ludcmp.o: $(STDH) ludcmp.c
	$(CC) $(COPTS)  -c ludcmp.c

 lubksb.o: $(STDH) lubksb.c
	$(CC) $(COPTS)  -c lubksb.c

 kernel_vector.o: $(STDH) kernel_vector.c
	$(CC) $(COPTS)  -c kernel_vector.c

 make_kernel.o: $(STDH) make_kernel.c
	$(CC) $(COPTS)  -c make_kernel.c

 kernel_convolve.o: $(STDH) kernel_convolve.c
	$(CC) $(COPTS)  -c kernel_convolve.c

 get_background.o: $(STDH) get_background.c
	$(CC) $(COPTS)  -c get_background.c

 cut_stamp.o: $(STDH) cut_stamp.c
	$(CC) $(COPTS)  -c cut_stamp.c

 defect_map.o: $(STDH) defect_map.c
	$(CC) $(COPTS)  -c defect_map.c

 check_stamps.o: $(STDH) check_stamps.c
	$(CC) $(COPTS)  -c check_stamps.c

 build_matrix0.o: $(STDH) build_matrix0.c
	$(CC) $(COPTS)  -c build_matrix0.c

 build_scprod0.o: $(STDH) build_scprod0.c
	$(CC) $(COPTS)  -c build_scprod0.c

 make_model.o: $(STDH) make_model.c
	$(CC) $(COPTS)  -c make_model.c

 check_again.o: $(STDH) check_again.c
	$(CC) $(COPTS)  -c check_again.c

 write_header.o: $(STDH) write_header.c
	$(CC) $(COPTS)  -c write_header.c

 write_image.o: $(STDH) write_image.c
	$(CC) $(COPTS)  -c write_image.c

 write_image2.o: $(STDH) write_image2.c
	$(CC) $(COPTS)  -c write_image2.c

 read_config.o: $(STDH) read_config.c
	$(CC) $(COPTS)  -c read_config.c

 read_header.o: $(STDH) read_header.c
	$(CC) $(COPTS)  -c read_header.c

 quick_sort.o: $(STDH) quick_sort.c
	$(CC) $(COPTS)  -c quick_sort.c

 get_index.o: $(STDH) get_index.c
	$(CC) $(COPTS)  -c get_index.c
	$(CC) $(COPTS)  -c get_index.c
	$(CC) $(COPTS)  -c get_index.c
	$(CC) $(COPTS)  -c get_index.c
