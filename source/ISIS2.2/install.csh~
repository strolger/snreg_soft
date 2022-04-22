#! /bin/csh -f

 setenv CC gcc
 setenv COPTS -O2
 

 rm -r bin
 mkdir bin

 cp   short.h extract 
 cp   short.h subtract
 cp   short.h interp  
 cp   short.h Bphot 
 cp   short.h Cphot
 cp   short.h stack 
 cp   short.h phot_ref 
 cp   short.h utils
 cp   short.h detect
 cp   short.h abs
  
 cd   fit2d;     rm *.o; make
 cd ../phot_ref; rm *.o; make
 cd ../extract;  rm *.o; make
 cd ../cross;    rm *.o; make
 cd ../subtract; rm *.o; make
 cd ../interp;   rm *.o; make
 cd ../Bphot;    rm *.o; make
 cd ../Cphot;    rm *.o; make
 cd ../stack;    rm *.o; make
 cd ../utils;    rm *.o; $CC $COPTS hist2.c -lm;mv a.out ../bin/hist2
 cd ../abs  ; rm *.o; make
 cd ../detect  ; rm *.o; make
 cd ../czerny  ; rm *.o; make

