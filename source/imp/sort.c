#include <stdio.h>

#define MAXSTACK 256
#define NSTOP 15

/* Sorting program that uses a quicksort algorithm 1980 (JT) */
qsort4(n, x)
     int n;
     float *x;
{
  float key, kl, kr, km, temp;
  int l, r, m, lstack[MAXSTACK], rstack[MAXSTACK], sp;
  register int i, j, k;
  int mgtl, lgtr, rgtm;

  sp = 0;
  lstack[sp] = 0;
  rstack[sp] = n-1;

  while(n > NSTOP && sp >= 0) {
/* Sort a subrecord off the stack */
    l = lstack[sp];
    r = rstack[sp];
    sp--;
    m = (l + r) / 2;
/* Set KEY = median of X(L), X(M), X(R) */
    kl = x[l];
    km = x[m];
    kr = x[r];
    mgtl = km > kl;
    rgtm = kr > km;
    lgtr = kl > kr;
    /*
    if(mgtl ^ rgtm) {
      if(mgtl ^ lgtr) key = kr;
      else            key = kl;
    } else {
      key = km;
    }
    */
/* Curiously enough, this non-median seems to work as well or better */
    if(mgtl ^ rgtm) {
      key = km;
    } else {
      if(mgtl ^ lgtr) key = kl;
      else            key = kr;
    }

    i = l;
    j = r;
    while(1) {
/* Find a big record on the left */
      while(x[i] < key) i++;

/* Find a small record on the right */
      while(x[j] > key) j--;

      if(i >= j) break;
/* Exchange records */
      temp = x[i];
      x[i] = x[j];
      x[j] = temp;
      i++;
      j--;
    }

/* Subfile is partitioned into two halves, left .le. right */
/* Push the two halves on the stack */
    if(j-l+1 > NSTOP) {
      sp = sp + 1;
      lstack[sp] = l;
      rstack[sp] = j;
    }
    if(r-j > NSTOP) {
      sp = sp + 1;
      lstack[sp] = j+1;
      rstack[sp] = r;
    }
    if(sp >= MAXSTACK) {
      fprintf(stderr,"QSORT4: Fatal error from stack overflow\n");
      fprintf(stderr,"Fall back on sort by insertion\n");
      break;
    }
  }

/* Sorting routine that sorts the N elements of single precision */
/* array X by straight insertion between previously sorted numbers */
  for(j=n-2; j>=0; j--) {
    k = j;
    for(i=j+1; i<n; i++) {
      if(x[j] <= x[i]) break;
      k = i;
    }
    if(k != j) {
      temp = x[j];
      for(i=j+1; i<=k; i++) x[i-1] = x[i];
      x[k] = temp;
    }
  }
}
