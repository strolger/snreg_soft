void median(float x[], int n, float *xmed);

/* This routine hacked to expect a normal C zero-indexed array */
void median(x,n,xmed)
float x[],*xmed;
int n;
{
	int n2,n2p;
	void sort(int n, float ra[]);

	x = &x[-1];   /* convert to NR style 1-indexed array */
	sort(n,x);
	n2p=(n2=n/2)+1;
	*xmed=(n % 2 ? x[n2p] : 0.5*(x[n2]+x[n2p]));
}

/* This routine expects an NR 1-indexed array, first element ra[1] */
void sort(n,ra)
int n;
float ra[];
{
	int l,j,ir,i;
	float rra;

	l=(n >> 1)+1;
	ir=n;
	for (;;) {
		if (l > 1)
			rra=ra[--l];
		else {
			rra=ra[ir];
			ra[ir]=ra[1];
			if (--ir == 1) {
				ra[1]=rra;
				return;
			}
		}
		i=l;
		j=l << 1;
		while (j <= ir) {
			if (j < ir && ra[j] < ra[j+1]) ++j;
			if (rra < ra[j]) {
				ra[i]=ra[j];
				j += (i=j);
			}
			else j=ir+1;
		}
		ra[i]=rra;
	}
}

