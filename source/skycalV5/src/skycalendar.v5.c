/*
   This is a self-contained c-language program to print a nighttime
   astronomical calendar for use in planning observations.  
   It prints to standard output (usually the terminal); the
   operator should capture this output (e. g., using redirection
   in UNIX or the /out= switch in VMS) and then print it on an
   appropriate output device.  The table which is printed is some
   125 columns wide, so a wide device is required (either a line
   printer or a laserprinter in LANDSCAPE mode.)  It is assumed that
   the ASCII form-feed character will actually begin a new page. 
   The original program was to run on VMS, but it should be very 
   transportable.  Non-vms users will probably want to change 
   'unixio.h' to 'stdio.h' in the first line.
   An explanatory text is printed at the beginning of the output, which
   includes the appropriate CAUTIONS regarding accuracy and applicability.

   A number of 'canned site' parameters have been included.  Be
   careful of time zones, DST etc. for foreign sites.  
   To customize to your own site, install an option in the
   routine 'load_site'.  The code is very straightforward; just do
   it exactly the same as the others.  You might also want to erase
   some seldom-used choices.  One can also specify new site parameters
   at run time.

   This program may be used freely by anyone for scientific or educational
   purposes.  If you use it for profit, I want a cut, and claim 
   a copyright herewith.  In any case please acknowledge the source:

			John Thorstensen
			Dept. of Physics and Astronomy
			Dartmouth College
			Hanover, NH 03755
			John.Thorstensen@dartmouth.edu

			May 26, 1993.

*/

#include <stdio.h>
#include <math.h>

double pi = 3.14159265358979;

#define PI                3.14159265358979
#define ARCSEC_IN_RADIAN  206264.8062
#define DEG_IN_RADIAN     57.2957795130823
#define HRS_IN_RADIAN     3.819718634
#define J2000             2451545.
#define SEC_IN_DAY        86400.
#define FLATTEN           0.003352813 
#define EQUAT_RAD         6378137.
#define TWILIGHT_ALT      -18.



struct coord
   {
     short sign;  /* carry sign explicitly since -0 not neg. */
     double hh;
     double mm;
     double ss;
   };

struct coord_pair
   {
      struct coord RA;
      struct coord dec;
   };

double bab_to_dec (bab)
   struct coord bab;
   {
   double x;
   x = bab.sign * (bab.hh + bab.mm / 60. + bab.ss / 3600.);
   return(x);
   }

void dec_to_bab (deci,bab)

   /* function for converting decimal to babylonian hh mm ss.ss */

double deci;
struct coord *bab;

{
   int hr_int, min_int;

   if (deci >= 0.) bab->sign = 1; 
   else {
      bab->sign = -1;
      deci = -1. * deci;
   }
   hr_int = deci;   /* use conversion conventions to truncate */
   bab->hh = hr_int;
   min_int = 60. * (deci - bab->hh);
   bab->mm = min_int;
   bab->ss = 3600. * (deci - bab->hh - bab->mm / 60.);
}

double get_coord()

/* Reads a string from the terminal and converts it into
   a double-precision coordinate.  This is trickier than 
   it appeared at first, since a -00 tests as non-negative; 
   the sign has to be picked out and handled explicitly. */
/* Prompt for input in the calling routine.*/
{
   short sign;
   double hrs, mins, secs;
   char hh_string[6];  /* string with the first coord (hh) */
   char hh1[1];
   short i = 0;

   /* read and handle the hour (or degree) part with sign */

   scanf("%s",hh_string);
   hh1[0] = hh_string[i];

   while(hh1[0] == ' ') {
       /* discard leading blanks */
       i++;
       hh1[0] = hh_string[i];
   }

   if(hh1[0] == '-') sign = -1;

     else sign = 1;

   sscanf(hh_string,"%lf", &hrs);
   if(sign == -1) hrs = -1. * hrs;

   /* read in the minutes and seconds normally */
   scanf("%lf %lf",&mins,&secs);

   return(sign * (hrs + mins / 60. + secs / 3600.));
}

void put_coords(deci, precision)

   double deci;  /* decimal coordinate */
   short precision;

/* prints out a struct coord in a nice format; precision
   is a code for how accurate you want it.  The options are:
     precision = 0;   minutes rounded to the nearest minute
     precision = 1;   minutes rounded to the nearest tenth.
     precision = 2;   seconds rounded to the nearest second
     precision = 3;   seconds given to the tenth
     precision = 4;   seconds given to the hundredth
   The program assumes that the line is ready for the coord
   to be printed and does NOT deliver a new line at the end
   of the output. */
 
{
   
   double minutes;  /* for rounding off if necess. */
   struct coord out_coord, coords;
   char out_string[20];  /* for checking for nasty 60's */

   dec_to_bab(deci,&coords);  /* internally convert to coords*/

   if(coords.sign == -1) printf("-");
	else printf(" "); /* to preserve alignment */

   if(precision == 0) {   /* round to nearest minute */
      minutes = coords.mm + coords.ss / 60.;
           /* check to be sure minutes aren't 60 */
      sprintf(out_string,"%.0f %02.0f",coords.hh,minutes);
      sscanf(out_string,"%lf %lf",&out_coord.hh,&out_coord.mm);
      if(fabs(out_coord.mm - 60.) < 1.0e-7) {
         out_coord.mm = 0.;
         out_coord.hh = out_coord.hh + 1.;
      }
      printf("%2.0f %02.0f",out_coord.hh,out_coord.mm);
   }

   else if(precision == 1) {    /* keep nearest tenth of a minute */
      minutes = coords.mm + coords.ss / 60.;
           /* check to be sure minutes are not 60 */
      sprintf(out_string,"%.0f %04.1f",coords.hh,minutes);
      sscanf(out_string,"%lf %lf",&out_coord.hh, &out_coord.mm);
      if(fabs(out_coord.mm - 60.) < 1.0e-7) {
         out_coord.mm = 0.;
         out_coord.hh = out_coord.hh + 1.;
      }
      printf("%2.0f %04.1f", out_coord.hh, out_coord.mm);
   }
   else if(precision == 2) {
          /* check to be sure seconds are not 60 */
      sprintf(out_string,"%.0f %02.0f %02.0f",coords.hh,coords.mm,coords.ss);
      sscanf(out_string,"%lf %lf %lf",&out_coord.hh,&out_coord.mm,
           &out_coord.ss);
      if(fabs(out_coord.ss - 60.) < 1.0e-7) {
          out_coord.mm = out_coord.mm + 1.;
          out_coord.ss = 0.;
          if(fabs(out_coord.mm - 60.) < 1.0e-7) {
              out_coord.hh = out_coord.hh + 1.;
              out_coord.mm = 0.;
          }
      }
      printf("%2.0f %02.0f %02.0f",out_coord.hh,out_coord.mm,out_coord.ss);
   }
   else if(precision == 3) {
          /* the usual shuffle to check for 60's */
      sprintf(out_string,"%.0f %02.0f %04.1f",coords.hh, coords.mm, coords.ss);
      sscanf(out_string,"%lf %lf %lf",&out_coord.hh,&out_coord.mm,
           &out_coord.ss);
      if(fabs(out_coord.ss - 60.) < 1.0e-7) {
          out_coord.mm = out_coord.mm + 1.;
          out_coord.ss = 0.;
          if(fabs(out_coord.mm - 60.) < 1.0e-7) {
             out_coord.hh = out_coord.hh + 1.;
             out_coord.mm = 0.;
          }
      }
      printf("%2.0f %02.0f %04.1f",out_coord.hh,out_coord.mm,out_coord.ss);
   }
   else {
      sprintf(out_string,"%.0f %02.0f %05.2f",coords.hh,coords.mm,coords.ss);
      sscanf(out_string,"%lf %lf %lf",&out_coord.hh,&out_coord.mm,
           &out_coord.ss);
      if(fabs(out_coord.ss - 60.) < 1.0e-6) {
         out_coord.mm = out_coord.mm + 1.;
         out_coord.ss = 0.;
         if(fabs(out_coord.mm - 60.) < 1.0e-6) {
            out_coord.hh = out_coord.hh + 1.;
            out_coord.mm = 0.;
         }
      }
      printf("%2.0f %02.0f %05.2f",out_coord.hh, out_coord.mm, out_coord.ss);
   }
}


double atan_circ(x,y)

double x, y;

{
	/* returns radian angle 0 to 2pi for coords x, y */

	double theta;

	if(x == 0) {
		if(y > 0.) theta = PI/2.;
		else if(y < 0.) theta = 3.*PI/2.;
		else theta = 0.;   /* x and y zero */
	}
	else theta = atan(y/x);
	if(x < 0.) theta = theta + PI;
	if(theta < 0.) theta = theta + 2.*PI;
	return(theta);
}

double altit(dec,ha,lat)

double dec,ha,lat;  /* dec deg, dec hrs, dec deg */

{
	double x;
	dec = dec / DEG_IN_RADIAN;
	ha = ha / HRS_IN_RADIAN;
	lat = lat / DEG_IN_RADIAN;  /* thank heavens for pass-by-value */
	x = DEG_IN_RADIAN * asin(cos(dec)*cos(ha)*cos(lat) + sin(dec)*sin(lat));
	return(x);
}

void min_max_alt(lat,dec,min,max)

double lat, dec, *min, *max;

{
	/* computes minimum and maximum altitude for a given dec and
            latitude. */

	double x;
	lat = lat / DEG_IN_RADIAN; /* pass by value! */
	dec = dec / DEG_IN_RADIAN;
	x = cos(dec)*cos(lat) + sin(dec)*sin(lat);
	if(fabs(x) <= 1.) {
		*max = asin(x) * DEG_IN_RADIAN;
	}
	else printf("Error in min_max_alt -- arcsin(>1)\n");
	x = sin(dec)*sin(lat) - cos(dec)*cos(lat);
	if(fabs(x) <= 1.) {
		*min = asin(x) * DEG_IN_RADIAN;
	}
	else printf("Error in min_max_alt -- arcsin(>1)\n");
}

double ha_alt(dec,lat,alt)

double dec,lat,alt;  /* dec deg */

{
	/* returns hour angle at which object at dec is at altitude alt */
	
	double x,coalt,min,max;
	
	min_max_alt(lat,dec,&min,&max);
	if(alt < min) 
		return(1000.);  /* flag value - always higher than asked */
	if(alt > max)
		return(-1000.); /* flag for object always lower than asked */
	dec = (0.5*PI) - dec / DEG_IN_RADIAN;
	lat = (0.5*PI) - lat / DEG_IN_RADIAN;
	coalt = (0.5*PI) - alt / DEG_IN_RADIAN;
	x = (cos(coalt) - cos(dec)*cos(lat)) / (sin(dec)*sin(lat));
	if(fabs(x) <= 1.) return(acos(x) * HRS_IN_RADIAN);
	else {
		printf("Error in ha_alt ... acos(>1).\n");
		return(1000.);
	}
}

double subtend(ra1,dec1,ra2,dec2)

double ra1, dec1, ra2, dec2;  /* dec hrs and dec degrees */

{
	/* angle subtended by two directions in the sky. */

	double x1, y1, z1, x2, y2, z2;
	double theta;

	ra1 = ra1 / HRS_IN_RADIAN;
	dec1 = dec1 / DEG_IN_RADIAN;
	ra2 = ra2 / HRS_IN_RADIAN;
	dec2 = dec2 / DEG_IN_RADIAN;
	x1 = cos(ra1)*cos(dec1);
	y1 = sin(ra1)*cos(dec1);
	z1 = sin(dec1);
	x2 = cos(ra2)*cos(dec2);
	y2 = sin(ra2)*cos(dec2);
	z2 = sin(dec2);
   	theta = acos(x1*x2+y1*y2+z1*z2);
	return(theta);
}

struct date_time
   {
	short y;
	short mo;
	short d;
	short h;
	short mn;
	float s;
   };

double date_to_jd(date)

struct date_time date;

{
	short yr1=0, mo1=1;	
	long jdzpt = 1720982, jdint, inter;
	double jd,jdfrac;


	if((date.y <= 1900) | (date.y >= 2100)) {
		printf("Date out of range.  1900 - 2100 only.\n");
		return(0.);
	}
	
	if(date.mo <= 2) {
		yr1 = -1;
		mo1 = 13;
	}

	jdint = 365.25*(date.y+yr1);  /* truncates */
        inter = 30.6001*(date.mo+mo1);
	jdint = jdint+inter+date.d+jdzpt;
	jd = jdint;
	jdfrac=date.h/24.+date.mn/1440.+date.s/86400;
	if(jdfrac < 0.5) {
		jdint--;
		jdfrac=jdfrac+0.5;
	}
	else jdfrac=jdfrac-0.5;			
	jd=jdint+jdfrac;
	return(jd);
}

void caldat(jdin,date)

#define IGREG 2299161

double jdin;
struct date_time *date;

{ 
	/* Returns date and time for a given julian date;
	   Adapted from Numerical Recipes, p. 12. */

	int mm, id, iyyy;  /* their notation */
	long ja, jdint, jalpha, jb, jc, jd, je;
	float jdfrac;

	jdin = jdin + 0.5;  /* adjust for 1/2 day */
	jdint = jdin;
	jdfrac = jdin - jdint;
	date->h = jdfrac * 24; /* truncate */
	date->mn = (jdfrac - ((float) date->h)/24.) * 1440.;
	date->s = (jdfrac - ((float) date->h)/24. - 
			((float) date->mn)/1440.) * 86400;
	
	if(jdint > IGREG) {
		jalpha=((float) (jdint-1867216)-0.25)/36524.25;
		ja=jdint+1+jalpha-(long)(0.25*jalpha);
	}
	else
		ja=jdint;
	jb=ja+1524;
	jc=6680.0+((float) (jb-2439870)-122.1)/365.25;
	jd=365*jc+(0.25*jc);
	je=(jb-jd)/30.6001;
	id=jb-jd-(int) (30.6001*je);
	mm=je-1;
	if(mm > 12) mm -= 12;
	iyyy=jc-4715;
	if(mm > 2) --iyyy;
	if (iyyy <= 0) --iyyy;
	date->y = iyyy;
	date->mo = mm;
	date->d = id;
}

void print_calendar(jdin,length)

double jdin;
short length;

{
	struct date_time date;
	char *months = "JanFebMarAprMayJunJulAugSepOctNovDec";
	char mo_out[4];

	caldat(jdin,&date);	
	mo_out[0] = *(months + 3*(date.mo - 1));
	mo_out[1] = *(months + 3*(date.mo - 1) + 1);
	mo_out[2] = *(months + 3*(date.mo - 1) + 2);
	mo_out[3] = '\0';
	if(length == 1) 
		printf("%d %s %d",date.y,mo_out,date.d);
	else printf("%s %02d",mo_out,date.d); /*no year */
}

void print_time(jdin,prec)

double jdin;
short prec;

{
	/* prints time only */
	struct date_time date;
	double temptime;
	
	caldat(jdin,&date);

	temptime = date.h + date.mn/60. + date.s/3600.;

	put_coords(temptime,prec);
}	

short day_of_week(jd)

double jd;

{ 
	/* returns day of week, 0 = Mon, 6 = Sun. */

	double x,y;
	long i;
	short d;
	
	x = (jd+0.5)/7.;
	d = 7.*(x - (long) x);   /* truncate */
	return(d);
}

void print_day(d)

short d;

{
	char *days = "MonTueWedThuFriSatSun";
	char day_out[4];

	day_out[0] = *(days+3*d);
	day_out[1] = *(days+3*d+1);
	day_out[2] = *(days+3*d+2);
	day_out[3] = '\0';  /* terminate with null char */
	
	printf("%s",day_out);
}


double lst(jd,longit) 

double jd,longit;

{
	/* returns the local MEAN sidereal time (dec hrs) at julian date jd
	   at west longitude long (decimal hours).  Follows
           definitions in 1992 Astronomical Almanac, pp. B7 and L2. 
           Expression for GMST at 0h ut referenced to Aoki et al, A&A 105,
	   p.359, 1982. */

	double t, ut, jdmid, jdint, jdfrac, sid_g, sid;
	double jdnoon2000jan1 = 2451545.;
	long jdin, sid_int;

	jdin = jd;         /* fossil code from earlier package which 
			split jd into integer and fractional parts ... */
	jdint = jdin;
	jdfrac = jd - jdint;
	if(jdfrac < 0.5) {
		jdmid = jdint - 0.5;
		ut = jdfrac + 0.5;
	}
	else {
		jdmid = jdint + 0.5;
		ut = jdfrac - 0.5;
	}
	t = (jdmid - jdnoon2000jan1)/36525;
	sid_g = (24110.54841+8640184.812866*t+0.093104*t*t-6.2e-6*t*t*t)/86400.;
	sid_int = sid_g;
	sid_g = sid_g - (double) sid_int;
	sid_g = sid_g + 1.0027379093 * ut - longit/24.;
	sid_int = sid_g;
	sid_g = (sid_g - (double) sid_int) * 24.;
	if(sid_g < 0.) sid_g = sid_g + 24.;
	return(sid_g);
}

double adj_time(x)

double x;

{
	/* adjusts a time (decimal hours) to be between -12 and 12. */

	if(fabs(x) < 100000.) {  /* too inefficient for this! */
		while(x > 12.) {
			x = x - 24.;
		}
		while(x < -12.) {
			x = x + 24.;
		}
	}
	else printf("Out of bounds in adj_time!\n");
	return(x);
}


double circulo(x)
	
double x;

{
	/* assuming x is an angle in degrees, returns 
	   modulo 360 degrees. */

	int n;

	n = (int)(x / 360.);
	return(x - 360. * n);
}	


void geocent(geolong, geolat, height, x_geo, y_geo, z_geo)

double geolong, geolat, height, *x_geo, *y_geo, *z_geo;

/* computes the geocentric coordinates from the geodetic 
(standard map-type) longitude, latitude, and height. 
These are assumed to be in decimal hours, decimal degrees, and
meters respectively.  Notation generally follows 1992 Astr Almanac, 
p. K11 */


{
	
	double denom, C_geo, S_geo;

	geolat = geolat / DEG_IN_RADIAN;
	geolong = geolong / HRS_IN_RADIAN;      
	denom = (1. - FLATTEN) * sin(geolat);
	denom = cos(geolat) * cos(geolat) + denom*denom;
	C_geo = 1. / sqrt(denom);
	S_geo = (1. - FLATTEN) * (1. - FLATTEN) * C_geo;
	C_geo = C_geo + height / EQUAT_RAD;  /* deviation from almanac
                       notation -- include height here. */
	S_geo = S_geo + height / EQUAT_RAD;
	*x_geo = C_geo * cos(geolat) * cos(geolong);
	*y_geo = C_geo * cos(geolat) * sin(geolong);
	*z_geo = S_geo * sin(geolat);
}

void eclrot(jd, x, y, z)

/* rotates ecliptic rectangular coords x, y, z to
   equatorial (all assumed of date.) */

double jd, *x, *y, *z;

{
	double incl;
	double xpr,ypr,zpr;
	double T;

	T = (jd - J2000) / 36525;  /* centuries since J2000 */
	
	incl = (23.439291 + T * (-0.0130042 - 0.00000016 * T))/DEG_IN_RADIAN; 
		/* 1992 Astron Almanac, p. B18, dropping the 
                   cubic term, which is 2 milli-arcsec! */
	ypr = cos(incl) * *y - sin(incl) * *z;
	zpr = sin(incl) * *y + cos(incl) * *z;
	*y = ypr;
	*z = zpr;
	/* x remains the same. */	
}

double etcorr(jd)

double jd;

{

	/* Given a julian date in 1900-2100, returns the jd corrected
           for delta t; delta t  is 
		TDT - UT (after 1983 and before 1994)
		ET - UT (before 1983)
		an extrapolated guess  (after 1994). 

	For dates in the past (<= 1994 and after 1900) the value is linearly
        interpolated on 5-year intervals; for dates after the present,
        an extrapolation is used, because the true value of delta t
	cannot be predicted precisely.  Note that TDT is essentially the
	modern version of ephemeris time with a slightly cleaner 
	definition.  

	Where the algorithm shifts there is an approximately 0.1 second
        discontinuity.  Also, the 5-year linear interpolation scheme can 
        lead to errors as large as 0.5 seconds in some cases, though
 	usually rather smaller. */

	double jd1900 = 2415019.5;
	double dates[20];
	double delts[20];  /* can't initialize this look-up table
            with stupid old sun compiler .... */
	double year, delt;
	int i;
	
	/* this stupid patch for primitive sun C compilers .... 
		do not allow automatic initialization of arrays! */

	for(i = 0; i <= 18; i++) dates[i] = 1900 + (double) i * 5.;
	dates[19] = 1994;


	delts[0] = -2.72;  delts[1] = 3.86; delts[2] = 10.46;
	delts[3] = 17.20;  delts[4] = 21.16; delts[5] = 23.62;
	delts[6] = 24.02;  delts[7] = 23.93; delts[8] = 24.33;
	delts[9] = 26.77;  delts[10] = 29.15; delts[11] = 31.07;
	delts[12] = 33.15;  delts[13] = 35.73; delts[14] = 40.18;
	delts[15] = 45.48;  delts[16] = 50.54; delts[17] = 54.34;
	delts[18] = 56.86;  delts[19] = 59.98;

	year = 1900. + (jd - 2415019.5) / 365.25;

	if(year < 1994.0 && year >= 1900.) {
		i = (year - 1900) / 5;
		delt = delts[i] + 
		 ((delts[i+1] - delts[i])/(dates[i+1] - dates[i])) * (year - dates[i]);
	}

	else if (year > 1994. && year < 2100.)
		delt = 33.15 + (2.164e-3) * (jd - 2436935.4);  /* rough extrapolation */

	else if (year < 1900) {
		printf("etcorr ... no ephemeris time data for < 1900.\n");
       		delt = 0.;
	}

	else if (year >= 2100.) {
		printf("etcorr .. very long extrapolation in delta T - inaccurate.\n");
		delt = 180.; /* who knows? */
	} 

	return(jd + delt/SEC_IN_DAY);
}
	
void accumoon(jd,geolat,lst,elevsea,topora,topodec,topodist)

double jd,geolat,lst,elevsea;  /* jd, dec. degr., dec. hrs., meters */
double *topora,*topodec,*topodist;

/* More accurate (but more elaborate and slower) lunar 
   ephemeris, from Jean Meeus' *Astronomical Formulae For Calculators*,
   pub. Willman-Bell.  Includes all the terms given there. */

{	
/*	double *eclatit,*eclongit, *pie,*ra,*dec,*dist; geocent quantities,
		formerly handed out but not in this version */
	double pie, dist;  /* horiz parallax */
	double Lpr,M,Mpr,D,F,Om,T,Tsq,Tcb;
	double e,lambda,B,beta,om1,om2;
	double sinx, x, y, z, l, m, n;
	double x_geo, y_geo, z_geo;  /* geocentric position of *observer* */	

	jd = etcorr(jd);   /* approximate correction to ephemeris time */
	T = (jd - 2415020.) / 36525.;   /* this based around 1900 ... */
	Tsq = T * T;
	Tcb = Tsq * T;

	Lpr = 270.434164 + 481267.8831 * T - 0.001133 * Tsq 
			+ 0.0000019 * Tcb;
	M = 358.475833 + 35999.0498*T - 0.000150*Tsq
			- 0.0000033*Tcb;
	Mpr = 296.104608 + 477198.8491*T + 0.009192*Tsq 
			+ 0.0000144*Tcb;
	D = 350.737486 + 445267.1142*T - 0.001436 * Tsq
			+ 0.0000019*Tcb;
	F = 11.250889 + 483202.0251*T -0.003211 * Tsq 
			- 0.0000003*Tcb;
	Om = 259.183275 - 1934.1420*T + 0.002078*Tsq 
			+ 0.0000022*Tcb;

	Lpr = circulo(Lpr);
	Mpr = circulo(Mpr);	
	M = circulo(M);
	D = circulo(D);
	F = circulo(F);
	Om = circulo(Om);

	
	sinx =  sin((51.2 + 20.2 * T)/DEG_IN_RADIAN);
	Lpr = Lpr + 0.000233 * sinx;
	M = M - 0.001778 * sinx;
	Mpr = Mpr + 0.000817 * sinx;
	D = D + 0.002011 * sinx;
	
	sinx = 0.003964 * sin((346.560+132.870*T -0.0091731*Tsq)/DEG_IN_RADIAN);

	Lpr = Lpr + sinx;
	Mpr = Mpr + sinx;
	D = D + sinx;
	F = F + sinx;

	sinx = sin(Om/DEG_IN_RADIAN);
	Lpr = Lpr + 0.001964 * sinx;
	Mpr = Mpr + 0.002541 * sinx;
	D = D + 0.001964 * sinx;
	F = F - 0.024691 * sinx;
	F = F - 0.004328 * sin((Om + 275.05 -2.30*T)/DEG_IN_RADIAN);

	e = 1 - 0.002495 * T - 0.00000752 * Tsq;

	M = M / DEG_IN_RADIAN;   /* these will all be arguments ... */
	Mpr = Mpr / DEG_IN_RADIAN;
	D = D / DEG_IN_RADIAN;
	F = F / DEG_IN_RADIAN;

	lambda = Lpr + 6.288750 * sin(Mpr)
		+ 1.274018 * sin(2*D - Mpr)
		+ 0.658309 * sin(2*D)
		+ 0.213616 * sin(2*Mpr)
		- e * 0.185596 * sin(M) 
		- 0.114336 * sin(2*F)
		+ 0.058793 * sin(2*D - 2*Mpr)
		+ e * 0.057212 * sin(2*D - M - Mpr)
		+ 0.053320 * sin(2*D + Mpr)
		+ e * 0.045874 * sin(2*D - M)
		+ e * 0.041024 * sin(Mpr - M)
		- 0.034718 * sin(D)
		- e * 0.030465 * sin(M+Mpr)
		+ 0.015326 * sin(2*D - 2*F)
		- 0.012528 * sin(2*F + Mpr)
		- 0.010980 * sin(2*F - Mpr)
		+ 0.010674 * sin(4*D - Mpr)
		+ 0.010034 * sin(3*Mpr)
		+ 0.008548 * sin(4*D - 2*Mpr)
		- e * 0.007910 * sin(M - Mpr + 2*D)
		- e * 0.006783 * sin(2*D + M)
		+ 0.005162 * sin(Mpr - D);

		/* And furthermore.....*/

	lambda = lambda + e * 0.005000 * sin(M + D)
		+ e * 0.004049 * sin(Mpr - M + 2*D)
		+ 0.003996 * sin(2*Mpr + 2*D)
		+ 0.003862 * sin(4*D)
		+ 0.003665 * sin(2*D - 3*Mpr)
		+ e * 0.002695 * sin(2*Mpr - M)
		+ 0.002602 * sin(Mpr - 2*F - 2*D)
		+ e * 0.002396 * sin(2*D - M - 2*Mpr)
		- 0.002349 * sin(Mpr + D)
		+ e * e * 0.002249 * sin(2*D - 2*M)
		- e * 0.002125 * sin(2*Mpr + M)
		- e * e * 0.002079 * sin(2*M)
		+ e * e * 0.002059 * sin(2*D - Mpr - 2*M)
		- 0.001773 * sin(Mpr + 2*D - 2*F)
		- 0.001595 * sin(2*F + 2*D)
		+ e * 0.001220 * sin(4*D - M - Mpr)
		- 0.001110 * sin(2*Mpr + 2*F)
		+ 0.000892 * sin(Mpr - 3*D)
		- e * 0.000811 * sin(M + Mpr + 2*D)
		+ e * 0.000761 * sin(4*D - M - 2*Mpr)
		+ e * e * 0.000717 * sin(Mpr - 2*M)
		+ e * e * 0.000704 * sin(Mpr - 2 * M - 2*D)
		+ e * 0.000693 * sin(M - 2*Mpr + 2*D)
		+ e * 0.000598 * sin(2*D - M - 2*F)
		+ 0.000550 * sin(Mpr + 4*D)
		+ 0.000538 * sin(4*Mpr)
		+ e * 0.000521 * sin(4*D - M)
		+ 0.000486 * sin(2*Mpr - D);
	
/*		*eclongit = lambda;  */

	B = 5.128189 * sin(F)
		+ 0.280606 * sin(Mpr + F)
		+ 0.277693 * sin(Mpr - F)
		+ 0.173238 * sin(2*D - F)
		+ 0.055413 * sin(2*D + F - Mpr)
		+ 0.046272 * sin(2*D - F - Mpr)
		+ 0.032573 * sin(2*D + F)
		+ 0.017198 * sin(2*Mpr + F)
		+ 0.009267 * sin(2*D + Mpr - F)
		+ 0.008823 * sin(2*Mpr - F)
		+ e * 0.008247 * sin(2*D - M - F) 
		+ 0.004323 * sin(2*D - F - 2*Mpr)
		+ 0.004200 * sin(2*D + F + Mpr)
		+ e * 0.003372 * sin(F - M - 2*D)
		+ 0.002472 * sin(2*D + F - M - Mpr)
		+ e * 0.002222 * sin(2*D + F - M)
		+ e * 0.002072 * sin(2*D - F - M - Mpr)
		+ e * 0.001877 * sin(F - M + Mpr)
		+ 0.001828 * sin(4*D - F - Mpr)
		- e * 0.001803 * sin(F + M)
		- 0.001750 * sin(3*F)
		+ e * 0.001570 * sin(Mpr - M - F)
		- 0.001487 * sin(F + D)
		- e * 0.001481 * sin(F + M + Mpr)
		+ e * 0.001417 * sin(F - M - Mpr)
		+ e * 0.001350 * sin(F - M)
		+ 0.001330 * sin(F - D)
		+ 0.001106 * sin(F + 3*Mpr)
		+ 0.001020 * sin(4*D - F)
		+ 0.000833 * sin(F + 4*D - Mpr);
     /* not only that, but */
	B = B + 0.000781 * sin(Mpr - 3*F)
		+ 0.000670 * sin(F + 4*D - 2*Mpr)
		+ 0.000606 * sin(2*D - 3*F)
		+ 0.000597 * sin(2*D + 2*Mpr - F)
		+ e * 0.000492 * sin(2*D + Mpr - M - F)
		+ 0.000450 * sin(2*Mpr - F - 2*D)
		+ 0.000439 * sin(3*Mpr - F)
		+ 0.000423 * sin(F + 2*D + 2*Mpr)
		+ 0.000422 * sin(2*D - F - 3*Mpr)
		- e * 0.000367 * sin(M + F + 2*D - Mpr)
		- e * 0.000353 * sin(M + F + 2*D)
		+ 0.000331 * sin(F + 4*D)
		+ e * 0.000317 * sin(2*D + F - M + Mpr)
		+ e * e * 0.000306 * sin(2*D - 2*M - F)
		- 0.000283 * sin(Mpr + 3*F);
	
	om1 = 0.0004664 * cos(Om/DEG_IN_RADIAN);	
	om2 = 0.0000754 * cos((Om + 275.05 - 2.30*T)/DEG_IN_RADIAN);
	
	beta = B * (1. - om1 - om2);
/*      *eclatit = beta; */
	
	pie = 0.950724 
		+ 0.051818 * cos(Mpr)
		+ 0.009531 * cos(2*D - Mpr)
		+ 0.007843 * cos(2*D)
		+ 0.002824 * cos(2*Mpr)
		+ 0.000857 * cos(2*D + Mpr)
		+ e * 0.000533 * cos(2*D - M)
		+ e * 0.000401 * cos(2*D - M - Mpr)
		+ e * 0.000320 * cos(Mpr - M)
		- 0.000271 * cos(D)
		- e * 0.000264 * cos(M + Mpr)
		- 0.000198 * cos(2*F - Mpr)
		+ 0.000173 * cos(3*Mpr)
		+ 0.000167 * cos(4*D - Mpr)
		- e * 0.000111 * cos(M)
		+ 0.000103 * cos(4*D - 2*Mpr)
		- 0.000084 * cos(2*Mpr - 2*D)
		- e * 0.000083 * cos(2*D + M)
		+ 0.000079 * cos(2*D + 2*Mpr)
		+ 0.000072 * cos(4*D)
		+ e * 0.000064 * cos(2*D - M + Mpr)
		- e * 0.000063 * cos(2*D + M - Mpr)
		+ e * 0.000041 * cos(M + D)
		+ e * 0.000035 * cos(2*Mpr - M)
		- 0.000033 * cos(3*Mpr - 2*D)
		- 0.000030 * cos(Mpr + D)
		- 0.000029 * cos(2*F - 2*D)
		- e * 0.000029 * cos(2*Mpr + M)
		+ e * e * 0.000026 * cos(2*D - 2*M)
		- 0.000023 * cos(2*F - 2*D + Mpr)
		+ e * 0.000019 * cos(4*D - M - Mpr);

	beta = beta/DEG_IN_RADIAN;
	lambda = lambda/DEG_IN_RADIAN;
	l = cos(lambda) * cos(beta);	
	m = sin(lambda) * cos(beta);
	n = sin(beta);
	eclrot(jd,&l,&m,&n);
	
	dist = 1/sin((pie)/DEG_IN_RADIAN);
	x = l * dist;
	y = m * dist;
	z = n * dist;

/*	*ra = atan_circ(l,m) * DEG_IN_RADIAN;
	*dec = asin(n) * DEG_IN_RADIAN;        */

	geocent(lst,geolat,elevsea,&x_geo,&y_geo,&z_geo);
	
	x = x - x_geo;  /* topocentric correction using elliptical earth fig. */
	y = y - y_geo;
	z = z - z_geo;

	*topodist = sqrt(x*x + y*y + z*z);
	
	l = x / (*topodist);
	m = y / (*topodist);
	n = z / (*topodist);

	*topora = atan_circ(l,m) * HRS_IN_RADIAN;
	*topodec = asin(n) * DEG_IN_RADIAN; 

}


lpsun(jd,ra,dec)

double jd, *ra, *dec;

/* Low precision formulae for the sun, from Almanac p. C24 (1990) */
/* ra and dec are returned as decimal hours and decimal degrees. */

{
	double n, L, g, lambda,epsilon,alpha,delta,x,y,z;

	n = jd - 2451545.0;
	L = 280.460 + 0.9856474 * n;
	g = (357.528 + 0.9856003 * n)/DEG_IN_RADIAN;
	lambda = (L + 1.915 * sin(g) + 0.020 * sin(2. * g))/DEG_IN_RADIAN;
	epsilon = (23.439 - 0.0000004 * n)/DEG_IN_RADIAN;

	x = cos(lambda); 
	y = cos(epsilon) * sin(lambda); 
	z = sin(epsilon)*sin(lambda);

	*ra = (atan_circ(x,y))*HRS_IN_RADIAN;
	*dec = (asin(z))*DEG_IN_RADIAN;
}

double jd_moon_alt(alt,jdguess,lat,longit)

double alt, jdguess, lat, longit;

{
	/* returns jd at which moon is at a given 
		altitude, given jdguess as a starting point. */

	double jdout;
	double deriv, err, del = 0.002;
	double ra,dec,dist,sid,ha,alt2,alt3;
	short i = 0;

	/* first guess */
	
	sid=lst(jdguess,longit);
	accumoon(jdguess,lat,sid,0.,&ra,&dec,&dist);
	ha = lst(jdguess,longit) - ra;
	alt2 = altit(dec,ha,lat);
	jdguess = jdguess + del;
	sid = lst(jdguess,longit);
	accumoon(jdguess,lat,sid,0.,&ra,&dec,&dist);
	alt3 = altit(dec,(sid - ra),lat);
	err = alt3 - alt;
	deriv = (alt3 - alt2) / del;
	while((fabs(err) > 0.01) && (i < 10)) {
		jdguess = jdguess - err/deriv;
		sid=lst(jdguess,longit);
		accumoon(jdguess,lat,sid,0.,&ra,&dec,&dist);
		alt3 = altit(dec,(sid - ra),lat);
		err = alt3 - alt;
		i++;
		if(i == 9) return (-1.0e10); /* bad status flag */
	}	
	jdout = jdguess;
	return(jdout);
}

double jd_sun_alt(alt,jdguess,lat,longit)

double alt, jdguess, lat, longit;

{
	/* returns jd at which sun is at a given 
		altitude, given jdguess as a starting point. */

	double jdout;
	double deriv, err, del = 0.002;
	double ra,dec,ha,alt2,alt3;
	short i = 0;

	/* first guess */
	
	lpsun(jdguess,&ra,&dec);
	ha = lst(jdguess,longit) - ra;
	alt2 = altit(dec,ha,lat);
	jdguess = jdguess + del;
	lpsun(jdguess,&ra,&dec);
	alt3 = altit(dec,(lst(jdguess,longit) - ra),lat);
	err = alt3 - alt;
	deriv = (alt3 - alt2) / del;
	while((fabs(err) > 0.1) && (i < 10)) {
		jdguess = jdguess - err/deriv;
		lpsun(jdguess,&ra,&dec);
		alt3 = altit(dec,(lst(jdguess,longit) - ra),lat);
		err = alt3 - alt;
		i++;
		if(i == 9) return(-1.0e10); /* bad status flag */
	}	
	jdout = jdguess;
	return(jdout);
}

void find_dst_bounds(yr,stdz,use_dst,jdb,jde)

short yr, use_dst;
double 	stdz,*jdb,*jde;

{
	/* finds jd's at which daylight savings time begins 
	    and ends.  The parameter use_dst allows for a number
            of conventions, namely:
		0 = don't use it at all (standard time all the time)
		1 = use USA post-1986 convention (1st Sun in April to
                     last Sun in Oct)
		2 = use Spanish convention (for Canary Islands)
		-1 = use Chilean convention (CTIO).
                -2 = use Australian convention.
	    Negative numbers denote sites in the southern hemisphere,
            where jdb and jde are beginning and end of STANDARD time for
            the year. */

	struct date_time trial;

	if((use_dst == 1) || (use_dst == 0)) { 
	    /* USA Convention, and including no DST to be defensive */
		trial.y = yr;
		trial.mo = 4;
		trial.d = 1; 
		trial.h = 2;
		trial.mn = 0;
		trial.s = 0;

		while(day_of_week(date_to_jd(trial)) != 6) {
			trial.d++;
		}
		*jdb = date_to_jd(trial) + stdz/24.;	
		trial.mo = 10;
		trial.d = 31;
		while(day_of_week(date_to_jd(trial)) != 6) {
			trial.d--;
		}
		*jde = date_to_jd(trial) + (stdz - 1.)/24.; 		
	}
	else if (use_dst == 2) {  /* Spanish, for Canaries */
		trial.y = yr;
		trial.mo = 3;
		trial.d = 31; 
		trial.h = 2;
		trial.mn = 0;
		trial.s = 0;

		while(day_of_week(date_to_jd(trial)) != 6) {
			trial.d--;
		}
		*jdb = date_to_jd(trial) + stdz/24.;	
		trial.mo = 9;
		trial.d = 30;
		while(day_of_week(date_to_jd(trial)) != 6) {
			trial.d--;
		}
		*jde = date_to_jd(trial) + (stdz - 1.)/24.; 		
	}		
	else if (use_dst == -1) {  /* Chilean, for CTIO, etc.  */
	   /* off daylight 2nd Sun in March, onto daylight 2nd Sun in October */
		trial.y = yr;
		trial.mo = 3;
		trial.d = 8;  /* earliest possible 2nd Sunday */
		trial.h = 2;
		trial.mn = 0;
		trial.s = 0;

		while(day_of_week(date_to_jd(trial)) != 6) {
			trial.d++;
		}
		*jdb = date_to_jd(trial) + (stdz - 1.)/24.;
			/* note jdb is beginning of STANDARD time in south,
				hence use stdz - 1. */	
		trial.mo = 10;
		trial.d = 8;
		while(day_of_week(date_to_jd(trial)) != 6) {
			trial.d++;
		}
		*jde = date_to_jd(trial) + stdz /24.; 		
	}		
        else if (use_dst == -2) {  /* For Anglo-Australian Telescope  */
           /* off daylight 1st Sun in March, onto daylight last Sun in October */
                trial.y = yr;
                trial.mo = 3;
                trial.d = 1;  /* earliest possible 1st Sunday */
                trial.h = 2;
                trial.mn = 0;
                trial.s = 0;

                while(day_of_week(date_to_jd(trial)) != 6) {
                        trial.d++;
                }
                *jdb = date_to_jd(trial) + (stdz - 1.)/24.;
                        /* note jdb is beginning of STANDARD time in south,
                                hence use stdz - 1. */
                trial.mo = 10;
                trial.d = 31;
                while(day_of_week(date_to_jd(trial)) != 6) {
                        trial.d--;
                }
                *jde = date_to_jd(trial) + stdz /24.;
        }
}



double zone(use_dst,stdz,jd,jdb,jde) 

short use_dst;
double stdz,jd,jdb,jde;

{
	/* Returns zone time offset when standard time zone is stdz,
	   when the first daylight time change date (for the year) is 
	   jdb, and the last is jde.  The parameter 
		use_dst = 0  (don't use it)
		use_dst = 1  (USA convention, post-1986)
		use_dst = 2  (Canary Islands convention)
		use_dst = -1 (Chilean convention)
		use_dst > 0  (other Northern Hemisphere sites)
		use_dst < 0  (reserved for Southern hemisphere sites)
	*/

	if(use_dst == 0) return(stdz);
	else if(((jd > jdb) && (jd < jde)) && (use_dst > 0)) return(stdz-1.);
	else if(((jd < jdb) || (jd > jde)) && (use_dst < 0)) return(stdz-1.);
	else return(stdz);
}
	
short get_line(s)
/* gets a line terminated by end-of-line and returns number of characters. */
char s[];

{	
	char c;
	short i = 0;

	c = getchar(); /* get the first character */
	/* chew through until you hit non white space */
	while((c == '\n') || (c == ' ') || (c == '\t')) c = getchar();

	s[i]=c;
	i++;

	/* keep going til the next newline */
	while((c=getchar()) != '\n') {
		s[i]=c;
		i++;
	}
	s[i]='\0';  /* terminate with null */
	return(i);
}

void load_site(longit,lat,stdz,use_dst,zone_name,elev,horiz,site_name)

double *lat, *longit, *stdz, *elev, *horiz;
short *use_dst;
char *site_name;
char *zone_name;


/* Here are all the site-specific quantities; they are:

		- west longitude, in decimal hours
		- north latitude, in decimal degrees
		- standard time zone offset, hours west of greenwich
		- flag for whether daylight time is to be used, 1=yes, 0=no
		- name of time zone
		- name of site
		- elev = site elevation above its horizon
		- horiz (derived from elev) = correction to rise/set zenith dist
	    Be sure to change all of them at once. Watch signs and units.  */

{
	short nch;
	char obs_code[3];  /* need only one char, but why not? */

	printf("*SELECT SITE* - Enter single-character code:\n");
	printf("   n .. new site (enter all parameters).\n");
	printf("   x .. exit without change (current: %s)\n",site_name);
	printf("   k .. Kitt Peak\n");
	printf("   s .. Shattuck Observatory\n");
	printf("   c .. Cambridge, MA - Harvard Coll. Obs.\n");
	printf("   h .. Mt. Hopkins, AZ (MMT, FLWO)\n");
	printf("   p .. Palomar Observatory\n");
	printf("   t .. Tololo (Cerro Tololo Interamerican Obs.)\n");
	printf("   r .. Roque de los Muchachos, La Palma, Canary Is.\n");
	printf("   b .. Black Moshannon Obs., Penn State U.\n");
	printf("   d .. Dominion Astrophysical Obs., Victoria, BC\n");
	printf("   m .. Mauna Kea, Hawaii\n");
	printf("   l .. Lick Observatory\n");
	printf("   Any other char .. OTHER (You'll be prompted for params.)\n");
	printf("Your answer --> ");
	scanf("%s",obs_code);
	if(obs_code[0] == 'x') {
		printf("No action taken.\n");
		return;
	}
	if(obs_code[0] == 'k') {	
		strcpy(site_name,"Kitt Peak");
		strcpy(zone_name,"Mountain");
		*use_dst = 0;
		*longit = 7.44111; /* decimal hours */
		*lat = 31.9533;    /* decimal degrees */
		*stdz = 7.;
		*elev = 500.;  /* approximate height above horizon */
	}
	else if (obs_code[0] == 's') {
		strcpy(site_name,"Shattuck Observatory");
		strcpy(zone_name,"Eastern");
		*use_dst = 1;
                *longit = 4.81907;  /* from GPS */
                *lat = 43.705;
		*stdz = 5.;
		*elev = 0.; /* no clear horizon */
	}
	else if (obs_code[0] == 'p') {
		strcpy(site_name,"Palomar Observatory");
		strcpy(zone_name,"Pacific");
		*use_dst = 1;
		*longit = 7.79089;
		*lat = 33.35667;
		*elev = 0.;  /* lots of trees, I think ... */
		*stdz = 8.;
	}
	else if (obs_code[0] == 't') {
		strcpy(site_name,"Cerro Tololo");
		strcpy(zone_name,"Chilean");
		*use_dst = -1; 
		*longit = 4.721;
		*lat = -30.165;
		*elev = 2215.;  /* for ocean horizon, not Andes! */
		*stdz = 4.;
	}
	else if (obs_code[0] == 'h') {
		strcpy(site_name,"Mount Hopkins, Arizona");
		strcpy(zone_name,"Mountain");
		*use_dst = 0;
		*longit = 7.39233;
		*lat = 31.6883;
		*elev = 500.;  /* approximate elevation above horizon mtns */
		*stdz = 7.;
	}
	else if (obs_code[0] == 'c') {
		strcpy(site_name,"Harvard College Observatory");
		strcpy(zone_name,"Eastern");
		*use_dst = 1;
		*longit = 4.742;
		*lat = 42.38;
		*elev = 0.;  /* city */
		*stdz = 5.;
	}
	else if (obs_code[0] == 'b') {
		strcpy(site_name,"Black Moshannon Observatory");
		strcpy(zone_name,"Eastern");
		*use_dst = 1;
		*longit = 5.20033;
		*lat = 40.92167;
		*elev = 0.;  /* not set ... who knows? */
		*stdz = 5.;
	}
	else if (obs_code[0] == 'd') {
		strcpy(site_name,"DAO, Victoria, BC");
		strcpy(zone_name,"Pacific");
		*use_dst = 1;
		printf("\n\nWARNING: United States conventions for DST assumed.\n\n");
		*longit = 8.22778;
		*lat = 48.52;
		*elev = 74.;  /* not that it makes much difference .. */
		*stdz = 8.;
	}
	else if (obs_code[0] == 'm') {
		strcpy(site_name,"Mauna Kea, Hawaii");
		strcpy(zone_name,"Hawaiian");
		*use_dst = 0;
		*longit = 10.36478;
		*lat = 19.8267;
		*elev = 4215.;  /* gasp!  Ocean horizon ... */
		*stdz = 10.;
	}
	else if (obs_code[0] == 'l') {
		strcpy(site_name,"Lick Observatory");
		strcpy(zone_name,"Pacific");
		*use_dst = 1;
		*longit = 8.10911;
		*lat = 37.3433;
		*elev = 1290.;  /* for Pacific sunsets, not Sierra sunrise! */
		*stdz = 8.;
	}
	else if (obs_code[0] == 'r') {
		strcpy(site_name,"Roque de los Muchachos Obs.");
		strcpy(zone_name,"pseudo-Greenwich");
		*use_dst = 2;
		*longit = 1.192;
		*lat = 28.75833;
		*elev = 2326.;  /* ocean horizons? */
		*stdz = 0;
	}
	else {	
		printf("West longitude, (h m s); current value ");
		put_coords(*longit,3);
		printf(": ");
		*longit=get_coord();
		printf("Latitude, (d m s); current value ");
		put_coords(*lat,2);
		printf(": ");
		*lat=get_coord();
		printf("Elevation above horizon, meters (for rise/set times):");
		printf("Current value = %5.0f:",*elev);
		scanf("%lf",elev);
		printf("\nStd time zone, hours W; currently %3.0f :",*stdz);
		scanf("%lf",stdz);
		printf("Site name (< 30 char): ");
		nch=get_line(site_name);
		printf("Time zone name, e. g., Central: ");
		nch = get_line(zone_name);
		printf("Enter daylight time convention:\n");
		printf("  0 ... do not use daylight time.\n");
		printf("  1 ... USA, 1986 and later.\n");
		printf("  2 ... Spanish convention (Canaries)\n");
		printf(" -1 ... Chilean convention.\n");
		printf("  your answer ---> ");
		scanf("%hd",use_dst);
	}
	/* now compute derived quantity "horiz" = additional depression
			 of horizon */
	*horiz = pow((2. * *elev / 6378140.),0.5) * DEG_IN_RADIAN;
}

void flmoon(n,nph,jdout) 

/* Gives jd (+- 2 min) of phase nph on lunation n; replaces
less accurate Numerical Recipes routine.  This routine 
implements formulae found in Jean Meeus' *Astronomical Formulae
for Calculators*, 2nd edition, Willman-Bell.  A very useful
book!! */

int n, nph;  /* lunation and phase; nph = 0 new, 1 1st, 2 full, 3 last */
double *jdout;  /* jd of requested phase */


{
	double jd, cor;
	double M, Mpr, F;
	double T;
	double lun;

	lun = (double) n + (double) nph / 4.;
	T = lun / 1236.85;
	jd = 2415020.75933 + 29.53058868 * lun 	
		+ 0.0001178 * T * T 
		- 0.000000155 * T * T * T
		+ 0.00033 * sin((166.56 + 132.87 * T - 0.009173 * T * T)/DEG_IN_RADIAN);
	M = 359.2242 + 29.10535608 * lun - 0.0000333 * T * T - 0.00000347 * T * T * T;
	M = M / DEG_IN_RADIAN;
	Mpr = 306.0253 + 385.81691806 * lun + 0.0107306 * T * T + 0.00001236 * T * T * T;
	Mpr = Mpr / DEG_IN_RADIAN;
	F = 21.2964 + 390.67050646 * lun - 0.0016528 * T * T - 0.00000239 * T * T * T;
	F = F / DEG_IN_RADIAN;
	if((nph == 0) || (nph == 2)) {/* new or full */
		cor =   (0.1734 - 0.000393*T) * sin(M)
			+ 0.0021 * sin(2*M)
			- 0.4068 * sin(Mpr)
			+ 0.0161 * sin(2*Mpr)
			- 0.0004 * sin(3*Mpr)
			+ 0.0104 * sin(2*F)
			- 0.0051 * sin(M + Mpr)
			- 0.0074 * sin(M - Mpr)
			+ 0.0004 * sin(2*F+M)
			- 0.0004 * sin(2*F-M)
			- 0.0006 * sin(2*F+Mpr)
			+ 0.0010 * sin(2*F-Mpr)
			+ 0.0005 * sin(M+2*Mpr);
		jd = jd + cor;
	}
	else {
		cor = (0.1721 - 0.0004*T) * sin(M)
			+ 0.0021 * sin(2 * M)
			- 0.6280 * sin(Mpr)
			+ 0.0089 * sin(2 * Mpr)
			- 0.0004 * sin(3 * Mpr)
			+ 0.0079 * sin(2*F)
			- 0.0119 * sin(M + Mpr)
			- 0.0047 * sin(M - Mpr)
			+ 0.0003 * sin(2 * F + M)
			- 0.0004 * sin(2 * F - M)
			- 0.0006 * sin(2 * F + Mpr)
			+ 0.0021 * sin(2 * F - Mpr)
			+ 0.0003 * sin(M + 2 * Mpr)
			+ 0.0004 * sin(M - 2 * Mpr)
			- 0.0003 * sin(2*M + Mpr);
		if(nph == 1) cor = cor + 0.0028 - 
				0.0004 * cos(M) + 0.0003 * cos(Mpr);
		if(nph == 3) cor = cor - 0.0028 +
				0.0004 * cos(M) - 0.0003 * cos(Mpr);
		jd = jd + cor;

	}
	*jdout = jd;
}

void setupTeX(option) 

short option;   /* 1 = 1 month per page, 2 = 2 months per page, 3 = landscape */

{

  printf("\n\n%%  ---- CUT HERE --- delete this line and everything above it.---- \n");
  printf("%% see the TeXBook, D. Knuth, Addison-Wesley, ISBN 0-201-13448-9, p.382\n");
  printf("%% The following sizing parameters may need to be tweaked for your setup:\n");
  if(option < 3) {
     printf("\\magnification=835\n");
     printf("\\hsize 7.6truein\n");
     printf("\\hoffset -0.7truein \n");
     if(option == 1) 
        printf("\\baselineskip=9.8pt \n");
     else {
        printf("\\baselineskip=9.8pt \n");
        printf("\\voffset -0.55truein\n");
        printf("\\vsize 10.0truein\n");
     }
  }
  else {
        printf("\\magnification=1040\n");
        printf("\\hsize 9.5truein\n");
        printf("\\baselineskip=10.2pt \n");
        printf("\\voffset -0.3truein\n");
        printf("\\vsize 7.2truein\n");
  }
  printf("%% The rest of this should be system-independent:\n");
  printf("\\nopagenumbers\n");
  printf("\\def\\uncatcodespecials{\\def\\do##1{\\catcode`##1=12 }\\dospecials} \n");
  printf("\\def\\doverbatim#1{\\def\\next##1#1{##1\\endgroup}\\next} \n");
  printf("\\def\\setupverbatim{\\tt \n");
  printf("  \\def\\par{\\leavevmode\\endgraf} \\catcode`\\`=\\active \n");
  printf("  \\obeylines \\uncatcodespecials \\obeyspaces} \n");
  printf("  {\\obeyspaces\\global\\let =\\ }\n");
  printf("\\def\\listing#1{\\par\\begingroup\\setupverbatim\\input#1 \\endgroup} \n");
  printf("\\def\\verbatim{\\begingroup\\setupverbatim\\doverbatim} \n");
  printf("\\verbatim$\n");
}

void page_top(site_name,longit,lat,zone_name,use_dst,elev,stdz) 

char *site_name, *zone_name;
double longit, lat, elev, stdz;
short use_dst;

{

	printf("Calendar for %s, west longitude (h.m.s) = ",site_name);
	put_coords(longit,2);
	printf(", latitude (d.m) = ");
	put_coords(lat,1);
	printf("\n");
	printf("Note that each line lists events of one night,");
	printf(" spanning two calendar dates.  Rise/set times are given\n");
	printf("in %s time (%3.0f hr W),",zone_name,stdz);
	if(elev == 0.) printf(" uncorrected for elevation, ");
	else printf(" for %4.0f m above surroundings, ",elev);
	if(use_dst == 0) printf("in standard time all year.\n");
	else printf("DAYLIGHT time used, * shows night clocks are reset.\n");
	printf("Moon coords. and illum. are for local midnight,");
	printf(" even if moon is down.  Program: John Thorstensen, Dartmouth College.\n");
}

void month_banner(y,m) 

short y, m;

{
	char moname[20];
	
	switch(m) {
	  case 1:  strcpy(moname,"JANUARY");
	           break;
	  case 2:  strcpy(moname,"FEBRUARY");
	           break;
	  case 3:  strcpy(moname,"MARCH");
	           break;
	  case 4:  strcpy(moname,"APRIL");
	           break;
	  case 5:  strcpy(moname,"MAY");
	           break;
	  case 6:  strcpy(moname,"JUNE");
	           break;
	  case 7:  strcpy(moname,"JULY");
	           break;
	  case 8:  strcpy(moname,"AUGUST");
	           break;
	  case 9:  strcpy(moname,"SEPTEMBER");
	           break;
	  case 10:  strcpy(moname,"OCTOBER");
	           break;
	  case 11:  strcpy(moname,"NOVEMBER");
	           break;
	  case 12:  strcpy(moname,"DECEMBER");
	           break;
	  default: strcpy(moname,"Month err,not 1-12!");
        }
	printf("                                           ***** %4d %s *****\n",
	   y,moname);
}

void column_heads(jdroot,year) 

long jdroot;
short year;

{
	printf("  Date (eve/morn)      JDmid    LMSTmidn   ---------- Sun: --------- ");
	printf("  LST twilight:  ------------- Moon: --------------\n");
	printf("  (%4d at start)    (-%7d)       ",year,jdroot);
	printf("     set  twi.end twi.beg rise");
	printf("    eve    morn    rise   set  %%illum   RA      Dec\n\n");

}

void info_page(year,site_name,use_dst)

short year, use_dst;
char *site_name;

{
	printf("\n         ***** %d Night-time Astronomical Calendar for %s ***** \n\n",year,site_name);
	printf("               By John Thorstensen, Dartmouth College\n\n");

	printf("   This calendar is designed to provide information useful for ");
	printf("the planning of nighttime observations.\n");
	printf("The format should minimize confusion;  each line gives the ");
	printf("phenomena for a single (local!) night,\n");
	printf("and each line is");
	printf(" labeled with both evening and morning (local) day and date.\n");
	printf("Note that all times given are LOCAL CIVIL (zone) times.  ");
	if(use_dst == 1) {
		printf("DAYLIGHT SAVINGS time is used from the\nfirst Sunday ");
		printf("in April to last in October; this is the present");
		printf("(1986+) convention in the U.S.A.");
	}
	else if (use_dst != 0) 
	  printf("DAYLIGHT SAVINGS time is used, with\na non-USA convention.\n");
	printf("\n\n");  
	printf("   The rise/set times printed are the times at which the center of ");
	printf("the object is 50 arcminutes below\nthe geometrical horizon.");
	printf("  At the given twilight, the center of the sun is %5.1f degrees");
	printf(" below the geometrical horizon.\n\n",(-1. * TWILIGHT_ALT));
	printf("   The moon positions (and rise/set times) are generated by");
	printf(" an implementation of the Low-Precision formulae");
	printf("\nin the Astronomical Almanac.  The Almanac states that");
	printf(" the error seldom exceeds 0.3 degrees.");
	printf("\nTopocentric corrections are included.  Comparisons with ");
	printf("tables for Kitt Peak in the NOAO Newsletter");
	printf("\nindicate that the rise-set times are good to +- 2 min");
	printf(" or so.  The moon's RA, Dec, and illuminated fraction");
	printf("\nare given for local midnight, regardless of whether the");
	printf(" moon is actually up at that time.");
	printf("\nNote that the moonrise and moonset times are not printed ");
	printf("if they occur near mid-day.");
	printf("\n\n   The LST at evening and morning twilight are");
	printf(" tabulated.  This gives an accurate idea of the");
	printf("\nrange of RA's accessible during the night.");
	printf("\n\n   The JD is given (severely rounded off) for local");
	printf(" midnight.  Again, this avoids any ambiguity.\n");
	printf("\n   Some credits:  The sidereal time and Julian date");
	printf(" routines were originally coded in PL/I by");
	printf("\nSteve Maker of Dartmouth College.  The algorithms");
	printf(" originated in the old American Ephemeris.");
	printf("\nThe routine to convert JD back to calendar date is");
	printf(" adapted from Numerical Recipes in C, by Press et al.");
	printf("\n\n    CAUTIONS:  I believe that the program which ");
	printf("generates these tables is reasonably accurate.  However,");
	printf("\nit has not been exhaustively tested, so you should be");
	printf(" sure to run 'sanity checks' on the results.  Also,");
	printf("\nin view of the approximations used, the results should");
	printf(" not be used when high precision is needed.");
	printf("\nExtension to dates far from the present (1990) ");
	printf("should be done with great caution.  The code has");
	printf("\nnot been tested for the eastern ");
	printf("or southern hemishpheres.  Rise/set times are slightly ");
	printf("incaccurate and\nrather confusing at circumpolar ");
	printf("latitudes, where the concept of a 'night' is blurry.\n");
	printf("\nThe daylight savings time conventions (if used) are quite");
	printf(" specific (to U. S., post-1986) and subject to change.");
	printf("\nI know that the code has many infelicities; ");
	printf("if you should find actual");
	printf(" errors,");
	printf("please notify\n");
	printf("   John.Thorstensen@dartmouth.edu\n\n\n");
	
	printf("  [This output comes from a (hopefully) portable, completely ");
	printf("self-contained program in the c language.  It is available\n");
	printf("from the author and may be used freely for scientific or ");
	printf("educational purposes.  If you use it for profit, \n");
	printf("please contact the author to arrange a (modest!) fee.");
	printf("\nSource code is copyright John Thorstensen, 1990.]\n");


}

main()

{
	struct date_time date,dateback, date1,date2;
	double jd, jdmid, jdeve, jdmorn, jdsunset, jdsunrise;
	double jdmoonrise,jdmoonset,jdetwilight,jdmtwilight,jd1,jd2,jdout;
	double jdbdst, jdedst;  /* jd at beg and end of dst this year */
	   /* for southern sites, interpretation of these numbers reversed. */
	double sid;
	double stmid;
	struct coord sidereal;
	short option;
	short TeX_out = 0;
	short done = 0, optdone = 0;
	double ramoon, decmoon, distmoon;
	double hamoonset, hasunset, hatwilight, tmoonset, tmoonrise;
	double rasun, decsun;
	float ill_frac, temp;
	short day, i, j, k, nights, jdtr;
	long jdroot;   /* for banner */
	short year, daymo[13];
	
	char ff = 12;  /* ascii form feed for pagination */
	char *months = "JanFebMarAprMayJunJulAugSepOctNovDec";
	char mo_out[4];
	char site_code[3];  /* only using one character, actually */
	short nchars, tries;

	int lunation, mphase;  /* for lunar calendar */
	struct date_time scr_date;
	double  jdjan0, jddec32;


/* Here are all the site-specific quantities; they are:

		- west longitude, in decimal hours
		- north latitude, in decimal degrees
		- standard time zone offset, hours west of greenwich
		- flag for whether daylight time is to be used, 1=yes, 0=no
		- hours from midnight to print moon rise or set (8 is good
                    for temperate latitudes; more at very high lat.)	
		- name of time zone
		- elevation of site above effective horizon, for rise/set
		- zenith angle correction due to elevation, dec. degr.
		- name of site
	    Be sure to change all of them at once. Watch signs and units.  */

/* Example .. for Kitt Peak, Arizona, MDM Observatory. */

	double longit = 7.44111; 
	double lat = 31.9533;    
	double stdz = 7.;       
	short use_dst = 0;
	double moon_pr = 7.5;     
	double elev;
	double horiz;
/* .... These are initialized later with strcpy for portability..... */
	char zone_name[25];
	char site_name[40];

	strcpy(site_name,"Kitt Peak");
	strcpy(zone_name,"MOUNTAIN STANDARD");

/* .. end of site-specific quantities.  To install a new 'canned' site,
      read and change routine 'load_site' above.  It's not hard.  */

/* also load contents of daymo[] explicitly here, to keep some
    compilers happy... */
	daymo[0] = 0;
	daymo[1] = 31;
	daymo[2] = 28;
	daymo[3] = 31;
	daymo[4] = 30;
	daymo[5] = 31;
	daymo[6] = 30;
	daymo[7] = 31;
	daymo[8] = 31;
	daymo[9] = 30;
	daymo[10] = 31;
	daymo[11] = 30;
	daymo[12] = 31;

/*	Introductory banner for user, with site information and 
	opportunity to escape foreground execution. */

	printf("Nighttime astronomical calendar program.\n");
	printf("Select a site: \n");

	load_site(&longit,&lat,&stdz,&use_dst,zone_name,&elev,
			&horiz,site_name);

	printf("The site you've selected is:  %s.\n",
			site_name);
	if(use_dst != 0) {
		if((lat < 0.) && (use_dst > 0)) {
			printf("You've done something ridiculous, namely\n");
			printf("selected DST (computed for NORTHERN (USA) summer)\n");
			printf("at a SOUTHERN site.  If you need DST, change\n");
			printf("the source code and change this error trap.\n");
			printf("Goodbye.\n");
			goto BLUNDER;
		}		
		printf("You've selected daylight savings time. If your site\n");
		printf("is in US it should switch at 2AM on 1st Sun in April and last Sun in October.\n");
		printf("This was mandated in 1986 for the U. S.; it was different before and\n");
		printf("always could change at a later date, and may not apply\n");
		printf("in other countries.  Foreign conventions are coded and\n");
		printf("added as needed, e. g. Spanish (for Canaries.), Chilean.\n");
	}	
	printf("Type 0 for ordinary text,\n");
	printf("1 for TeX-style output with one month per page, or\n");
	printf("2 for TeX-style output with two months per page, or\n");
	printf("3 for landscape TeX-style output with one month per page:\n");
	scanf("%hd",&TeX_out);
	if(TeX_out != 0) 
              printf("\nYou'll have to edit the TeX output ... look for 'CUT HERE'\n");
	printf("Run this code in background and capture output.\n");
        printf("See source code or documentation\n");
	printf("for details of what program will expect.\n");
	printf("To exit now, give a negative year at the next prompt.\n");
	printf("Year to print, negative to abort: ");
	scanf("%hd",&year);

	if(year < 0) goto BLUNDER;
	if(year < 100) {
		year = year + 1900;
		printf("Your answer taken to be %d ! \n",year);
	}
	if((year < 1901) || (year > 2099)) {  /* filter input a bit more */
		printf("Year %d out of range -- exiting.\n",year);
		goto BLUNDER;
	}
	if(use_dst != 0) find_dst_bounds(year,stdz,use_dst,&jdbdst,&jdedst);
 

	temp = ((float) year) / 4.;
	if(fabs(temp - (int) temp ) < 0.01) daymo[2] = 29; /* leapyr -
                      ok for 2000, which is a leap year */ 
	date1.y = year;
	date1.mo = 1;
	date1.d = 1;
	date1.h = 18;  /* local afternoon */
	date1.mn = 0;
	date1.s = 0;  /* afternoon */
	date2.y = year;
	date2.mo = 12;
	date2.d = 31;
	date2.h = 18;  /* local afternoon */
	date2.mn = 0;
	date2.s = 0;  /* afternoon */
	jd1 = date_to_jd(date1);
	jd2 = date_to_jd(date2);
	nights = (short)(jd2 - jd1);
	if((jd2-jd1) > 500) {
		printf("excessively long time span - exiting.\n");
		goto BLUNDER;
	}
	if(jd1 < 2400000.) {
		printf("Dates too early.  Use all four digits for year.\n");
		goto BLUNDER;
	}


	/* print some info about the program */

	if(TeX_out == 0) printf("%c",ff); /* form feed */
	else setupTeX(TeX_out);

	info_page(year,site_name,use_dst);  /* "throwaway page" of info */

	if(TeX_out == 0) printf("%c",ff); /* separate page for lunar calendar */
	else {
		printf("$\n");
		printf("\\par\\vfill\\eject\n");
		printf("\\verbatim$\n");
	}
	printf("      MOON PHASES FOR %d, at %s\n\n",year,site_name);
	printf("Times and dates are given in local time, zone = %3.0f hr West.\n",
			stdz); 
	printf("They are generally better than +- 2 minutes.\n");
	if(use_dst != 0) printf("Daylight savings time used.\n");
	printf("\n   The end of the previous year and the beginning of the next\n");
	printf("are included for continuity.\n\n");
	
	scr_date.y = year;
	scr_date.mo = 1;
	scr_date.d = 0;
	scr_date.h = 0;
	scr_date.mn = 0;
	scr_date.s = 0;
	jdjan0 = date_to_jd(scr_date);
	scr_date.mo = 12;
	scr_date.d = 32;
	jddec32 = date_to_jd(scr_date);

	lunation = (year - 1900) * 12.3;  /* puts it before beginning */
	jd = 0.;
	/* find first lunation of yr */
	while((jd < jdjan0) && (lunation < 2000)) { 
		flmoon(lunation,0,&jd);
		lunation++;
	}
	if(lunation < 1995) { /* skip lunar calendar if results make no sense */
		lunation = lunation - 2; /* back up */
	printf("    NEW             1ST             FULL            LAST\n\n");
		while(jd < jddec32) {
			for(mphase=0; mphase <= 3; mphase++) {
				flmoon(lunation, mphase, &jd);
				print_calendar((jd -
				    zone(use_dst,stdz,jd,jdbdst,jdedst)/24.),0);
				printf(" ");
				print_time((jd -
 				   zone(use_dst,stdz,jd,jdbdst,jdedst)/24.),0);
				printf("   ");
			}
			printf("\n\n");
			lunation++;
		}
	}
	else printf("Some error in lunations ... no lunar calendar printed.\n");

	for(j = 1; j <= 12; j++) {    /*   unindented loop! */

	if(TeX_out == 0) {
		printf("%c",ff);  /*  form feed for new month*/
	}
	else if(TeX_out == 1 || TeX_out == 3) {
		printf("$\n");
		if(j < 13) {
			printf("\\par\\vfill\\eject\n");
			printf("\\verbatim$\n");
		}
	}
	else if(TeX_out == 2) {
		switch(j) {  /* awkward but it works .... */
			case 1:
			case 3:
			case 5:
			case 7:
			case 9:
			case 11:   		
				printf("$\n");
				printf("\\par\\vfill\\eject\n");
				printf("\\verbatim$\n");
				page_top(site_name,longit,lat,zone_name,
					use_dst,elev,stdz);  
				break;
			default: 
				break;
		}
	}
	printf("\n");

	month_banner(year,j);

	if(TeX_out < 2) {
		printf("\n");
		page_top(site_name,longit,lat,zone_name,use_dst,
				elev,stdz);
	}
	printf("\n");
	date.y = year;
	date.mo = j;
	date.d = 1;  
	date.h = 18;  /* evening of first night of month */
	date.mn = 0;
	date.s = 0;
	jd = date_to_jd(date); /* not really jd; local equivalent */
	jdtr = jd / 10000.;  /* type short -- truncate */
	jdroot = jdtr * 10000; /* type long */
        column_heads(jdroot,date.y);
	for(i=1;i<=daymo[j];i++) {  
		date.d = i;
		jd = date_to_jd(date); /* not really jd; local equivalent */
		if(day_of_week(jd) == 6) printf("\n"); /* blank line at Sunday */
		print_day(day_of_week(jd));
		printf(" ");
		print_calendar(jd,0); /* translate back e.g. 11/31-12/1 */
		printf("/");
		jd = jd + .5;  /* local morning */
		print_day(day_of_week(jd));
		printf(" ");
		print_calendar(jd,0);
		jd = jd - 0.25;  /* local midnight */
		jdmid = jd + zone(use_dst,stdz,jd,jdbdst,jdedst) / 24.;  /* corresponding ut */
		if(use_dst != 0) {
			if((fabs(jdmid-jdbdst)<0.49)||(fabs(jdmid-jdedst)<0.49)) 
				printf("*");  /* print a star if changing times */
			else printf(" ");
		}
		else printf(" ");  
		jdout = jdmid - jdroot;
		printf("%7.1f  ",jdout);	
		stmid = lst(jdmid,longit);
		put_coords(stmid,2);
		printf("  ");
		accumoon(jdmid,lat,stmid,0.,&ramoon,&decmoon,&distmoon);
		lpsun(jdmid,&rasun,&decsun);
		hasunset = ha_alt(decsun,lat,-(0.83+horiz));
/* there follows some awful flow of control.  This is an artifact of the
   history of the program - it started as a program for mid-latitudes, and
   then I changed it to handle polar latitudes reasonably well.  This led
   to a lot of really weird branching.  Sorry.  It appears to work, though. */
		if(hasunset > 900.) {  /* sun isn't going to set; no twilight */
			printf(" .....  .....   ....  ....   .....  .....");
		        moon_pr = 10.;   /* print moonrise/set if within 10 hr of midn. */
			goto MIDNIGHT_SUN; /* horrible flow of control */
		}
		if(hasunset < -900.) {  /* sun ain't gonna rise, but may be twilight */
			printf(" .....");
		        moon_pr = 10.;
			goto TWILIGHT;    /* more horrible flow of control */
		}	
		jdsunset = jdmid + adj_time(rasun+hasunset-stmid)/24.; /* initial guess */
		jdsunset = jd_sun_alt(-(0.83+horiz),
			jdsunset,lat,longit); /* refinement */
		if(jdsunset > 0.) /* convergence */
		  print_time((jdsunset-
		     zone(use_dst,stdz,jdsunset,jdbdst,jdedst)/24.),0);
		else printf(" .....");
	      TWILIGHT: hatwilight = ha_alt(decsun,lat,TWILIGHT_ALT);
		if((hatwilight > 900.) | (hatwilight < -900.)) {
			printf("  .....   ....");
			jdetwilight = -1.0e10;
			jdmtwilight = -1.0e10;
			goto SUNRISE; /* even more horrible */
		}
		jdetwilight = jdmid + adj_time(rasun+hatwilight-stmid)/24.;
		jdetwilight = jd_sun_alt(TWILIGHT_ALT,jdetwilight,lat,longit);
		printf(" ");
		if(jdetwilight > 0.) 
		  print_time((jdetwilight-
			zone(use_dst,stdz,jdetwilight,jdbdst,jdedst)/24.),0);
		else printf(" ....."); /* no convergence */
		jdmtwilight = jdmid + adj_time(rasun-hatwilight-stmid)/24.;
		jdmtwilight = jd_sun_alt(TWILIGHT_ALT,jdmtwilight,lat,longit);
		printf(" ");
		if(jdmtwilight > 0.) 
		  print_time((jdmtwilight-
			zone(use_dst,stdz,jdmtwilight,jdbdst,jdedst)/24.),0);
		else printf(" ....."); /* no convergence */
	      SUNRISE: if(hasunset < -900.) printf(" .....");
	      else {
		jdsunrise = jdmid + adj_time(rasun-hasunset-stmid)/24.;
		jdsunrise = jd_sun_alt(-(0.83+horiz),jdsunrise,lat,longit);
		if(jdsunrise > 0.) 
		  print_time((jdsunrise-
			zone(use_dst,stdz,jdsunrise,jdbdst,jdedst)/24.),0);
		else printf(" .....");
	      }
                if((jdsunrise > 0.) && (jdsunset > 0.)) moon_pr = 
			1.5 + 0.5 * (jdsunrise - jdsunset) * 24.;
	        if(moon_pr < 6.5) moon_pr = 6.5;  		
		printf("  ");
	      if((jdetwilight > 0.) && (jdmtwilight > 0.)) {
		sid = lst(jdetwilight,longit);
		put_coords(sid,0);
		printf(" ");		
		sid = lst(jdmtwilight,longit);
		put_coords(sid,0);
	      }
	      else printf(" .....  .....");
	       MIDNIGHT_SUN: printf("  ");
		hamoonset = ha_alt(decmoon,lat,-(0.83+horiz)); /* rough approx. */
	    if(fabs(hamoonset) < 100.) { /* if it's really gonna set */
		tmoonrise = adj_time(ramoon-hamoonset-stmid);
		tmoonset = adj_time(ramoon+hamoonset-stmid);
		jdmoonrise = jdmid + tmoonrise / 24.;
		jdmoonrise = jd_moon_alt(-(0.83+horiz),jdmoonrise,lat,longit);
	       /* moonrise/set printed only if 
                             within moon_pr hours of local midnight */
		if(fabs((jdmoonrise - jdmid) * 24.) < moon_pr) { 
			/* if convergent, and  more-or-less at night */
			print_time((jdmoonrise-
			   zone(use_dst,stdz,jdmoonrise,jdbdst,jdedst)/24.),0);
			printf(" ");
		}
		else printf(" ..... ");
		jdmoonset = jdmid + tmoonset / 24.;
		jdmoonset = jd_moon_alt(-(0.83+horiz),jdmoonset,lat,longit);
		if(fabs((jdmoonset - jdmid) * 24.) < moon_pr) {
			print_time((jdmoonset-
			  zone(use_dst,stdz,jdmoonset,jdbdst,jdedst)/24.),0);
			printf("  ");
		}
		else printf(" .....  ");
	     }  
	     else {
		printf(" .....  .....  "); /* no rise or set likely */
	     }
		ill_frac=0.5*(1.-cos(subtend(ramoon,decmoon,rasun,decsun)));
		printf("%4.0f ",
			(100.*ill_frac));
		put_coords(ramoon,1);
		printf(" ");
		put_coords(decmoon,0);
		printf("\n");
	   }
	   if((TeX_out != 0) && (j == 12)) {
	      printf("$\n");
	      printf("\\par\\vfill\\supereject\\end\n");
	   }

	   }  /* end of the UNINDENTED LOOP */

	   BLUNDER: ;  /* dummy statement - dump here if bad input. */
}
