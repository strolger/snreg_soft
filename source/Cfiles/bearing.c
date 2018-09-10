/* BEARING - calculate the shifts and rotations needed to be
              applied to HST cameras using the jitter data.

Original version: Andy Fruchter, August 1998

Revised by Richard Hook, August 1998. 

Revised by Andy Fruchter, August 1998

Added preliminary NICMOS and STIS parameters, Richard Hook, September 1998.

Added STIS NUV and FUV MAMAs, Richard Hook, October 1998
      Note that the signs of the offsets for the STISCCD were also
      changed to fit in with the convention for NICMOS.

Distance subroutine revised by Andy Fruchter, July 1999, to more accurately
      handle small angles.

*****
Note that this script assumes that subsequent Drizzle runs
            use shft_un=input, shft_fr=output and align=center.

*/

#include <stdio.h> 
#include <math.h>
#define SHIFT 2.0362e-3
#define DEGREE (3.1415926535898/180.)
#define HOUR (3.1415926535898/12.)
#define TWOPI 6.283185307179
#define PI 3.1415926535898
#define SQR(A) ((A)*(A))
#define gmax(A,B) ((A)>(B)?(A):(B))

/* Instrument IDs */
#define PC1 0
#define WF2 1
#define WF3 2
#define WF4 3
#define NIC1 4
#define NIC2 5
#define NIC3 6
#define STISCCD 7
#define STISFUV 8
#define STISNUV 9
#define NIN 10

/* Routines for calculating bearing and angular separation */
double bearing(double x,double y,double x_0,double y_0);
double distance(double x,double y,double x_0,double y_0);

/* These numbers refer to the X axis angles from the aperture location
   files on the STScI Web site */

static double orient_inst[NIN] = {134.908,224.388,314.698,45.258,
                                  225.312,224.507,224.851,
                                  45.0559,45.0,45.0};

static double scale_x[NIN] = {0.04557,0.09961,0.09958,0.09965,
                              0.04320,0.0760242,0.203859,
                              0.050763,0.012371,0.012375};

static double scale_y[NIN] = {0.04557,0.09961,0.09958,0.09965,
                              0.043026,0.0753412,0.203113,
                              0.050763,0.012371,0.012375};

static double offset_x[NIN] = {578.085,453.150,266.752,253.831,
                               33.0,20.0,11.0,
                               22.384,527.65,527.65};

static double offset_y[NIN] = {571.715,274.320,250.987,454.497,
                               -29.0,31.0,6.0,
                               3.67,528.675,528.675};

main( argc, argv )
  int     argc;
  char    *argv[];
  {
	int counter,CamNum;
        char image[100], CamName[4];
        double ra,dec,orient,ra_0,dec_0,orient_0,bear,dist;
	double angle_x,bear_x,delta_or,x_dsh,y_dsh,x_sh,y_sh;

        /* Check that a camera has been specified */
        if ( argc < 2 || argc > 4 )
        {
        printf("usage: bearing pc1|wf2|wf3|wf4|nic1|nic2|nic3|stisccd|stisfuv|stisnuv\n");
        exit(2);
        }

        /* Extract the camera name from the command line and set the number */
        sscanf( argv[1], "%s", CamName);
        if (strcmp(CamName,"pc1") == 0) CamNum=0;
        else if (strcmp(CamName,"wf2") == 0) CamNum=1;
        else if (strcmp(CamName,"wf3") == 0) CamNum=2;
        else if (strcmp(CamName,"wf4") == 0) CamNum=3;
        else if (strcmp(CamName,"nic1") == 0) CamNum=4;
        else if (strcmp(CamName,"nic2") == 0) CamNum=5;
        else if (strcmp(CamName,"nic3") == 0) CamNum=6;
        else if (strcmp(CamName,"stisccd") == 0) CamNum=7;
        else if (strcmp(CamName,"stisfuv") == 0) CamNum=8;
        else if (strcmp(CamName,"stisnuv") == 0) CamNum=9;
        else { 
             printf("! Invalid camera name \n"); 
             printf("Options are: pc1|wf2|wf3|wf4|nic1|nic2|nic3|stisccd|stisfuv|stisnuv \n");
             exit(2); }

        /* Read lines from the input until EOF */
	counter = 0;
	while (scanf("%s %lf %lf %lf",image,&ra,&dec,&orient) > 0 ) 
          {
	   counter++;

           /* Convert to radians */
           ra = ra * DEGREE;
           dec = dec * DEGREE;
           orient = orient * DEGREE;

           /* The first line defines the reference */
	   if (counter == 1) 
	     {
	      ra_0 = ra;
	      dec_0 = dec;
	      orient_0 = orient;
	      printf("# Setting origin using %s\n",image);
	     }

              /* Calculate the bearing and offset in celestial coordinates */
	      dist = distance(ra,dec,ra_0,dec_0);
	      bear = bearing(ra,dec,ra_0,dec_0);

              /* Derive relative rotation angles */
              angle_x =  orient_0 + orient_inst[CamNum]*DEGREE;
	      delta_or = orient - orient_0; 

              /* Force angles into the appropriate range */
              while (bear < 0) bear = bear + TWOPI;
              while (bear > TWOPI) bear  = bear  - TWOPI;

              while (angle_x < 0) angle_x = angle_x + TWOPI;
              while (angle_x > TWOPI) angle_x = angle_x - TWOPI;

              while (delta_or < (-1*PI)) delta_or = delta_or + TWOPI;
              while (delta_or > PI) delta_or = delta_or - TWOPI;
	  
              bear_x = bear - angle_x;

              /* Apply correction to the centre of the chip as opposed
                 to the reference point */
              rotate_detector(delta_or,&x_dsh,&y_dsh,
                              offset_x[CamNum],offset_y[CamNum]);

              /* Convert from radians to pixels */
              x_sh = cos(bear_x)*dist*3600.0/scale_x[CamNum]/DEGREE;
              y_sh = sin(bear_x)*dist*3600.0/scale_y[CamNum]/DEGREE;

              /* Apply the corrections; sign conversion due to convention */
              x_sh = -x_sh + x_dsh;
              y_sh = -y_sh + y_dsh;

              /* Write out the shifts etc ready for Drizzle */
              printf (" %s %d %.10f %.10f %.10f \n",
                      image,CamNum+1,x_sh,y_sh,delta_or/DEGREE);

              }
	  }

/* Calculate the bearing in a general way from one point to another on
   the sky */
double bearing (double ra_1,double dec_1,double ra_2,double dec_2)
{
	double x,y,d_ra,bear;
	d_ra = ra_2 - ra_1;
	y = sin(d_ra)*cos(dec_2);
	x = sin(dec_2)*cos(dec_1)-cos(dec_2)*sin(dec_1)*cos(d_ra);
	if ((x != 0.0) || (y != 0.0)) 
	   bear = atan2(y,x);
	else
	    bear = 0.0;
        while (bear < 0)
          bear = bear + TWOPI;
	return(bear);
}

/* Calculate the separation in a general way from one point to another on
   the sky */
double distance (double ra_1,double dec_1,double ra_2,double dec_2)
{
	double x_1,y_1,z_1,x_2,y_2,z_2,dist,cross,w;

	x_1 = cos(ra_1)*cos(dec_1);
	y_1 = sin(ra_1)*cos(dec_1);
	z_1 = sin(dec_1);

	x_2 = cos(ra_2)*cos(dec_2);
	y_2 = sin(ra_2)*cos(dec_2);
	z_2 = sin(dec_2);

        cross = x_1*x_2 + y_1*y_2 + z_1*z_2;
	dist = acos(cross);

	/* Modulus squared of half the difference vector */
	w = SQR(x_1 - x_2) + SQR(y_1 - y_2) + SQR(z_1 - z_2);
	w /= 4.0;

	/* Angle between the vectors */
	/* This trick for an accurate conversion from P.T. Wallace's Starlink */

        return 2.0 * atan2 ( sqrt ( w ), sqrt ( gmax ( 0.0, 1.0 - w )));

    
}

/* Calculate the corrections in the offsets from the fact that the
   reference position is not at the centre of the detectors (which is
   where drizzle uses as rotational reference position) */
rotate_detector (double angle, double * x_dsh, double * y_dsh,
                 double xoff, double yoff)
{
        *x_dsh = cos(angle)*xoff - sin(angle)*yoff
                 - xoff;
        *y_dsh = sin(angle)*xoff + cos(angle)*yoff          
                 - yoff;
}
