/*
 				check.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SExtractor
*
*	Author:		E.BERTIN, IAP/Leiden
*
*	Contents:	handling of "check-images".
*
*	Last modify:	05/09/97
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*--------------------------------- structures ------------------------------*/
/* Check-image parameters */
typedef struct structcheck
  {
  char		filename[MAXCHAR];	/* ptr to check-image filename */
  FILE		*file;			/* ptr to check-image file structure */
  char		*fitshead;		/* ptr to check-image FITS header */
  int		fitsheadsize;		/* size of check-image FITS header */
  void		*pix;			/* ptr to check-image pixmap */
  int		width, height;		/* size of check-image */
  size_t	npix;			/* number of pixels in check-image */
  int		y;			/* current line in check-image */
  PIXTYPE	overlay;		/* intensity of the overlayed plots */
  PIXTYPE	*line;			/* buffered image line */
  checkenum	type;			/* CHECKIMAGE_TYPE */
  }	checkstruct;

/*------------------------------- functions ---------------------------------*/

checkstruct	*initcheck(picstruct *, char *, checkenum);

void		addcheck(checkstruct *, float *, int,int, int,int, float),
		blankcheck(checkstruct *, PIXTYPE *, int,int,int,int,PIXTYPE),
		endcheck(picstruct *field, checkstruct *),
		writecheck(checkstruct *, PIXTYPE *, int);
