# snreg_soft
This will eventually be a complete package for common SNSEARCH tools, updated from when they were used nearly a decade ago. Here's a write of what it currenly does.

Gathering photometry and other measures
Here’s a quick writeup of what I did to measure the surface brightness, magnitudes, and effective radii for the host galaxy data we have so far. This assumes the data (from APO, SOAR, and Gemini) have been fully reduced (bias and flat-field corrected) and are essentially science ready. There’s a couple additional steps for SOAR.

First steps
The package assumes you’re running this under astroconda
some additional packages may need to be installed, which can generally be obtained using pip. This is compiled on macosx, but could be recompiled for any unix/linux system. There source files are provided, and there are some notes in the appropriate makefiles.

Fixing SOAR Data
The Goodman spectrograph in imaging mode (a) is highly vignetted and (b) has an incomplete image header.
1. Fix Goodman header problems, and no WCS
    1. use “fix_goodman_header.py [image_name].fits" to fix header info
    2. use "wcs.py [image_name].fits” to find an astrometric solution, and apply to image. 
2. Create a mask image for source extractor
    1. use "mask_masking.py [flat_field].fits" on a bias-subtracted and combined flat-field image. The output image is called “flagimg.fits”.

Determining zero points
As long as the header WCS is mostly correct, we can use a simple routine to download the relevant part of the USNO catalog, match sources, then compare SExtractor magnitudes with USNO magnitudes to determine the zeropoint offset. This should be true for Gemini and APO data, but SOAR Goodman data needs the fixes above. 

If working on data that requires a source mask, use:
“get_zeropoint.py —flagimg flagimg.fits [image_name].fits”
otherwise, just
“get_zeropoint.py [image_name].fits”

It should return an average difference value, mode difference value (preferred), and the standard deviation. 

The mode is preferred as there tends to be a tail of poor photometry in USNO for the faintest sources. This is seen as a positive tail in the distribution of difference magnitudes. A histogram image titled [image_name]_zmag.png is created to illustrate any zernpoint biases.

Generating a final photometry catalog
The last (or second to last) step is to generate a new source catalog with the appropriate magnitudes, kron radii, and magnitude errors.

the script:
“get_photometry.py —flagimg flagimg.fits —zp zeropoint [image_name].fits”
or without a mask image,
“get_photometry.py —zp zeropoint [image_name].fits”
will do the trick. zernpoint is the mode (or average) zeropoint derived in the last step. You have to specify the value.

A last script:
find_targ.py [image_name].fits x y 
extracts the magnitude and effective radius for the source closest to the specified x and y position. 


