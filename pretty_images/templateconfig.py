#!/usr/bin/env python
"""
Input Configurations for plotting.py which pretty plots continuum and moment maps
"""

# overwrite parameter 
# if true will overwrite output files
# else will just add an iterator to make non destructive
overwrite = False

# Filename of the input continuum to overlay (must not be blank)
# This also decides how many images to make (i.e. this is the iterator for the rest of the inputs)
# also input the sigma (RMS) level of the image
inputContinuum = ("")
continuumParam = (0.)

# boolean, True if show contours, false otherwise
# red blue respectively
showContours = ((True,True))

# filenames of the input moment maps
# red and blue respectively
inputMoments=(("",""))

# defines the contour level (RMS), the starting contourlvl and the interval of contour levels
redMomentParam = ((0.,0.,0))
blueMomentParam = ((0.,0.,0))

# ouflow parameters
# define if to show outflow vectors and if true give the position angles
# Allow for defining of outflows both red and blue respectively
showOutflows=(False)
outflowParams = ((0.,0.))

# center of the output image (must be within input image range)
# form can be hh:mm:ss or deg.deg values
centerRaDec = (("",""))

# lower and upper bounds (respectively) for the flux
# This just controls the output
# suggest setting lower to noise and upper to 80% max
boundFlux = ((0.,0.))

# image size (arcseconds)
# must be within input images at the centerRaDec +- imsize defined
imSize = ((0.,0.))

# scaling to provide for the image
# supports sqrt, log, none, linear
imScale = ("")

# size of the scalebar (arcseconds)
scaleBar = (0.)

# Distance to the source in parsec
sourceDist = (0.)

# first and secondary plotting labels
plotLabels = (("",""))

# boolean expression string (True,False,Both)
# output Color (True False) or (Both) color and grayscale
# if color is being output, can specify color map from file
outputColor = ("")
colorMap = ("")

# sub image to a smaller region
# just redefine a center and a size with a new scale bar
# everything else stays the same
# center,"(imsize)",bar scale size
subImage = (("","()",""))
