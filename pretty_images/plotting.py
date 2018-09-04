#!/usr/bin/env python

# import modules
import colormap as colorMap
import aplpy
import numpy as np
import matplotlib.pyplot as mpl

# look at bottom for main calls and where to set stuff

def parseRaDec(inputString):
    """
    """
    # parse RA (03:25:36.3005 +030.45.14.4112) to deg
    ra,dec = inputString.split(" ")
    hh,mm,ss = map(float,ra.split(":"))
    deg = (hh+(mm/60.)+(ss/3600.))*15.
    fin = "{}".format(deg)
    dd,mm,ss = map(float,dec.split(".")[0:3])
    ss += float("0." + dec.split(".")[3])
    deg = (dd+(mm/60.)+(ss/3600.))
    fin += " {}".format(deg)
    return fin

class mappingContours():
    """
    """
    def __init__(self,cont,ra,dec,minpixval=0.,maxpixval=0.,\
                size=(0.0),scalebar=0.,distance=0.,name='',imagestretch='sqrt',\
                colororgray='false',colormap='ds9',plotlabel='',misc='',\
                textcolor='',extension=('pdf',),showLabels=False):
        self.cont = cont
        self.ra = ra
        self.dec = dec
        self.minpixval = minpixval
        self.maxpixval = maxpixval
        self.xsize,self.ysize = size
        print("Size:",10,int(10*self.ysize/self.xsize))
        self.fig = mpl.figure(figsize=(10,int(10*self.ysize/self.xsize)))
        self.scalebar = scalebar
        self.distance = distance
        self.name = name
        self.imagestretch = imagestretch
        self.colororgray = colororgray
        if 'ds9' in cmap:
            colorMap.main('ds9')
            self.colormap = colormap
        self.plotlabel = plotlabel
        if textcolor == '':
            self.auto = True
            if colororgray == 'true':
                textcolor = 'white'
            else:
                textcolor = 'black'
        else:
            self.auto = False

        self.textcolor = textcolor
        self.extension = extension
        self.misc = misc
        self.showLabels = showLabels

    def drawContinuum(self):
        """
        """
        if self.cont != "":
            self.gc1 = aplpy.FITSFigure(self.cont, figure=self.fig)
            if (self.colororgray.lower() == 'true'):
                self.gc1.show_colorscale(vmin=self.minpixval,vmax=self.maxpixval,\
                    stretch=self.imagestretch,cmap=self.colormap) 
            if (self.colororgray.lower() == 'false'):
                self.gc1.show_grayscale(vmin=self.minpixval,vmax=self.maxpixval,\
                    stretch=self.imagestretch,invert=True)

    def addBeam(self,imageName):
        """
        """
        dele = mpl.figure(figsize=(1,1))
        self.ContBeam = aplpy.FITSFigure(imageName, figure=dele)
        self.ContBeam.add_beam()   
        dele = None


    def setupContours(self,color,setup,show,imageName,neg,contourStart,contourInterval,contourNoise,label):
        """
        """
        if setup:
            self.addBeam(imageName)
            self.gc1.add_beam()
            for x in self.ContBeam.beam._base_settings:
                self.gc1.beam._base_settings[x] = self.ContBeam.beam._base_settings[x]
        if show:
            contours = self.drawContours(neg,contourStart,contourInterval,contourNoise)
            self.gc1.show_contour(imageName,levels=contours,colors=color,linewidths=1.0)
            self.addLabel(color," contour: "+label)


    def drawContours(self,neg=False,contstart=None,continterval=None,contnoise=None):
        """
        """
        contours1=np.arange(contstart,contstart*10.0,continterval,dtype='float32')
        contours2=np.arange(contstart*10.0,contstart*40.0,continterval*3.0,dtype='float32')
        contours3=np.arange(contstart*40.0,contstart*100.0,continterval*10.0,dtype='float32')
        poscontours=np.concatenate((contours1,contours2,contours3))
        negcontours=poscontours[::-1]*(-1.0)
        contours = poscontours * contnoise
        if neg:
            contours=np.concatenate((poscontours,negcontours))*contnoise
        return contours

    def plotFormat(self):
        """
        """
        # centering
        self.gc1.recenter(self.ra,self.dec,width=(self.xsize/3600.0),height=(self.ysize/3600.0))

        # axis labels
        self.gc1.axis_labels.set_font(size='x-large')
        self.gc1.tick_labels.set_style('colons')
        self.gc1.tick_labels.set_xformat('hh:mm:ss.s')
        self.gc1.tick_labels.set_yformat('dd:mm:ss.s')
        self.gc1.tick_labels.set_font(size='small')
        self.gc1.ticks.set_color('black')
        self.gc1.ticks.show()
        self.gc1.axis_labels.set_ypad(-20)

        # scalebar
        self.gc1.add_scalebar(self.scalebar/3600.0,color=self.textcolor)
        self.gc1.scalebar.set_corner('bottom left')
        self.gc1.scalebar.set_label(str(self.scalebar)+'" ('+str(self.scalebar*self.distance)+' AU)')
        self.gc1.scalebar.set_linewidth(6)
        self.gc1.scalebar.set_font_size('xx-large')

        # beam
        self.gc1.beam.set_corner('bottom right')
        self.gc1.beam.set_color(self.textcolor)
        self.gc1.beam.set_hatch('+')

        # title/misc labels
        self.gc1.add_label(0.5, 0.95, r'{}'.format(self.name), \
            relative=True,size='x-large',color=self.textcolor,weight='heavy')
        self.gc1.add_label(0.1, 0.95, r'{}'.format(self.plotlabel), \
            relative=True,size='x-large',color=self.textcolor,weight='heavy')
        self.gc1.add_label(0.5, 0.875, r'{}'.format(self.misc), relative=True,\
            size='large',color=self.textcolor,weight='heavy')


    def drawOutflow(self,ra,dec,paR,paB):
        """draws outflows
        @input : xlength draws an outflow of this xlength
        @input : ylength "" ylength
        @input : ra places outflow at this pos
        @input : dec ""
        @input : paR and paB position angles of Red and blue outflow
        """
        xlength = self.xsize/8.
        ylength = self.ysize/8.
        if paR != None:
            dxred=xlength*np.cos(paR)
            dyred=ylength*np.sin(paR)
            self.gc1.show_arrows(ra,dec, dxred, dyred, color='red')
        if paB != None:
            dxblue=xlength*np.cos(paB)
            dyblue=ylength*np.sin(paB)
            self.gc1.show_arrows(ra,dec, dxblue, dyblue, color='blue')

    def showSources(self,ra,dec,color='white'):
        """
        """
        self.gc1.show_markers(ra,dec, c=color,marker='+',zorder=20)

    def showColorbar(self):
        """
        """
        self.gc1.add_colorbar()
        self.gc1.colorbar.set_width(0.1)
        self.gc1.colorbar.set_location('right')

    def addLabel(self,color,label):
        """
        """
        try:
            self.labelNum += 0.02
            labelNum = self.labelNum
        except:
            labelNum = 0.02
            self.labelNum = 0.02
        self.gc1.add_label(0.13, 0.98-labelNum, r'{} {}'.format(color,label), \
            relative=True,size='large',color=color)

    def save(self,outfilename,dpi=400): 
        """plotting module
        @input extension tuple of extension types
        @input dpi, 400 is pretty high quality
        """ 
        self.gc1.list_layers()
        #mpl.tight_layout()
        for ext in self.extension:
            self.fig.savefig(outfilename+'.'+ext,dpi=dpi,adjust_bbox=True,format=ext)

# vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
if __name__ == "__main__":
    # general file structure. Use -M- for molecule, 
    # -D- for doppler (r==red b==blue), and -G- for general qualifier(medium,low,etc)

    fileStructure = "../../L1448IRS3B_-M-_image_taper1500k_-D-_-G-.fits"
    distance = 293 # distance to source in pc
    colorBars = False # plot color bars?

    # if supplied will overlay continuum
    # also can set color True,False, or Both
    # also set the min/max flux values for pixels
    continuumOverlay = ("../../L1448IRS3B_cont_robust0.5.fits","false") 
    continuumParams = (0.001,0.04)
    fluxStretch = 'sqrt' # ways to stretch the flux for faint things
    cmap = 'ds9bb' # color map to use (only works if color mode is specified)
    txtColor = 'black' # overrides, can be left blank for auto generation

    # molecules the length of this determines most of the other params.
    # All the plot params need to have a length on axis 0 of the same length as mols
    mols = ("C17O","H13COp")

    # plotting
    imSize = ((7,7),(7,7)) # size of the output image in arcseconds
    title = (r"L1448N C$^{17}$O",r"L1448N H$^{13}$CO+")
    miscLabel1 = ('','')
    miscLabel2 = ('','')
    center=("03:25:36.32 30.45.14.8","03:25:36.32 30.45.14.8")
    scalebar = (1.,1.) # scalebar to add in arcsec
    show = ((True,True),(True,True)) # show
    #show =((False,False),(False,False))
    allowNegContours = (False,False)
    # Plot source markers. [True/False,((color,"RA DEC",label))] doesnt need length of mols
    markers = [True,(("red","03:25:36.3290 30.45.15.042",'+ Kinematic'),\
                     ("yellow","03:25:36.382 30.45.14.715",'+ tertiary'),\
                     ("black","03:25:36.32 30.45.14.8",'+ max PV'))]

    # draw outflows? and position angle of red blue respectively
    # this must match the # of markers in markers[2]
    outflows = ([False,0,0],[False,0,0],[False,0,0])

    # Specific vel ranges or "qualifiers"
    # for every set of quals it will plot them concurrently
    # only supports up to total num colors/2
    # labels correspond to the labels that will be plotted on the graph
    quals = (('ultralow','low','med','high','ultrahigh-2x'),
             ('lowbin','highbin'))
    quals = (('low','ultrahigh-2x'),
             ('lowbin','highbin'))
    # labelsfollow the major qualifier then red blue then the minor qualifier
    showLabels = True # show the labels?
    labels = ((('ultralow','low','med','high','ultrahigh-2x'),('ultralow','low','med','high','ultrahigh-2x')),
             (('lowbin','highbin'),('lowbin','highbin')))
    labels = ((('low','high','ultrahigh-2x'),('low','high','ultrahigh-2x')),
             (('lowbin','highbin'),('lowbin','highbin')))
    output = ("splitMoments","binnedMoments") # output formatter. fileStructure.replace('.fits',output+ext)
    extensions = ('pdf',) # extensions for formatted output
    colors = (("firebrick","red","magenta","tomato","indianred"),# reds 
              ("royalblue","C","blue","slateblue","mediumslateblue")) # blues
    # now there might be a lot of contours, you will need 1 per qual per RB side
    # so in total 2 * #totalquals
    # defined as sets with start,interval,rms
    contours = ((# molecule 1
                 (# qual Major 1
                  # these iterations are now the minor quals with each row being Red and Ble
                  #((7,3,0.0013),
                  # (7,3,0.0012)),
                  ((7,3,0.0012),
                   (5,3,0.0017)),
                  #((7,3,0.0013),
                  # (7,3,0.0012)), 
                  #((7,3,0.0013),
                  # (7,3,0.0015)),
                  ((7,3,0.0011),
                   (7,3,0.0012))),
                 (# qual Major 2
                  ((7,3,0.0029),
                   (7,3,0.0029)),
                  ((7,3,0.0024),
                   (7,3,0.0023))
                 )
                ),
                (# molecule 2
                 (# qual Major 1
                  # these iterations are now the minor quals with each row being Red and Ble
                  #((7,3,0.0013),
                  # (7,3,0.0012)),
                  ((4,3,0.0012),
                   (4,3,0.0017)),
                  #((7,3,0.0013),
                  # (7,3,0.0012)), 
                  #((7,3,0.0013),
                  # (7,3,0.0015)),
                  ((4,3,0.0011),
                   (4,3,0.0012))),
                 (# qual Major 2
                  ((4,3,0.0038),
                   (4,3,0.0019)),
                  ((4,3,0.0039),
                   (4,3,0.0039))
                 )
               ))

# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for molI,molV in enumerate(mols):
        for qualI1,qualV1 in enumerate(quals):
            ra,dec = map(float,parseRaDec(center[molI]).split(' '))
            if continuumOverlay[1] == 'both':
                for colorOverlay in ('true','false'):
                    # setup plot
                    plot = mappingContours(continuumOverlay[0],ra,dec,minpixval=continuumParams[0],\
                        maxpixval=continuumParams[1], size=imSize[molI],scalebar=scalebar[molI],\
                        distance=distance,name=title[molI],imagestretch=fluxStretch,\
                        colororgray=colorOverlay,colormap=cmap,plotlabel=miscLabel1[molI],\
                        misc=miscLabel1[molI],textcolor=txtColor,extension=extensions,showLabels=showLabels)
                    if colorOverlay == 'true':
                        ite = 'color'
                        if plot.auto:
                            plot.textcolor = 'white'
                    else:
                        ite = 'bw'
                        if plot.auto:
                            plot.textcolor = 'black'
                    plot.drawContinuum()
                    plot.plotFormat()
                    # plotting rb channels for all qualifiers specified
                    start = True
                    for qualI2,qualV2 in enumerate(qualV1):
                        for dopplerI,dopplerV in enumerate(('r','b')):
                            inputFname = fileStructure.replace('-M-',molV).replace('-D-',dopplerV)\
                                                      .replace('-G-',qualV2)
                            print(colors[dopplerI][qualI2],start,\
                                show[molI][dopplerI],inputFname,allowNegContours[molI],*contours[molI][qualI1][qualI2][dopplerI],labels[qualI1][dopplerI][qualI2])
                            plot.setupContours(colors[dopplerI][qualI2],qualI2 == 0,\
                                show[molI][dopplerI],inputFname,allowNegContours[molI],*contours[molI][qualI1][qualI2][dopplerI],labels[qualI1][dopplerI][qualI2])
                            start = False
                    # plotting markers and outflows
                    # plotting markers and outflows
                    if markers[0]:
                        for sourceI,sourceV in enumerate(markers[1]):
                            # plot marker
                            mRa,mDec = map(float,parseRaDec(sourceV[1]).split(' '))
                            plot.showSources(mRa,mDec,sourceV[0])
                            plot.addLabel(sourceV[0],sourceV[2])
                            # plot outflow
                            if outflows[sourceI][0]:
                                plot.drawOutflow(mRa,mDec,outflows[sourceI][1],outflows[sourceI][2])

                    # now plotting final steps and save
                    plot.plotFormat()
                    if colorBars:
                        plot.showColorbar()

                    outputFname = fileStructure.split('/')[-1].replace('-M-',molV).replace('-D-','')\
                                                      .replace('-G-','')\
                                                      .replace('.fits',output[qualI1]+'_'+ite)
                    plot.save(outputFname,200)
                    plot = None

            else:
                # setup plot
                colorOverlay = continuumOverlay[1]
                plot = mappingContours(continuumOverlay[0],ra,dec,minpixval=continuumParams[0],\
                    maxpixval=continuumParams[1], size=imSize[molI],scalebar=scalebar[molI],\
                    distance=distance,name=title[molI],imagestretch=fluxStretch,\
                    colororgray=colorOverlay,colormap=cmap,plotlabel=miscLabel1[molI],\
                    misc=miscLabel1[molI],textcolor=txtColor,extension=extensions,showLabels=showLabels)
                if plot.auto:
                    plot.textcolor = 'black'
                plot.drawContinuum()
                # plotting rb channels for all qualifiers specified
                start = True
                for qualI2,qualV2 in enumerate(qualV1):
                    for dopplerI,dopplerV in enumerate(('r','b')):
                        inputFname = fileStructure.replace('-M-',molV).replace('-D-',dopplerV)\
                                                  .replace('-G-',qualV2)
                        print((colors[dopplerI][qualI2],start,\
                                show[molI][dopplerI],inputFname,allowNegContours[molI],*contours[molI][qualI1][qualI2][dopplerI],labels[qualI1][dopplerI][qualI2]))
                        plot.setupContours(colors[dopplerI][qualI2], start,\
                            show[molI][dopplerI],inputFname,allowNegContours[molI],*contours[molI][qualI1][qualI2][dopplerI],labels[qualI1][dopplerI][qualI2])
                        start = False
                # plotting markers and outflows
                if markers[0]:
                    for sourceI,sourceV in enumerate(markers[1]):
                        # plot marker
                        mRa,mDec = map(float,parseRaDec(sourceV[1]).split(' '))
                        plot.showSources(mRa,mDec,sourceV[0])
                        plot.addLabel(sourceV[0],sourceV[2])
                        # plot outflow
                        if outflows[sourceI][0]:
                            plot.drawOutflow(mRa,mDec,outflows[sourceI][1],outflows[sourceI][2])

                # now plotting final steps and save
                plot.plotFormat()
                if colorBars:
                    plot.showColorbar()

                outputFname = "./"+fileStructure.split('/')[-1].replace('-M-',molV).replace('-D-','')\
                                                  .replace('-G-','')\
                                                  .replace('.fits',output[qualI1])
                plot.save(outputFname,200)
                plot = None
