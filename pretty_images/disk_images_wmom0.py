import aplpy
import numpy as np
import pylab as P
import matplotlib.pyplot as mpl
import math
import sys
from readcol import *
#from tycholib import *
import os
#import ds9cmap
#import viridis
from astropy import wcs
from astropy.io import fits

def main(image,outfilename,ra,dec,minpixval,maxpixval,\
   size,scalebar,distance,name,pa,showoutflow,sigma,showcontours,showsources,\
   imagestretch,redimage,redstart,redinterval,rednoise,\
   showredimage,blueimage,bluestart,blueinterval,bluenoise,showblueimage,\
   colororgray,colormap,plotlabel,textcolor,extension):
   def standardStuff():
      gc1.axis_labels.set_font(size='x-large')
      gc1.tick_labels.set_style('colons')
      gc1.tick_labels.set_xformat('hh:mm:ss.s')
      gc1.tick_labels.set_yformat('dd:mm:ss.s')
      gc1.tick_labels.set_font(size='small')
      gc1.ticks.set_color('black')
      gc1.ticks.show()
      gc1.axis_labels.set_ypad(-20)
      #gc1.tick_labels.hide()
      #gc1.ticks.hide()

      gc1.add_scalebar(scalebar/3600.0,color=textcolor)
      gc1.scalebar.set_corner('bottom left')
      gc1.scalebar.set_label(str(scalebar)+'" ('+str(scalebar*distance)+' AU)')
      gc1.scalebar.set_linewidth(6)
      gc1.scalebar.set_font_size('xx-large')

      gc1.add_beam()
      gc1.beam.set_corner('bottom right')
      gc1.beam.set_color(textcolor)
      gc1.beam.set_hatch('+')

   def drawOutflow():
      dxred=0.0
      dyred=0.0
      dxblue=0.0
      dyblue=0.0
      length=size/8.0/3600.0
      if size > 5.0:
         length=size/8.0/3600.0
      else:
         length=size/4.0/3600.0

      if pa > 180.0:
        panew=pa-180.0 
      else:
        panew=pa
      angle=abs(abs(panew)-90.0)
      xlength=length*math.cos(angle*math.pi/180.0)
      ylength=length*math.sin(angle*math.pi/180.0)
     
      xpt1=xlength
      ypt1=ylength
      xpt2=-xlength
      ypt2=-ylength

      if pa == 180.0 or pa == 360.0:
         dxred=0.0
         dxblue=0.0
      if pa == 90.0 or pa == 270.0:
         dyred=0.0 
         dyblue=0.0
      if pa > 180.0 and pa <= 270.0:
         dxblue=xpt2
         dyblue=ypt2
         dxred=xpt1 
         dyred=ypt1
      if pa > 180.0 and pa > 270.0:
         dxblue=xpt2
         dyblue=ypt1
         dxred=xpt1
         dyred=ypt2
      if pa < 180.0 and pa > 90.0:
         dxblue=xpt1
         dyblue=ypt2
         dxred=xpt2
         dyred=ypt1
      if pa < 180.0 and pa <=90.0:
         dxblue=xpt1
         dyblue=ypt1
         dxred=xpt2
         dyred=ypt2
      #print dxblue*3600.0,dyblue*3600.0,dxred*3600.0,dyred*3600.0
   #   gc1.show_arrows(ra+dxblue/4.0,dec+dyblue/4.0, dxblue, dyblue, color='cyan')
   #   gc1.show_arrows(ra+dxred/4.0,dec+dyred/4.0, dxred, dyred, color='red')

      if size > 5.0:
         ranew=ra-size/6.0/3600.0
         decnew=dec-size/6.0/3600.0
      else:
         ranew=ra-size/3.0/3600.0
         decnew=dec-size/3.0/3600.0
      gc1.show_arrows(ranew,decnew, dxblue, dyblue, color='blue')
      gc1.show_arrows(ranew,decnew, dxred, dyred, color='red')

   def drawContours(contstart=None,continterval=None,contnoise=None):
      contours1=np.arange(contstart,contstart*10.0,continterval,dtype='float32')
      contours2=np.arange(contstart*10.0,contstart*40.0,continterval*3.0,dtype='float32')
      contours3=np.arange(contstart*40.0,contstart*100.0,continterval*10.0,dtype='float32')
      poscontours=np.concatenate((contours1,contours2,contours3))
      negcontours=poscontours[::-1]*(-1.0)
      contours = poscontours * contnoise
      #contours=np.concatenate((poscontours,negcontours))*contnoise
      
      return contours

   fig = mpl.figure(figsize=(16,12))

   dx=0.0
   dy=0.0
   #os.system('rm -rf '+ image+'cut.fits')
   #fitscut2d(image,image+'cut.fits',ra,dec,300)
   #image=image+'cut.fits'

   #print name, image, ra,dec
   gc1 = aplpy.FITSFigure(image, figure=fig, subplot=[0.15+dx,0.1+dy,0.7,0.9])

   if colororgray == 'color':
      gc1.show_colorscale(vmin=minpixval,vmax=maxpixval,stretch=imagestretch,cmap=colormap)
   if colororgray == 'gray':
      gc1.show_grayscale(vmin=minpixval,vmax=maxpixval,stretch=imagestretch,invert=True)

   gc1.recenter(ra,dec,width=(size/3600.0),height=(size/3600.0))

   if showsources == 'y':
      mainsource,sourcename,ra_src,dra_junk,dec_src,ddec_junk,freq9_junk,intflux9_junk,eintflux9_junk,pflux9_junk,rms9_junk,freq8_junk,intflux8_junk,eintflux8_junk,pflux8_junk,rms8_junk,freq10_junk,intflux10_junk,eintflux10_junk,pflux10_junk,rms10_junk,spindex_junk,espindex_junk,pspindex_junk,epspindex_junk=readcol('/data2/EVLA/Perseus/FluxTable/all-fitresults-multiples.txt',twod=False)
      if showblueimage == 'y':
         #gc1.show_markers(ra_src[:]-0.06/3600.0,dec_src[:]-0.035/3600.0, size/3600.0*0.03, c='white',marker='+')
         #gc1.show_markers(ra_src[:]-0.06/3600.0,dec_src[:]-0.035/3600.0,c='white')
         gc1.show_markers(ra_src[:]-0.06/3600.0,dec_src[:]-0.06/3600.0, c='white',marker='+',zorder=20)
      else:
         gc1.show_markers(ra_src[:],dec_src[:],c='red',marker='+',zorder=20)


   gc1.add_label(0.5, 0.95, name, relative=True,size='x-large',color=textcolor,weight='heavy')
   gc1.add_label(0.1, 0.95, plotlabel, relative=True,size='x-large',color=textcolor,weight='heavy')
   #gc1.add_label(0.5, 0.875, r'$\Delta$'+separation, relative=True,size='large',color='white',weight='heavy')

   #gc1.add_colorbar()
   #gc1.colorbar.set_width(0.1)
   #gc1.colorbar.set_location('right')

   standardStuff()

   if showredimage == 'y':
      contours=drawContours(redstart,redinterval,rednoise)
      gc1.show_contour(redimage,levels=contours,colors='red',linewidths=1.0)

   if showblueimage == 'y':
      contours=drawContours(bluestart,blueinterval,bluenoise)
      gc1.show_contour(blueimage,levels=contours,colors='blue',linewidths=1.0)

   if pa < 360.0 and showoutflow == 'y':
      drawOutflow()

   gc1.list_layers()
   #os.system('rm -rf '+ image)
   #os.system('rm -rf '+ redimage)
   #os.system('rm -rf '+ blueimage)

   fig.savefig(outfilename,dpi=400,adjust_bbox=True,format=extension)

