import aplpy
#print aplpy.version.version
import numpy as np
import pylab as P
import matplotlib.pyplot as mpl
import math
import sys
from readcol import *
#from tycholib import *
import os
import ds9cmap
#import viridis
from astropy import wcs
from astropy.io import fits

def main(image,outfilename,ra,dec,minpixval,maxpixval,size,scalebar,\
   distance,name,pa,showoutflow,sigma,showcontours,showsources,imagestretch,\
   redimage,redstart,redinterval,rednoise,showredimage,blueimage,bluestart,\
   blueinterval,bluenoise,showblueimage,colororgray,colormap,plotlabel,textcolor,\
   ra1,dec1,ra2,dec2,pa1,pa2,extension):
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
      dxred1=0.0
      dyred1=0.0
      dxblue1=0.0
      dyblue1=0.0
      length1=size/8.0/3600.0
      if size > 5.0:
         length1=size/8.0/3600.0
      else:
         length1=size/4.0/3600.0

      if pa1 > 180.0:
        panew1=pa1-180.0 
      else:
        panew1=pa1
      angle1 =abs(abs(panew1)-90.0)
      xlength1 =length1*math.cos(angle1*math.pi/180.0)
      ylength1 =length1*math.sin(angle1*math.pi/180.0)
     
      xpt11=xlength1
      ypt11=ylength1
      xpt21=-xlength1
      ypt21=-ylength1

      if pa1 == 180.0 or pa1 == 360.0:
         dxred1 =0.0
         dxblue1 =0.0
      if pa1 == 90.0 or pa1 == 270.0:
         dyred1 =0.0 
         dyblue1=0.0
      if pa1 > 180.0 and pa1 <= 270.0:
         dxblue1=xpt21
         dyblue1=ypt21
         dxred1=xpt11
         dyred1=ypt11
      if pa1 > 180.0 and pa1 > 270.0:
         dxblue1=xpt21
         dyblue1=ypt11
         dxred1=xpt11
         dyred1=ypt21
      if pa1 < 180.0 and pa1 > 90.0:
         dxblue1=xpt11
         dyblue1=ypt21
         dxred1=xpt21
         dyred1=ypt11
      if pa1 < 180.0 and pa1 <=90.0:
         dxblue1=xpt11
         dyblue1=ypt11
         dxred1=xpt21
         dyred1=ypt21
      #print dxblue1*3600.0,dyblue1*3600.0,dxred1*3600.0,dyred1*3600.0
   #   gc1.show_arrows(ra+dxblue/4.0,dec+dyblue/4.0, dxblue, dyblue, color='cyan')
   #   gc1.show_arrows(ra+dxred/4.0,dec+dyred/4.0, dxred, dyred, color='red')

      if size > 5.0:
         ranew1=ra1-size/6.0/3600.0
         decnew1=dec1-size/6.0/3600.0
      else:
         ranew1=ra1-size/3.0/3600.0
         decnew1=dec1-size/3.0/3600.0

      dxred2=0.0
      dyred2=0.0
      dxblue2=0.0
      dyblue2=0.0
      length2=size/8.0/3600.0
      if size > 5.0:
         length2=size/8.0/3600.0
      else:
         length2=size/4.0/3600.0

      if pa2 > 180.0:
        panew2=pa2-180.0 
      else:
        panew2=pa2
      angle2 =abs(abs(panew2)-90.0)
      xlength2 =length2*math.cos(angle2*math.pi/180.0)
      ylength2 =length2*math.sin(angle2*math.pi/180.0)
     
      xpt12=xlength2
      ypt12=ylength2
      xpt22=-xlength2
      ypt22=-ylength2

      if pa2 == 180.0 or pa2 == 360.0:
         dxred2 =0.0
         dxblue2 =0.0
      if pa2 == 90.0 or pa2 == 270.0:
         dyred2 =0.0 
         dyblue2=0.0
      if pa2 > 180.0 and pa2 <= 270.0:
         dxblue2=xpt22
         dyblue2=ypt22
         dxred2=xpt12
         dyred2=ypt12
      if pa2 > 180.0 and pa2 > 270.0:
         dxblue2=xpt22
         dyblue2=ypt12
         dxred2=xpt12
         dyred2=ypt22
      if pa2 < 180.0 and pa2 > 90.0:
         dxblue2=xpt12
         dyblue2=ypt22
         dxred2=xpt22
         dyred2=ypt12
      if pa2 < 180.0 and pa2 <=90.0:
         dxblue2=xpt12
         dyblue2=ypt12
         dxred2=xpt22
         dyred2=ypt22
      #print dxblue2*3600.0,dyblue2*3600.0,dxred2*3600.0,dyred2*3600.0

      if size > 5.0:
         ranew2=ra2-size/6.0/3600.0
         decnew2=dec2-size/6.0/3600.0
      else:
         ranew2=ra2-size/3.0/3600.0
         decnew2=dec2-size/3.0/3600.0
   
      gc1.show_arrows(ranew1,decnew1, dxblue1, dyblue1, color='blue')
      gc1.show_arrows(ranew1,decnew1, dxred1, dyred1, color='red')
      gc1.show_arrows(ranew2,decnew2, dxblue2, dyblue2, color='blue')
      gc1.show_arrows(ranew2,decnew2, dxred2, dyred2, color='magenta')   


   def drawContours(contstart=None,continterval=None,contnoise=None):
      contours1=np.arange(contstart,contstart*10.0,continterval,dtype='float32')
      contours2=np.arange(contstart*10.0,contstart*40.0,continterval*3.0,dtype='float32')
      contours3=np.arange(contstart*40.0,contstart*100.0,continterval*10.0,dtype='float32')
      poscontours=np.concatenate((contours1,contours2,contours3))
      negcontours=poscontours[::-1]*(-1.0)
      contours = poscontours * contnoise
      #contours=np.concatenate((poscontours,negcontours))*contnoise
      
      return contours

   fig = mpl.figure(figsize=(7,7))

   dx=0.0
   dy=0.0
   #os.system('rm -rf '+ image+'cut.fits')
   #fitscut2d(image,image+'cut.fits',ra,dec,300)
   #image=image+'cut.fits'

   #print name, image, ra,dec
   gc1 = aplpy.FITSFigure(image, figure=fig)

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
      contours.sort()
      gc1.show_contour(redimage,levels=contours,colors='red',linewidths=1.0)

   if showblueimage == 'y':
      contours=drawContours(bluestart,blueinterval,bluenoise)
      contours.sort()
      gc1.show_contour(blueimage,levels=contours,colors='blue',linewidths=1.0)

   if pa < 360.0 and showoutflow == 'y':
      drawOutflow()

   gc1.list_layers()
   #os.system('rm -rf '+ image)
   #os.system('rm -rf '+ redimage)
   #os.system('rm -rf '+ blueimage)

   fig.savefig(outfilename,dpi=400,adjust_bbox=True,format='pdf')

