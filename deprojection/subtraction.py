_TEMP_='L1448IRS3B_cont_robust2.deproj'
imagename=_TEMP_+'.image'
region=_TEMP_+'.rgn'
frequ='335.500GHz'
newmodel="L1448IRS3B_cont_robust2.deproj.model"
newvis='L1448IRS3B.cont.ms.apselfcal.concat.deproj'

# values are read from estimate file from imfit

cl.done()
#cl.open('L1527.A.FF.cl')
direction='J2000 03h25m36.37176s +030d45m14.74540s'
direction1='J2000 03h25m36.382s +030d45m14.720s'
direction2='J2000 03h25m36.380s +030d45m14.716s'
cl.addcomponent(dir=direction1, flux=0.0651926, fluxunit='Jy', freq=frequ, shape="Gaussian", 
                majoraxis="0.327832arcsec", minoraxis='0.230896arcsec', positionangle='113.673deg')
cl.addcomponent(dir=direction2, flux=0.0249433, fluxunit='Jy', freq=frequ, shape="Point", 
                majoraxis="0.67arcsec", minoraxis='0.54arcsec', positionangle='29.5deg')

ia.fromshape(newmodel,[2048,2048,1,1],overwrite=True)
cs=ia.coordsys()
cs.setunits(['rad','rad','','Hz'])
cell_rad=qa.convert(qa.quantity("0.01arcsec"),"rad")['value']
cs.setincrement([-cell_rad,cell_rad],'direction')
cs.setreferencevalue([qa.convert("3.426755556h",'rad')['value'],qa.convert("30.754147222deg",'rad')['value']],type="direction")

cs.setreferencevalue(frequ,'spectral')
cs.setincrement('3.0GHz','spectral')
ia.setcoordsys(cs.torecord())
ia.setbrightnessunit("Jy/pixel")
ia.modify(cl.torecord(),subtract=False)

os.system('rm -rf ' + newvis) 
os.system('cp -rf ./deproj/' + newvis + ' ./')
setjy(vis=newvis,model=newmodel,standard='manual',fluxdensity=0.15,reffreq=frequ,spix=0.0)
uvsub(vis=newvis)

contvis = vis
field='0' 
cell='0.01arcsec' # cell size for imaging.
imsize = [2048,2048] # size of image in pixels.
outframe='lsrk' # velocity reference frame. See science goals.
veltype='radio' # velocity type. See note below.
weighting = 'briggs'
niter=1000
threshold = '0.0mJy'
imagermode='csclean'

contvis=newvis

source='L1448IRS3B'
contimagename=source+'_cont_dirty_image.deproj.subclump'
os.system('rm -rf '+contimagename+'*')
clean(vis=contvis,
      imagename=contimagename,
      field=field,
      #phasecenter=phasecenter, # uncomment if mosaic.      
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      niter = niter, 
      threshold = threshold, 
      interactive = True,
      imagermode = imagermode)


# 7200 iterations
contimagename=source+'_cont_robust2.deproj.subclump'
robust=2.0
os.system('rm -rf '+contimagename+'*')
clean(vis=contvis,
      imagename=contimagename,
      field=field,
      #phasecenter=phasecenter, # uncomment if mosaic.      
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      interactive = True,
      imagermode = imagermode)
