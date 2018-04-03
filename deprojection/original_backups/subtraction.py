# values are read from estimate file from imfit
'''
0.0300588, 1027.17, 1023.47, 0.4763 arcsec, 0.4147 arcsec, 29.5 deg
0.05599, 1020.42, 1025.35, 0.301785 arcsec, 0.203298 arcsec, 112.683 deg

'''
frequ='335.500GHz'
newmodel="L1448IRS3B_cont_robust2.deproj.presubclump.model"
newvis='L1448IRS3B.cont.ms.apselfcal.concat.deproj.subclump'

cl.purge()
cl.done()
#cl.open('L1527.A.FF.cl')
direction='J2000 03h25m36.383s +030d45m14.729s'
directionG='J2000 03h25m36.379s +030d45m14.704s'
directionG1='J2000 03h25m36.38s +030d45m14.727s'
directionG2='J2000 03h25m36.38s +030d45m14.700s'
directionP='J2000 03h25m36.384s +030d45m14.723s'
directionP2='J2000 03h25m36.389s +030d45m14.770s'
directionG3='J2000 03h25m36.390s +030d45m14.765s'

cl.addcomponent(dir=directionP, flux=0.010505, fluxunit='Jy', freq=frequ, shape="Point",positionangle='29.5deg')#, majoraxis="0.5763arcsec", minoraxis='0.4147arcsec', )7
cl.addcomponent(dir=directionP2, flux=0.04, fluxunit='Jy', freq=frequ, shape="Gaussian", majoraxis="0.8763arcsec", minoraxis='0.8447arcsec', positionangle='112.654deg')
cl.addcomponent(dir=directionG, flux=0.04056, fluxunit='Jy', freq=frequ, shape="Gaussian", majoraxis="1.01785arcsec", minoraxis='0.7503298arcsec', positionangle='112.654deg')
cl.addcomponent(dir=directionG1, flux=0.045, fluxunit='Jy', freq=frequ, shape="Gaussian", majoraxis="1.51785arcsec", minoraxis='1.23298arcsec', positionangle='112.654deg')
cl.addcomponent(dir=directionG2, flux=0.06, fluxunit='Jy', freq=frequ, shape="Gaussian", majoraxis="0.25785arcsec", minoraxis='0.23298arcsec', positionangle='112.654deg')
cl.addcomponent(dir=directionG3, flux=0.04, fluxunit='Jy', freq=frequ, shape="Gaussian", majoraxis="0.25785arcsec", minoraxis='0.23298arcsec', positionangle='112.654deg')



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
os.system('cp -rf ./' + newvis + '.bak ./'+newvis)
setjy(vis=newvis,model=newmodel,standard='manual',fluxdensity=0.0051,reffreq=frequ,spix=0.0)
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
      niter = 0, 
      threshold = '1000Jy', 
      interactive =False,
      imagermode = imagermode)


# 7200 iterations
contimagename=source+'_cont_robust2_image.deproj.subclump'
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
