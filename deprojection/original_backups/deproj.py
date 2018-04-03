vis='L1448IRS3B.cont.ms.apselfcal.concat'
fitsfile='L1448IRS3B.cont.ms.apselfcal.concat.uv'
field='0'
spw='5,11'

exportuvfits(vis=vis,spw=spw,field=field,fitsfile=fitsfile,datacolumn='corrected')

###########

# run modify_uv.py now
# python3 modify_uv.py --input L1448IRS3B.cont.ms.apselfcal.concat.uv -p -29.5 -i 45 -o L1448IRS3B.cont.ms.apselfcal.concat.

###########

fitsfile='L1448IRS3B.cont.ms.apselfcal.conc.deproj.fits'
vis='L1448IRS3B.cont.ms.apselfcal.concat.deproj'

importuvfits(fitsfile=fitsfile,vis=vis)

###########

source="L1448IRS3B"
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


#############################################
# Imaging the Continuuum
# Set the ms and continuum image name.

# If necessary, run the following commands to get rid of older clean
# data.
#clearcal(vis=contvis)
#delmod(vis=contvis)

contimagename=source+'_cont_dirty_image.deproj'
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
contimagename=source+'_cont_robust2.deproj'
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

###################################################

vis='L1448IRS3B.C17O.ms.selfcal.contsub.concat'
fitsfile='L1448IRS3B.C17O.ms.selfcal.contsub.concat.uv'
field='0'
spw='0,1,2'
os.system('rm -rf ' + fitsfile)

exportuvfits(vis=vis,spw=spw,field=field,fitsfile=fitsfile,datacolumn='corrected')

###########

# run modify_uv.py now
# python3 modify_uv.py --input L1448IRS3B.C17O.ms.selfcal.contsub.concat.uv -p -29.5 -i 45 -o L1448IRS3B.C17O.ms.selfcal.contsub.concat

###########

fitsfile='L1448IRS3B.cont.ms.apselfcal.concat.deproj.fits'
vis='L1448IRS3B.cont.ms.apselfcal.concat.deproj'
os.system('rm -rf ' + vis)
importuvfits(fitsfile=fitsfile,vis=vis)

source='L1448IRS3B'
field='0' 
imagermode='csclean' 
cell='0.02arcsec' # cell size for imaging.
imsize = [2048,2048] # size of image in pixels.
outframe='lsrk' # velocity reference frame. See science goals.
veltype='radio' # velocity type. See note below.
weighting = 'natural'
robust=0.5
niter=10000
threshold = '0mJy'

##############################################
# Image line emission [REPEAT AS NECESSARY]

linevis = source+'.C17O.ms.selfcal.contsub.concat'
linename = 'C17O' # name of transition (see science goals in OT for name) 
restfreq='337.06121GHz' # Typically the rest frequency of the line of      
spw='' # uncomment and replace with appropriate spw if necessary.

#plotms(vis=linevis,averagedata=True,avgtime='1e8',avgscan=True,xaxis='velocity',restfreq=restfreq)

start='-5km/s' # start velocity. See science goals for appropriate value.
width='0.11km/s' # velocity width. See science goals.
nchan = 200  # number of channels. See science goals for appropriate value.
mode='velocity'


##########################
lineimagename = source+'_'+linename+'_dirty_image'
os.system('rm -rf '+lineimagename+'.*')
clean(vis=linevis,
      imagename=lineimagename, 
      field=field,
      spw=spw,
#      phasecenter=phasecenter, # uncomment if mosaic.      
      mode=mode,
      start=start,
      width=width,
      nchan=nchan, 
      outframe=outframe, 
      veltype=veltype, 
      restfreq=restfreq, 
      niter=niter,  
      threshold=threshold, 
      interactive=True,
      cell=cell,
      imsize=imsize, 
      weighting=weighting, 
      imagermode=imagermode)



