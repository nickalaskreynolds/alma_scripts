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

linevis = source+'.H13COp.ms.selfcal.contsub.concat'
linename = 'H13COp' # name of transition (see science goals in OT for name) 
restfreq='346.998347GHz' # Typically the rest frequency of the line of      
spw='' # uncomment and replace with appropriate spw if necessary.

#plotms(vis=linevis,averagedata=True,avgtime='1e8',avgscan=True,xaxis='velocity',restfreq=restfreq)

start='-5km/s' # start velocity. See science goals for appropriate value.
width='0.2km/s'# velocity width. See science goals.
nchan = 100  # number of channels. See science goals for appropriate value.
mode='velocity'
uvrange='>100klambda'


##########################
# 500 iterations
lineimagename = source+'_'+linename+'_dirty_image_2_binned'
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
      imagermode=imagermode,
      uvrange=uvrange)


##########################
# 500 iterations
lineimagename = source+'_'+linename+'_clean_image_2_binned'
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
      imagermode=imagermode,
      uvrange=uvrange)


##########################
uvtaper=True
outertaper='1500klambda'
# 400 iterations
lineimagename = source+'_'+linename+'_image_2_binned_taper1500k'
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
      imagermode=imagermode,
      uvtaper=uvtaper,
      outertaper=outertaper,
      uvrange=uvrange)

##########################
uvtaper=True
outertaper='1000klambda'

# 500 iterations
lineimagename = source+'_'+linename+'_image_2_binned_taper1000k'
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
      imagermode=imagermode,
      uvtaper=uvtaper,
      outertaper=outertaper,
      uvrange=uvrange)


#!rm -rf *.model *.mask *.psf *.residual *p0* *p1* *p3*


