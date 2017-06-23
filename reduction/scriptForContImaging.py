source='L1448IRS3B'
contvis = source+'_cont.ms.apselfcal'         

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

# 7000 iterations
contimagename=source+'_cont_robust2'
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


robust=0.5
contimagename=source+'_cont_robust0.5'
os.system('rm -rf '+contimagename+'*')
# 5000 iterations
clean(vis=contvis,
      imagename=contimagename,
      field=field,
#      #phasecenter=phasecenter, # uncomment if mosaic.      
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

robust=0.5
contimagename=source+'_cont_robust0.5_gt25klambda'
rangelimit='>25klambda'
os.system('rm -rf '+contimagename+'*')
# 810 iterations
clean(vis=contvis,
      imagename=contimagename,
      field=field,
#      #phasecenter=phasecenter, # uncomment if mosaic.      
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      spw='',
      interactive = True,
      imagermode = imagermode,uvrange=rangelimit)

####

robust=-0.5
rangelimit=''
contimagename=source+'_cont_robust-0.5'
os.system('rm -rf '+contimagename+'*')
# 290 iterations
clean(vis=contvis,
      imagename=contimagename,
      field=field,
#      #phasecenter=phasecenter, # uncomment if mosaic.      
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      spw='',
      interactive = True,
      imagermode = imagermode)


robust=-0.5
contimagename=source+'_cont_robust-0.5_gt25klambda'
rangelimit='>25klambda'
os.system('rm -rf '+contimagename+'*')
# 530 iterations
clean(vis=contvis,
      imagename=contimagename,
      field=field,
#      #phasecenter=phasecenter, # uncomment if mosaic.      
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      spw='',
      interactive = True,
      imagermode = imagermode,uvrange=rangelimit)


robust=-1
contimagename=source+'_cont_robust-1'
os.system('rm -rf '+contimagename+'*')
# 5200 Iterations
clean(vis=contvis,
      imagename=contimagename,
      field=field,
#      #phasecenter=phasecenter, # uncomment if mosaic.      
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      spw='',
      interactive = True,
      imagermode = imagermode,uvrange='')


robust=-2
contimagename=source+'_cont_robust-2'
os.system('rm -rf '+contimagename+'*')
# 14500 iterations
clean(vis=contvis,
      imagename=contimagename,
      field=field,
#      #phasecenter=phasecenter, # uncomment if mosaic.      
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      spw='',
      interactive = True,
      imagermode = imagermode)


weighting='superuniform'
robust=-2
contimagename=source+'_cont_superuniform'
os.system('rm -rf '+contimagename+'*')
# 20000 iterations
clean(vis=contvis,
      imagename=contimagename,
      field=field,
#      #phasecenter=phasecenter, # uncomment if mosaic.      
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      spw='',
      interactive = True,
      imagermode = imagermode)


weighting='briggs'
robust=2
contimagename=source+'_cont_robust2_taper500kl'
os.system('rm -rf '+contimagename+'*')
# 10000 iterations
clean(vis=contvis,
      imagename=contimagename,
      field=field,
#      #phasecenter=phasecenter, # uncomment if mosaic.      
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      niter = niter, 
      threshold = threshold, 
      spw='',
      interactive = True,
      imagermode = imagermode,uvtaper=True,outertaper='500klambda',uvrange='')

weighting='superuniform'
robust=-2
contimagename=source+'_cont_superuniform_robust-2_taper300kl'
os.system('rm -rf '+contimagename+'*')
# 20000 iterations
clean(vis=contvis,
      imagename=contimagename,
      field=field,
#      #phasecenter=phasecenter, # uncomment if mosaic.      
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust = robust,
      niter = niter, 
      threshold = threshold, 
      spw='',
      interactive = True,
      imagermode = imagermode,outertaper='500klambda')


