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

# for each line and continuum
linename = 'H13COp' # name of transition (see science goals in OT for name) 
restfreq='346.998347GHz' # Typically the rest frequency of the line of      
spw='' # uncomment and replace with appropriate spw if necessary.
start='-0km/s' # start velocity. See science goals for appropriate value.
width='0.1km/s'# velocity width. See science goals.
nchan = 100  # number of channels. See science goals for appropriate value.
mode='velocity'

sourcems = source + "." + linename + ".ms"
concatms = source + "." + linename + "_comb.ms"

#plotms(vis=linevis,averagedata=True,avgtime='1e8',avgscan=True,xaxis='velocity',restfreq=restfreq)
concat(vis=['../../ALMA-Per33-hr/'+sourcems,'../../ALMA-Per33-lr/'+sourcems],concatvis=concatms,visweightscale=[1.0,0.3])
