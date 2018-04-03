import glob
import os

_TEMP_='L1448IRS3B_cont_robust0.5'

# both point and gauss
default(imfit)
imagename='L1448IRS3B_cont_robust0.5.image'
region='L1448IRS3B_cont_robust0.5.rgn'
logfile='L1448IRS3B_cont_robust0.5.sub.imfitlog'
estimates='L1448IRS3B_cont_robust0.5.sub.imfitest'
summary='L1448IRS3B_cont_robust0.5.sub.imfitsummary'
model='L1448IRS3B_cont_robust0.5.sub.imfitmodel'
residual='L1448IRS3B_cont_robust0.5.sub.imfitresidual'
newestimates='L1448IRS3B_cont_robust0.5.sub.imfitnewest'
dooff=True
excludepix=[-1e10,0]
rms=2E-02
overwrite=True
go()

viewer('L1448IRS3B_cont_robust0.5.sub.imfitresidual')

_TEMP_='L1448IRS3B_cont_robust0.5'

outfile='L1448IRS3B_cont_robust0.5.subclump.image'
os.system('rm -rf ' + outfile)
immath(imagename=['L1448IRS3B_cont_robust0.5.image','L1448IRS3B_cont_robust0.5.sub.imfitmodel'],expr='IM0-IM1', outfile=outfile)

exportfits(outfile,outfile[:-6]+'.fits',overwrite=True)

