source=['L1448IRS3B_cont_robust2.deproj']#['continuum-concat','continuum-tobin']
offset=10.31E-06 #[7.034589E-03,6.188E-03]
import glob
import os


# summary of above commands

dooff=True
fixoffset=True
stretch=True

_TEMP_='L1448IRS3B_cont_robust2.deproj'
imagename=_TEMP_+'.image'
region=_TEMP_+'.rgn'
estimates='continuum-concat-deproj_0_point_output.estimates'
newvis='L1448IRS3B.cont.ms.apselfcal.concat.deproj'
frequ='335.500GHz'

'''
# peak intensity must be in map units
# 'f' (peak intensity), 'x' (peak x position), 'y' (peak y position), 'a' (major axis), 'b' (minor axis), 'p' (position angle)  
0.09, 1022, 1026, 0.67arcsec, 0.54arcsec, c, 120deg, abp 
'''
# sing source

imfit(imagename=imagename,region=region,model='continuum-concat-deproj_0_gauss_output.model',residual='continuum-concat-deproj_0_gauss_output.residual',logfile='continuum-concat-deproj_0_gauss_output.logfile',summary='continuum-concat-deproj_0_gauss_output.summary',newestimates='continuum-concat-deproj_0_gauss_output.newestimates',dooff=dooff,offset=offset)

imfit(imagename=imagename,region=region,model='continuum-concat-deproj_0_point_output.model',residual='continuum-concat-deproj_0_point_output.residual',logfile='continuum-concat-deproj_0_point_output.logfile',summary='continuum-concat-deproj_0_point_output.summary',estimates=estimates,newestimates='continuum-concat-deproj_0_point_output.newestimates',dooff=dooff,fixoffset=fixoffset,offset=offset)


# both point and gauss

estimates='continuum-concat-deproj_0_point_output.estimates'
imfit(imagename=imagename,region=region,model='continuum-concat-deproj_0_point_gauss_output.model',residual='continuum-concat-deproj_0_point_gauss_output.residual',logfile='continuum-concat-deproj_0_point_gauss_output.logfile',summary='continuum-concat-deproj_0_point_gauss_output.summary',estimates=estimates,newestimates='continuum-concat-deproj_0_point_gauss_output.newestimates',dooff=dooff,fixoffset=fixoffset,offset=offset)





immath(imagename=['L1448IRS3B_cont_robust2.deproj.image','continuum-concat-deproj_0_gauss_output.model'],expr='IM0-IM1',  outfile='L1448IRS3B_cont_robust2.deproj.subclump_gauss.image') 

immath(imagename=['L1448IRS3B_cont_robust2.deproj.image','continuum-concat-deproj_0_point_output.model.bak'],expr='IM0-IM1',  outfile='L1448IRS3B_cont_robust2.deproj.subclump_point.image') 