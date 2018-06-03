# Run commands in order, otherwise might have conflicts unresolved
source='Per33_L1448_IRS3B'
sourcename='L1448IRS3B'
sourcems='L1448IRS3B.ms'

plotms(vis='uid___A002_Xbbe66a_X1954.ms',field='4',xaxis='time',yaxis='amplitude',coloraxis='field',spw='20~40',averagedata=True,avgscan=True)

#Preserving the original file


#Listobs output
  SpwID  Name                                       #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs  
  0      BB_1#SQLD                                       1   TOPO  282994.575   2000000.000   2000000.0 282994.5750        1  XX  YY
  1      BB_2#SQLD                                       1   TOPO  284932.075   2000000.000   2000000.0 284932.0750        2  XX  YY
  2      BB_3#SQLD                                       1   TOPO  294994.575   2000000.000   2000000.0 294994.5750        3  XX  YY
  3      BB_4#SQLD                                       1   TOPO  296994.575   2000000.000   2000000.0 296994.5750        4  XX  YY
  4      WVR#NOMINAL                                     4   TOPO  184550.000   1500000.000   7500000.0 187550.0000        0  XX
  5      X0000000000#ALMA_RB_07#BB_1#SW-01#FULL_RES    128   TOPO  283986.763    -15625.000   2000000.0 282994.5750        1  XX  YY
  6      X0000000000#ALMA_RB_07#BB_1#SW-01#CH_AVG        1   TOPO  282971.138   1781250.000   1781250.0 282971.1375        1  XX  YY
  7      X0000000000#ALMA_RB_07#BB_2#SW-01#FULL_RES    128   TOPO  285924.263    -15625.000   2000000.0 284932.0750        2  XX  YY
  8      X0000000000#ALMA_RB_07#BB_2#SW-01#CH_AVG        1   TOPO  284908.638   1781250.000   1781250.0 284908.6375        2  XX  YY
  9      X0000000000#ALMA_RB_07#BB_3#SW-01#FULL_RES    128   TOPO  294002.388     15625.000   2000000.0 294994.5750        3  XX  YY
  10     X0000000000#ALMA_RB_07#BB_3#SW-01#CH_AVG        1   TOPO  294971.138   1781250.000   1781250.0 294971.1375        3  XX  YY
  11     X0000000000#ALMA_RB_07#BB_4#SW-01#FULL_RES    128   TOPO  296002.388     15625.000   2000000.0 296994.5750        4  XX  YY
  12     X0000000000#ALMA_RB_07#BB_4#SW-01#CH_AVG        1   TOPO  296971.138   1781250.000   1781250.0 296971.1375        4  XX  YY
  13     BB_1#SQLD                                       1   TOPO  346525.932   2000000.000   2000000.0 346525.9316        1  XX  YY
  14     BB_2#SQLD                                       1   TOPO  346213.432   2000000.000   2000000.0 346213.4316        2  XX  YY
  15     BB_3#SQLD                                       1   TOPO  336171.432   2000000.000   2000000.0 336171.4316        3  XX  YY
  16     BB_4#SQLD                                       1   TOPO  335470.830   2000000.000   2000000.0 335470.8304        4  XX  YY
  17     X258738026#ALMA_RB_07#BB_1#SW-01#FULL_RES     128   TOPO  345533.744     15625.000   2000000.0 346525.9316        1  XX  YY
  18     X258738026#ALMA_RB_07#BB_1#SW-01#CH_AVG         1   TOPO  346518.119   1875000.000   1875000.0 346518.1191        1  XX  YY
  19     X258738026#ALMA_RB_07#BB_2#SW-01#FULL_RES     128   TOPO  345221.244     15625.000   2000000.0 346213.4316        2  XX  YY
  20     X258738026#ALMA_RB_07#BB_2#SW-01#CH_AVG         1   TOPO  346205.619   1875000.000   1875000.0 346205.6191        2  XX  YY
  21     X258738026#ALMA_RB_07#BB_3#SW-01#FULL_RES     128   TOPO  337163.619    -15625.000   2000000.0 336171.4316        3  XX  YY
  22     X258738026#ALMA_RB_07#BB_3#SW-01#CH_AVG         1   TOPO  336163.619   1875000.000   1875000.0 336163.6191        3  XX  YY
  23     X258738026#ALMA_RB_07#BB_4#SW-01#FULL_RES     128   TOPO  336463.018    -15625.000   2000000.0 335470.8304        4  XX  YY
  24     X258738026#ALMA_RB_07#BB_4#SW-01#CH_AVG         1   TOPO  335463.018   1875000.000   1875000.0 335463.0179        4  XX  YY
  25     X258738026#ALMA_RB_07#BB_1#SW-01#FULL_RES    1920   TOPO  345541.679       244.141    468750.0 345775.9316        1  XX  YY
  26     X258738026#ALMA_RB_07#BB_1#SW-01#CH_AVG         1   TOPO  345775.810    468750.000    468750.0 345775.8096        1  XX  YY
  27     X258738026#ALMA_RB_07#BB_1#SW-02#FULL_RES    1920   TOPO  347035.606       244.141    468750.0 347269.8586        1  XX  YY
  28     X258738026#ALMA_RB_07#BB_1#SW-02#CH_AVG         1   TOPO  347269.737    468750.000    468750.0 347269.7366        1  XX  YY
  29     X258738026#ALMA_RB_07#BB_2#SW-01#FULL_RES    1920   TOPO  345251.121        61.035    117187.5 345309.6841        2  XX  YY
  30     X258738026#ALMA_RB_07#BB_2#SW-01#CH_AVG         1   TOPO  345309.654    117187.500    117187.5 345309.6536        2  XX  YY
  31     X258738026#ALMA_RB_07#BB_2#SW-02#FULL_RES    1920   TOPO  346909.629        61.035    117187.5 346968.1924        2  XX  YY
  32     X258738026#ALMA_RB_07#BB_2#SW-02#CH_AVG         1   TOPO  346968.162    117187.500    117187.5 346968.1619        2  XX  YY
  33     X258738026#ALMA_RB_07#BB_3#SW-01#FULL_RES    3840   TOPO  337149.001       -61.035    234375.0 337031.8442        3  XX  YY
  34     X258738026#ALMA_RB_07#BB_3#SW-01#CH_AVG         1   TOPO  337031.814    234375.000    234375.0 337031.8137        3  XX  YY
  35     X258738026#ALMA_RB_07#BB_4#SW-01#FULL_RES    1920   TOPO  336407.842      -976.562   1875000.0 335470.8304        4  XX  YY
  36     X258738026#ALMA_RB_07#BB_4#SW-01#CH_AVG         1   TOPO  335470.586   1875000.000   1875000.0 335470.5863        4  XX  YY


split(vis=sourcems,field=source,spw='33',datacolumn='all',keepflags=False,outputvis=sourcename+'.C17O.ms')
split(vis=sourcems,field=source,spw='31',datacolumn='all',keepflags=False,outputvis=sourcename+'.H13COp.ms')
split(vis=sourcems,field=source,spw='29',datacolumn='all',keepflags=False,outputvis=sourcename+'.H13CN.ms')
split(vis=sourcems,field=source,spw='25',datacolumn='all',keepflags=False,outputvis=sourcename+'.CO.ms')
split(vis=sourcems,field=source,spw='27',datacolumn='all',keepflags=False,outputvis=sourcename+'.SiO.ms')
split(vis=sourcems,field=source,spw='35',datacolumn='all',keepflags=False,outputvis=sourcename+'.335.5GHz.ms')

plotms(vis='L1448IRS3B.ms',xaxis='channel',spw='33',avgtime='1e8',averagedata=True,avgscan=True) 
#check each spectral window

#C17O 337.06121 spw 33
#H13CO+ 346.99834 spw 31
#H13CN 345.33977 spw 29
#CO 345.79599 - spw 2
#SiO 347.330  - spw 27

flagdata(vis=sourcems,mode='manual',spw='33:1820~2020')
flagdata(vis=sourcems,mode='manual',spw='31:850~1050')
flagdata(vis=sourcems,mode='manual',spw='29:580~1050')
flagdata(vis=sourcems,mode='manual',spw='25:750~1300')
flagdata(vis=sourcems,mode='manual',spw='27:850~1050')
flagdata(vis=sourcems,mode='manual',spw='35:850~1050')
#all of 5 and 11 look good 
#flagdata(vis=sourcems,mode='manual',spw=5:xxx~xx,11:xxx~xxx')

split(vis=sourcems,width=16,datacolumn='all',field='4',spw='33,31,29,25,27,35',outputvis='L1448IRS3B_cont.ms',keepflags=False)

#source parameters
# ------------------

field='0' # science field(s). For a mosaic, select all mosaic fields. DO NOT LEAVE BLANK ('') OR YOU WILL TRIGGER A BUG IN CLEAN THAT WILL PUT THE WRONG COORDINATE SYSTEM ON YOUR FINAL IMAGE.

# image parameters.
# ----------------
imagermode='csclean' 
cell='0.05arcsec' # cell size for imaging.
imsize = [512,512] # size of image in pixels.
outframe='lsrk' # velocity reference frame. See science goals.
veltype='radio' # velocity type. See note below.
weighting = 'briggs'
robust=0.5
niter=1000
threshold = '0.0mJy'


contvis = 'L1448IRS3B_cont.ms'         
contimagename = sourcename+'_cont_image'
refant = 'DV03,DV15' 
spwmap = [0,0,0,0,0,0,0]


os.system('rm -rf ' +contimagename + '_nosc*')
clean(vis=contvis,
      imagename=contimagename + '_nosc',
      field=field,
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=robust,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      usescratch=True, 
      imagermode=imagermode)


#100 iterations
os.system('rm -rf ' +contimagename + '_p0*')
clean(vis=contvis,
      imagename=contimagename + '_p0',
      field=field,
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=robust,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      usescratch=True, 
      imagermode=imagermode)


# per scan solution
rmtables('pcal1')
gaincal(vis=contvis,
        caltable='pcal1',
        field=field,
        gaintype='T',
        refant=refant, 
        calmode='p',
        combine='spw', 
        solint='inf',
        minsnr=3.0,
        uvrange='>100klambda',
        minblperant=6)

# Check the solution
plotcal(caltable='pcal1',
        xaxis='time',
        yaxis='phase',
        iteration='antenna',
        subplot=421,
	figfile='pcal1.png',
        plotrange=[0,0,-180,180],timerange='2016/12/19~2016/12/20')
# 19-Dec-2016/02:27:35.4   to   19-Dec-2016/03:22:31.0
# made a backup contvis file
# apply the calibration to the data for next round of imaging
applycal(vis=contvis,
         field=field,
         spwmap=spwmap, 
         gaintable=['pcal1'],
         gainfield='',
         calwt=F, 
         flagbackup=F,
         applymode='calonly')

# clean deeper
#for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
#    rmtables(contimagename + '_p1'+ ext)
#240 iterations
os.system('rm -rf ' +contimagename + '_p1*')
clean(vis=contvis,
      field=field,
      imagename=contimagename + '_p1',
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=robust,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      usescratch=True, # needed for 4.3 and 4.4 (and maybe 4.5)
      imagermode=imagermode)


# Note number of iterations performed.

 #shorter solution
rmtables('pcal2')
gaincal(vis=contvis,
        field=field,
        caltable='pcal2',
        gaintype='T',
        refant=refant, 
        calmode='p',
        combine='spw', 
        solint='30.25', # solint=30.25s gets you five 12m integrations, while solint=50.5s gets you five 7m integration
        minsnr=3.0,
	uvrange='>100klambda',
        minblperant=6)

# Check the solution
plotcal(caltable='pcal2',
        xaxis='time',
        yaxis='phase',
        iteration='antenna',
        subplot=421,
	figfile='pcal2.png',
        plotrange=[0,0,-180,180],timerange='2016/09/20~2016/10/02')

os.system('rm -rf ' +contimagename + '_p1_full*')
clean(vis=contvis,
      field=field,
      imagename=contimagename + '_p1_full',
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=robust,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      imagermode=imagermode)



# apply the calibration to the data for next round of imaging
applycal(vis=contvis,
         field=field,
         spwmap=spwmap, 
         gaintable=['pcal2'],
         gainfield='',
         calwt=F, 
         flagbackup=F,
         applymode='calonly')


# clean deeper
#for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
#    rmtables(contimagename + '_p2'+ ext)
#600 iterations
os.system('rm -rf ' +contimagename + '_p2*')
clean(vis=contvis,
      imagename=contimagename + '_p2',
      field=field,
      phasecenter=phasecenter, # uncomment if mosaic.            
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=robust,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      usescratch=True, # needed for 4.3 and 4.4 (and maybe 4.5)
      imagermode=imagermode)

#shorter solution
rmtables('pcal3')
gaincal(vis=contvis,
        field=field,
        caltable='pcal3',
        gaintype='T',
        refant=refant, 
        calmode='p',
        combine='spw', 
        solint='12.1', # solint=30.25s gets you five 12m integrations, while solint=50.5s gets you five 7m integration
        minsnr=3.0,uvrange='>100klambda',
        minblperant=6)

# Check the solution
plotcal(caltable='pcal3',
        xaxis='time',
        yaxis='phase',
        iteration='antenna',
        subplot=421,
	figfile='pcal3.png',
        plotrange=[0,0,-180,180],timerange='2016/09/20~2016/10/02')


os.system('rm -rf ' +contimagename + '_p2_full*')
clean(vis=contvis,
      imagename=contimagename + '_p2_full',
      field=field,
      phasecenter=phasecenter, # uncomment if mosaic.            
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=robust,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      imagermode=imagermode)

# apply the calibration to the data for next round of imaging
applycal(vis=contvis,
         field=field,
         spwmap=spwmap, 
         gaintable=['pcal3'],
         gainfield='',
         calwt=F, 
         flagbackup=F,
         applymode='')


# clean deeper
#for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
#    rmtables(contimagename + '_p2'+ ext)
#1500 iterations
os.system('rm -rf ' +contimagename + '_p3*')
clean(vis=contvis,
      imagename=contimagename + '_p3',
      field=field,
      phasecenter=phasecenter, # uncomment if mosaic.            
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=robust,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      usescratch=True, # needed for 4.3 and 4.4 (and maybe 4.5)
      imagermode=imagermode)

#shorter solution
# do the amplitude self-calibration.

rmtables('apcal')
gaincal(vis=contvis,
        field=field,
        caltable='apcal',
        gaintype='T',
        refant=refant,
        calmode='ap',
        combine='spw',
        solint='inf',
        minsnr=3.0,
        minblperant=4,
        gaintable='pcal3',
        spwmap=spwmap,
        uvrange='>100klambda',
        solnorm=True)

plotcal(caltable='apcal',
        xaxis='time',
        yaxis='amp',
        timerange='',
        iteration='antenna',
        subplot=421,
        plotrange=[0,0,0.2,1.8],timerange='2016/09/20~2016/10/02')

plotcal(caltable='apcal',
        xaxis='time',
        yaxis='amp',
        iteration='antenna',
        subplot=421,
        plotrange=[0,0,0.2,1.8],timerange='2016/10/03~2016/10/05')

applycal(vis=contvis,
         spwmap=[spwmap,spwmap], # select which spws to apply the solutions for each table
         field=field,
         gaintable=['pcal3','apcal'],
         gainfield='',
         calwt=F,
         flagbackup=F,
         applymode='calonly')

# Make amplitude and phase self-calibrated image.
#for ext in ['.flux','.image','.mask','.model','.pbcor','.psf','.residual','.flux.pbcoverage']:
#    rmtables(contimagename + '_ap'+ ext)
os.system('rm -rf ' +contimagename + '_ap*')
clean(vis=contvis,
      imagename=contimagename + '_ap',
      field=field, 
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=robust,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      imagermode=imagermode,uvrange='>50klambda')

cell='0.01arcsec' # cell size for imaging.
imsize = [2048,2048] # size of image in pixels.
outframe='lsrk' # velocity reference frame. See science goals.

os.system('rm -rf ' +contimagename + '_robust-0.5*')
clean(vis=contvis,
      imagename=contimagename + '_robust-0.5',
      field=field, 
      mode='mfs',
      psfmode='clark',
      imsize = imsize, 
      cell= cell, 
      weighting = weighting, 
      robust=-0.5,
      niter=niter, 
      threshold=threshold, 
      interactive=True,
      imagermode=imagermode,uvrange='>50klambda')

# Save results of self-cal in a new ms
os.system('rm -rf '+contvis+'.selfcal')
split(vis=contvis,
      outputvis=contvis+'.selfcal',
      datacolumn='corrected')
os.system('rm -rf '+contvis)


 cd ../

applycal(vis=contvis,
         spwmap=[spwmap,spwmap], # select which spws to apply the solutions for each table
         field=field,
         gaintable=['cont-only/pcal3','cont-only/apcal'],
         gainfield='',
         calwt=F,
         flagbackup=F,
         applymode='calonly')

os.system('rm -rf '+contvis+'.selfcal')
split(vis=contvis,
      outputvis=contvis+'.selfcal',
      datacolumn='corrected')
os.system('rm -rf '+contvis)


