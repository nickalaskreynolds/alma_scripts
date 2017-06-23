# Run commands in order, otherwise might have conflicts unresolved
source='Per33_L1448_IRS3B'
sourcename='L1448IRS3B'
sourcems='L1448IRS3B.ms'
contms='L1448IRS3B_cont.ms'
origms1='uid___A002_Xb8fee5_X3055.ms.split.cal'
origms2='uid___A002_Xb8e115_X3798.ms.split.cal'
concatorig='uuid__concat.ms'
#plotms(vis='uid___A002_Xbbe66a_X1954.ms',field='4',xaxis='time',yaxis='amplitude',coloraxis='field',spw='20~40',averagedata=True,avgscan=True)
#depending if original ms file has calibrated data or not will determine what <datacolumn> modifier you will use.
#Preserving the original file
os.system('rm -rf ' +concatorig + '*')
concat(vis=[origms1,origms2],concatvis=concatorig)
os.system('rm -rf ' +sourcems + '*')
split(vis=concatorig,datacolumn='all',field='4',outputvis=sourcems,keepflags=False)

#Listobs output
2017-01-25 16:21:19 INFO listobs	Antennas: 46:
2017-01-25 16:21:19 INFO listobs	  ID   Name  Station   Diam.    Long.         Lat.                Offset from array center (m)                ITRF Geocentric coordinates (m)        
2017-01-25 16:21:19 INFO listobs	                                                                     East         North     Elevation               x               y               z
2017-01-26 21:03:24 INFO listobs	  SpwID  Name                                      #Chans   Frame   Ch0(MHz)  ChanWid(kHz)  TotBW(kHz) CtrFreq(MHz) BBC Num  Corrs  
2017-01-26 21:03:24 INFO listobs	  0      X276863869#ALMA_RB_07#BB_1#SW-01#FULL_RES   1920   TOPO  345778.059       244.141    468750.0 346012.3123        1  XX  YY
2017-01-26 21:03:24 INFO listobs	  1      X276863869#ALMA_RB_07#BB_1#SW-02#FULL_RES   1920   TOPO  346778.059       244.141    468750.0 347012.3123        1  XX  YY
2017-01-26 21:03:24 INFO listobs	  2      X276863869#ALMA_RB_07#BB_2#SW-01#FULL_RES   1920   TOPO  345293.511        61.035    117187.5 345352.0738        2  XX  YY
2017-01-26 21:03:24 INFO listobs	  3      X276863869#ALMA_RB_07#BB_2#SW-02#FULL_RES   1920   TOPO  346952.019        61.035    117187.5 347010.5821        2  XX  YY
2017-01-26 21:03:24 INFO listobs	  4      X276863869#ALMA_RB_07#BB_3#SW-01#FULL_RES   3840   TOPO  337190.170       -61.035    234375.0 337073.0133        3  XX  YY
2017-01-26 21:03:24 INFO listobs	  5      X276863869#ALMA_RB_07#BB_4#SW-01#FULL_RES   1920   TOPO  336448.906      -976.562   1875000.0 335511.8943        4  XX  YY
2017-01-26 21:03:24 INFO listobs	  6      X276863869#ALMA_RB_07#BB_1#SW-01#FULL_RES   1920   TOPO  345776.946       244.141    468750.0 346011.1993        1  XX  YY
2017-01-26 21:03:24 INFO listobs	  7      X276863869#ALMA_RB_07#BB_1#SW-02#FULL_RES   1920   TOPO  346776.946       244.141    468750.0 347011.1993        1  XX  YY
2017-01-26 21:03:24 INFO listobs	  8      X276863869#ALMA_RB_07#BB_2#SW-01#FULL_RES   1920   TOPO  345292.398        61.035    117187.5 345350.9609        2  XX  YY
2017-01-26 21:03:24 INFO listobs	  9      X276863869#ALMA_RB_07#BB_2#SW-02#FULL_RES   1920   TOPO  346950.906        61.035    117187.5 347009.4692        2  XX  YY
2017-01-26 21:03:24 INFO listobs	  10     X276863869#ALMA_RB_07#BB_3#SW-01#FULL_RES   3840   TOPO  337189.088       -61.035    234375.0 337071.9309        3  XX  YY
2017-01-26 21:03:24 INFO listobs	  11     X276863869#ALMA_RB_07#BB_4#SW-01#FULL_RES   1920   TOPO  336447.841      -976.562   1875000.0 335510.8297        4  XX  YY



os.system('rm -rf ' +sourcems + '.*' + '*.ms')
split(vis=sourcems,field=source,spw='4,10',datacolumn='all',keepflags=False,outputvis=sourcename+'.C17O.ms')
split(vis=sourcems,field=source,spw='3,9',datacolumn='all',keepflags=False,outputvis=sourcename+'.H13COp.ms')
split(vis=sourcems,field=source,spw='2,8',datacolumn='all',keepflags=False,outputvis=sourcename+'.H13CN.ms')
split(vis=sourcems,field=source,spw='0,6',datacolumn='all',keepflags=False,outputvis=sourcename+'.CO.ms')
split(vis=sourcems,field=source,spw='1,7',datacolumn='all',keepflags=False,outputvis=sourcename+'.SiO.ms')
split(vis=sourcems,field=source,spw='5,11',datacolumn='all',keepflags=False,outputvis=sourcename+'.335.5GHz.ms')

#plotms(vis='L1448IRS3B.ms',xaxis='channel',spw='3',avgtime='1e8',averagedata=True,avgscan=True) 
#check each spectral window

#C17O 337.06121 spw  4 10
#H13CO+ 346.99834 spw 3 9
#H13CN 345.33977 spw 2 8
#CO 345.79599 - spw 0 6
#SiO 347.330  - spw 1 7

flagdata(vis=sourcems,mode='manual',spw='0:0~400,6:10~400')
flagdata(vis=sourcems,mode='manual',spw='3:900~1000,9:900~1000')
flagdata(vis=sourcems,mode='manual',spw='1:925~975,7:925~975')
flagdata(vis=sourcems,mode='manual',spw='2:875~1000,8:875~1000')
flagdata(vis=sourcems,mode='manual',spw='4:1850~2000,11:1850~2000')

#flagdata(vis=sourcems,mode='manual',spw=5:xxx~xx,11:xxx~xxx')

split(vis=sourcems,width=16,datacolumn='all',field='0',spw='0~11',outputvis='L1448IRS3B_cont.ms',keepflags=False)

#source parameters
# ------------------

field='0' # science field(s). For a mosaic, select all mosaic fields. DO NOT LEAVE BLANK ('') OR YOU WILL TRIGGER A BUG IN CLEAN THAT WILL PUT THE WRONG COORDINATE SYSTEM ON YOUR FINAL IMAGE.

# image parameters.
# ----------------
imagermode='csclean' 
cell='0.01arcsec' # cell size for imaging.
imsize = [2048,2048] # size of image in pixels.
outframe='lsrk' # velocity reference frame. See science goals.
veltype='radio' # velocity type. See note below.
weighting = 'briggs'
robust=0.5
niter=1000
threshold = '0.0mJy'


contvis = 'L1448IRS3B_cont.ms'         
contimagename = sourcename+'_cont_image'
refant = 'DA49,DV18,DV09,DA46' 
spwmap = [0,0,0,0,0,0,0,0,0,0,0,0]

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
        minblperant=6)

# Check the solution
plotcal(caltable='pcal1',
        xaxis='time',
        yaxis='phase',
        iteration='antenna',
        subplot=421,
	figfile='pcal1.png',
        plotrange=[0,0,-180,180])

 
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
#100 iterations
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

# shorter solution
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
        minblperant=6)

# Check the solution
plotcal(caltable='pcal2',
        xaxis='time',
        yaxis='phase',
        iteration='antenna',
        subplot=421,
	figfile='pcal2.png',
        plotrange=[0,0,-180,180])

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
#300 iterations
os.system('rm -rf ' +contimagename + '_p2*')
clean(vis=contvis,
      imagename=contimagename + '_p2',
      field=field,
      #phasecenter=phasecenter, # uncomment if mosaic.            
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
        minsnr=3.0,
        minblperant=6)

# Check the solution
plotcal(caltable='pcal3',
        xaxis='time',
        yaxis='phase',
        iteration='antenna',
        subplot=421,
	figfile='pcal3.png',
        plotrange=[0,0,-180,180])


os.system('rm -rf ' +contimagename + '_p2_full*')
clean(vis=contvis,
      imagename=contimagename + '_p2_full',
      field=field,
      #phasecenter=phasecenter, # uncomment if mosaic.            
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
#3000 iterations
os.system('rm -rf ' +contimagename + '_p3*')
clean(vis=contvis,
      imagename=contimagename + '_p3',
      field=field,
      #phasecenter=phasecenter, # uncomment if mosaic.            
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
        solnorm=True)

plotcal(caltable='apcal',
        xaxis='time',
        yaxis='amp',
        timerange='',
        iteration='antenna',
        subplot=421,
        plotrange=[0,0,0.2,1.8])
####################################################################################################################current progress
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
# 2500 iterations
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
      imagermode=imagermode)

cell='0.05arcsec' # cell size for imaging.
imsize = [512,512] # size of image in pixels.
outframe='lsrk' # velocity reference frame. See science goals.

# 220 iterations
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
      imagermode=imagermode)

# Save results of self-cal in a new ms
os.system('rm -rf '+contvis+'.selfcal')
split(vis=contvis,
      outputvis=contvis+'.selfcal',
      datacolumn='corrected')


applycal(vis=contvis,
         spwmap=[spwmap,spwmap], # select which spws to apply the solutions for each table
         field=field,
         gaintable=['pcal3','apcal'],
         gainfield='',
         calwt=F,
         flagbackup=F,
         applymode='calonly')

os.system('rm -rf '+contvis+'.apselfcal')
split(vis=contvis,
      outputvis=contvis+'.apselfcal',
      datacolumn='corrected')


