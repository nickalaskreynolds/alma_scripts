source='L1448IRS3B'


#########################################################
# Apply continuum self-calibration to line data [OPTIONAL]


# save original flags in case you don't like the self-cal
#flagmanager(vis=linevis,mode='save',versionname='before_selfcal',merge='replace')
sourcevis=source+'.C17O.ms'
restfreq='337.06121GHz'
spwmap_line = [0] # Mapping self-calibration solution to the individual line spectral windows.
applycal(vis=sourcevis,
         spwmap=[spwmap_line], # select which spws to apply the solutions for each table
         gaintable=['apcal','pcal3'],
         gainfield='',
         calwt=F,
         flagbackup=F,
         applymode='calonly')

# Save results of self-cal in a new ms and reset the image name.
os.system("rm -rvf " + sourcevis+'.selfcal')
split(vis=sourcevis,
      outputvis=sourcevis+'.selfcal',
      datacolumn='corrected')

# #plotms(vis='L1448IRS3B.C17O.ms',xaxis='channel',restfreq=restfreq,spw='0',avgtime='1e8',averagedata=True,avgscan=True)
#plotms(vis=sourcevis+'.selfcal',averagedata=True,avgtime='1e8',avgscan=True,xaxis='channel',restfreq=restfreq,spw='0,1')
#C17O 337.06121 spw  4 10
#H13CO+ 346.99834 spw 3 9
#H13CN 345.33977 spw 2 8
#CO 345.79599 - spw 0 6
#flagdata(vis=sourcems,mode='manual',spw='0:0~400,6:10~400')
#flagdata(vis=sourcems,mode='manual',spw='3:900~1000,9:900~1000')
#flagdata(vis=sourcems,mode='manual',spw='2:875~1000,8:875~1000')
#flagdata(vis=sourcems,mode='manual',spw='4:1850~2000,11:1850~2000')

linevis=source+'.C17O.ms.selfcal'
fitspw = '0:1820~2020,1:1820~2020' # line-free channel for fitting continuum
linespw = '0,1'
os.system("rm -rvf " + linevis+'.contsub')
uvcontsub(vis=linevis,
          spw=linespw, # spw to do continuum subtraction on
          fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
          combine='spw', 
          solint='int',
          fitorder=1,
          want_cont=False,excludechans=True) # This value should not be changed.
linevis = linevis+'.contsub'

#os.system('rm -rvf '+source+'.C17O.ms.selfcal')
#os.system('rm -rvf '+source+'.C17O.ms')

##########################
sourcevis=source+'.H13COp.ms'
spwmap_line = [0] # Mapping self-calibration solution to the individual line spectral windows.
applycal(vis=sourcevis,
         spwmap=[spwmap_line], # select which spws to apply the solutions for each table
         gaintable=['apcal','pcal3'],
         gainfield='',
         calwt=F,
         flagbackup=F,
         applymode='calonly')

# Save results of self-cal in a new ms and reset the image name.
os.system("rm -rvf " + sourcevis+'.selfcal')
split(vis=sourcevis,
      outputvis=sourcevis+'.selfcal',
      datacolumn='corrected')

# #plotms(vis=sourcevis+'.selfcal',averagedata=True,avgtime='1e8',avgscan=True,xaxis='channel',restfreq=restfreq,spw='0,1,2',coloraxis='spw')
#plotms(vis=sourcevis+'.selfcal',averagedata=True,avgtime='1e8',avgscan=True,xaxis='channel',spw='0,1')

linevis=source+'.H13COp.ms.selfcal'
fitspw = '0:875~1000,1:875~1000' # line-free channel for fitting continuum
linespw = '0,1'
os.system("rm -rvf " + linevis+'.contsub')
uvcontsub(vis=linevis,
          spw=linespw, # spw to do continuum subtraction on
          fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
          combine='spw', 
          solint='int',
          fitorder=1,
          want_cont=False,excludechans=True)# This value should not be changed.
linevis = linevis+'.contsub'
#os.system('rm -rvf '+source+'.H13COp.ms.selfcal')
#os.system('rm -rvf '+source+'.H13COp.ms')

##########################
sourcevis=source+'.H13CN.ms'
spwmap_line = [0] # Mapping self-calibration solution to the individual line spectral windows.
applycal(vis=sourcevis,
         spwmap=[spwmap_line,spwmap_line], # select which spws to apply the solutions for each table
         gaintable=['apcal','pcal3'],
         gainfield='',
         calwt=F,
         flagbackup=F,
         applymode='calonly')

# Save results of self-cal in a new ms and reset the image name.
#os.system('rm -rvf ' + sourcevis + '.selfcal')
os.system("rm -rvf " + sourcevis+'.selfcal')
split(vis=sourcevis,
      outputvis=sourcevis+'.selfcal',
      datacolumn='corrected')

##plotms(vis=sourcevis+'.selfcal',averagedata=True,avgtime='1e8',avgscan=True,xaxis='channel',restfreq=restfreq,spw='0,1')

linevis=source+'.H13CN.ms.selfcal'
fitspw = '0:875~1000,1:875~1000' # line-free channel for fitting continuum
linespw = '0,1'
os.system("rm -rvf " + linevis+'.contsub')
uvcontsub(vis=linevis,
          spw=linespw, # spw to do continuum subtraction on
          fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
          combine='spw', 
          solint='int',
          fitorder=1,
          want_cont=False,excludechans=True)# This value should not be changed.
linevis = linevis+'.contsub'

#os.system('rm -rvf '+source+'.H13CN.ms.selfcal')
#os.system('rm -rvf '+source+'.H13CN.ms')


##########################========
sourcevis=source+'.12CO.ms'
tsourcevis=source+'.12CO.ms'
spwmap_line = [0] # Mapping self-calibration solution to the individual line spectral windows.
applycal(vis=tsourcevis,
         spwmap=[spwmap_line,spwmap_line], # select which spws to apply the solutions for each table
         gaintable=['apcal','pcal3'],
         gainfield='',
         calwt=F,
         flagbackup=F,
         applymode='calonly')

# Save results of self-cal in a new ms and reset the image name.
os.system("rm -rvf " + sourcevis+'.selfcal')
split(vis=tsourcevis,
      outputvis=sourcevis+'.selfcal',
      datacolumn='corrected')

#plotms(vis=sourcevis+'.selfcal',averagedata=True,avgtime='1e8',avgscan=True,xaxis='channel',restfreq=restfreq,spw='0,1')

linevis=source+'.12CO.ms.selfcal'
fitspw = '0:0~400,1:0~400' # line-free channel for fitting continuum
linespw = '0,1'
os.system("rm -rvf " + linevis+'.contsub')
uvcontsub(vis=linevis,
          spw=linespw, # spw to do continuum subtraction on
          fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
          combine='spw', 
          solint='int',
          fitorder=1,
          want_cont=False,excludechans=True) # This value should not be changed.
linevis = linevis+'.contsub'

#os.system('rm -rvf '+source+'.CO.ms.selfcal')
#os.system('rm -rvf '+source+'.CO.ms')


##########################
sourcevis=source+'.335.5GHz.ms'
spwmap_line = [0] # Mapping self-calibration solution to the individual line spectral windows.
applycal(vis=sourcevis,
         spwmap=[spwmap_line,spwmap_line], # select which spws to apply the solutions for each table
         gaintable=['apcal','pcal3'],
         gainfield='',
         calwt=F,
         flagbackup=F,
         applymode='calonly')

# Save results of self-cal in a new ms and reset the image name.
os.system("rm -rvf " + sourcevis+'.selfcal')
split(vis=sourcevis,
      outputvis=sourcevis+'.selfcal',
      datacolumn='corrected')

#plotms(vis=sourcevis+'.selfcal',averagedata=True,avgtime='1e8',avgscan=True,xaxis='channel',restfreq=restfreq,spw='0,1')

linevis=source+'.335.5GHz.ms.selfcal'
fitspw = '0:0~4000,1:0~4000' # line-free channel for fitting continuum
linespw = '0,1'
os.system("rm -rvf " + linevis+'.contsub')
uvcontsub(vis=linevis,
          spw=linespw, # spw to do continuum subtraction on
          fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
          combine='spw', 
          solint='int',
          fitorder=1,
          want_cont=False,excludechans=True) # This value should not be changed.
linevis = linevis+'.contsub'

#os.system('rm -rvf '+source+'.335.5GHz.ms.selfcal')
#os.system('rm -rvf '+source+'.335.5GHz.ms')