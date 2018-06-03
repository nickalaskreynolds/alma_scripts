source='L1448IRS3B'


#########################################################
# Apply continuum self-calibration to line data [OPTIONAL]


# save original flags in case you don't like the self-cal
#flagmanager(vis=linevis,mode='save',versionname='before_selfcal',merge='replace')
sourcevis=source+'.C17O.ms'
spwmap_line = [0,0] # Mapping self-calibration solution to the individual line spectral windows.
applycal(vis=sourcevis,
         spwmap=[spwmap_line,spwmap_line], # select which spws to apply the solutions for each table
         gaintable=['cont-only/apcal','cont-only/pcal3'],
         gainfield='',
         calwt=F,
         flagbackup=F,
         applymode='calonly')

# Save results of self-cal in a new ms and reset the image name.
split(vis=sourcevis,
      outputvis=sourcevis+'.selfcal',
      datacolumn='corrected')



linevis=source+'.C17O.ms.selfcal'
fitspw = '0:0~1800,0:2100~3840,1:0~1800,1:2100~3840' # line-free channel for fitting continuum
linespw = '0,1'
uvcontsub(vis=linevis,
          spw=linespw, # spw to do continuum subtraction on
          fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
          combine='spw', 
          solint='int',
          fitorder=1,
          want_cont=False) # This value should not be changed.
linevis = linevis+'.contsub'

os.system('rm -rf '+source+'.C17O.ms.selfcal')
os.system('rm -rf '+source+'.C17O.ms')


sourcevis=source+'.H13COp.ms'
spwmap_line = [0,0] # Mapping self-calibration solution to the individual line spectral windows.
applycal(vis=sourcevis,
         spwmap=[spwmap_line,spwmap_line], # select which spws to apply the solutions for each table
         gaintable=['cont-only/apcal','cont-only/pcal3'],
         gainfield='',
         calwt=F,
         flagbackup=F,
         applymode='calonly')

# Save results of self-cal in a new ms and reset the image name.
split(vis=sourcevis,
      outputvis=sourcevis+'.selfcal',
      datacolumn='corrected')

plotms(vis=sourcevis+'.selfcal',averagedata=True,avgtime='1e8',avgscan=True,xaxis='channel',spw='0,1,2',coloraxis='spw')


linevis=source+'.H13COp.ms.selfcal'
fitspw = '0:0~800;1100~1919,1:0~800;1100~1919' # line-free channel for fitting continuum
linespw = '0,1'
uvcontsub(vis=linevis,
          spw=linespw, # spw to do continuum subtraction on
          fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
          combine='spw', 
          solint='int',
          fitorder=1,
          want_cont=False) # This value should not be changed.
linevis = linevis+'.contsub'

os.system('rm -rf '+source+'.H13COp.ms.selfcal')
os.system('rm -rf '+source+'.H13COp.ms')


sourcevis=source+'.H13CN.ms'
spwmap_line = [0,0] # Mapping self-calibration solution to the individual line spectral windows.
applycal(vis=sourcevis,
         spwmap=[spwmap_line,spwmap_line], # select which spws to apply the solutions for each table
         gaintable=['cont-only/apcal','cont-only/pcal3'],
         gainfield='',
         calwt=F,
         flagbackup=F,
         applymode='calonly')

# Save results of self-cal in a new ms and reset the image name.
split(vis=sourcevis,
      outputvis=sourcevis+'.selfcal',
      datacolumn='corrected')

plotms(vis=sourcevis+'.selfcal',averagedata=True,avgtime='1e8',avgscan=True,xaxis='channel',spw='0,1',coloraxis='spw')


linevis=source+'.H13CN.ms.selfcal'
fitspw = '0:1100~1919,1:1100~1919' # line-free channel for fitting continuum
linespw = '0,1'
uvcontsub(vis=linevis,
          spw=linespw, # spw to do continuum subtraction on
          fitspw=fitspw, # select spws to fit continuum. exclude regions with strong lines.
          combine='spw', 
          solint='int',
          fitorder=1,
          want_cont=False) # This value should not be changed.
linevis = linevis+'.contsub'

os.system('rm -rf '+source+'.H13COp.ms.selfcal')
os.system('rm -rf '+source+'.H13COp.ms')

