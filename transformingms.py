# recentering
# L1448IRS3B.H13CN.ms.selfcal.contsub.concat.wide
# L1448IRS3B.C17O.ms.selfcal.contsub.concat.triplet
# L1448IRS3B.cont.ms.apselfcal.concat.wide
# L1448IRS3B.cont.ms.apselfcal.concat.triplet

default(mstransform)
vis                = 'L1448IRS3B.C17O.ms.selfcal.contsub.concat'
keepflags          =  False
timebin            =  "10s"
timeaverage        =  True
mode               =  "velocity"
nchan              =  55
start              =  '-0.5km/s'
width              =  '0.163km/s'
combinespws        =  True
regridms           =  True
datacolumn         =  "data"
restfreq           =  "337061.104MHz"
outframe           =  "lsrk"
outputvis          =  'L1448IRS3B.C17O.ms.selfcal.contsub.concat.lsrk.tavg.chavg'
go()

default(mstransform)
vis                = 'L1448IRS3B.H13CN.ms.selfcal.contsub.concat'
keepflags          =  False
timebin            =  "30s"
timeaverage        =  True
mode               =  "velocity"
nchan              =  50
start              =  '2km/s'
width              =  '0.2651km/s'
combinespws        =  True
regridms           =  True
datacolumn         =  "data"
restfreq           =  "345.339756GHz"
outframe           =  "lsrk"
outputvis          = 'L1448IRS3B.H13CN.ms.selfcal.contsub.concat.lsrk.tavg.chavg'
go()
