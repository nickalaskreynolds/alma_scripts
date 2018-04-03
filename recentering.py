# recentering
# L1448IRS3B.H13CN.ms.selfcal.contsub.concat.wide
# L1448IRS3B.C17O.ms.selfcal.contsub.concat.triplet
# L1448IRS3B.cont.ms.apselfcal.concat.wide
# L1448IRS3B.cont.ms.apselfcal.concat.triplet
# EXAMPLE.....  phasecenter = 'J2000 19h30m00 -40d00m00'

default(fixvis)
vis='L1448IRS3B.H13CN.ms.selfcal.contsub.concat.lsrk.tavg.chavg'
outputvis='L1448IRS3B.H13CN.ms.selfcal.contsub.concat.wide'
phasecenter='J2000 03h25m36.502 30d45m21.866'
go()

default(fixvis)
vis='L1448IRS3B.C17O.ms.selfcal.contsub.concat.lsrk.tavg.chavg'
outputvis='L1448IRS3B.C17O.ms.selfcal.contsub.concat.triplet'
phasecenter='J2000 03h25m36.324 30d45m14.941'
go()

os.system('rm -rf L1448IRS3B.cont.ms.apselfcal.concat.wide')
default(fixvis)
vis='L1448IRS3B.cont.ms.apselfcal.concat'
outputvis='L1448IRS3B.cont.ms.apselfcal.concat.wide'
phasecenter='J2000 03h25m36.502 30d45m21.866'
go()

os.system('rm -rf L1448IRS3B.cont.ms.apselfcal.concat.binary')
default(fixvis)
vis='L1448IRS3B.cont.ms.apselfcal.concat'
outputvis='L1448IRS3B.cont.ms.apselfcal.concat.binary'
phasecenter='J2000 03h25m36.321 30d45m15.049'
go()

os.system("rm -rf L1448IRS3B.cont.ms.apselfcal.concat.tertiary")
default(fixvis)
vis='L1448IRS3B.cont.ms.apselfcal.concat'
outputvis='L1448IRS3B.cont.ms.apselfcal.concat.tertiary'
phasecenter='J2000 03h25m36.382 30d45m14.715'
go()
