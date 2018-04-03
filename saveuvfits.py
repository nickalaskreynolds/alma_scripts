sourcems = 'L1448IRS3B.cont.ms.apselfcal.concat.triplet'

spws=[['4','10','16'],
['3','9','15'],
['2','8','14'],
['0','6','12'],
['1','7','13'],
['5','11','17']]

names=['C17O',
'H13COp',
'H13CN',
'CO',
'SiO',
'335.5GHz']

dire = 'conttrip'

for x in range(len(spws)):
    for y in spws[x]:
        findir="{}/{}.{}".format(dire,names[x],y)
        if not os.path.isdir(findir):
            os.system('mkdir -p {}'.format(findir))
        else:
            os.system('rm -rf {}/*'.format(findir))
        split(vis=sourcems,spw=y,datacolumn='data',width=999,keepflags=False,outputvis='{}/{}.{}.{}.ms'.format(findir,sourcems,names[x],y))
        exportuvfits(vis='{}/{}.{}.{}.ms'.format(findir,sourcems,names[x],y),fitsfile='{}/{}.{}.{}.uv.fits'.format(findir,sourcems,names[x],y))


sourcems = 'L1448IRS3B.cont.ms.apselfcal.concat.wide'

dire = 'contwide'

for x in range(len(spws)):
    for y in spws[x]:
        findir="{}/{}.{}".format(dire,names[x],y)
        if not os.path.isdir(findir):
            os.system('mkdir -p {}'.format(findir))
        else:
            os.system('rm -rf {}/*'.format(findir))
        split(vis=sourcems,spw=y,datacolumn='data',width=999,keepflags=False,outputvis='{}/{}.{}.{}.ms'.format(findir,sourcems,names[x],y))
        exportuvfits(vis='{}/{}.{}.{}.ms'.format(findir,sourcems,names[x],y),fitsfile='{}/{}.{}.{}.uv.fits'.format(findir,sourcems,names[x],y))


sourcems = 'L1448IRS3B.cont.ms.apselfcal.concat.tertiary'

dire = 'conttert'

for x in range(len(spws)):
    for y in spws[x]:
        findir="{}/{}.{}".format(dire,names[x],y)
        if not os.path.isdir(findir):
            os.system('mkdir -p {}'.format(findir))
        else:
            os.system('rm -rf {}/*'.format(findir))
        split(vis=sourcems,spw=y,datacolumn='data',width=999,keepflags=False,\
            outputvis='{}/{}.{}.{}.ms'.format(findir,sourcems,names[x],y))
        exportuvfits(vis='{}/{}.{}.{}.ms'.format(findir,sourcems,names[x],y),\
            fitsfile='{}/{}.{}.{}.uv.fits'.format(findir,sourcems,names[x],y))



sourcems = 'L1448IRS3B.cont.ms.apselfcal.concat.binary'

dire = 'contbinary'

for x in range(len(spws)):
    for y in spws[x]:
        findir="{}/{}.{}".format(dire,names[x],y)
        if not os.path.isdir(findir):
            os.system('mkdir -p {}'.format(findir))
        else:
            os.system('rm -rf {}/*'.format(findir))
        split(vis=sourcems,spw=y,datacolumn='data',width=999,keepflags=False,\
            outputvis='{}/{}.{}.{}.ms'.format(findir,sourcems,names[x],y))
        exportuvfits(vis='{}/{}.{}.{}.ms'.format(findir,sourcems,names[x],y),\
            fitsfile='{}/{}.{}.{}.uv.fits'.format(findir,sourcems,names[x],y))

