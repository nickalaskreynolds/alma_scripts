source='L1448IRS3B'
lines=['SiO','H13COp','H13CN','cont','12CO','C17O','335.5GHz']
#lines = ['335.5GHz']
lines = ['C17O']
for num in range(len(lines)):
    print('Working on: {}'.format(lines[num]))

    if 'cont' in lines[num]:
        sourcems = source + '.' + lines[num] + '.ms.selfcal'
        sourceconcatms = source + '.' + lines[num] + '.hlr' + '.ms.selfcal'
        os.system('rm -rf {}'.format(sourceconcatms))
        concat(vis=['../ALMA-Per33-hr/'+sourcems,'../ALMA-Per33-lr/'+sourcems],concatvis=sourceconcatms)

        output=source+'.'+lines[num]+'.ms.apselfcal.concat'
        os.system('rm -rf {}'.format(output))
        default(split)
        vis=sourceconcatms
        outputvis=output
        go()
    else:
        sourcems = source + '.' + lines[num] + '.ms.selfcal.contsub'
        sourceconcatms = source + '.' + lines[num] + '.hlr' + '.ms.selfcal.contsub'
        os.system('rm -rf {}'.format(sourceconcatms))
        concat(vis=['../ALMA-Per33-hr/'+sourcems,'../ALMA-Per33-lr/'+sourcems],concatvis=sourceconcatms)

        output=source+'.'+lines[num]+'.ms.apselfcal.contsub.concat'
        os.system('rm -rf {}'.format(output))
        default(split)
        vis=sourceconcatms
        outputvis=output
        go()


######
