import os

source = 'L1448IRS3B'
moments = [-1,0]
tapers = ['1500k','1000k'] # in lambda
linenames=['12CO','C17O','H13CN','H13COp']
#         r        b sequence
chanrb=['38~186','207~321','57~90','90~115','25~52','53~99','28~47','48~67']

print("Do not cancel this task. Let it finish, high chance of corruption if cancelled.")

counter = 0
countern = 0
for j in linenames:
    for i in tapers:
        counter = countern
        for k in range(2):
            lineimagename = source+'_'+j+'_image_taper' + i
            if ((counter%2)==0):
                outname=lineimagename+"_r"
                noutname = outname+'.image'
            elif ((counter%2)==1):
                outname=lineimagename+"_b"
                noutname = outname+'.image'
            os.system("rm -rf " + outname + "*integrated " + outname + "*average " + outname + "*fits")
            immoments(imagename = lineimagename + '.image', moments = moments,axis='spectral',mask = lineimagename + '.mask',chans=chanrb[counter],outfile=noutname)
            exportfits(imagename = noutname+'.integrated', fitsimage=outname+'.fits')
            counter+=1
    countern+=2
