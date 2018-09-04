'''
Pretty plots files from ALMA data sets
Original Author John Tobin
Adapted by Nickalas Reynolds
April 2017
'''

# import libraries
import os
import glob
import aplpy
import sys
import ds9cmap

if len(sys.argv) == 1:
    filesourcename = 'plot_sources.txt'
else:
    filesourcename = sys.argv[1]

while True:
    try:
        answer = input('Which source position is your entire FoV (if NA just [RET],starting from 1): ')
        if answer == '':
            break
        else:
            answer = int(answer) -1
    except ValueError: 
        continue
    break

import disk_images_panel as prettyplot
import disk_images_wmom0 as prettymom
import disk_images_wmom0_2s as prettymom2
with open(filesourcename,'r') as f:
    print('Opened file')
    raw = f.read().splitlines()
    print('Found ',len(raw),' lines')
for i in range(int(len(raw)/40)):
    tmp0=raw[40*i]
    tmp1=raw[40*i+1]
    tmp2=float(raw[40*i+2])
    tmp3=float(raw[40*i+3])
    tmp4=float(raw[40*i+4])
    tmp5=float(raw[40*i+5])
    tmp6=[float(x) for x in raw[40*i+6].split(',')]
    tmp7=float(raw[40*i+7])
    tmp8=float(raw[40*i+8])
    tmp9=r'{}'.format(raw[40*i+9])
    tmp10=float(raw[40*i+10])
    tmp11=raw[40*i+11]
    tmp12=float(raw[40*i+12])
    tmp13=raw[40*i+13]
    tmp14=raw[40*i+14]
    tmp15=raw[40*i+15]
    tmp16=raw[40*i+16]
    tmp17=raw[40*i+17]
    tmp18=raw[40*i+18]
    tmp19=raw[40*i+19]
    tmp100=raw[40*i+20]
    tmp101=float(raw[40*i+21])
    tmp102=float(raw[40*i+22])
    tmp103=float(raw[40*i+23])
    tmp104=raw[40*i+24]
    tmp105=raw[40*i+25]
    tmp106=float(raw[40*i+26])
    tmp107=float(raw[40*i+27])
    tmp108=float(raw[40*i+28])
    tmp109=raw[40*i+29]
    tmp110=raw[40*i+30]
    tmp111=r'{}'.format(raw[40*i+31])
    tmp112=raw[40*i+32]
    tmp113=raw[40*i+33]
    tmp114=float(raw[40*i+34])
    tmp115=float(raw[40*i+35])
    tmp116=float(raw[40*i+36])
    tmp117=float(raw[40*i+37])
    tmp118=float(raw[40*i+38])
    tmp119=float(raw[40*i+39])
    os.system("rm -vf " + tmp1 + " " + tmp110)
    print('Answer: ',answer)
    if i == answer:
        tmp11 = 'n'
    prettyplot.main(tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,'pdf')
    if tmp114 == '':
        prettymom.main(tmp0,tmp110,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp111,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15,tmp100,tmp101,tmp102,tmp103,tmp104,tmp105,tmp106,tmp107,tmp108,tmp109,tmp112,tmp17,tmp18,tmp113,'pdf')
    if tmp114 != '':
        prettymom2.main(tmp0,tmp110,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp111,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15,tmp100,tmp101,tmp102,tmp103,tmp104,tmp105,tmp106,tmp107,tmp108,tmp109,tmp112,tmp17,tmp18,tmp113,tmp114,tmp115,tmp116,tmp117,tmp118,tmp119,'pdf')
