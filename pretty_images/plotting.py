
'''
Pretty plots files from ALMA data sets
Original Author John Tobin
Adapted by Nickalas Reynolds
April 2017
'''

# import libraries
from __future__ import print_function
import os
import glob
import aplpy
import sys
import ds9cmap

# make example
print("Making example file.")
with open("example_plot_sources.txt",'w') as f:
    f.write("image name,outfilename,ra,dec,minpixval,maxpixval,size,scalebar,distance,name,pa,showoutflow,sigma,showcontours,showsources,imagestretch,colororgray,colormap,plotlabel,textcolor,moment outfilename , contourplotlabel,contour color scale, contour textcolor\n")
    f.write('Will plot a pretty image and then plot the moment map.\n')
    f.write("# capable of handling multiple sources, just begin the next source on the next line with no line break inbetween.\n")
    f.write("#^^^^^^^^ exclude these first few lines with # ^^^^^^^^\n")
    f.write("L1448IRS3B_cont_superuniform.fits\n") 
    f.write("L1448IRS3B_cont_superuniform_C17O_all.png\n")
    f.write("51.4016208333\n")
    f.write("30.7550666667\n")
    f.write("0.00014571\n")
    f.write("0.0361\n") 
    f.write("15.0\n")
    f.write("0.1\n")
    f.write("230.0\n")
    f.write("L1448 IRS3B Cont 870micron\n")
    f.write("25.0\n")
    f.write("n\n")
    f.write("00.00014571\n")
    f.write("n\n")
    f.write("n\n")
    f.write("linear\n")
    f.write("color\n")
    f.write("jet\n")
    f.write("\n")
    f.write("white\n")
    f.write("L1448IRS3B_C17O_image_taper1500k_r.fits\n")
    f.write("3\n")
    f.write("2\n")
    f.write("0.00231\n")
    f.write("y\n")
    f.write("L1448IRS3B_C17O_image_taper1500k_b.fits\n")
    f.write("3\n")
    f.write("2\n")
    f.write("0.00675\n")
    f.write("y\n")
    f.write("L1448IRS3B_cont_superuniform_C17O_mom0_all.png\n")
    f.write("L1448 IRS3B C17O\n")
    f.write("gray\n")
    f.write("black")
    f.close()

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
    raw = f.read().splitlines()
for i in range(int(len(raw)/40)):
    tmp0=raw[40*i]
    tmp1=raw[40*i+1]
    tmp2=float(raw[40*i+2])
    tmp3=float(raw[40*i+3])
    tmp4=float(raw[40*i+4])
    tmp5=float(raw[40*i+5])
    tmp6=float(raw[40*i+6])
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
    if i == answer:
        tmp11 = 'n'
    prettyplot.main(tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15,tmp16,tmp17,tmp18,tmp19,'pdf')
    if tmp114 == '':
        prettymom.main(tmp0,tmp110,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp111,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15,tmp100,tmp101,tmp102,tmp103,tmp104,tmp105,tmp106,tmp107,tmp108,tmp109,tmp112,tmp17,tmp18,tmp113,'pdf')
    if tmp114 != '':
        prettymom2.main(tmp0,tmp110,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp111,tmp10,tmp11,tmp12,tmp13,tmp14,tmp15,tmp100,tmp101,tmp102,tmp103,tmp104,tmp105,tmp106,tmp107,tmp108,tmp109,tmp112,tmp17,tmp18,tmp113,tmp114,tmp115,tmp116,tmp117,tmp118,tmp119,'pdf')
