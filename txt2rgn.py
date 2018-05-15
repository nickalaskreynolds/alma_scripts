#!/usr/bin/env python
from argparse import ArgumentParser as ap

parser = ap()
parser.add_argument('-i', '--input',default='cores.txt')
parser.add_argument('-o', '--output',default='cores.reg')
parser.add_argument('-op', '--options',default=0,type=int)
parser.add_argument('-c', '--colour',default='green')
args = parser.parse_args()

option  = args.options
infile  = args.input
outfile = args.output

ds9="box($0,$1,10\",10\") # text={$2}\n"

if option == 0:
    with open(infile,'r') as f:
        allines = [x.strip('\n').split() for x in f.readlines() if (x.strip('\n') != '') and (x[0] != '#')]


    with open(outfile,'w') as f:
        f.write('global color={} dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1' \
            ' highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n'.format(args.colour))
        for x in allines:
            nstr = ds9.replace('$0',x[1]).replace('$1',x[2]).replace('$2',x[0])
            print(nstr)
            f.write(nstr)

#
if option == 1: # for cox 2018
    with open(infile,'r') as f:
        allines = [x.strip('\n') for x in f.readlines() if x.strip('\n') != ''] 
    stripped = []
    i = 0
    while (i < (len(allines)-1)):
        cur = allines[i]
        if len(cur.split(':')) > 1:
            stripped.append([allines[i-1],cur,allines[i+1]])
            i+=2
        else:
            i+=1
        
    with open(outfile,'w') as f:
        f.write('global color={} dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1' \
            ' highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n'.format(args.colour))
        for x in stripped:
            f.write('{}  {} {}\n'.format(*x))


if option == 2: # for Young et al 2015
    with open(infile,'r') as f:
        allines = [x.strip('\n') for x in f.readlines() if x.strip('\n') != ''] 

    with open(outfile,'w') as f:
        f.write('global color={} dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1' \
            ' highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n'.format(args.colour))
        nstrip = []
        for y in allines:
            if y[0] == '#':
                f.write(y+'\n')
            else:
                x = y.split(' ')
                nstr = ds9.replace('$0',x[1]).replace('$1',x[2]).replace('$2',x[0])
                f.write(nstr)

if option == 3: # for Sadavoy et al 2014
    with open(infile,'r') as f:
        allines = [x.strip('\n') for x in f.readlines() if x.strip('\n') != ''] 
    stripped = []
    i = 0
    while (i < (len(allines)-1)):
        cur = allines[i]
        if len(cur.split(':')) > 1:
            a = allines[i-2].strip(' ')
            b = allines[i-1].strip(' ')
            c = cur.split('   ')[0]
            d = cur.split('   ')[1]
            stripped.append(["{}_{}".format(a,b),c,d])
            i+=2
        else:
            i+=1
        
    with open(outfile,'w') as f:
        f.write('global color={} dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1' \
            ' highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n'.format(args.colour))
        for x in stripped:
            f.write('{}  {} {}\n'.format(*x))

if option == 4: # for Vandam multiple
    with open(infile,'r') as f:
        allines = [x.strip('\n') for x in f.readlines() if x.strip('\n') != ''] 
    stripped = []
    for x in allines:
        y = x.split('  ')
        stripped.append([y[1],":".join(y[2].split(' ')),":".join(y[4].split(' '))])
        
    with open(outfile,'w') as f:
        f.write('global color={} dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1' \
            ' highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n'.format(args.colour))
        for x in stripped:
            f.write('{}  {} {}\n'.format(*x))

if option == 5: # for possible
    with open(infile,'r') as f:
        allines = [x.strip('\n').split() for x in f.readlines() if (x.strip('\n') != '') and (x[0] != '#')]


    with open(outfile,'w') as f:
        f.write('global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1' \
            ' highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')
        for x in allines:
            nstr = ds9.replace('$0',x[1]).replace('$1',x[2]).replace('$2',x[0])
            print(nstr)
            f.write(nstr)
        f.write('box(3:32:51.919,+31:01:50.099,648.629",616.319",359.04385) # text={Reynolds18}\n')
        f.write('box(3:28:38.343,+31:03:24.742,659.451",593.740",359.04385) # text={Reynolds18}\n')
        f.write('box(3:28:36.058,+30:20:48.138,501.642",424.634",359.04385) # text={Reynolds18}\n')
        f.write('box(3:26:36.119,+30:16:03.969,453.426",343.537",359.04385) # text={Reynolds18}\n')
        f.write('box(3:28:49.964,+30:44:59.815,386.266",329.523",359.04385) # text={Reynolds18}\n')
        f.write('box(3:30:27.003,+30:27:03.403,571.192",673.606",359.04385) # text={Reynolds18}\n')

if option == 6: # generate perseus region boxes
    with open(outfile,'w') as f:
        f.write('global color=red dashlist=8 3 width=1 '+\
            'font="helvetica 10 normal roman" select=1 '+\
            'highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        f.write('fk5\n')
        f.write('box(3:28:52.753,+31:19:08.248,2285.198",2419.549",359.04385) # text={NGC1333}\n')
        f.write('box(3:25:23.991,+30:44:44.486,1719.280",977.603",359.04385)  # text={L1448}\n')
        f.write('box(3:32:54.951,+31:05:53.169,2148.804",2363.770",359.04385) # text={B1}\n')
        f.write('box(3:35:56.349,+31:09:07.653,1576.778",1715.353",359.04385) # text={B1E}\n')
        f.write('box(3:39:23.335,+31:20:26.009,3243.804",1881.050",359.04385) # text={L1468}\n')
        f.write('box(3:43:48.117,+32:04:03.412,3667.516",2918.979",359.04385) # text={IC348}\n')
        f.write('box(3:47:25.352,+32:50:28.866,1617.085",1811.439",359.04385) # text={B5}\n')
        f.write('box(3:27:57.091,+30:07:20.606,2109.465",2325.722",359.04385) # text={L1455}\n')
        f.write('box(3:24:42.689,+30:16:17.297,2089.241",1360.326",359.04385) # text={L1451}\n')
