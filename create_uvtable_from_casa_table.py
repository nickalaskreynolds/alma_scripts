#! /usr/bin/env python

"""
create_uvtable_from_casa_table: 
                Program used to create a table of visibi-
                lities (u,v,Re,Im,Weight) from an input
                visibility file created by the CASA table
                tools. This program computes the visibility
                value for I (coming from XX,YY dual pol
                visibilities) and creates a simple nice
                uv-table that can be used elsewhere.

VERY IMPORTANT: XX and YY polarizations coming from the dual-pol 
                receivers are calibrated separately (for historical
                reasons), hence they separately provide the total 
                intensity I. This means that each visibility coming
                from XX and YY are like two independent visibility 
                measurements of the total intensity, and you will
                need to average them together to get 'I' instead of 
                adding them together (as presented in TMS eq. 4.58)

"""

import sys
from math import pi
from math import sin
from math import cos

usage = "usage: create_uvtable_from_tb filename_input filename_output"

if len(sys.argv) != 3:
    print(usage)
    sys.exit(1)

in_file  = sys.argv[1]
out_file = sys.argv[2]

# Open the input file (Formated by CASA table tools)
# casapy> tb.toasciifmt(asciifile='output.txt', columns=['UVW','DATA','WEIGHT'])
f_in = open(in_file, "r")

# Open the output file (Format will be: u,v,Re,Im,We)
f_out = open(out_file, "w")

# Read one line at a time and write output file:
i = 0
for line in f_in:
    i += 1
    if i > 2.0:
        
        variable = line.replace('[',' ').replace(']',' ').replace('(',' ').replace(')',' ').replace(',',' ').split()

        u = float(variable[0])
        v = float(variable[1])
        Re_xx = float(variable[3])
        Im_xx = float(variable[4])
        Re_yy = float(variable[5])
        Im_yy = float(variable[6])
        wei_xx = float(variable[7])
        wei_yy = float(variable[8])
        
        # Add XX and YY components
        Re = (Re_xx + Re_yy) / 2.0
        Im = (Im_xx + Im_yy) / 2.0
        We = (wei_xx**2 + wei_yy**2)**0.5 / 2.0
        
        # Print visibilities to file
        f_out.write("%f \t %f \t %f \t %f \t %f" % (u,v,Re,Im,We))
        
print(variable)

# Close files    
f_in.close()
f_out.close()

# Print output message
print("UV-table creation is completed!")
print("Visibilities written in: %s" % (out_file))

