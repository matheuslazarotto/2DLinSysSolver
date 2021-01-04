# Function plot 
#
#   Usage:
#   python3.6 plot.py
#
#   Matheus J. Lazarotto (30/10/2020)

# Import auxiliary file
import sys
sys.path.append('../python')
from aux_par import *
from aux_plt import *

import os
import argparse

parser = argparse.ArgumentParser()
# Arg: minimum x value
parser.add_argument('-x_min', help = 'Minimum x value.',
                    action = 'store', type = float, dest = 'x_min')
# Arg: maximum x value
parser.add_argument('-x_max', help = 'Maximum x value.',
                    action = 'store', type = float, dest = 'x_max')
# Arg: minimum y value
parser.add_argument('-y_min', help = 'Minimum y value.',
                    action = 'store', type = float, dest = 'y_min')
# Arg: maximum x value
parser.add_argument('-y_max', help = 'Maximum y value.',
                    action = 'store', type = float, dest = 'y_max')
args = parser.parse_args()

# Data parse function
def parse_data(flname):
    with open(flname, 'r') as inpt:
        dat = []

        for line in inpt:
            if line[:3] == '# [' in line:
                break
            elif 'a11 =' in line:
                a11 = float(line.split()[-1])
            elif 'a12 =' in line:
                a12 = float(line.split()[-1])
            elif 'a21 =' in line:
                a21 = float(line.split()[-1])
            elif 'a22 =' in line:
                a22 = float(line.split()[-1])

        dat_x, dat_y = [], []
        dat_dxdt, dat_dydt = [], []
            
        while (line[:3] == '# [' and line != ''):
            
            tmp_dxdt, tmp_dydt = [], []
            tmp_x, tmp_y = [], []
            line = inpt.readline()

            while (line[:3] != '# [' and line != ''):
                tmp_x.append(float(line.split()[1]))
                tmp_y.append(float(line.split()[2]))
                line = inpt.readline()

            dat_x.append(tmp_x)
            dat_y.append(tmp_y)
            
        dat.append(dat_x)
        dat.append(dat_y)

        return dat

####################################
#-------------- Main --------------#
####################################

# Read files from working directory.
wd = os.getcwd()
for flname in os.listdir(wd):
    if flname.endswith('.dat') and 'out_points' in flname:
        dat = parse_data(flname)
        plot_streamline(dat, args.x_min, args.x_max, 
                             args.y_min, args.y_max)

