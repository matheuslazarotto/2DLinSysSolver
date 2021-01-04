# Auxiliary plot routines
#
# Made by: Matheus J. Lazarotto (30/09/2020)

# Import auxiliary file
import sys
sys.path.append('../python')
from aux_par import *

from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FormatStrFormatter

# Plot routine 
def plot_streamline(dat, Xmin, Xmax, Ymin, Ymax):
    # Edit text font
    font = {'family' : 'serif',
            'weight' : 'normal',
            'size'   : 55}
    plt.rc('font', **font)
    
    # plt.rc('text', usetex = True)  # LaTeX text
    # plt.rc('font', family = 'serif', size = 80)
    fig, ax = plt.subplots(figsize = (30, 25))

    # Set X axes #
    ax.set_xlim([Xmin, Xmax])
    ax.set_xlabel("x", rotation = 0, labelpad = 15, fontsize = 70)

    # Set Y axes #
    ax.set_ylim([Ymin, Ymax])
    ax.set_ylabel("y", rotation = 0, labelpad = 25, fontsize = 70)

    # Set global features #
    ax.tick_params(which = 'major', direction = 'in', length = 18, 
                   width = 2.2, bottom = True, top = True, left = True, 
                   right = True, pad = 20)
    ax.tick_params(which = 'minor', direction = 'in', length = 10, 
                   width = 1.8, bottom = True, top = True, left = True, 
                   right = True)
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(4.0)

    # Loop over data dictionary
    color = 'black'
    color_arrow = 'blue'
    for i in range(len(dat[0])):
        plt.plot(dat[0][i], dat[1][i], linestyle = '-', linewidth = 4, 
                                       color = color, zorder = 1)
        for pt in range(len(dat[0][i])):
           if ((pt % 100 == 0) and (pt+1 != len(dat[0][i])) and
              dat[0][i][pt] != 0 and dat[0][i][pt] != 0 and 
              (dat[0][i][pt+1] - dat[0][i][pt]) != 0 and
              (dat[0][i][pt+1] - dat[0][i][pt]) != 0):
               
               plt.arrow(dat[0][i][pt], dat[1][i][pt], 
                         (dat[0][i][pt+1] - dat[0][i][pt]), 
                         (dat[1][i][pt+1] - dat[1][i][pt]), 
                         shape = 'full', lw = 0, length_includes_head = True, 
                         head_width=0.03, color = color_arrow, zorder = 2)
    
    plt.grid()
    # plt.legend(loc = 'best', ncol = 2, prop = {'size': 40})

    # Save and close #
    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(1.0)
    ax.patch.set_facecolor('white')
    ax.patch.set_alpha(1.0)
    plt.savefig('functions.png', bbox_inches='tight', 
                 facecolor=fig.get_facecolor(), edgecolor='none')
    plt.clf()
    plt.close()
