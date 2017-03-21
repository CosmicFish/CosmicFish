#----------------------------------------------------------------------------------------
#
# This file is part of CosmicFish.
#
# Copyright (C) 2015-2017 by the CosmicFish authors
#
# The CosmicFish code is free software;
# You can use it, redistribute it, and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# The full text of the license can be found in the file LICENSE at
# the top level of the CosmicFish distribution.
#
#----------------------------------------------------------------------------------------

"""

Simple Python code to plot a set of SN mocks.

The CosmicFish SN Fisher matrix calculator outputs SN Mock data sets files and this Python
script plots them.

The ouput will be a pdf file with the plot.

Invoking the help option ``compare_SN_mocks.py -h`` will result in::

    usage: compare_SN_mocks.py [-h] [-o OUTROOT] [-v] [-q] files [files ...]
    
    SN mock plotter
    
    positional arguments:
      files                 a list of files with SN mocks
    
    optional arguments:
      -h, --help            show this help message and exit
      -o OUTROOT, --outroot OUTROOT
                            path and name of the output file
      -v, --version         show program's version number and exit
      -q, --quiet           decides wether something gets printed to screen or not


Developed by Marco Raveri (mraveri@uchicago.edu) for the CosmicFish code.

"""

# ***************************************************************************************

__version__ = '1.1' #: version of the application

# ***************************************************************************************

""" Hard coded options """

x_spacing     = 1.0           #: x spacing between plots on the grid. In cm.
y_spacing     = 1.0           #: y spacing between plots on the grid. In cm.
x_size        = 10.0           #: x dimension of the single subplot. In cm.
y_size        = 5.0           #: y dimension of the single subplot. In cm.
main_fontsize = 10.0          #: fontsize for the plots
errorevery    = 5             #: plot error bars on a sub sample of points
dpi           = 300           #: dpi for non-vector figure export

# ***************************************************************************************

# import first dependencies:
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines    as mlines
import numpy as np
import argparse
import math
import sys
import os

matplotlib.rcParams['savefig.dpi'] = dpi

from matplotlib.patches import Rectangle

# get the path of the application and the CosmicFish library:
here = os.path.dirname(os.path.abspath(__file__))
cosmicfish_pylib_path = here+'/..'
sys.path.insert(0, os.path.normpath(cosmicfish_pylib_path))

# import the CosmicFish pylib
import cosmicfish_pylib.utilities as fu
import cosmicfish_pylib.colors    as fc

# ***************************************************************************************

# protection against importing:
if __name__ == "__main__":
    
    # parse command line arguments:
    parser = argparse.ArgumentParser(description='SN mock plotter')
    # parse file names:
    parser.add_argument('files', metavar='files', type=str, nargs='+',
                         help='a list of files with SN mocks')
    # parse the output root:
    parser.add_argument('-o','--outroot', dest='outroot', type=str,
                         help='path and name of the output file')
    # version:
    parser.add_argument('-v','--version', action='version', version='%(prog)s '+__version__)
    # quiet mode:
    parser.add_argument('-q','--quiet', dest='quiet', action='store_true', 
                        help='decides wether something gets printed to screen or not')
    # do the parsing:
    args = parser.parse_args()
    # print the CosmicFish header:
    if not args.quiet:
        fu.CosmicFish_write_header(' SN Mock plotter v'+__version__)
    # process input arguments:
    files          = args.files
    number_SN_mock = len(files)
    if args.outroot is not None:
        outroot    = args.outroot
    else:
        outroot    = os.path.join( os.path.splitext(files[0])[0] )
    # import the data:
    data = []
    for i in files:
        data.append( np.loadtxt(i) )
    # extract the names:
    names = []
    for i in files:
        names.append( os.path.splitext(os.path.basename(i))[0] )
    # set up the plots:
    fig = plt.gcf()
    num_plots_x = 2
    num_plots_y = 4
    # convert into inches:
    x_spacing = x_spacing/2.54
    y_spacing = y_spacing/2.54
    x_size  = (x_size+x_spacing)/2.54*num_plots_x
    y_size  = (y_size+y_spacing)/2.54*num_plots_y
    # setup the plot:
    fig.set_size_inches( x_size, y_size )
    plot_grid = gridspec.GridSpec( num_plots_y, num_plots_x, wspace=x_spacing, hspace=y_spacing )

    # print to screen number of SN:
    if not args.quiet:
        for k in xrange(0, number_SN_mock):
            print 'Total number of SN in '+names[k]+' :', len(data[k][:,1])
        print
    
    # plot SN redshift distribution and luminosity distance:
    ax_1  = plt.subplot(plot_grid[1,0]) # distribution
    ax_2  = plt.subplot(plot_grid[1,1]) # d_L
    # cycle over files:
    for k in xrange(0, number_SN_mock):
        # grab the x values:
        x_values = data[k][:,1]
        y_values = data[k][:,2]
        # do the instogram
        ax_1.hist(x_values, 50, histtype='step', color=fc.nice_colors(k) ) 
        # do the luminosity plot:
        y_err    = np.sqrt( data[k][:,5] )
        # do the plot:
        ax_2.errorbar(x_values, y_values, yerr=y_err, 
                    fmt='o', color=fc.nice_colors(k), 
                    markersize=1, markeredgecolor = 'none', 
                    markerfacecolor=fc.nice_colors(k), errorevery=errorevery, capsize=0)
        
    # finalise the plot:
    ax_1.set_xlabel( '$z$' )
    ax_1.set_ylabel( '$n_{\\rm SN}(z)$' )
    ax_2.set_xlabel( '$z$' )
    ax_2.set_ylabel( '$\mu_{\\rm obs}(z)$' )

    # plot SN mock color:
    ax_1  = plt.subplot(plot_grid[2,0])
    ax_2  = plt.subplot(plot_grid[2,1])
    # cycle over files:
    for k in xrange(0, number_SN_mock):
        # grab the x values:
        x_values = data[k][:,1]
        y_values = data[k][:,3]
        # do the plot:
        ax_1.hist(y_values, 50, histtype='step', color=fc.nice_colors(k) )
        ax_2.plot(x_values, y_values,
                    'o', color=fc.nice_colors(k), 
                    markersize=1, markeredgecolor = 'none', 
                    markerfacecolor=fc.nice_colors(k) )
        
    # finalise the first plot:
    ax_1.set_xlabel( '$C$' )
    ax_1.set_ylabel( '$n_{C}$' )
    ax_2.set_xlabel( '$z$' )
    ax_2.set_ylabel( '$C(z)$' )
    
    # plot SN mock stretch:
    ax_1  = plt.subplot(plot_grid[3,0])
    ax_2  = plt.subplot(plot_grid[3,1])
    # cycle over files:
    for k in xrange(0, number_SN_mock):
        # grab the x values:
        x_values = data[k][:,1]
        y_values = data[k][:,4]
        # do the plot:
        ax_1.hist(y_values, 50, histtype='step', color=fc.nice_colors(k) )
        ax_2.plot(x_values, y_values,
                    'o', color=fc.nice_colors(k), 
                    markersize=1, markeredgecolor = 'none', 
                    markerfacecolor=fc.nice_colors(k) )
        
    # finalise the first plot:
    ax_1.set_xlabel( '$X_1$' )
    ax_1.set_ylabel( '$n_{X_1}$' )
    ax_2.set_xlabel( '$z$' )
    ax_2.set_ylabel( '$X_1(z)$' )

    # create legend handlers
    leg_handlers = []
    for k in xrange(0, number_SN_mock):
        leg_handlers.append( mlines.Line2D([], [], color=fc.nice_colors(k)) )
    # create the legend
    legend_anchor =  plot_grid[0,0].get_position(fig)
    legend_anchor.set_points( np.array( [[0.0,legend_anchor.y0], [1.0,legend_anchor.y0+legend_anchor.height]] ) )
    lgd = fig.legend( handles=leg_handlers, labels=names, fancybox=True,
                      bbox_to_anchor=legend_anchor, borderaxespad=0.0, fontsize=main_fontsize, ncol=2, loc='center')
    # global title of the plot
    tit = plt.suptitle('SN Mocks', fontsize=1.5*main_fontsize)
    # apply tight layout:
    plot_grid.tight_layout( fig )
    # save the figure and close
    plt.savefig(outroot, bbox_extra_artists=(lgd,tit))
    plt.clf()
    # print some final feedback:
    if not args.quiet:
        print 'Done. Saved results in: ', outroot
    # exit without error:
    exit(0)
