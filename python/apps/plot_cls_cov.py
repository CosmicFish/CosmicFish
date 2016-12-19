#----------------------------------------------------------------------------------------
#
# This file is part of CosmicFish.
#
# Copyright (C) 2015-2016 by the CosmicFish authors
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

Simple Python code to plot a set of cls covariance.

The CosmicFish Cls Fisher matrix calculator outputs Cls covariance files and this Python
script plots them.

The ouput will be a pdf file with the plot.

Invoking the help option ``plot_cls_cov.py -h`` will result in::

    usage: plot_cls_cov.py [-h] [-o OUTROOT] [-l LABELS [LABELS ...]] [-v] [-q]
                           files [files ...]

    Cls covariance plotter

    positional arguments:
      files                 a list of files with Cls covariances

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTROOT, --outroot OUTROOT
                            path and name of the output file
      -l LABELS [LABELS ...], --label LABELS [LABELS ...]
                            set the labels for the legend. If empty will use the
                            file names.
      -v, --version         show program's version number and exit
      -q, --quiet           decides wether something gets printed to screen or not

Developed by Marco Raveri (mraveri@uchicago.edu) for the CosmicFish code.

"""

# ***************************************************************************************

__version__ = '1.1' #: version of the application

# ***************************************************************************************

""" Hard coded options """

do_abs        = True          #: wether to plot the absolute value of the Cls or not
cutoff        = 10.0**(-16)   #: lower bound where the values of the Cls is flattened
x_spacing     = 1.0           #: x spacing between plots on the grid. In cm.
y_spacing     = 1.0           #: y spacing between plots on the grid. In cm.
x_size        = 4.0           #: x dimension of the single subplot. In cm.
y_size        = 4.0           #: y dimension of the single subplot. In cm.
adjust_l      = True          #: wether to plot l(l+1)Cls/2pi or Cls.
do_lin        = False         #: wether to do linear or log plots.
main_fontsize = 10.0          #: fontsize for the plots
dpi           = 300               #: dpi for non-vector figure export

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
    parser = argparse.ArgumentParser(description='Cls covariance plotter')
    # parse file names:
    parser.add_argument('files', metavar='files', type=str, nargs='+',
                         help='a list of files with Cls covariances')
    # parse the output root:
    parser.add_argument('-o','--outroot', dest='outroot', type=str,
                         help='path and name of the output file')
    # parse the plot legend labels:
    parser.add_argument('-l','--label', dest='labels', type=str, nargs='+',
                         help='set the labels for the legend. If empty will use the file names.')
    # version:
    parser.add_argument('-v','--version', action='version', version='%(prog)s '+__version__)
    # quiet mode:
    parser.add_argument('-q','--quiet', dest='quiet', action='store_true', 
                        help='decides wether something gets printed to screen or not')
    # do the parsing:
    args = parser.parse_args()
    # print the CosmicFish header:
    if not args.quiet:
        fu.CosmicFish_write_header(' Cls covariance plotter v'+__version__)
    # process input arguments:
    files          = args.files
    number_cls_cov = len(files)
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
    if args.labels is None:
        # extract from file names:
        for i in files:
            names.append( os.path.splitext(os.path.basename(i))[0] )
    else:
        names = args.labels
    # check wether all the cls covariances have the same size:
    for i in xrange(0, number_cls_cov):
        if ( data[i].shape[1] != data[0].shape[1] ):
            raise ValueError( 'The cls covariance '+names[i]+' does not have the same size as '+names[0] )
    # number of Cls:
    num     = int(math.sqrt(data[0].shape[1]-1))
    # set up the plots:
    fig = plt.gcf()
    # convert into inches:
    x_spacing = x_spacing/2.54
    y_spacing = y_spacing/2.54
    x_size  = (x_size+x_spacing)/2.54*num
    y_size  = (y_size+y_spacing)/2.54*num
    # protection against very small figures:
    x_size = max( x_size, 21.0/2.54 )
    y_size = max( y_size, 21.0/2.54 )
    # setup the plot:
    fig.set_size_inches( x_size, y_size )
    plot_grid = gridspec.GridSpec( num, num, wspace=x_spacing, hspace=y_spacing )
    # loop over the Cls grid:
    for ind in xrange(1, num+1):
        for ind2 in xrange(1, ind+1):
    
            ax  = plt.subplot(plot_grid[ind-1,ind2-1])
            col  = ind + num*(ind2-1)
    
            plotrange_x = []
            plotrange_y = []
            for k in xrange(0, number_cls_cov):
    
                # grab the x values:
                x_values = data[k][:,0]
                # convert the bare Cls to l(l+1)/2piCls
                if adjust_l:
                    factor   = (x_values*(x_values+1))/(2.0*math.pi)
                else:
                    factor   = 1.0
                # get the Cls values:
                if do_abs:
                    y_values = np.abs( factor*data[k][:,col] )
                else:
                    y_values = factor*data[k][:,col]
                # protect against zeroes:
                np.place(y_values, y_values==0.0, [cutoff])
                # get the plot bounds:
                try: lower_y = 0.9*np.amin(y_values[np.abs(y_values)>cutoff])
                except: lower_y = cutoff
                try: upper_y = 1.1*np.amax(y_values[np.abs(y_values)>cutoff])
                except: upper_y = cutoff
                plotrange_x.append( [ np.amin(x_values), np.amax(x_values) ] )
                plotrange_y.append( [ lower_y,upper_y ] )
                # do the plot:
                ax.plot( x_values, y_values, color=fc.nice_colors(k) )
            # now set the appearence of the plot:
            plotrange_x = np.array( plotrange_x )
            plotrange_y = np.array( plotrange_y )
            # log scale:
            if not do_lin:
                ax.set_xscale('log')
                ax.set_yscale('log')
            # ranges:
            x_min = fu.v_nice_number( np.amin(plotrange_x[:,0]),1)
            x_max = fu.v_nice_number( np.amax(plotrange_x[:,1]),1)
            if ( x_min != x_max ):
                ax.set_xlim( [ x_min,x_max ] )
            y_min = fu.v_nice_number( np.amin(plotrange_y[:,0]),2)
            y_max = fu.v_nice_number( np.amax(plotrange_y[:,1]),0)
            if ( y_min != y_max ):
                ax.set_ylim( [ y_min,y_max ] )
            # labels
            if ( ind != num ):
                ax.set_xticklabels( [] )
            if ( ind == num ):
                ax.set_xlabel('$\ell$')
            if ( ind2 == 1):
                if adjust_l:
                    ax.set_ylabel( '$\ell(\ell+1) \, C_{\ell} \,\, / \,\, 2\pi $' )
                else:
                    ax.set_ylabel( '$C_{\ell}$' )
            
    # create legend handlers
    leg_handlers = []
    for k in xrange(0, len(names)):
        leg_handlers.append( mlines.Line2D([], [], color=fc.nice_colors(k)) )
    # create the legend
    legend_anchor =  plot_grid[0,num-1].get_position(fig)
    lgd = fig.legend( handles=leg_handlers, labels=names, fancybox=True,
                      bbox_to_anchor=legend_anchor, borderaxespad=0.0, fontsize=main_fontsize)
    
    # apply tight layout:
    plot_grid.tight_layout( fig )
    # global title of the plot
    tit = plt.suptitle('Cls matrix', 
                        fontsize=1.5*main_fontsize,
                        x=legend_anchor.corners()[3][0], 
                        y=legend_anchor.corners()[3][1],
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        )
    # save the figure and close
    plt.savefig(outroot, bbox_extra_artists=(lgd,tit))
    plt.clf()
    # print some final feedback:
    if not args.quiet:
        print 'Done. Saved results in: ', outroot
    # exit without error:
    exit(0)
