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

Simple Python code to compare two cls covariance.

The CosmicFish Cls Fisher matrix calculator outputs Cls covariance files and this Python
script compares two of them.

The ouput will be a pdf file with the plot.

Invoking the help option ``compare_cls_cov.py -h`` will result in::

    usage: compare_cls_cov.py [-h] [-o OUTROOT] [-l LABELS [LABELS ...]] [-v] [-q]
                              files files

    Cls covariance comparison plotter

    positional arguments:
      files                 two files with the Cls covariances to compare

    optional arguments:
      -h, --help            show this help message and exit
      -o OUTROOT, --outroot OUTROOT
                            path and name of the output file
      -l LABELS [LABELS ...], --label LABELS [LABELS ...]
                            set the labels for the legend. If empty will use the
                            file names.
      -f, --filling         whether to fill the curves. Default does not fill.
      -v, --version         show program's version number and exit
      -q, --quiet           decides wether something gets printed to screen or not

Developed by Marco Raveri (mraveri@sissa.it) for the CosmicFish code.

"""

# ***************************************************************************************

__version__ = '1.0' #: version of the application

# ***************************************************************************************

""" Hard coded options """

do_abs        = True              #: wether to plot the absolute value of the Cls or not
cutoff        = 10.0**(-16)       #: lower bound where the values of the Cls is flattened
x_spacing     = 1.0               #: x spacing between plots on the grid. In cm.
y_spacing     = 1.0               #: y spacing between plots on the grid. In cm.
x_size        = 4.0               #: x dimension of the single subplot. In cm.
y_size        = 4.0               #: y dimension of the single subplot. In cm.
do_lin        = False             #: wether to do linear or log plots.
main_fontsize = 10.0              #: fontsize for the plots
color         = (42.0/255.0, 46.0/255.0, 139.0/255.0) #: color to use for the plot

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
from scipy import interpolate

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
    parser = argparse.ArgumentParser(description='Cls covariance comparison plotter')
    # parse file names:
    parser.add_argument('files', metavar='files', type=str, nargs=2,
                         help='two files with the Cls covariances to compare')
    # parse the output root:
    parser.add_argument('-o','--outroot', dest='outroot', type=str,
                         help='path and name of the output file')
    # parse the plot legend labels:
    parser.add_argument('-l','--label', dest='labels', type=str, nargs='+',
                         help='set the labels for the legend. If empty will use the file names.')
    # whether we want filling:
    parser.add_argument('-f','--filling', dest='filling', action='store_true', 
                         help='whether to fill the curves. Default does not fill.')
    # version:
    parser.add_argument('-v','--version', action='version', version='%(prog)s '+__version__)
    # quiet mode:
    parser.add_argument('-q','--quiet', dest='quiet', action='store_true', 
                        help='decides wether something gets printed to screen or not')
    # do the parsing:
    args = parser.parse_args()
    # print the CosmicFish header:
    if not args.quiet:
        fu.CosmicFish_write_header(' Cls covariance comparison plotter')
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
    
    for ind in xrange(1, num+1):
        for ind2 in xrange(1, ind+1):
    
            ax  = plt.subplot(plot_grid[ind-1,ind2-1])
            col  = ind + num*(ind2-1)
    
            plotrange_x = []
            plotrange_y = []
            
            # perform interpolation of the Cls:
            Cls_1  = interpolate.InterpolatedUnivariateSpline(data[0][:,0], data[0][:,col])
            Cls_2  = interpolate.InterpolatedUnivariateSpline(data[1][:,0], data[1][:,col])
            
            # get the x minimum and maximum values for the plots:
            x_min = np.amax( np.array([np.amin(data[0][:,0]), np.amin(data[1][:,0])]) )
            x_max = np.amin( np.array([np.amax(data[0][:,0]), np.amax(data[1][:,0])]) )
            
            # get the x values:
            x_values = np.linspace( x_min, x_max, x_max-x_min )
            
            # get the comparison values:
            yvalues_1 = Cls_1(x_values)
            yvalues_2 = Cls_2(x_values)
            
            # protect against zeroes:
            yvalues_1_temp = np.abs( yvalues_1 )
            yvalues_2_temp = np.abs( yvalues_2 )
            try:
                min2val_1      = np.amin(yvalues_1_temp[np.nonzero(yvalues_1_temp)])
                min2val_2      = np.amin(yvalues_2_temp[np.nonzero(yvalues_2_temp)])
            except:
                min2val_1      = cutoff
                min2val_2      = cutoff
            np.place(yvalues_1, yvalues_1 == 0.0, min2val_1)
            np.place(yvalues_2, yvalues_2 == 0.0, min2val_2)
            # computation of the percentual comparison:
            yvalues_comp_1 = (yvalues_1 - yvalues_2)/np.abs(yvalues_1)*100.0
            yvalues_comp_2 = (yvalues_1 - yvalues_2)/np.abs(yvalues_2)*100.0
            # protection against values too small:
            np.place(yvalues_comp_1, abs(yvalues_comp_1)<cutoff, [cutoff])
            np.place(yvalues_comp_2, abs(yvalues_comp_2)<cutoff, [cutoff])
            # get the Cls values:
            if do_abs:
                yvalues_comp_1 = np.abs( yvalues_comp_1 )
                yvalues_comp_2 = np.abs( yvalues_comp_2 )
            
            # get the plot bounds for comp 1:
            try: lower_y = 0.9*np.amin(yvalues_comp_1[np.abs(yvalues_comp_1)>cutoff])
            except: lower_y = cutoff
            try: upper_y = 1.1*np.amax(yvalues_comp_1[np.abs(yvalues_comp_1)>cutoff])
            except: upper_y = cutoff
            plotrange_x.append( [ x_min  , x_max   ] )
            plotrange_y.append( [ lower_y, upper_y ] )
            # get the plot bounds for comp 2:
            try: lower_y = 0.9*np.amin(yvalues_comp_2[np.abs(yvalues_comp_2)>cutoff])
            except: lower_y = cutoff
            try: upper_y = 1.1*np.amax(yvalues_comp_2[np.abs(yvalues_comp_2)>cutoff])
            except: upper_y = cutoff
            plotrange_x.append( [ x_min  , x_max   ] )
            plotrange_y.append( [ lower_y, upper_y ] )
                
            # do the plot:
            ax.plot( x_values, yvalues_comp_1, color=color )
            ax.plot( x_values, yvalues_comp_2, color=color )
            if args.filling:
                ax.fill_between(x_values, yvalues_comp_1, yvalues_comp_2, facecolor=color, interpolate=True)
            
            # plot cosmic variance:
            cosmic_variance = np.array( [ math.sqrt(2.0/(2.0*l + 1.0))*100.0 for l in x_values ] )
            ax.plot( x_values, cosmic_variance, color='k' )
            # get the plot bounds for cosmic variance:
            try: lower_y = 0.9*np.amin(cosmic_variance[np.abs(cosmic_variance)>cutoff])
            except: lower_y = cutoff
            try: upper_y = 1.1*np.amax(cosmic_variance[np.abs(cosmic_variance)>cutoff])
            except: upper_y = cutoff
            plotrange_y.append( [ lower_y,upper_y ] )
            
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
                ax.set_ylabel( '$ \Delta C_{\ell} \,\, / \,\, C_{\ell} \,\, \%$' )
            
    # create legend handlers
    leg_handlers = []
    # add the comparison
    names = [ names[0]+'\n VS \n'+names[1] ]
    leg_handlers.append( mlines.Line2D([], [], color=color ) )
    # add cosmic variance to the legend:
    names.append('Cosmic Variance')
    leg_handlers.append( mlines.Line2D([], [], color='k' ) )
    # create the legend
    legend_anchor =  plot_grid[0,num-1].get_position(fig)
    lgd = fig.legend( handles=leg_handlers, labels=names, fancybox=True,
                      bbox_to_anchor=legend_anchor, borderaxespad=0.0, fontsize=main_fontsize)
    # apply tight layout:
    plot_grid.tight_layout( fig )
    # global title of the plot
    tit = plt.suptitle('Cls comparison matrix', 
                        fontsize=1.5*main_fontsize,
                        x=legend_anchor.corners()[3][0], 
                        y=legend_anchor.corners()[3][1],
                        horizontalalignment='right',
                        verticalalignment='bottom',
                        )
    # save the figure and close
    plt.savefig(outroot+'_compCovCls.pdf', bbox_extra_artists=(lgd,tit))
    plt.clf()
    # print some final feedback:
    if not args.quiet:
        print 'Done. Saved results in: ', outroot+'_compCovCls.pdf'
    # exit without error:
    exit(0)
