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

Simple Python code to compare a set of Fisher matrices.

If fed with more than two Fisher matrices the code will compare all of them two by two.

Invoking the help option ``compare_fisher_matrices.py -h`` will result in::

    usage: compare_fisher_matrices.py [-h] [-o OUTROOT] [-p PARAMS [PARAMS ...]]
                                      [-e] [-f FORMAT] [-d] [-v] [-q]
                                      files [files ...]
    
    Fisher matrix comparison plotter
    
    positional arguments:
      files                 a list of files with Fisher matrices
    
    optional arguments:
      -h, --help            show this help message and exit
      -o OUTROOT, --outroot OUTROOT
                            path and name of the output file. Will be added in
                            front of all files.
      -p PARAMS [PARAMS ...], --params PARAMS [PARAMS ...]
                            names of the parameters to plot. If this option is not
                            present the plot will contain all parameters.
      -e, --eliminate       if parameters are passed from the command line this
                            option avoids marginalization over the others
      -f FORMAT, --format FORMAT
                            format for the output plot.
      -d, --derived         decides wether to look for derived parameters when
                            importing the Fisher matrix
      -v, --version         show program's version number and exit
      -q, --quiet           decides wether something gets printed to screen or not

Developed by Marco Raveri (mraveri@uchicago.edu) for the CosmicFish code.

"""

# ***************************************************************************************

__version__ = '1.0' #: version of the application

# ***************************************************************************************

""" Hard coded options """

x_size        = 21.0               #: x and y dimension for 10 elements of the plot. The plot will be squared. In cm.
main_fontsize = 10.0              #: fontsize for the plots
colormap      = 'Blues'
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
import itertools
from mpl_toolkits.axes_grid1 import make_axes_locatable

# get the path of the application and the CosmicFish library:
here = os.path.dirname(os.path.abspath(__file__))
cosmicfish_pylib_path = here+'/..'
sys.path.insert(0, os.path.normpath(cosmicfish_pylib_path))

matplotlib.rcParams['savefig.dpi'] = dpi

# import the CosmicFish pylib
import cosmicfish_pylib.utilities            as fu
import cosmicfish_pylib.colors               as fc
import cosmicfish_pylib.fisher_matrix        as fm
import cosmicfish_pylib.fisher_derived       as fd
import cosmicfish_pylib.fisher_operations    as fo
import cosmicfish_pylib.fisher_plot_settings as fps
import cosmicfish_pylib.fisher_plot_analysis as fpa
import cosmicfish_pylib.fisher_plot          as fp

# ***************************************************************************************

# protection against importing:
if __name__ == "__main__":
    
    # parse command line arguments:
    parser = argparse.ArgumentParser(description='Fisher matrix comparison plotter')
    # parse file names:
    parser.add_argument('files', metavar='files', type=str, nargs='+',
                         help='a list of files with Fisher matrices')
    # parse the output root:
    parser.add_argument('-o','--outroot', dest='outroot', type=str,
                         help='path and name of the output file. Will be added in front of all files.')
    # parse parameters if wanted:
    parser.add_argument('-p','--params', dest='params', type=str, nargs='+',
                         help='names of the parameters to plot. If this option is not present the plot will contain all parameters.')
    # parse an optional argument to avoid marginalization:
    parser.add_argument('-e','--eliminate', dest='eliminate', action='store_true', 
                        help='if parameters are passed from the command line this option avoids marginalization over the others')
    # parse format argument:
    parser.add_argument('-f','--format', dest='format', type=str,
                         help='format for the output plot.')
    # parse a sum argument to sum all Fisher matrices:
    parser.add_argument('-d','--derived', dest='derived', action='store_true', 
                        help='decides wether to look for derived parameters when importing the Fisher matrix')
    # version:
    parser.add_argument('-v','--version', action='version', version='%(prog)s '+__version__)
    # quiet mode:
    parser.add_argument('-q','--quiet', dest='quiet', action='store_true', 
                        help='decides wether something gets printed to screen or not')
    
    # do the parsing:
    args = parser.parse_args()
    
    # print the CosmicFish header:
    if not args.quiet:
        fu.CosmicFish_write_header(' Fisher matrix comparison plotter v'+__version__)
        
    # process input arguments:
    files          = args.files
    number_fish    = len(files)
    if args.outroot is not None:
        outroot    = args.outroot
    else:
        outroot    = os.path.dirname(files[0])+'/'

    params = args.params
    
    output_format = args.format
    if output_format is None:
        output_format = 'pdf'
    
    # import the Fisher matrices:
    if args.derived:
        fishers = fpa.CosmicFish_FisherAnalysis(fisher_path=files, with_derived=True)
    else:
        fishers = fpa.CosmicFish_FisherAnalysis(fisher_path=files, with_derived=False)
    
    # eliminate parameters that are not wanted:
    if args.eliminate:
        if params is None:
            raise ValueError('Avoiding marginalization works only if parameters are specified')
        fishers = fishers.reshuffle( params=params )
    
    # some feedback:
    if not args.quiet:
        print
        
    # get the couples of Fisher matrices:
    fisher_names = list(itertools.combinations( fishers.get_fisher_name_list(), 2))
    
    # print some warning feedback it it's the case:
    if not args.quiet:
        if len(fisher_names) > 10:
            print ' WARNING: comparison will result in '+str(len(fisher_names))+' plots.'
            print ' That is a lot and will take some time...'
            print
    
    # cycle over the couples:
    for fisher_couples in fisher_names:
        
        # print some feedback:
        if not args.quiet:
            print ' Comparing '+fisher_couples[0]+' and '+fisher_couples[1]+': ',
        
        # get the names:
        fisher_name_1 = fisher_couples[0]
        fisher_name_2 = fisher_couples[1]
        
        # get the fisher matrices:
        fisher_1 = fishers.get_fisher_matrix( names=fisher_name_1 )[0]
        fisher_2 = fishers.get_fisher_matrix( names=fisher_name_2 )[0]
        
        # test for equality:
        if fisher_1 == fisher_2:
            # print feedback:
            if not args.quiet:
                print 'equal.'
        else:
            # print feedback:
            if not args.quiet:
                print 'different. Plotting results.'
            
            # make sure that the two matrices have the same parameters:
            fisher_1_paramnames = fisher_1.get_param_names()
            fisher_2_paramnames = fisher_2.get_param_names()
            paramnames = [ name for name in fisher_1_paramnames if name in fisher_2_paramnames ]
            
            # check if there's parameters left:
            if len(paramnames)==0:
                print 'the two Fisher matrices do not have common parameters.'
                continue
            
            # now reshuffle:
            fisher_1 = fo.reshuffle( fisher_1, paramnames )
            fisher_2 = fo.reshuffle( fisher_2, paramnames )
        
            # compare the two Fisher matrices:
            fisher_comparison = np.abs( fisher_1.get_fisher_matrix()-fisher_2.get_fisher_matrix())/np.abs(fisher_1.get_fisher_matrix())*100.0
        
            # setup the structure of the plot:
            fig, ax = plt.subplots()
            plot_grid = gridspec.GridSpec(1, 1)
            ax1 = plt.subplot(plot_grid[0, 0])
            
            plot_size_temp = float(len(fisher_1.get_param_names()))/10.0*x_size
            fig.set_size_inches(plot_size_temp/2.54 , plot_size_temp/2.54)
            
            # plot:
            fisher_plot  = ax1.imshow( fisher_comparison,
                           vmin=np.floor( 0.9*np.amin( fisher_comparison )),
                           vmax=np.ceil( 1.1*np.amax( fisher_comparison )),
                           interpolation='nearest', cmap=colormap )

            # put the numbers on top:
            min_val, max_val = 0, len(fisher_1.get_param_names())
            ind_array = np.arange(min_val, max_val, 1.0)
            x, y = np.meshgrid(ind_array, ind_array)
                        
            for i, (x_val, y_val) in enumerate(zip(x.flatten(), y.flatten())):
                c = '$'+str( round( fisher_comparison[ int(x_val), int(y_val) ],1 ) )+'\\%$'
                ax1.text(x_val, y_val, c, va='center', ha='center', fontsize=0.8*main_fontsize)
                
            # get numbers and parameter names:
            numbers = np.arange(len(fisher_comparison[0]))
            labels  = fisher_1.get_param_names_latex()
            
            # refactor names:
            labels_temp = []
            for lab in labels:
                labels_temp.append( '$'+lab+'$' )
            labels = labels_temp
    
            # arrange the ticks:
            ax1.xaxis.tick_top()
            ax1.set_xticks( numbers )
            ax1.set_yticks( numbers )
            ax1.set_xticklabels( labels , rotation = 90, fontsize=0.9*main_fontsize)
            ax1.set_yticklabels( labels , fontsize=0.9*main_fontsize )

            # set the color bar:
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes("right", size="5%", pad=0.05)

            cbar = plt.colorbar( fisher_plot, cax=cax, orientation='vertical')

            bar_label = np.round( np.linspace( np.floor( 0.9*np.amin( fisher_comparison )),
                                               np.ceil( 1.1*np.amax( fisher_comparison )),
                                               10 ),
                                               1)
            
            cbar.set_ticks( bar_label )
            cbar.ax.set_yticklabels([ "$"+str(round(i,2))+"\%$" for i in bar_label], fontsize=0.9*main_fontsize )
            
            # align the bar labels:
            fig.canvas.draw()
            cbar.ax.get_yticklabels()[0].set_verticalalignment('bottom')
            cbar.ax.get_yticklabels()[-1].set_verticalalignment('top')
            
            # get the dimensions of the labels:
            fig.canvas.draw()
            # get the default renderer:
            renderer = matplotlib.backend_bases.RendererBase()
            # default dpi is 72.0:
            default_dpi = 72.0
            # get the size of the figure:
            figure_x_size = fig.get_size_inches()[0] #: in inches
            figure_y_size = fig.get_size_inches()[1] #: in inches
            # get the maximum dimensions of the labels:
            max_x_tick_label = [0.0,0.0]
            for xlabel in ax1.get_xticklabels():
                x_dimension = xlabel.get_window_extent(renderer).width/default_dpi   #: in inches
                y_dimension = xlabel.get_window_extent(renderer).height/default_dpi  #: in inches
                if x_dimension>max_x_tick_label[0]: max_x_tick_label[0]=x_dimension
                if y_dimension>max_x_tick_label[1]: max_x_tick_label[1]=y_dimension
        
            # plot title:
            title = 'Comparison of Fisher matrices: \n '+fisher_name_1+' and '+fisher_name_2
            tit   = fig.suptitle( title, fontsize=main_fontsize, x=0.5, y=(matplotlib.rcParams['axes.labelpad']/default_dpi)/figure_y_size, 
                          horizontalalignment='center',
                          verticalalignment='bottom')
            
            # get the size of the title:
            title_y_size = tit.get_window_extent(renderer).height/default_dpi
            
            # get the size of the colorbar:
            max_x_tick_label_cbar = [0.0,0.0]
            for xlabel in cbar.ax.get_yticklabels():
                x_dimension = xlabel.get_window_extent(renderer).width/default_dpi   #: in inches
                y_dimension = xlabel.get_window_extent(renderer).height/default_dpi  #: in inches
                if x_dimension>max_x_tick_label_cbar[0]: max_x_tick_label_cbar[0]=x_dimension
                if y_dimension>max_x_tick_label_cbar[1]: max_x_tick_label_cbar[1]=y_dimension
            
            # update dimensions:
            bottom = ( title_y_size+ 2.0*matplotlib.rcParams['axes.labelpad']/default_dpi)/figure_y_size
            top    = 1.0 -(max_x_tick_label[1]+ 2.0*matplotlib.rcParams['axes.labelpad']/default_dpi)/figure_y_size
            left   = (max_x_tick_label[1]+ 2.0*matplotlib.rcParams['axes.labelpad']/default_dpi)/figure_x_size
            right  = 1.0 -(max_x_tick_label_cbar[0]+ 2.0*matplotlib.rcParams['axes.labelpad']/default_dpi)/figure_y_size
            plot_grid.update( bottom= bottom, top=top, left=left, right=right)
            
            # save the figure and close:
            plt.savefig( outroot+fisher_name_1+'_vs_'+fisher_name_2+'.'+output_format )   
            plt.clf()
            plt.cla()  
            plt.close()
            
            # print some final feedback:
            if not args.quiet:
                print '  Done. Saved results in: ', outroot+fisher_name_1+'_vs_'+fisher_name_2+'.'+output_format
        
    # exit without error:
    exit(0)