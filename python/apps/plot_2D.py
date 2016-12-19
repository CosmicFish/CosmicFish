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

Simple Python code to do 2D plot of a set of Fisher matrices.

The ouput will be a pdf file with the plot.

Invoking the help option ``plot_2D.py -h`` will result in::

    usage: plot_2D.py [-h] [-o OUTROOT] [-p PARAMS [PARAMS ...]] [-s SUM_FISH]
                  [-e] [-f FORMAT] [-d] [-i INI_FILE] [-v] [-q]
                  files [files ...]

    2D confidence level plotter
    
    positional arguments:
      files                 a list of files with Fisher matrices
    
    optional arguments:
      -h, --help            show this help message and exit
      -o OUTROOT, --outroot OUTROOT
                            path and name of the output file. If this option is
                            not present the name and path of the first Fisher
                            matrix will be used.
      -p PARAMS [PARAMS ...], --params PARAMS [PARAMS ...]
                            names of the parameters to plot. If this option is not
                            present all parameters will be used.
      -s SUM_FISH, --sum SUM_FISH
                            decides wether to sum all the input Fisher matrices.
                            If selected the argument will be the label of the
                            Fisher matrices combination.
      -e, --eliminate       If parameters are passed from the command line, this
                            option avoids marginalization over the others.
      -f FORMAT, --format FORMAT
                            format for the output plot.
      -d, --derived         decides wether to look for derived parameters when
                            importing the Fisher matrix
      -i INI_FILE, --ini INI_FILE
                            path and name of file with options
      -v, --version         show program's version number and exit
      -q, --quiet           decides wether something gets printed to screen or not

Developed by Matteo Martinelli (m.martinelli@thphys.uni-heidelberg.de)
and Marco Raveri (mraveri@uchicago.edu) for the CosmicFish code.

"""

# ***************************************************************************************

__version__ = '1.1' #: version of the application

# ***************************************************************************************

""" Hard coded options """

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
import itertools as it

# get the path of the application and the CosmicFish library:
here = os.path.dirname(os.path.abspath(__file__))
cosmicfish_pylib_path = here+'/..'
sys.path.insert(0, os.path.normpath(cosmicfish_pylib_path))

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
    parser = argparse.ArgumentParser(description='2D confidence level plotter')
    # parse file names:
    parser.add_argument('files', metavar='files', type=str, nargs='+',
                         help='a list of files with Fisher matrices')
    # parse the output root:
    parser.add_argument('-o','--outroot', dest='outroot', type=str,
                         help='path and name of the output file. If this option is not present the name and path of the first Fisher matrix will be used.')
    # parse parameters if wanted
    parser.add_argument('-p','--params', dest='params', type=str, nargs='+',
                         help='names of the parameters to plot. If this option is not present all parameters will be used.')
    # parse a sum argument to sum all Fisher matrices:
    parser.add_argument('-s','--sum', dest='sum_fish', type=str,
                        help='decides wether to sum all the input Fisher matrices. If selected the argument will be the label of the Fisher matrices combination.')
    # parse an optional argument to avoid marginalization
    parser.add_argument('-e','--eliminate', dest='eliminate', action='store_true',
                        help='If parameters are passed from the command line, this option avoids marginalization over the others.')
    # parse format argument:
    parser.add_argument('-f','--format', dest='format', type=str,
                         help='format for the output plot.')
    # parse to look for derived parameters:
    parser.add_argument('-d','--derived', dest='derived', action='store_true',
                        help='decides wether to look for derived parameters when importing the Fisher matrix')
    # parse the name of a parameter input file:
    parser.add_argument('-i','--ini', dest='ini_file', type=str,
                        help='path and name of file with options')
    # version:
    parser.add_argument('-v','--version', action='version', version='%(prog)s '+__version__)
    # quiet mode:
    parser.add_argument('-q','--quiet', dest='quiet', action='store_true', 
                        help='decides wether something gets printed to screen or not')
    # do the parsing:
    args = parser.parse_args()
 
    # print the CosmicFish header:
    if not args.quiet:
        fu.CosmicFish_write_header(' 2D confidence level Fisher matrix plotter v'+__version__)

    # process input arguments:
    files          = args.files
    number_fish    = len(files)
    if args.outroot is not None:
        outroot    = args.outroot
    else:
        outroot    = os.path.join( os.path.splitext(files[0])[0] )
    params = args.params

    # warning for the ini file:
    if args.ini_file is not None:
        raise ValueError('Not yet implemented')

    # import Fisher matrices:
    if args.derived:
       fishers = fpa.CosmicFish_FisherAnalysis(fisher_path=files, with_derived=True)
    else:
       fishers = fpa.CosmicFish_FisherAnalysis(fisher_path=files, with_derived=False)

    # manipulating Fisher:
    if args.sum_fish is not None:
        fishers_temp = fpa.CosmicFish_FisherAnalysis()
        fisher_list = fishers.get_fisher_matrix()
        for fish in fisher_list[1:]:
            fisher_list[0] = fisher_list[0]+fish
        fisher_list[0].name = args.sum_fish
        fishers_temp.add_fisher_matrix( fisher_list[0] )
        fishers = fishers_temp

    if args.eliminate:
        if params is None:
            raise ValueError('Avoiding marginalization only works if parameters are specified')
        fishers = fishers.reshuffle( params=params )

    # set up the plot settings
    plot_settings = fps.CosmicFish_PlotSettings()

    # set up the plotter
    plotter = fp.CosmicFishPlotter( settings=plot_settings, fishers=fishers)

    if params is not None:    
        params = [ list(i) for i in it.combinations(params, 2)]
        if len(params)==0:
           raise ValueError('Not enough parameters for 2D plot.')
    
    plotter.new_plot()
    plotter.plot2D( params=params )

    plotter.export( outroot+'.'+args.format )
    plotter.close_plot()
    
    # print some final feedback:
    if not args.quiet:
        print 'Done. Saved results in: ', outroot+'.'+args.format

    # exit without error:
    exit(0)
