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

Simple Python code to compute the narginal bounds of a set of Fisher matrices.

The ouput will be printed on screen if an output root is not provided.

Invoking the help option ``bounds.py -h`` will result in::

    usage: bounds.py [-h] [-o OUTROOT] [-p PARAMS [PARAMS ...]] [-s SUM_FISH] [-e]
                 [-d] [-l] [-i INI_FILE] [-v] [-q]
                 files [files ...]

    confidence bounds calculator
    
    positional arguments:
      files                 a list of files with Fisher matrices
    
    optional arguments:
      -h, --help            show this help message and exit
      -o OUTROOT, --outroot OUTROOT
                            path and name of the output file. If this option is
                            not present the name and path of the first Fisher
                            matrix will be used
      -p PARAMS [PARAMS ...], --params PARAMS [PARAMS ...]
                            names of the parameters to plot. If this option is not
                            present the plot will contain all parameters.
      -s SUM_FISH, --sum SUM_FISH
                            decides wether to sum all the input Fisher matrices.
                            If selected the argument will be the label of the
                            Fisher matrices combination
      -e, --eliminate       if parameters are passed from the command line this
                            option avoids marginalization over the others
      -d, --derived         decides wether to look for derived parameters when
                            importing the Fisher matrix
      -l, --latex           decides wether output the bounds in LaTeX format
      -i INI_FILE, --ini INI_FILE
                            name and path of an input file with options
      -v, --version         show program's version number and exit
      -q, --quiet           decides wether something gets printed to screen or not

Developed by Marco Raveri (mraveri@sissa.it) 
and Matteo Martinelli (m.martinelli@thphys.uni-heidelberg.de) for the CosmicFish code.

"""

# ***************************************************************************************

__version__ = '1.0' #: version of the application

# ***************************************************************************************

""" Hard coded options """

latex_num_col = 3

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
    parser = argparse.ArgumentParser(description='confidence bounds calculator')
    # parse file names:
    parser.add_argument('files', metavar='files', type=str, nargs='+',
                         help='a list of files with Fisher matrices')
    # parse the output root:
    parser.add_argument('-o','--outroot', dest='outroot', type=str,
                         help='path and name of the output file. If this option is not present the name and path of the first Fisher matrix will be used')
    # parse parameters if wanted:
    parser.add_argument('-p','--params', dest='params', type=str, nargs='+',
                         help='names of the parameters to plot. If this option is not present the plot will contain all parameters.')
    # parse a sum argument to sum all Fisher matrices:
    parser.add_argument('-s','--sum', dest='sum_fish', type=str,
                        help='decides wether to sum all the input Fisher matrices. If selected the argument will be the label of the Fisher matrices combination')
    # parse an optional argument to avoid marginalization:
    parser.add_argument('-e','--eliminate', dest='eliminate', action='store_true', 
                        help='if parameters are passed from the command line this option avoids marginalization over the others')
    # parse a sum argument to sum all Fisher matrices:
    parser.add_argument('-d','--derived', dest='derived', action='store_true', 
                        help='decides wether to look for derived parameters when importing the Fisher matrix')
    # parse an argument to output the table in LaTeX format:
    parser.add_argument('-l','--latex', dest='latex', action='store_true', 
                        help='decides wether output the bounds in LaTeX format')
    
    # parse the name of a parameter input file:
    parser.add_argument('-i','--ini', dest='ini_file', type=str,
                        help='name and path of an input file with options')
    # version:
    parser.add_argument('-v','--version', action='version', version='%(prog)s '+__version__)
    # quiet mode:
    parser.add_argument('-q','--quiet', dest='quiet', action='store_true', 
                        help='decides wether something gets printed to screen or not')
    # do the parsing:
    args = parser.parse_args()
    
    # print the CosmicFish header:
    if not args.quiet:
        fu.CosmicFish_write_header(' confidence bounds calculator version '+__version__)
        
    # process input arguments:
    files          = args.files
    number_fish    = len(files)
    if args.outroot is not None:
        outroot    = args.outroot
    else:
        outroot    = None
    params = args.params
    # warning for the ini file:
    if args.ini_file is not None:
        raise ValueError('Not yet implemented.')
    
    # import the Fisher matrices:
    if args.derived:
        fishers = fpa.CosmicFish_FisherAnalysis(fisher_path=files, with_derived=True)
    else:
        fishers = fpa.CosmicFish_FisherAnalysis(fisher_path=files, with_derived=False)
    # sum all the fishers if wanted:
    if args.sum_fish is not None:
        fishers_temp = fpa.CosmicFish_FisherAnalysis()
        fisher_list  = fishers.get_fisher_matrix()
        for fish in fisher_list[1:]:
            fisher_list[0] = fisher_list[0]+fish
        fisher_list[0].name = args.sum_fish
        fishers_temp.add_fisher_matrix( fisher_list[0] )
        fishers = fishers_temp
    
    # eliminate parameters that are not wanted:
    if args.eliminate:
        if params is None:
            raise ValueError('Avoiding marginalization works only if parameters are specified')
        fishers = fishers.reshuffle( params=params )
    elif params is not None:
        fishers = fishers.marginalise( params=params )
    
    # open the file if wanted:
    if outroot is not None:
        out_file = open(outroot,"w")
    
    # do some first printing for Latex:
    if args.latex:
        if outroot is not None:
            out_file.write( '\\begin{tabular}{ |'+''.join(['l|' for i in range(latex_num_col) ])+' }\n' )
        else:
            print '\\begin{tabular}{ |'+''.join(['l|' for i in range(latex_num_col) ])+' }'
            
    # print the confidence bounds on the Fisher matrices:
    for num, fish in enumerate(fishers.get_fisher_list()):
        # get the bounds:
        Bounds_68 = list( fu.v_nice_number( fish.get_confidence_bounds( 0.680 ), mode=1 ) )
        Bounds_95 = list( fu.v_nice_number( fish.get_confidence_bounds( 0.950 ), mode=1 ) )
        Bounds_99 = list( fu.v_nice_number( fish.get_confidence_bounds( 0.997 ), mode=1 ) )
        # get the fiducial:
        fiducial = []
        for num, par in enumerate(fish.get_param_fiducial()):
            fiducial.append( fu.significant_digits( (par, Bounds_68[num]), mode=1 ) )
        # get the names:
        if args.latex:
            parameter_names_latex = fish.get_param_names_latex()
        else:
            parameter_names       = fish.get_param_names()
        # do the printing:
        if args.latex:
            print_table = []
            for par,fid,bound in zip( parameter_names_latex,fiducial,Bounds_68 ):
                print_table.append( '$'+str(par)+' = '+str(fid)+' \pm '+str(bound)+'$' )
            if len(print_table)%latex_num_col == 0:
                table_length = len(print_table)/latex_num_col
            else:
                table_length = len(print_table)/latex_num_col +1
            print_table = fu.grouper( table_length, print_table,fillvalue='' )
            print_table = [ list(i) for i in print_table]      
            print_table = map(list, zip(*print_table))
            col_width = [max(len(str(x)) for x in col) for col in zip(*print_table)]
            
            if outroot is not None:
                out_file.write( '\hline'+'\n' )
                out_file.write( '\multicolumn{'+str(latex_num_col)+'}{|c|}{'+fish.name.replace('_',' ')+'} \\\[1mm]'+'\n' )
                out_file.write( '\hline'+'\n' )
                for line in print_table:
                    out_file.write( " " + " & ".join("{:{}}".format(x, col_width[i]) for i, x in enumerate(line)) + " \\\[1mm]\n" )
            else:
                print '\hline'
                print '\multicolumn{'+str(latex_num_col)+'}{|c|}{'+fish.name.replace('_',' ')+'} \\\[1mm]'
                print '\hline'
                for line in print_table:
                    print " " + " & ".join("{:{}}".format(x, col_width[i]) for i, x in enumerate(line)) + " \\\[1mm]"
                
        else:
            # put on top the labels of the columns:
            parameter_names.insert(0,' Parameter ')
            fiducial.insert(0, ' fiducial')
            Bounds_68.insert(0,' 68% C.L.')
            Bounds_95.insert(0,' 95% C.L.')
            Bounds_99.insert(0,' 99.7% C.L.')
            # put a white space:
            parameter_names.insert(1,' ')
            fiducial.insert(1, ' ')
            Bounds_68.insert(1,' ')
            Bounds_95.insert(1,' ')
            Bounds_99.insert(1,' ')
            #
            print_table = [parameter_names,fiducial,Bounds_68, Bounds_95, Bounds_99]
            if outroot is not None:
                out_file.write( ''.join([ '*' for i in xrange(len('Parameter bounds for the Fisher matrix: '+fish.name)+1)])+'\n' )
                out_file.write( 'Parameter bounds for the Fisher matrix: '+fish.name+'\n' )
                out_file.write( ''.join([ '*' for i in xrange(len('Parameter bounds for the Fisher matrix: '+fish.name)+1)])+'\n' )
                out_file.write( '\n' )
                print_table = map(list, zip(*print_table))
                col_width = [max(len(str(x)) for x in col) for col in zip(*print_table)]
                # print it to file:
                for line in print_table:
                    out_file.write( "| " + " | ".join("{:{}}".format(x, col_width[i]) for i, x in enumerate(line)) + " |"+'\n' )
                out_file.write( '\n' )
            else:
                # print the table:
                # header to screen
                print ''.join([ '*' for i in xrange(len('Parameter bounds for the Fisher matrix: '+fish.name)+1)])
                print 'Parameter bounds for the Fisher matrix: ', fish.name
                print ''.join([ '*' for i in xrange(len('Parameter bounds for the Fisher matrix: '+fish.name)+1)])
                fu.print_table( print_table )
        
    # finalize and close file:
    if args.latex:
        if outroot is not None:
            out_file.write( '\hline\n' )
            out_file.write( '\end{tabular}' )
        else:
            print '\hline'
            print '\end{tabular}'
        
    if outroot is not None:
        out_file.close()    
    
    # print some final feedback:
    if not args.quiet:
        if outroot is not None:
            print 'Done. Saved results in: ', outroot+'.txt'
        else:
            print 'Done.'
    # exit without error:
    exit(0)