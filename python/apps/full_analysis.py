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

Simple Python code to perform analysis of Fisher matrices (plot_1D, plot_2D, bounds...)

The ouput will be a set of pdf with 1D, 2D, triangular plots and a file with bounds


Invoking the help option ``full_analysis.py -h`` will result in::

    usage: full_analysis.py [-h] [-v] [-q] inifile [inifile ...]

    Analysis tool for plot and bounds

    positional arguments:
      inifile        file with a list of instructions

    optional arguments:
      -h, --help     show this help message and exit
      -v, --version  show program's version number and exit
      -q, --quiet    decides wether something gets printed to screen or not

Developed by Matteo Martinelli (martinelli@lorentz.leidenuniv.nl)
and Marco Raveri (mraveri@sissa.it) for the CosmicFish code.

"""

# ***************************************************************************************

__version__ = '1.0' #: version of the application

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
import copy
import itertools as it
import ConfigParser

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
    parser = argparse.ArgumentParser(description='Analysis tool for plot and bounds')
    # parse file names:
    parser.add_argument('inifile', metavar='inifile', type=str, nargs='+',
                         help='file with a list of instructions')
    # version:
    parser.add_argument('-v','--version', action='version', version='%(prog)s '+__version__)
    # quiet mode:
    parser.add_argument('-q','--quiet', dest='quiet', action='store_true', 
                        help='decides wether something gets printed to screen or not')
    # do the parsing:
    args = parser.parse_args()
 

    # print the CosmicFish header:
    if not args.quiet:
        fu.CosmicFish_write_header('Global analysis app version '+__version__)

    # process input arguments:
    inifile          = args.inifile

    #function used to deal with ini sections
    def ConfigSectionMap(section):
        dict1 = {}
        options = Config.options(section)
        for option in options:
            try:
                dict1[option] = Config.get(section, option)
                if dict1[option] == -1:
                   DebugPrint("skip: %s" % option)
            except:
                print("exception on %s!" % option)
                dict1[option] = None
        return dict1
    
    #initializing and reading the config file
    Config = ConfigParser.ConfigParser()
    Config.read(inifile)

    #Reading general options
    outroot   = ConfigSectionMap("General Options")['outroot']
    files     = Config.get("General Options", "fishers").split("\n")
    derived   = Config.getboolean("General Options", "derived")
    sum_fish  = Config.get("General Options", "sum_fish").split("\n")
    eliminate = Config.getboolean("General Options", "eliminate")
    fishnames = Config.get("General Options", "names").split("\n")

    #General screen output
    if not args.quiet:
       print 'GENERAL OPTIONS:'
       print ' Output root='+outroot
       print ' Using derived parameters='+str(derived)
       print ' Eliminate rather than marginalize='+str(eliminate)
       print ' ---------------------------------'
       print ' Bounds from these matrices will be computed:'
       for elem in files:
           print elem
       if sum_fish[0]:
           print 'Also the combination of these will be computed:'
           for elem in sum_fish:
               print elem
       print ' ---------------------------------'
       print
       print

    if not files[0]:
       if not sum_fish[0]:
          print 'NO MATRICES TO WORK WITH!'
          exit()
       else:
          files = sum_fish
          print 'No fishers to plot, using only the combined one'


    #MOD: too much putput here!
    if derived is not False:
       fishers = fpa.CosmicFish_FisherAnalysis(fisher_path=files, with_derived=True)
       if sum_fish[0]:
          print 'NOT HERE'
          summing = fpa.CosmicFish_FisherAnalysis(fisher_path=sum_fish, with_derived=True)
    else:
       fishers = fpa.CosmicFish_FisherAnalysis(fisher_path=files, with_derived=False)
       if sum_fish[0]:
          summing = fpa.CosmicFish_FisherAnalysis(fisher_path=sum_fish, with_derived=False)
   


    fishers_temp = fpa.CosmicFish_FisherAnalysis()
    fisher_list = fishers.get_fisher_matrix()
    if sum_fish[0]:
       summing_list = summing.get_fisher_matrix()
       for fish in summing_list[1:]:
          summing_list[0] = summing_list[0]+fish
       fisher_list.append(summing_list[0])

    for i in range(len(fisher_list)):
       fisher_list[i].name = fishnames[i]


    fishers_temp.add_fisher_matrix( fisher_list[:] )
    fishers = fishers_temp


    #producing 1D plots
    num1D = Config.items( "1Dplot" )
        
    if not args.quiet and len(num1D)>0:
        print
        print 'Producing 1D plots:'
    for key, params in num1D:
        params = Config.get("1Dplot", key).split(",")

        fishers_temp = fishers
        if eliminate is not False:
           fishers_temp = fishers_temp.reshuffle( params=params )

        plot_settings = fps.CosmicFish_PlotSettings()
        plotter = fp.CosmicFishPlotter( settings=plot_settings, fishers=fishers_temp)

        plotter.new_plot()
        plotter.plot1D( params=params )

        plotter.export( outroot+'_1Dplot_'+str(key)+'.pdf' )
        plotter.close_plot()
        if not args.quiet:
           print ' 1D plots done for parameters '+str(params)
           print ' Saved results in: ', outroot+'_1Dplot_'+str(key)+'.pdf'

    if not args.quiet and len(num1D)>0:
        print '1D plots done!'

    #Producing 2D plots
    num2D = Config.items( "2Dplot" )
    if not args.quiet and len(num2D)>0:
        print
        print 'Producing 2D plots:'
    for key, params in num2D:
        params = Config.get("2Dplot", key).split(",")

        fishers_temp = fishers
        if eliminate is not False:
           fishers_temp = fishers_temp.reshuffle( params=params )

        plot_settings = fps.CosmicFish_PlotSettings()
        plotter = fp.CosmicFishPlotter( settings=plot_settings, fishers=fishers_temp)

        if params is not None:
           params = [ list(i) for i in it.combinations(params, 2)]
           if len(params)==0:
              raise ValueError('Not enough parameters for 2D plot.')

        plotter.new_plot()
        plotter.plot2D( params=params )

        plotter.export( outroot+'_2Dplot_'+str(key)+'.pdf' )
        plotter.close_plot()
        if not args.quiet:
           print ' 2D plots done for parameters '+str(params)
           print ' Saved results in: ', outroot+'_2Dplot_'+str(key)+'.pdf'

    if not args.quiet and len(num2D)>0:
        print '2D plots done!'

    #Producing triangular plots
    numtri = Config.items( "triplot" )
    if not args.quiet and len(numtri)>0:
        print 
        print 'Producing triangular plots:'
    for key, params in numtri:
        params = Config.get("triplot", key).split(",")

        fishers_temp = fishers
        if eliminate is not False:
           fishers_temp = fishers_temp.reshuffle( params=params )

        plot_settings = fps.CosmicFish_PlotSettings()
        plotter = fp.CosmicFishPlotter( settings=plot_settings, fishers=fishers_temp)

        plotter.new_plot()
        plotter.plot_tri( params=params )

        plotter.export( outroot+'_triplot_'+str(key)+'.pdf' )
        plotter.close_plot()
        if not args.quiet:
           print ' Triangular plots done for parameters '+str(params)
           print ' Saved results in: ', outroot+'_triplot_'+str(key)+'.pdf'

    if not args.quiet and len(numtri)>0:
        print 'Triangular plots done!'

    #Producing bounds files:
    # get the parameters:
    numbounds     = [ i for i in Config.items( "bounds" ) if "params" in i[0] ]
    use_latex     = Config.getboolean( "bounds",'use_latex')
    latex_num_col = Config.getint( "bounds",'latex_num_col')

    if len(numbounds)>0:
    
        if not args.quiet:
            print
            print 'Producing bounds:'
    
        # open the file if wanted:
        if outroot is not None:
            out_file = open(outroot+'_bounds.txt',"w")
    
        # do some first printing for Latex:
        if use_latex:
            out_file.write( '\\begin{tabular}{ |'+''.join(['l|' for i in range(latex_num_col) ])+' }\n' )
                
        for key, params in numbounds:
            params = Config.get("bounds", key).split(",")
            
            fishers_temp = fishers
            if eliminate is not False:
               fishers_temp = fishers_temp.reshuffle( params=params )
            elif params is not None:
               fishers_temp = fishers_temp.marginalise( params=params )
            
            for num, fish in enumerate(fishers_temp.get_fisher_list()):
                # get the bounds:
                Bounds_68 = list( fu.v_nice_number( fish.get_confidence_bounds( 0.680 ), mode=1 ) )
                Bounds_95 = list( fu.v_nice_number( fish.get_confidence_bounds( 0.950 ), mode=1 ) )
                Bounds_99 = list( fu.v_nice_number( fish.get_confidence_bounds( 0.997 ), mode=1 ) )
                # get the fiducial:
                fiducial = []
                for num, par in enumerate(fish.get_param_fiducial()):
                    fiducial.append( fu.significant_digits( (par, Bounds_68[num]), mode=1 ) )
                # get the names:
                if use_latex:
                    parameter_names_latex = copy.deepcopy(fish.get_param_names_latex())
                else:
                    parameter_names       = copy.deepcopy(fish.get_param_names())
                # do the printing:
                if use_latex:
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
                
            if not args.quiet:
               print ' Bounds computed for parameters '+str( fishers_temp.get_parameter_list() )
        
        # finalize the latex part:
        if use_latex:
            out_file.write( '\hline\n' )
            out_file.write( '\end{tabular}' )
        # close the file:
            out_file.close()  
        
        if not args.quiet:
            print ' Saved results in: ', outroot+'_bounds.txt'
            print 'bounds done!'

    # finalize:
    if not args.quiet:
       print
       print 'It seems everything is done...'
       print 'It was nice working with you. See you soon!'
    # everything's fine, exit without error:
    exit(0)
