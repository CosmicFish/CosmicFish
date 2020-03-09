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

# import first packages
import os, sys, copy
# define path to the executable and add all the relevant folders to the path where python looks for files.
here = os.path.dirname(os.path.abspath(__file__))
cosmicfish_pylib_path = here+'/..'
test_input = here+'/test_input'
sys.path.insert(0, os.path.normpath(cosmicfish_pylib_path))

import numpy as np
import cosmicfish_pylib.fisher_matrix as fm
import cosmicfish_pylib.fisher_plot_analysis as fpa
import cosmicfish_pylib.fisher_plot as fp
import cosmicfish_pylib.fisher_plot_settings as fps

from nose.tools import with_setup
from nose.tools import assert_raises

import cosmicfish_pylib.colors as fc
color_print = fc.bash_colors()

plot_dump = here+'/test_plot/'

# ***************************************************************************************
    
class test_setup():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_setup.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_setup.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        num_param = 6
        num_fish  = 3
        matrix = np.identity(num_param)
        for j in xrange(num_fish):
            for i in xrange(num_param):
                #matrix[i,i] = i+j+1
                matrix[i,i] = (i+1)/float(j+1)
            param_names_latex = [ 'm'+str(i) for i in xrange(num_param) ] 
            fiducial = [ float(i) for i in xrange(num_param) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test
        
    def test_init(self):
        set = fps.CosmicFish_PlotSettings(do_legend=False)
        test_plotter_1 = fp.CosmicFishPlotter( fishers=self.fisher_list_test, settings=set )
        assert test_plotter_1.plot_settings.do_legend == False
        test_plotter_2 = fp.CosmicFishPlotter( fishers=self.fisher_list_test, settings={'do_legend':False} )
        assert test_plotter_2.plot_settings.do_legend == False
        
    def test_init_raises(self):
        assert_raises( ValueError, fp.CosmicFishPlotter, fishers=self.fisher_list_test, settings='test' )
        assert_raises( ValueError, fp.CosmicFishPlotter, fishers='test' )
    
    def test_figure_1D_raises(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        place = test_plotter.plot_grid[0,0]
        assert_raises( ValueError, test_plotter.figure_1D, subplot=place, param='test', names=['fisher1'] )
        assert_raises( ValueError, test_plotter.figure_1D, subplot=place, param='p1', names=['fisher1','test'] )
        test_plotter.close_plot()
    
    def test_figure_2D_raises(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        test_plotter.new_plot()
        test_plotter.plot2D([['p1','p2']])
        place = test_plotter.plot_grid[0,0]
        assert_raises( ValueError, test_plotter.figure_2D, subplot=place, param1='test', param2='p2', names=['fisher1'] )
        assert_raises( ValueError, test_plotter.figure_2D, subplot=place, param1='p1', param2='test', names=['fisher1'] )
        assert_raises( ValueError, test_plotter.figure_2D, subplot=place, param1='p1', param2='p2', names=['test'] )
        test_plotter.close_plot()

# ***************************************************************************************
    
class test_utilities():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_utilities.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_utilities.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        num_param = 6
        num_fish  = 3
        matrix = np.identity(num_param)
        for j in xrange(num_fish):
            for i in xrange(num_param):
                #matrix[i,i] = i+j+1
                matrix[i,i] = (i+1)/float(j+1)
            param_names_latex = [ 'm'+str(i) for i in xrange(num_param) ] 
            fiducial = [ float(i) for i in xrange(num_param) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test
        
    def test_setting_setter_raises(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        assert_raises( ValueError, test_plotter.setting_setter, 'test')
        
    def test_bind_plot_settings_to_names(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        test_plotter.bind_plot_settings_to_names( names='fisher1', 
                                                  line_colors=['k'], 
                                                  solid_colors=['k'],
                                                  labels = ['test'],
                                                  linestyle = ['-'] )
        assert test_plotter.bind_linestyle == {'fisher1': '-'}
        assert test_plotter.bind_labels == {'fisher1': 'test'}
        assert test_plotter.bind_solid_colors == {'fisher1': 'k'}
        assert test_plotter.bind_line_colors == {'fisher1': 'k'}
        
# ***************************************************************************************
    
class test_plot1D():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_plot1D.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_plot1D.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        num_param = 6
        num_fish  = 3
        matrix = np.identity(num_param)
        for j in xrange(num_fish):
            for i in xrange(num_param):
                #matrix[i,i] = i+j+1
                matrix[i,i] = (i+1)/float(j+1)
            param_names_latex = [ 'm'+str(i) for i in xrange(num_param) ] 
            fiducial = [ float(i) for i in xrange(num_param) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test
        
    def test_init_plot_single_1D(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 1 parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        test_plotter.export(plot_dump+'test_init_plot_single_1D_1.pdf')
        test_plotter.close_plot()
        # get plot 1D of 1 parameter and some fishers:
        test_plotter.new_plot()
        test_plotter.plot1D('p1',['fisher1','fisher3'])
        test_plotter.export(plot_dump+'test_init_plot_single_1D_2.pdf')
        test_plotter.close_plot()
        
    def test_init_plot_multiple_1D(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 2 parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot1D(['p1','p2','par1'])
        test_plotter.export(plot_dump+'test_init_plot_multiple_1D_1.pdf')
        test_plotter.close_plot()
        # get plot 1D of 2 parameter and some fishers:
        test_plotter.new_plot()
        test_plotter.plot1D(['p1','p3'],['fisher1','fisher3'])
        test_plotter.export(plot_dump+'test_init_plot_multiple_1D_2.pdf')
        test_plotter.close_plot()
        # get plot 1D of all parameters and some fishers:
        test_plotter.new_plot()
        test_plotter.plot1D(names=['fisher1','fisher3'])
        test_plotter.export(plot_dump+'test_init_plot_multiple_1D_3.pdf')
        test_plotter.close_plot()
        # get plot 1D of all parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot1D(num_plots_per_line=2)
        test_plotter.export(plot_dump+'test_init_plot_multiple_1D_4.pdf')
        test_plotter.close_plot()
        # assert raises:
        test_plotter.new_plot()
        assert_raises( ValueError, test_plotter.plot1D, ['p1'],names=['fisher_3'] )
        test_plotter.close_plot()
        test_plotter.new_plot()
        assert_raises( ValueError, test_plotter.plot1D, ['p_1'],names=['fisher3'] )
        test_plotter.close_plot()
        
    def test_plot_1D_options(self):
        # get plot 1D of 1 parameter and all fishers:
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, D1_filled=False, D1_show_x_ticks=False, D1_show_y_ticks=False)
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        test_plotter.export(plot_dump+'test_plot_1D_options_1.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, D1_show_x_ticks_labels=False, D1_show_y_ticks_labels=False )
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        test_plotter.export(plot_dump+'test_plot_1D_options_2.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, D1_show_xaxis_label=False, D1_show_yaxis_label=False )
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        test_plotter.export(plot_dump+'test_plot_1D_options_3.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, D1_norm_prob=True )
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        test_plotter.export(plot_dump+'test_plot_1D_options_4.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, do_legend=False )
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        test_plotter.export(plot_dump+'test_plot_1D_options_5.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, use_fixed_figure_width=True )
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        test_plotter.export(plot_dump+'test_plot_1D_options_6.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, use_fixed_figure_height=True )
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        test_plotter.export(plot_dump+'test_plot_1D_options_7.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, use_fixed_figure_width=True, use_fixed_figure_height=True )
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        test_plotter.export(plot_dump+'test_plot_1D_options_8.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, D1_ylabel_right=True, D1_xlabel_up=True )
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        test_plotter.export(plot_dump+'test_plot_1D_options_9.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, legend_takes_place_plot=True )
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        test_plotter.export(plot_dump+'test_plot_1D_options_10.pdf')
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, legend_takes_place_plot=False )
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        test_plotter.export(plot_dump+'test_plot_1D_options_11.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, tight_layout=True )
        test_plotter.new_plot()
        test_plotter.plot1D('p1')
        test_plotter.export(plot_dump+'test_plot_1D_options_12.pdf')
        test_plotter.close_plot()

        # test for extra_vlines (single)
        pnames = self.fisher_list_test.get_fisher_list()[0].get_param_names()
        pvalues = self.fisher_list_test.get_fisher_list()[0].get_param_fiducial()
        extra_vlines = {pname : pvalue * 1.5 for pname, pvalue in zip(pnames, pvalues)}
        extra_vlines[None] = {
            'linestyle' : '--',
            'color' : 'green',
            'linewidth' : 3,
        }
        test_plotter = fp.CosmicFishPlotter(
            fishers=self.fisher_list_test,
            D1_extra_vlines=extra_vlines,
        )
        test_plotter.new_plot()
        test_plotter.plot1D('p2')
        test_plotter.export(plot_dump+'test_plot_1D_options_13.pdf')
        test_plotter.close_plot()

        # test for extra_vlines (multiple)
        extra_vlines = [
            {pname : pvalue * (1 + index / 3.) for pname, pvalue in zip(pnames, pvalues)} for index in range(3)
        ]
        extra_vlines[0][None] = {
            'linestyle' : '--',
            'color' : 'green',
            'linewidth' : 2,
        }
        extra_vlines[1][None] = {
            'linestyle' : ':',
            'color' : 'orange',
            'linewidth' : 2,
        }
        extra_vlines[2][None] = {
            'linestyle' : '-.',
            'color' : 'magenta',
            'linewidth' : 2,
        }

        test_plotter = fp.CosmicFishPlotter(
            fishers=self.fisher_list_test,
            D1_extra_vlines=extra_vlines,
        )
        test_plotter.new_plot()
        test_plotter.plot1D('p2')
        test_plotter.export(plot_dump+'test_plot_1D_options_14.pdf')
        test_plotter.close_plot()

        
# ***************************************************************************************
    
class test_plot2D():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_plot2D.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_plot2D.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        num_param = 3
        num_fish  = 3
        matrix = np.identity(num_param)
        for j in xrange(num_fish):
            for i in xrange(num_param):
                #matrix[i,i] = i+j+1
                matrix[i,i] = (i+1)/float(j+1)
            param_names_latex = [ 'm'+str(i) for i in xrange(num_param) ] 
            fiducial = [ float(i) for i in xrange(num_param) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test
        
    def test_init_plot_single_2D(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 1 parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot2D([['p1','p2']])
        test_plotter.export(plot_dump+'test_init_plot_single_2D_1.pdf')
        test_plotter.close_plot()
        # get plot 1D of 1 parameter and some fishers:
        test_plotter.new_plot()
        test_plotter.plot2D([['p1','p2']],['fisher1','fisher3'])
        test_plotter.export(plot_dump+'test_init_plot_single_2D_2.pdf')
        test_plotter.close_plot()
        
    def test_init_plot_multiple_2D(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 2 parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot2D([['p1','p2'],['par1','p2'],['p3','p2']])
        test_plotter.export(plot_dump+'test_init_plot_multiple_2D_1.pdf')
        test_plotter.close_plot()
        # get plot 1D of 2 parameter and some fishers:
        test_plotter.new_plot()
        test_plotter.plot2D([['p1','p2'],['par1','p2'],['p3','p2']],['fisher1','fisher3'])
        test_plotter.export(plot_dump+'test_init_plot_multiple_2D_2.pdf')
        test_plotter.close_plot()
        # get plot 1D of 2 parameter and some fishers:
        test_plotter.new_plot()
        test_plotter.plot2D(names=['fisher1'])
        test_plotter.export(plot_dump+'test_init_plot_multiple_2D_3.pdf')
        test_plotter.close_plot()
        # get plot 1D of all parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot2D()
        test_plotter.export(plot_dump+'test_init_plot_multiple_2D_4.pdf')
        test_plotter.close_plot()
        # assert raises:
        test_plotter.new_plot()
        assert_raises( ValueError, test_plotter.plot2D, [['p1','p2']],names=['fisher_3'] )
        test_plotter.close_plot()
        test_plotter.new_plot()
        assert_raises( ValueError, test_plotter.plot2D, [['p_1','p_2']],names=['fisher3'] )
        test_plotter.close_plot()
    
    def test_plot_2D_options(self):
        # get plot 1D of 1 parameter and all fishers:
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, D2_show_x_ticks=False, D2_show_y_ticks=False )
        test_plotter.new_plot()
        test_plotter.plot2D([['p1','p2']])
        test_plotter.export(plot_dump+'test_plot_2D_options_1.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, D2_show_x_ticks_labels=False, D2_show_y_ticks_labels=False )
        test_plotter.new_plot()
        test_plotter.plot2D([['p1','p2']])
        test_plotter.export(plot_dump+'test_plot_2D_options_2.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, D2_show_xaxis_label=False, D2_show_yaxis_label=False )
        test_plotter.new_plot()
        test_plotter.plot2D([['p1','p2']])
        test_plotter.export(plot_dump+'test_plot_2D_options_3.pdf')
        test_plotter.close_plot()
    
    def test_plot_2D_legend(self):
        # get plot 1D of 1 parameter and all fishers:
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test)
        legend_locs = ['upper right','upper left','lower left','lower right','right','center left','center right','lower center','upper center','center'] 
        for num,place in enumerate(legend_locs):
            test_plotter.new_plot()
            test_plotter.plot2D([['p1','p2']], legend_loc=place)
            test_plotter.export(plot_dump+'test_plot_2D_legend_'+str(num+1)+'.pdf')
            test_plotter.close_plot()
        
# ***************************************************************************************
    
class test_plot3D():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_plot3D.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_plot3D.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        num_param = 4
        num_fish  = 1
        matrix = np.identity(num_param)
        for j in xrange(num_fish):
            for i in xrange(num_param):
                #matrix[i,i] = i+j+1
                matrix[i,i] = (i+1)/float(j+1)
            param_names_latex = [ 'm'+str(i) for i in xrange(num_param) ] 
            fiducial = [ float(i) for i in xrange(num_param) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test
        
    def test_init_plot_single_3D(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 1 parameter and all fishers:
        test_plotter.new_plot()
        assert_raises( ValueError, test_plotter.plot3D, [['p1','p2','p3']] )
        test_plotter.export(plot_dump+'test_init_plot_single_3D_1.pdf')
        test_plotter.close_plot()
        # get plot 1D of 1 parameter and some fishers:
        test_plotter.new_plot()
        assert_raises( ValueError, test_plotter.plot3D, [['p1','p2','p3']], ['fisher1','fisher3'] )
        test_plotter.export(plot_dump+'test_init_plot_single_3D_2.pdf')
        test_plotter.close_plot()
    
    def test_init_plot_multiple_3D(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 2 parameter and all fishers:
        test_plotter.new_plot()
        assert_raises( ValueError, test_plotter.plot3D, [['p1','p2','p3'],['par1','p2','p3'],['p3','p2','p4']] )
        test_plotter.export(plot_dump+'test_init_plot_multiple_3D_1.pdf')
        test_plotter.close_plot()

# ***************************************************************************************
    
class test_plot_tri():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_plot_tri.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_plot_tri.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        num_param = 4
        num_fish  = 3
        matrix = np.identity(num_param)
        for j in xrange(num_fish):
            for i in xrange(num_param):
                #matrix[i,i] = i+j+1
                matrix[i,i] = (i+1)/float(j+1)
            param_names_latex = [ 'm'+str(i) for i in xrange(num_param) ] 
            fiducial = [ float(i) for i in xrange(num_param) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test
        
    def test_init_plot_tri(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 2 parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p2','par1'])
        test_plotter.export(plot_dump+'test_init_plot_tri_1.pdf')
        test_plotter.close_plot()
        # get plot 1D of 2 parameter and some fishers:
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p3'],['fisher1','fisher3'])
        test_plotter.export(plot_dump+'test_init_plot_tri_2.pdf')
        test_plotter.close_plot()
        # get plot 1D of all parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot_tri(names=['fisher1','fisher3'])
        test_plotter.export(plot_dump+'test_init_plot_tri_3.pdf')
        test_plotter.close_plot()
        # get plot 1D of all parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot_tri()
        test_plotter.export(plot_dump+'test_init_plot_tri_4.pdf')
        test_plotter.close_plot()
        # assert raises:
        test_plotter.new_plot()
        assert_raises( ValueError, test_plotter.plot_tri, ['p1','p2'],names=['fisher_3'] )
        test_plotter.close_plot()
        test_plotter.new_plot()
        assert_raises( ValueError, test_plotter.plot_tri, ['p_1','p_2'],names=['fisher3'] )
        test_plotter.close_plot()
    
    def test_plot_tri_options(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, legend_takes_place_plot_tri=False )
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p2'])
        test_plotter.export(plot_dump+'test_plot_tri_options_1.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, use_fixed_figure_height=True )
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p2'])
        test_plotter.export(plot_dump+'test_plot_tri_options_2.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, use_fixed_figure_width=True )
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p2'])
        test_plotter.export(plot_dump+'test_plot_tri_options_3.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, use_fixed_figure_width=True, use_fixed_figure_height=True )
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p2'])
        test_plotter.export(plot_dump+'test_plot_tri_options_4.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p2'])
        for subplot in test_plotter.plot_dict.values():
            subplot.yaxis.set_label_position('right')
            subplot.xaxis.set_label_position('top')
        test_plotter.set_triplot_dimensions( num_col=2, num_rows=2 )
        test_plotter.export(plot_dump+'test_plot_tri_options_5.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p2'])
        for subplot in test_plotter.plot_dict.values():
            subplot.yaxis.set_label_position('left')
            subplot.xaxis.set_label_position('bottom')
        test_plotter.set_triplot_dimensions( num_col=2, num_rows=2 )
        test_plotter.export(plot_dump+'test_plot_tri_options_6.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, do_legend=False )
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p2'])
        test_plotter.export(plot_dump+'test_plot_tri_options_7.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p2'])
        del test_plotter.legend
        test_plotter.set_legend()
        test_plotter.export(plot_dump+'test_plot_tri_options_8.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, legend_filled=False )
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p2'])
        test_plotter.export(plot_dump+'test_plot_tri_options_9.pdf')
        test_plotter.close_plot()
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, tight_layout=True )
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p2'])
        test_plotter.export(plot_dump+'test_plot_tri_options_10.pdf')
        test_plotter.close_plot()

        # test for extra_vlines (single)
        pnames = self.fisher_list_test.get_fisher_list()[0].get_param_names()
        pvalues = self.fisher_list_test.get_fisher_list()[0].get_param_fiducial()
        extra_vlines = {pname : pvalue * 1.5 for pname, pvalue in zip(pnames, pvalues)}
        extra_vlines[None] = {
            'linestyle' : '--',
            'color' : 'green',
            'linewidth' : 3,
        }
        test_plotter = fp.CosmicFishPlotter(
            fishers=self.fisher_list_test,
            D1_extra_vlines=extra_vlines,
        )
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p2'])
        test_plotter.export(plot_dump+'test_plot_tri_options_11.pdf')
        test_plotter.close_plot()

        # test for extra_vlines (multiple)
        extra_vlines = [
            {pname : pvalue * (1 + index / 3.) for pname, pvalue in zip(pnames, pvalues)} for index in range(3)
        ]
        extra_vlines[0][None] = {
            'linestyle' : '--',
            'color' : 'green',
            'linewidth' : 2,
        }
        extra_vlines[1][None] = {
            'linestyle' : ':',
            'color' : 'orange',
            'linewidth' : 2,
        }
        extra_vlines[2][None] = {
            'linestyle' : '-.',
            'color' : 'magenta',
            'linewidth' : 2,
        }

        test_plotter = fp.CosmicFishPlotter(
            fishers=self.fisher_list_test,
            D1_extra_vlines=extra_vlines,
        )
        test_plotter.new_plot()
        test_plotter.plot_tri(['p1','p2'])
        test_plotter.export(plot_dump+'test_plot_tri_options_12.pdf')
        test_plotter.close_plot()
    
    def test_plot_tri_legend(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test, legend_takes_place_plot_tri=False )
        legend_locs = ['upper right','upper left','lower left','lower right','right','center left','center right','lower center','upper center','center'] 
        for num,place in enumerate(legend_locs):
            test_plotter.new_plot()
            test_plotter.plot_tri(['p1','p2'],legend_loc=place)
            test_plotter.export(plot_dump+'test_plot_tri_legend_'+str(num+1)+'.pdf')
            test_plotter.close_plot()
 
# ***************************************************************************************
    
class test_plot_mixed():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_plot_mixed.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_plot_mixed.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        num_param = 4
        num_fish  = 3
        matrix = np.identity(num_param)
        for j in xrange(num_fish):
            for i in xrange(num_param):
                #matrix[i,i] = i+j+1
                matrix[i,i] = (i+1)/float(j+1)
            param_names_latex = [ 'm'+str(i) for i in xrange(num_param) ] 
            fiducial = [ float(i) for i in xrange(num_param) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test
        
    def test_init_plot_mixed(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 2 parameter and all fishers:
        test_plotter.new_plot()
        assert_raises( ValueError, test_plotter.plot_mixed, [['p1','p2'],['p1'],['p1','p2','p3']] )
        test_plotter.export(plot_dump+'test_init_plot_mixed_1.pdf')
        test_plotter.close_plot()
        
# ***************************************************************************************

class test_plot_realistic():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_plot_realistic.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_plot_realistic.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_1      = fm.fisher_matrix( file_name=test_input+'/dummy_fisher_matrix.dat'  )
        fisher_1.name = 'realistic fisher'
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_1)
    
    def teardown(self):
        del self.fisher_list_test 
        
    def test_init_plot_single_1D(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 1 parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot1D('omegabh2', num_plots_per_line=2, title="Wonderfull plot" )
        test_plotter.export(plot_dump+'test_realistic_1.png')
        test_plotter.export(plot_dump+'test_realistic_1.pdf')
        test_plotter.close_plot()
        
    def test_init_plot_multiple_1D(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 2 parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot1D(['omegach2','omeganuh2','h'])
        test_plotter.export(plot_dump+'test_realistic_2.pdf')
        test_plotter.close_plot()
        # get plot 1D of all parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot1D()
        test_plotter.export(plot_dump+'test_realistic_3.pdf')
        test_plotter.close_plot()
    
    def test_init_plot_single_2D(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 1 parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot2D([['yhe','logA']])
        test_plotter.export(plot_dump+'test_realistic_4.pdf')
        test_plotter.close_plot()
        
    def test_init_plot_multiple_2D(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 2 parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot2D([['ns','nrun'],['nt','r'],['tau','Bias_W_1']],  title="Wonderfull plot")
        test_plotter.export(plot_dump+'test_realistic_5.pdf')
        test_plotter.close_plot()
        
    def test_init_plot_tri(self):
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 2 parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot_tri(['omegabh2','omegach2','omeganuh2'], title="Wonderfull plot")
        test_plotter.export(plot_dump+'test_realistic_6.pdf')
        test_plotter.export(plot_dump+'test_realistic_6.png')
        test_plotter.close_plot()
    
    def test_init_plot_tri_2(self):
        # do spectral protection a couple of times:
        for i in xrange(4):
            fish = copy.deepcopy( self.fisher_list_test.get_fisher_list()[-1] )
            fish.protect_degenerate( cache=True )
            fish.name = str(i+1)
            self.fisher_list_test.add_fisher_matrix( fish )
        test_plotter = fp.CosmicFishPlotter( fishers=self.fisher_list_test )
        # get plot 1D of 2 parameter and all fishers:
        test_plotter.new_plot()
        test_plotter.plot_tri(['omegabh2','omegach2','omeganuh2','h','logA','ns'], title="Wonderfull plot")
        test_plotter.export(plot_dump+'test_realistic_7.png')
        test_plotter.close_plot()
        
        
# ***************************************************************************************