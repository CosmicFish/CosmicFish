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
import os, sys
# define path to the executable and add all the relevant folders to the path where python looks for files.
here = os.path.dirname(os.path.abspath(__file__))
cosmicfish_pylib_path = here+'/..'
test_input = here+'/test_input'
sys.path.insert(0, os.path.normpath(cosmicfish_pylib_path))

import numpy as np

import cosmicfish_pylib.fisher_matrix as fm
import cosmicfish_pylib.fisher_plot_analysis as fpa

from nose.tools import with_setup
from nose.tools import assert_raises

import cosmicfish_pylib.colors as fc
color_print = fc.bash_colors()


# ***************************************************************************************
    
class test_class_getters():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_class_getters.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_class_getters.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        matrix = np.identity(3)
        for j in xrange(3):
            for i in xrange(3):
                matrix[i,i] = i+j+1
            param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
            fiducial = [ float(i) for i in xrange(3) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test 
        
    # test getters:
    def test_get_fisher_name_list(self):
        assert self.fisher_list_test.get_fisher_name_list() == ['fisher1', 'fisher2', 'fisher3']
    
    def test_get_fisher_list(self):
        assert isinstance( self.fisher_list_test.get_fisher_list(), list )
        for i in self.fisher_list_test.get_fisher_list():
            assert isinstance( i, fm.fisher_matrix )
        
    def test_get_fisher_matrix(self):
        assert  [ i.name for i in self.fisher_list_test.get_fisher_matrix() ] == ['fisher1', 'fisher2', 'fisher3']
        assert [ i.name for i in self.fisher_list_test.get_fisher_matrix('fisher1') ] == ['fisher1']
        assert [ i.name for i in self.fisher_list_test.get_fisher_matrix(['fisher1','fisher2']) ] == ['fisher1', 'fisher2']
    
    def test_delete_fisher_matrix_1(self):
        self.fisher_list_test.delete_fisher_matrix(['fisher1','fisher2'])
        assert self.fisher_list_test.get_fisher_name_list() == ['fisher3'] 
    
    def test_delete_fisher_matrix_2(self):
        self.fisher_list_test.delete_fisher_matrix(['fisher0','fisher2'])
        assert self.fisher_list_test.get_fisher_name_list() == ['fisher1','fisher3'] 
        
    def test_delete_fisher_matrix_3(self):
        self.fisher_list_test.delete_fisher_matrix()
        assert self.fisher_list_test.get_fisher_name_list() == [] 

    def test_delete_fisher_matrix_4(self):
        self.fisher_list_test.delete_fisher_matrix(['fisher2','fisher2'])
        assert self.fisher_list_test.get_fisher_name_list() == ['fisher1','fisher3'] 

    def test_delete_fisher_matrix_5(self):
        self.fisher_list_test.delete_fisher_matrix(['fisher0','fisher0'])
        assert self.fisher_list_test.get_fisher_name_list() == ['fisher1','fisher2','fisher3'] 
    
    def test_get_parameter_list(self):
        assert self.fisher_list_test.get_parameter_list() == ['p1', 'p2', 'p3']

# ***************************************************************************************
    
class test_init():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_init.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_init.teardown_class() -------')

    def setup(self):
        pass
    
    def teardown(self):
        pass
        
    def test_init_from_python(self):
        # create the list of Fisher matrices:
        fisher_list = []
        matrix = np.identity(3)
        for j in xrange(3):
            for i in xrange(3):
                matrix[i,i] = i+j+1
            param_names = [ 'pa'+str(i+j) for i in xrange(3) ] 
            param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
            fiducial = [ float(i) for i in xrange(3) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
        # quality check:
        assert fisher_list_test.fisher_name_list == [ fish.name for fish in fisher_list_test.fisher_list ]
        assert fisher_list_test.fisher_name_list == [ fish.name for fish in fisher_list ]
    
    def test_init_from_python_with_add(self):
        # create the list of Fisher matrices:
        fisher_list = []
        matrix = np.identity(3)
        for j in xrange(3):
            for i in xrange(3):
                matrix[i,i] = i+j+1
            param_names = [ 'pa'+str(i+j) for i in xrange(3) ] 
            param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
            fiducial = [ float(i) for i in xrange(3) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        fisher_list_test = fpa.CosmicFish_FisherAnalysis()
        fisher_list_test.add_fisher_matrix( fisher_list )
        # quality check:
        assert fisher_list_test.fisher_name_list == [ fish.name for fish in fisher_list_test.fisher_list ]
        assert fisher_list_test.fisher_name_list == [ fish.name for fish in fisher_list ]
        
    def test_init_with_add_raises(self):
        # init the list:
        fisher_list_test = fpa.CosmicFish_FisherAnalysis()
        # check first raise:
        assert_raises( ValueError, fisher_list_test.add_fisher_matrix, [0.0] )
        # check second raise:
        matrix = np.identity(3)
        for i in xrange(3):
            matrix[i,i] = i+1
        fisher = fm.fisher_matrix( fisher_matrix=matrix )
        # add a fisher without a name:
        fisher_list_test.add_fisher_matrix( fisher )
        # add another fisher with the same (empty) name:
        assert_raises( ValueError, fisher_list_test.add_fisher_matrix, fisher )
        
    def test_init_from_path_search_false(self):
        # create the list of Fisher matrices:
        fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_path=test_input, search_fisher_guess=False)
        # quality check:
        assert fisher_list_test.fisher_name_list == [ fish.name for fish in fisher_list_test.fisher_list ]
        
    def test_init_from_path_search_true(self):
        # create the list of Fisher matrices:
        fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_path=test_input, search_fisher_guess=True)
        # quality check:
        assert fisher_list_test.fisher_name_list == [ fish.name for fish in fisher_list_test.fisher_list ]
            
    def test_init_search(self):
        # create the two list with different methods:
        fisher_list_test_1 = fpa.CosmicFish_FisherAnalysis(fisher_path=test_input, search_fisher_guess=False)
        fisher_list_test_2 = fpa.CosmicFish_FisherAnalysis(fisher_path=test_input, search_fisher_guess=True)
        fisher_list_test_3 = fpa.CosmicFish_FisherAnalysis(fisher_path=[test_input+'/dummy_fisher_matrix.dat', 
                                                                       test_input+'/dummy_fisher_matrix_2.dat',
                                                                       test_input+'/dummy_fisher_matrix_derived.dat'], search_fisher_guess=True)
        # check wether they catch the same fishers:
        for fish_1 in fisher_list_test_1.fisher_list:
            assert fish_1 in fisher_list_test_2.fisher_list
        for fish_2 in fisher_list_test_2.fisher_list:
            assert fish_2 in fisher_list_test_3.fisher_list
        
    def test_init_without_derived(self):
        # create the two list with different methods:
        fisher_list_test_1 = fpa.CosmicFish_FisherAnalysis(fisher_path=test_input, search_fisher_guess=False, with_derived=True)
        fisher_list_test_2 = fpa.CosmicFish_FisherAnalysis(fisher_path=test_input, search_fisher_guess=True , with_derived=False)
        fisher_list_test_3 = fpa.CosmicFish_FisherAnalysis(fisher_path=[test_input+'/dummy_fisher_matrix.dat', 
                                                                       test_input+'/dummy_fisher_matrix_2.dat', ], search_fisher_guess=True, with_derived=True)
        # check wether they catch the same fishers:
        for fish_1 in fisher_list_test_1.fisher_list:
            assert fish_1 not in fisher_list_test_2.fisher_list
        for fish_2 in fisher_list_test_2.fisher_list:
            assert fish_2 in fisher_list_test_3.fisher_list        
        
# ***************************************************************************************
    
class test_get_fisher_matrix():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_get_fisher_matrix.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_get_fisher_matrix.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        matrix = np.identity(3)
        for j in xrange(3):
            for i in xrange(3):
                matrix[i,i] = i+j+1
            param_names = [ 'pa'+str(i+j) for i in xrange(3) ] 
            param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
            fiducial = [ float(i) for i in xrange(3) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test 
        
    def test_get_fisher_matrix(self):
        
        fisher_list_get = self.fisher_list_test.get_fisher_matrix('fisher1')
        for fish in fisher_list_get:
            assert isinstance( fish, fm.fisher_matrix )
            assert fish.name=='fisher1'
            
# ***************************************************************************************
    
class test_delete_fisher_matrix():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_delete_fisher_matrix.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_delete_fisher_matrix.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        matrix = np.identity(3)
        for j in xrange(3):
            for i in xrange(3):
                matrix[i,i] = i+j+1
            param_names = [ 'pa'+str(i+j) for i in xrange(3) ] 
            param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
            fiducial = [ float(i) for i in xrange(3) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test 
        
    def test_delete_fisher_matrix(self):
        
        fisher_to_delete = self.fisher_list_test.get_fisher_matrix('fisher1')[0]
        self.fisher_list_test.delete_fisher_matrix('fisher1')
        assert fisher_to_delete not in self.fisher_list_test.fisher_list

# ***************************************************************************************
    
class test_get_parameter_list():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_get_parameter_list.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_get_parameter_list.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        matrix = np.identity(3)
        for j in xrange(3):
            for i in xrange(3):
                matrix[i,i] = i+j+1
            param_names = [ 'pa'+str(i+j) for i in xrange(3) ] 
            param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
            fiducial = [ float(i) for i in xrange(3) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test 
        
    def test_get_parameter_list(self):
        
        assert self.fisher_list_test.get_parameter_list() == ['pa0', 'pa1', 'pa2', 'pa3', 'pa4']
        assert self.fisher_list_test.get_parameter_list('fisher1') == ['pa0', 'pa1', 'pa2']
        assert self.fisher_list_test.get_parameter_list(['fisher1','fisher2']) == ['pa0', 'pa1', 'pa2', 'pa3']

# ***************************************************************************************
    
class test_get_parameter_latex_names():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_get_parameter_latex_names.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_get_parameter_latex_names.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        matrix = np.identity(3)
        for j in xrange(3):
            for i in xrange(3):
                matrix[i,i] = i+j+1
            param_names = [ 'pa'+str(i+j) for i in xrange(3) ] 
            param_names_latex = [ 'm'+str(i+j) for i in xrange(3) ] 
            fiducial = [ float(i) for i in xrange(3) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test 
        
    def test_get_parameter_latex_names(self):
        
        assert self.fisher_list_test.get_parameter_latex_names() == {'pa0': 'm0', 'pa1': 'm1', 'pa2': 'm2', 'pa3': 'm3', 'pa4': 'm4'}
        assert self.fisher_list_test.get_parameter_latex_names('fisher1') == {'pa0': 'm0', 'pa1': 'm1', 'pa2': 'm2'}
        assert self.fisher_list_test.get_parameter_latex_names(['fisher1','fisher2']) == {'pa0': 'm0', 'pa1': 'm1', 'pa2': 'm2', 'pa3': 'm3'}     

# ***************************************************************************************
    
class test_reshuffle():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_reshuffle.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_reshuffle.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        matrix = np.identity(3)
        for j in xrange(3):
            for i in xrange(3):
                matrix[i,i] = i+j+1
            param_names = [ 'pa'+str(i+j) for i in xrange(3) ] 
            param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
            fiducial = [ float(i) for i in xrange(3) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test 
        
    def test_reshuffle_1(self):
        temp = self.fisher_list_test.reshuffle(params=['pa1','pa2'])
        assert temp.get_parameter_list() == ['pa1', 'pa2']
        temp = self.fisher_list_test.reshuffle(params=['pa1','pa2'], names=['fisher1'])
        assert temp.get_parameter_list() == ['pa1', 'pa2']
        assert temp.get_fisher_name_list() == ['fisher1_reshuffled']
        
# ***************************************************************************************
    
class test_marginalise():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_marginalise.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_marginalise.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        matrix = np.identity(3)
        for j in xrange(3):
            for i in xrange(3):
                matrix[i,i] = i+j+1
            param_names = [ 'pa'+str(i+j) for i in xrange(3) ] 
            param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
            fiducial = [ float(i) for i in xrange(3) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test 
        
    def test_reshuffle_1(self):
        temp = self.fisher_list_test.marginalise(params=['pa1','pa2'])
        assert temp.get_parameter_list() == ['pa1', 'pa2']
        temp = self.fisher_list_test.marginalise(params=['pa1','pa2'], names=['fisher1'])
        assert temp.get_parameter_list() == ['pa1', 'pa2']
        assert temp.get_fisher_name_list() == ['fisher1_marginal']
        
# ***************************************************************************************
    
class test_plot_range():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_plot_range.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_plot_range.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        matrix = np.identity(3)
        for j in xrange(3):
            for i in xrange(3):
                matrix[i,i] = i+j+1
            param_names = [ 'pa'+str(i+j) for i in xrange(3) ] 
            param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
            fiducial = [ float(i) for i in xrange(3) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test 
        
    def test_compute_plot_range_nice(self):
        assert self.fisher_list_test.compute_plot_range() == {'pa0': [-1.0, 1.0], 'pa1': [-0.8, 1.8], 'pa2': [-0.6, 2.6], 'pa3': [0.5, 2.5], 'pa4': [1.5, 2.5]}
        
    def test_compute_plot_range_ugly(self):
        dict = self.fisher_list_test.compute_plot_range( nice=False )
        test = [ dict[i] for i in self.fisher_list_test.get_parameter_list() ]
        given = [[-0.99445788321, 0.99445788321], [-0.703187912822, 1.70318791282], [-0.574150526569, 2.57415052657], [0.502771058395, 2.4972289416], [1.55526491448, 2.44473508552]]
        assert np.allclose( test,given )
    
    def test_compute_plot_range_parameters(self):
        assert self.fisher_list_test.compute_plot_range('pa1') == {'pa1': [-0.8, 1.8]}
        assert self.fisher_list_test.compute_plot_range(['pa1','pa2']) == {'pa1': [-0.8, 1.8], 'pa2': [-0.6, 2.6]}
        
    def test_compute_plot_range_names(self):
        assert self.fisher_list_test.compute_plot_range('pa1') == {'pa1': [-0.8, 1.8]}
        assert self.fisher_list_test.compute_plot_range('pa1', names='fisher1') == {'pa1': [0.2, 1.8]}
        assert self.fisher_list_test.compute_plot_range('pa1', names=['fisher2']) == {'pa1': [-0.8, 0.8]}
        assert self.fisher_list_test.compute_plot_range('pa1', names=['fisher3']) == {}
        
# ***************************************************************************************
    
class test_compute_gaussian():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_compute_gaussian.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_compute_gaussian.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        matrix = np.identity(3)
        for j in xrange(3):
            for i in xrange(3):
                matrix[i,i] = i+j+1
            param_names = [ 'pa'+str(i+j) for i in xrange(3) ] 
            param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
            fiducial = [ float(i) for i in xrange(3) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test 
        
    def test_compute_gaussian_all(self):
        
        results_1 = self.fisher_list_test.compute_gaussian(num_points=10)
        
        assert results_1.keys() == ['pa0', 'pa1', 'pa2', 'pa3', 'pa4']
        assert results_1.values()[0].keys() == ['fisher2', 'fisher3', 'fisher1']
        
    def test_compute_gaussian_one_param(self):
        
        results_1 = self.fisher_list_test.compute_gaussian( params='pa1', num_points=10)
        
        assert results_1.keys() == ['pa1']
        assert results_1.values()[0].keys() == ['fisher2', 'fisher3', 'fisher1']
        
    def test_compute_gaussian_one_fisher(self):
        
        results_1 = self.fisher_list_test.compute_gaussian( names='fisher1', num_points=10)
        
        assert results_1.keys() == ['pa0', 'pa1', 'pa2']
        assert results_1.values()[0].keys() == ['fisher1']

    def test_compute_gaussian_without_nice(self):
        
        results_1 = self.fisher_list_test.compute_gaussian( names='fisher1', num_points=10, nice_bounds=False )
        
    def test_compute_gaussian_with_normalized(self):
        
        results_1 = self.fisher_list_test.compute_gaussian( names='fisher1', num_points=10, normalized=True )
        
# ***************************************************************************************

class test_compute_compute_ellipse():

    @classmethod
    def setup_class(cls):
        print color_print.header(__name__+': test_compute_compute_ellipse.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold(__name__+': test_compute_compute_ellipse.teardown_class() -------')

    def setup(self):
        # create a list of fisher matrices:
        fisher_list = []
        matrix = np.identity(3)
        for j in xrange(3):
            for i in xrange(3):
                matrix[i,i] = i+j+1
            param_names = [ 'pa'+str(i+j) for i in xrange(3) ] 
            param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
            fiducial = [ float(i) for i in xrange(3) ]
            fisher = fm.fisher_matrix( fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
            fisher.name = 'fisher'+str(j+1)
            fisher_list.append( fisher )
        self.fisher_list_test = fpa.CosmicFish_FisherAnalysis(fisher_list=fisher_list)
    
    def teardown(self):
        del self.fisher_list_test 
        
    def test_compute_ellipse_all(self):
        
        results_1 = self.fisher_list_test.compute_ellipse(num_points=10)
        
        assert results_1.keys() == ['pa0', 'pa1', 'pa2', 'pa3', 'pa4']
        assert results_1.values()[0].keys() == ['pa0', 'pa1', 'pa2', 'pa3', 'pa4']
        assert results_1.values()[0].values()[0].keys() == ['fisher2', 'fisher3', 'fisher1']
        
    def test_compute_ellipse_one_param_multiple_param(self):
        
        results_1 = self.fisher_list_test.compute_ellipse( params1='pa1', num_points=10)
        
        assert results_1.keys() == ['pa1']
        assert results_1.values()[0].keys() == ['pa0', 'pa1', 'pa2', 'pa3', 'pa4']
        assert results_1.values()[0].values()[0].keys() == ['fisher2', 'fisher3', 'fisher1']
    
    def test_compute_ellipse_one_param_one_param(self):
        
        results_1 = self.fisher_list_test.compute_ellipse( params1='pa1',params2='pa2', num_points=10)
        
        assert results_1.keys() == ['pa1']
        assert results_1.values()[0].keys() == ['pa2']
        assert results_1.values()[0].values()[0].keys() == ['fisher2', 'fisher3', 'fisher1']
        
    def test_compute_ellipse_one_fisher(self):
        
        results_1 = self.fisher_list_test.compute_ellipse( names='fisher1', num_points=10)
        
        assert results_1.keys() == ['pa0', 'pa1', 'pa2']
        assert results_1.values()[0].keys() == ['pa0', 'pa1', 'pa2']
        assert results_1.values()[0].values()[0].keys() == ['fisher1']

        
# ***************************************************************************************
        