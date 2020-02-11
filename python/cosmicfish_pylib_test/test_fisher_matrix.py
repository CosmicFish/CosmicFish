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
from dircache import cache
# define path to the executable and add all the relevant folders to the path where python looks for files.
here = os.path.dirname(os.path.abspath(__file__))
cosmicfish_pylib_path = here+'/..'
test_input  = here+'/test_input'
test_output = here+'/test_output'
sys.path.insert(0, os.path.normpath(cosmicfish_pylib_path))

import numpy as np
import cosmicfish_pylib.fisher_matrix as fm

from nose.tools import with_setup
from nose.tools import assert_raises

import cosmicfish_pylib.colors as fc
color_print = fc.bash_colors()

# ***************************************************************************************
    
class test_getters():

    @classmethod
    def setup_class(cls):
        print(color_print.header( __name__+': test_getters.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold( __name__+': test_getters.teardown_class() -------'))

    def setup(self):
        # generate the Fisher matrix. In this case a simple diagonal matrix.
        matrix = np.identity(3)
        for i in range(3):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in range(3) ] 
        fiducial = [ float(i) for i in range(3) ]
        # initialize the Fisher type:
        self.fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )

    # test against known output:
    def test_get_fisher_matrix(self):
        #  1- fisher matrix:
        assert np.allclose( self.fisher_1.get_fisher_matrix(), np.array([[1.0, 0.0, 0.0],[0.0,2.0,0.0],[0.0,0.0,3.0]]) )
    def test_get_fisher_eigenvalues(self):
        #  2- eigenvalues:
        assert np.allclose( self.fisher_1.get_fisher_eigenvalues(), np.array([1.0,2.0,3.0]))
    def test_get_fisher_eigenvectors(self):
        #  3- eigenvectors:
        assert np.allclose( self.fisher_1.get_fisher_eigenvectors(), np.array([[1.0, 0.0, 0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]) )
    def test_get_fisher_inverse(self):
        #  4- inverse:
        assert np.allclose( self.fisher_1.get_fisher_inverse(), np.array([[1.0, 0.0, 0.0],[0.0,0.5,0.0],[0.0,0.0,1.0/3.0]]))
    def test_get_param_names(self):
        #  5- param names:
        assert self.fisher_1.get_param_names()==['p1','p2','p3']
    def test_get_param_names_latex(self):
        #  6- latex param names:
        assert self.fisher_1.get_param_names_latex()==['m0','m1','m2']
    def test_get_param_fiducial(self):
        #  7- fiducial parameters:
        assert np.allclose( self.fisher_1.get_param_fiducial(), np.array([0.0,1.0,2.0]))
    
    # test against class members:
    def test_get_fisher_matrix_self(self):
        #  1- fisher matrix:
        assert np.allclose( self.fisher_1.get_fisher_matrix(), self.fisher_1.fisher_matrix )
    def test_get_fisher_eigenvalues_self(self):
        #  2- eigenvalues:
        assert np.allclose( self.fisher_1.get_fisher_eigenvalues(), self.fisher_1.fisher_eigenvalues )
    def test_get_fisher_eigenvectors_self(self):
        #  3- eigenvectors:
        assert np.allclose( self.fisher_1.get_fisher_eigenvectors(), self.fisher_1.fisher_eigenvectors )
    def test_get_fisher_inverse_self(self):
        #  4- inverse:
        assert np.allclose( self.fisher_1.get_fisher_inverse(), self.fisher_1.fisher_matrix_inv )
    def test_get_param_names_self(self):
        #  5- param names:
        assert self.fisher_1.get_param_names()==self.fisher_1.param_names
    def test_get_param_names_latex_self(self):
        #  6- latex param names:
        assert self.fisher_1.get_param_names_latex()==self.fisher_1.param_names_latex
    def test_get_param_fiducial_self(self):
        #  7- fiducial parameters:
        assert np.allclose( self.fisher_1.get_param_fiducial(), self.fisher_1.param_fiducial )

# ***************************************************************************************
    
class test_advanced_getters():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_advanced_getters.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_advanced_getters.teardown_class() -------'))

    def setup(self):
        # generate the Fisher matrix. In this case a simple diagonal matrix.
        matrix = np.identity(3)
        for i in range(3):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in range(3) ] 
        fiducial = [ float(i) for i in range(3) ]
        # initialize the Fisher type:
        self.fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )

    # test advanced getters:
    def test_get_param_name(self):
        assert [ self.fisher_1.get_param_name(1),self.fisher_1.get_param_name(2),self.fisher_1.get_param_name(3)] == ['p1','p2','p3']
    
    def test_get_param_name_vector(self):
        assert self.fisher_1.get_param_name([1,2,3]) == ['p1','p2','p3']
                
    def test_get_param_ind(self):
        assert [ self.fisher_1.get_param_index('p1'),self.fisher_1.get_param_index('p2'),self.fisher_1.get_param_index('p3')] == [0,1,2]
    
    def test_get_param_ind_vector(self):
        assert self.fisher_1.get_param_index( ['p1','p2','p3'] ) == [0,1,2]
    
    def test_get_param_number(self):
        assert [ self.fisher_1.get_param_number('p1'),self.fisher_1.get_param_number('p2'),self.fisher_1.get_param_number('p3')] == [1,2,3]

    def test_get_param_number_vector(self):
        assert self.fisher_1.get_param_number(['p1','p2','p3']) == [1,2,3]

    def test_get_param_name_latex(self):
        assert [self.fisher_1.get_param_name_latex('p1'),self.fisher_1.get_param_name_latex('p2'),self.fisher_1.get_param_name_latex('p3')] == ['m0','m1','m2']

    def test_get_param_name_latex_vector(self):
        assert self.fisher_1.get_param_name_latex(['p1','p2','p3']) == ['m0','m1','m2']

    def test_get_param_fiducial(self):
        assert [ self.fisher_1.get_fiducial('p1'), self.fisher_1.get_fiducial('p2'), self.fisher_1.get_fiducial('p3')] == [0.0,1.0,2.0]
    
    def test_get_param_fiducial_vector(self):
        assert self.fisher_1.get_fiducial(['p1','p2','p3']) == [0.0,1.0,2.0]
        
# ***************************************************************************************
    
class test_fisher_init():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_fisher_init.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_fisher_init.teardown_class() -------'))

    def setup(self):
        pass
    
    # test the class init:
    def test_init_fisher_invalid(self):
        assert_raises( ValueError, fm.fisher_matrix ) 
    
    def test_init_fisher_nonsymm_1(self):
        fisher_asymm = []
        for i in range(2):
            fisher_asymm.append([])
            for j in range(3):
                fisher_asymm[i].append(0.0)
        assert_raises( ValueError, fm.fisher_matrix, fisher_matrix=fisher_asymm )
    
    def test_init_fisher_nonsymm_2(self):
        matrix = np.identity(10)
        matrix[1,2] = +1.0
        matrix[2,1] = -1.0
        assert_raises( ValueError, fm.fisher_matrix, fisher_matrix=matrix )
    
    def test_init_from_file(self):
        fisher_1 = fm.fisher_matrix( file_name=test_input+'/dummy_fisher_matrix.dat' )
        assert np.isclose( fisher_1.fisher_cutoff , 0.100488387718 )
        assert fisher_1.path  == test_input+'/dummy_fisher_matrix.dat'
        assert fisher_1.name  == 'dummy_fisher_matrix'
        assert fisher_1.indir == test_input
        assert fisher_1.num_params == 28
        assert fisher_1.param_names == ['omegabh2', 'omegach2', 'omeganuh2', 'h', 'yhe', 'logA', 'ns', 'nrun', 'nt', 'r', 'tau', 'Bias_W_1', 'Bias_W_2', 'Bias_W_3', 'Bias_W_4', 'Bias_W_5', 'Bias_W_6', 'Bias_W_7', 'Bias_W_8', 'Bias_W_9', 'Bias_W_10', 'Bias_W_11', 'Bias_W_12', 'Bias_W_13', 'Bias_W_14', 'alpha_SN', 'beta_SN', 'M0_SN']
        assert fisher_1.param_names_latex == ['\\Omega_b h^2', '\\Omega_c h^2', '\\Omega_\\nu h^2', 'h', 'Y_{He}', '{\\rm{ln}}(10^{10} A_s)', 'n_s', 'n_{\\rm run}', 'n_t', 'r', '\\tau', 'b_1', 'b_2', 'b_3', 'b_4', 'b_5', 'b_6', 'b_7', 'b_8', 'b_9', 'b_10', 'b_11', 'b_12', 'b_13', 'b_14', '\\alpha_{\\rm SN}', '\\beta_{\\rm SN}', 'M_0^{\\rm SN}']
      
    def test_init_from_python(self):
        matrix = np.identity(10)
        for i in range(10):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in range(10) ] 
        fiducial = [ float(i) for i in range(10) ]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        assert fisher_1.fisher_cutoff == 1e-09
        assert np.allclose( fisher_1.fisher_matrix, matrix )
        assert fisher_1.path  == ''
        assert fisher_1.name  == ''
        assert fisher_1.indir == ''
        assert fisher_1.num_params == 10
        assert np.allclose( fisher_1.fisher_eigenvalues, [ 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0 ] )
        assert fisher_1.param_names == ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10']
        assert fisher_1.param_names_latex == ['m0', 'm1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9']
        assert np.allclose( fisher_1.param_fiducial, [ 0., 1., 2., 3., 4., 5., 6., 7., 8., 9.] )
        assert fisher_1.param_names_dict == {'p2': 2, 1: 'p1', 2: 'p2', 'p1': 1, 4: 'p4', 5: 'p5', 'p4': 4, 'p3': 3, 8: 'p8', 9: 'p9', 'p8': 8, 'p5': 5, 'p10': 10, 3: 'p3', 7: 'p7', 'p6': 6, 6: 'p6', 'p9': 9, 10: 'p10', 'p7': 7}
    
        
    def test_init_0D_fisher(self):
        matrix = 10
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix )
        assert fisher_1.fisher_cutoff == 1e-09
        assert np.allclose( fisher_1.fisher_matrix, matrix )
        assert fisher_1.path  == ''
        assert fisher_1.name  == ''
        assert fisher_1.indir == ''
        assert fisher_1.num_params == 1
        assert np.allclose( fisher_1.fisher_eigenvalues, [ 10.0 ] )
        assert np.allclose( fisher_1.fisher_eigenvectors, [[ 1.]] )
        assert np.allclose( fisher_1.fisher_matrix_inv, [[ 0.1]] )
        assert fisher_1.param_names == ['p1']
        assert fisher_1.param_names_latex == fisher_1.param_names
        assert np.allclose( fisher_1.param_fiducial, [ 0.0 ] )
        assert fisher_1.param_names_dict == {1: 'p1', 'p1': 1}

    def test_init_1D_fisher(self):
        matrix = [10]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix )
        # after calling init all the objects of the class have to be initialized properly:
        assert fisher_1.fisher_cutoff == 1e-09
        assert np.allclose( fisher_1.fisher_matrix, matrix )
        assert fisher_1.path  == ''
        assert fisher_1.name  == ''
        assert fisher_1.indir == ''
        assert fisher_1.num_params == 1
        assert np.allclose( fisher_1.fisher_eigenvalues, [ 10.0 ] )
        assert np.allclose( fisher_1.fisher_eigenvectors, [[ 1.]] )
        assert np.allclose( fisher_1.fisher_matrix_inv, [[ 0.1]] )
        assert fisher_1.param_names == ['p1']
        assert fisher_1.param_names_latex == fisher_1.param_names
        assert np.allclose( fisher_1.param_fiducial, [ 0.0 ] )
        assert fisher_1.param_names_dict == {1: 'p1', 'p1': 1}
        
    def test_init_from_python_only_fisher(self):
        matrix = np.identity(10)
        for i in range(10):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in range(10) ] 
        fiducial = [ float(i) for i in range(10) ]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix )
        # after calling init all the objects of the class have to be initialized properly:
        assert fisher_1.fisher_cutoff == 1e-09
        assert np.allclose( fisher_1.fisher_matrix, matrix )
        assert fisher_1.path  == ''
        assert fisher_1.name  == ''
        assert fisher_1.indir == ''
        assert fisher_1.num_params == 10
        assert np.allclose( fisher_1.fisher_eigenvalues, [ 1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0 ] )
        assert fisher_1.param_names == ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10']
        assert fisher_1.param_names_latex == fisher_1.param_names
        assert np.allclose( fisher_1.param_fiducial, [ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.] )
        assert fisher_1.param_names_dict == {'p2': 2, 1: 'p1', 2: 'p2', 'p1': 1, 4: 'p4', 5: 'p5', 'p4': 4, 'p3': 3, 8: 'p8', 9: 'p9', 'p8': 8, 'p5': 5, 'p10': 10, 3: 'p3', 7: 'p7', 'p6': 6, 6: 'p6', 'p9': 9, 10: 'p10', 'p7': 7}

    def test_init_from_python_invalid_param_names(self):
        matrix = np.identity(10)
        for i in range(10):
            matrix[i,i] = i+1
        param_names = [ 'q'+str(i) for i in range(9) ]     
        param_names_latex = [ 'm'+str(i) for i in range(10) ] 
        fiducial = [ float(i) for i in range(10) ]
        assert_raises( ValueError, fm.fisher_matrix, fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
    
    def test_init_from_python_invalid_param_names_latex(self):
        matrix = np.identity(10)
        for i in range(10):
            matrix[i,i] = i+1
        param_names = [ 'q'+str(i) for i in range(10) ]     
        param_names_latex = [ 'm'+str(i) for i in range(9) ] 
        fiducial = [ float(i) for i in range(10) ]
        assert_raises( ValueError, fm.fisher_matrix, fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
    
    def test_init_from_python_invalid_fiducial(self):
        matrix = np.identity(10)
        for i in range(10):
            matrix[i,i] = i+1
        param_names = [ 'q'+str(i) for i in range(10) ]     
        param_names_latex = [ 'm'+str(i) for i in range(10) ] 
        fiducial = [ float(i) for i in range(9) ]
        assert_raises( ValueError, fm.fisher_matrix, fisher_matrix=matrix, param_names=param_names, param_names_latex=param_names_latex, fiducial=fiducial )
         
# ***************************************************************************************
    
class test_fisher_load_paramnames_from_file():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_fisher_load_paramnames_from_file.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_fisher_load_paramnames_from_file.teardown_class() -------'))

    def setup(self):
        pass
    
    # test the loading of the parameter names from file:
    def test_load_paramnames_from_file_1(self):
        matrix = np.identity(5)
        for i in range(5):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in range(5) ] 
        fiducial = [ float(i) for i in range(5) ]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix )
        fisher_1.load_paramnames_from_file( file_name=test_input+'/dummy_paramnames_1.paramnames')
        # test if everything is properly initialized:
        assert fisher_1.param_names == ['p1', 'p2', 'p3', 'p4', 'p5']
        assert fisher_1.param_names_latex == ['p1', 'p2', 'p3', 'p4', 'p5']
        assert np.allclose( fisher_1.param_fiducial, [ 0.0,0.0,0.0,0.0,0.0 ] )

    def test_load_paramnames_from_file_2(self):
        matrix = np.identity(5)
        for i in range(5):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in range(5) ] 
        fiducial = [ float(i) for i in range(5) ]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix )
        fisher_1.load_paramnames_from_file( file_name=test_input+'/dummy_paramnames_2.paramnames')
        # test if everything is properly initialized:
        assert fisher_1.param_names == ['p1', 'p2', 'p3', 'p4', 'p5']
        assert fisher_1.param_names_latex == ['p_1', 'p_2', 'p_3', 'p_4', 'p_5']
        assert np.allclose( fisher_1.param_fiducial, [ 0.0,0.0,0.0,0.0,0.0 ] )

    def test_load_paramnames_from_file_3(self):
        matrix = np.identity(5)
        for i in range(5):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in range(5) ] 
        fiducial = [ float(i) for i in range(5) ]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix )
        fisher_1.load_paramnames_from_file( file_name=test_input+'/dummy_paramnames_3.paramnames')
        # test if everything is properly initialized:
        assert fisher_1.param_names == ['p1', 'p2', 'p3', 'p4', 'p5']
        assert fisher_1.param_names_latex == ['p_1', 'p_2', 'p_3', 'p_4', 'p_5']
        assert np.allclose( fisher_1.param_fiducial, [ 1.0,2.0,3.0,4.0,5.0 ] )
        
    def test_load_paramnames_from_file_4(self):
        matrix = np.identity(5)
        for i in range(5):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in range(5) ] 
        fiducial = [ float(i) for i in range(5) ]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix )
        fisher_1.load_paramnames_from_file( file_name=test_input+'/dummy_paramnames_4.paramnames')
        # test if everything is properly initialized:
        assert fisher_1.param_names == ['p1', 'p2', 'p3', 'p4', 'p5']
        assert fisher_1.param_names_latex == ['p1', 'p2', 'p3', 'p4', 'p5']
        assert np.allclose( fisher_1.param_fiducial, [ 1.0,2.0,3.0,4.0,5.0 ] )
    
    def test_load_paramnames_from_file_invalid_num_params(self):
        matrix = np.identity(5)
        for i in range(5):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in range(5) ] 
        fiducial = [ float(i) for i in range(5) ]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix )
        assert_raises( ValueError, fisher_1.load_paramnames_from_file, file_name=test_input+'/dummy_paramnames_5.paramnames' )

# ***************************************************************************************
    
class test_fisher_save_paramnames_to_file():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_fisher_save_paramnames_to_file.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_fisher_save_paramnames_to_file.teardown_class() -------'))

    def setup(self):
        pass
    
    # test the loading of the parameter names from file:
    def test_fisher_save_paramnames_to_file_1(self):
        matrix = np.identity(5)
        for i in range(5):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in range(5) ] 
        fiducial = [ float(i) for i in range(5) ]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        fisher_1.save_paramnames_to_file(file_name=test_output+'/dummy_paramnames_out.paramnames')
        fisher_1.indir = './dont_exist'
        fisher_1.name  = 'dont_exist'
        assert_raises( IOError, fisher_1.save_paramnames_to_file )

# ***************************************************************************************
    
class test_save_to_file():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_save_to_file.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_save_to_file.teardown_class() -------'))

    def setup(self):
        pass
    
    # test the loading of the parameter names from file:
    def test_save_to_file_1(self):
        # load the fisher matrix:
        fisher_1 = fm.fisher_matrix( file_name=test_input+'/dummy_fisher_matrix_2.dat' )
        # save the fisher matrix back to file:
        fisher_1.save_to_file( file_name=test_output+'/dummy_fisher_matrix_out' )        
        # import it as a second fisher:
        fisher_2 = fm.fisher_matrix( file_name=test_output+'/dummy_fisher_matrix_out.dat' )
        # test equality, to do so we have to set a couple of things to be the same...
        fisher_1.path  = fisher_2.path
        fisher_1.name  = fisher_2.name
        fisher_1.indir = fisher_2.indir
        fisher_1.protect_degenerate( cache=False ) # we have to do it twice otherwise eigenvalues are slightly different
        assert fisher_1==fisher_2
                               
# ***************************************************************************************
    
class test_fisher_overload_operations():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_fisher_overload_operations.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_fisher_overload_operations.teardown_class() -------'))

    def setup(self):
        pass
    
    # test the fisher matrix addition:
    def test_add_same_params(self):
        matrix = np.identity(3)
        for i in range(3):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in range(3) ] 
        fiducial = [ float(i) for i in range(3) ]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        fisher_2 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        fisher_3 = fisher_1+fisher_2
        assert np.allclose( 2.0*fisher_1.fisher_matrix, fisher_3.fisher_matrix )
        assert fisher_3.path  == fisher_3.path
        assert fisher_3.name  == "_"
        assert fisher_3.indir == fisher_3.indir
        assert fisher_3.num_params == 3
        assert np.allclose( fisher_3.fisher_eigenvalues, [ 2., 4., 6.] )
        assert fisher_3.param_names == ['p1', 'p2', 'p3']
        assert fisher_3.param_names_latex == fisher_1.param_names_latex
        assert np.allclose( fisher_3.param_fiducial, [ 0., 1., 2.] )
    
    # test the loading of the parameter names from file:
    def test_add_different_fiducial(self):
        matrix = np.identity(3)
        for i in range(3):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in range(3) ] 
        fiducial_1 = [ float(i) for i in range(3) ]
        fiducial_2 = [ float(i+1) for i in range(3) ]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial_1 )
        fisher_2 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial_2 )
        
        assert_raises( ValueError, lambda: fisher_1+fisher_2 )
    
    def test_add_different_params(self):
        matrix = np.identity(2)
        for i in range(2):
            matrix[i,i] = i+1
        param_names_1 = [ 'm'+str(i) for i in range(2) ] 
        param_names_2 = [ 'b'+str(i) for i in range(2) ] 
        fiducial = [ float(i) for i in range(2) ]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names = param_names_1, fiducial=fiducial )
        fisher_2 = fm.fisher_matrix( fisher_matrix=matrix, param_names = param_names_2, fiducial=fiducial )
        fisher_3 = fisher_1+fisher_2
        
        assert np.allclose( fisher_3.fisher_matrix, [[1.0,0.0,0.0,0.0],[0.0,2.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,2.0]] )
        assert fisher_3.num_params == 4
        assert np.allclose( fisher_3.fisher_eigenvalues, [ 1.0,1.0,2.0,2.0] )
        assert fisher_3.param_names == ['m0', 'm1', 'b0', 'b1']
        assert np.allclose( fisher_3.param_fiducial, [ 0., 1., 0., 1.] ) 
        
    def test_equality(self):
        matrix = np.identity(2)
        for i in range(2):
            matrix[i,i] = i+1
        param_names = [ 'm'+str(i) for i in range(2) ] 
        fiducial = [ float(i) for i in range(2) ]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names = param_names, fiducial=fiducial )
        
        assert fisher_1 != [0.0]
    
    def test_equality_except(self):
        matrix = np.identity(2)
        for i in range(2):
            matrix[i,i] = i+1
        param_names = [ 'm'+str(i) for i in range(2) ] 
        fiducial = [ float(i) for i in range(2) ]
        fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names = param_names, fiducial=fiducial )
        
        matrix = np.identity(3)
        for i in range(3):
            matrix[i,i] = i+1
        param_names = [ 'm'+str(i) for i in range(3) ] 
        fiducial = [ float(i) for i in range(3) ]
        fisher_2 = fm.fisher_matrix( fisher_matrix=matrix, param_names = param_names, fiducial=fiducial )
        
        assert fisher_1 != fisher_2
    
# ***************************************************************************************
    
class test_fisher_std_operations():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_fisher_load_paramnames_from_file.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_fisher_load_paramnames_from_file.teardown_class() -------'))

    def setup(self):
        matrix = np.identity(3)
        param_names_latex = [ 'm'+str(i) for i in range(3) ] 
        fiducial = [ float(i) for i in range(3) ]
        self.fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        pass
    
    def test_determinant(self):
        assert np.isclose( self.fisher_1.determinant(), 1.0 )
    
    def test_protect_degenerate(self):
        self.fisher_1.protect_degenerate(cache=False)
    
    def test_confidence_bound_error(self):
        assert_raises( ValueError, lambda: self.fisher_1.get_confidence_bounds( confidence_level=-1.0 ) )
        
    def test_confidence_bound_cache(self):
        assert np.allclose( self.fisher_1.get_confidence_bounds( ), self.fisher_1.get_confidence_bounds( cache=True ) )

# ***************************************************************************************
    
class test_fisher_setters():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_fisher_setters.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_fisher_setters.teardown_class() -------'))

    def setup(self):
        matrix = np.identity(3)
        param_names_latex = [ 'm'+str(i) for i in range(3) ] 
        fiducial = [ float(i) for i in range(3) ]
        self.fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        pass
    
    def test_set_fisher_matrix_0D(self):
        matrix = 10.0
        self.fisher_1.set_fisher_matrix( fisher_matrix=matrix )
        
    def test_set_fisher_matrix_1D(self):
        matrix = [10.0]
        self.fisher_1.set_fisher_matrix( fisher_matrix=matrix )
        
    def test_set_fisher_matrix_2D(self):
        matrix = 2.0*np.identity(3)
        self.fisher_1.set_fisher_matrix( fisher_matrix=matrix )
        
    def test_set_fisher_matrix_nonsymm(self):
        matrix = np.identity(10)
        matrix[1,2] = +1.0
        matrix[2,1] = -1.0
        assert_raises( ValueError, self.fisher_1.set_fisher_matrix, fisher_matrix=matrix )
        
    def test_set_param_names_illegal(self):
        list = [ 'm'+str(i) for i in range(4) ] 
        assert_raises( ValueError, self.fisher_1.set_param_names, list )
    
    def test_set_param_names(self):
        list = [ 'new'+str(i) for i in range(3) ] 
        self.fisher_1.set_param_names(list)
        assert self.fisher_1.param_names == ['new0', 'new1', 'new2']
        assert self.fisher_1.param_names == self.fisher_1.param_names_latex
    
    def test_set_param_names_latex_illegal(self):
        list = [ 'm'+str(i) for i in range(4) ] 
        assert_raises( ValueError, self.fisher_1.set_param_names_latex, list )
    
    def test_set_param_names_latex(self):
        list = [ 'new'+str(i) for i in range(3) ] 
        self.fisher_1.set_param_names_latex(list)
        assert self.fisher_1.param_names_latex == ['new0', 'new1', 'new2']
    
    def test_set_fiducial_illegal(self):
        fiducial = [ float(i) for i in range(4) ]
        assert_raises( ValueError, self.fisher_1.set_fiducial, fiducial )
    
    def test_set_fiducial(self):
        fiducial = [ float(i+2) for i in range(3) ]
        self.fisher_1.set_fiducial(fiducial)
        assert np.allclose( self.fisher_1.param_fiducial, [2.0,3.0,4.0] )
        
# ***************************************************************************************
        