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
import cosmicfish_pylib.fisher_derived as fd

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
        # Jacobian matrix:
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        # parameter names:
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names_latex = [ 'md'+str(i) for i in range(2) ]
        # fiducial:
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(2) ]
        # initialize the Fisher type:
        self.fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        # initialize the derived matrix:
        self.fisher_der = fd.fisher_derived( derived_matrix=matrix_derived, 
                                             param_names_latex=param_names_latex, 
                                             derived_param_names_latex=derived_param_names_latex, 
                                             fiducial=fiducial, fiducial_derived=fiducial_derived)
        
    # test against known output:
    def test_get_derived_matrix(self):
        #  1- Jacobian matrix:
        assert np.allclose( self.fisher_der.get_derived_matrix(), self.fisher_der.derived_matrix )
    def test_get_param_names(self):
        #  2- parameter names:
        assert self.fisher_der.get_param_names() == self.fisher_der.param_names
    def test_get_param_names_latex(self):
        #  3- latex parameter names:
        assert self.fisher_der.get_param_names_latex() == self.fisher_der.param_names_latex
    def test_get_param_fiducial(self):
        #  4- parameter fiducial:
        assert np.allclose( self.fisher_der.get_param_fiducial(), self.fisher_der.param_fiducial )
    def test_get_derived_param_names(self):
        #  5- derived parameter names:
        assert self.fisher_der.get_derived_param_names() == self.fisher_der.derived_param_names 
    def test_get_derived_param_names_latex(self):
        #  6- latex derived parameters names:
        assert self.fisher_der.get_derived_param_names_latex() == self.fisher_der.derived_param_names_latex 
    def test_get_derived_param_fiducial(self):
        #  7- derived parameters fiducial:
        assert np.allclose( self.fisher_der.get_derived_param_fiducial(), self.fisher_der.derived_param_fiducial ) 

# ***************************************************************************************
    
class test_fisher_derived_init():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_fisher_derived_init.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_fisher_derived_init.teardown_class() -------'))

    def setup(self):
        pass
    
    # test the class init:
    def test_init_fisher_derived_invalid(self):
        assert_raises( ValueError, fd.fisher_derived )
    
    def test_init_from_file(self):
        fisher_1 = fd.fisher_derived( file_name=test_input+'/dummy_fisher_matrix_derived.dat' )
        assert fisher_1.path == test_input+'/dummy_fisher_matrix_derived.dat'
        assert fisher_1.name == 'dummy_fisher_matrix_derived'
        assert fisher_1.indir == test_input
        assert fisher_1.num_params == 28
        assert fisher_1.num_derived == 9
        assert fisher_1.param_names == ['omegabh2', 'omegach2', 'omeganuh2', 'h', 'yhe', 'logA', 'ns', 'nrun', 'nt', 'r', 'tau', 'Bias_W_1', 'Bias_W_2', 'Bias_W_3', 'Bias_W_4', 'Bias_W_5', 'Bias_W_6', 'Bias_W_7', 'Bias_W_8', 'Bias_W_9', 'Bias_W_10', 'Bias_W_11', 'Bias_W_12', 'Bias_W_13', 'Bias_W_14', 'alpha_SN', 'beta_SN', 'M0_SN']
        assert fisher_1.param_names_latex == ['\\Omega_b h^2', '\\Omega_c h^2', '\\Omega_\\nu h^2', 'h', 'Y_{He}', '{\\rm{ln}}(10^{10} A_s)', 'n_s', 'n_{\\rm run}', 'n_t', 'r', '\\tau', 'b_1', 'b_2', 'b_3', 'b_4', 'b_5', 'b_6', 'b_7', 'b_8', 'b_9', 'b_10', 'b_11', 'b_12', 'b_13', 'b_14', '\\alpha_{\\rm SN}', '\\beta_{\\rm SN}', 'M_0^{\\rm SN}']
        assert fisher_1.derived_param_names == ['omegab', 'omegac', 'omeganu', 'omegav', 'omegak', 'omegam', 'theta', 'mnu', 'z_re']
        assert fisher_1.derived_param_names_latex == ['\\Omega_b', '\\Omega_c', '\\Omega_{\\nu}', '\\Omega_{\\Lambda}', '\\Omega_{K}', '\\Omega_{m}', '100\\theta_{MC}', '\\Sigma m_\\nu', 'z_{\\rm re}']

    def test_init_0D_fisher_derived(self):
        matrix = 10
        fisher_1 = fd.fisher_derived( derived_matrix=matrix )
        assert np.allclose( fisher_1.derived_matrix, matrix )           
        assert fisher_1.path == ''
        assert fisher_1.name == ''
        assert fisher_1.indir == ''
        assert fisher_1.num_params == 1
        assert fisher_1.num_derived == 1
        assert fisher_1.param_names == ['p1']
        assert fisher_1.param_names_latex == fisher_1.param_names
        assert np.allclose( fisher_1.param_fiducial, [0.0] )
        assert fisher_1.derived_param_names == ['p2']
        assert fisher_1.derived_param_names_latex == fisher_1.derived_param_names
        assert np.allclose( fisher_1.derived_param_fiducial, [0.0] )
        
    def test_init_1D_fisher_derived(self):
        matrix = [10]
        fisher_1 = fd.fisher_derived( derived_matrix=matrix )
        assert np.allclose( fisher_1.derived_matrix, matrix )           
        assert fisher_1.path == ''
        assert fisher_1.name == ''
        assert fisher_1.indir == ''
        assert fisher_1.num_params == 1
        assert fisher_1.num_derived == 1
        assert fisher_1.param_names == ['p1']
        assert fisher_1.param_names_latex == fisher_1.param_names
        assert np.allclose( fisher_1.param_fiducial, [0.0] )
        assert fisher_1.derived_param_names == ['p2']
        assert fisher_1.derived_param_names_latex == fisher_1.derived_param_names
        assert np.allclose( fisher_1.derived_param_fiducial, [0.0] )
        
    def test_init_from_python_only_fisher_derived(self):
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        # parameter names:
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names_latex = [ 'md'+str(i) for i in range(2) ]
        # fiducial:
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(2) ]
        # initialize the derived matrix:
        fisher_1 = fd.fisher_derived( derived_matrix=matrix_derived )
        # after calling init all the objects of the class have to be initialized properly:
        assert np.allclose( fisher_1.derived_matrix, matrix_derived )           
        assert fisher_1.path == ''
        assert fisher_1.name == ''
        assert fisher_1.indir == ''
        assert fisher_1.num_params == 3
        assert fisher_1.num_derived == 2
        assert fisher_1.param_names == ['p1','p2','p3']
        assert fisher_1.param_names_latex == fisher_1.param_names
        assert np.allclose( fisher_1.param_fiducial, [0.0,0.0,0.0] )
        assert fisher_1.derived_param_names == ['p4','p5']
        assert fisher_1.derived_param_names_latex == fisher_1.derived_param_names
        assert np.allclose( fisher_1.derived_param_fiducial, [0.0,0.0] )
    
    def test_init_from_python_invalid_param_names(self):
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        # parameter names:
        param_names = [ 'q'+str(i) for i in range(4) ]   
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names = [ 'qd'+str(i) for i in range(2) ]   
        derived_param_names_latex = [ 'md'+str(i) for i in range(2) ]
        # fiducial:
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(2) ]
        # initialize the derived matrix:
        fisher_1 = fd.fisher_derived( derived_matrix=matrix_derived )
        assert_raises( ValueError,  fd.fisher_derived, derived_matrix=matrix_derived,
                                             param_names = param_names,
                                             derived_param_names = derived_param_names,
                                             param_names_latex=param_names_latex, 
                                             derived_param_names_latex=derived_param_names_latex, 
                                             fiducial=fiducial, fiducial_derived=fiducial_derived)
    
    def test_init_from_python_invalid_param_names_latex(self):
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        # parameter names:
        param_names = [ 'q'+str(i) for i in range(3) ]   
        param_names_latex = [ 'm'+str(i) for i in range(4) ]
        derived_param_names = [ 'qd'+str(i) for i in range(2) ]   
        derived_param_names_latex = [ 'md'+str(i) for i in range(2) ]
        # fiducial:
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(2) ]
        # initialize the derived matrix:
        fisher_1 = fd.fisher_derived( derived_matrix=matrix_derived )
        assert_raises( ValueError,  fd.fisher_derived, derived_matrix=matrix_derived,
                                             param_names = param_names,
                                             derived_param_names = derived_param_names,
                                             param_names_latex=param_names_latex, 
                                             derived_param_names_latex=derived_param_names_latex, 
                                             fiducial=fiducial, fiducial_derived=fiducial_derived)
    
    def test_init_from_python_invalid_derived_param_names(self):
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        # parameter names:
        param_names = [ 'q'+str(i) for i in range(3) ]   
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names = [ 'qd'+str(i) for i in range(3) ]   
        derived_param_names_latex = [ 'md'+str(i) for i in range(2) ]
        # fiducial:
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(2) ]
        # initialize the derived matrix:
        fisher_1 = fd.fisher_derived( derived_matrix=matrix_derived )
        assert_raises( ValueError,  fd.fisher_derived, derived_matrix=matrix_derived,
                                             param_names = param_names,
                                             derived_param_names = derived_param_names,
                                             param_names_latex=param_names_latex, 
                                             derived_param_names_latex=derived_param_names_latex, 
                                             fiducial=fiducial, fiducial_derived=fiducial_derived)
    
    def test_init_from_python_invalid_derived_param_names_latex(self):
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        # parameter names:
        param_names = [ 'q'+str(i) for i in range(3) ]   
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names = [ 'qd'+str(i) for i in range(2) ]   
        derived_param_names_latex = [ 'md'+str(i) for i in range(3) ]
        # fiducial:
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(2) ]
        # initialize the derived matrix:
        fisher_1 = fd.fisher_derived( derived_matrix=matrix_derived )
        assert_raises( ValueError,  fd.fisher_derived, derived_matrix=matrix_derived,
                                             param_names = param_names,
                                             derived_param_names = derived_param_names,
                                             param_names_latex=param_names_latex, 
                                             derived_param_names_latex=derived_param_names_latex, 
                                             fiducial=fiducial, fiducial_derived=fiducial_derived) 
    
    def test_init_from_python_invalid_fiducial(self):
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        # parameter names:
        param_names = [ 'q'+str(i) for i in range(3) ]   
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names = [ 'qd'+str(i) for i in range(2) ]   
        derived_param_names_latex = [ 'md'+str(i) for i in range(2) ]
        # fiducial:
        fiducial = [ float(i) for i in range(4) ]
        fiducial_derived = [ float(i) for i in range(2) ]
        # initialize the derived matrix:
        fisher_1 = fd.fisher_derived( derived_matrix=matrix_derived )
        assert_raises( ValueError,  fd.fisher_derived, derived_matrix=matrix_derived,
                                             param_names = param_names,
                                             derived_param_names = derived_param_names,
                                             param_names_latex=param_names_latex, 
                                             derived_param_names_latex=derived_param_names_latex, 
                                             fiducial=fiducial, fiducial_derived=fiducial_derived)
    
    def test_init_from_python_invalid_fiducial_derived(self):
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        # parameter names:
        param_names = [ 'q'+str(i) for i in range(3) ]   
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names = [ 'qd'+str(i) for i in range(2) ]   
        derived_param_names_latex = [ 'md'+str(i) for i in range(2) ]
        # fiducial:
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(3) ]
        # initialize the derived matrix:
        fisher_1 = fd.fisher_derived( derived_matrix=matrix_derived )
        assert_raises( ValueError,  fd.fisher_derived, derived_matrix=matrix_derived,
                                             param_names = param_names,
                                             derived_param_names = derived_param_names,
                                             param_names_latex=param_names_latex, 
                                             derived_param_names_latex=derived_param_names_latex, 
                                             fiducial=fiducial, fiducial_derived=fiducial_derived)  

# ***************************************************************************************
    
class test_fisher_derived_load_paramnames_from_file():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_fisher_derived_load_paramnames_from_file.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_fisher_derived_load_paramnames_from_file.teardown_class() -------'))

    def setup(self):
        pass
    
    # test the loading of the parameter names from file:
    def test_load_paramnames_from_file_1(self):
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        param_names = [ 'q'+str(i) for i in range(3) ]   
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names = [ 'qd'+str(i) for i in range(2) ]   
        derived_param_names_latex = [ 'md'+str(i) for i in range(2) ]
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(2) ]
        fisher_1 = fd.fisher_derived( derived_matrix=matrix_derived )
        fisher_1.load_paramnames_from_file( file_name=test_input+'/dummy_paramnames_1_derived.paramnames')
        # test if everything is properly initialized:
        assert fisher_1.param_names == ['p1', 'p2', 'p3']     
        assert fisher_1.param_names_latex == ['p1', 'p2', 'p3']
        assert np.allclose( fisher_1.param_fiducial, [0.0, 0.0, 0.0] )
        assert fisher_1.derived_param_names == ['p4', 'p5']
        assert fisher_1.derived_param_names_latex == ['p4', 'p5']
        assert np.allclose( fisher_1.derived_param_fiducial, [0.0, 0.0] ) 
    
    # test the loading of the parameter names from file:
    def test_load_paramnames_from_file_2(self):
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        param_names = [ 'q'+str(i) for i in range(3) ]   
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names = [ 'qd'+str(i) for i in range(2) ]   
        derived_param_names_latex = [ 'md'+str(i) for i in range(2) ]
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(2) ]
        fisher_1 = fd.fisher_derived( derived_matrix=matrix_derived )
        fisher_1.load_paramnames_from_file( file_name=test_input+'/dummy_paramnames_2_derived.paramnames')
        # test if everything is properly initialized:
        assert fisher_1.param_names == ['p1', 'p2', 'p3']     
        assert fisher_1.param_names_latex == ['p_1', 'p_2', 'p_3']
        assert np.allclose( fisher_1.param_fiducial, [0.0, 0.0, 0.0] )
        assert fisher_1.derived_param_names == ['p4', 'p5']
        assert fisher_1.derived_param_names_latex == ['p_4', 'p_5']
        assert np.allclose( fisher_1.derived_param_fiducial, [0.0, 0.0] ) 
    
    # test the loading of the parameter names from file:
    def test_load_paramnames_from_file_3(self):
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        param_names = [ 'q'+str(i) for i in range(3) ]   
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names = [ 'qd'+str(i) for i in range(2) ]   
        derived_param_names_latex = [ 'md'+str(i) for i in range(2) ]
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(2) ]
        fisher_1 = fd.fisher_derived( derived_matrix=matrix_derived )
        fisher_1.load_paramnames_from_file( file_name=test_input+'/dummy_paramnames_3_derived.paramnames')
        # test if everything is properly initialized:
        assert fisher_1.param_names == ['p1', 'p2', 'p3']     
        assert fisher_1.param_names_latex == ['p_1', 'p_2', 'p_3']
        assert np.allclose( fisher_1.param_fiducial, [1.0, 2.0, 3.0] )
        assert fisher_1.derived_param_names == ['p4', 'p5']
        assert fisher_1.derived_param_names_latex == ['p_4', 'p_5']
        assert np.allclose( fisher_1.derived_param_fiducial, [4.0, 5.0] ) 
 
    # test the loading of the parameter names from file:
    def test_load_paramnames_from_file_4(self):
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        param_names = [ 'q'+str(i) for i in range(3) ]   
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names = [ 'qd'+str(i) for i in range(2) ]   
        derived_param_names_latex = [ 'md'+str(i) for i in range(2) ]
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(2) ]
        fisher_1 = fd.fisher_derived( derived_matrix=matrix_derived )
        fisher_1.load_paramnames_from_file( file_name=test_input+'/dummy_paramnames_4_derived.paramnames')
        # test if everything is properly initialized:
        assert fisher_1.param_names == ['p1', 'p2', 'p3']     
        assert fisher_1.param_names_latex == ['p1', 'p2', 'p3']
        assert np.allclose( fisher_1.param_fiducial, [1.0, 2.0, 3.0] )
        assert fisher_1.derived_param_names == ['p4', 'p5']
        assert fisher_1.derived_param_names_latex == ['p4', 'p5']
        assert np.allclose( fisher_1.derived_param_fiducial, [4.0, 5.0] ) 
 
    # test the loading of the parameter names from file:
    def test_load_paramnames_from_file_5(self):
        matrix_derived = np.zeros( (3, 4) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        param_names = [ 'q'+str(i) for i in range(3) ]   
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names = [ 'qd'+str(i) for i in range(4) ]   
        derived_param_names_latex = [ 'md'+str(i) for i in range(4) ]
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(4) ]
        fisher_1 = fd.fisher_derived( derived_matrix=matrix_derived )
        fisher_1.load_paramnames_from_file( file_name=test_input+'/dummy_paramnames_5_derived.paramnames')
        # test if everything is properly initialized:
        assert fisher_1.param_names == ['p1', 'p2', 'p3']     
        assert fisher_1.param_names_latex == ['p_1', 'p_2', 'p_3']
        assert np.allclose( fisher_1.param_fiducial, [1.0, 2.0, 3.0] )
        assert fisher_1.derived_param_names == ['p4', 'p5', 'p6', 'p7']
        assert fisher_1.derived_param_names_latex == ['p_4', 'p_5', 'p_6', 'p7']
        assert np.allclose( fisher_1.derived_param_fiducial, [4.0, 5.0,0.0,0.0] ) 
    
    # test the loading of the parameter names from file:
    def test_load_paramnames_invalid_num(self):
        matrix_derived = np.zeros( (3, 3) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        param_names = [ 'q'+str(i) for i in range(3) ]   
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names = [ 'qd'+str(i) for i in range(4) ]   
        derived_param_names_latex = [ 'md'+str(i) for i in range(4) ]
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(4) ]
        fisher_1 = fd.fisher_derived( derived_matrix=matrix_derived )
        assert_raises( ValueError, fisher_1.load_paramnames_from_file, file_name=test_input+'/dummy_paramnames_5_derived.paramnames' )

# ***************************************************************************************

class test_add_derived():

    @classmethod
    def setup_class(cls):
        print(color_print.header( __name__+': test_add_derived.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold( __name__+': test_add_derived.teardown_class() -------'))

    def setup(self):
        # generate the Fisher matrix. In this case a simple diagonal matrix.
        matrix = np.identity(3)
        for i in range(3):
            matrix[i,i] = i+1
        # Jacobian matrix:
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.0
        matrix_derived[1,1] = 1.0
        matrix_derived[2,0] = 1.0
        # parameter names:
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names_latex = [ 'md'+str(i) for i in range(2) ]
        # fiducial:
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(2) ]
        # initialize the Fisher type:
        self.fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        # initialize the derived matrix:
        self.fisher_der = fd.fisher_derived( derived_matrix=matrix_derived, 
                                             param_names_latex=param_names_latex, 
                                             derived_param_names_latex=derived_param_names_latex, 
                                             fiducial=fiducial, fiducial_derived=fiducial_derived)
        
    def test_invalid_matrix(self):
        assert_raises( ValueError, self.fisher_der.add_derived, fisher_matrix=[0.0] )
    
    def test_invalid_names(self):
        fake_names = [ 'm'+str(i) for i in range(3) ]
        self.fisher_1.set_param_names( fake_names )
        assert_raises( ValueError, self.fisher_der.add_derived, fisher_matrix=self.fisher_1 )
        
    def test_invalid_fiducial(self):
        fake_fiducial = [ float(i+1) for i in range(3) ]
        self.fisher_1.set_fiducial( fake_fiducial )
        assert_raises( ValueError, self.fisher_der.add_derived, fisher_matrix=self.fisher_1 )
        
    def test_preserve_input_false(self):
        fisher_2 = self.fisher_der.add_derived( fisher_matrix=self.fisher_1, preserve_input=False)
        assert np.allclose( fisher_2.fisher_matrix, [[ 0.75, 0.],[ 0., 2.]] )
        assert fisher_2.path == ''
        assert fisher_2.name == '_derived'
        assert fisher_2.indir == ''
        assert fisher_2.num_params == 2
        assert np.allclose( fisher_2.fisher_eigenvalues, [ 0.75, 2.0  ] )
        assert fisher_2.param_names == ['p4', 'p5']
        assert fisher_2.param_names_latex == ['md0', 'md1']
        assert np.allclose( fisher_2.param_fiducial, [ 0., 1.] )
    
    def test_preserve_input_true(self):
        fisher_2 = self.fisher_der.add_derived( fisher_matrix=self.fisher_1, preserve_input=True )
        assert fisher_2.path == ''
        assert fisher_2.name == ''
        assert fisher_2.indir == ''
        assert fisher_2.num_params == 5
        assert fisher_2.param_names == ['p1','p2','p3','p4', 'p5']
        assert fisher_2.param_names_latex == ['m0', 'm1', 'm2', 'md0', 'md1']
        assert np.allclose( fisher_2.param_fiducial, [ 0., 1., 2., 0., 1.] )  
    
    def test_spectral_change(self):
        # redefine the Jacobian matrix:
        matrix_derived = np.zeros( (3, 2) )
        matrix_derived[0,0] = 1.e10
        matrix_derived[2,0] = 1.e-10
        # parameter names:
        param_names_latex = [ 'm'+str(i) for i in range(3) ]
        derived_param_names_latex = [ 'md'+str(i) for i in range(2) ]
        # fiducial:
        fiducial = [ float(i) for i in range(3) ]
        fiducial_derived = [ float(i) for i in range(2) ]
        # initialize the derived matrix:
        fisher_der = fd.fisher_derived( derived_matrix=matrix_derived, 
                                             param_names_latex=param_names_latex, 
                                             derived_param_names_latex=derived_param_names_latex, 
                                             fiducial=fiducial, fiducial_derived=fiducial_derived)
        
        fisher_2 = fisher_der.add_derived( fisher_matrix=self.fisher_1, preserve_input=True )
        
    def test_realistic_add(self):
        # get the original real Fisher:
        fisher   = fm.fisher_matrix( file_name=test_input+'/dummy_fisher_matrix.dat' )
        # get the derived Fisher:
        fisher_1 = fd.fisher_derived( file_name=test_input+'/dummy_fisher_matrix_derived.dat' )
        # add them:
        new_fisher = fisher_1.add_derived( fisher_matrix=fisher, preserve_input=True )
        # check wether the initial spectrum is contained in the final spectrum:
        old_eigen = [ i for i in fisher.fisher_eigenvalues if i> 1.1*fisher.fisher_cutoff ] # remove the cutoff modes
        assert np.amin(new_fisher.fisher_eigenvalues) < np.amin(old_eigen)
        assert np.amax(new_fisher.fisher_eigenvalues) > np.amax(old_eigen)
        
# ***************************************************************************************

