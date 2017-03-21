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
import math
# define path to the executable and add all the relevant folders to the path where python looks for files.
here = os.path.dirname(os.path.abspath(__file__))
cosmicfish_pylib_path = here+'/..'
test_input = here+'/test_input'
sys.path.insert(0, os.path.normpath(cosmicfish_pylib_path))

import numpy as np
import cosmicfish_pylib.fisher_matrix as fm
import cosmicfish_pylib.fisher_derived as fd
import cosmicfish_pylib.utilities as fu
import cosmicfish_pylib.fisher_operations as fo

from nose.tools import with_setup
from nose.tools import assert_raises

import cosmicfish_pylib.colors as fc
color_print = fc.bash_colors()

# ***************************************************************************************
    
class test_eliminate_columns_rows():

    @classmethod
    def setup_class(cls):
        print color_print.header( __name__+': test_eliminate_columns_rows.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold( __name__+': test_eliminate_columns_rows.teardown_class() -------')

    def setup(self):
        # generate the Fisher matrix. In this case a simple diagonal matrix.
        matrix = np.identity(4)
        for i in xrange(4):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i+1) for i in xrange(4) ] 
        fiducial = [ float(i) for i in xrange(4) ]
        # initialize the Fisher type:
        self.fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )

    def test_invalid_input(self):
        assert_raises( ValueError, fo.eliminate_columns_rows, 0.0, [0,1] )

    def test_eliminate_column_rows_1(self):
        fisher_2 = fo.eliminate_columns_rows( self.fisher_1, [0,1] )
        # quality check:
        assert np.allclose( fisher_2.fisher_matrix, [[ 3., 0.],[ 0., 4.]] )
        assert fisher_2.path == ''
        assert fisher_2.name == '_reduced'
        assert fisher_2.indir == ''
        assert fisher_2.num_params == 2
        assert np.allclose( fisher_2.fisher_eigenvalues, [ 3., 4.] )
        assert fisher_2.param_names == ['p3', 'p4']
        assert fisher_2.param_names_latex == ['m3', 'm4']
        assert np.allclose( fisher_2.param_fiducial, [ 2., 3.] )

# ***************************************************************************************
    
class test_eliminate_parameters():

    @classmethod
    def setup_class(cls):
        print color_print.header( __name__+': test_eliminate_parameters.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold( __name__+': test_eliminate_parameters.teardown_class() -------')

    def setup(self):
        # generate the Fisher matrix. In this case a simple diagonal matrix.
        matrix = np.identity(3)
        for i in xrange(3):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
        fiducial = [ float(i) for i in xrange(3) ]
        # initialize the Fisher type:
        self.fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )

    def test_invalid_input(self):
        assert_raises( ValueError, fo.eliminate_parameters, 0.0, ['p1'] )
    
    def test_invalid_param(self):
        assert_raises( ValueError, fo.eliminate_parameters, self.fisher_1, ['pfake'] )

    def test_eliminate_parameters_1(self):
        fisher_2 = fo.eliminate_parameters( self.fisher_1, ['p1'] )
        # quality check:
        assert np.allclose( fisher_2.fisher_matrix, [[ 2., 0.],[ 0., 3.]] )
        assert fisher_2.path == ''
        assert fisher_2.name == '_reduced'
        assert fisher_2.indir == ''
        assert fisher_2.num_params == 2
        assert np.allclose( fisher_2.fisher_eigenvalues, [ 2., 3.] )
        assert fisher_2.param_names == ['p2', 'p3']
        assert fisher_2.param_names_latex == ['m1', 'm2']
        assert np.allclose( fisher_2.param_fiducial, [ 1., 2.] )

# ***************************************************************************************
    
class test_reshuffle():

    @classmethod
    def setup_class(cls):
        print color_print.header( __name__+': test_reshuffle.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold( __name__+': test_reshuffle.teardown_class() -------')

    def setup(self):
        # generate the Fisher matrix. In this case a simple diagonal matrix.
        matrix = np.identity(3)
        for i in xrange(3):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
        fiducial = [ float(i) for i in xrange(3) ]
        # initialize the Fisher type:
        self.fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        
    def test_invalid_input(self):
        assert_raises( ValueError, fo.reshuffle, 0.0, ['p1'] )
        
    def test_invalid_param(self):
        assert_raises( ValueError, fo.reshuffle, self.fisher_1, ['pfake'] )
    
    def test_eliminate_parameters_1(self):
        fisher_2 = fo.reshuffle( self.fisher_1, ['p1','p2'] )
        # quality check:
        assert np.allclose( fisher_2.fisher_matrix, [[ 1., 0.],[ 0., 2.]] )
        assert fisher_2.path == ''
        assert fisher_2.name == '_reshuffled'
        assert fisher_2.indir == ''
        assert fisher_2.num_params == 2
        assert np.allclose( fisher_2.fisher_eigenvalues, [ 1., 2.] )
        assert fisher_2.param_names == ['p1', 'p2']
        assert fisher_2.param_names_latex == ['m0', 'm1']
        assert np.allclose( fisher_2.param_fiducial, [ 0., 1.] )
        
# ***************************************************************************************

class test_marginalise():

    @classmethod
    def setup_class(cls):
        print color_print.header( __name__+': test_marginalise.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold( __name__+': test_marginalise.teardown_class() -------')

    def setup(self):
        # generate the Fisher matrix. In this case a simple diagonal matrix.
        matrix = np.identity(3)
        for i in xrange(3):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
        fiducial = [ float(i) for i in xrange(3) ]
        # initialize the Fisher type:
        self.fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        
    def test_invalid_input(self):
        assert_raises( ValueError, fo.marginalise, 0.0, ['p1'] )
        
    def test_invalid_param(self):
        assert_raises( ValueError, fo.marginalise, self.fisher_1, ['pfake'] )
    
    def test_eliminate_parameters_1(self):
        fisher_2 = fo.marginalise( self.fisher_1, ['p1','p2'] )
        # quality check:
        assert np.allclose( fisher_2.fisher_matrix, [[ 1., 0.],[ 0., 2.]] )
        assert fisher_2.path == ''
        assert fisher_2.name == '_marginal'
        assert fisher_2.indir == ''
        assert fisher_2.num_params == 2
        assert np.allclose( fisher_2.fisher_eigenvalues, [ 1., 2.] )
        assert fisher_2.param_names == ['p1', 'p2']
        assert fisher_2.param_names_latex == ['m0', 'm1']
        assert np.allclose( fisher_2.param_fiducial, [ 0., 1.] )
        
# ***************************************************************************************

class test_marginalise_over():

    @classmethod
    def setup_class(cls):
        print color_print.header( __name__+': test_marginalise_over.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold( __name__+': test_marginalise_over.teardown_class() -------')

    def setup(self):
        # generate the Fisher matrix. In this case a simple diagonal matrix.
        matrix = np.identity(3)
        for i in xrange(3):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
        fiducial = [ float(i) for i in xrange(3) ]
        # initialize the Fisher type:
        self.fisher_1 = fm.fisher_matrix( fisher_matrix=matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        
    def test_invalid_input(self):
        assert_raises( ValueError, fo.marginalise_over, 0.0, ['p1'] )
        
    def test_invalid_param(self):
        assert_raises( ValueError, fo.marginalise_over, self.fisher_1, ['pfake'] )
    
    def test_eliminate_parameters_1(self):
        fisher_2 = fo.marginalise_over( self.fisher_1, ['p1','p2'] )
        # quality check:
        assert np.allclose( fisher_2.fisher_matrix, [[3.]] )
        assert fisher_2.path == ''
        assert fisher_2.name == '_marginal'
        assert fisher_2.indir == ''
        assert fisher_2.num_params == 1
        assert np.allclose( fisher_2.fisher_eigenvalues, [ 3.] )
        assert fisher_2.param_names == ['p3']
        assert fisher_2.param_names_latex == ['m2']
        assert np.allclose( fisher_2.param_fiducial, [ 2.] )

# ***************************************************************************************

class test_information_gain():

    @classmethod
    def setup_class(cls):
        print color_print.header( __name__+': test_information_gain.setup_class() ----------')
       
    @classmethod
    def teardown_class(cls):
        print color_print.bold( __name__+': test_information_gain.teardown_class() -------')

    def setup(self):
        # generate the Fisher matrix. In this case a simple diagonal matrix.
        matrix = np.identity(3)
        for i in xrange(3):
            matrix[i,i] = i+1
        param_names_latex = [ 'm'+str(i) for i in xrange(3) ] 
        fiducial = [ float(i) for i in xrange(3) ]
        # initialize the Fisher type:
        self.fisher_1 = fm.fisher_matrix( fisher_matrix=1.0*matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        self.fisher_2 = fm.fisher_matrix( fisher_matrix=10.0*matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        self.fisher_3 = fm.fisher_matrix( fisher_matrix=0.0*matrix, param_names_latex=param_names_latex, fiducial=fiducial )
        
    def test_same_distro(self):
        assert np.isclose( fo.information_gain( self.fisher_1, self.fisher_1, self.fisher_1, stat=True  ), 3.0/4.0*3.0/2.0/math.log(2) )
        assert np.isclose( fo.information_gain( self.fisher_1, self.fisher_1, self.fisher_1, stat=False ), 0.0 )

    def test_different_distro(self):
        """
        print fo.information_gain( self.fisher_1, self.fisher_2, self.fisher_3, stat=True  )
        print fo.information_gain( self.fisher_1, self.fisher_2, self.fisher_3, stat=False )
        print
        print fo.information_gain( self.fisher_2, self.fisher_1, self.fisher_3, stat=True  )
        print fo.information_gain( self.fisher_2, self.fisher_1, self.fisher_3, stat=False )
        """
        pass
        
# ***************************************************************************************