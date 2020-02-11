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
import cosmicfish_pylib.utilities as fu

from nose.tools import with_setup
from nose.tools import assert_raises

import cosmicfish_pylib.colors as fc
color_print = fc.bash_colors()

# ***************************************************************************************

class test_nice_number():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_nice_number.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_nice_number.teardown_class() -------'))

    def setup(self):
        pass
    
    def teardown(self):
        pass 
        
    def test_nice_number_ceil(self):
        assert fu.nice_number(0.1234, mode=0) == 0.2
        assert fu.nice_number(-0.1234, mode=0) == -0.1
        
    def test_nice_number_round(self):
        assert fu.nice_number(0.1234, mode=1) == 0.1
        assert fu.nice_number(-0.1234, mode=1) == -0.1

    def test_nice_number_floor(self):
        assert fu.nice_number(0.1234, mode=2) == 0.1
        assert fu.nice_number(-0.1234, mode=2) == -0.2
    
    def test_nice_number_invalid(self):
        assert_raises( ValueError, fu.nice_number, -0.1234, mode=4 )  

# ***************************************************************************************

class test_significant_digits():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_nice_number.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_nice_number.teardown_class() -------'))

    def setup(self):
        pass
    
    def teardown(self):
        pass 
        
    def test_test_significant_digits_ceil_plus(self):
        assert fu.significant_digits( (1.234,0.012), mode=0 ) == 1.24
        assert fu.significant_digits( (1.234,10.23), mode=0 ) == 10.0
        assert fu.significant_digits( (1.234,1.234), mode=0 ) == 2.0
        assert fu.significant_digits( (0.000145,0.000123), mode=0 ) == 0.0002
        assert fu.significant_digits( (0.000145,0.000000000123), mode=0 ) == 0.000145
    
    def test_test_significant_digits_ceil_minus(self):
        assert fu.significant_digits( (-1.234,0.012), mode=0 ) == -1.23
        assert fu.significant_digits( (-1.234,10.23), mode=0 ) == -0.0
        assert fu.significant_digits( (-1.234,1.234), mode=0 ) == -1.0
        assert fu.significant_digits( (-0.000145,0.000123), mode=0 ) == -0.0001
        assert fu.significant_digits( (-0.000145,0.000000000123), mode=0 ) == -0.000145
    
    def test_test_significant_digits_round_plus(self):
        assert fu.significant_digits( (1.234,0.012), mode=1 ) == 1.23
        assert fu.significant_digits( (1.234,10.23), mode=1 ) == 0.0
        assert fu.significant_digits( (1.234,1.234), mode=1 ) == 1.0
        assert fu.significant_digits( (0.000145,0.000123), mode=1 ) == 0.0001
        assert fu.significant_digits( (0.000145,0.000000000123), mode=1 ) == 0.000145
    
    def test_test_significant_digits_round_minus(self):
        assert fu.significant_digits( (-1.234,0.012), mode=1 ) == -1.23
        assert fu.significant_digits( (-1.234,10.23), mode=1 ) == -0.0
        assert fu.significant_digits( (-1.234,1.234), mode=1 ) == -1.0
        assert fu.significant_digits( (-0.000145,0.000123), mode=1 ) == -0.0001
        assert fu.significant_digits( (-0.000145,0.000000000123), mode=1 ) == -0.000145
        
    def test_test_significant_digits_floor_plus(self):
        assert fu.significant_digits( (1.234,0.012), mode=2 ) == 1.23
        assert fu.significant_digits( (1.234,10.23), mode=2 ) == 0.0
        assert fu.significant_digits( (1.234,1.234), mode=2 ) == 1.0
        assert fu.significant_digits( (0.000145,0.000123), mode=2 ) == 0.0001
        assert fu.significant_digits( (0.000145,0.000000000123), mode=2 ) == 0.000145
        
    def test_test_significant_digits_floor_minus(self):
        assert fu.significant_digits( (-1.234,0.012), mode=2 ) == -1.24
        assert fu.significant_digits( (-1.234,10.23), mode=2 ) == -10.0
        assert fu.significant_digits( (-1.234,1.234), mode=2 ) == -2.0
        assert fu.significant_digits( (-0.000145,0.000123), mode=2 ) == -0.0002
        assert fu.significant_digits( (-0.000145,0.000000000123), mode=2 ) == -0.000145
        
    def test_test_significant_digits_error(self):
        assert_raises( ValueError, fu.significant_digits, (1.234,0.012), mode=3 )  

# ***************************************************************************************

class test_make_list():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_make_list.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_make_list.teardown_class() -------'))

    def setup(self):
        pass
    
    def teardown(self):
        pass 
        
    def test_make_list_with_list(self):
        output = fu.make_list( [0,1,2] )
        assert isinstance( output, list )
        assert np.allclose( output, [0,1,2] )
        
    def test_make_list_without_list(self):
        output = fu.make_list( 0 )
        assert isinstance( output, list )
        assert np.allclose( output, [0] )

# ***************************************************************************************

class test_print_table():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_print_table.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_print_table.teardown_class() -------'))

    def setup(self):
        pass
    
    def teardown(self):
        pass 
        
    def test_print_table(self):
        fu.print_table( [['test','values','of','table'],[0.0,1.0,2.0,3.0]] )

# ***************************************************************************************

class test_confidence_coefficient():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_confidence_coefficient.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_confidence_coefficient.teardown_class() -------'))

    def setup(self):
        pass
    
    def teardown(self):
        pass 
        
    def test_confidence_coefficient_1(self):
        assert np.allclose( fu.confidence_coefficient( 0.0 ), 0.0 )
    
# ***************************************************************************************

class test_num_to_mant_exp():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_num_to_mant_exp.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_num_to_mant_exp.teardown_class() -------'))

    def setup(self):
        pass
    
    def teardown(self):
        pass 
        
    def test_num_to_mant_exp_raises(self):
        assert fu.num_to_mant_exp(0.0) == (0,0)
    
# ***************************************************************************************

class test_grouper():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_grouper.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_grouper.teardown_class() -------'))

    def setup(self):
        pass
    
    def teardown(self):
        pass 
        
    def test_grouper_1(self):
        list = [i for i in range(9)]
        assert fu.grouper(2, list) == [(0, 1), (2, 3), (4, 5), (6, 7), (8, None)]
        assert fu.grouper(2, list, fillvalue=0) == [(0, 1), (2, 3), (4, 5), (6, 7), (8, 0)]
        assert fu.grouper(3, list, fillvalue=0) == [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
         
# ***************************************************************************************

class test_CosmicFish_write_header():

    @classmethod
    def setup_class(cls):
        print(color_print.header(__name__+': test_CosmicFish_write_header.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold(__name__+': test_CosmicFish_write_header.teardown_class() -------'))

    def setup(self):
        pass
    
    def teardown(self):
        pass 
        
    def test_CosmicFish_write_header(self):
        fu.CosmicFish_write_header(' Test header.')
         
# ***************************************************************************************
