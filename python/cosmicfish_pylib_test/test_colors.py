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

# This file contains the test for cosmicfish_pylib.colors 

# import first packages
import os, sys
# define path to the executable and add all the relevant folders to the path where python looks for files.
here = os.path.dirname(os.path.abspath(__file__))
cosmicfish_pylib_path = here+'/..'
test_input = here+'/test_input'
sys.path.insert(0, os.path.normpath(cosmicfish_pylib_path))

import cosmicfish_pylib.colors as col

# ***************************************************************************************

def test_nice_colors():
    for i in range(10):
        # check that the entrances are 3:
        assert len(col.nice_colors(i))==3
        # check that they are positive and smaller than 1, i.e. it is legal RGB:
        for j in col.nice_colors(i):
            assert j>=0.0
            assert j<=1.0
            
# ***************************************************************************************

def test_bash_colors():
    color_print = col.bash_colors()
    print() 
    print(color_print.header(__name__+'.test_bash_colors: Header color'))
    print(color_print.blue(__name__+'.test_bash_colors: Blue color'))
    print(color_print.green(__name__+'.test_bash_colors: Green color'))
    print(color_print.warning(__name__+'.test_bash_colors: Warning color'))
    print(color_print.fail(__name__+'.test_bash_colors: Fail color'))
    print(color_print.bold(__name__+'.test_bash_colors: Bold color'))
    print(color_print.underline(__name__+'.test_bash_colors: Underline color'))
    
# ***************************************************************************************
