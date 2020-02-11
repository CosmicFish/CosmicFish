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
import cosmicfish_pylib.fisher_plot_settings as fps

from nose.tools import with_setup
from nose.tools import assert_raises

import cosmicfish_pylib.colors as fc
color_print = fc.bash_colors()

# ***************************************************************************************
    
class test_init():

    @classmethod
    def setup_class(cls):
        print(color_print.header( __name__+': test_init.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold( __name__+': test_init.teardown_class() -------'))

    def test_init_settings_default(self):
        settings = fps.CosmicFish_PlotSettings()
        
    def test_init_settings_invalid(self):
        assert_raises( ValueError, fps.CosmicFish_PlotSettings, 'test' )
    
    def test_init_settings_with_dict(self):
        settings = fps.CosmicFish_PlotSettings({'do_legend':False})
        assert settings.do_legend==False
        
    def test_init_settings_with_kwargs(self):
        settings = fps.CosmicFish_PlotSettings( do_legend=False, D2_show_xaxis_label=False )
        assert settings.do_legend==False
        assert settings.D2_show_xaxis_label == False

# ***************************************************************************************
    
class test_update():

    @classmethod
    def setup_class(cls):
        print(color_print.header( __name__+': test_update.setup_class() ----------'))
       
    @classmethod
    def teardown_class(cls):
        print(color_print.bold( __name__+': test_update.teardown_class() -------'))

    def test_update_dict(self):
        settings = fps.CosmicFish_PlotSettings()
        settings.update({'do_legend':False})
        assert settings.do_legend==False
    
    def test_update_kwargs(self):
        settings = fps.CosmicFish_PlotSettings()
        settings.update( do_legend=False, D2_show_xaxis_label=False )
        assert settings.do_legend==False
        assert settings.D2_show_xaxis_label == False
    
    def test_update_raises(self):
        settings = fps.CosmicFish_PlotSettings()
        assert_raises( ValueError, settings.update , 'test' )
        
# ***************************************************************************************

        