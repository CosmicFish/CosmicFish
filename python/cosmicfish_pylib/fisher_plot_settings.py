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

"""
.. module:: fisher_plot_settings
   :platform: Unix
   :synopsis: Module that contains the class for plot settings and their default values. 
   
.. moduleauthor:: Marco Raveri <mraveri@uchicago.edu> for the CosmicFish code.

"""

# ***************************************************************************************

import matplotlib

# ***************************************************************************************

""" 
Some global parameters 
"""

paper_width = 21.0                #: width of an A4 paper in cm

inch_to_cm  = 2.54                #: conversion factor from inches to centimeter
cm_to_inch  = 1.0/inch_to_cm      #: conversion factor from centimeter to inches 

points_to_mm = 0.352778           #: conversion factor from points to mm 
mm_to_points = 1.0/points_to_mm   #: conversion factor from mm to points 

points_width_to_mm = 0.18
points_heights_to_mm = 0.30


# ***************************************************************************************

class CosmicFish_PlotSettings():
    """
    Class containing the plot settings for the CosmicFish plotter.
    
    :ivar use_nice_colors: wether to use nice_colors for the plot. Not yet implemented.
    :ivar nice_colors_palette: nice_color palette to use. Default:
    :ivar line_colors: list of color to use for lines. Default:
    :ivar solid_colors: list of color to use for filled stuff. Default:
    :ivar linestyle: list of line style to use. Default:
        
    :ivar D1_confidence_level: list with the confidence levels for the 1D plot. Default: [ 0.95, 0.68 ]
    :ivar D1_num_points: number of points per line in the plot. Default: 100
    :ivar D1_line_thickness: thickness of the line. Default: 1.0
    :ivar D1_norm_prob: wether to normalize the distribution or set the maximum to one. Default: False
    :ivar D1_use_nice_numbers: wether to use nice numbers when doing 1D plots. Default: True
    :ivar D1_filled: wether to fill the gaussians highlighting the confidence levels. Default: True
    :ivar D1_alpha: alpha used for the 1D filling. Each layer is filled with the same alpha and the superposition of the layers result in the highlihting effect. Default: 0.1
    :ivar D1_ylabel_right: wether to put y labels on the right side of the plot. Default: False
    :ivar D1_xlabel_up: wether to put x labels on top of the plot. Default: False

    :ivar D1_show_x_ticks: wether to show x axis thicks in 1D plot. Default: True
    :ivar D1_number_x_ticks: number of x thicks to show. Default: 3
    :ivar D1_show_y_ticks: wether to show y axis thicks in 1D plot. Default: True
    :ivar D1_number_y_ticks: number of y thicks to show. Default: 2
    :ivar D1_show_x_ticks_labels: wether to show x axis tick label. Default: True
    :ivar D1_show_y_ticks_labels: wether to show y axis tick label. Default: True
    
    :ivar D1_show_xaxis_label: wether to show x axis label. Default: True
    :ivar D1_show_yaxis_label: wether to show y axis label. Default: True
    :ivar D1_x_label_rotation: rotation of the x axis label. Default: 0
    :ivar D1_y_label_rotation: rotation of the y axis label. Default: 90    
    :ivar D1_prob_label: y axis label to use when P is not normalized. Default: u'$P/P_{\\rm max}$' 
    :ivar D1_norm_prob_label: y axis label to use when P is normalized. Default: u'$P$'  
    :ivar D1_main_fontsize: main fontsize for the 1D plot. Default: 10.0   
    :ivar D1_secondary_fontsize: secondary fontsize for the 1D plot. Default:
    :ivar D1_show_best_fit: show a vertical line on the fiducial. Default: False

    :ivar D2_confidence_levels: list with confidence levels for the 2D plot. Default: [ 0.95, 0.68 ]      
    :ivar D2_num_points: number of points per line in the plot. Default: 100 
    :ivar D2_line_thickness: thickness of the lines. Default: 1.0 
    :ivar D2_use_nice_numbers: wether to use nice numbers when doing 2D plots. Default: True
    :ivar D2_filled: wether to fill the contours highlighting the confidence levels. Default: True
    :ivar D2_alphas: alpha for the filling of the confidence regions. Default: [ 0.40, 0.85 ]      
    :ivar D2_ylabel_right: wether to put y labels on the right side of the plot. Default: False
    :ivar D2_xlabel_up: wether to put x labels on top of the plot. Default: False

    :ivar D2_show_x_ticks: wether to show x axis thicks in 1D plot. Default: True
    :ivar D2_number_x_ticks: number of x thicks to show. Default: 3
    :ivar D2_show_y_ticks: wether to show y axis thicks in 1D plot. Default: True
    :ivar D2_number_y_ticks: number of y thicks to show. Default: 3
    :ivar D2_show_x_ticks_labels: wether to show x axis tick label. Default: True
    :ivar D2_show_y_ticks_labels: wether to show y axis tick label. Default: True

    :ivar D2_show_xaxis_label: wether to show x axis label. Default: True
    :ivar D2_show_yaxis_label: wether to show y axis label. Default: True
    :ivar D2_x_label_rotation: rotation of the x axis label. Default: 0
    :ivar D2_y_label_rotation: rotation of the y axis label. Default: 90 
    :ivar D2_main_fontsize: main fontsize. Default: 10.0
    :ivar D2_secondary_fontsize: secondary fontsize. Default: 9.0
    :ivar D2_show_best_fit: show an intersecting line on the fiducial. Default: False
        
    :ivar num_plots_per_line: number of plots per line. Used in all plots but triangulars. Default: 3
    :ivar figure_width: width of the image in cm. Default: paper_width
    :ivar figure_height: height of the image in cm. Default: paper_width
    :ivar tight_layout: wether to use tight_layout at the end of plotting. Default: False
    :ivar subplot_x_size: x size of the subplot in cm. Default: 5.0 cm
    :ivar subplot_y_size: y size of the subplot in cm. Default: 5.0 cm
    :ivar use_fixed_figure_width: wether to use a fixed size for the figure width. If true the width of the subplots is determined automatically, if false is determined by subplot_x_size. Default: False
    :ivar use_fixed_figure_height: wether to use a fixed size for the figure height. If true the height of the subplots is determined automatically, if false is determined by subplot_y_size. Default: False

    :ivar do_legend: wether to put the legend or not. Default: True
    :ivar legend_filled: wether to use a legend with filled patches. Default: True
    :ivar legend_fontsize: font size to be used in the legend. Default: 10.0     
    :ivar legend_frame: wether to put the frame on the legend. Default: True
    :ivar legend_fancybox: wether to use a fancy box for the legend. Default: False
    :ivar legend_shadow: wether to put the shadow below the legend. Default: False
    :ivar legend_markerfirst: wether to put the marker first and the text after. Default: True
    :ivar legend_ncol: number of columns of the legend. Default: 1
    :ivar legend_loc: location of the legend, relative to the global figure if legend_takes_place_plot is False or relative to the legend sub plot if legend_takes_place_plot is True. Default: 'lower center'
    :ivar legend_takes_place_plot: wether the legend should take the place of an additional subplot in 1D and 2D plots. Default: False
    :ivar legend_takes_place_plot_tri: wether the legend should take the place of an additional subplot in a triangular plot. Default: True
    
    :ivar title_fontsize: fontsize of the title of the plot. Default: 10.0  
    
    """

    # -----------------------------------------------------------------------------------
    
    def __init__( self, dictionary=None, **kwargs ):
        """
        **CosmicFish_PlotSettings class constructor**. The constructor initializes the plot 
        settings to the default values.
        By default everything (font sizes, ecc..) is chosen so that every plot looks good to us...
        The settings can be initialized to default by passing nothing.
        
        :param dictionary: dictionary containing the users settings. Settings that are not passed
            are initialized to their default values. If a setting is not in __accepted_settings__ it
            will simply be ignored.
        :type dictionary: :class:`dict`
        :param kwargs: keyword arguments containing settings. Settings that are not passed in this way
            are either initialized to default or by dictionary. If a setting is not in __accepted_settings__ it
            will simply be ignored. Notice that keyword arguments have higher priority with respect to 
            dictionary settings.
        
        """
        
        # colors:
        self.use_nice_colors     = False
        self.nice_colors_palette = None
        self.line_colors         = ['#006FED', '#E03424', 'black', '#009966', '#000866', '#336600', '#006633', 'm','r']  
        self.solid_colors        = ['#006FED', '#E03424', 'gray', '#009966', '#000866', '#336600', '#006633', 'm','r']  
        self.linestyle           = ['-','-','-','-','-','-','-','-','-'] 
        
        # options for 1D plots:
        self.D1_confidence_level     = [ 0.95, 0.68 ]      
        self.D1_num_points           = 100                 
        self.D1_line_thickness       = 1.0                 
        self.D1_norm_prob            = False               
        self.D1_use_nice_numbers     = True                
        self.D1_filled               = True                
        self.D1_alpha                = 0.1                 
        self.D1_ylabel_right         = False               
        self.D1_xlabel_up            = False               
        # - thicks
        self.D1_show_x_ticks         = True                
        self.D1_number_x_ticks       = 3                   
        self.D1_show_y_ticks         = True                
        self.D1_number_y_ticks       = 2                   
        self.D1_show_x_ticks_labels  = True                
        self.D1_show_y_ticks_labels  = True                
        # - labels
        self.D1_show_xaxis_label     = True                
        self.D1_show_yaxis_label     = True                
        self.D1_x_label_rotation     = 0                   
        self.D1_y_label_rotation     = 90                  
        self.D1_prob_label           = '$P/P_{\\rm max}$' 
        self.D1_norm_prob_label      = '$P$'              
        self.D1_main_fontsize        = 10.0                
        self.D1_secondary_fontsize   = 9.0                 
        self.D1_show_best_fit        = False               
    
        # options for 2D plots:
        self.D2_confidence_levels    = [ 0.95, 0.68 ]      
        self.D2_num_points           = 100                 
        self.D2_line_thickness       = 1.0                 
        self.D2_use_nice_numbers     = True                
        self.D2_filled               = True                
        self.D2_alphas               = [ 0.40, 0.85 ]      
        self.D2_ylabel_right         = False               
        self.D2_xlabel_up            = False               
        # - thicks
        self.D2_show_x_ticks         = True                
        self.D2_number_x_ticks       = 3                   
        self.D2_show_y_ticks         = True                
        self.D2_number_y_ticks       = 3                   
        self.D2_show_x_ticks_labels  = True                
        self.D2_show_y_ticks_labels  = True                
        # - labels
        self.D2_show_xaxis_label     = True               
        self.D2_show_yaxis_label     = True                
        self.D2_x_label_rotation     = 0                   
        self.D2_y_label_rotation     = 90                  
        self.D2_main_fontsize        = 10.0                
        self.D2_secondary_fontsize   = 9.0                 
        self.D2_show_best_fit        = False               
        
        # options for the global figure:
        self.num_plots_per_line      = 3                   
        self.figure_width            = paper_width         
        self.figure_height           = paper_width         
        self.tight_layout            = False               
        self.subplot_x_size          = 5.0                 
        self.subplot_y_size          = 5.0                 
        self.use_fixed_figure_width  = False               
        self.use_fixed_figure_height = False               
        
        # options for the legend:
        self.do_legend               = True                
        self.legend_filled           = True                
        self.legend_fontsize         = 10.0                
        self.legend_frame            = True                
        self.legend_fancybox         = False               
        self.legend_shadow           = False               
        self.legend_markerfirst      = True                
        self.legend_ncol             = 1                   
        self.legend_loc              = 'lower center'      
        self.legend_takes_place_plot = False               
        self.legend_takes_place_plot_tri = True
        
        # options for the title:
        self.title_fontsize          = 10.0                

        # matplotlib relevant options:
        matplotlib.rcParams['savefig.transparent'] = True
        matplotlib.rcParams['figure.autolayout']   = False
        matplotlib.rcParams['xtick.major.pad']     = 4.0
        matplotlib.rcParams['ytick.major.pad']     = 4.0
        matplotlib.rcParams['axes.labelpad']       = 5.0
        matplotlib.rcParams['savefig.dpi']         = 300
        matplotlib.rcParams['savefig.pad_inches']  = 0.0
        
        # latex options:
        matplotlib.rcParams['text.usetex'] = True
        
        # replace settings with those passed from dictionary:
        if dictionary is not None:
            if not isinstance(dictionary, dict):
                raise ValueError('Invalid input. Dictionary shall be a dictionary of settings.')
            else:
                temp = [ key for key in list(dictionary.keys()) if key in __accepted_settings__ ]
                for key in temp:
                    setattr(self, key, dictionary[key])
        # now parse kwargs. Esplicit kwarg takes precedence wrt dictionary:
        if kwargs is not None:
            temp = [ key for key in list(kwargs.keys()) if key in __accepted_settings__ ]
            for key in temp:
                setattr(self, key, kwargs[key])
            
    # -----------------------------------------------------------------------------------
    
    def update(self, dictionary=None, **kwargs):
        """
        Update settings. Works very similarly to the class constructor with the difference
        being that settings that are not passed are not initialized to default.
        
        :param dictionary: dictionary containing the users settings. Settings that are not passed
            are not changed. If a setting is not in __accepted_settings__ it
            will simply be ignored.
        :type dictionary: :class:`dict`
        :param kwargs: keyword arguments containing settings. Settings that are not passed are not changed.
            If a setting is not in __accepted_settings__ it will simply be ignored. 
            Notice that keyword arguments have higher priority with respect to 
            dictionary settings.
            
        """
        # replace settings with those passed from dictionary:
        if dictionary is not None:
            if not isinstance(dictionary, dict):
                raise ValueError('Invalid input. Dictionary shall be a dictionary of settings.')
            else:
                temp = [ key for key in list(dictionary.keys()) if key in __accepted_settings__ ]
                for key in temp:
                    setattr(self, key, dictionary[key])
        # now parse kwargs. Esplicit kwarg takes precedence wrt dictionary:
        if kwargs is not None:
            temp = [ key for key in list(kwargs.keys()) if key in __accepted_settings__ ]
            for key in temp:
                setattr(self, key, kwargs[key])
        
    # -----------------------------------------------------------------------------------

# ***************************************************************************************

__accepted_settings__=[ 
    'use_nice_colors',
    'nice_colors_palette', 
    'line_colors',    
    'solid_colors',
    'linestyle',  
    'D1_confidence_level',      
    'D1_num_points',       
    'D1_line_thickness',    
    'D1_norm_prob',
    'D1_use_nice_numbers',
    'D1_filled',
    'D1_alpha',
    'D1_ylabel_right',
    'D1_xlabel_up',
    'D1_show_x_ticks',
    'D1_number_x_ticks',
    'D1_show_y_ticks',
    'D1_number_y_ticks',
    'D1_show_x_ticks_labels',
    'D1_show_y_ticks_labels',
    'D1_show_xaxis_label',
    'D1_show_yaxis_label',
    'D1_x_label_rotation',
    'D1_y_label_rotation',
    'D1_prob_label',
    'D1_norm_prob_label',
    'D1_main_fontsize',
    'D1_secondary_fontsize',
    'D1_show_best_fit',
    'D2_confidence_levels',
    'D2_num_points',
    'D2_line_thickness',
    'D2_use_nice_numbers',
    'D2_filled',
    'D2_alphas',
    'D2_ylabel_right',
    'D2_xlabel_up',
    'D2_show_x_ticks',      
    'D2_number_x_ticks',     
    'D2_show_y_ticks',       
    'D2_number_y_ticks',      
    'D2_show_x_ticks_labels',
    'D2_show_y_ticks_labels',
    'D2_show_xaxis_label',
    'D2_show_yaxis_label',
    'D2_x_label_rotation',
    'D2_y_label_rotation',
    'D2_main_fontsize',
    'D2_secondary_fontsize',
    'D2_show_best_fit',
    'num_plots_per_line',    
    'figure_width',
    'figure_height',
    'tight_layout',
    'subplot_x_size',
    'subplot_y_size',
    'use_fixed_figure_width',
    'use_fixed_figure_height',
    'do_legend',        
    'legend_filled',
    'legend_fontsize',
    'legend_frame',     
    'legend_fancybox',       
    'legend_shadow',
    'legend_markerfirst',
    'legend_ncol',
    'legend_loc',
    'legend_takes_place_plot',
    'legend_takes_place_plot_tri',
    'title_fontsize',
    ]

# ***************************************************************************************
