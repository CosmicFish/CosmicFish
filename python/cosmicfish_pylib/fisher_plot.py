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
.. module:: fisher_plot
   :platform: Unix
   :synopsis: Module that contains a set of plotting tools. 

.. moduleauthor:: Marco Raveri <mraveri@uchicago.edu> for the CosmicFish code.

"""

# ***************************************************************************************

import os
import math
import copy
import itertools as it
import numpy     as np
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot   as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches  as mpatches
import matplotlib.lines    as mlines

import utilities     as fu
import fisher_matrix as fm
import fisher_plot_analysis as fpa

from fisher_plot_settings import *

# ***************************************************************************************

class CosmicFishPlotter():
    """
    Main class for making plots from one or more Fisher matrices. All the functions create
    matplotlib objects that can later be accessed and modified.
    
    :ivar plot_fishers: :class:`cosmicfish_pylib.fisher_plot_analysis.CosmicFish_FisherAnalysis` containing the list of Fisher matrices that are being plotted.
    :ivar plot_settings: :class:`cosmicfish_pylib.fisher_plot_settings.CosmicFish_PlotSettings` containing the plot settings that are used if not overwritten when plotting.
    :ivar figure: :class:`matplotlib.figure.Figure` main figure object.
    :ivar legend: :class:`matplotlib.legend.Legend` legend object.
    :ivar title: :class:`matplotlib.text.Text` title object.
    :ivar plot_grid: :class:`matplotlib.gridspec.GridSpec` grid with the subplots. Created at initialization as empty and then arranged by plotting functions.
    :ivar plot_dict: :class:`dict` dictionary containing the mapping of parameters and subplots :class:`matplotlib.axes._subplots.AxesSubplot`.
    :ivar plot_number: :class:`int` number of plots.
    :ivar bind_line_colors: :class:`dict` dictionary mapping the names of the Fisher matrices to line colors. This is used to ensure consistency of settings and Fishers across plots.
    :ivar bind_solid_colors: :class:`dict` dictionary mapping the names of the Fisher matrices to solid colors. This is used to ensure consistency of settings and Fishers across plots.
    :ivar bind_labels: :class:`dict` dictionary mapping the names of the Fisher matrices to labels. This is used to ensure consistency of settings and Fishers across plots.
    :ivar bind_linestyle: :class:`dict` dictionary mapping the names of the Fisher matrices to line styles. This is used to ensure consistency of settings and Fishers across plots.
    
    """
    
    # -----------------------------------------------------------------------------------

    def __init__(self, settings=None, fishers=None, **kwargs):
        """
        **CosmicFishPlotter class constructor**. The constructor will create the plotting
        objects based on the passed settings and Fisher matrices.
        
        :param settings: an instance of CosmicFish_PlotSettings or a dictionary containing the plot settings.
        :type settings: :class:`cosmicfish_pylib.fisher_plot_settings.CosmicFish_PlotSettings` or :class:`dict`
        :param fishers:  an instance of fisher_plot.analysis.CosmicFish_FisherAnalysis containing a series of
            fisher matrices.
        :type fishers: :class:`cosmicfish_pylib.fisher_plot_analysis.CosmicFish_FisherAnalysis`
        :param kwargs: these contain the possibility of passing additional options or other options that would
            override the ones in settings.

        """
        # initialize the objects of the class:
        self.plot_fishers  = fpa.CosmicFish_FisherAnalysis()
        self.plot_settings = None
        self.figure        = None
        self.legend        = None
        self.title         = None
        self.plot_grid     = None
        self.plot_dict     = {}
        self.plot_number   = 0        
        # get settings:
        if settings is None:
            self.plot_settings = CosmicFish_PlotSettings()
        else:
            if isinstance(settings, CosmicFish_PlotSettings):
                self.plot_settings = settings
            elif isinstance(settings, dict):
                self.plot_settings = CosmicFish_PlotSettings( settings )
            else:
                raise ValueError('Error in initializing the CosmicFishPlotter: settings is not type CosmicFish_PlotSettings nor a dictionary.')
        # update settings with kwargs:
        self.plot_settings.update( **kwargs )
        # get the fisher list:
        if fishers is not None:
            if not isinstance(fishers, fpa.CosmicFish_FisherAnalysis):
                raise ValueError('CosmicFishPlotter error: fishers is not type CosmicFish_FisherAnalysis')
            else:
                self.plot_fishers = copy.deepcopy( fishers )
        
        # bind settings to names initially:
        self.bind_plot_settings_to_names( **kwargs )
    
    # -----------------------------------------------------------------------------------

    def __del__(self):
        """
        CosmicFishPlotter class destructor.
        """
        self.close_plot()
        del self.plot_fishers
        self.plot_settings = None
        self.figure        = None
        self.plot_dict     = {}
    
    # -----------------------------------------------------------------------------------

    def new_plot(self):
        """
        Creates a new plot erasing everything that happened before.
        No need to call close_plot if you want to do another plot.
        """
        self.close_plot()
        self.figure = plt.gcf()

    # -----------------------------------------------------------------------------------

    def close_plot(self):
        """
        Closes the plot erasing everything that happened before.
        """
        plt.cla()
        plt.clf()
        plt.close("all")
        self.plot_dict = {}
        
    # ----------------------------------------------------------------------------------- 
    
    def export(self, filename, **kwargs):
        """
        Export the plot to file. Almost just a wrapper to :class:`matplotlib.savefig`
        
        :param filename: filename and path of the output file.
        :type filename: :class:`string`
        :param kwargs: optional keyword arguments directly fed to plt.savefig.
        
        """
        # if the user does not supply options use ours:
        
        # call savefig:
        plt.savefig( filename, **kwargs )
        
    # ----------------------------------------------------------------------------------- 
    
    def plot1D(self, params=None, names=None, title=None, **kwargs):
        """
        Main 1D plotting function. This takes a set of parameter names, a set of fisher
        names and produces a figure with 1D plots.
        Uses settings in plot_settings but the resulting figure (figure) can be edited
        just like all matplotlib.pyplot figures. Optionally a title can be supplied and
        any kind of settings can be passed by kwargs.
        
        :param params: a name of a parameter or a list of parameter names. If None does all.
        :type params: :class:`string` or a :class:`list` of :class:`string`
        :param names: a name of a fisher matrix or a list of names of fisher matrices. If None does all.
        :type names: :class:`string` or a :class:`list` of :class:`string`
        :param title: the optional title for the plot.
        :type title: :class:`string`
        :param kwargs: optional keyword settings.

        """
        # get params:
        all_names  = copy.deepcopy( self.plot_fishers.fisher_name_list )
        all_params = self.plot_fishers.get_parameter_list()
        if params==None and names==None:
            names_temp  = copy.deepcopy( all_names )
            params_temp = copy.deepcopy( all_params )
        elif params==None and names!=None:
            names_temp  = [ i for i in fu.make_list(names) if i in all_names ]
            params_temp = self.plot_fishers.get_parameter_list( names=names_temp )
        elif params!=None and names==None:
            params_temp = [ i for i in fu.make_list(params) if i in all_params ]
            names_temp  = copy.deepcopy( all_names )
        else:
            names_temp  = [ i for i in fu.make_list(names) if i in all_names ]
            params_temp = [ i for i in fu.make_list(params) if i in self.plot_fishers.get_parameter_list(names_temp) ]
        # check input:
        if len(names_temp) == 0:
            raise ValueError('ERROR PLOT 1D: no valid input fisher matrix.\nPossible values are: '+str(all_names))
        if len(params_temp) == 0:
            raise ValueError('ERROR PLOT 1D: no valid input parameters.\nPossible values are: '+str(all_params))
        # override settings:
        plot_per_line           = self.setting_setter('num_plots_per_line', **kwargs)
        legend_takes_place_plot = self.setting_setter('legend_takes_place_plot', **kwargs)
        # override name bind settings:
        self.bind_plot_settings_to_names( names=names_temp, **kwargs )
        # number of parameters:
        num_plots = len(params_temp)
        if legend_takes_place_plot: num_plots = num_plots +1
        # get the number of rows:
        if num_plots%plot_per_line == 0: 
            num_row = max(1,num_plots/plot_per_line)
        else: 
            num_row = max(1,num_plots/plot_per_line+1)
        # get the number of columns:
        num_col = min( num_plots, plot_per_line )
        # create the layout of the figure:
        self.plot_grid = gridspec.GridSpec( num_row, num_col )
        self.plot_number = num_plots
        # fill the single plotters:
        self.plot_dict = {}
        for i,j in enumerate(params_temp):
            self.plot_dict[j] = plt.subplot( self.plot_grid[i/plot_per_line,i%plot_per_line] )
            self.figure_1D( subplot=self.plot_dict[j], param=j, names=names_temp, **kwargs)
        # do the legend:
        self.set_legend( names=names_temp, **kwargs )
        # do the title:
        self.set_title( title=title, **kwargs)
        # set the global appearence of the plot:
        self.set_plot_dimensions( num_col=num_col, num_rows=num_row, **kwargs )
    
    # -----------------------------------------------------------------------------------
    
    def plot2D(self, params=None, names=None, title=None, **kwargs):
        """
        Main 2D plotting function. This takes a set of parameter names couples, a set of fisher
        names and produces a figure with 2D plots.
        Uses settings in plot_settings but the resulting figure (figure) can be edited
        just like all matplotlib.pyplot figures. Optionally a title can be supplied and
        any kind of settings can be passed by kwargs.
       
        :param params: a couple of names of a parameter or a list of parameter names couples. If None does all.
        :type params: :class:`string` or a :class:`list` of :class:`string`
        :param names: a name of a fisher matrix or a list of names of fisher matrices. If None does all.
        :type names: :class:`string` or a :class:`list` of :class:`string`
        :param title: the optional title for the plot.
        :type title: :class:`string`
        :param kwargs: optional keyword settings.
        
        """
        # get params:
        all_names  = copy.deepcopy( self.plot_fishers.fisher_name_list )
        all_params = self.plot_fishers.get_parameter_list()
        all_legal_params_couples = [ list(i) for i in it.permutations(all_params, 2) ]
        if params==None and names==None:
            names_temp  = copy.deepcopy( all_names )
            params_temp = [ list(i) for i in it.combinations(all_params, 2) ]
        elif params==None and names!=None:
            names_temp  = [ i for i in fu.make_list(names) if i in all_names ]
            params_temp = self.plot_fishers.get_parameter_list( names=names_temp )
            params_temp = [ list(i) for i in it.combinations(params_temp, 2) ]
        elif params!=None and names==None:
            params_temp = [ i for i in fu.make_list(params) if i in all_legal_params_couples ]
            names_temp  = copy.deepcopy( all_names )
        else:
            params_temp = [ i for i in fu.make_list(params) if i in all_legal_params_couples ]
            names_temp  = [ i for i in fu.make_list(names) if i in all_names ]
        # check input:
        if len(params_temp) == 0:
            raise ValueError('ERROR PLOT 2D: no valid input parameters.\nPossible values are: '+str(all_params))
        if len(names_temp) == 0:
            raise ValueError('ERROR PLOT 2D: no valid input fisher matrix.\nPossible values are: '+str(all_names))
        # override settings:
        plot_per_line = self.setting_setter('num_plots_per_line', **kwargs)
        legend_takes_place_plot = self.setting_setter('legend_takes_place_plot', **kwargs)
        # override name bind settings:
        self.bind_plot_settings_to_names( names=names_temp, **kwargs )
        # number of parameters:
        num_plots = len(params_temp)
        if legend_takes_place_plot: num_plots = num_plots +1
        # get the number of rows:
        if num_plots%plot_per_line == 0: 
            num_row = max(1,num_plots/plot_per_line)
        else: 
            num_row = max(1,num_plots/plot_per_line+1)
        # get the number of columns:
        num_col = min( num_plots, plot_per_line )
        # create the layout of the figure:
        self.plot_grid = gridspec.GridSpec( num_row, num_col )
        self.plot_number = num_plots
        # fill the single plotters:
        self.plot_dict = {}
        for i,j in enumerate(params_temp):
            # generate the dictionary: we want symmetric keys just to be sure.
            self.plot_dict['['+str(j[0])+','+str(j[1])+']'] = plt.subplot(self.plot_grid[i/plot_per_line,i%plot_per_line])
            self.plot_dict['['+str(j[1])+','+str(j[0])+']'] = plt.subplot(self.plot_grid[i/plot_per_line,i%plot_per_line])
            # create the figures:
            self.figure_2D( subplot=self.plot_dict['['+str(j[0])+','+str(j[1])+']'], param1=j[0], param2=j[1], names=names_temp, **kwargs )
        # do the legend:
        self.set_legend( names=names_temp, **kwargs )
        # do the title:
        self.set_title( title=title, **kwargs)
        # set the global appearence of the plot:
        self.set_plot_dimensions( num_col=num_col, num_rows=num_row, **kwargs )
        
    # -----------------------------------------------------------------------------------
    
    def plot3D(self, params=None, names=None, **kwargs):
        """
        Main 3D plotting function. This takes a set of parameter names triplets, a set of fisher
        names and produces a figure with 3D plots.
        Uses settings in plot_settings but the resulting figure (figure) can be edited
        just like all matplotlib.pyplot figures. Optionally a title can be supplied and
        any kind of settings can be passed by kwargs.
        
        :param params: a triplet of names of a parameter or a list of parameter names triplets. If None does all.
        :type params: :class:`string` or a :class:`list` of :class:`string`
        :param names: a name of a fisher matrix or a list of names of fisher matrices. If None does all.
        :type names: :class:`string` or a :class:`list` of :class:`string`
        :param title: the optional title for the plot.
        :type title: :class:`string`
        :param kwargs: optional keyword settings.
        
        """
        raise ValueError('plot3D: not yet implemented.')   
        
    # -----------------------------------------------------------------------------------
    
    def plot_tri(self, params=None, names=None, title=None, **kwargs):
        """
        Main triangle plotting function. This takes a set of parameter names, a set of fisher
        names and produces a figure with triangular plots.
        Uses settings in plot_settings but the resulting figure (figure) can be edited
        just like all matplotlib.pyplot figures. Optionally a title can be supplied and
        any kind of settings can be passed by kwargs.
        
        :param params: a list of names of a parameter. If None does all.
        :type params: :class:`string` or a :class:`list` of :class:`string`
        :param names: a name of a fisher matrix or a list of names of fisher matrices. If None does all.
        :type names: :class:`string` or a :class:`list` of :class:`string`
        :param title: the optional title for the plot.
        :type title: :class:`string`
        :param kwargs: optional keyword settings.
        
        """
        # get params:
        all_names  = copy.deepcopy( self.plot_fishers.fisher_name_list )
        all_params = self.plot_fishers.get_parameter_list()
        if params==None and names==None:
            names_temp  = copy.deepcopy( all_names )
            params_temp = copy.deepcopy( all_params )
        elif params==None and names!=None:
            names_temp  = [ i for i in fu.make_list(names) if i in all_names ]
            params_temp = self.plot_fishers.get_parameter_list( names=names_temp )
        elif params!=None and names==None:
            params_temp = [ i for i in fu.make_list(params) if i in all_params ]
            names_temp  = copy.deepcopy( all_names )
        else:
            names_temp  = [ i for i in fu.make_list(names) if i in all_names ]
            params_temp = [ i for i in fu.make_list(params) if i in self.plot_fishers.get_parameter_list(names_temp) ]
        # check input:
        if len(names_temp) == 0:
            raise ValueError('ERROR PLOT tri: no valid input fisher matrix.\nPossible values are: '+str(all_names))
        if len(params_temp) == 0:
            raise ValueError('ERROR PLOT tri: no valid input parameters.\nPossible values are: '+str(all_params))
        # override settings:
        
        # override name bind settings:
        self.bind_plot_settings_to_names( names=names_temp, **kwargs )
        # create the layout of the figure:
        self.plot_grid = gridspec.GridSpec( len(params_temp), len(params_temp) )
        self.plot_number = len(params_temp)*(len(params_temp)-1)/2
        # fill the single plotters:
        self.plot_dict = {}
        for i,name1 in enumerate(params_temp):
            for j,name2 in enumerate(params_temp):
                # get only one diagonal:
                if i<j: continue
                # generate the dictionary: we want symmetric keys just to be sure.
                self.plot_dict['['+str(name1)+','+str(name2)+']'] = plt.subplot(self.plot_grid[i,j])
                self.plot_dict['['+str(name2)+','+str(name1)+']'] = plt.subplot(self.plot_grid[i,j])
                # create the figures:
                if i==j:
                    # setting up the label exclusion for the triangular:
                    if i!=0:
                        # y label on the other side:
                        kwargs['D1_ylabel_right'] = True
                    else:   
                        kwargs['D1_ylabel_right'] = False
                        
                    if i!=len(params_temp)-1:
                        kwargs['D1_show_x_ticks_labels'] = False
                        kwargs['D1_show_xaxis_label']    = False
                    else:   
                        kwargs['D1_show_x_ticks_labels'] = True
                        kwargs['D1_show_xaxis_label']    = True
                        
                    self.figure_1D( subplot=self.plot_dict['['+str(name1)+','+str(name2)+']'], param=name1, names=names_temp, **kwargs)
                else:
                    # setting up the label exclusion for the triangular:
                    if j!=0:
                        kwargs['D2_show_y_ticks_labels'] = False
                        kwargs['D2_show_yaxis_label']    = False
                    else:   
                        kwargs['D2_show_y_ticks_labels'] = True
                        kwargs['D2_show_yaxis_label']    = True
                        
                    if i!=len(params_temp)-1:
                        kwargs['D2_show_x_ticks_labels'] = False
                        kwargs['D2_show_xaxis_label']    = False
                    else:   
                        kwargs['D2_show_x_ticks_labels'] = True
                        kwargs['D2_show_xaxis_label']    = True
                    # doing the plot:
                    self.figure_2D( subplot=self.plot_dict['['+str(name1)+','+str(name2)+']'], param1=name2, param2=name1, names=names_temp, **kwargs)
        # do the legend:
        self.set_legend( names=names_temp, **kwargs )
        # do the title:
        self.set_title( title=title, **kwargs)
        # set the global appearence of the plot:
        self.set_triplot_dimensions( num_col=len(params_temp), num_rows=len(params_temp), **kwargs )
        
    # -----------------------------------------------------------------------------------
    
    def plot_mixed(self, params=None, names=None, **kwargs):
        """
        Does mixed plot. This takes a set of parameter names, a set of fisher
        names and produces a figure with combined plots.
        Uses settings in plot_settings but the resulting figure (figure) can be edited
        just like all matplotlib.pyplot figures. Optionally a title can be supplied and
        any kind of settings can be passed by kwargs.
        
        :param params: a name of a parameter or a list of parameter names. If None does nothing.
        :type params: :class:`string` or a :class:`list` of :class:`string`
        :param names: a name of a fisher matrix or a list of names of fisher matrices. If None does all.
        :type names: :class:`string` or a :class:`list` of :class:`string`
        :param title: the optional title for the plot.
        :type title: :class:`string`
        :param kwargs: optional keyword settings.
        
        """
        raise ValueError('plot_mixed: not yet implemented.')     
    
    # -----------------------------------------------------------------------------------
    
    def figure_1D( self, subplot, param, names, **kwargs ):
        """
        Creates a 1D gaussian plot of one parameter and with all the Fisher in names.
        Notice this is NOT one of the main plotting functions but rather used by other
        functions. 
        
        :param subplot: a subplot. The plot will be done here.
        :type subplot: :class:`matplotlib.axes._subplots.AxesSubplot`
        :param param: the name of the parameter to plot.
        :type param: :class:`string`
        :param names: a name of a Fisher matrix or a list of names of Fisher matrices. If None does none.
        :type names: :class:`string` or a :class:`list` of :class:`string`
        :param kwargs: optional keyword settings.    
        
        """
        # check the input:
        if param not in self.plot_fishers.get_parameter_list():
            raise ValueError('The input parameter '+str(param)+' is not in the parameters of the list.')
        names_temp  = copy.deepcopy( fu.make_list(names)  )
        for name in names_temp:
            if name not in self.plot_fishers.fisher_name_list:
                raise ValueError('The input names are not in the list of fisher names')
        # check the settings:
        confidence_level = self.setting_setter('D1_confidence_level', **kwargs)
        num_points = self.setting_setter('D1_num_points', **kwargs)
        linewidth  = self.setting_setter('D1_line_thickness', **kwargs)
        normalized = self.setting_setter('D1_norm_prob', **kwargs)
        nice       = self.setting_setter('D1_use_nice_numbers', **kwargs)
        show_x_ticks = self.setting_setter('D1_show_x_ticks', **kwargs)
        number_x_ticks = self.setting_setter('D1_number_x_ticks', **kwargs)
        show_y_ticks = self.setting_setter('D1_show_y_ticks', **kwargs)
        number_y_ticks = self.setting_setter('D1_number_y_ticks', **kwargs)
        show_x_ticks_labels = self.setting_setter('D1_show_x_ticks_labels', **kwargs)
        show_y_ticks_labels = self.setting_setter('D1_show_y_ticks_labels', **kwargs)
        show_xaxis_label = self.setting_setter('D1_show_xaxis_label', **kwargs)
        show_yaxis_label = self.setting_setter('D1_show_yaxis_label', **kwargs)
        x_label_rotation = self.setting_setter('D1_x_label_rotation', **kwargs)
        y_label_rotation = self.setting_setter('D1_y_label_rotation', **kwargs)
        prob_label = self.setting_setter('D1_prob_label', **kwargs)
        norm_prob_label = self.setting_setter('D1_norm_prob_label', **kwargs)
        main_fontsize = self.setting_setter('D1_main_fontsize', **kwargs)
        secondary_fontsize = self.setting_setter('D1_secondary_fontsize', **kwargs)
        filled = self.setting_setter('D1_filled', **kwargs)
        filling_alpha = self.setting_setter('D1_alpha', **kwargs)
        ylabel_right = self.setting_setter('D1_ylabel_right', **kwargs)
        xlabel_up = self.setting_setter('D1_xlabel_up', **kwargs)
        # do the plotting:
        if filled:
            for i,confl in enumerate(confidence_level):
                # get the data:
                data = self.plot_fishers.compute_gaussian( params=param, confidence_level=confl, names=names_temp, num_points=num_points/len(confidence_level), normalized=normalized, nice_bounds=False )
                # do the plotting:
                for i,name in enumerate(names):
                    subplot.fill_between( data[param][name][0], 0.0, data[param][name][1],
                                          color     = self.bind_solid_colors[name],
                                          linewidth = 0.0,
                                          alpha     = filling_alpha )
        # get the data:
        data = self.plot_fishers.compute_gaussian( params=param, confidence_level=max(confidence_level), names=names_temp, num_points=num_points, normalized=normalized, nice_bounds=nice )
        # do the plotting:
        for i,name in enumerate(names):
            subplot.plot( data[param][name][0], data[param][name][1], 
                          color     = self.bind_line_colors[name], 
                          linestyle = self.bind_linestyle[name],
                          linewidth = linewidth, 
                          label     = self.bind_labels[name] )
        # set the plot range:
        x_limits = self.plot_fishers.compute_plot_range(params=param, confidence_level=max(confidence_level), names=names_temp, nice=nice)[param]
        subplot.set_xlim( x_limits )
        subplot.set_ylim( [0.0, 1.05] )
        # set the ticks and ticks labels:
        if show_x_ticks:
            subplot.set_xticks( np.linspace(x_limits[0], x_limits[1], number_x_ticks)  )
        else:
            subplot.set_xticks([])
        if show_y_ticks:
            subplot.set_yticks( np.linspace(0,1, number_y_ticks)  )
        else:
            subplot.set_yticks([])
        # xticks:
        if not show_x_ticks_labels:
            subplot.set_xticklabels( [] )
        else:
            xticks = np.linspace(x_limits[0], x_limits[1], number_x_ticks)
            xticks = [ u'$'+str(i)+'$' for i in xticks ]
            subplot.set_xticklabels( xticks, fontsize=secondary_fontsize )
            if xlabel_up: subplot.xaxis.tick_top()
        # yticks:
        if not show_y_ticks_labels:
            subplot.set_yticklabels( [] )
        else:
            yticks = np.linspace(0,1, number_y_ticks) 
            yticks = [ u'$'+str(i)+'$' for i in yticks ]
            subplot.set_yticklabels( yticks, fontsize=secondary_fontsize )
            if ylabel_right: subplot.yaxis.tick_right()
        # align the left and right tick labels:
        try:
            subplot.xaxis.get_majorticklabels()[0].set_horizontalalignment('left')
            subplot.xaxis.get_majorticklabels()[-1].set_horizontalalignment('right')
            subplot.yaxis.get_majorticklabels()[0].set_verticalalignment('bottom')
            subplot.yaxis.get_majorticklabels()[-1].set_verticalalignment('top')
        except:
            pass
        # set axis labels:
        if show_xaxis_label:
            subplot.set_xlabel( u'$'+self.plot_fishers.get_parameter_latex_names()[param]+'$', fontsize=main_fontsize, rotation=x_label_rotation )
            if xlabel_up: subplot.xaxis.set_label_position("top")
        if show_yaxis_label:
            if normalized:
                subplot.set_ylabel( norm_prob_label, fontsize=main_fontsize, rotation=y_label_rotation )
            else:
                subplot.set_ylabel( prob_label, fontsize=main_fontsize, rotation=y_label_rotation )
            if ylabel_right: subplot.yaxis.set_label_position("right")
            
    # -----------------------------------------------------------------------------------
    
    def figure_2D( self, subplot, param1, param2, names, **kwargs ):
        """
        Creates a 2D gaussian ellipse plot of param1 and param2 and with all the Fisher in names.
        Notice this is NOT one of the main plotting functions but rather used by other
        functions. 
        
        :param subplot: a subplot. The plot will be done here.
        :type subplot: :class:`matplotlib.axes._subplots.AxesSubplot`
        :param param1: the name of the first parameter of the plot.
        :type param1: :class:`string`
        :param param2: the name of the second parameter of the plot.
        :type param2: :class:`string`
        :param names: a name of a Fisher matrix or a list of names of Fisher matrices. If None does none.
        :type names: :class:`string` or a :class:`list` of :class:`string`
        :param kwargs: optional keyword settings.    
        
        """
        # check the input:
        if param1 not in self.plot_fishers.get_parameter_list():
            raise ValueError('The input parameter '+str(param1)+' is not in the parameters of the list.')
        if param2 not in self.plot_fishers.get_parameter_list():
            raise ValueError('The input parameter '+str(param2)+' is not in the parameters of the list.')
        names_temp  = copy.deepcopy( fu.make_list(names)  )
        for name in names_temp:
            if name not in self.plot_fishers.fisher_name_list:
                raise ValueError('The input names are not in the list of fisher names')
        # check the settings:
        confidence_level = self.setting_setter('D2_confidence_levels', **kwargs)
        num_points = self.setting_setter('D2_num_points', **kwargs)
        linewidth  = self.setting_setter('D2_line_thickness', **kwargs)
        nice       = self.setting_setter('D2_use_nice_numbers', **kwargs)
        show_x_ticks = self.setting_setter('D2_show_x_ticks', **kwargs)
        number_x_ticks = self.setting_setter('D2_number_x_ticks', **kwargs)
        show_y_ticks = self.setting_setter('D2_show_y_ticks', **kwargs)
        number_y_ticks = self.setting_setter('D2_number_y_ticks', **kwargs)
        show_x_ticks_labels = self.setting_setter('D2_show_x_ticks_labels', **kwargs)
        show_y_ticks_labels = self.setting_setter('D2_show_y_ticks_labels', **kwargs)
        show_xaxis_label = self.setting_setter('D2_show_xaxis_label', **kwargs)
        show_yaxis_label = self.setting_setter('D2_show_yaxis_label', **kwargs)
        x_label_rotation = self.setting_setter('D2_x_label_rotation', **kwargs)
        y_label_rotation = self.setting_setter('D2_y_label_rotation', **kwargs)
        main_fontsize = self.setting_setter('D2_main_fontsize', **kwargs)
        secondary_fontsize = self.setting_setter('D2_secondary_fontsize', **kwargs)
        filled = self.setting_setter('D2_filled', **kwargs)
        filling_alpha = self.setting_setter('D2_alphas', **kwargs)
        ylabel_right = self.setting_setter('D2_ylabel_right', **kwargs)
        xlabel_up = self.setting_setter('D2_xlabel_up', **kwargs)
        # get the data:
        for j,confidence in enumerate(confidence_level):
            # get the data from the fisher_plot_analysis:
            data = self.plot_fishers.compute_ellipse( params1=param1, params2=param2, confidence_level=confidence, names=names_temp, num_points=num_points )
            # get the area of the ellipses:
            areas = [ math.pi*data[param1][param2][name][2][2]*data[param1][param2][name][2][3] for name in names_temp ] 
            # reorder the names to plot big and then small:
            names_reorder = [ y for (x,y) in sorted( zip( areas, names_temp), reverse=True) ]
            for i,name in enumerate(names_reorder):
                # filled contour:
                if filled:
                    subplot.fill(data[param1][param2][name][0], data[param1][param2][name][1],
                                 alpha=filling_alpha[j], 
                                 facecolor=self.bind_solid_colors[name],
                                 linewidth=0.0, 
                                 )
                # contours border:
                subplot.plot( data[param1][param2][name][0], data[param1][param2][name][1],
                              color     = self.bind_line_colors[name], 
                              linestyle = self.bind_linestyle[name],
                              linewidth = linewidth, 
                              label     = self.bind_labels[name] )
                
        # set the plot range:
        ranges = self.plot_fishers.compute_plot_range(params=[param1,param2], confidence_level=max(confidence_level), names=names_temp, nice=nice)
        subplot.set_xlim( ranges[param1] )
        subplot.set_ylim( ranges[param2] )
        # set the ticks and ticks labels:
        if show_x_ticks:
            subplot.set_xticks( np.linspace(ranges[param1][0], ranges[param1][1], number_x_ticks)  )
        else:
            subplot.set_xticks([])
        if show_y_ticks:
            subplot.set_yticks( np.linspace( ranges[param2][0], ranges[param2][1], number_y_ticks)  )
        else:
            subplot.set_yticks([])
        # xticks:
        if not show_x_ticks_labels:
            subplot.set_xticklabels( [] )
        else:
            xticks = np.linspace(ranges[param1][0], ranges[param1][1], number_x_ticks)
            xticks = [ u'$'+str(i)+'$' for i in xticks ]
            subplot.set_xticklabels( xticks, fontsize=secondary_fontsize )
            if xlabel_up: subplot.xaxis.tick_top()
        # yticks:
        if not show_y_ticks_labels:
            subplot.set_yticklabels( [] )
        else:
            yticks = np.linspace( ranges[param2][0], ranges[param2][1], number_y_ticks)
            yticks = [ u'$'+str(i)+'$' for i in yticks ]
            subplot.set_yticklabels( yticks, fontsize=secondary_fontsize )
            if ylabel_right: subplot.yaxis.tick_right()
        # align the left and right tick labels:
        try:
            subplot.xaxis.get_majorticklabels()[0].set_horizontalalignment('left')
            subplot.xaxis.get_majorticklabels()[-1].set_horizontalalignment('right')
            subplot.yaxis.get_majorticklabels()[0].set_verticalalignment('bottom')
            subplot.yaxis.get_majorticklabels()[-1].set_verticalalignment('top')
        except:
            pass
        # set axis labels:
        if show_xaxis_label:
            subplot.set_xlabel( u'$'+self.plot_fishers.get_parameter_latex_names()[param1]+'$', fontsize=main_fontsize, rotation=x_label_rotation )
            if xlabel_up: subplot.xaxis.set_label_position("top")
        if show_yaxis_label:
            subplot.set_ylabel( u'$'+self.plot_fishers.get_parameter_latex_names()[param2]+'$', fontsize=main_fontsize, rotation=y_label_rotation )
            if ylabel_right: subplot.yaxis.set_label_position("right")    
                
    # -----------------------------------------------------------------------------------
    
    def setting_setter(self, key, **kwargs):
        """
        Small utility to check wether a specific setting is in plot_settings or is passed
        as a keyword argument (in kwargs).

        :param key: input setting keyword
        :type key: :class:`string`
        :param kwargs: keyword arguments to check
        :returns: the value of the setting. If the setting is present in kwargs then the value in there
            otherwise the value in :class:`cosmicfish_pylib.fisher_plot.CosmicFishPlotter.plot_settings`
        
        """
        if kwargs.has_key(key): 
            return kwargs[key]
        else:
            try:
                return getattr(self.plot_settings, key) 
            except:
                raise ValueError('Error in setting_setter: '+key+' is nor in **kwargs nor in plot_setting')
    
    # -----------------------------------------------------------------------------------
    
    def bind_plot_settings_to_names( self, names=None, **kwargs ):
        """
        This function binds settings to names using the default if something else is not 
        supplied as a keyword argument.
        
        :param names: a name of a fisher matrix or a list of names of fisher matrices to bind to settings. If None does all.
        :type names: :class:`string` or a :class:`list` of :class:`string`
        :param kwargs: optional keyword settings to bind to names.  
        
        """
        # get the names:
        all_names  = copy.deepcopy( self.plot_fishers.fisher_name_list )
        if names is not None:
            names_temp  = [ i for i in fu.make_list(names) if i in all_names ]
        else:
            names_temp  = all_names
        # let's start with colors. We want them binded to names so there is consistency of colors through plots.
        # line colors
        self.bind_line_colors = {}
        if kwargs.has_key('line_colors'): 
            for i, name in enumerate(names_temp):
                self.bind_line_colors[name] = kwargs['line_colors'][i]
        else:
            for i, name in enumerate(all_names):
                self.bind_line_colors[name] = self.plot_settings.line_colors[i]
        # solid colors
        self.bind_solid_colors = {}
        if kwargs.has_key('solid_colors'): 
            for i, name in enumerate(names_temp):
                self.bind_solid_colors[name] = kwargs['solid_colors'][i]
        else:
            for i, name in enumerate(all_names):
                self.bind_solid_colors[name] = self.plot_settings.solid_colors[i]
        # labels:
        self.bind_labels = {}
        if kwargs.has_key('labels'): 
            for i, name in enumerate(names_temp):
                self.bind_labels[name] = kwargs['labels'][i]
        else:
            for i, name in enumerate(all_names):
                self.bind_labels[name] = name
        # linestyle:
        self.bind_linestyle = {}
        if kwargs.has_key('linestyle'): 
            for i, name in enumerate(names_temp):
                self.bind_linestyle[name] = kwargs['linestyle'][i]
        else:
            for i, name in enumerate(all_names):
                self.bind_linestyle[name] = self.plot_settings.linestyle[i]
    
    # -----------------------------------------------------------------------------------
    
    def get_dimensions_plot_obj(self):
        """
        This function computes the sizes in inches of all the objects in the figure to 
        allow the computation of the figure spacings and disposition.
        Considers only the objects already present in the figure.
        
        :returns: a dictionary with all the spacings.
        :rtype: :class:`dict`
        
        """
        # initialize the dictionary of sizes:
        sizes_dictionary = {}
        # draw the canvas with the default renderer. In this way we can read the rendered dimensions of the objects.
        self.figure.canvas.draw()
        # get the default renderer:
        renderer = matplotlib.backend_bases.RendererBase()
        # default dpi is 72.0:
        default_dpi = 72.0
        # get the size of the figure:
        figure_x_size = self.figure.get_size_inches()[0] #: in inches
        figure_y_size = self.figure.get_size_inches()[1] #: in inches
        sizes_dictionary['figure'] = self.figure.get_size_inches() #: in inches
        # get max (x,y) size of the y ticks:
        max_y_tick_label = [0.0,0.0]
        for plot in self.plot_dict.values():
            for ylabel in plot.get_yticklabels():
                x_dimension = ylabel.get_window_extent(renderer).width/default_dpi  #: in inches
                y_dimension = ylabel.get_window_extent(renderer).height/default_dpi  #: in inches
                if x_dimension>max_y_tick_label[0]: max_y_tick_label[0]=x_dimension
                if y_dimension>max_y_tick_label[1]: max_y_tick_label[1]=y_dimension
        sizes_dictionary['y_tick_labels'] = max_y_tick_label #: in inches
        # get max (x,y) size of the x ticks:
        max_x_tick_label = [0.0,0.0]
        for plot in self.plot_dict.values():
            for xlabel in plot.get_xticklabels():
                x_dimension = xlabel.get_window_extent(renderer).width/default_dpi  #: in inches
                y_dimension = xlabel.get_window_extent(renderer).height/default_dpi  #: in inches
                if x_dimension>max_x_tick_label[0]: max_x_tick_label[0]=x_dimension
                if y_dimension>max_x_tick_label[1]: max_x_tick_label[1]=y_dimension
        sizes_dictionary['x_tick_labels'] = max_x_tick_label #: in inches
        # get max (x,y) size of the y label:
        max_y_label = [0.0,0.0]
        for plot in self.plot_dict.values():    
            x_dimension = plot.yaxis.get_label().get_window_extent(renderer).width/default_dpi  #: in inches
            y_dimension = plot.yaxis.get_label().get_window_extent(renderer).height/default_dpi  #: in inches
            if x_dimension>max_y_label[0]: max_y_label[0]=x_dimension
            if y_dimension>max_y_label[1]: max_y_label[1]=y_dimension
        sizes_dictionary['y_labels'] = max_y_label #: in inches 
        # get max (x,y) size of the x label:
        max_x_label = [0.0,0.0]
        for plot in self.plot_dict.values():     
            x_dimension = plot.xaxis.get_label().get_window_extent(renderer).width/default_dpi  #: in inches
            y_dimension = plot.xaxis.get_label().get_window_extent(renderer).height/default_dpi  #: in inches
            if x_dimension>max_x_label[0]: max_x_label[0]=x_dimension
            if y_dimension>max_x_label[1]: max_x_label[1]=y_dimension
        sizes_dictionary['x_labels'] = max_x_label #: in inches 
        # get the padding:
        sizes_dictionary['xtick_pad'] = max( matplotlib.rcParams['xtick.major.pad'], matplotlib.rcParams['xtick.minor.pad'] )/default_dpi #: in inches
        sizes_dictionary['ytick_pad'] = max( matplotlib.rcParams['ytick.major.pad'], matplotlib.rcParams['ytick.minor.pad'] )/default_dpi #: in inches
        sizes_dictionary['label_pad'] = matplotlib.rcParams['axes.labelpad']/default_dpi #: in inches
        # get the size of the title:
        if self.title is not None:
            sizes_dictionary['title'] = [ self.title.get_window_extent(renderer).width/default_dpi,
                                          self.title.get_window_extent(renderer).height/default_dpi ]
        else:
            sizes_dictionary['title'] = [0.0, 0.0]
        # get the size of the legend:
        if self.legend is not None:
            sizes_dictionary['legend'] = [ self.legend.get_window_extent(renderer).width/default_dpi,
                                          self.legend.get_window_extent(renderer).height/default_dpi ]
        else:
            sizes_dictionary['legend'] = [0.0, 0.0]
        # return the dictionary:
        return sizes_dictionary
                
    # -----------------------------------------------------------------------------------
    
    def set_plot_dimensions( self, num_col, num_rows, **kwargs ):
        """
        Sets the dimensions of the plot based on the settings passed.
        
        :param num_col: number of colums in the plot grid.
        :type num_col: :class:`int`
        :param num_rows: number of rows in the plot grid.
        :type num_rows: :class:`int`
        :param kwargs: optional keyword settings
        
        """
        # override settings:
        tight_layout            = self.setting_setter('tight_layout', **kwargs)
        x_size_in               = self.setting_setter('figure_width', **kwargs)
        y_size_in               = self.setting_setter('figure_height', **kwargs)
        subplot_x_size          = self.setting_setter('subplot_x_size', **kwargs)
        subplot_y_size          = self.setting_setter('subplot_y_size', **kwargs)
        use_fixed_figure_width  = self.setting_setter('use_fixed_figure_width', **kwargs)
        use_fixed_figure_height = self.setting_setter('use_fixed_figure_height', **kwargs)
        legend_takes_place_plot = self.setting_setter('legend_takes_place_plot', **kwargs)
        legend_loc = self.setting_setter('legend_loc', **kwargs)
        # do the magic: first get the sizes of everything
        dimensions  = self.get_dimensions_plot_obj()
        small_space = 1.5*points_to_mm/10.0/inch_to_cm
        # compute x_size and y_size:
        x_size = num_col*( subplot_x_size*cm_to_inch +dimensions['ytick_pad'] +2.0*dimensions['label_pad'] +dimensions['y_tick_labels'][0] +dimensions['y_labels'][0] )
        y_size = num_rows*( subplot_y_size*cm_to_inch +dimensions['xtick_pad'] +2.0*dimensions['label_pad'] +dimensions['x_tick_labels'][1] +dimensions['x_labels'][1] )
        # add the title:
        if self.title is not None:
            x_size = max( x_size, dimensions['title'][0] )
            y_size = y_size +dimensions['title'][1] +dimensions['label_pad']
        # add the legend:
        if legend_takes_place_plot or self.legend is None:
            pass # no need to add. The legend is inside a plot.
        else:
            legend_code = self.legend.codes[legend_loc]
            # add space on the x dimension
            if legend_code in [1,2,3,4,5,6,7]:
                x_size = x_size +dimensions['legend'][0]
            # add space on the y dimension
            if legend_code in [1,2,3,4,8,9]:
                y_size = y_size +dimensions['legend'][1]
            # get extra size if the legend is bigger than the plot:
            if legend_code in [8,9,10]:
                x_size = max(x_size, dimensions['legend'][0])
            if legend_code in [5,6,7,10]:
                y_size = max(y_size, dimensions['legend'][1])
        # now fix the x_size and y_size of the plot:
        if use_fixed_figure_width and use_fixed_figure_height:
            # both sizes fixed. Do nothing, just copy the setting.
            x_size = x_size_in/inch_to_cm
            y_size = y_size_in/inch_to_cm
        elif use_fixed_figure_width and not use_fixed_figure_height:
            x_size = x_size_in/inch_to_cm
        elif not use_fixed_figure_width and use_fixed_figure_height:
            y_size = y_size_in/inch_to_cm
        elif not use_fixed_figure_width and not use_fixed_figure_height:
            pass
        # set figure size:
        self.figure.set_size_inches( x_size, y_size )
        # compute the global paddings, and the space between sub plots
        wspace = ( +dimensions['ytick_pad'] +2.0*dimensions['label_pad'] +dimensions['y_tick_labels'][0] +dimensions['y_labels'][0] )/2.0
        hspace = ( +dimensions['xtick_pad'] +2.0*dimensions['label_pad'] +dimensions['x_tick_labels'][1] +dimensions['x_labels'][1] )/2.0
        if 'left' in [ subplot.yaxis.get_label_position() for subplot in self.plot_dict.values() ]:
            left   = 2.0*wspace
        else:
            left   = small_space
        if 'right' in [ subplot.yaxis.get_label_position() for subplot in self.plot_dict.values() ]:
            right  = x_size - 2.0*wspace
        else:
            right  = x_size - small_space
        if 'bottom' in [ subplot.xaxis.get_label_position() for subplot in self.plot_dict.values() ]:
            bottom  = 2.0*hspace
        else:
            bottom  = small_space
        if 'top' in [ subplot.xaxis.get_label_position() for subplot in self.plot_dict.values() ]:
            top  = y_size - 2.0*hspace
        else:
            top  = y_size - small_space
        # add the title:
        if self.title is not None: top = top -dimensions['title'][1] -dimensions['label_pad']
        # add the legend:
        if legend_takes_place_plot or self.legend is None:
            pass # no need to add. The legend is inside a plot.
        else:
            legend_code = self.legend.codes[legend_loc]
            # add space on the x dimension
            if legend_code in [1,4,5,7]:
                right = right -dimensions['legend'][0]
            if legend_code in [2,3,6]:
                left  = left  +dimensions['legend'][0]
            if legend_code in [1,2,9]:
                top = top -dimensions['legend'][1]
            if legend_code in [3,4,8]:
                bottom = bottom +dimensions['legend'][1]
        # rearrange the plot grid:
        self.plot_grid.update( left   = left/x_size,
                               right  = right/x_size, 
                               bottom = bottom/y_size,
                               top    = top/y_size,
                               wspace = wspace, 
                               hspace = hspace )

        # position the legend:
        if self.legend is not None:
            if legend_takes_place_plot:
                num_plots = self.plot_number
                plot_per_line = self.plot_grid.get_geometry()[1]
                bbox     = self.plot_grid[(num_plots-1)/plot_per_line,(num_plots-1)%plot_per_line].get_position(self.figure)
                self.legend.set_bbox_to_anchor( bbox=bbox )
            
        # position the title:
        if self.title is not None:
            self.title.set_position((0.5, 1.0-dimensions['ytick_pad']/y_size))
        
        # use tightplot if wanted: 
        if tight_layout:
            self.plot_grid.tight_layout( self.figure )
        
        # draw the canvas, just to be sure...
        self.figure.canvas.draw()
        
    # -----------------------------------------------------------------------------------
    
    def set_triplot_dimensions( self, num_col, num_rows, **kwargs ):
        """
        Sets the dimensions of the triangular plot based on the settings passed.
        This is different from the other dimension setter because triplots need some
        personalization.
        
        :param num_col: number of colums in the plot grid.
        :type num_col: :class:`int`
        :param num_rows: number of rows in the plot grid.
        :type num_rows: :class:`int`
        :param kwargs: optional keyword settings
        
        """
        # override settings:
        tight_layout            = self.setting_setter('tight_layout', **kwargs)
        x_size_in               = self.setting_setter('figure_width', **kwargs)
        y_size_in               = self.setting_setter('figure_height', **kwargs)
        subplot_x_size          = self.setting_setter('subplot_x_size', **kwargs)
        subplot_y_size          = self.setting_setter('subplot_y_size', **kwargs)
        use_fixed_figure_width  = self.setting_setter('use_fixed_figure_width', **kwargs)
        use_fixed_figure_height = self.setting_setter('use_fixed_figure_height', **kwargs)
        legend_takes_place_plot = self.setting_setter('legend_takes_place_plot_tri', **kwargs)
        legend_loc = self.setting_setter('legend_loc', **kwargs)
        # do the magic: first get the sizes of everything
        dimensions  = self.get_dimensions_plot_obj()
        small_space = 1.5*points_to_mm/10.0/inch_to_cm
        # spacing between sub plots:
        tri_subspace = 0.1 # cm
        # compute x_size and y_size:
        x_size = num_col*subplot_x_size/inch_to_cm +2.0*( dimensions['ytick_pad'] +2.0*dimensions['label_pad'] +dimensions['y_tick_labels'][0] +dimensions['y_labels'][0] ) +(num_col-1.0)*tri_subspace/inch_to_cm
        y_size = num_rows*subplot_y_size/inch_to_cm +1.0*( dimensions['xtick_pad'] +2.0*dimensions['label_pad'] +dimensions['x_tick_labels'][1] +dimensions['x_labels'][1] ) +(num_rows-1.0)*tri_subspace/inch_to_cm
        # add the title:
        if self.title is not None:
            x_size = max( x_size, dimensions['title'][0] )
            y_size = y_size +dimensions['title'][1] +dimensions['label_pad']
        # add the legend:
        if legend_takes_place_plot or self.legend is None: pass # no need to add. The legend is inside a plot.
        else:
            legend_code = self.legend.codes[legend_loc]
            # add space on the x dimension
            if legend_code in [1,2,3,4,5,6,7]:
                x_size = x_size +dimensions['legend'][0]
            # add space on the y dimension
            if legend_code in [1,2,3,4,8,9]:
                y_size = y_size +dimensions['legend'][1]
        # now fix the x_size and y_size of the plot:
        if use_fixed_figure_width and use_fixed_figure_height:
            # both sizes fixed. Do nothing, just copy the setting.
            x_size = x_size_in/inch_to_cm
            y_size = y_size_in/inch_to_cm
        elif use_fixed_figure_width and not use_fixed_figure_height:
            x_size = x_size_in/inch_to_cm
        elif not use_fixed_figure_width and use_fixed_figure_height:
            y_size = y_size_in/inch_to_cm
        elif not use_fixed_figure_width and not use_fixed_figure_height:
            pass
        # set figure size:
        self.figure.set_size_inches( x_size, y_size )
        # compute the global paddings, and the space between sub plots
        wspace = ( +dimensions['ytick_pad'] +2.0*dimensions['label_pad'] +dimensions['y_tick_labels'][0] +dimensions['y_labels'][0] )/2.0
        hspace = ( +dimensions['xtick_pad'] +2.0*dimensions['label_pad'] +dimensions['x_tick_labels'][1] +dimensions['x_labels'][1] )/2.0
        if 'left' in [ subplot.yaxis.get_label_position() for subplot in self.plot_dict.values() ]:
            left   = 2.0*wspace
        else:
            left   = small_space
        if 'right' in [ subplot.yaxis.get_label_position() for subplot in self.plot_dict.values() ]:
            right  = x_size - 2.0*wspace
        else:
            right  = x_size - small_space
        if 'bottom' in [ subplot.xaxis.get_label_position() for subplot in self.plot_dict.values() ]:
            bottom  = 2.0*hspace
        else:
            bottom  = small_space
        if 'top' in [ subplot.xaxis.get_label_position() for subplot in self.plot_dict.values() ]:
            top  = y_size - 2.0*hspace
        else:
            top  = y_size - small_space
        # add the title:
        if self.title is not None: top = top -dimensions['title'][1] -dimensions['label_pad']
        # add the legend:
        if not legend_takes_place_plot and self.legend is not None:
            legend_code = self.legend.codes[legend_loc]
            # add space on the x dimension
            if legend_code in [1,4,5,7]:
                right = right -dimensions['legend'][0]
            if legend_code in [2,3,6]:
                left  = left  +dimensions['legend'][0]
            if legend_code in [1,2,9]:
                top = top -dimensions['legend'][1]
            if legend_code in [3,4,8]:
                bottom = bottom +dimensions['legend'][1]
        # rearrange the plot grid:
        self.plot_grid.update( left   = left/x_size,
                               right  = right/x_size, 
                               bottom = bottom/y_size,
                               top    = top/y_size,
                               wspace = tri_subspace*cm_to_inch, 
                               hspace = tri_subspace*cm_to_inch)
        
        # position the legend:
        if self.legend is not None:
            if legend_takes_place_plot:
                plot_per_line = self.plot_grid.get_geometry()[1]
                bbox     = self.plot_grid[0,plot_per_line-1].get_position(self.figure)
                self.legend.set_bbox_to_anchor( bbox=bbox )
            
        # position the title:
        if self.title is not None:
            self.title.set_position((0.5, 1.0-dimensions['ytick_pad']/y_size))
        
        # use tightplot if wanted: 
        if tight_layout:
            self.plot_grid.tight_layout( self.figure )
        
        # draw the canvas, just to be sure...
        self.figure.canvas.draw()
        
    # -----------------------------------------------------------------------------------
    
    def set_legend( self, names=None, **kwargs ):
        """
        Creates the legend of the plot.
        
        :param names: a name of a fisher matrix or a list of names of fisher matrices. If None does all.
        :type names: :class:`string` or a :class:`list` of :class:`string`
        :param kwargs: optional keyword settings.
        
        """
        # get the names:
        all_names  = copy.deepcopy( self.plot_fishers.fisher_name_list )
        if names==None:
            names_temp  = copy.deepcopy( all_names )
        else:
            names_temp  = [ i for i in fu.make_list(names) if i in all_names ]
        # override settings:
        do_legend = self.setting_setter('do_legend', **kwargs)
        filled    = self.setting_setter('legend_filled', **kwargs)
        main_fontsize  = self.setting_setter('legend_fontsize', **kwargs)
        frameon   = self.setting_setter('legend_frame', **kwargs)
        fancybox  = self.setting_setter('legend_fancybox', **kwargs)
        shadow    = self.setting_setter('legend_shadow', **kwargs)
        markerfirst  = self.setting_setter('legend_markerfirst', **kwargs)
        ncol  = self.setting_setter('legend_ncol', **kwargs)
        legend_loc  = self.setting_setter('legend_loc', **kwargs)
        # check:
        if not do_legend: return
        # create legend handlers
        leg_handlers = []
        for i, name in enumerate(names_temp):
            if (filled):
                leg_handlers.append( mpatches.Patch(facecolor = self.bind_solid_colors[name],
                                                    edgecolor = self.bind_line_colors[name],
                                                    linestyle = self.bind_linestyle[name], ) )
            else:
                leg_handlers.append( mlines.Line2D([], [], color = self.bind_line_colors[name],
                                                   linestyle = self.bind_linestyle[name], ) )
        # process names:
        names_legend = [ u'$\\mathrm{'+str(i).replace(" ", "\ ").replace('_', '\ ')+'}$' for i in names_temp]
        self.legend = self.figure.legend( handles  = leg_handlers, 
                                          labels   = names_legend, 
                                          fontsize = main_fontsize,
                                          frameon  = frameon,
                                          fancybox = fancybox,
                                          shadow   = shadow,
                                          markerfirst = markerfirst,
                                          ncol     = ncol,
                                          borderaxespad = 0.0,
                                          loc = legend_loc,
                                          )
            
    # -----------------------------------------------------------------------------------
    
    def set_title( self, title=None, **kwargs ):
        """
        Creates the title of the plot.
        
        :param names: a name of a fisher matrix or a list of names of fisher matrices. If None does all.
        :type names: :class:`string` or a :class:`list` of :class:`string`
        :param kwargs: optional keyword settings.
        
        """
        # check input:
        if title is None: return
        # override settings:
        fontsize = self.setting_setter('title_fontsize', **kwargs)
        # create the title:
        self.title = self.figure.suptitle( u'$\\mathrm{'+str(title).replace(" ", "\ ")+'}$', fontsize=fontsize )
        
    # -----------------------------------------------------------------------------------
    
# ***************************************************************************************
