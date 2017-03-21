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
.. module:: fisher_plot_analysis
   :platform: Unix
   :synopsis: Module that contains a set of tools to perform the analysis of a set of Fisher matrices. 
   
.. moduleauthor:: Marco Raveri <mraveri@uchicago.edu> for the CosmicFish code.

"""

# ***************************************************************************************

import os
import math
import copy
import numpy as np
import utilities         as fu
import fisher_matrix     as fm
import fisher_derived    as fd
import fisher_operations as fo
from __builtin__ import dict

import __init__ as CosmicFishPyLib

# ***************************************************************************************

class CosmicFish_FisherAnalysis():
    """
    This class takes care of handeling a set of Fisher matrices with plotting in mind.
    This class is meant to hold a list of Fisher matrices and have defined on this list
    a set of vectorized operations.
    For now no caching is implemented.
    
    :ivar fisher_list: list of Fisher matrices :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`
    :ivar fisher_name_list: list of names of the Fisher matrices. The names are used as the unique 
        identifier of the Fisher matrix. No double name is enforced.
    
    """
    
    # -----------------------------------------------------------------------------------
    
    # class getters:
    
    def get_fisher_list(self):
        """ :returns: the list fisher matrices. """
        return self.fisher_list
    
    def get_fisher_name_list(self):
        """ :returns: the list fisher matrices names. These are the unique identifiers of the list."""
        return self.fisher_name_list
        
    # -----------------------------------------------------------------------------------

    def __init__( self, fisher_list=None, fisher_path=None, search_fisher_guess=False, with_derived=True ):
        """
        **CosmicFish_FisherAnalysis class constructor**. The constructor builds the Fisher matrices list.
        If all the arguments of the constructor are None then the Fisher list will be created and will be empty.
        If this is initialized with fisher_list then this will constitute the base Fisher set.
        if this is initialized with fisher_path then the code will search the desired path for Fisher matrices.
        If initialized with both fisher_list and fisher_path then both will be added.

        :param fisher_list: Fisher matrix or list of Fisher matrices.
        :type fisher_list: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix` or :class:`list` of :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`
        :param fisher_path: path or list of paths to search for Fisher matrices.
        :type fisher_path: :class:`string` or :class:`list` of :class:`string`
        :param search_fisher_guess: wether to guess the Fisher matrix name or not. 
            Guessing assumes that the Fisher matrices have 'fisher_matrix' and '.dat' in the name as happens
            with Fisher matrices produced with CosmicFish.
        :type search_fisher_guess: :class:`bool`
        :param with_derived: wether to search for derived Fisher matrices to add to the base Fisher that are found.
        :type with_derived: :class:`bool`

        """
        # create class instances:
        self.fisher_list      = []
        self.fisher_name_list = []
        # check wether we have to do something:
        if not ( fisher_list==None and fisher_path==None ):
            # do the initialization from fisher_list:
            if not fisher_list==None:
                self.add_fisher_matrix(fisher_list)
            # do the initialization from path:
            if not fisher_path==None:
                self.search_fisher_path( fisher_path=fisher_path, search_fisher_guess=search_fisher_guess, with_derived=with_derived )
        
    # -----------------------------------------------------------------------------------

    def __del__(self):
        """
        CosmicFish_FisherAnalysis class destructor. 
        Makes sure everything is gone when called.
        """
        self.fisher_list[:] = []
        self.fisher_name_list[:] = []
    
    # -----------------------------------------------------------------------------------

    def search_fisher_path(self, fisher_path, search_fisher_guess=False, with_derived=True ):
        """
        Searches a path for fisher matrices. 
        Will detect wether fisher_path contains directly the paths to the Fisher files or folder.
        If a list of folders is passed all the folders will be searched, first for Fisher matrices
        then for derived Fisher matrices.

        :param fisher_path: path or list of paths to search for Fisher matrices. If this contains Fisher matrices files those are imported as well.
        :type fisher_path: :class:`string` or :class:`list` of :class:`string`
        :param search_fisher_guess: wether to guess the Fisher matrix name or not. 
            Guessing assumes that the Fisher matrices have 'fisher_matrix' and '.dat' in the name as happens
            with Fisher matrices produced with CosmicFish.
        :type search_fisher_guess: :class:`bool`
        :param with_derived: wether to search for derived Fisher matrices to add to the base Fisher that are found.
        :type with_derived: :class:`bool`
        
        """
        fisher_path = copy.deepcopy( fu.make_list(fisher_path) )
        # explore all paths:
        fisher_files = []
        derived_fisher_files = []
        for path in fisher_path:
            if os.path.isfile(path):
                fisher_files.append(path)
                derived_fisher_files.append(path)
                continue
            if search_fisher_guess:
                # search for all files with fisher_matrix in the name and fisher_matrix_derived in the name
                # this is the standard cosmicfish output and should be safe if you are using cosmicfish.
                # get the file list:
                files = [ os.path.join(path, f) for f in os.listdir(path) if os.path.isfile( os.path.join(path, f) )]
                # filter the list by guessing the names to get Fisher matrices:
                for f in files:
                    f_temp = os.path.basename(f)
                    if 'fisher_matrix' in f_temp and '.dat' in f_temp and '_derived' not in f_temp:
                        fisher_files.append( f )
                # filter the list by guessing the names to get derived Fisher matrices:
                for f in files:
                    f_temp = os.path.basename(f)
                    if 'fisher_matrix' in f_temp and '_derived' in f_temp and '.dat' in f_temp:
                        derived_fisher_files.append( f )  
            else:
                # search for all the files and try to import them as fisher matrices then as derived fishers
                # get the file list:
                files = [ os.path.join(path, f) for f in os.listdir(path) if os.path.isfile( os.path.join(path, f) )]
                for f in files:
                    fisher_files.append( f )
                    derived_fisher_files.append( f ) 
                    
        # import fisher matrices:
        for file in fisher_files:
            if CosmicFishPyLib.__feedback__ > 1: print 'Trying to import: '+file+' as a Fisher matrix: ',
            try:    
                self.fisher_list.append( fm.fisher_matrix( file_name = file ))
                self.fisher_name_list.append( self.fisher_list[-1].name )
                if CosmicFishPyLib.__feedback__ > 1: print 'SUCCESS'
                if CosmicFishPyLib.__feedback__ == 1: print 'Imported as a Fisher matrix: '+file
                # remove from the derived fishers if necessary:
                derived_fisher_files.remove(file)
            except:
                if CosmicFishPyLib.__feedback__ > 1: print 'FAIL'
        # import derived fisher matrices:
        if with_derived:
            for file in derived_fisher_files:
                # get the derived matrix:
                if CosmicFishPyLib.__feedback__ > 1: print 'Trying to import: '+file+' as a derived Fisher matrix: ',
                try:    
                    fisher_derived = fd.fisher_derived( file_name = file )
                    if CosmicFishPyLib.__feedback__ > 1: print 'SUCCESS'
                except: 
                    fisher_derived = None
                    if CosmicFishPyLib.__feedback__ > 1: print 'FAIL'
                # add derived to the fisher matrix when possible:
                if fisher_derived is not None:
                    for i in xrange( len( self.fisher_name_list) ):
                        if CosmicFishPyLib.__feedback__ > 1: print 'Trying to add derived from: '+fisher_derived.name+' to the Fisher matrix: '+self.fisher_list[i].name,
                        try: 
                            self.fisher_list[i] = fisher_derived.add_derived( fisher_matrix=self.fisher_list[i], preserve_input=True )
                            self.fisher_name_list[i] = self.fisher_list[i].name 
                            if CosmicFishPyLib.__feedback__ > 1: print 'SUCCESS'
                            if CosmicFishPyLib.__feedback__ == 1: print 'Added derived parameters from: '+fisher_derived.name+' to Fisher matrix: '+self.fisher_list[i].name
                        except:
                            if CosmicFishPyLib.__feedback__ > 1: print 'FAIL'
                
    # -----------------------------------------------------------------------------------
    
    def add_fisher_matrix(self, fisher):
        """
        Add a set of Fisher matrices to the already existing set.
        Rejects Fisher matrices if the name is double i.e. the name is the unique 
        identifier of the Fisher matrix.
        Checks wether the elements that are passed are really Fisher matrices.

        :param fisher_list: Fisher matrix or list of Fisher matrices.
        :type fisher_list: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix` or :class:`list` of :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`

        """
        fisher_in = copy.deepcopy(fu.make_list(fisher))
        #
        for i in fisher_in:
            # check wether all the input elements are Fisher matrices:
            if not isinstance( i, fm.fisher_matrix ):
                raise ValueError('The input element is not a Fisher matrix of type fisher_matrix.')
            else:
                # get the name of the Fisher matrix:
                name_temp = copy.deepcopy(i.name)
                # check wether it is already an element of the list:
                if name_temp in self.fisher_name_list:
                    raise ValueError('A Fisher matrix with name '+name_temp+' already exists.')
                else:
                    # add it to the Fisher list:
                    self.fisher_name_list.append(name_temp)
                    self.fisher_list.append( copy.deepcopy(i) )
        
    # -----------------------------------------------------------------------------------
    
    def get_fisher_matrix(self, names=None):
        """
        Returns the list of Fisher matrices corresponding to the given names.

        :param names: names of the Fisher matrices.
        :type names: a :class:`string` or a :class:`list` of :class:`string`
        :returns: a list containing the desired Fisher matrices. Notice if the name is not found
            in the list no error is risen and and the entrance is just ignored.
        :rtype: a :class:`list` of :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix` 

        """
        # process names:
        if names==None:
            names_temp = self.fisher_name_list
        else:
            names_temp = [ i for i in fu.make_list(names) if i in self.fisher_name_list ]
        # get the fishers:
        fisher_list = []
        for name in names_temp:
            for fish in self.fisher_list:
                if fish.name == name:
                    fisher_list.append(fish)
            
        return fisher_list
    
    # -----------------------------------------------------------------------------------
    
    def delete_fisher_matrix(self, names=None):
        """
        Delete the fisher matrix or the fisher matrices in names from the Fisher list.

        :param names: names of the Fisher matrices to delete.
        :type names: a :class:`string` or a :class:`list` of :class:`string`

        """
        # process names:
        if names==None:
            names_temp = self.fisher_name_list
        else:
            names_temp = [ i for i in fu.make_list(names) if i in self.fisher_name_list ]
        # get names to delete:
        names_to_delete = [ i for i in names_temp if i in self.fisher_name_list ]
        # get indeces to delete:
        indexes_to_delete = [ i for i, s in enumerate(self.fisher_name_list) if s in names_to_delete ]
        # delete them. We have to use this way because fisher_matrix might not have a delete method and we do not want to move in memory the arrays.
        for i in xrange(len(indexes_to_delete)):
            self.fisher_name_list.pop(indexes_to_delete[i]-i)
            self.fisher_list.pop(indexes_to_delete[i]-i)
        
    # -----------------------------------------------------------------------------------
    
    def get_parameter_list(self, names=None):
        """
        Returns the list of parameter names of all the matrices identified in names.
        
        :param names: names of the Fisher matrices.
        :type names: a :class:`string` or a :class:`list` of :class:`string`
        :returns: a list containing the parameter names. 
        :rtype: a :class:`list` of :class:`string`

        """
        # process names:
        if names==None:
            names_temp = self.fisher_name_list
        else:
            names_temp = [ i for i in fu.make_list(names) if i in self.fisher_name_list ]
        # get the parameter names:
        parameter_names = []
        for i in self.get_fisher_matrix(names_temp):
            for j in i.param_names:
                if j not in parameter_names: parameter_names.append(j)
        
        return parameter_names
    
    # -----------------------------------------------------------------------------------
    
    def reshuffle(self, params, names=None):
        """
        Reshuffles all the Fisher matrices. 
        
        :param params: parameters to reshuffle.
        :type params: a :class:`string` or a :class:`list` of :class:`string`
        :param names: names of the Fisher matrices.
        :type names: a :class:`string` or a :class:`list` of :class:`string`
        :returns: a new Fisher list with the reshuffled Fishers
        :rtype: a :class:`cosmicfish_pylib.fisher_plot_analysis.CosmicFish_FisherAnalysis`

        """
        # process names:
        if names==None:
            names_temp = self.fisher_name_list
        else:
            names_temp = [ i for i in fu.make_list(names) if i in self.fisher_name_list ]
        # create the empty Fisher list:
        temp = CosmicFish_FisherAnalysis()
        # get the parameter names:
        for fish in self.get_fisher_matrix(names_temp):
            parameters = [ par for par in params if par in fish.get_param_names() ]
            if len(parameters) == 0: continue # skip the element if it has no parameters
            temp.add_fisher_matrix( fisher = fo.reshuffle( fish, parameters ) )
        
        return temp
    
    # -----------------------------------------------------------------------------------
    
    def marginalise(self, params, names=None):
        """
        Marginalise all the Fisher matrices over all the parameters that are not in params. 
        
        :param params: list of names of the parameters of the output Fisher matrix,
            in the order that will appear in the output Fisher matrix. All other parameters 
            will be marginalized over.
        :type params: a :class:`string` or a :class:`list` of :class:`string`
        :param names: names of the Fisher matrices.
        :type names: a :class:`string` or a :class:`list` of :class:`string`
        :returns: a new Fisher list with the reshuffled Fishers
        :rtype: a :class:`cosmicfish_pylib.fisher_plot_analysis.CosmicFish_FisherAnalysis`

        """
        # process names:
        if names==None:
            names_temp = copy.deepcopy( self.fisher_name_list )
        else:
            names_temp = copy.deepcopy( [ i for i in fu.make_list(names) if i in self.fisher_name_list ] )
        
        # create the empty Fisher list:
        temp = CosmicFish_FisherAnalysis()
        # get the parameter names:
        for fish in self.get_fisher_matrix(names_temp):
            parameters = [ par for par in params if par in fish.get_param_names() ]
            if len(parameters)>0:
                temp.add_fisher_matrix( fisher = fo.marginalise( fish, parameters ) )

        return temp
    
    # -----------------------------------------------------------------------------------
    
    def get_parameter_latex_names(self, names=None):
        """
        Returns a dictionary mapping parameter names and latex parameter names.
        
        :param names: names of the Fisher matrices.
        :type names: a :class:`string` or a :class:`list` of :class:`string`
        :returns: a dictionary containing parameter names and the LaTeX parameter names. 
        :rtype: :class:`dict`

        """
        # process names:
        if names==None:
            names_temp = self.fisher_name_list
        else:
            names_temp = [ i for i in fu.make_list(names) if i in self.fisher_name_list ]
        # get the parameter names:
        
        # get the names in latex:
        latex_names = {}
        for name in self.get_parameter_list(names_temp):
            for fish in self.get_fisher_matrix(names_temp):
                try:
                    latex_names[name] = fish.get_param_name_latex(name)
                except:
                    continue
                if latex_names.has_key(name): continue
        
        return latex_names
    
    # -----------------------------------------------------------------------------------
    
    def compute_plot_range(self, params=None, confidence_level=0.68, names=None, nice=True ):
        """
        Function that computes a meaningfull plot range for the plots involving the specified parameters
        and the specified Fisher names.

        :param params: name of the parameter or list of names of parameters.
        :type params: a :class:`string` or a :class:`list` of :class:`string`
        :param confidence_level: (optional) Confidence Level of the bounds. Default 68%.
        :type confidence_level: :class:`float`
        :param names: names of the Fisher matrices.
        :type names: a :class:`string` or a :class:`list` of :class:`string`
        :param nice: wether the number is properly rounded to be nice.
        :type nice: :class:`bool`
        :returns: a dictionary of name and bounds
        :rtype: :class:`dict`

        """
        # process names:
        if names==None:
            names_temp = self.fisher_name_list
        else:
            names_temp = [ i for i in fu.make_list(names) if i in self.fisher_name_list ]
        # process parameters:
        total_paramnames_list = self.get_parameter_list(names_temp)
        if params==None:
            params_temp = total_paramnames_list
        else:
            params_temp = [ i for i in fu.make_list(params) if i in total_paramnames_list ]
        # get the fishers:
        fisher_temp_list = self.get_fisher_matrix(names_temp)
        # compute the confidence coefficient:
        confidence_coefficient = fu.confidence_coefficient( confidence_level )
        # get the ranges:
        range = []
        for i in params_temp:
            lower_bound = []
            upper_bound  = []
            for j in fisher_temp_list:
                # get the index of the parameter:
                try:
                    index = j.get_param_index(i)
                except: 
                    continue
                # get sigma and fiducial:
                sigma    = j.get_fisher_inverse()[index,index]
                sigma    = math.sqrt( sigma )
                fiducial = j.get_fiducial(i)
                # apply a coefficient to safeguard:
                if not nice:
                    sigma = confidence_coefficient*sigma
                # store the values:
                if nice:
                    lower_bound.append(fu.significant_digits( (fiducial-confidence_coefficient*sigma, sigma), mode=2 ) )
                    upper_bound.append(fu.significant_digits( (fiducial+confidence_coefficient*sigma, sigma), mode=0 ) )
                else:
                    lower_bound.append(fiducial-sigma)
                    upper_bound.append(fiducial+sigma)
            # decide what to use:
            upper_bound = np.array(upper_bound)
            lower_bound = np.array(lower_bound)
            range.append([ float(str(np.amin(lower_bound))), float(str(np.amax(upper_bound))) ])
        
        dict = {}
        for i,j in zip(params_temp,xrange(len(params_temp))):
            dict[i] = range[j]
            
        return dict
            
    # -----------------------------------------------------------------------------------
    
    def compute_gaussian( self, params=None, confidence_level=0.68, names=None, num_points=100, normalized=False, nice_bounds=True ):
        """
        Function that computes the (1D) gaussian distribution of a given parameter.
        Returns a dictionary with all the meaningul information about the gaussian.

        :param params: name of the parameter or list of names of parameters.
        :type params: a :class:`string` or a :class:`list` of :class:`string`
        :param confidence_level: (optional) Confidence Level of the bounds. Default 68%.
        :type confidence_level: :class:`float`
        :param names: names of the Fisher matrices.
        :type names: a :class:`string` or a :class:`list` of :class:`string`
        :param num_points: number of (x,y) points. 
        :type num_points: :class:`int`
        :param normalized: wether the distribution is normalized or not.
        :type normalized: :class:`bool`
        :param nice_bounds: wether the x range is properly rounded to be nice or not.
        :type nice_bounds: :class:`bool`
        :returns: a dictionary mapping name and parameter to a tuple of: [x, y, [fiducial,sigma]]
        :rtype: :class:`dict`
        
        """
        # process names:
        if names==None:
            names_temp = self.fisher_name_list
        else:
            names_temp = [ i for i in fu.make_list(names) if i in self.fisher_name_list ]
        # process parameters:
        total_paramnames_list = self.get_parameter_list(names_temp)
        if params==None:
            params_temp = total_paramnames_list
        else:
            params_temp = [ i for i in fu.make_list(params) if i in total_paramnames_list ]
        # get the fishers:
        fisher_temp_list = self.get_fisher_matrix(names_temp)
        # compute the confidence coefficient:
        confidence_coefficient = fu.confidence_coefficient( confidence_level )
        # get the plot ranges:
        plot_ranges = self.compute_plot_range( params=params_temp, confidence_level=confidence_level, names=names_temp, nice=nice_bounds )
        # get the distributions:
        gaussian_distro = {}
        for par1 in params_temp:
            dict_names = {}
            for mat,name in zip(fisher_temp_list,names_temp):
                # get the index of the parameter:
                try:
                    index = mat.get_param_index(par1)
                except:
                    dict_names[name] = [np.array([0.0]), np.array([0.0]), [0.0,0.0]]
                    continue
                # get sigma and fiducial:
                sigma    = math.sqrt( mat.get_fisher_inverse()[index,index] )
                fiducial = mat.get_fiducial(par1)
                # get optimal x:
                if nice_bounds:
                    lower_bound = plot_ranges[par1][0]
                    upper_bound = plot_ranges[par1][1]
                else: 
                    lower_bound = fiducial -confidence_coefficient*sigma
                    upper_bound = fiducial +confidence_coefficient*sigma
                x_points = np.linspace( lower_bound, upper_bound, num_points )
                # get y:
                if normalized:
                    y_points = np.array([ np.exp( -( x-fiducial )**2/(2.0*sigma**2) )/( sigma*np.sqrt(2.0*math.pi) ) for x in x_points ])
                else:
                    y_points = np.array([ np.exp( -( x-fiducial )**2/(2.0*sigma**2) ) for x in x_points ])
                    
                dict_names[name] = [x_points, y_points, [fiducial,sigma]]
            gaussian_distro[par1] = dict_names
        
        return gaussian_distro
    
    # -----------------------------------------------------------------------------------
    
    def compute_ellipse(self, params1=None, params2=None, confidence_level=0.68, names=None, num_points=100):
        """
        Function that computes the (2D) ellipses for a given parameters combination.
        Returns a dictionary with all the meaningul information about the ellipses.

        :param params1: name of the first parameter or list of names of parameters.
        :type params1: a :class:`string` or a :class:`list` of :class:`string`
        :param params2: name of the second parameter or list of names of parameters.
        :type params3: a :class:`string` or a :class:`list` of :class:`string`
        :param confidence_level: (optional) Confidence Level of the bounds. Default 68%.
        :type confidence_level: :class:`float`
        :param names: names of the Fisher matrices.
        :type names: a :class:`string` or a :class:`list` of :class:`string`
        :param num_points: number of (x,y) points. 
        :type num_points: :class:`int`
        :returns: a dictionary mapping name and parameters to a tuple of: [x, y, [fiducial_x, fiducial_y, coeff_a, coeff_b, theta_0]]
        :rtype: :class:`dict`    
        
        """
        # process names:
        if names==None:
            names_temp = self.fisher_name_list
        else:
            names_temp = [ i for i in fu.make_list(names) if i in self.fisher_name_list ]
        # process parameters:
        total_paramnames_list = self.get_parameter_list(names_temp)
        if params1==None:
            params_temp_1 = total_paramnames_list
        else:
            params_temp_1 = [ i for i in fu.make_list(params1) if i in total_paramnames_list ]
        if params2==None:
            params_temp_2 = total_paramnames_list
        else:
            params_temp_2 = [ i for i in fu.make_list(params2) if i in total_paramnames_list ]    
        # get the fishers:
        fisher_temp_list = self.get_fisher_matrix(names_temp)
        # compute the confidence coefficient:
        confidence_coefficient = fu.confidence_coefficient( confidence_level )
        # get the distributions:
        gaussian_distro = {}
        for par1 in params_temp_1:
            dict_names_1 = {}
            for par2 in params_temp_2:
                dict_names_2 = {}
                for mat,name in zip(fisher_temp_list,names_temp):
                    # get the index of the parameter:
                    try:
                        index_1 = mat.get_param_index(par1)
                        index_2 = mat.get_param_index(par2)
                    except:
                        dict_names_2[name] = [np.array([0.0]), np.array([0.0]), [0.0, 0.0, 0.0, 0.0, 0.0]]
                        continue
                    # get sigmas:
                    sigma_x  = mat.get_fisher_inverse()[index_1,index_1]
                    sigma_y  = mat.get_fisher_inverse()[index_2,index_2]
                    sigma_xy = mat.get_fisher_inverse()[index_1,index_2]
                    # get fiducial:
                    fiducial_x = mat.get_fiducial(par1)
                    fiducial_y = mat.get_fiducial(par2)
                    # compute the ellipse coefficients:
                    coeff_a = confidence_coefficient*math.sqrt((sigma_x + sigma_y)/2.0 + math.sqrt( (sigma_x - sigma_y)**2/4.0 + sigma_xy**2 ))
                    coeff_b = confidence_coefficient*math.sqrt((sigma_x + sigma_y)/2.0 - math.sqrt( (sigma_x - sigma_y)**2/4.0 + sigma_xy**2 ))
                    theta_0 = math.atan2( (2.0*sigma_xy), (sigma_x - sigma_y) )/2.0
                    # generate the ellipses
                    angles  = np.linspace( 0, 2.0*math.pi, num_points )
                    x_points = np.array( [+coeff_a*math.cos(theta)*math.cos(theta_0)
                                          -coeff_b*math.sin(theta)*math.sin(theta_0) + fiducial_x for theta in angles ] )
                    y_points = np.array( [+coeff_a*math.cos(theta)*math.sin(theta_0)
                                          +coeff_b*math.sin(theta)*math.cos(theta_0) + fiducial_y for theta in angles ] )
                    # save the result:
                    dict_names_2[name] = [x_points, y_points, [fiducial_x, fiducial_y, coeff_a, coeff_b, theta_0]]
                dict_names_1[par2] = dict_names_2
            gaussian_distro[par1] = dict_names_1

        return gaussian_distro

    # -----------------------------------------------------------------------------------

# ***************************************************************************************