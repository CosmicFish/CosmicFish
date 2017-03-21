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
.. module:: fisher_derived
   :platform: Unix
   :synopsis: Module that contains the fisher_derived class and the operations defined on it.
              This is meant to handle efficiently and safely Jacobian matrices transforming a 
              Fisher matrix for some parameters to a Fisher matrix for some other parameters.

.. moduleauthor:: Marco Raveri <mraveri@uchicago.edu> for the CosmicFish code.

"""

# ***************************************************************************************

import os
import math
import copy
import numpy as np
import utilities as fu
import fisher_matrix as fm

# ***************************************************************************************

class fisher_derived():
    """
    This class contains the relevant code to define a matrix that
    contains the relevant information to reparametrize a Fisher matrix.
    Generally this is a rectangular matrix containing the Jacobian of the transformation
    from the original Fisher to the derived one.
    
    :ivar derived_matrix: Numpy array containing the Jacobian of the transformation between the Fisher matrix and the derived Fisher matrix. Passed to the constructor of by file.
    :ivar path: Absolute path of the input Jacobian matrix. Computed at initialization if passing a file.
    :ivar name: Name of the input Jacobian matrix. Computed at initialization if passing a file.
    :ivar indir: Absolute path of the directory containing the input Jacobian matrix. Computed at initialization if passing a file.
    :ivar num_params: Number of base parameters of the Jacobian matrix.
    :ivar num_derived: Number of derived parameters.
    :ivar param_names: Names of the base parameters. Used as the identifier of the parameters. Initialized, if possible, through a .paramnames file.
    :ivar param_names_latex: LaTeX names of the base parameters.
    :ivar param_fiducial: Numpy array with the values of the fiducial of the base parameters. Passed to the constructor or by file.
    :ivar derived_param_names: Names of the derived parameters. Used as the identifier of the parameters. Initialized, if possible, through a .paramnames file.
    :ivar derived_param_names_latex: LaTeX names of the derived parameters.
    :ivar derived_param_fiducial: Numpy array with the values of the fiducial of the derived parameters. Passed to the constructor or by file.
              
    """

    # -----------------------------------------------------------------------------------
    
    # class getters:
    
    def get_derived_matrix(self):
        """ :returns: the derived Jacobian matrix. """
        return self.derived_matrix
    
    def get_param_names(self):
        """ :returns: the base parameter names. """
        return self.param_names
    
    def get_param_names_latex(self):
        """ :returns: the LaTeX version of the base parameter names. """
        return self.param_names_latex
    
    def get_param_fiducial(self):
        """ :returns: the base parameter fiducial values. """
        return self.param_fiducial
    
    def get_derived_param_names(self):
        """ :returns: the derived parameters names. """
        return self.derived_param_names
    
    def get_derived_param_names_latex(self):
        """ :returns: the LaTeX version of the derived parameters names. """
        return self.derived_param_names_latex
    
    def get_derived_param_fiducial(self):
        """ :returns: the derived parameter fiducial values. """
        return self.derived_param_fiducial
    
    # -----------------------------------------------------------------------------------

    def __init__( self, derived_matrix=None, 
                  param_names=None, derived_param_names=None, 
                  param_names_latex=None, derived_param_names_latex=None, 
                  fiducial=None, fiducial_derived=None, file_name=None ):
        """
        **fisher_derived class constructor**. The class constructor will read from file the 
        Fisher derived matrix if it is initialized with the name of a file (and the file exists). 
        Otherwise it will read the matrix and the parameter names as passed by the user.

        :param derived_matrix: array containing the input Jacobian matrix.
        :type derived_matrix: 2D :class:`list` or :class:`numpy.array`
        :param  param_names: names of the parameters of the derived parameters. If initialized from file it will
            read them if a file file_name.paramnames is found. If it is none when itialized from
            python it will just use some defaults names (p1, p2,...).
        :type param_names: :class:`list` of :class:`string`
        :param derived_param_names: names of the derived parameters. If initialized from file it will
            read them if a file file_name.paramnames and expects them to have a * appened. 
            If it is none when itialized from python it will just use some defaults names.
        :type derived_param_names: :class:`list` of :class:`string`
        :param param_names_latex: LaTeX names of the parameters of the Jacobian matrix also appearing
            in the Fisher matrix.
        :type param_names_latex: :class:`list` of :class:`string`
        :param derived_param_names_latex: LaTeX names of the parameters of the Jacobian matrix that are 
            derived parameters.
        :type derived_param_names_latex: :class:`list` of :class:`string`
        :param fiducial: values of the fiducial parameters of the Fisher matrix. If initialized from file it will have the value found in .paramnames.
            If not found on a file it will be zero. Passing it to the constructor overwrites the values.
        :type fiducial: :class:`list` of :class:`float` or :class:`numpy.array`
        :param fiducial_derived: values of the fiducial derived parameters. If initialized from file it will have the value found in .paramnames.
            If not found on a file it will be zero. Passing it to the constructor overwrites the values.
        :type fiducial_derived: :class:`list` of :class:`float` or :class:`numpy.array`
        :param file_name: name of the file (and path) of the input Jacobian matrix.
        :type file_name: :class:`string`
        
        """
        # check that the input is legal:
        if derived_matrix is None and file_name is None:
            raise ValueError('Error in initializing the Fisher Jacobian matrix: derived_matrix and file_name are both None.')
        # initialize class members:
        self.derived_matrix            = np.array([]) 
        self.path                      = ''           
        self.name                      = ''           
        self.indir                     = ''           
        self.num_params                = 0            
        self.num_derived               = 0            
        self.param_names               = []           
        self.param_names_latex         = []           
        self.param_fiducial            = np.array([]) 
        self.derived_param_names       = []           
        self.derived_param_names_latex = []           
        self.derived_param_fiducial    = np.array([]) 
        # initialize the Jacobian matrix:
        if derived_matrix is None:
            # initialize from file:
            self.derived_matrix = np.loadtxt(file_name)
            # get the file name and path:
            self.path = os.path.abspath(file_name)
            self.name = os.path.splitext(os.path.basename(file_name))[0]
            self.indir = os.path.dirname(self.path)
        else:
            # read the fisher matrix from input:
            self.derived_matrix = np.array(derived_matrix)
        # protection against 1D Fisher matrices
        if np.ndim(self.derived_matrix) == 0:
            self.derived_matrix = np.array([[self.derived_matrix]])
        elif np.ndim(self.derived_matrix) == 1:
            self.derived_matrix = np.array([self.derived_matrix])
        # get the number of parameters:
        self.num_params  = self.derived_matrix.shape[0]
        self.num_derived = self.derived_matrix.shape[1]
        # load the parameter names:
        if param_names is None or derived_param_names is None:
            try:
                self.load_paramnames_from_file()
            except ValueError:
                self.param_names               = [ 'p'+str(i+1) for i in xrange(self.num_params) ]
                self.param_names_latex         = [ 'p'+str(i+1) for i in xrange(self.num_params) ]
                self.param_fiducial            = np.array( [0.0 for i in self.param_names] ) 
                self.derived_param_names       = [ 'p'+str(i+1) for i in xrange(self.num_params, self.num_derived+self.num_params) ]
                self.derived_param_names_latex = [ 'p'+str(i+1) for i in xrange(self.num_params, self.num_derived+self.num_params) ]
                self.derived_param_fiducial    = np.array( [0.0 for i in self.derived_param_names] ) 
        else:
            self.param_names                 = copy.deepcopy(param_names)
            self.param_names_latex           = copy.deepcopy(param_names)
            self.param_fiducial              = np.array( [0.0 for i in self.param_names] ) 
            self.derived_param_names         = copy.deepcopy(derived_param_names)
            self.derived_param_names_latex   = copy.deepcopy(derived_param_names)
            self.derived_param_fiducial      = np.array( [0.0 for i in self.derived_param_names] ) 
        # over write what is explicitly given:
        if param_names is not None:
            self.param_names       = copy.deepcopy( param_names )
        if param_names_latex is not None:
            self.param_names_latex = copy.deepcopy(param_names_latex)
        if fiducial is not None:
            self.param_fiducial    = np.array( copy.deepcopy(fiducial) )
        if derived_param_names is not None:
            self.derived_param_names = copy.deepcopy( derived_param_names )
        if derived_param_names_latex is not None:
            self.derived_param_names_latex = copy.deepcopy( derived_param_names_latex )
        if fiducial_derived is not None:
            self.derived_param_fiducial = np.array( copy.deepcopy(fiducial_derived) )
        # check the validity:
        if len(self.param_names) != self.num_params:
            raise ValueError('The input param_names has not '+str(self.num_params)+' elements.')
        if len(self.param_names_latex) != self.num_params:
            raise ValueError('The input param_names_latex has not '+str(self.num_params)+' elements.')
        if len(self.param_fiducial) != self.num_params:
            raise ValueError('The input param_fiducial has not '+str(self.num_params)+' elements.')
        if len(self.derived_param_names) != self.num_derived:
            raise ValueError('The input derived_param_names has not '+str(self.num_derived)+' elements.')
        if len(self.derived_param_names_latex) != self.num_derived:
            raise ValueError('The input derived_param_names_latex has not '+str(self.num_derived)+' elements.')
        if len(self.derived_param_fiducial) != self.num_derived:
            raise ValueError('The input derived_param_fiducial has not '+str(self.num_derived)+' elements.')
        
    # -----------------------------------------------------------------------------------       
        
    def load_paramnames_from_file( self, file_name=None ):
        """
        Loads the paramnames array, of a derived Fisher matrix, from a file

        :param file_name: (optional) file name and path of the parameter names file.
            If file_name is None this reads the file self.name+.paramnames.
        :type file_name: :class:`string`
        
        """
        if file_name is None:
            name = self.indir+'/'+self.name+'.paramnames'
        else:
            name = file_name
        # check wether the file exists:
        if not os.path.isfile(name):
            raise ValueError('No .paramnames file found at: '+name)
        # init:
        param_names               = []
        param_names_latex         = []
        param_fiducial            = []
        derived_param_names       = []
        derived_param_names_latex = []
        derived_param_fiducial    = []
        # parse:
        with open(name) as f:
            for line in f:
                if line[0]!='#' and line[1]!='#':
                    split_line = [ i.strip() for i in line.split('    ') ]
                    split_line = [ i for i in split_line if i is not '' ]
                    temp_line=split_line[0].strip()
                    if temp_line[-1]=='*':
                        derived_param_names.append(temp_line[0:-1])
                        if len(split_line) == 1:
                            # latex and fiducial missing:
                            derived_param_names_latex.append(temp_line[0:-1])
                            derived_param_fiducial.append(0.0)
                        elif len(split_line) == 2:
                            # one of the two is missing:
                            try:
                                temp = float(split_line[1].strip())
                                derived_param_fiducial.append(temp)
                                derived_param_names_latex.append(temp_line[0:-1])
                            except:
                                temp = 0.0
                                derived_param_fiducial.append(temp)
                                derived_param_names_latex.append(split_line[1].strip())
                        elif len(split_line) >= 3:
                            # nothing is missing:
                            derived_param_names_latex.append(split_line[1].strip())
                            derived_param_fiducial.append(float(split_line[2].strip()))
                    else:
                        param_names.append(temp_line)
                        if len(split_line) == 1:
                            # latex and fiducial missing:
                            param_names_latex.append(split_line[0].strip())
                            param_fiducial.append(0.0)
                        elif len(split_line) == 2:
                            # one of the two is missing:
                            try:
                                temp = float(split_line[1].strip())
                                param_fiducial.append(temp)
                                param_names_latex.append(split_line[0].strip())
                            except:
                                temp = 0.0
                                param_fiducial.append(temp)
                                param_names_latex.append(split_line[1].strip())
                        elif len(split_line) >= 3:
                            # nothing is missing:
                            param_names_latex.append(split_line[1].strip())
                            param_fiducial.append(float(split_line[2].strip()))
                        
        # check the valitidy of the result:
        if len(derived_param_names) != self.num_derived:
            raise ValueError('Wrong number of derived parameters in the .paramnames file')
        
        # save the result:
        self.param_names         = param_names
        self.param_names_latex   = param_names_latex
        self.param_fiducial      = param_fiducial
        
        self.derived_param_names       = derived_param_names
        self.derived_param_names_latex = derived_param_names_latex
        self.derived_param_fiducial    = derived_param_fiducial
        
    # -----------------------------------------------------------------------------------       
        
    def add_derived( self, fisher_matrix, preserve_input=False ):
        """
        This function computes the derived fisher_matrix given an input Fisher matrix based on the Jacobian
        contained in derived_matrix.

        :param fisher_matrix: input Fisher matrix that will be used as a base for the derived Fisher matrix.
        :type fisher_matrix: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`
        :param preserve_input: wether to preserve input parameters in the output Fisher. 
                Default to false because it might lead to strage results if not used properly.
        :type preserve_input: :class:`bool`        
        :return: output Fisher matrix with derived parameters. 
        :rtype: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`

        """
        # check the type of the input Fisher matrix:
        if ( not isinstance(fisher_matrix, fm.fisher_matrix) ):
            raise ValueError('Error, input fisher_matrix is not a fisher_matrix')
        # check if the fisher matrix is compatible with the fisher_derived matrix:
        if ( fisher_matrix.param_names != self.param_names):
            raise ValueError('Error, paramnames of derived matrix and fisher matrix do not match')
        # check wether the two fiducials are the same:
        if not np.allclose(self.param_fiducial, fisher_matrix.param_fiducial ):
            raise ValueError('Error, fiducial of derived matrix and fisher matrix do not match')
        # create the new parameter names:
        if preserve_input:
            temp_param_names = fisher_matrix.param_names + self.derived_param_names
            temp_param_names_latex = fisher_matrix.param_names_latex + self.derived_param_names_latex
            temp_param_fiducial = np.append(fisher_matrix.param_fiducial, self.derived_param_fiducial)
        else:
            temp_param_names = self.derived_param_names
            temp_param_names_latex = self.derived_param_names_latex
            temp_param_fiducial = self.derived_param_fiducial
        # prepare the derived matrix to preserve the original prameters:
        if preserve_input:
            temp_derived_matrix = np.hstack( (np.identity( self.num_params ), self.derived_matrix) )
        else:
            temp_derived_matrix = self.derived_matrix
        # compute the derived inverse Fisher matrix:
        A = temp_derived_matrix
        AT = np.transpose( temp_derived_matrix )
        original_fisher = fisher_matrix.get_fisher_inverse()
        fisher_new_inverse = np.dot( np.dot(AT, fisher_matrix.get_fisher_inverse() ), A )
        # get initial PCA and spectral cutoff:
        (w,v)              = fisher_matrix.PCA()
        initial_cutoff     = fisher_matrix.fisher_cutoff
        constrained_eigen  = [ i for i in w if i> 1.1*initial_cutoff ]
        initial_best_mode  = 1.0/np.amax(constrained_eigen)
        initial_worse_mode = 1.0/np.amin(constrained_eigen)
        # get the new matrix spectrum:
        (w_new, v_new) = np.linalg.eigh( fisher_new_inverse )
        # establish the cutoff. This is chosen to leave unaltered the original eigenvalues of the Fisher matrix as much as we can.
        spectral_width = ( np.log10( fisher_matrix.fisher_spectrum ) -(np.log10( initial_worse_mode) - np.log10( initial_best_mode )) ) /2.0
        upper_cutoff = 10.0**(np.log10(initial_worse_mode)+spectral_width)
        lower_cutoff = 10.0**(np.log10(initial_best_mode)-spectral_width)
        # check the cutoff:
        if np.amin(np.abs(w_new))<lower_cutoff or np.amax(np.abs(w_new))>upper_cutoff:
            print 'WARNING: in add_derived name:', self.name, ' fisher:', fisher_matrix.name
            print '** derived parameters are strongly degenerate and might alter the quality of the original Fisher matrix.'
            print '** Try removing degenerate parameters from the Fisher matrix and the derived matrix to get rid of this warning.'
        # apply the cutoff:
        if preserve_input:
            temp = np.zeros((self.num_params + self.num_derived, self.num_params + self.num_derived), float)
        else:
            temp = np.zeros(( self.num_derived, self.num_derived), float)
        for i, el in enumerate(temp):
            if w_new[i] < lower_cutoff:
                temp[i,i] = 1.0/( lower_cutoff )
            elif np.abs(w_new[i]) > upper_cutoff:
                temp[i,i] = 1.0/( upper_cutoff )
            else:
                temp[i,i] = 1.0/w_new[i]
        # rebuild the Fisher matrix:
        fisher_new_inverse =  np.dot( v_new, np.dot(temp, np.transpose(v_new) ))
        # return the new Fisher matrix:
        fisher_matrix_new = fm.fisher_matrix( fisher_matrix=fisher_new_inverse, param_names=temp_param_names, param_names_latex=temp_param_names_latex, fiducial=temp_param_fiducial )
        fisher_matrix_new.path  = fisher_matrix.path
        if preserve_input:
            fisher_matrix_new.name  = fisher_matrix.name
        else:
            fisher_matrix_new.name  = fisher_matrix.name+'_derived'
        fisher_matrix_new.indir = fisher_matrix.indir
        # return the new Fisher matrix with the derived parameters
        return fisher_matrix_new
    
    # -----------------------------------------------------------------------------------
        
# ***************************************************************************************    
