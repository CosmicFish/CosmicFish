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
.. module:: fisher_matrix
   :platform: Unix
   :synopsis: Module that contains the fisher_matrix class and the operations defined on it.
              The idea is that of creating a Fisher matrix safely, with all operations being
              well defined and safe-guarded once the input matrix is accepted as a Fisher
              matrix.

.. moduleauthor:: Marco Raveri <mraveri@uchicago.edu> for the CosmicFish code.

"""

# ***************************************************************************************

import os
import math
import copy
import numpy as np
from . import utilities as fu

# ***************************************************************************************

class fisher_matrix():
    """
    This class contains the relevant code to define a fisher matrix
    and basic operations on it.

    :ivar fisher_cutoff: cutoff for the spectrum of the Fisher matrix. Parameter of the class. Starts at 10**(-9) but, if needed, is fixed during computations.
    :ivar fisher_spectrum: maximum condition number allowed for the Fisher matrix. Worse constrained modes that go above this value will be flattened to this value.
    :ivar fisher_matrix: numpy array with the fisher matrix. Passed to the constructor of by file.
    :ivar path: absolute path of the input Fisher matrix. Computed at initialization if passing a file or just an empty string.
    :ivar name: name of the input Fisher matrix. Computed at initialization if passing a file or just an empty string.
    :ivar indir: absoulte path of the directory of the input Fisher matrix. Computed at initialization if passing a file or just an empty string.
    :ivar num_params: number of parameters of the input Fisher matrix. Computed at initialization.
    :ivar fisher_eigenvalues: eigenvalues of the input Fisher matrix. Computed at initialization or by the PCA function.
    :ivar fisher_eigenvectors: eigenvectors of the input Fisher matrix. Computed at initialization or by the PCA function.
    :ivar fisher_matrix_inv: inverse of the input Fisher matrix. Computed at initialization or by inverse_fisher_matrix.
    :ivar param_names: name of the parameters of the input Fisher matrix. Passed to the constructor of by file.
    :ivar param_names_latex: LaTeX name of the parameters of the Fisher matrix. Passed to the constructor of by file.
    :ivar param_fiducial: numpy array with the values of the fiducial parameters. Passed to the constructor of by file.
    :ivar param_names_dict: a dictionary that maps parameter names to numbers and vice versa.

    .. automethod:: cosmicfish_pylib.fisher_matrix.fisher_matrix.__add__
    .. automethod:: cosmicfish_pylib.fisher_matrix.fisher_matrix.__eq__
    .. automethod:: cosmicfish_pylib.fisher_matrix.fisher_matrix.__ne__

    """

    # -----------------------------------------------------------------------------------

    # class getters:

    def get_fisher_matrix(self):
        """ :returns: the fisher matrix as a numpy array. """
        return self.fisher_matrix

    def get_fisher_eigenvalues(self):
        """ :returns: the eigenvalues of the Fisher matrix as a numpy array. """
        return self.fisher_eigenvalues

    def get_fisher_eigenvectors(self):
        """ :returns: the eigenvectors of the Fisher matrix. """
        return self.fisher_eigenvectors

    def get_fisher_inverse(self):
        """ :returns: the inverse of the Fisher matrix. """
        return self.fisher_matrix_inv

    def get_param_names(self):
        """ :returns: the parameter names of the Fisher matrix. """
        return self.param_names

    def get_param_names_latex(self):
        """ :returns: the parameter names, in LaTeX format of the Fisher matrix. """
        return self.param_names_latex

    def get_param_fiducial(self):
        """ :returns: the fiducial values of the parameters of the Fisher matrix. """
        return self.param_fiducial

    # -----------------------------------------------------------------------------------

    # advanced class getters:

    def get_param_name(self, number):
        """
        Returns the name of the parameter corresponding to the given number.

        :param number: number of the parameter or list of numbers. Notice that parameters are numbered starting from 1.
        :type number: :class:`int` or a :class:`list` of :class:`int`
        :returns: the name or a list of names of the parameters.
        :rtype: :class:`string` or a :class:`list` of :class:`string`.

        """
        if isinstance(number, int):
            return self.param_names_dict[number]
        elif len(number)>1:
            return [ self.param_names_dict[i] for i in number ]

    def get_param_index(self,name):
        """
        Returns the index of a parameter as specified by his name. Notice that indices starts at 0.

        :param name: input name or list of names of the parameters.
        :type name: :class:`string` or a :class:`list` of :class:`string`
        :returns: the index of the parameter or a list of numbers.
        :rtype: :class:`int` or a :class:`list` of :class:`int`

        """
        if isinstance(name, str):
            return self.param_names_dict[name]-1
        elif len(name)>1:
            return [ self.param_names_dict[i]-1 for i in name ]

    def get_param_number(self,name):
        """
        Returns the number of a parameter as specified by his name.
        Notice this differs from get_param_index becasue number = index+1.

        :param name: input name or list of names of the parameters.
        :type name: :class:`string` or a :class:`list` of :class:`string`
        :returns: the index of the parameter or a list of numbers.
        :rtype: :class:`int` or a :class:`list` of :class:`int`

        """
        if isinstance(name, str):
            return self.param_names_dict[name]
        elif len(name)>1:
            return [ self.param_names_dict[i] for i in name ]

    def get_param_name_latex(self, name):
        """
        Returns the Latex name of the parameter called name.

        :param name: input name or list of names of the parameters.
        :type name: :class:`string` or a :class:`list` of :class:`string`
        :returns: the Latex name or a list of Latex names corresponding to the parameter names.
        :rtype: :class:`string` or a :class:`list` of :class:`string`.

        """
        if isinstance(name, str):
            return self.param_names_latex[self.get_param_index(name)]
        elif len(name)>1:
            return [ self.param_names_latex[i] for i in self.get_param_index(name) ]

    def get_fiducial(self, name):
        """
        Returns the fiducial of the parameter called name.

        :param name: input name or list of names of the parameters.
        :type name: :class:`string` or a :class:`list` of :class:`string`
        :returns: the fiducial or a list of fiducials corresponding to the parameter names.
        :rtype: :class:`float` or a :class:`list` of :class:`float`.

        """
        if isinstance(name, str):
            return self.param_fiducial[self.get_param_index(name)]
        elif len(name)>1:
            return [ self.param_fiducial[i] for i in self.get_param_index(name) ]

    # -----------------------------------------------------------------------------------

    def __init__( self, fisher_matrix=None, param_names=None, param_names_latex=None, fiducial=None, file_name=None ):
        """
        **fisher_matrix class constructor**. The constructor will read from file the Fisher matrix if it is initialized
        with the name of a file (and the file exists). Otherwise it will read the matrix and the
        parameter names as passed by the user.

        :param fisher_matrix: array containing the input Fisher matrix. You can initialize a Fisher matrix
            passing a matrix or by filename.
        :type fisher_matrix: 2D :class:`list` or :class:`numpy.array`
        :param param_names: names of the parameters of the Fisher matrix. If initialized from file it will
            read them if a file file_name.paramnames is found. If a file is found and param_names are given to
            the constructor the latter will be used. If param_names is None when itialized from
            python it will just use some defaults names (p1, p2,...).
        :type param_names: :class:`list` of :class:`string`
        :param param_names_latex: LaTeX names of the parameters of the Fisher matrix. If initialized from file it will
            read them if a file file_name.paramnames is found. If a file is found and param_names_latex are given to
            the constructor the latter will be used. If it is none when itialized from
            python it will just use some defaults names (p1, p2,...).
        :type param_names_latex: :class:`list` of :class:`string`
        :param fiducial: values of the fiducial parameters. If initialized from file it will have the value found in .paramnames.
            If not found on a file it will be zero. Passing it to the constructor overwrites the values.
        :type fiducial: :class:`list` of :class:`float` or :class:`numpy.array`
        :param file_name: name of the file (and path) of the input Fisher matrix.
        :type file_name: :class:`string`

        """
        # check that the input is legal:
        if fisher_matrix is None and file_name is None:
            raise ValueError('Error in initializing the Fisher matrix: fisher_matrix and file_name are both None.')
        # initialize all the objects:
        self.fisher_cutoff   = 10**(-9)
        self.fisher_spectrum = 10**(10)
        self.fisher_matrix = np.array([])
        self.path  = ''
        self.name  = ''
        self.indir = ''
        self.num_params  = 0
        self.fisher_eigenvalues  = np.array([])
        self.fisher_eigenvectors = np.array([[]])
        self.fisher_matrix_inv   = np.array([[]])
        self.param_names         = ['']
        self.param_names_latex   = ['']
        self.param_fiducial      = np.array([])
        self.param_names_dict = {}
        # initialize the Fisher matrix:
        if fisher_matrix is None:
            # initialize from file:
            self.fisher_matrix = np.loadtxt(file_name)
            # get the file name and path:
            self.path  = os.path.abspath(file_name)
            self.name  = os.path.splitext(os.path.basename(file_name))[0]
            self.indir = os.path.dirname(self.path)
        else:
            # read the fisher matrix from input:
            self.fisher_matrix = np.array(fisher_matrix)

        # protection against 1D Fisher matrices
        if np.ndim(self.fisher_matrix) == 0:
            self.fisher_matrix = np.array([[self.fisher_matrix]])
        elif np.ndim(self.fisher_matrix) == 1:
            self.fisher_matrix = np.array([self.fisher_matrix])
        # check if the Fisher matrix is symetric:
        if not np.allclose(self.fisher_matrix, np.transpose(self.fisher_matrix), rtol=1e-03, atol=1e-06 ):
            raise ValueError('The input Fisher matrix is not equal to its transpose')
        # get the number of parameters:
        self.num_params = self.fisher_matrix.shape[0]
        # load the parameter names:
        if param_names is None:
            try:
                self.load_paramnames_from_file()
            except ValueError:
                self.param_names  = [ 'p'+str(i+1) for i in range(self.num_params) ]
                self.param_names_latex = [ 'p'+str(i+1) for i in range(self.num_params) ]
                self.param_fiducial    = np.array( [0.0 for i in self.param_names] )
        else:
            self.param_names = copy.deepcopy(param_names)
            self.param_names_latex = copy.deepcopy(param_names)
            self.param_fiducial    = np.array( [0.0 for i in self.param_names] )
        # over write the fiducial if it is given:
        if param_names_latex is not None:
            self.param_names_latex = copy.deepcopy(param_names_latex)
        if fiducial is not None:
            self.param_fiducial    = copy.deepcopy(fiducial)
        self.param_fiducial = np.array( self.param_fiducial )
        # check the validity:
        if len(self.param_names) != self.num_params:
            raise ValueError('The input param_names has not '+str(self.num_params)+' elements.')
        if len(self.param_names_latex) != self.num_params:
            raise ValueError('The input param_names_latex has not '+str(self.num_params)+' elements.')
        if len(self.param_fiducial) != self.num_params:
            raise ValueError('The input param_fiducial has not '+str(self.num_params)+' elements.')
        # create a dictionary of param names:
        self.param_names_dict = {}
        for i in range(len(self.param_names)):
            self.param_names_dict[i+1] = self.param_names[i]
            self.param_names_dict[self.param_names[i]] = i+1
        # do PCA and store results:
        (self.fisher_eigenvalues,self.fisher_eigenvectors) = self.PCA()
        # protect against degenerate parameters:
        self.protect_degenerate()
        # invert the Fisher matrix and store the result:
        self.fisher_matrix_inv = self.inverse_fisher_matrix()

    # -----------------------------------------------------------------------------------

    def load_paramnames_from_file( self, file_name=None ):
        """
        Loads the paramnames array from a file.

        :param file_name: (optional) file name and path of the parameter names file.
            If file_name is None this reads the file self.name+.paramnames.

        """
        if file_name is None:
            name = self.indir+'/'+self.name+'.paramnames'
        else:
            name = file_name
        # check input:
        if not os.path.isfile(name):
            raise ValueError('No .paramnames file found at: '+name)
        # parse the param names:
        self.param_names       = []
        self.param_names_latex = []
        self.param_fiducial    = []
        with open(name) as f:
            for line in f:
                if line[0]!='#' and line[1]!='#':
                    split_line = [ i.strip() for i in line.split('    ') ]
                    split_line = [ i for i in split_line if i is not '' ]
                    self.param_names.append(split_line[0].strip())
                    if len(split_line) == 1:
                        # latex and fiducial missing:
                        self.param_names_latex.append(split_line[0].strip())
                        self.param_fiducial.append(0.0)
                    elif len(split_line) == 2:
                        # one of the two is missing:
                        try:
                            temp = float(split_line[1].strip())
                            self.param_fiducial.append(temp)
                            self.param_names_latex.append(split_line[0].strip())
                        except:
                            temp = 0.0
                            self.param_fiducial.append(temp)
                            self.param_names_latex.append(split_line[1].strip())
                    elif len(split_line) >= 3:
                        # nothing is missing:
                        self.param_names_latex.append(split_line[1].strip())
                        self.param_fiducial.append(float(split_line[2].strip()))
        self.param_fiducial = np.array( self.param_fiducial )
        # check the validity of the param names:
        if len(self.param_names) != self.num_params:
            raise ValueError('Error in load_paramnames_from_file: wrong number of parameters in the .paramnames file')

    # -----------------------------------------------------------------------------------

    def save_paramnames_to_file( self, file_name=None ):
        """
        Saves the paramnames to a file.

        :param file_name: (optional) file name and path of the parameter names file.
            If file_name is None this saves the file self.name+.paramnames.

        """
        if file_name is None:
            name = self.indir+'/'+self.name+'.paramnames'
        else:
            name = file_name
        # open the output file:
        out_file = open( name, 'w' )
        # write the header:
        out_file.write( '#\n' )
        out_file.write( '# This file contains the parameter names for a Fisher matrix.\n' )
        out_file.write( '#\n' )
        # write the parameters:
        for ind in range( self.num_params ):
            param_name = self.get_param_name( ind+1 )
            out_file.write( str(param_name)+'    '+ \
                            str(self.get_param_name_latex( param_name ))+'    '+ \
                            str(self.get_fiducial( param_name ))+'\n'
                            )
        # close the output file:
        out_file.close()

    # -----------------------------------------------------------------------------------

    def save_to_file( self, file_name ):
        """
        Saves the fisher matrix to a file. Notice that the file name has to be specified
        to avoid overwriting an existing fisher matrix.

        :param file_name: file name and path of the output fisher matrix.
        The file extension gets automatically added as is not needed.

        """
        # save the param name file:
        self.save_paramnames_to_file( file_name=file_name+'.paramnames' )
        # open the output file:
        out_file = open( file_name+'.dat', 'w' )
        # write the header:
        out_file.write( '#\n' )
        out_file.write( '# This file contains a Fisher matrix created with the CosmicFish code.\n' )
        out_file.write( '#\n' )
        out_file.write( '# The parameters of this Fisher matrix are:\n' )
        out_file.write( '#\n' )
        # write the param names commented:
        for ind in range( self.num_params ):
            param_name = self.get_param_name( ind+1 )
            out_file.write( '#'+'    '+str(ind+1)+'    '+ \
                            str(param_name)+'    '+ \
                            str(self.get_param_name_latex( param_name ))+'    '+ \
                            str(self.get_fiducial( param_name ))+'\n'
                            )
        out_file.write( '#\n' )
        # write the fisher matrix:
        fisher_matrix = self.get_fisher_matrix()
        for i in range( self.num_params ):
            for j in range( self.num_params ):
                out_file.write( str( format(fisher_matrix[i,j],'.16E') )+'     ' )
            out_file.write( '\n' )


        # close the output file:
        out_file.close()

    # -----------------------------------------------------------------------------------

    def __add__(self, other):
        """
        Addition operator (+). Safeguarded agains adding Fisher matrices with different parameters
        and different fiducials.
        The addition will add parameters with the same name and append parameters with different names.
        Notice that if a parameter is in one of the two Fisher matrices but not in the other it will be
        assumed independent from the other.
        """
        # get the new param names and fiducial:
        param_names_new       = copy.deepcopy( self.param_names )
        param_names_latex_new = copy.deepcopy( self.param_names_latex )
        param_fiducial_new    = copy.deepcopy( self.param_fiducial )
        # check the parameters of the other Fisher matrix:
        for i in other.param_names:
            # if the parameter is in common check that they have the same fiducial value:
            if i in self.param_names_dict:
                ind1 = self.param_names_dict[i]-1
                ind2 = other.param_names_dict[i]-1
                if not np.allclose( self.param_fiducial[ind1], other.param_fiducial[ind2]):
                    raise ValueError('Error in addition: parameter '+str(i)+' has different fiducials: '+str(self.param_fiducial[ind1])+' and '+str(other.param_fiducial[ind2]))
            # otherwise we add them:
            else:
                ind = other.param_names_dict[i]-1
                param_names_new.append(i)
                param_names_latex_new.append(other.param_names_latex[ind])
                param_fiducial_new = np.append( param_fiducial_new, other.param_fiducial[ind])
        # initialize an empty matrix:
        num_param_new = len(param_names_new)
        new_matrix = np.zeros([num_param_new,num_param_new])
        # fill the new matrix:
        for i in range(num_param_new):
            for j in range(num_param_new):
                # get the new parameters name at each entry:
                x = param_names_new[i]
                y = param_names_new[j]
                # try to see if we have a corresponding entry on the first matrix:
                try:
                    x1 = self.param_names_dict[x]-1
                    y1 = self.param_names_dict[y]-1
                    fact_1 = self.fisher_matrix[x1,y1]
                except:
                    fact_1 = 0.0
                # try to see if we have a corresponding entry on the second matrix:
                try:
                    x1 = other.param_names_dict[x]-1
                    y1 = other.param_names_dict[y]-1
                    fact_2 = other.fisher_matrix[x1,y1]
                except:
                    fact_2 = 0.0
                # write the entrance of the new matrix:
                new_matrix[i,j] = fact_1 + fact_2
        # create the new Fisher matrix:
        fisher_new = fisher_matrix( fisher_matrix=new_matrix, param_names=param_names_new, param_names_latex=param_names_latex_new, fiducial=param_fiducial_new )
        fisher_new.name  = self.name + '_' + other.name
        fisher_new.path  = self.path
        fisher_new.indir = self.indir

        return fisher_new

    # -----------------------------------------------------------------------------------

    def __eq__(self, other):
        """
        Equality check operator (==). Ensures equality in all properties of the Fisher matrix.
        Notice that also name, path and indir are checked.
        """
        try:
            return_value = isinstance( other, fisher_matrix ) and \
                           np.allclose(self.fisher_cutoff, other.fisher_cutoff) and \
                           np.allclose(self.fisher_matrix, other.fisher_matrix) and \
                           self.path == other.path and \
                           self.name == other.name and \
                           self.indir == other.indir and \
                           self.num_params == other.num_params and \
                           np.allclose(self.fisher_eigenvalues, other.fisher_eigenvalues) and \
                           np.allclose(self.fisher_eigenvectors,other.fisher_eigenvectors) and \
                           np.allclose(self.fisher_matrix_inv, other.fisher_matrix_inv) and \
                           self.param_names == other.param_names and \
                           self.param_names_latex == other.param_names_latex and \
                           np.allclose( self.param_fiducial, other.param_fiducial ) and \
                           self.param_names_dict == other.param_names_dict
        except:
            return_value = False

        return return_value

    # -----------------------------------------------------------------------------------

    def __ne__(self, other):
        """
        Non-equality operator (!=). Simply implemented as the inverse of the equality operator.
        """
        return not self.__eq__(other)

    # -----------------------------------------------------------------------------------

    def inverse_fisher_matrix(self):
        """
        Invert the Fisher matrix.

        :returns: a matrix containing the inverse of the Fisher matrix.
        :rtype: :class:`numpy.array`

        """
        return np.linalg.inv( np.array( self.fisher_matrix ) )

    # -----------------------------------------------------------------------------------

    def PCA( self ):
        """
        This function performs the principal component analysis of the Fisher matrix returning its eigenvalues and its eigenvectors.
        As of now it just works just as a wrapper for numpy.

        :returns: a :class:`list` with (eigenvalues, eigenvectors) as :class:`numpy.array`.
        :rtype: :class:`list` of :class:`numpy.array`

        """
        w, v = np.linalg.eigh( self.fisher_matrix )
        return ( w, v )

    # -----------------------------------------------------------------------------------

    def determinant( self ):
        """
        This function returns the determinant of the Fisher matrix.

        :return: a :class:`float` with the determinant of the Fisher matrix.

        """
        return np.linalg.det( self.fisher_matrix )

    # -----------------------------------------------------------------------------------

    def protect_degenerate(self, cache=True):
        """
        Protects the Fisher matrix against degeneracies. Modifies the spectrum to ensure that
        the absolute value of the eigenvalues is bounded. This will make the Fisher matrix
        strictly positive definite. It will modify the magnitude of the worst constrained
        parameter combinations without modifying the degeneracies directions.

        :param cache: (optional) wether to use cached results or compute everything again
        :type cache: bool

        """
        # check cache:
        if not cache:
            (self.fisher_eigenvalues,self.fisher_eigenvectors) = self.PCA()
        # get the eigenvalues and eigenvactors:
        fisher_eigenvectors    = self.fisher_eigenvectors
        fisher_eigenvectors_m1 = np.transpose( np.array( self.fisher_eigenvectors ) )
        #
        redo_PCA = False
        # sometime the cutoff has to be adjusted if the condition number is too large:
        minimum_spectrum = np.amin( np.abs( self.fisher_eigenvalues ) )
        maximum_spectrum = np.amax( np.abs( self.fisher_eigenvalues ) )

        if np.isclose( minimum_spectrum, 0.0 ):
            condition_number = maximum_spectrum/self.fisher_cutoff
        else:
            condition_number = maximum_spectrum/minimum_spectrum

        if np.abs(condition_number) > self.fisher_spectrum:
            self.fisher_cutoff = maximum_spectrum/self.fisher_spectrum
        # cycle through the eigenvalues
        for i in range( len(self.fisher_eigenvalues) ):
            # detect very small eigenvalues
            if self.fisher_eigenvalues[i] < self.fisher_cutoff:
                # tell the code to redo PCA:
                redo_PCA = True
                # define the perturbation to the eigenvalue
                self.fisher_eigenvalues[i] = self.fisher_cutoff
        # modify the Fisher matrix if necessary:
        if redo_PCA:
            # create a zero matrix:
            temp = np.zeros(( len(self.fisher_eigenvalues), len(self.fisher_eigenvalues)), float )
            # write:
            for ind, eigh in enumerate( self.fisher_eigenvalues ):
                temp[ind,ind] = eigh
            # get the Fisher:
            self.fisher_matrix = np.dot( np.dot( fisher_eigenvectors, temp ), fisher_eigenvectors_m1 )
            # redo PCA:
            (self.fisher_eigenvalues,self.fisher_eigenvectors) = self.PCA()

    # -----------------------------------------------------------------------------------

    def get_confidence_bounds(self, confidence_level=0.68, cache=False ):
        """
        Computes the marginal 1D confidence bounds on the Fisher parameters

        :param confidence_level: (optional) C.L. of the bounds. Default 68%.
        :type confidence_level: :class:`float` in [0,1]
        :param cache: (optional) wether to use cached results or compute everything again
        :type cache: bool

        """
        # check input:
        if confidence_level<0.0 or confidence_level>1.0:
            raise ValueError('Invalid confidence level. Legal input is between 0 and 1.')
        # invert the Fisher matrix:
        if cache:
            fisher_matrix_inv_temp = self.fisher_matrix_inv
        else:
            fisher_matrix_inv_temp = self.inverse_fisher_matrix()
        # compute the coefficient that correspond to the desired confidence level:
        coefficient = fu.confidence_coefficient( confidence_level )
        # compute the result:
        return coefficient*np.sqrt( np.diagonal( fisher_matrix_inv_temp ) )

    # -----------------------------------------------------------------------------------

    # class setters:

    def set_fisher_matrix(self, fisher_matrix):
        """
        Function sets a new fisher matrix substituting the old one. Notice that
        this will reset parameter names, latex parameter names and fiducial values.

        :param fisher_matrix: :class:`numpy.array` containing the input Fisher matrix.
        :type fisher_matrix: :class:`numpy.array`

        """
        # initialize the Fisher matrix:
        temp_fisher_matrix = np.array(fisher_matrix)
        # protection against 1D Fisher matrices
        if np.ndim(temp_fisher_matrix) == 0:
            temp_fisher_matrix = np.array([[temp_fisher_matrix]])
        elif np.ndim(temp_fisher_matrix) == 1:
            temp_fisher_matrix = np.array([temp_fisher_matrix])
        # check if the Fisher matrix is symetric:
        if not np.allclose(temp_fisher_matrix, np.transpose(temp_fisher_matrix)):
            raise ValueError('The input Fisher matrix is not equal to its transpose')
        # accept the fisher matrix and start doing computations:
        self.fisher_matrix = np.array(temp_fisher_matrix)
        # get the number of parameters:
        self.num_params = self.fisher_matrix.shape[0]
        # load the parameter names:
        self.param_names       = [ 'p'+str(i+1) for i in range(self.num_params) ]
        self.param_names_latex = [ 'p'+str(i+1) for i in range(self.num_params) ]
        self.param_fiducial    = np.array( [0.0 for i in self.param_names] )
        # re-create a dictionary of param names:
        self.param_names_dict = {}
        for i in range(len(self.param_names)):
            self.param_names_dict[i+1] = self.param_names[i]
            self.param_names_dict[self.param_names[i]] = i+1
        # do PCA and store results:
        (self.fisher_eigenvalues,self.fisher_eigenvectors) = self.PCA()
        # protect against degenerate parameters:
        self.protect_degenerate( )
        # invert the Fisher matrix and store the result:
        self.fisher_matrix_inv = self.inverse_fisher_matrix()

    def set_param_names(self, param_names):
        """
        Function sets a new list of param names substituting the old one. Notice that
        latex parameter names will be reset.

        :param param_names: list containing the new parameter names.
        :type param_names: :class:`list` of :class:`string`

        """
        # check the input:
        if len(param_names) != self.num_params:
            raise ValueError('set_param_names: the input param_names has not '+str(self.num_params)+' elements.')
        # accept input:
        self.param_names = copy.deepcopy(param_names)
        # overwrite the latex parameter names:
        self.param_names_latex = copy.deepcopy(param_names)
        # re-create a dictionary of param names:
        self.param_names_dict = {}
        for i in range(len(self.param_names)):
            self.param_names_dict[i+1] = self.param_names[i]
            self.param_names_dict[self.param_names[i]] = i+1

    def set_param_names_latex(self, param_names_latex):
        """
        Function sets a new list of LaTeX param names substituting the old one.

        :param param_names_latex: list containing the new LaTeX parameter names.
        :type param_names: :class:`list` of :class:`string`

        """
        # check the input:
        if len(param_names_latex) != self.num_params:
            raise ValueError('set_param_names_latex: the input param_names_latex has not '+str(self.num_params)+' elements.')
        # accept input:
        self.param_names_latex = copy.deepcopy(param_names_latex)

    def set_fiducial(self, fiducial):
        """
        Function sets a new fiducial substituting the old one.

        :param fiducial: list containing the new fiducial.
        :type param_names: :class:`list` of :class:`float`

        """
        # check the input:
        if len(fiducial) != self.num_params:
            raise ValueError('set_fiducial: the input fiducial has not '+str(self.num_params)+' elements.')
        # accept input:
        self.param_fiducial = np.array( fiducial )

    # -----------------------------------------------------------------------------------

# ***************************************************************************************
