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
.. module:: fisher_operations
   :platform: Unix
   :synopsis: Module that contains operations that can be performed on Fisher matrices.
              All of them are safeguarded against non-Fisher input.

.. moduleauthor:: Marco Raveri <mraveri@uchicago.edu> for the CosmicFish code.

"""

# ***************************************************************************************

import numpy as np
import fisher_matrix as fm
import math

# ***************************************************************************************

def eliminate_columns_rows( fisher_matrix, indexes ):
    """
    This function eliminates the row and columns corresponding to the given indexes
    from the Fisher matrix. It also deletes all the other informations like the names
    of the parameters. Notice that the index corresponding to the first parameter
    is zero.

    :param fisher_matrix: input Fisher matrix
    :type fisher_matrix: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`
    :param indexes: list of integers with the indexes to delete from the Fisher matrix
    :type indexes: :class:`list` of :class:`int`        
    :returns: A Fisher matrix with the columns and rows deleted
    :rtype: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`
    
    """
    # check validity of the input:
    if ( not isinstance(fisher_matrix, fm.fisher_matrix) ):
            raise ValueError('Error, input fisher_matrix is not a fisher_matrix')
    # write the param names:
    new_param_names       = []
    new_param_names_latex = []
    new_param_fiducial    = []
    for i in xrange( fisher_matrix.num_params ):
        if i not in indexes:
            new_param_names.append( fisher_matrix.param_names[i] )
            new_param_names_latex.append( fisher_matrix.param_names_latex[i] )
            new_param_fiducial.append( fisher_matrix.param_fiducial[i] )
    # write the Fisher matrix:
    fisher_temp = np.delete ( np.delete( fisher_matrix.fisher_matrix, indexes , 0 ), indexes , 1 )
    # initialize the new Fisher matrix:
    fisher_new = fm.fisher_matrix(fisher_matrix=fisher_temp, param_names=new_param_names, param_names_latex=new_param_names_latex, fiducial=new_param_fiducial )        
    fisher_new.name  = fisher_matrix.name + '_reduced'
    fisher_new.path  = fisher_matrix.path
    fisher_new.indir = fisher_matrix.indir

    return fisher_new

# ***************************************************************************************

def eliminate_parameters( fisher_matrix, names ):
    """
    This function eliminates the row and columns corresponding to the given parameter
    name from the Fisher matrix. It also deletes all the other informations like the names
    of the parameters.
      
    :param fisher_matrix: input Fisher matrix
    :type fisher_matrix: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`
    :param names: list of names of the parameters to delete from the Fisher matrix
    :type names: :class:`list` of :class:`string`        
    :returns: A Fisher matrix with the parameters deleted
    :rtype: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`    
    
    """
    # check validity of the input:
    if ( not isinstance(fisher_matrix, fm.fisher_matrix) ):
            raise ValueError('Error, input fisher_matrix is not a fisher_matrix')
    # get the indexes of the parameters:
    index_list = []
    for i in names:
        if not fisher_matrix.param_names_dict.has_key(i):
            raise ValueError('Error, parameter '+str(i)+' is not in a parameter of fisher_matrix')
        index_list.append(fisher_matrix.param_names_dict[i]-1)
    # elminate them from the list and return:
    return eliminate_columns_rows( fisher_matrix, index_list )

# ***************************************************************************************

def reshuffle( fisher_matrix, names ):
    """
    This function reshuffles a Fisher matrix. The new Fisher matrix will have the 
    parameters specified in names, in the order specified by names.
    Can be used to delete parameters, change their order or extract the Fisher 
    for some parameters without marginalizing over the others.
        
    :param fisher_matrix: input Fisher matrix
    :type fisher_matrix: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`
    :param names: list of names of the parameters that are desired in the output Fisher 
        matrix, in the desired order.
    :type names: :class:`list` of :class:`string`        
    :returns: A Fisher matrix with the new parameters
    :rtype: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`
    
    """
    # check validity of the input:
    if ( not isinstance(fisher_matrix, fm.fisher_matrix) ):
            raise ValueError('Error, input fisher_matrix is not a fisher_matrix')
    # check wether the names required are inside the Fisher matrix:
    for i in names:
        if not fisher_matrix.param_names_dict.has_key(i):
            raise ValueError('Error, parameter '+str(i)+' is not in a parameter of fisher_matrix')
    # get the new latex names and fiducial:
    new_param_names_latex = []
    new_param_fiducial    = []
    for i in names:
        ind = fisher_matrix.param_names_dict[i] -1 
        new_param_names_latex.append(fisher_matrix.param_names_latex[ind])
        new_param_fiducial.append(fisher_matrix.param_fiducial[ind])
    # initialize an empty matrix:
    num_param_new = len(names)
    new_matrix = np.zeros([num_param_new,num_param_new])    
    # fill the new matrix:
    for i in xrange(num_param_new):
        for j in xrange(num_param_new):
            # get the name:
            x = names[i]
            y = names[j]
            # get the parameter name:
            x1 = fisher_matrix.param_names_dict[x]-1
            y1 = fisher_matrix.param_names_dict[y]-1
            # get the entrance of the new matrix:
            new_matrix[i,j] = fisher_matrix.fisher_matrix[x1,y1]
    # create the new Fisher matrix:
    fisher_new = fm.fisher_matrix(fisher_matrix=new_matrix, param_names=names, param_names_latex=new_param_names_latex, fiducial=new_param_fiducial)        
    fisher_new.name  = fisher_matrix.name + '_reshuffled'
    fisher_new.path  = fisher_matrix.path
    fisher_new.indir = fisher_matrix.indir        
    
    return fisher_new    
    
# ***************************************************************************************

def marginalise( fisher_matrix, names ):
    """
    This function marginalises a Fisher matrix over all parameters but the ones in names. 
    The new Fisher matrix will have the parameters specified in names, in the order specified by names.
    The calculation is performed in the numerically stable way.
        
    :param fisher_matrix: input Fisher matrix
    :type fisher_matrix: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`
    :param names: list of names of the parameters of the output Fisher matrix,
        in the order that will appear in the output Fisher matrix. All other parameters 
        will be marginalized over.
    :type names: :class:`list` of :class:`string`        
    :returns: A Fisher matrix with the marginalized parameters
    :rtype: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`    
    
    """
    # check validity of the input:
    if ( not isinstance(fisher_matrix, fm.fisher_matrix) ):
            raise ValueError('Error, input fisher_matrix is not a fisher_matrix')
    # check wether the names required are inside the Fisher matrix:
    for i in names:
        if not fisher_matrix.param_names_dict.has_key(i):
            raise ValueError('Error, parameter '+str(i)+' is not in a parameter of fisher_matrix')
    # get the new latex names and fiducial:
    new_param_names_latex = []
    new_param_fiducial    = []
    for i in names:
        ind = fisher_matrix.param_names_dict[i] -1 
        new_param_names_latex.append(fisher_matrix.param_names_latex[ind])
        new_param_fiducial.append(fisher_matrix.param_fiducial[ind])
    # initialize an empty matrix:
    num_param_new = len(names)
    new_matrix = np.zeros([num_param_new,num_param_new])    
    # fill the new inverse matrix:
    for i in xrange(num_param_new):
        for j in xrange(num_param_new):
            # get the name:
            x = names[i]
            y = names[j]
            # get the parameter name:
            x1 = fisher_matrix.param_names_dict[x]-1
            y1 = fisher_matrix.param_names_dict[y]-1
            # get the entrance of the new matrix:
            new_matrix[i,j] = fisher_matrix.get_fisher_inverse()[x1,y1]
    
    fisher_temp = np.linalg.inv( new_matrix )
    # create the new Fisher matrix:
    fisher_new = fm.fisher_matrix(fisher_matrix=fisher_temp, param_names=names, param_names_latex=new_param_names_latex, fiducial=new_param_fiducial)        
    fisher_new.name  = fisher_matrix.name + '_marginal'
    fisher_new.path  = fisher_matrix.path
    fisher_new.indir = fisher_matrix.indir        
    
    return fisher_new     

# ***************************************************************************************

def marginalise_over( fisher_matrix, names ):
    """
    This function marginalises a Fisher matrix over the parameters in names. 
    The new Fisher matrix will not have the parameters specified in names.
    The calculation is performed in the numerically stable way.
        
    :param fisher_matrix: input Fisher matrix
    :type fisher_matrix: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`
    :param names: list of names of the parameters over which the Fisher will be marginalised.
    :type names: :class:`list` of :class:`string`        
    :returns: A Fisher matrix with the names parameters marginalized.
    :rtype: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`      
        
    """
    # check validity of the input:
    if ( not isinstance(fisher_matrix, fm.fisher_matrix) ):
            raise ValueError('Error, input fisher_matrix is not a fisher_matrix')
    # check wether the names required are inside the Fisher matrix:
    for i in names:
        if not fisher_matrix.param_names_dict.has_key(i):
            raise ValueError('Error, parameter '+str(i)+' is not in a parameter of fisher_matrix')
    # get the indexes:
    new_names = [ i for i in fisher_matrix.param_names if i not in names ]
        
    return marginalise( fisher_matrix, new_names )

# ***************************************************************************************

def information_gain( fisher_1, fisher_2, fisher_prior, units=math.log(2.0), stat=True ):
    """
    This function computes the Fisher approximation of Kullback-Leibler information gain.
    For the details of the formula we refer to the CosmicFish notes.
        
    :param fisher_1: first input Fisher matrix
    :type fisher_1: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`
    :param fisher_2: second input Fisher matrix
    :type fisher_2: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`
    :param fisher_prior: input Fisher matrix with the prior information.
    :type fisher_prior: :class:`cosmicfish_pylib.fisher_matrix.fisher_matrix`
    :param units: Units of information gain. Optional by default in Bits.
    :type units: :class:`float`   
    :param stat: wether to output the expected value and variance
    :type stat: :class:`logical`   
    :returns: a :class:`float` with the information gain.
    :rtype: :class:`float`   
    """
    info_gain = 0.0
    # first computations:
    F1p = fisher_1 + fisher_prior
    F2p = fisher_2 + fisher_prior
    # get common parameter names:
    param_names = [ name for name in F1p.get_param_names() if name in F2p.get_param_names() ]
    # reshuffle the second matrix:
    F1p = reshuffle( F1p, param_names )
    F2p = reshuffle( F2p, param_names )
    # define a dummy Fisher matrix with empty entrances and with the same parameters as the others:
    fisher_temp = fm.fisher_matrix( fisher_matrix=0.0*F2p.get_fisher_matrix(),
                                    param_names=F2p.get_param_names(),
                                    param_names_latex=F2p.get_param_names_latex(),
                                    fiducial=F2p.get_param_fiducial() )
    fisher_temp = fisher_2 + fisher_temp
    # the first term:
    info_gain = info_gain -math.log( F1p.determinant()/F2p.determinant() ) 
    info_gain = info_gain -F1p.get_fisher_matrix().shape[0]
    # the second trace term:
    info_gain = info_gain + np.trace( np.dot( F2p.get_fisher_inverse() , F1p.get_fisher_matrix() ) )
    # add additional term if statistical average over data is wanted
    if stat:
        # we break down the third term into two pieces:
        temp = np.dot( np.dot( np.dot( fisher_temp.get_fisher_matrix(), F2p.get_fisher_inverse() ),F1p.get_fisher_matrix() ), F2p.get_fisher_inverse() )
        temp = temp + np.dot( np.dot( temp,fisher_temp.get_fisher_matrix() ), F1p.get_fisher_inverse() )
        info_gain = info_gain + np.trace( temp )
        # compute variance:
        temp          = np.dot( temp, temp )
        info_variance = np.trace( temp )
    # output
    info_gain     = info_gain/2.0/units
    return info_gain

# ***************************************************************************************
    