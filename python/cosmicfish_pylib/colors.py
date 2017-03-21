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
.. module:: colors
   :platform: Unix
   :synopsis: A collection of color utilities.

.. moduleauthor:: Marco Raveri <mraveri@uchicago.edu> for the CosmicFish code.
"""

# ***************************************************************************************

def nice_colors( num ):
    """
    This function returns a color from a colormap defined below according to the
    number entered.

    :param num: input number. Can be an integer or float. Notice 
        that the colormap contains only a small numbers of colors. Even if the input is a float
        the output will still be one of the few colors in the colormap.
    :type num: :class:`int` or :class:`float`
    :return: tuple of :class:`float` containing the three RGB coordinates of the color.
    :rtype: tuple
        
    """
    # default colormap
    colormap = { 0: (0.0, 0.75, 0.75),
                 1: (0.0, 0.0, 1.0),
                 2: (1.0, 0.0, 0.0),
                 3: (0.0, 0.5, 0.0),
                 4: (0.75, 0.75, 0),
                 5: (0.0, 0.0, 0.0),
                 6: (0.75, 0, 0.75),
               }
    colormap = { 0: (203.0/255.0, 15.0/255.0, 40.0/255.0),
                 1: (255.0/255.0, 165.0/255.0, 0.0),
                 2: (42.0/255.0, 46.0/255.0, 139.0/255.0),
                 3: (0.0/255.0, 153.0/255.0, 204.0/255.0),
                 4: (0.0/255.0, 221.0/255.0, 52.0/255.0),
                 5: (0.0, 0.0, 0.0),
                 6: (0.0, 0.75, 0.75),
               }
    index = int( round( num%7 ))
    return colormap[index]

# ***************************************************************************************

class bash_colors:
    """
    This class contains the necessary definitions to print to bash screen with colors.
    Sometimes it can be useful and nice!
    
    :ivar HEADER: ANSI color for light purple.
    :ivar OKBLUE: ANSI color for blue.
    :ivar OKGREEN: ANSI color for green.
    :ivar WARNING: ANSI color for yellow.
    :ivar FAIL: ANSI color for red.
    :ivar BOLD: ANSI code for bold text.
    :ivar UNDERLINE: ANSI code for underlined text.
    :ivar ENDC: ANSI code to restore the bash default.
    
    """
    
    # -----------------------------------------------------------------------------------
    
    HEADER    = '\033[95m' #: ANSI color for light purple.
    OKBLUE    = '\033[94m' #: ANSI color for blue.
    OKGREEN   = '\033[92m' #: ANSI color for green.
    WARNING   = '\033[93m' #: ANSI color for yellow.
    FAIL      = '\033[91m' #: ANSI color for red.
    BOLD      = '\033[1m'  #: ANSI code for bold text.
    UNDERLINE = '\033[4m'  #: ANSI code for underlined text.
    ENDC      = '\033[0m'  #: ANSI code to restore the bash default.
    
    # -----------------------------------------------------------------------------------
    
    def __init__(self):
        pass
    
    # -----------------------------------------------------------------------------------
    
    def header(self, string):
        """ 
        Function that returns a string that can be printed to bash in :class:`cosmicfish_pylib.colors.bash_colors.HEADER` color.
        
        :param string: input string.
        :type string: string
        :return: the input string with the relevant ANSI code at the beginning and at the end.
        :rtype: string 
           
        """
        return self.HEADER+str(string)+self.ENDC
    
    # -----------------------------------------------------------------------------------
    
    def blue(self,string):
        """ 
        Function that returns a string that can be printed to bash in :class:`cosmicfish_pylib.colors.bash_colors.OKBLUE` color.
        
        :param string: input string.
        :type string: string
        :return: the input string with the relevant ANSI code at the beginning and at the end.
        :rtype: string 
           
        """
        return self.OKBLUE+str(string)+self.ENDC
    
    # -----------------------------------------------------------------------------------
    
    def green(self,string):
        """ 
        Function that returns a string that can be printed to bash in :class:`cosmicfish_pylib.colors.bash_colors.OKGREEN` color.
        
        :param string: input string.
        :type string: string
        :return: the input string with the relevant ANSI code at the beginning and at the end.
        :rtype: string 
           
        """
        return self.OKGREEN+str(string)+self.ENDC
    
    # -----------------------------------------------------------------------------------
    
    def warning(self,string):
        """ 
        Function that returns a string that can be printed to bash in :class:`cosmicfish_pylib.colors.bash_colors.WARNING` color.
        
        :param string: input string.
        :type string: string
        :return: the input string with the relevant ANSI code at the beginning and at the end.
        :rtype: string 
           
        """
        return self.WARNING+str(string)+self.ENDC
    
    # -----------------------------------------------------------------------------------
    
    def fail(self,string):
        """ 
        Function that returns a string that can be printed to bash in :class:`cosmicfish_pylib.colors.bash_colors.FAIL` color.
        
        :param string: input string.
        :type string: string
        :return: the input string with the relevant ANSI code at the beginning and at the end.
        :rtype: string 
           
        """
        return self.FAIL+str(string)+self.ENDC
    
    # -----------------------------------------------------------------------------------
    
    def bold(self,string):
        """ 
        Function that returns a string that can be printed to bash in :class:`cosmicfish_pylib.colors.bash_colors.BOLD` color.
        
        :param string: input string.
        :type string: string
        :return: the input string with the relevant ANSI code at the beginning and at the end.
        :rtype: string 
           
        """
        return self.BOLD+str(string)+self.ENDC
    
    # -----------------------------------------------------------------------------------
    
    def underline(self,string):
        """ 
        Function that returns a string that can be printed to bash in :class:`cosmicfish_pylib.colors.bash_colors.UNDERLINE` color.
        
        :param string: input string.
        :type string: string
        :return: the input string with the relevant ANSI code at the beginning and at the end.
        :rtype: string
        
        """
        return self.UNDERLINE+str(string)+self.ENDC
    
    # -----------------------------------------------------------------------------------
    
# ***************************************************************************************
