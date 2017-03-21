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

#
# This file contains a script that compiles in several ways the CosmicFish library to
# test wether it can be safely deployed, at least in the tested cases.
#

# get the path of the script:
SCRIPT_PATH="`dirname \"$0\"`"                  # relative
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$SCRIPT_PATH" ] ; then exit 1 ; fi     # error; for some reason, the path is not accessible to the script

# get the path of the user while it's executing the script:
USER_PATH="`dirname \"$1\"`"                    # relative
USER_PATH="`( cd \"$USER_PATH\" && pwd )`"      # absolutized and normalized
if [ -z "$USER_PATH" ] ; then exit 1 ; fi       # error; for some reason, the path is not accessible to the script

# go to the script directory:
cd $SCRIPT_PATH

# go to the root of the CosmicFish dir:
cd ..

# start the test:
make distclean && \
make all HAS_IFORT=0 DEBUG=0 && \
make distclean && \
make all HAS_IFORT=0 DEBUG=1 && \
make distclean && \
make all HAS_IFORT=1 DEBUG=0 && \
make distclean && \
make all HAS_IFORT=1 DEBUG=1 \
make distclean

# if the code arrives here then everything is fine. Return:
exit 0