#----------------------------------------------------------------------------------------
#
# This file is part of CosmicFish.
#
# Copyright (C) 2015-2016 by the CosmicFish authors
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
# This file contains a script that packs the CosmicFish release.
#

# get the path of the script:
SCRIPT_PATH="`dirname \"$0\"`"                  # relative
SCRIPT_PATH="`( cd \"$SCRIPT_PATH\" && pwd )`"  # absolutized and normalized
if [ -z "$SCRIPT_PATH" ] ; then exit 1 ; fi     # error; for some reason, the path is not accessible to the script

# get the path of the user while it's executing the script:
USER_PATH="`dirname \"$1\"`"                  # relative
USER_PATH="`( cd \"$USER_PATH\" && pwd )`"    # absolutized and normalized
if [ -z "$USER_PATH" ] ; then exit 1 ; fi     # error; for some reason, the path is not accessible to the script

# go to the script directory:
cd $SCRIPT_PATH

# get colors:
source colors.sh

# some feedback:
echo -e ${BRed} 'Producing the CosmicFish distribution tar file' ${Color_Off}

# do the job:
cd ..
make distclean && \
#rsync -av --progress ./* ./CosmicFish --exclude CosmicFish --exclude .git && \
#tar cvzf --delete CosmicFish.tar.gz CosmicFish
cd .. && \
tar cvzf cosmicfish.tar.gz cosmicfish \
    --exclude=.git* \
    --exclude=*.DS_Store* \
    --exclude=.project

# if the code arrives here then everything is fine. Return:
exit 0