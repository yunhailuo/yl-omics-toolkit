#!/bin/bash

# Get the full directory name of the script no matter where it is being called from.
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-which-directory-it-is-stored-in/246128#246128
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

python $DIR/../../yl_sam_tools.py -s input.sam
