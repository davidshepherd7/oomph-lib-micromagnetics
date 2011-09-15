#!/bin/bash

# Script for all the actions required when creating a new folder in user_drivers.
# Input is the name of the new folder, a driver file will also be created with _driver.cc appended to the folder name.

# Input Options
######################################################################

# Change these to the path to your oomph-lib install and user_drivers directory name respectively
OOMPHPATH=$HOME/oomph-lib/trunk
USERDIRNAME=david_shepherd

# Change this to the folder containg the closest code to yours and the name of the driver code (without the extension!) respectively.
OLDFOLDER=$OOMPHPATH/user_drivers/$USERDIRNAME/demag
OLDDRIVER=one_d_micromag_driver



# Main code:
######################################################################
NEWFOLDER=$OOMPHPATH/user_drivers/$USERDIRNAME/$1

# Create the folder
mkdir $NEWFOLDER

# Create sub-folders to hold results
mkdir $NEWFOLDER/RESLT $NEWFOLDER/RESLT_old

# Add to oomph-lib user_src directory list (if not already added)
# Check if already added (options: ignore case, only match entire lines, quiet - no output and exit with 0 status if match is found)
if ! grep -Fxq "user_drivers/$USERDIRNAME/$1" $OOMPHPATH/config/configure.ac_scripts/user_drivers.dir_list
then
    #echo "user_drivers/$USERDIRNAME/$1" >> $OOMPHPATH/config/configure.ac_scripts/user_drivers.dir_list
    echo "Wrote: user_drivers/$USERDIRNAME/$1 into user_drivers list"
fi

# Copy makerun.sh and plot.sh ??ds should probably make these work more generally and add to path
cp $OOMPHPATH/user_drivers/$USERDIRNAME/demag/makerun.sh $OOMPHPATH/user_drivers/$USERDIRNAME/demag/plot.sh $NEWFOLDER

# Copy over closest existing code and rename
cp $OLDFOLDER/$OLDDRIVER.cc $NEWFOLDER/$1_driver.cc

# Create Makefile.am from file for old code
cp $OLDFOLDER/Makefile.am $NEWFOLDER/Makefile.am
# Find all references to the old driver name and replace with the new one
sed s/$OLDDRIVER/$1_driver/ $NEWFOLDER/Makefile.am > $NEWFOLDER/Makefile.am

# Run quickautogen.sh (or autogen.sh if does not exist) to generate Makefile
# (the script requires us to be in it's directory to run properly)
# cd $OOMPHPATH
# if [ -e quickautogen.sh ]; then
#     $OOMPHPATH/quickautogen.sh
# else
#     $OOMPHPATH/autogen.sh
# fi

# Add new folder to git and push (cd back to new folder first)
cd $NEWFOLDER
git add $1_driver.cc Makefile.am #--dry-run
git commit --message="Added new folder: $1" #--dry-run
git push
