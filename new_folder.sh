#!/bin/bash

# Script for all the actions required when creating a new folder in user_drivers.
# Input is the name of the new folder, a driver file will also be created with _driver.cc appended to the folder name.

# Exit if any errors occur
set -o errexit

# Exit if any variables are undefined
set -o nounset

# Input Options
######################################################################

# Change these to the path to your oomph-lib install and user_drivers directory name respectively
OOMPHPATH="$HOME/oomph-lib/trunk"
USERDIRNAME="david_shepherd"

# Change this to the folder containg the closest code to yours and the name of the driver code (without the extension!) respectively.
OLDFOLDER="$OOMPHPATH/user_drivers/$USERDIRNAME/one_d_micromag"
OLDDRIVER="one_d_micromag" # do not include _driver part



# Main code:
######################################################################
NEWFOLDER="$OOMPHPATH/user_drivers/$USERDIRNAME/$1"

# Create the folder (-p gives no error if already exists, creates parents if needed)
mkdir -p $NEWFOLDER

# Create sub-folders to hold results
mkdir -p $NEWFOLDER/RESLT $NEWFOLDER/RESLT_old

# Add to oomph-lib user_src directory list (if not already added)
# Check if already added (options: ignore regex characters, only match entire lines, quiet - no output and exit with 0 status if match is found)
if ! grep -Fxq "user_drivers/$USERDIRNAME/$1" $OOMPHPATH/config/configure.ac_scripts/user_drivers.dir_list
then
    echo "user_drivers/$USERDIRNAME/$1" >> $OOMPHPATH/config/configure.ac_scripts/user_drivers.dir_list
fi

# Copy makerun.sh and give it permisson to run as a script
cp $OOMPHPATH/user_drivers/$USERDIRNAME/makerun-general.sh $NEWFOLDER/makerun.sh
chmod +x $NEWFOLDER/makerun.sh

# Copy over closest existing code and rename
cp $OLDFOLDER/${OLDDRIVER}_driver.cc $NEWFOLDER/$1_driver.cc

# Create Makefile.am from file for old code
cp $OLDFOLDER/Makefile.am $NEWFOLDER/Makefile.am
# Find all references to the old driver name and replace with the new one
sed s/$OLDDRIVER/$1_driver/ $OLDFOLDER/Makefile.am > $NEWFOLDER/Makefile.am

# Run quickautogen.sh (or autogen.sh if does not exist) to generate Makefile
# (the script requires us to be in it's directory to run properly)
cd $OOMPHPATH
if [ -e quickautogen.sh ]; then
    $OOMPHPATH/quickautogen.sh
else
    $OOMPHPATH/autogen.sh --rebuild
fi

# Add new folder to git and push (cd back to new folder first)
cd $NEWFOLDER
git add $1_driver.cc Makefile.am
git commit --message="Added new folder: $1"
git push
