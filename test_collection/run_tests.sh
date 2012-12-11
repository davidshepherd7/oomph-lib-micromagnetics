#!/bin/bash

# Set this to an exact directory if needed
test_folder=`pwd`

# Recursively find a list of folders containing Makefiles
dir_list=`find "$test_folder" -type f -name "Makefile.am" -print0 | xargs -0 -l1 dirname`

# In each folder run make then all excutables
for dir in $dir_list; do
    cd $test_folder
    cd $dir
    # echo `pwd`

    # make (with output recorded in a trace)
    if make > "make_trace"
    then
        # if there is a test.sh script run it
        if [ -e test.sh ]
        then
            if ./test.sh > test_trace
            then echo `basename $dir` "passed"
            else echo "*** $dir failed!"
            fi
        else
  	    # otherwise run all the executables in the folder (ignore backups)
	    exec_files=`find "./" -executable -type f | grep -v ".~/"`
	    for executable in $exec_files; do

	    # Record the output and report the result:
	        if ${executable} > "${executable}_trace"
	        then
		    echo "${executable} passed"
	        else
		    echo "*** ${executable} FAILED!"
	        fi
	    done
        fi
	# if make failed then say so
    else
	echo "*** Make FAILED in $dir!"
    fi
done