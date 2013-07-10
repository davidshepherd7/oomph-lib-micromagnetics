#!/bin/bash

set -o errexit
set -o nounset

# If no new commits since we last wrote this file (and file exists) then
# just exit to avoid rebuilds.
a=$(git log --since=`stat -c%Y git_version.h`)
if [[ $a == "" ]]; then
    if [[ -e git_version.h ]]; then
        exit 0
    fi
fi

# Generate a header file containing current git version hash and time
echo "" > git_version.h
echo "// AUTOMATICALLY GENERATED FILE, DO NOT EDIT" >> git_version.h
echo "#ifndef OOMPH_GIT_VERSION_H" >> git_version.h
echo "#define OOMPH_GIT_VERSION_H" >> git_version.h
echo "" >> git_version.h
echo "namespace GitVersion" >> git_version.h
echo "{" >> git_version.h
echo "static const char* const VERSION = \"$(git rev-parse HEAD)\";" >> git_version.h
echo "}" >> git_version.h
echo "" >> git_version.h
echo "#endif" >> git_version.h


# ??ds deal with systems that don't have git?
# ??ds use git describe to get a better description?
