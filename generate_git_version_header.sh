#!/bin/sh

set -o errexit
set -o nounset

# Generate a header file containing current git version hash and time
echo "" > git_version.h
echo "// AUTOMATICALLY GENERATED FILE, DO NO EDIT" >> git_version.h
echo "#ifndef OOMPH_GIT_VERSION_H" >> git_version.h
echo "#define OOMPH__GIT_VERSION_H" >> git_version.h
echo "" >> git_version.h
echo "namespace GitVersion" >> git_version.h
echo "{" >> git_version.h
echo "extern const char[] VERSION = \"$(git rev-parse HEAD)\";" >> git_version.h
echo "extern const char[] BUILD_TIME = \"$(date)\";" >> git_version.h
echo "}" >> git_version.h
echo "" >> git_version.h
echo "#endif" >> git_version.h


# ??ds deal with systems that don't have git?
# ??ds use git describe to get a better description?
