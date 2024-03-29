#!/bin/sh

# 
# This script is used to generate ExternCMangler.H 
# files which are needed for multidimensional builds.
#

bootstrapdir=`dirname $0`
if ! test $bootstrapdir = "."; then
    echo "*****************************************************"
    echo "Error: you must run bootstrap from its own directory."
    echo "*****************************************************"
    exit 1
fi


smallbuildGrep=[A-Za-z]

grepfile=/tmp/temp.grep.$$
rm -f $grepfile
for g in $smallbuildGrep; do
  echo $g >> $grepfile  # Can't seem to pass grep a shell variable...
done

min_spacedim=1
max_spacedim=6

#
# Generate *_ExternCMangler.H files.
#
echo 'Generating *_ExternCMangler.H files'
dirlist=`find ./lib/src -maxdepth 2 -mindepth 2  -type d -name multidim | grep -f $grepfile`
for d in $dirlist; do
    if test -f ${d}/extern.list; then
        (cd $d
         cd ..
         lib_src_dir=`basename $PWD`
         mangle_header=${lib_src_dir}_ExternC_Mangler.H
         cd multidim
         rm -f ../$mangle_header
         ../../../util/multidim/mangle_externs.sh $lib_src_dir $min_spacedim $max_spacedim > ../$mangle_header
        )
    fi
done

rm $grepfile

set -x



