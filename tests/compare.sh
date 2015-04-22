#!/bin/bash

# compare.sh - compare files and dsiplays if they match
#
# Author: Max Sagebaum <max.sagebaum@scicomp.uni-kl.de>
# Date:   2015-04-22
# Category: File Comparison, TXT

# Colored output for 'ok' and 'failure'
ok="$(tput setaf 2)OK$(tput sgr 0)"
failure="$(tput setaf 1)FAILURE$(tput sgr 0)"

Usage () {
    echo >&2 "$0 - compare files with a base file and print out if they match
usage: $0 [-n name] -b baseFile files ..."

    exit 0
}

baseFileName=
testName=
while getopts n:b: opt
do
    case "$opt" in
      n)  testName=" $OPTARG";;
      b)  baseFileName="$OPTARG";;
      \?)   # unknown flag
          Usage;;
    esac
done
shift `expr $OPTIND - 1`

if [ -z $baseFileName ];
then Usage;
fi;

# arguments have been read now iterate over the files and compare them
res=
while [ $# -gt 0 ]
do
    res+=" "
    cmp $baseFileName $1
    if [ $? -eq 0 ];
    then res+=$ok;
    else res+=$failure;
    fi;
    shift
done
echo Test$testName:$res
