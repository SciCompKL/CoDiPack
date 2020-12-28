#
# CoDiPack, a Code Differentiation Package
#
# Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
# Homepage: http://www.scicomp.uni-kl.de
# Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
#
# Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
#
# This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
#
# CoDiPack is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# CoDiPack is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU General Public License for more details.
# You should have received a copy of the GNU
# General Public License along with CoDiPack.
# If not, see <http://www.gnu.org/licenses/>.
#
# Authors:
#  - SciComp, TU Kaiserslautern:
#     Max Sagebaum
#     Tim Albring
#     Johannes Blühdorn
#

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
usage: $0 [-u] [-n name] -b baseFile files ..."

    exit 0
}

baseFileName=
testName=
updateResult=
while getopts n:b:u opt
do
    case "$opt" in
      n)  testName=" $OPTARG";;
      b)  baseFileName="$OPTARG";;
      u)  updateResult=1;;
      \?)   # unknown flag
          Usage;;
    esac
done
shift `expr $OPTIND - 1`

if [ -z $baseFileName ];
then Usage;
fi;

fail=0
if [ -z $updateResult ];
then
  # arguments have been read now iterate over the files and compare them
  res=
  if [[ 0 == $# ]];
  then
      res+=" $failure no drivers run this test."
      fail=-1
  else
      while [ $# -gt 0 ]
      do
          TEST_NAME_PATTERN='_([^/]+)\.out'
          res+=" "
          [[ $1 =~ $TEST_NAME_PATTERN ]] &&
              res+=${BASH_REMATCH[1]}:
          build/compare.exe -t 1e-14 $baseFileName $1
          if [ $? -eq 0 ];
          then res+=$ok;
          else
              res+=$failure
              fail=-1
          fi;
          shift
      done
  fi;
  echo Test$testName:$res
else

  if [[ 0 == $# ]];
  then
      res+=" $failure No file for update."
      fail=-1
  else
      # just update the results file
      echo "Test$testName: updating $1 --> $baseFileName"
      mkdir -p $(dirname $baseFileName)
      cp $1 $baseFileName
  fi;
fi;

exit $fail
