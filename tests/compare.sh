#!/bin/bash

ok="$(tput setaf 2)OK$(tput sgr 0)"
failure="$(tput setaf 1)FAILURE$(tput sgr 0)"
cmp $1 $2
if [ $? -eq 0 ];
then res1=$ok;
else res1=$failure;
fi;
cmp $1 $3
if [ $? -eq 0 ];
then res2=$ok;
else res2=$failure;
fi;
echo Test: $4 $res1 $res2
