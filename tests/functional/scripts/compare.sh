#!/bin/bash

NAME=$1
shift
TESTS=$@

if [ "ALL" = "$TESTS" ]; then
	diff -rq results/${NAME} build/src/${NAME}_run
else
	for value in $TESTS; do
	  diff -rq results/${value}_run build/src/$value
  done
fi
