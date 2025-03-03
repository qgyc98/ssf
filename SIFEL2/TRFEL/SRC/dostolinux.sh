#!/bin/bash



all=''
all1=`find . -name '*.cpp'`
all2=`find . -name '*.h'`
all3=`find . -name '*.in'`
all4=`find . -name '*.tex'`

all="$all1 $all2 $all3 $all4"

echo 'checked files:'
echo $all
echo 'changed files:'

for file in $all; do
    tmpfile=`mktemp` || exit 1
    cat $file | tr -d '\r' > $tmpfile
    if `cmp -s $file $tmpfile`
    then
	rm -f $tmpfile
    else
        mv $tmpfile $file
	echo $file
    fi
done