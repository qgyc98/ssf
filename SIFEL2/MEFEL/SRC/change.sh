#!/bin/bash

# execution path: *
# parameter: files (file.1 file.2, file.* , *.in *.cpp , *.* , ...)
# function: v aktualnim adresari prohleda vsechny specifikovane soubory
#           a vsechny retezce 'old' nahradi retezcem 'new'
# request: v retezci 'old' musi byt "regularni vyraz"
#          v retezci 'new' musi byt "obycejny vyraz"


# Set this strings, please :

old='\/\* adaptivity \*\/'
new='\/\* termitovo \*\/'


if [ -z $old ]; then echo "Set string 'old' within file 'change.sh' !!!"; exit 1; fi
if [ -z $@ ]; then echo "Set parameter, please !!!"; exit 1; fi


all=''
all1=''

for i in $@ ; do
    all1=`find . \! -type d -name $i`
    all="$all $all1"
done


echo 'checked files:'
echo $all
echo 'changed files:'

for file in $all; do
    sed s/"$old"/"$new"/g $file > aux.txt
    if `cmp -s $file aux.txt`
    then
	rm -f aux.txt
    else
        mv aux.txt $file
	echo $file
    fi
done


exit 0

