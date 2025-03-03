#!/bin/bash

# This script aimed to install user-defined objective function
# in "file.cpp" into desired "directory". 

#  

if [ "$1" = "" ] || [ "$2" = "" ]; then
    echo
    echo "Usage: optim_install.sh file.cpp directory"
    echo
    exit 1
fi
if [ ! -f "$1" ]; then
    echo "File $1 doesn't exist." ; exit 1
fi
if [ ! -d "$2" ]; then
    echo "Directory $2 doesn't exist." ; exit 1
fi

target="${2}/obj_funct.cpp"
if [ ! -f "$target" ]; then
    echo "Target file obj_funct.cpp is not in directory $2." ; exit 1
fi
tmpfile=`mktemp ${2}/temp.XXXXXX`
if [ ! -f $tmpfile ]; then
    echo "Temporary files cannot be created." ; exit 1
fi

mv $target $tmpfile 
cp $1 $target
cd $2 ; make -k ; cd -
mv $tmpfile $target 

echo "Everythink should be OK."

exit 0

