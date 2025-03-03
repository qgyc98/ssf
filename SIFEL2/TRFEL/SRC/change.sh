#!/bin/bash

# execution path: *
# parameter: files (file.1 file.2, file.* , *.in *.cpp , *.* , ...)
# function: v aktualnim adresari prohleda vsechny specifikovane soubory
#           a vsechny retezce 'old' nahradi retezcem 'new'
# request: v retezci 'old' musi byt "regularni vyraz"
#          v retezci 'new' musi byt "obycejny vyraz"


# Set this strings, please :

old='transmition'
new='transmission'


[ "$old" = "" ] && echo "Set string 'old' within file 'change.sh' !!!" && exit 1
[ "$1"   = "" ] && echo "Set parameter, please !!!"                    && exit 1


all=''
all1=''

for i in $@ ; do
    i=${i#./}
    all1=`find . -type f -name "$i"`
    all="$all $all1"
done


printf "\n checked files:\n\n"
echo $all
printf "\n changed files:\n\n"

for file in $all; do
    [ "$file" = "*change.sh*" ] && continue
    
    cp $file $file~
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

