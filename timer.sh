#!/bin/bash

usage() {
  cat <<EOF
timer.sh
========
Timer script

Usage:
./timer.sh X "Y"
where the command Y is repeated X number of times.
Quotation marks are needed.
EOF
}

# Test if first parameter is a positive number, and that there's a command
# http://www.unix.com/shell-programming-scripting/21668-how-check-whether-string-number-not.html
if [ `echo $1 | sed 's/^[0-9][0-9]*//' | wc -c` -ne 1 ] || [ -z "$2" ]; then
  usage
  exit 0
fi

echo "Repeat the \"$2\" command exactly $1 times"
a=1
while [ $a -le "$1" ]
do
    echo " "
    echo "Repeat $a:"
    result=`time $2`
    a=$((a+1))
done

exit 0
