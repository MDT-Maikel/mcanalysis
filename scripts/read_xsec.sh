#!/bin/sh
#
# script to read the cross section from a madgraph file


# print usage, if no argument is specified
if [ $# -eq 0 ] 
then
    echo "Usage: $0 <file>"
    exit 1
fi

# get file and check whether it exists
xsec_file=$1
if [ ! -e $xsec_file ]; then
	echo "$0: Error opening file $3" 
	exit 2
fi   

# get the cross section using the identifier and strip useless characters.
xsec=$( awk '/^#  Integrated weight/ { $1=""; print $0; exit}' $xsec_file )
printf -v xsec '%s' $xsec
echo "$xsec" | tr -d [:alpha:] | tr -d '(' | tr -d ')' | tr -d ':'

#succeeded
exit
