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

# get the cross section using the identifier and then take second argument after ':'
xsec=`sed -n "/\<Integrated weight\>/p" $xsec_file  | cut -d":" -f2`
echo $xsec

#succeeded
exit
