#!/bin/sh
#
# script to change a madgraph parameter card entry


# print usage, if no argument is specified
if [ $# -lt 3 ] 
then
    echo "Usage: $0 <param card> <param> <new value>"
    exit 1
fi

# get the script parameters
param_card=$1
par=$2
new_val=$3

# replace the parameter and store in tmp_file
tmp=`sed -n "/\<"$par"\>/p" $param_card  | cut -d"#" -f1`
old_val=`echo $tmp | cut -d" " -f3`
if [ -z $old_val ]
then
	old_val=`echo $tmp | cut -d" " -f2`
fi
sed "s/\<"$old_val" # "$par"\>/"$new_val" # "$par"/" $param_card > tmp0

# replace original card with tmp_file
cp tmp0 $param_card
rm tmp0

#succeeded
exit
