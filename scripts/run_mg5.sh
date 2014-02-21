#!/bin/sh
#
# script to run a madgraph instance


# print usage, if no argument is specified
if [ $# -lt 3 ] 
then
    echo "Usage: $0 <process dir> <event dir> <run name>"
    exit 1
fi

# get the program parameters
process_dir=$1
event_dir=$2
run_name=$3
wkdir=`pwd`

# copy process dir into run dir
run_dir=$wkdir"/run_"$run_name
if [ -d $run_dir ]
then
    echo "Warning: $run_dir already exists"
    exit
fi
cp -r $process_dir $run_dir

# generate events with madgraph
$run_dir/bin/madevent < $wkdir/mg5_script.dat

# store the results in the event dir
if [ ! -d $event_dir ]
then
	mkdir $event_dir
fi
mv $run_dir/Events/run_01/run_01_tag_1_banner.txt $event_dir/$run_name.txt
mv $run_dir/Events/run_01/unweighted_events.lhe.gz $event_dir/$run_name.lhe.gz
mv $run_dir/Events/run_01/tag_1_delphes_events.lhco.gz $event_dir/$run_name.lhco.gz
mv $run_dir/HTML/run_01/results.html $event_dir/$run_name.html

# change the rights for the results
chmod g=rwx $event_dir/$run_name.txt
chmod g=rwx $event_dir/$run_name.lhe.gz
chmod g=rwx $event_dir/$run_name.lhco.gz
chmod g=rwx $event_dir/$run_name.html

# remove the run dir
rm -rf $run_dir

#succeeded
exit



