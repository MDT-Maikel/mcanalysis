#!/bin/sh
#
# script to run a madgraph instance


# print usage, if no argument is specified
if [ $# -lt 3 ] 
then
    echo "Usage: $0 <process dir> <event dir> <run name>"
    exit 1
fi

# get the script parameters
process_dir=$1
event_dir=$2
run_name=$3

# generate events with madgraph
cat <<- EOF > $process_dir/mg5_script
	generate_events
	delphes
	0
EOF
$process_dir/bin/madevent < $process_dir/mg5_script
rm $process_dir/mg5_script

# store the results in the event dir
if [ ! -d $event_dir ]
then
	mkdir $event_dir
fi
mv $process_dir/Events/run_01/run_01_tag_1_banner.txt $event_dir/$run_name.txt
mv $process_dir/Events/run_01/unweighted_events.lhe.gz $event_dir/$run_name.lhe.gz
mv $process_dir/Events/run_01/tag_1_delphes_events.lhco.gz $event_dir/$run_name.lhco.gz
mv $process_dir/HTML/run_01/results.html $event_dir/$run_name.html

# change the rights for the results
chmod g=rwx $event_dir/$run_name.txt
chmod g=rwx $event_dir/$run_name.lhe.gz
chmod g=rwx $event_dir/$run_name.lhco.gz
chmod g=rwx $event_dir/$run_name.html

#succeeded
exit



