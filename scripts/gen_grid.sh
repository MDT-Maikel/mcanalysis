#!/bin/bash
#
# generates a grid of mc data and possible analyzes it


# print usage, if no argument is specified
if [ $# -lt 2 ] 
then
    echo "Usage: $0 <process dir> <event dir> <result file>"
    exit 1
fi

# get the script parameters, specifying result file is optional
process_dir=$1
event_dir=$2
if [ $# -ge 3 ] 
then
	results_txt=$3
else
	results_txt=`pwd`"/results.txt"
fi

# initiate or clear the contents of the result_txt
> $results_txt

# construct loops over the independent parameters to form the grid
# example loop over an integer parameter
for pi in `seq 0 10 100`; do
	
	# example loop over an floating point parameter
	for pf in $(echo "for(pf=0.0; pf<=1.0; pf+=0.1) pf"|bc); do

		# create a unique run name depending on the parameters
		run_name="pi"$pi"_pf"$pf

		# copy process dir into a unique run dir, based on the run name
		run_dir=`pwd`"/run_"$run_name
		if [ -d $run_dir ]
		then
			echo "Warning: $run_dir already exists"
			exit
		fi
		cp -r $process_dir $run_dir

		# edit the parameter card using a C++ program based on ./calc_val
		# ./calc_val is used to calculate dependent parameters and to make
		# sure the correct formatting is inserted into the param card
		param_card=$run_dir/Cards/param_card.dat
		# edit first independent parameters, because of formatting
		read par_pi <<< $(./calc_val "PI" $pi $pf)
		./change_par.sh $param_card "PI" $par_pi
		read par_pf <<< $(./calc_val "PF" $pi $pf)	
		./change_par.sh $param_card "PF" $par_pf
		# then edit the dependent paramaters
		read par_dep <<< $(./calc_val "PD" $pi $pf)	
		./change_par.sh $param_card "PD" $par_dep
		
		# run mg5 with given process dir, event dir and event name 
		./run_mg5.sh $run_dir $event_dir $run_name
		
		# remove the run dir
		rm -rf $run_dir
	
		# read the cross section from the event sample
		read cross_section <<< $(./read_xsec.sh $event_dir"/"$run_name".txt")

		# possible run an C++ analysis on the produced event sample
		# in this case ./analysis takes an lhco file and prints
		# the efficiency to the screen
		read analysis_result <<< $(./analysis $event_dir$"/"run_name".lhco.gz")

		# add the results to a file
		echo $pi $pf $cross_section $analysis_result >> $results_txt

	# exit the floating point parameter loop 
	done
# exit the integer parameter loop
done

#succeeded
exit
