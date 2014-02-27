#!/bin/sh
#
# generates lhco data from lhe data using pythia and delphes


##### initialisation #####

# print warnings
if [[ $3 == "" ]]; then
    echo [lhe input] [lhco output] [delphes card] [optional: --fast]"\n"[optional: tmp folder] [optional: merge process] [optional: merge njmax] [optional: merge scale] [optional: merge njadd]
    exit
fi

# set the directories
script_dir=`pwd`
mcanalysis_dir=$script_dir/..
build_dir=$mcanalysis_dir/build
delphes_dir=$SOFTDIR/Delphes-3.0.12
tmp_dir=$5
if [[ $5 == "" ]]; then
    tmp_dir=$mcanalysis_dir/tmp
fi
mkdir $tmp_dir

# set .lhe input
lhe_input=$1

# set .lhco output
lhco_output=$2.lhco.gz

# set Delphes card
delphes_card=$3

# fast showering options
opt_fast=""
if [[ $4 != "" ]]; then
	opt_fast="--fast"
fi

# merging options
opt_proc=""
if [[ $6 != "" ]]; then
	opt_proc="--merge_process="$6
fi

opt_njmax=""
if [[ $7 != "" ]]; then
	opt_njmax="--merge_njmax="$7
fi

opt_scale=""
if [[ $8 != "" ]]; then
	opt_scale="--merge_scale="$8
fi

opt_add=""
if [[ $9 != "" ]]; then
	opt_add="--merge_njadd="$9
fi

# compile gen_hepmc
cd $build_dir
make gen_hepmc -j4
cd $mcanalysis_dir/scripts


##### core processes #####

# run gen_hepmc to create .hepmc file from input .lhe
$build_dir/tools/gen_hepmc $opt_fast $opt_proc $opt_njmax $opt_scale $opt_njadd $lhe_input $tmp_dir/hepmc_output.hepmc
echo "\n===== hepmc file created =====\n"

# run Delphes to create .root file from input .hepmc
$delphes_dir/DelphesHepMC $delphes_card $tmp_dir/root_output.root $tmp_dir/hepmc_output.hepmc
echo "\n===== root file created =====\n"

# run root2lhco to create .lhco file from input .root
$delphes_dir/root2lhco $tmp_dir/root_output.root $tmp_dir/lhco_output_tmp.lhco

# remove the first line of the generated .lhco file
tail -n +2 $tmp_dir/lhco_output_tmp.lhco > $tmp_dir/lhco_output.lhco

# compress the .lhco file
gzip $tmp_dir/lhco_output.lhco
chmod g=rwx $tmp_dir/lhco_output.lhco.gz
echo "\n===== lhco.gz file created =====\n"

# copy .lhco.gz file in correct folder and remove tmp files
cp $tmp_dir/lhco_output.lhco.gz $lhco_output
rm -rf $tmp_dir
echo "===== end of script ====="

exit
