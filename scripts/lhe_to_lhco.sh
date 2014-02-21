#!/bin/sh
#
# generates lhco data from lhe data using pythia and delphes


##### initialisation #####

# print warnings
if [[ $1 == "" || $2 == "" || $3 == "" ]]; then
    echo [lhe input] [lhco output] [delphes card] [optional: tmp folder]
    exit
fi

# set the directories
script_dir=`pwd`
mcanalysis_dir=$script_dir/..
build_dir=$mcanalysis_dir/build
delphes_dir=$SOFTDIR/Delphes-3.0.12
tmp_dir=$4
if [[ $4 == "" ]]; then
    tmp_dir=$mcanalysis_dir/tmp
fi
mkdir $tmp_dir

# set .lhe input
lhe_name=$1
lhe_input=$mcanalysis_dir/files/$lhe_name.lhe

# set .lhco output
lhco_name=$2
lhco_output=$mcanalysis_dir/files/$lhco_name.lhco.gz

# set Delphes card
delphes_card=$3

# compile gen_hepmc
cd $build_dir
make gen_hepmc -j4
cd $tools_dir


##### core processes #####

# run gen_hepmc to create .hepmc file from input .lhe
$build_dir/tools/gen_hepmc $lhe_input $tmp_dir/hepmc_output.hepmc
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
