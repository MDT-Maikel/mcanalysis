#!/bin/bash

##### basic setup
samp_init=1
samp_fin=10
nsamp=10

prod="th_incl"
mass=1000
gstar=0.1
#prod="ttz"


##### average cross section
tot_xsec=0
for irun in `seq $samp_init $samp_fin`; do

    run_name=run$irun"_mass"$mass"_g"$gstar
    xsec=`sed -n "/\<Integrated weight\>/p" /data/btag/events/signal/$prod/mg5/mass$mass"_g"$gstar/$run_name"_mg5_output.txt" | cut -d":" -f2`
    if [[ "$OSTYPE" = "darwin"* ]]
    then
        xsec=`sed -n "/[[:<:]]Integrated weight[[:>:]]/p" /data/btag/events/signal/$prod/mg5/mass$mass"_g"$gstar/$run_name"_mg5_output.txt" | cut -d":" -f2`
    fi

    #xsec=`sed -n "/\<Integrated weight\>/p" /data/btag/events/background/$prod/mg5/run$irun"_merge_output.txt" | cut -d":" -f2`
    #if [[ "$OSTYPE" = "darwin"* ]]
    #then
    #    xsec=`sed -n "/[[:<:]]Integrated weight[[:>:]]/p" /data/btag/events/background/$prod/mg5/run$irun"_merge_output.txt" | cut -d":" -f2`
    #fi

    xsec=`echo ${xsec} | sed -e 's/[eE]+*/*10^/g'`
    tot_xsec=`echo "$tot_xsec + $xsec" | bc -l`
done

average_xsec=`echo "$tot_xsec / $nsamp" | bc -l`

echo "process "$prod" averaged cross section: "$average_xsec

##### succeeded
echo "===== done ====="
exit
