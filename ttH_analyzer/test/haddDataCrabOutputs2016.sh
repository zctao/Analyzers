#!/bin/bash

directory=${1:-/eos/uscms/store/user/ztao/ttH_80X/DoubleEG/}
datasets=${2:-[b,c,d]}
batch=${3:-jan2017}
date=${4:-170119}

cd $directory

sample=$(basename "$directory")

case "$sample" in
	DoubleEG) s="_dieg_";;
	DoubleMuon) s="_dimu_";;
	MuonEG) s="_mueg_";;
	SingleElectron) s="_e_";;
	SingleMuon) s="_mu_";;
esac

eval hadd ./"$batch"/output_data_obs.root crab_data_obs"$s"2016*/"$date"_*/*/output_data_obs_*.root
eval hadd ./"$batch"/output_fakes_data.root crab_fakes_data"$s"2016*/"$date"_*/*/output_fakes_data_*.root
eval hadd ./"$batch"/output_flips_data.root crab_flips_data"$s"2016*/"$date"_*/*/output_flips_data_*.root

cd -
