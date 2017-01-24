#!/bin/bash

infile=${1:-./task_to_resubmit.txt}
outfile=${2:-./task_to_resubmit_updated.txt}

rm $outfile
touch $outfile

while read dir; do
	echo "Checking job status: " $dir
	crab status -d "$dir" > tmp
	if !(grep -q COMPLETED tmp); then
		echo $dir | cat >> "${outfile}"
	fi
	rm tmp
done <${infile}