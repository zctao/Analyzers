#!/bin/bash

outfile=${1:-./task_to_purge.txt}

rm $outfile
touch $outfile

for dir in ~/nobackup/crab/*/; do
	echo "Checking job status: " $dir
	crab status -d "$dir" > tmp_p
	if (grep -q COMPLETED tmp_p); then
		echo $dir | cat >> "${outfile}"
	fi
	rm tmp_p
done