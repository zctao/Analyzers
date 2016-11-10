#!/bin/bash
infile=${1:-./task_to_resubmit.txt}

while read line; 
do
	#echo "$line"
	crab resubmit -d "$line"
done < $infile