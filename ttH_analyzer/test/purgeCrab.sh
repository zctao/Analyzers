#!/bin/bash
infile=${1:-./task_to_purge.txt}

while read line; 
do
	#echo "$line"
	crab purge -d "$line"
done < $infile