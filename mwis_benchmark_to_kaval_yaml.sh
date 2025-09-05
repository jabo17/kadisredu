#!/bin/bash

# input
instances=$1 # "mwis_benchmark_instances_weighted_parhip.txt"

while read instance; do
	instance_name=$(basename $instance)
	instance_name="${instance_name%.*}"
	echo "  - generator: dummy"
	echo "    name: ${instance_name}"
	echo "    nokey: $instance"
	#echo "    kagen_option_string: \"file;filename=$instance;distribution=balance-edges\""
done < "$instances"
