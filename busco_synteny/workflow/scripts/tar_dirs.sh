#!/bin/bash
START=$(pwd)
echo $START
for d in $( find . -type d -name 'run_metazoa_odb10' ) ; do
    echo $d
    cd $START
    cd ${d}
    echo $PWD
    tar -czf metauk_output.tar.gz metauk_output --remove-files
    tar -czf hmmer_output.tar.gz hmmer_output --remove-files
    tar -czf busco_sequences.tar.gz busco_sequences --remove-files
done
