#!/usr/bin/sh

while [ True ]
do
    output_var=$(ls /mnt/heintz-bambi1/processing/RNASEQ_jobs | grep -E '^[0-9].*.txt$')
    if [ -n "$output_var" ]; then
        # parse jobs and feed to pipeline
        nohup pipelineParser $output_var > `date +%Y-%m-%d`.$$.out 2>&1 &
    fi
    sleep 12h
done

