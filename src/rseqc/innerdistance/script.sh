#!/bin/bash

set -exo pipefail 

prefix=$(openssl rand -hex 8)

if $par_paired; then
    inner_distance.py \
        -i $par_input \
        -r $par_refgene \
        -o $prefix \
        -k $par_sample_size \
        -l $par_lower_bound_size \
        -u $par_upper_bound_size \
        -s $par_step_size \
        -q $par_map_qual \
    > stdout.txt

    head -n 2 stdout.txt > $par_output_stats
    
    mv $prefix.inner_distance.txt $par_output_dist
    mv $prefix.inner_distance_plot.pdf $par_output_plot
    mv $prefix.inner_distance_plot.r $par_output_plot_r
    mv $prefix.inner_distance_freq.txt $par_output_freq

fi

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    rseqc: \$(inner_distance.py --version | sed -e "s/inner_distance.py //g")
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi