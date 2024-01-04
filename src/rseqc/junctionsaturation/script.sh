#!/bin/bash

set -eo pipefail 

prefix=$(openssl rand -hex 8)

junction_saturation.py \
    -i $par_input \
    -r $par_refgene \
    -o $prefix \
    -l $par_sampling_percentile_lower_bound \
    -u $par_sampling_percentile_upper_bound \
    -s $par_sampling_percentile_step \
    -m $par_min_intron \
    -v $par_min_splice_read \
    -q $par_map_qual

mv $prefix.junctionSaturation_plot.pdf $par_output_plot
mv $prefix.junctionSaturation_plot.r $par_output_plot_r

# Version
read -r -d '' text <<- END_VERSIONS
"${meta_functionality_name}":
    rseqc: \$(junction_saturation.py --version | sed -e "s/junction_saturation.py //g")
END_VERSIONS

if [ -e "$par_versions" ]; then
    echo "$text" >> "$par_versions"
else
    echo "$text" > "$par_versions"
fi