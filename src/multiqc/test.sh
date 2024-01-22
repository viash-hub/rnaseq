#!/bin/bash

unzip $meta_resources_dir/data.zip

echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
  --input "$meta_resources_dir/data/" \
  --report multiqc_report.html \
  --data multiqc_data \
  --plots multiqc_plots

echo ">>> Check whether output exists"
[ ! -f "multiqc_report.htm"l ] && echo "MultiQC report does not exist!" && exit 1
[ ! -s "multiqc_report.html" ] && echo "MultiQC report is empty!" && exit 1
[ ! -d "multiqc_data" ] && echo "MultiQC data directory does not exist!" && exit 1
[ -z "$(ls -A 'multiqc_data')" ] && echo "MultiQC data directory is empty!" && exit 1
[ ! -d "multiqc_plots" ] && echo "MultiQC plots directory does not exist!" && exit 1
[ -z "$(ls -A 'multiqc_plots')" ] && echo "MultiQC plots directory is empty!" && exit 1

echo "All tests succeeded!"
exit 0