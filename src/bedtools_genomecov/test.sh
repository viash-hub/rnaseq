#!/bin/bash

id="SRR6357070"
echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
    --strandedness unstranded \
    --bam $meta_resources_dir/chr19.bam \
    --bedgraph_forward chr19_forward.bedgraph \
    --bedgraph_reverse chr19_reverse.bedgraph

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

# check whether output exists
[ ! -f "chr19_forward.bedgraph" ] && echo "File 'chr19_forward.bedgraph' does not exist!" && exit 1
[ ! -s "chr19_forward.bedgraph" ] && echo "File 'chr19_forward.bedgraph' is empty!" && exit 1
[ ! -f "chr19_reverse.bedgraph" ] && echo "File 'chr19_reverse.bedgraph' does not exist!" && exit 1
[ ! -s "chr19_reverse.bedgraph" ] && echo "File 'chr19_reverse.bedgraph' is empty!" && exit 1

echo "All tests succeeded!"
exit 0