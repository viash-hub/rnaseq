#!/bin/bash

id="SRR6357070"
echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
    --strandedness unstranded \
    --bam $meta_resources_dir/test.paired_end.sorted.bam \
    --bedgraph_forward forward.bedgraph \
    --bedgraph_reverse reverse.bedgraph

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

# check whether output exists
[ ! -f "forward.bedgraph" ] && echo "File 'forward.bedgraph' does not exist!" && exit 1
[ ! -s "forward.bedgraph" ] && echo "File 'forward.bedgraph' is empty!" && exit 1
[ ! -f "reverse.bedgraph" ] && echo "File 'reverse.bedgraph' does not exist!" && exit 1
[ ! -s "reverse.bedgraph" ] && echo "File 'reverse.bedgraph' is empty!" && exit 1

echo "All tests succeeded!"
exit 0