#!/bin/bash

id="SRR6357070"
echo ">>> Testing $meta_functionality_name"

"$meta_executable" \
    --id $id \
    --strandedness unstranded \
    --bam $meta_resources_dir/$id.sorted.bam \
    --bedgraph_forward $id.forward.bedgraph \
    --bedgraph_reverse $id.reverse.bedgraph

exit_code=$?
[[ $exit_code != 0 ]] && echo "Non zero exit code: $exit_code" && exit 1

# check whether output exists
[ ! -f $id.forward.bedgraph ] && echo "File '$id.forward.bedgraph' does not exist!" && exit 1
[ ! -f $id.reverse.bedgraph ] && echo "File '$id.reverse.bedgraph' does not exist!" && exit 1

echo "All tests succeeded!"
exit 0