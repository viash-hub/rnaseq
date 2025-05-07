# rnaseq v0.2.1

## Minor changes

* Pin biobox version to v0.3.1

## Bug fixes

* Fix `summarizedexperiment` build PR (#42).
 
* Fix an issue with the `deseq2_qc` component not being able to create the DESeq2 object (PR #41).

## Known issues

The following caveats are known and will be addressed in future releases:

- [`bbmap_bbsplit` input file logic requires revision](https://github.com/viash-hub/rnaseq/issues/30)
- [Setting `--skip_deseq2_qc=false` results in an error](https://github.com/viash-hub/rnaseq/issues/31)
