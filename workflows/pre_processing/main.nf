workflowDir = params.rootDir + "/workflows"
targetDir = params.rootDir + "/target/nextflow"

openPipelineDir = params.rootDir + "/target/dependencies/git/github.com/openpipelines-bio/openpipeline/0.9.0/nextflow"

include { fastqc } from openPipelineDir + "/qc/fastqc/main.nf"

include { readConfig; helpMessage; readCsv; preprocessInputs; channelFromParams } from workflowDir + "/utils/WorkflowHelper.nf"
include {  setWorkflowArguments; getWorkflowArguments; passthroughMap as pmap; passthroughFlatMap as pFlatMap; strictMap as smap } from workflowDir + "/utils/DataflowHelper.nf"

config = readConfig("$workflowDir/pre_processing/config.vsh.yaml")

workflow {
  helpMessage(config)

  channelFromParams(params, config)
    | run_wf

}

workflow run_wf {
  take:
    input_ch

  main:
    output_ch = input_ch
      | view { "Input: $it" }
      | fastqc.run( auto: [ publish: true ] )
      | view { "Output: $it" }
      

  emit:
    output_ch
}

// ===============================
// === start of test workflows ===
// ===============================

workflow test_wf {
  helpMessage(config)

  // allow changing the resources_test dir
  params.resources_test = params.rootDir + "/resources_test"

  // or when running from s3: params.resources_test = "s3://openpipelines-data/"
  testParams = [
    param_list: [
    ]
  ]

  output_ch =
    channelFromParams(testParams, config)
      | view { "Input: $it" }
      | run_wf
      | view { output ->
        assert output.size() == 2 : "outputs should contain two elements; [id, file]"
        // ...
        "Output: $output"
      }
      | toSortedList()
      | map { output_list ->
        assert output_list.size() == 1 : "output channel should contain one event"
        // ...
      }
  
}
