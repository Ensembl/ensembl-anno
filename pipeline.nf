#!/usr/bin/env nextflow


/*
See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/


nextflow.enable.dsl=2


/*
Usage:
   nextflow run pipeline.nf --run_simple_features --genome_file <genome_file> --output_dir <output_dir> --num_threads <number_of_threads>
*/


// workflow parameters
params.run_simple_features = false
params.output_dir = "."
params.genome_file = ""
params.num_threads = 1


// default workflow
workflow {
    // input data channel
    input_channel = Channel.fromPath(params.genome_file)

    run_anno(input_channel)

    // print process output to the terminal
    run_anno.out.view()
}


process run_anno {
    input:
    path genome_file

    output:
    stdout

    script:
    """
    python ensembl_anno.py --output_dir "${params.output_dir}" --genome_file "${genome_file}" --run_simple_features --num_threads "${params.num_threads}"
    """
}
