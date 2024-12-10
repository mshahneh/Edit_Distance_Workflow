#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.batch_x = 150
params.batch_y = 100

process CalcProcess {
    errorStrategy 'ignore'
    // conda "${baseDir}/environment.yml"
    conda '/home/user/mambaforge/envs/wout_data_analysis'
    maxForks 90
    publishDir "/home/user/LabData/Reza/data/Wout/120924/nf_output/", mode: 'copy'
    // publishDir "/home/user/LabData/Reza/Code/Wout_data/nf_output", mode: 'copy'
    // time '15m'

    input:
    val batch_index

    output:
    file "res/*" optional true

    script:
    """
    mkdir res
    python $baseDir/pairs.py \
    '${baseDir}/data/mols.csv' \
    '${baseDir}/data/mols.pkl' \
    '${baseDir}/data/motifs_mols.pkl' \
    'res/res_${batch_index}.csv' \
    --batch_x ${params.batch_x} --batch_y ${params.batch_y} --batch_index ${batch_index-1}\
    --top_limit 4 \
    """
    // """
    // mkdir res
    // python $baseDir/pairs.py \
    // '${baseDir}/DataForReza/ALL_GNPS_cleaned_unique.csv' \
    //     'res/res_${batch_index}.csv' --batch_x ${params.batch_x} --batch_y ${params.batch_y} --batch_index ${batch_index-1} \
    //     --cached_dir  "/home/user/LabData/Reza/data/Wout/nf_output/res/"
    // """
}


workflow {
    // using batch_x and batch_y to define the number of batches
    def batches = Channel.from(1..params.batch_x*params.batch_y)
    all_csvs = CalcProcess(batches)
    all_excels = all_csvs.collectFile(name: '/home/user/LabData/Reza/data/Wout/120924/nf_output/res_nf.csv', newLine: false, keepHeader: true, skip: 1)
    // all_excels = all_csvs.collectFile(name: '/home/user/LabData/Reza/Code/ModiFinder_Network.csv', newLine: false, keepHeader: true, skip: 1)
}