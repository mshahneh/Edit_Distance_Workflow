#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.batch_x = 20
params.batch_y = 20
params.calculate_mces = 0
params.calculate_motif_based = 0
params.mols_csv = "$baseDir/data/mols.csv"
params.motifs_csv = "$baseDir/data/motifs.csv"
params.output_dir = 'nf_output'

// /home/user/LabData/Reza/data/Wout/120924/nf_output/

process CacheMols{
    conda "$baseDir/environment.yml"
    // conda '/home/user/mambaforge/envs/wout_data_analysis'
    input:
    file mols_csv
    file motifs_csv

    output:
    file 'cahced_output/mols.pkl'
    file 'cahced_output/motifs_mols.pkl'

    script:
    """
    mkdir cahced_output
    python $baseDir/create_cache_files.py $mols_csv $motifs_csv cahced_output
    """
}

process CalcProcess {
    errorStrategy 'ignore'
    conda "${baseDir}/environment.yml"
    // conda '/home/user/mambaforge/envs/wout_data_analysis'
    // maxForks 60
    publishDir "${params.output_dir}", mode: 'copy'
    // time '15m' 

    input:
    val batch_index
    file mols_csv
    file mols_pkl
    file motifs_mols_pkl
    val argument_string

    output:
    file "res/*" optional true

    script:
    """
    mkdir res
    python $baseDir/main_edit_distance.py \
    $mols_csv \
    $mols_pkl \
    $motifs_mols_pkl \
    'res/res_${batch_index}.csv' \
    --batch_x ${params.batch_x} --batch_y ${params.batch_y} --batch_index ${batch_index-1}\
    --cached_dir  "${params.output_dir}/res/" \
    $argument_string \
    --top_limit 4 \
    """
}


workflow {
    def mols = Channel.fromPath(params.mols_csv)
    def motifs = Channel.fromPath(params.motifs_csv)
    (cached_mols, cached_motifs) = CacheMols(mols, motifs)
    
    // convert to value channel
    cached_mols = cached_mols.collect()
    cached_motifs = cached_motifs.collect()
    mols = mols.collect()

    // using batch_x and batch_y to define the number of batches
    def batches = Channel.from(1..params.batch_x*params.batch_y)

    // create the boolean argument string
    argument_string = ""
    if (params.calculate_mces == 1 || params.calculate_mces == "1")
    {
        argument_string = argument_string + " --calculate_mces"
    }
    if (params.calculate_motif_based == 1 || params.calculate_motif_based == "1")
    {
        argument_string = argument_string + " --calculate_motif_based"
    }
    
    all_csvs = CalcProcess(batches, mols, cached_mols, cached_motifs, argument_string)
    
    // collect all the csvs to a single file
    all_excels = all_csvs.collectFile(name: "${params.output_dir}/res_nf.csv", newLine: false, keepHeader: true, skip: 1)
}