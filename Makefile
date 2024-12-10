run:
	nextflow run ./main.nf -resume -c nextflow.config --batch_x 150 --batch_y 100 

run_slurm:
	nextflow run ./main.nf -resume -c nextflow_slurm.config

run_docker:
	nextflow run ./main.nf -resume -with-docker <CONTAINER NAME>

run_sample:
	nextflow run ./main.nf -resume -c nextflow.config --mols_csv 'data/mols_small.csv' --output_dir 'nf_output_sample/'