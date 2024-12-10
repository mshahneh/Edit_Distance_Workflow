run:
	nextflow run ./main.nf -resume -c nextflow.config

run_slurm:
	nextflow run ./main.nf -resume -c nextflow_slurm.config

run_docker:
	nextflow run ./main.nf -resume -with-docker <CONTAINER NAME>