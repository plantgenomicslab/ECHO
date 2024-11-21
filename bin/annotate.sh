#!/usr/bin/env bash

set -x
set -e
#--rerun-triggers mtime \
snakemake -p \
	--rerun-incomplete \
	--use-singularity \
	--use-conda \
	--keep-going \
	--keep-incomplete \
	--stats annotation_runtime_stats.json \
	--cluster-config config_annotate.yml \
	--snakefile bin/Snakefile_annotate \
	--groups Sam2Transfrag=group0 --group-components group0=100 \
	--max-jobs-per-second 50 \
	--max-status-checks-per-second 50 \
	--jobs 150 \
	--latency-wait 30 \
	--notemp \
	--cluster "sbatch --mem={cluster.memory} --cpus-per-task={cluster.ncpus} \
			-J {cluster.name} \
			--parsable -A {cluster.account} -p {cluster.partition} \
			-t {cluster.time} -o {cluster.output} -e {cluster.error}" \
		"$@"
