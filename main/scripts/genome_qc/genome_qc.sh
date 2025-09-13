# requires set up of genome_qc
# https://github.com/hallamlab/genome_qc
# modified gtdb to use only 1 core for pplacer
# still required 110 GB RAM

GTDB=/path/to/gtdb_r226
cpus=14
apptainer exec -B ./ws:/app/ws,./data/gunc:/app/data/gunc,${GTDB}:/app/data/gtdb,./cache:/app/.snakemake,./cache2:/home/$USER/.cache docker://quay.io/hallamlab/genome_qc:1.0 \
    bash -c "cd /app && snakemake -p -s Snakefile.genome_qc --use-conda --cores ${cpus}"
