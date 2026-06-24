#!/bin/bash -l
#SBATCH --output=/scratch_tmp/users/%u/%j.out


module load cellranger/7.0.0-gcc-13.2.0

cellranger  count --id 5A  --transcriptome /rds/prj/id_iris/IRIS/IRIS_Raw_data/refdata-gex-GRCh38-2024-A/ --fastqs ./ --sample 5A --expect-cells 10000 
