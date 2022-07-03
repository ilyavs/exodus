# Exodus
Exodus is a bioinformatics pipeline aimed to quantify relative abundance of highly similar sequences.
# Dependencies
The software is distributed as a self contained docker image biomx/exodus:1.0 on Docker Hub.
# Input data preparation
## Illumina short reads
The short reads need to have their adapters and low quality bases trimmed. Cutadapt is highly recommended with the appropriate adapter sequences. For quality trimming we recommend minimum quality score of 30 and minimum read length of 140.
Furthermore, when the sequenced sample has host DNA (for example host bacteria for phages) we recommend discarding reads that map to the host genome. One can run bowtie2 in very-fast-local mode and retain read pairs that did not align concordantly to the host genome.

## Reference set
The quality of the reference set is of the highest importance for successful analysis. The reference set must contain the references for all the genomes that are expected in the sample. It must not contain duplicate reference genomes (even if they are shifted in sequence). The reference sequences must be complete and not miss any parts of the reference genome. If the reference genome is circular, we recommend extending the reference by the read length (for example 250bp) to allow for read mapping to the area where the genome has been linearized. This will reduce false contamination by other closely related genomes but it will also have a minor effect on the relative abundance calculation (which accounts for genome size).

## Input config file
The pipeline is implemented as a snakemake pipeline executed in a docker container. The input for the pipeline is a config json or yaml file containing samples with their R1 and R2 paths and references with the path to the fasta file. Please see the sample config file in this repository for the format.

# Usage
To run the pipeline, place the config file in a new working directory and change to that directory. The pipeline will be executed within a docker container and therefore the paths containing all the input files specified in the config file must be mounted to the docker container. For example, if the files are in `/home/ilyav` then make sure to mount this path when executing the pipeline using `-v /home/ilyav:/home/ilyav`. An example command when you are in the working directory with the config file:
```docker container run --rm --init -v $(pwd):/workdir --env host_path_to_workdir=$(pwd) -v /home/ilyav:/home/ilyav -v /var/run/docker.sock:/var/run/docker.sock --workdir /workdir biomx/exodus:1.0 snakemake -s /exodus/exodus.smk --cores all --configfile /workdir/config.yaml -d /workdir```

# Testing
To run the test that is included in this repository, change to a working directory for saving output files and run:
```docker container run --rm --init -v $(pwd):/workdir --env host_path_to_workdir=$(pwd) -v /home/ilyav:/home/ilyav -v /var/run/docker.sock:/var/run/docker.sock --workdir /workdir biomx/exodus:1.0 snakemake -s /exodus/exodus.smk --cores all --configfile /exodus/config.yaml -d /workdir```
The relative abundances are in `Summary/Combined_relative_abundances.csv` and should match the results in `test_data/Combined_relative_abundances.csv` in this repository.

# License
The software is distributed under the GPLv3 license.

# Citation
If using exodus, please cite:
Vainberg-Slutskin, I. et al. Exodus: sequencing-based pipeline for quantification of pooled variants. Bioinformatics 38, 3288–3290 (2022).
