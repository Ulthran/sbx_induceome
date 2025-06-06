<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_induceome

<!-- Badges start -->
[![Tests](https://github.com/Ulthran/sbx_induceome/actions/workflows/tests.yml/badge.svg)](https://github.com/Ulthran/sbx_induceome/actions/workflows/tests.yml)
![Condabot](https://img.shields.io/badge/condabot-active-purple)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_induceome)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_induceome/)
<!-- Badges end -->

## Introduction

sbx_induceome is a [sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for identifying induced, (previously) integrated phages. This pipeline uses [bwa](https://bio-bwa.sourceforge.net/) for alignment to the host and [samtools](https://www.htslib.org/) to find the induced sites.

NOTE: sbx_induceome makes a couple assumptions about your reference genomes and reads. It assumes that 1) the reference files are named `cd{number}.fas` (under the directory you specify in the config) and 2) the read files are formatted like `...T{number}_....fastq.gz` (the `...` here can be whatever other characters you have in your sample name, although be careful not to include other instances of `T{number}_`).

## Options for config.yml

  - threads: Number of cores to use when running bwa
  - reference_fp: Directory with reference genomes
  - mapping_fp: Path to the file mapping samples to strains (.csv)
  - blast_db: Path to the viral blast db (full path including .faa file name)
  - phold_db: Path to the phold db (full path to the directory that will hold the tsv and other db files)
  - min_width: The sliding window size for peak finding
  - smoothing_factor: The number of below-threshold sliding windows to tolerate as being part of one peak
    
## Docs

More [docs](https://sunbeam.readthedocs.io/en/stable/extensions.html).
