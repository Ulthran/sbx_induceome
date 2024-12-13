<img src="https://github.com/sunbeam-labs/sunbeam/blob/stable/docs/images/sunbeam_logo.gif" width=120, height=120 align="left" />

# sbx_induceome

<!-- Badges start -->
[![Tests](https://github.com/sunbeam-labs/sbx_induceome/actions/workflows/tests.yml/badge.svg)](https://github.com/sunbeam-labs/sbx_induceome/actions/workflows/tests.yml)
[![DockerHub](https://img.shields.io/docker/pulls/sunbeamlabs/sbx_induceome)](https://hub.docker.com/repository/docker/sunbeamlabs/sbx_induceome/)
<!-- Badges end -->

## Introduction

sbx_induceome is a [sunbeam](https://github.com/sunbeam-labs/sunbeam) extension for identifying induced, (previously) integrated phages. This pipeline uses [bwa](https://bio-bwa.sourceforge.net/) for alignment to the host and [samtools](https://www.htslib.org/) to find the induced sites.

## Installation

Extension install is as simple as passing the extension's URL on GitHub to `sunbeam extend`:

    sunbeam extend https://github.com/Ulthran/sbx_induceome

Any user-modifiable parameters specified in `config.yml` are automatically added on `sunbeam init`. If you're installing an extension in a project where you already have a config file, run the following to add the options for your newly added extension to your config (the `-i` flag means in-place config file modification; remove the `-i` flag to see the new config in stdout):

    sunbeam config update -i /path/to/project/sunbeam_config.yml

Installation instructions for older versions of Sunbeam are included at the end of this README.

## Running

To run an extension, simply run Sunbeam as usual with your extension's target rule specified:

    sunbeam run --profile /path/to/project/ all_induceome

NOTE: sbx_induceome makes a couple assumptions about your reference genomes and reads. It assumes that 1) the reference files are named `cd{number}.fas` (under the directory you specify in the config) and 2) the read files are formatted like `...T{number}_....fastq.gz` (the `...` here can be whatever other characters you have in your sample name, although be careful not to include other instances of `T{number}_`).

### Options for config.yml

  - threads: Number of cores to use when running bwa
  - reference_fp: Directory with reference genomes
    
## Installing an extension (legacy instructions for sunbeam <3.0)

Installing an extension is as simple as cloning (or moving) your extension directory into the sunbeam/extensions/ folder, installing requirements through Conda, and adding the new options to your existing configuration file: 

    git clone https://github.com/sunbeam-labs/sbx_induceome/ sunbeam/extensions/sbx_induceome
    cat sunbeam/extensions/sbx_induceome/config.yml >> sunbeam_config.yml

## Issues with pipeline

Please post any issues with this extension [here](https://github.com/Ulthran/sbx_induceome/issues).
