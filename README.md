# ONT-adaptive-sampling

## Overview

This script is designed to gather information on the decision made during adaptive sampling for each read in a nanopore sequencing experiment, along with read length and quality information.

It will output a `.tsv` file with the read ID, decision, target type (i.e. on-target or off-target), length, Phred quality score, mapping quality score, and average alignment score per base for each read. This script will also generate a `.tsv` file containing summary statistics, and optionally a collection of plots to visualize the distribution of read lengths and quality scores for each subset of reads.

The script can also be run on non-adaptive sampling experiments and will output ‘not_applicable’ in the decision field.

## Dependencies

- samtools
- pandas
- matplotlib
- numpy
- pysam

## Pre-processing

Basecalling should be performed on raw nanopore data (FAST5 or POD5 format) using ONT’s guppy or dorado software and output in FASTQ format.

Example:

    guppy_basecaller -i /path/to/raw/data -s /path/to/output/directory -c guppy_config_file

Or:

    dorado basecaller /path/to/model /path/to/raw/data --emit-fastq > /path/to/output/file

Alignment should be performed using minimap2 with the following command:

    minimap2 -ax splice -ub -t [threads] -k 14 -w 4 /path/to/reference/FASTA/file /path/to/input/FASTQ/file > /path/to/output/SAM/file

Optionally, users can provide a reference BED file to prioritize annotated junctions and aid in alignment:

    --junc-bed /path/to/reference/BED/file

Optionally, users can choose to exclude secondary alignments:

    --secondary=no

## Usage


    get_info.py [-h]
      -t {AS,nAS}
      -q /path/to/FASTQ/file
      -b /path/to/BAM/file
      -d /path/to/BED/file
      -o /path/to/output/directory
      [-n SAMPLE_NAME]
      [-s /path/to/sequencing/summary/file]
      [-p /path/to/adaptive/sampling/report]
      [--split_bam_by_decision]
      [--chunk_size CHUNK_SIZE]
      [--output_plots]
      [--min_length MIN_LENGTH]
      [--max_length MAX_LENGTH]


### Arguments:

- `-t`: Specify whether the experiment uses adaptive sampling (`AS`) or non-adaptive sampling (`nAS`).
- `-q`: Path to the input FASTQ file, basecalled and ideally untrimmed and unfiltered.
- `-b`: Path to the BAM file, aligned and ideally unfiltered.
- `-d`: Path to the BED file containing regions of interest (optional for improved junction prioritization).
- `-o`: Output directory where results will be stored.
- `-n`: (Optional) Sample name for labeling output files (default: "sample").
- `-s`: Path to the sequencing summary file. If `experiment_type=AS`, either `-s` or `-p` must be provided.
- `-p`: Path to the adaptive sampling report. If `experiment_type=AS`, either `-s` or `-p` must be provided.
- `--split_bam_by_decision`: (Optional) If set, outputs separate BAM files for each decision category.
- `--chunk_size`: (Optional) Chunk size for processing large files (default: 1,000,000).
- `--output_plots`: (Optional) If set, generates plots for read length, Phred quality, alignment score, and mapping quality by decision (default: True).
- `--min_length`: (Optional) Minimum read length for plots (default: 0).
- `--max_length`: (Optional) Maximum read length for plots (default: none).
