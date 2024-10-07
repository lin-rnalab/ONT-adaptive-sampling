# ONT-adaptive-sampling

## Overview

This script is designed to gather information on the decision made during adaptive sampling for each read in a nanopore sequencing experiment, along with read length and quality information.

It will output a .tsv file with the read ID, decision, target type (i.e. on-target or off-target), length, Phred quality score, mapping quality score, and average alignment score per base for each read. This script will also generate a .tsv file containing summary statistics, and optionally a collection of plots to visualize the distribution of read lengths and quality scores for each subset of reads. 

The script can also be run on non-adaptive sampling experiments and will output ‘not_applicable’ in the decision field. 

## Table of Contents

* [Dependencies](#dependencies)
* [Pre-processing](#pre-processing)
* [Usage](#usage)
  + [Detect Differential Isoforms](#detect-differential-isoforms)
  + [Visualize Isoforms and Abundance](#visualize-isoforms-and-abundance)
  + [Classify Isoform Differences](#classify-isoform-differences)
  + [Example](#example)

## Dependencies

Dependencies can be installed by running [./install](./install). Alternatively, dependencies can be installed individually from conda.

samtools (v__)

pandas (v__)

matplotlib (v__)

numpy (v__)

pysam (v__)

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

process_sample.py [-h] \
    -t {AS,nAS} \
    -q /path/to/FASTQ/file \
    -b /path/to/BAM/file \
    -d /path/to/BED/file \
    -o /path/to/output/directory \
     [-n SAMPLE_NAME] \
     [-s /path/to/sequencing/summary/file] \
     [-p /path/to/adaptive/sampling/report] \
     [--split_bam_by_decision] \
     [--chunk_size CHUNK_SIZE] \
     [--output_plots] \
     [--min_length MIN_LENGTH] \
     [--max_length MAX_LENGTH] \
     [--report_supplementary] \
     [--report_secondary]

Arguments:

-t: Specify whether the experiment uses adaptive sampling (AS) or not (nAS).

-q: Path to the input FASTQ file.

-b: Path to the BAM file generated from alignment.

-d: Path to the BED file (optional for improved junction prioritization).

-o: Output directory to save results.

-n: (Optional) Sample name for labeling output.

-s: Path to the sequencing summary file. If experiment_type=AS, either -s or -p must be provided.

-p: Path to the adaptive sampling report. If experiment_type=AS, either -s or -p must be provided.

--split_bam_by_decision: (Optional) Output separate BAM files for each decision category.

--chunk_size: (Optional) Chunk size for reading in large files.

--output_plots: (Optional) Generate visualizations for read length and quality distributions.

--min_length/max_length: (Optional) Min/max read length for plots

--report_supplementary: (Optional) Report the maximum number of bases from reads with on-target supplementary alignments

--report_secondary: (Optional) Report the maximum number of bases from reads with on-target secondary alignments
