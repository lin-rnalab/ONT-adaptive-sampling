#!/usr/bin/env python3

# Author: Nicole DeBruyne (Lin Lab)
# Date: 2024.10.07

# This script processes data from a nanopore sequencing experiment (either adaptive sampling or non-adaptive sampling) 
# to obtain the final decision made during the adaptive sampling process ('accepted', 'rejected', or 'no_decision' for adaptive sampling; 'not_applicable' for non-adaptive sampling),
# target type (on-target, off-target, or unmapped), read length, Phred quality score, alignment score normalized by read length, and mapping quality score for each read.

# Optionally, the script can also:
# 1) split the input BAM file into separate BAM files for each decision category
# 2) generate plots of read lengths and quality scores by decision and target type

# The script requires the following input:
# - Experiment type (AS for adaptive sampling, nAS for non-adaptive sampling)
# - FASTQ file containing basecalled reads (ideally untrimmed and unfiltered)
# - BAM file containing aligned reads (ideally unfiltered)
# - BED file containing regions of interest
# - Sequencing summary file OR adaptive sampling report file (for adaptive sampling experiments only)
# - Output directory to store results

import argparse
import os
import pandas as pd
from collections import defaultdict
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import pysam
import math
import matplotlib.gridspec as gridspec
from matplotlib import rcParams

def parse_arguments():
    '''Parse command-line arguments'''

    parser = argparse.ArgumentParser(description="Obtain read depth, read lengths, and quality/alignment scores"
                                                 "for an adaptive sampling (AS) or non-adaptive sampling (nAS) experiment")
    parser.add_argument('-t', '--experiment_type', required=True, choices=['AS', 'nAS'],
                        help='Type of experiment ("AS" for adaptive sampling, "nAS" for non-adaptive sampling)')
    parser.add_argument('-q', '--fastq', required=True, help='Sample FASTQ file, basecalled and (ideally) untrimmed and unfiltered.')
    parser.add_argument('-b', '--bam', required=True, help='Sample BAM file, aligned and (ideally) unfiltered.')
    parser.add_argument('-d', '--bed', required=True, help='BED file containing regions of interest')
    parser.add_argument('-n', '--sample_name', default="sample", help='Name of the sample')   
    parser.add_argument('-s', '--seqsum', help='Sequencing summary file for AS experiments. Either --seqsum or --AS_report is required for AS experiments. Do not provide for nAS experiments.')
    parser.add_argument('-p', '--AS_report', help='Adaptive sampling report file for AS experiments. Either --seqsum or --AS_report is required for AS experiments. Do not provide for nAS experiments.')
    parser.add_argument('-o', '--output-directory', required=True, help='Output directory to store results')
    parser.add_argument('--split_bam_by_decision', action='store_true', help='Output BAM files for each decision category')
    parser.add_argument('--chunk_size', type=int, default=1000000, help='Chunk size for processing large files (default: 1,000,000)')
    parser.add_argument('--output_plots', action='store_true', help='Plot read length, Phred quality, alignment score, and mapping quality by decision (default: True)')
    parser.add_argument('--min_length', type=int, default=0, help='Minimum read length for plots (default: 0)')
    parser.add_argument('--max_length', type=int, help='Maximum read length for plots (default: none)')
    args = parser.parse_args()

    # Check experiment type. AS = adaptive sampling, nAS = non-adaptive sampling 
    # If experiment_type is AS, ensure either sequencing_summary.txt or adaptive_sampling.csv is provided, but not both.
    # If experiment_type is nAS, ensure sequencing_summary.txt and adaptive_sampling.csv are NOT provided, and split_bam_by_decision is NOT specified.
    if args.experiment_type == "AS":
        if not (args.seqsum and os.path.isfile(args.seqsum)) and not (args.AS_report and os.path.isfile(args.AS_report)):
            parser.error("Either '-s' ('--seqsum') or '-p' ('--AS_report') option must be provided for experiment_type 'AS'.")
        elif (args.seqsum and os.path.isfile(args.seqsum)) and (args.AS_report and os.path.isfile(args.AS_report)):
            parser.error("Conflicting options. Provide either '-s' ('--seqsum') or '-p' ('--AS_report'), but not both.")
    elif args.experiment_type == "nAS":
        if args.seqsum or args.AS_report:
            parser.error("Options '-s' ('--seqsum') and '-p' ('--AS_report') should not be provided for experiment_type 'nAS'.")
        if args.split_bam_by_decision:
            parser.error("Option '--split_bam_by_decision' should not be provided for experiment_type 'nAS'.")
    print(f"Arguments parsed successfully.\n")

    # Print arguments for user to review
    if args.sample_name:
        print(f"Sample name: {args.sample_name}")
    print(f"Experiment type: {args.experiment_type}")
    print(f"FASTQ file: {args.fastq}")
    print(f"BAM file: {args.bam}")
    print(f"BED file: {args.bed}")
    print(f"Output directory: {args.output_directory}")
    if args.seqsum:
        print(f"Sequencing summary file: {args.seqsum}")
    if args.AS_report:
        print(f"Adaptive sampling report file: {args.AS_report}")
    if args.split_bam_by_decision:
        print("Output BAM files will be generated for each decision category.")
    if args.output_plots:
        print("Read length, Phred quality, alignment score, and mapping quality plots will be generated.")
        if args.max_length:
            print(f"Maximum read length for plots: {args.max_length}")
        if args.min_length:
            print(f"Minimum read length for plots: {args.min_length}")
    return args

def get_decision(args):
    '''Get the end decisions made for each read from the sequencing summary or adaptive sampling report'''

    # Initialize sets to store read IDs for each decision category
    accepted_readIDs = set()
    rejected_readIDs = set()
    nodecision_readIDs = set()
    # Read in the appropriate file and add read IDs to the respective sets
    if args.seqsum:
        print("Collecting accepted and rejected read IDs from sequencing summary...")
        with open(args.seqsum, 'r') as f:
            header = f.readline().strip().split('\t')
            read_id_index = header.index('parent_read_id')
            for line in f:
                if 'signal_positive' in line:
                    read_id = line.strip().split('\t')[read_id_index]
                    accepted_readIDs.add(read_id)
                elif 'data_service_unblock_mux_change' in line:
                    read_id = line.strip().split('\t')[read_id_index]
                    rejected_readIDs.add(read_id)
    elif args.AS_report:
        print("Collecting accepted, rejected, and nodecision read IDs from adaptive sampling report...")
        with open(args.AS_report, 'r') as f:
            header = f.readline().strip().split('\t')
            read_id_index = header.index('read_id')
            for line in f:
                if 'stop_receiving' in line:
                    read_id = line.strip().split('\t')[read_id_index]
                    accepted_readIDs.add(read_id)
                elif 'unblock' in line:
                    read_id = line.strip().split('\t')[read_id_index]
                    rejected_readIDs.add(read_id)
                elif 'no_decision' in line:
                    read_id = line.strip().split('\t')[read_id_index]
                    nodecision_readIDs.add(read_id)
    return accepted_readIDs, rejected_readIDs, nodecision_readIDs

def get_target_type(args):
    '''Get on-target and off-target read IDs from a BAM file based on regions defined in a BED file'''
    
    # Initialize sets to store read IDs for on-target and off-target reads
    all_mapped_readIDs = set()
    on_target_readIDs = set()
    off_target_readIDs = set()

    try:
         # Use samtools to fetch all mapped, primary reads
        print("Collecting all mapped read IDs with 'samtools view -F 3844'...")
        # Process the output line by line (to reduce memory) and add read IDs to the set
        all_mapped_cmd = f"samtools view -F 3844 {args.bam_file}"
        with subprocess.Popen(all_mapped_cmd, shell=True, stdout=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                if line:
                    all_mapped_readIDs.add(line.split()[0])
        # Use samtools to fetch reads whose primary alignments overlap with regions in the BED file
        print("Collecting (primary) on-target read IDs with 'samtools view -F 3844 -L BED'...")
        # Process the output line by line (to reduce memory) and add read IDs to the set
        on_target_cmd = f"samtools view -F 3844 {args.bam_file} -L {args.bed_file}"
        with subprocess.Popen(on_target_cmd, shell=True, stdout=subprocess.PIPE, text=True) as proc:
            for line in proc.stdout:
                if line:
                    on_target_readIDs.add(line.split()[0])
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running samtools: {e}")
        return None, None
    
    # Calculate off-target read IDs
    print("Collecting (primary) off-target read IDs as (mapped read IDs - on-target read IDs)...")
    off_target_readIDs = all_mapped_readIDs - on_target_readIDs

    return on_target_readIDs, off_target_readIDs

def split_bam(args, accepted_readIDs, rejected_readIDs, nodecision_readIDs):
    '''Split the BAM file by decision category'''
    
    print("Splitting BAM file by decision with samtools...")
    bam_dir = os.path.join(args.output_directory, "split_bams")
    if not os.path.isdir(bam_dir):
        os.makedirs(bam_dir)
    
    # Function to write read IDs to a file
    def write_read_ids_to_file(read_ids, file_path):
        if read_ids:
            with open(file_path, 'w') as file:
                file.write("\n".join(read_ids) + "\n")

    # Write read IDs to files if the sets are not empty
    if accepted_readIDs:
        accepted_readIDs_file = os.path.join(bam_dir, f"{args.sample_name}_accepted_readIDs.txt")
        write_read_ids_to_file(accepted_readIDs, accepted_readIDs_file)
        accepted_bam = os.path.join(bam_dir, f'{args.sample_name}_accepted.bam')
        subprocess.run(f"samtools view -b -N {accepted_readIDs_file} -o {accepted_bam} {args.bam_file}", shell=True, check=True)
        subprocess.run(['samtools', 'index', accepted_bam])

    if rejected_readIDs:
        rejected_readIDs_file = os.path.join(bam_dir, f"{args.sample_name}_rejected_readIDs.txt")
        write_read_ids_to_file(rejected_readIDs, rejected_readIDs_file)
        rejected_bam = os.path.join(bam_dir, f'{args.sample_name}_rejected.bam')
        subprocess.run(f"samtools view -b -N {rejected_readIDs_file} -o {rejected_bam} {args.bam_file}", shell=True, check=True)
        subprocess.run(['samtools', 'index', rejected_bam])

    if nodecision_readIDs:
        nodecision_readIDs_file = os.path.join(bam_dir, f"{args.sample_name}_nodecision_readIDs.txt")
        write_read_ids_to_file(nodecision_readIDs, nodecision_readIDs_file)
        nodecision_bam = os.path.join(bam_dir, f'{args.sample_name}_nodecision.bam')
        subprocess.run(f"samtools view -b -N {nodecision_readIDs_file} -o {nodecision_bam} {args.bam_file}", shell=True, check=True)
        subprocess.run(['samtools', 'index', nodecision_bam])

def get_read_info(args, accepted_readIDs, rejected_readIDs, nodecision_readIDs, on_target_readIDs, off_target_readIDs):
    '''Extract read information (ID, decision, target_type, length, Phred_quality_score, alignment_score, mapping_quality_score) from the FASTQ file and BAM file'''
    
    # Initialize dictionary to store data for each read
    read_data = {}
    
    # Read in the FASTQ file line by line (to reduce memory usage)
    with open(args.fastq_file, 'r') as f:
        print("Extracting decision and target type, calculating read length, and obtaining Phred quality score for each read ID from the FASTQ file...")
        read_count = 0 # Initialize progress counter (to display progress every 1 million reads)
        for line in f:
            read_ID_line = line.strip()
            sequence_line = next(f).strip()
            placeholder = next(f).strip()
            Phred_quality_score_line = next(f).strip()
            # Extract read ID from the first line
            read_ID = read_ID_line.split()[0][1:]
            # Determine decision and target type for the current read
            decision = (
                'accepted' if read_ID in accepted_readIDs else
                'rejected' if read_ID in rejected_readIDs else
                'nodecision' if read_ID in nodecision_readIDs else
                'other'
            ) if args.exp_type == "AS" else 'not_applicable'
            target_type = (
                'on_target' if read_ID in on_target_readIDs else
                'off_target' if read_ID in off_target_readIDs else
                'unmapped'
            )
            # Calculate read length
            length = len(sequence_line)
            # Calculate the error rate for each base
            error_probabilities = [10 ** (-(ord(char) - 33) / 10.0) for char in Phred_quality_score_line]
            # Calculate the mean error probability
            mean_error_probability = sum(error_probabilities) / len(error_probabilities)
            # Convert mean error probability back to Phred score
            Phred_quality_score = -10 * math.log10(mean_error_probability)
            # Store the read data
            read_data[read_ID] = {}
            read_data[read_ID]['decision'] = decision
            read_data[read_ID]['target_type'] = target_type
            read_data[read_ID]['length'] = length
            read_data[read_ID]['Phred_quality_score'] = Phred_quality_score
            # Increment read count and display progress every 1 million reads
            read_count += 1
            if read_count % 1000000 == 0:
                print("Processed ", read_count, " reads from FASTQ...")
    
    # Extract alignment and mapping quality scores from the BAM file
    print("Extracting alignment and mapping quality scores from the BAM file...")
    read_count = 0  # Re-initialize progress counter (to display progress every 1 million reads)
    # Process the BAM file line by line (to reduce memory usage)
    with pysam.AlignmentFile(args.bam_file, 'rb') as bam:
        for read_alignment in bam.fetch():
            read_ID = read_alignment.query_name
            # Ensure read_ID is in read_data (i.e., it was present in the FASTQ file)
            if read_ID not in read_data:
                print(f"Warning: Read ID {read_ID} from BAM file not found in FASTQ file.")
                continue
            # If read is unmapped, set alignment and mapping quality scores to be empty and continue
            if read_alignment.is_unmapped:
                read_data[read_ID]['alignment_score_per_base'] = None
                read_data[read_ID]['mapping_quality_score'] = None
                continue
            # Extract the alignment and mapping quality scores for the primary alignment
            if not read_alignment.is_secondary and not read_alignment.is_supplementary:
                read_data[read_ID]['alignment_score_per_base'] = (read_alignment.get_tag('AS') / read_alignment.query_length) if read_alignment.has_tag('AS') else None
                read_data[read_ID]['mapping_quality_score'] = read_alignment.mapping_quality
            # Increment read count and display progress every 1 million reads
            read_count += 1
            if read_count % 1000000 == 0:
                print("Processed ", read_count, " reads from BAM...")

    return read_data

def calculate_n50(lengths):
    '''Calculate N50 from a list of read lengths'''
    sorted_lengths = sorted(lengths, reverse=True)
    total_bases = sum(sorted_lengths)  
    n50 = 0
    running_total = 0
    for length in sorted_lengths:
        running_total += length
        if running_total >= total_bases / 2:
            n50 = length
            break
    return n50

def calculate_stats(decision, tsv_file, chunk_size):
    '''Calculate statistics for a given decision category'''

    # Initialize dictionaries to store data
    data = {}
    data['length'] = {}
    data['Phred_quality_score'] = {}
    data['alignment_score_per_base'] = {}
    data['mapping_quality_score'] = {}

    # Filter the DataFrame based on the decision category
    if decision != 'all':
        read_df = read_df[read_df['decision'] == decision]

    # Filter the DataFrame based on the target type and append read lengths, Phred quality scores, alignment scores, and mapping quality scores to the respective lists
    for target_type in ['on_target', 'off_target', 'unmapped']:
        filtered_df = read_df[read_df['target_type'] == target_type]
        data['length'][target_type] = (filtered_df['length'])
        data['Phred_quality_score'][target_type] = (filtered_df['Phred_quality_score'])
        if target_type != 'unmapped':
            data['alignment_score_per_base'][target_type] = (filtered_df['alignment_score_per_base'])
            data['mapping_quality_score'][target_type] = (filtered_df['mapping_quality_score'])

    # Append metrics for all reads and mapped reads
    data['length']['all'] = data['length']['on_target'] + data['length']['off_target'] + data['length']['unmapped']
    data['Phred_quality_score']['all'] = data['Phred_quality_score']['on_target'] + data['Phred_quality_score']['off_target'] + data['Phred_quality_score']['unmapped']
    data['length']['mapped'] = data['length']['on_target'] + data['length']['off_target']
    data['Phred_quality_score']['mapped'] = data['Phred_quality_score']['on_target'] + data['Phred_quality_score']['off_target']
    data['alignment_score_per_base']['mapped'] = data['alignment_score_per_base']['on_target'] + data['alignment_score_per_base']['off_target']
    data['mapping_quality_score']['mapped'] = data['mapping_quality_score']['on_target'] + data['mapping_quality_score']['off_target']

    # Initialize dictionaries to store results
    stats = {}

    # Calculate statistics for each target type
    for target_type in ['all', 'unmapped', 'mapped', 'on_target', 'off_target']:
        # Extract reads and bases
        if target_type == 'all':
            stats['reads'] = len(data['length'][target_type])
            stats['bases'] = sum(data['length'][target_type])
        else:
            stats[f'{target_type}_reads'] = len(data['length'][target_type])
            stats[f'{target_type}_bases'] = sum(data['length'][target_type])
        # Iterate over metrics and extract other statistics
        for metric in ['length', 'Phred_quality_score', 'alignment_score_per_base', 'mapping_quality_score']:
            # Skip alignment and mapping quality scores for unmapped reads
            if target_type in ['all', 'unmapped'] and metric in ['alignment_score_per_base', 'mapping_quality_score']:
                continue
            # Set label for target type
            if target_type == 'all':
                label = ''
            else:
                label = f'{target_type}_'
            stats[f'min_{label}{metric}'] = np.min(data[metric][target_type])
            stats[f'max_{label}{metric}'] = np.max(data[metric][target_type])
            stats[f'avg_{label}{metric}'] = np.mean(data[metric][target_type])
            stats[f'std_{label}{metric}'] = np.std(data[metric][target_type])
            stats[f'q1_{label}{metric}'] = np.percentile(data[metric][target_type], 25)
            stats[f'median_{label}{metric}'] = np.median(data[metric][target_type])
            stats[f'q3_{label}{metric}'] = np.percentile(data[metric][target_type], 75)
            if metric == 'length':
                stats[f'n50_{label}{metric}'] = calculate_n50(data[metric][target_type])

    return stats
        
def calculate_all_stats(read_df, chunk_size, experiment_type):
    '''Calculate statistics for all reads and each decision category'''

    # Initialize dictionary to store statistics
    stats_dict = defaultdict(dict)
    # First calculate stats for each decision category if experiment_type is AS
    if experiment_type == "AS":
        for decision in read_df['decision'].unique():
            print(f"Calculating stats for {decision} reads...")
            stats_dict[decision] = calculate_stats(decision, read_df)
    # Then calculate stats for all reads
    print("Calculating stats for all reads...")
    stats_dict['total'] = calculate_stats('all', read_df)

    return stats_dict

def create_histogram(read_df, plot_output_directory, exp_type, plot_type, sample_name, min=None, max=None):
    '''Create histogram of read lengths or qualities by decision and target type'''

    # Define plot parameters based on plot_type
    if plot_type == 'length_reads':
        suptitle = 'Read Length Distribution'
        x_label = 'Read Length (bp)'
        y_label = 'Reads'
        file_suffix = 'length_histogram_reads'
        data_column = 'length'
    elif plot_type == 'length_bases':
        suptitle = 'Bases Sequenced by Read Length'
        x_label = 'Read Length (bp)'
        y_label = 'Bases'
        file_suffix = 'length_histogram_bases'
        data_column = 'length'
    elif plot_type == 'Phred_quality_score':
        suptitle = 'Phred Quality Score Distribution'
        x_label = 'Phred Quality Score'
        y_label = 'Reads'
        file_suffix = 'Phred_quality_histogram'
        data_column = 'Phred_quality_score'
    elif plot_type == 'alignment_score_per_base':
        suptitle = 'Alignment Score Distribution'
        x_label = 'Alignment Score per Base'
        y_label = 'Reads'
        file_suffix = 'alignment_score_histogram'
        data_column = 'alignment_score_per_base'
    elif plot_type == 'mapping_quality_score':
        suptitle = 'Mapping Quality Score Distribution'
        x_label = 'Mapping Quality Score'
        y_label = 'Reads'
        file_suffix = 'mapping_quality_histogram'
        data_column = 'mapping_quality_score'
    
    # Define the order of decisions and target types to be plotted
    decisions_desired_order = ['not_applicable', 'accepted', 'rejected', 'nodecision']
    decisions_unique = read_df['decision'].unique()
    decisions_order = [decision for decision in decisions_desired_order if decision in decisions_unique]
    target_types_desired_order = ['on_target', 'off_target', 'unmapped']
    target_types_unique = read_df['target_type'].unique()
    target_types_order = [target for target in target_types_desired_order if target in target_types_unique]

    # Define colors and labels for each decision and target type
    decision_colors = {'accepted': '#2e8b57', 'rejected': '#7a6097', 'nodecision': 'black', 'not_applicable': 'black'}
    target_type_colors = {'on_target': '#db4b62', 'off_target': '#007fb7', 'unmapped': '#999da0'}

    # Filter data and create bins
    min_value = int(np.floor(read_df[data_column].min())) if min is None else min
    max_value = int(np.ceil(read_df[data_column].max())) if max is None else max
    bin_size = 10 if max_value - min_value > 1000 else 1 if max_value - min_value > 10 else 0.01
    bins = np.arange(min_value, max_value + bin_size, bin_size)
    read_df = read_df[(read_df[data_column] >= min_value) & (read_df[data_column] <= max_value)]

    # Create figure background
    fig = plt.figure(figsize=(10, 6))

    # Create grid to separate data by decision category
    num_decisions = len(decisions_order)
    gs = gridspec.GridSpec(nrows=num_decisions, ncols=1, figure=fig, hspace=0.5)

    for i, decision in enumerate(decisions_order):
        decision_df = read_df[read_df['decision'] == decision]

        # Create subplot for this decision
        ax = fig.add_subplot(gs[i])

        # Initialize cumulative counts
        cumulative_count = np.zeros(len(bins) - 1)

        for target_type in target_types_order:
            if target_type in decision_df['target_type'].unique():
                target_data = decision_df[decision_df['target_type'] == target_type]
            else:
                continue

            # Calculate histogram
            count, _ = np.histogram(target_data[data_column], bins=bins, weights=None if plot_type != 'length_bases' else target_data['length'])

            # Plot
            ax.bar(bins[:-1], count, width=bin_size, bottom=cumulative_count, label=target_type, color=target_type_colors[target_type], align='edge')
            cumulative_count += count

        # Set plot labels and titles
        if exp_type == "AS":
            ax.set_title(f'Adaptive Sampling: {decision.capitalize()}', color=decision_colors[decision], fontsize=16)
        if i == num_decisions - 1:
            ax.set_xlabel(x_label, fontsize=14)
        
        ax.set_ylabel(y_label, fontsize=14)
        ax.set_yticks([])
        ax.set_yticklabels([])

    # Adjust layout
    fig.subplots_adjust(bottom=0.20, top=0.85)

    # Add suptitle
    fig.suptitle(suptitle, fontsize=16)

    # Create custom legend handles
    handles = [plt.Rectangle((0, 0), 1, 1, fc=target_type_colors['on_target'], edgecolor='none'),
               plt.Rectangle((0, 0), 1, 1, fc=target_type_colors['off_target'], edgecolor='none'),
               plt.Rectangle((0, 0), 1, 1, fc=target_type_colors['unmapped'], edgecolor='none')]
    labels = ['On-Target', 'Off-Target', 'Unmapped']
    fig.legend(handles, labels, loc='lower center', ncol=3, frameon=False, fontsize=14)

    # Save figure
    output_filename = f"{sample_name}_{file_suffix}"
    fig.savefig(os.path.join(plot_output_directory, f"{output_filename}.pdf"), bbox_inches='tight')
    fig.savefig(os.path.join(plot_output_directory, f"{output_filename}.png"), bbox_inches='tight', dpi=300)
    print(f"{output_filename} saved.")

def create_plot(df, plot_output_directory, sample_name, exp_type, plot_type, min_length=0, max_length=None):  
    '''Create violin or bar plots of read lengths and quality scores by decision and target type'''

    # Define plot parameters based on plot_type
    if 'length' in plot_type:
        suptitle='Read Length by Target Type'
        y_label = 'Read Length (bp)'
        file_suffix = 'length_violin_plot' if 'violin' in plot_type else 'length_barplot'
        data_column = 'length'
    elif 'Phred_quality_score' in plot_type:
        suptitle='Phred Quality Score by Target Type'
        y_label = 'Phred Quality Score'
        file_suffix = 'Phred_quality_violin_plot' if 'violin' in plot_type else 'Phred_quality_barplot'
        data_column = 'Phred_quality_score'
    elif 'alignment_score_per_base' in plot_type:
        suptitle='Alignment Score by Target Type'
        y_label = 'Alignment Score per Base'
        file_suffix = 'alignment_score_violin_plot' if 'violin' in plot_type else 'alignment_score_barplot'
        data_column = 'alignment_score_per_base'
    elif 'mapping_quality_score' in plot_type:
        suptitle='Mapping Quality Score by Target Type'
        y_label = 'Mapping Quality Score'
        file_suffix = 'mapping_quality_violin_plot' if 'violin' in plot_type else 'mapping_quality_barplot'
        data_column = 'mapping_quality_score'
    
    # Define the order of decisions and target types to be plotted
    decisions_desired_order = ['not_applicable', 'accepted', 'rejected', 'nodecision']
    decisions_unique = df['decision'].unique()
    decisions_order = [decision for decision in decisions_desired_order if decision in decisions_unique]
    if 'alignment_score_per_base' in plot_type or 'mapping_quality_score' in plot_type:
        target_types_desired_order = ['on_target', 'off_target']
    else:
        target_types_desired_order = ['on_target', 'off_target', 'unmapped']
    target_types_unique = df['target_type'].unique()
    target_types_order = [target for target in target_types_desired_order if target in target_types_unique]

    # Define colors and labels for each decision and target type
    decision_colors = {'accepted': '#2e8b57', 'rejected': '#7a6097', 'nodecision': 'black', 'not_applicable': 'black'}
    target_type_colors = {'on_target': '#db4b62', 'off_target': '#007fb7', 'unmapped': '#999da0'}

    # Define min and max values to be plotted
    if 'length' in plot_type:
        min_value = 0 if min_length is None else min_length
        max_value = df[data_column].max() if max_length is None else max_length
    else:
        min_value = min(df[data_column].min(), 0)
        max_value = df[data_column].max()

    # Create figure background
    fig = plt.figure(figsize=(10, 6))
    
    # Create grid to separate data by decision category
    num_decisions = len(decisions_order)
    gs = gridspec.GridSpec(nrows=1, ncols=num_decisions, figure=fig, hspace=0.5)

    # Iterate over decisions
    for i, decision in enumerate(decisions_order):
        decision_data = df[df['decision'] == decision]
        num_target_types = len(target_types_order)
        gs_decision = gridspec.GridSpecFromSubplotSpec(1, num_target_types, subplot_spec=gs[i], wspace=0.1)

        # Initialize a list to hold the center position of each subplot for this decision
        subplot_centers = []
        
        # Iterate over target types
        for j, target_type in enumerate(target_types_order):
            target_data = decision_data[decision_data['target_type'] == target_type]
            if target_data.empty:
                print(f"No data for {decision}, {target_type}")
                continue

            # Filter read length outliers for better plotting
            if 'length' in plot_type:
                Q1 = np.percentile(target_data[data_column], 25)
                Q3 = np.percentile(target_data[data_column], 75)
                IQR = Q3 - Q1
                lower_bound = Q1 - 1.5 * IQR
                upper_bound = Q3 + 1.5 * IQR
                outliers = target_data[(target_data[data_column] < lower_bound) | (target_data[data_column] > upper_bound)]
                target_data_filtered = target_data[(target_data[data_column] >= lower_bound) & (target_data[data_column] <= upper_bound)]
            else:
                target_data_filtered = target_data
            
            # Make plot
            ax = fig.add_subplot(gs_decision[j])
            if 'violin' in plot_type:
                # Plot violin
                parts = ax.violinplot(target_data_filtered[data_column], positions=[0], widths=0.8, showmeans=False, showextrema=False, showmedians=False, points=200, bw_method=0.15)
                parts['bodies'][0].set_facecolor(target_type_colors[target_type])
                parts['bodies'][0].set_edgecolor('black')
                parts['bodies'][0].set_alpha(0.7)
                # Overlay box plot
                ax.boxplot(target_data[data_column], positions=[0], showfliers=False, patch_artist=True, boxprops=dict(facecolor=target_type_colors[target_type], color='black'), whiskerprops=dict(color='black'), \
                            medianprops=dict(color='black'), capprops=dict(color='black'))
                # Plot outliers
                if plot_type == 'length_violin' and not outliers.empty:
                    x = np.random.normal(0, 0.04, size=len(outliers))
                    ax.scatter(x, outliers[data_column], color=target_type_colors[target_type], edgecolor=None, alpha=0.7, s=0.5)

            elif 'bar' in plot_type:
                # Calculate proportions
                bin_size = 10 if max_value - min_value > 100 else 1 if max_value - min_value > 10 else 0.01
                bin_edges = np.arange(min_value, max_value + bin_size, bin_size)
                hist_counts, _ = np.histogram(target_data_filtered[data_column], bins=bin_edges)
                total_count = len(target_data_filtered)
                proportions = hist_counts / total_count if total_count > 0 else np.zeros_like(hist_counts)
                # Plot histogram
                ax.barh(bin_edges[:-1], proportions, height=bin_size, color=target_type_colors[target_type], alpha=0.7, label=target_type, edgecolor=(target_type_colors[target_type] if plot_type == 'length_bar' else 'black'))


            # Set axis labels and limits
            ax.set_xticks([])
            ax.set_xticklabels([])
            if i == 0 and j == 0:
                ax.set_ylabel(y_label, fontsize=14)
            else:
                ax.set_yticks([])
                ax.set_yticklabels([])
            if 'length' in plot_type:
                ax.set_ylim(min_value, max_value)

            # Remove spines
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            
            # Get center position of the subplot
            x_position = ax.get_position().x0
            width = ax.get_position().width
            subplot_center = x_position + width / 2
            subplot_centers.append(subplot_center)
        
        # Add title for the decision
        if exp_type == "AS":
            decision_center = sum(subplot_centers) / len(subplot_centers)
            fig.text(decision_center, 0.17, f'Adaptive Sampling: {decision.capitalize()}', ha='center', va='top', color=decision_colors[decision], fontsize=14)
        
    # Adjust layout
    fig.subplots_adjust(bottom=0.20, top=0.9)
    
    # Add title
    fig.suptitle(suptitle, y=0.98, fontsize=16)

    # Create custom legend handles
    handles = [plt.Rectangle((0, 0), 1, 1, fc=target_type_colors['on_target'], edgecolor='none'),
                plt.Rectangle((0, 0), 1, 1, fc=target_type_colors['off_target'], edgecolor='none'),
                plt.Rectangle((0, 0), 1, 1, fc=target_type_colors['unmapped'], edgecolor='none')]
    labels = ['On-Target', 'Off-Target', 'Unmapped']
    fig.legend(handles, labels, loc='lower center', ncol=3, frameon=False, fontsize=14)

    # Save figures
    output_filename = f"{sample_name}_{file_suffix}"
    plt.savefig(os.path.join(plot_output_directory, f"{output_filename}.pdf"), bbox_inches='tight')
    plt.savefig(os.path.join(plot_output_directory, f"{output_filename}.png"), bbox_inches='tight', dpi=300)
    print(f"{output_filename} saved.")

def create_plots(read_df, plot_output_directory, args):
    ''' Create all plots'''

    print("Making plots...")
    create_histogram(read_df, plot_output_directory, 'length_reads', args.exp_type, args.sample_name, args.min_length, args.max_length)
    create_histogram(read_df, plot_output_directory, 'length_bases', args.exp_type, args.sample_name, args.min_length, args.max_length)
    create_histogram(read_df, plot_output_directory, 'Phred_quality_score', args.exp_type, args.sample_name)
    create_histogram(read_df, plot_output_directory, 'alignment_score_per_base', args.exp_type, args.sample_name)
    create_histogram(read_df, plot_output_directory, 'mapping_quality_score', args.exp_type, args.sample_name)
    create_plot(read_df, plot_output_directory, 'length_violin', args.exp_type, args.sample_name, args.min_length, args.max_length)
    create_plot(read_df, plot_output_directory, 'Phred_quality_score_violin', args.exp_type, args.sample_name)
    create_plot(read_df, plot_output_directory, 'alignment_score_per_base_violin', args.exp_type, args.sample_name)
    create_plot(read_df, plot_output_directory, 'mapping_quality_score_violin', args.exp_type, args.sample_name)
    create_plot(read_df, plot_output_directory, 'length_bar', args.exp_type, args.sample_name, args.min_length, args.max_length)
    create_plot(read_df, plot_output_directory, 'Phred_quality_score_bar', args.exp_type, args.sample_name)
    create_plot(read_df, plot_output_directory, 'alignment_score_per_base_bar', args.exp_type, args.sample_name)
    create_plot(read_df, plot_output_directory, 'mapping_quality_score_bar', args.exp_type, args.sample_name)

def main():
    '''Main function'''

    # Parse arguments
    args = parse_arguments()
    print()

    # Make output directory if it does not exist
    if args.output_plots and not os.path.isdir(plot_output_directory):
        plot_output_directory = os.path.join(args.output_directory, 'plots')
        os.makedirs(plot_output_directory)
        print(f"Created output directory.\n")
    elif not os.path.isdir(args.output_directory):
        os.makedirs(args.output_directory)
        print(f"Created output directory.\n")

    # Identify accepted, rejected, and nodecision read IDs (if applicable)
    if args.experiment_type == "AS":
        accepted_readIDs, rejected_readIDs, nodecision_readIDs = get_decision(args)
    else:
        accepted_readIDs, rejected_readIDs, nodecision_readIDs = set(), set(), set()
    # Determine on-target and off-target read IDs
    on_target_readIDs, off_target_readIDs = get_target_type(args.bam, args.bed)
    print(f"Read ID information collected.\n")    

    # Split BAM file by decision if experiment type is AS and --split_bam_by_decision flag is set
    if args.experiment_type == "AS" and args.split_bam_by_decision:
        split_bam(args, accepted_readIDs, rejected_readIDs, nodecision_readIDs)
        print(f"BAM files split by decision.\n")
    
    # Iterate over each read in the FASTQ file and determine its decision, target type, length, and quality score
    read_data = get_read_info(args, accepted_readIDs, rejected_readIDs, nodecision_readIDs, on_target_readIDs, off_target_readIDs)
    print(f"Read information, length, and quality scores collected.\n")
    
    # Create a dataframe from the collected read information and write to a TSV file
    print("Converting read information to DataFrame and writing to TSV file...")
    read_df = pd.DataFrame.from_dict(read_data, orient='index')
    read_df.reset_index(inplace=True)
    read_df.rename(columns={'index': 'Read_ID'}, inplace=True)
    read_info_file = os.path.join(args.output_directory, f"{args.sample_name}_read_information.tsv")
    read_df.to_csv(read_info_file, sep='\t', index=False)
    print(f"Read information written to {read_info_file}\n")

    # Collect statistics
    stats_dict = calculate_all_stats(read_df, args.chunk_size, args.experiment_type)
    stats_df = pd.DataFrame.from_dict(stats_dict, orient='index')
    print("Stats calculated successfully.")

    # Write the DataFrames to TSV files
    summary_stats_file = os.path.join(args.output_directory, f"{args.sample_name}_summary_statistics.tsv")
    stats_df.T.to_csv(summary_stats_file, sep='\t')
    print(f"Summary statistics written to {summary_stats_file}\n")

    # Create plots
    if args.output_plots:
        rcParams['pdf.fonttype'] = 42
        create_plots(read_df, plot_output_directory, args)

if __name__ == "__main__":
    main()



