import numpy as np
import csv, sys, os
import pandas as pd
import argparse
#from constants import *
# 
def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mutation_file_groups', type=str, required=True,
                        action='append', nargs='*')
    parser.add_argument('-ct', '--cancer_types', type=str, required=True, nargs='*')
    parser.add_argument('-pw', '--patient_whitelist', type=str, required=False)
    parser.add_argument('-hf', '--hypermutators_file', type=str, required=False, default=None)
    parser.add_argument('-ivc', '--ignored_variant_classes', type=str, required=False, nargs='*',
                        default=["Silent", "Intron", "3'UTR", "5'UTR", "IGR", "lincRNA", "RNA"])
    parser.add_argument('-ivt', '--ignored_variant_types', type=str, required=False, nargs='*',
                        default=['Germline'])
    parser.add_argument('-ivs', '--ignored_validation_statuses', type=str, required=False, nargs='*',
                        default=['Wildtype', 'Invalid'])
    parser.add_argument('-o', '--output_file', type=str, required=True)
    parser.add_argument('-v', '--verbose', type=int, default=1, required=False, choices=range(5))
    return parser

# Parse arguments
def get_parser():
    description = 'Script for processing a list of mutations across tumor samples in MAF format to a mutation matrix and to a format for '
    parser = argparse.ArgumentParser(description=description)
    # Files
    parser.add_argument('-maf', '--maf', required=True, help='MAF file name')
    parser.add_argument('-sf', '--sample_file', required=False, default=None, help='File with a list of samples to be included (optional).')
    parser.add_argument('-gf', '--gene_file', required=False, default=None, help='File with a list of genes to be included (optional).')
    # mutations
    parser.add_argument('-vce', '--variant_classes_to_exclude', type=str, required=False, nargs='*',
                        default=["Silent", "Intron", "3'UTR", "5'UTR", "5'Flank", "IGR", "lincRNA", "RNA"])
    # column names
    parser.add_argument('-G', '--GENE', required=False, default='Hugo_Symbol', help='Column name for gene')
    parser.add_argument('-MF', '--MUTATION_FEATURE', required=False, default='Hugo_Symbol', help='Column name for features')

    parser.add_argument('-S', '--SAMPLE', required=False, default='Tumor_Sample_Barcode', help='Column name for samples')
    parser.add_argument('-VC', '--VARIANT_CLASSIFICATION', required=False, 
                        default='Variant_Classification', help='Column name for variant class')
    # general
    parser.add_argument('-bo', '--binary_output_file', required=False, help='Name of output file (binary matrix)')
    parser.add_argument('-io', '--integer_output_file', required=False, help='Name of output file (integer matrix)')

    parser.add_argument('-v', '--verbosity', default=0, type=int, help='Flag verbose output')

    return parser


def run(args):
    # column names
    MAF_df = pd.read_csv(args.maf, sep='\t', index_col=False)
    if args.sample_file:
        if args.verbose > 0: print('* Loading samples to include...')
        with open(args.sample_file, 'r') as IN:
            sample_list = [l.rstrip().split()[0] for l in IN if not l.startswith("#")]
        print("Original number of cell lines: {}".format(len(MAF_df[args.SAMPLE].unique())))
        MAF_df = MAF_df[MAF_df[args.SAMPLE].isin(sample_list)]
        print("New number of cell lines: {}".format(len(MAF_df[args.SAMPLE].unique())))
    if args.gene_file:
        if args.verbose > 0: print('* Loading genes to include...')
        with open(args.gene_file, 'r') as IN:
            gene_list = [l.rstrip().split()[0] for l in IN if not l.startswith("#")]
        print("Original number of genes: {}".format(len(MAF_df[args.GENE].unique())))
        MAF_df = MAF_df[MAF_df[args.GENE].isin(gene_list)]
        print("New number of genes: {}".format(len(MAF_df[args.GENE].unique())))

    # Check what variant classes are in the data
    variant_classes_to_exclude = set(map(str.lower, args.variant_classes_to_exclude))
    original_variant_classes = MAF_df[args.VARIANT_CLASSIFICATION].unique()
    # remove 'conserving' mutations
    MAF_df = MAF_df[~MAF_df[args.VARIANT_CLASSIFICATION].isin(args.variant_classes_to_exclude)]
    # Check what variant classes have been removed
    new_variant_classes = MAF_df[args.VARIANT_CLASSIFICATION].unique()
    for variant in sorted(original_variant_classes):
        if variant not in new_variant_classes:
            print("CLASS  {:<16} has been removed".format(variant))

    # construct sample-feature matrices
    MAF_df['count'] = 1
    MAF_df_agg = MAF_df.groupby([args.MUTATION_FEATURE, args.SAMPLE]).agg(int_count=('count', 'sum'), bin_count=('count', 'min')).reset_index()
    if args.binary_output_file:
        bin_mat_df = MAF_df_agg.pivot_table(columns=args.MUTATION_FEATURE, index=args.SAMPLE, values="bin_count").reset_index()
        bin_mat_df = bin_mat_df.fillna(0)
        bin_mat_df.to_csv(args.binary_output_file, index=False, sep='\t')
    if args.integer_output_file:
        int_mat_df = MAF_df_agg.pivot_table(columns=args.MUTATION_FEATURE, index=args.SAMPLE, values="int_count").reset_index()
        int_mat_df = int_mat_df.fillna(0)
        int_mat_df.to_csv(args.integer_output_file, index=False, sep='\t')
    return

if __name__ == "__main__":
    run(get_parser().parse_args(sys.argv[1:]))
