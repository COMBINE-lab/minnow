import sys
import os
import numpy as np
import pandas as pd
import gzip
import argparse
from collections import defaultdict


def read_quants_csv_alevin(truth_dir):
    truth = pd.read_csv(os.path.join(truth_dir, 'quants_mat.csv'), header = None)
    # We have to take care of trailing comma
    truth.drop(truth.columns[[-1,]], axis=1, inplace=True)
    
    cell_names = pd.read_csv(os.path.join(truth_dir, 'quants_mat_rows.txt'), header = None)[0].values
    gene_names = pd.read_csv(os.path.join(truth_dir, 'quants_mat_cols.txt'), header = None)[0].values
    
    assert(len(truth.columns) == len(gene_names))
    assert(len(truth.index) == len(cell_names))
    
    truth.index = cell_names
    truth.columns = gene_names
    print('Done reading truth with shape ',truth.shape)
    return (truth)


def read_sim_count(filename, minnow_dir = ""):
    cell_count = {}
    cell_dict = {}
    if minnow_dir != "":
        cell_names_10X = pd.read_csv(
            os.path.join(minnow_dir, 'alevin','quants_mat_rows.txt'), header = None
        )[0].values
        
        original_cell_names = pd.read_csv(
            os.path.join(minnow_dir, 'alevin','true_cell_names.txt'), header = None
        )[0].values
        print(original_cell_names[:2])
        cell_dict = dict(zip(cell_names_10X, original_cell_names))
        
    with gzip.open(filename) as fp:
        num_of_cells = int(fp.readline().strip())
        for cell_id in range(num_of_cells):
            gene_count_dict = {}
            cell_name = fp.readline().strip()
            num_of_genes = int(fp.readline().strip())
            print(cell_id, end='\r')
            #print(cell_name, num_of_genes)
            for gene_id in range(num_of_genes):
                #print(gene_id)
                gene_name, gene_count = fp.readline().strip().decode('ascii').split('\t')
                gene_count_dict[gene_name] = int(gene_count)
            cell_count[cell_name] = gene_count_dict
        count_df = pd.DataFrame.from_dict(cell_count, orient='index')
        if len(cell_dict) != 0:
            #count_df.index = [cell_id.decode('ascii') for cell_id in count_df.index]
            count_df.index = [cell_dict[cell_id.decode('ascii')] for cell_id in count_df.index]
        else:
            count_df.index = [cell_id.decode('ascii') for cell_id in count_df.index]
        count_df.fillna(0, inplace = True)
        print("\nDone reading ",filename,' with shape ',count_df.shape)
    return count_df, cell_dict


def cell_level_gene_histogram(truth, count_df):
    cell_level_dict = {}
    bad_cells = 0
    good_cells = 0
    ind = 0
    for cell_id in truth.index:
        print(ind, end='\r')
        ind += 1
        df_joined = pd.DataFrame(truth.loc[cell_id]).join(
            pd.DataFrame(count_df.loc[cell_id]),
            lsuffix = '_truth', rsuffix = '_fastq',
            how = 'outer'
        ).fillna(0)
        df_joined.columns = ['truth', 'fastq']
        d = len(df_joined.loc[
            ((df_joined.truth != 0) & (df_joined.fastq == 0))|
            ((df_joined.truth == 0) & (df_joined.fastq != 0))
        ])
        if d != 0:
            bad_cells += 1
        else:
            good_cells += 1
        cell_level_dict[cell_id] = df_joined
    if(bad_cells == 0):
        print('Comparison test --- ',u'\u2713')
    else:
        print('Comparison test --- ',u'\u2717')


def main():
    truth_dir = sys.argv[1]
    exp_dir = sys.argv[2]
    truth = read_quants_csv_alevin(truth_dir)
    count_df, cell_dict = read_sim_count(exp_dir)
    cell_level_gene_histogram(truth, count_df)


if __name__ == '__main__':
    main()