import pandas as pd 
import glob
import numpy as np
import argparse 

parser = argparse.ArgumentParser(description='Merge ztfquery metadata tables')
parser.add_argument('--meta-dir', type=str, required=True,
                    help='Folder containing metadata tables in csv format')
parser.add_argument('--outname', type=str, required=True,
                    help='Filename where merged table will be written')
args = parser.parse_args()

all_tables = np.sort(glob.glob(f'{args.meta_dir}/*.csv'))
print(f'Found {all_tables.size} table in {args.meta_dir}')
tabs = [pd.read_csv(t, index_col=0) for t in all_tables]

big_tab = pd.concat(tabs, ignore_index=True)
print(f'Merged table contains {len(big_tab)} images')
big_tab.to_csv(f'{args.outname}')
print(f'Merged table written to {args.outname}')
