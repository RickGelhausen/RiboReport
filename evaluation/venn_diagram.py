import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from itertools import combinations
from simple_venn import venn2, venn3, venn4
import argparse
import os, sys
import numpy as np

def get_sets(a,b,c,d):
    aANDbANDcANDd =  len(list(set(a) & set(b) & set(c) & set(d)))
    aANDbANDc = len(list(set(a) & set(b) & set(c)))
    aANDbANDd = len(list(set(a) & set(b) & set(d)))
    aANDcANDd = len(list(set(a) & set(c) & set(d)))
    bANDcANDd = len(list(set(b) & set(c) & set(d)))

    aANDb = len(list(set(a) & set(b)))
    aANDc = len(list(set(a) & set(c)))
    aANDd = len(list(set(a) & set(d)))
    bANDc = len(list(set(b) & set(c)))
    bANDd = len(list(set(b) & set(d)))
    cANDd = len(list(set(c) & set(d)))
    
    
    a_remainder = len(list(set(a) - set(b)- set(c)- set(d)))
    b_remainder = len(list(set(b) - set(a)- set(c)- set(d)))
    c_remainder = len(list(set(c) - set(a)- set(b)- set(d)))
    d_remainder = len(list(set(d) - set(a)- set(b)- set(c)))
    subsets = [a_remainder, b_remainder, c_remainder, d_remainder, aANDb,aANDc,aANDd, bANDc, bANDd, cANDd, aANDbANDc, aANDbANDd, aANDcANDd, bANDcANDd, aANDbANDcANDd]
    return subsets




def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input1_df", action="store", dest="input_df_path", required=True, 
                                            help= "path to the reference data.")                                      
    parser.add_argument("-o", "--save_path", action="store", dest="save_path", required=True, 
                                       help= "path to save data to.")
    parser.add_argument("-n", "--name_folder", action="store", dest="name_folder", required=True, 
                                       help= "name of the folder containing all plots e.g. ncRNA.")
    parser.add_argument("-c", "--coverage_percent", action="store", dest="coverage_percent", required=True, 
                                       help= "filterd on percent overlab e.g. 05 or 001.")                                   

    args = parser.parse_args()

    plot_dir = args.save_path
    input_df = args.input_df_path
    name_folder = args.name_folder
    coverage_percent = args.coverage_percent
    
    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)
    
    # generate save path
    save_venn_diag = plot_dir + '/venn_diagram'+ coverage_percent + '_' + name_folder + '.pdf'
    
    
    df_stat = pd.read_csv(input_df,sep='\t')
    # NA is defined in the stat.py script. If this is changed that comand also dose not work    
    df_stat.replace('NA',np.NaN)
    
    # order of tools:
    tools = ('DeepRibo', 'Ribo-TISH', 'Reparation', 'IRSOM')
    # removes nan and convert into lists
    tool1_list = list(df_stat['deepribo'].dropna())
    tool2_list = list(df_stat['ribotish'].dropna())
    tool3_list = list(df_stat['reparation'].dropna())
    tool4_list = list(df_stat['irsom'].dropna())
    
    #compute the numbers for each subset
    subsets = get_sets(tool1_list,tool2_list,tool3_list,tool4_list)
    
    
    ##################################plot############################
    fig, ax = plt.subplots() # 1, 1, figsize=(24, 8)
    
    venn4(subsets, set_labels=tools, ax=ax)
    #ax.set_title('venn4', fontsize=24)
    #fig.suptitle('simple_venn Demo', fontsize=30)
    plt.savefig(save_venn_diag, format='pdf', dpi=300, bbox_inches='tight')



if __name__ == '__main__':
    main()
