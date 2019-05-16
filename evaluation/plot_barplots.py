#!/usr/bin/env python
import pandas as pd
#import mathe
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import argparse
import os
from operator import itemgetter



def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_df", action="store", dest="input_df", required=True, default="./data/all_filtered_massspec.gtf"
                                           , help= "path to the reference data.")
    parser.add_argument("-o", "--save_path", action="store", dest="save_path", required=True, default="./result_dfs/"
                                      , help= "path to save data to.")

    args = parser.parse_args()

    plot_dir = args.save_path

    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)
        
    plot_key_dir = plot_dir

    df_stat = pd.read_table(args.input_df,sep='\t')

    ##############################################################################################################
    ##############################################################################################################
    sns.set(style="whitegrid", font_scale=1)

    fig, ax = plt.subplots()

    sns.barplot(x="tool", y="FNR", data=df_stat, ax=ax)

    ax.set_xlabel('')
    #ax.set_ylabel('recovery')

    ax.set_ylim([0,1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    fig.show()
    fig.savefig(plot_key_dir + 'bar_FNR.pdf', format='pdf', dpi=300, bbox_inches='tight')

    ##############################################################################################################
    ##############################################################################################################
    sns.set(style="whitegrid", font_scale=1)

    fig, ax = plt.subplots()

    sns.barplot(x="tool", y="recall", data=df_stat, ax=ax)

    ax.set_xlabel('')
    #ax.set_ylabel('recovery')

    ax.set_ylim([0,1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    fig.show()
    fig.savefig(plot_key_dir + 'bar_recall.pdf', format='pdf', dpi=300, bbox_inches='tight')

    ##############################################################################################################
    ##############################################################################################################
    sns.set(style="whitegrid", font_scale=1)

    fig, ax = plt.subplots()

    sns.barplot(x="tool", y="presision", data=df_stat, ax=ax)

    ax.set_xlabel('')
    #ax.set_ylabel('recovery')

    ax.set_ylim([0,1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    fig.show()
    fig.savefig(plot_key_dir + 'bar_presision.pdf', format='pdf', dpi=300, bbox_inches='tight')

    ##############################################################################################################
    ##############################################################################################################
    sns.set(style="whitegrid", font_scale=1)

    fig, ax = plt.subplots()

    sns.barplot(x="tool", y="FDR", data=df_stat, ax=ax)

    ax.set_xlabel('')
    #ax.set_ylabel('recovery')

    ax.set_ylim([0,1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    fig.show()
    fig.savefig(plot_key_dir + 'bar_FDR.pdf', format='pdf', dpi=300, bbox_inches='tight')

    ##############################################################################################################
    ##############################################################################################################
    sns.set(style="whitegrid", font_scale=1)

    fig, ax = plt.subplots()

    sns.barplot(x="tool", y="F1", data=df_stat, ax=ax)

    ax.set_xlabel('')
    #ax.set_ylabel('recovery')

    ax.set_ylim([0,1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    fig.show()
    fig.savefig(plot_key_dir + 'bar_F1.pdf', format='pdf', dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    main()
