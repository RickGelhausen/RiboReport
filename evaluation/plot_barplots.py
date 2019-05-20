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
    parser.add_argument("-i1", "--input1_df", action="store", dest="input1_df", required=True
                                           , help= "path to the reference data.")
    parser.add_argument("-i2", "--input2_df", action="store", dest="input2_df", required=True
                                           , help= "path to the reference data.")
    parser.add_argument("-o", "--save_path", action="store", dest="save_path", required=True
                                      , help= "path to save data to.")

    args = parser.parse_args()

    plot_dir = args.save_path

    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)

    plot_key_dir = plot_dir


    dummy = pd.DataFrame(columns=["TP", "FP", "FN", "recall", "FNR", "precision", "FDR", "F1", "tool"], data=[[0,0,0,0,0,0,0,0,"reparation"],[0,0,0,0,0,0,0,0,"ribotish"],[0,0,0,0,0,0,0,0,"deepribo"], [0,0,0,0,0,0,0,0,"irsom"]])

    if os.path.isfile(args.input1_df):
        df_stat1 = pd.read_csv(args.input1_df,sep='\t')
    else:
        df_stat1 = dummy

    if os.path.isfile(args.input2_df):
        df_stat2 = pd.read_csv(args.input2_df,sep='\t')
    else:
        df_stat2 = dummy


    df_stat1['overlap'] = 0.01
    df_stat2['overlap'] = 0.5

    df_stat = pd.concat([df_stat1, df_stat2])


    ##############################################################################################################
    ##############################################################################################################
    sns.set(style="whitegrid", font_scale=1)

    fig, ax = plt.subplots()

    sns.barplot(x="tool", y="FNR", data=df_stat, hue="overlap", ax=ax)

    ax.set_xlabel('')
    #ax.set_ylabel('recovery')

    ax.set_ylim([0,1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    fig.savefig(plot_key_dir + 'bar_FNR.pdf', format='pdf', dpi=300, bbox_inches='tight')

    ##############################################################################################################
    ##############################################################################################################
    sns.set(style="whitegrid", font_scale=1)

    fig, ax = plt.subplots()

    sns.barplot(x="tool", y="recall", data=df_stat, hue="overlap", ax=ax)

    ax.set_xlabel('')
    #ax.set_ylabel('recovery')

    ax.set_ylim([0,1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    fig.savefig(plot_key_dir + 'bar_recall.pdf', format='pdf', dpi=300, bbox_inches='tight')

    ##############################################################################################################
    ##############################################################################################################
    sns.set(style="whitegrid", font_scale=1)

    fig, ax = plt.subplots()

    sns.barplot(x="tool", y="precision", data=df_stat, hue="overlap", ax=ax)

    ax.set_xlabel('')
    #ax.set_ylabel('recovery')

    ax.set_ylim([0,1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    fig.savefig(plot_key_dir + 'bar_precision.pdf', format='pdf', dpi=300, bbox_inches='tight')

    ##############################################################################################################
    ##############################################################################################################
    sns.set(style="whitegrid", font_scale=1)

    fig, ax = plt.subplots()

    sns.barplot(x="tool", y="FDR", data=df_stat, hue="overlap", ax=ax)

    ax.set_xlabel('')
    #ax.set_ylabel('recovery')

    ax.set_ylim([0,1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    fig.savefig(plot_key_dir + 'bar_FDR.pdf', format='pdf', dpi=300, bbox_inches='tight')

    ##############################################################################################################
    ##############################################################################################################
    sns.set(style="whitegrid", font_scale=1)

    fig, ax = plt.subplots()

    sns.barplot(x="tool", y="F1", hue="overlap", data=df_stat, ax=ax)

    ax.set_xlabel('')
    #ax.set_ylabel('recovery')

    ax.set_ylim([0,1])
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

    fig.savefig(plot_key_dir + 'bar_F1.pdf', format='pdf', dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    main()
