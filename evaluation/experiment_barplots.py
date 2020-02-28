#!/usr/bin/env python
import pandas as pd
#import mathe
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import argparse
import os
from operator import itemgetter
import pickle
import glob



def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-e", "--experiment_dict_path", action="store", dest="experiment_dict_path", required=True
                                           , help= "path to dict containing stats.cvs path ordered. Key=experiment Value=[(dataset_path, hue_value),(...),...].")
#    parser.add_argument("-i2", "--input2_df", action="store", dest="input2_df", required=True
#                                           , help= "path to the reference data.")
    parser.add_argument("-o", "--save_path", action="store", dest="save_path", required=True
                                      , help= "path to save data to.")

    args = parser.parse_args()
    experiment_dict_path = args.experiment_dict_path
    plot_dir = args.save_path

    with open(args.experiment_dict_path, 'rb') as handle:
        experiment_dict = pickle.load(handle)

    # generate empty dict
    #df_stat = pd.DataFrame(columns=["TP", 'TN', "FP", "FN", "recall", "specificity", "FNR", "precision", "FDR", "F1", "accuracy", "subopt_tp", "subopt_fp", "tool"])
    dummy = pd.DataFrame(columns=["TP", 'TN', "FP", "FN", "recall", "specificity", "FNR", "precision", "FDR", "F1", "accuracy", "subopt_tp", "subopt_fp", "no_genes", "tool"], data=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,"reparation"],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,"ribotish"],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,"deepribo"], [0,0,0,0,0,0,0,0,0,0,0,0,0,0,"irsom"]])

    for experiment in experiment_dict:
    	# print(experiment_dict[experiment])
        df_stat = pd.DataFrame(columns=["TP", 'TN', "FP", "FN", "recall", "specificity", "FNR", "precision", "FDR", "F1", "accuracy", "subopt_tp", "subopt_fp", "no_genes", "tool"], data=[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,"test"]])

        #df_stat = pd.DataFrame(columns=["TP", 'TN', "FP", "FN", "recall", "specificity", "FNR", "precision", "FDR", "F1", "accuracy", "subopt_tp", "subopt_fp", "tool"])
        no_entys = len(experiment_dict[experiment])
        #print(experiment_dict[experiment])
        #print('\n')
        #print(no_entys)
        #print('\n')
        counter = 0

        for path in experiment_dict[experiment]:
            path_stat = path[0]
            # print(path_stat)
            if os.path.isfile(path_stat):
                df_temp = pd.read_csv(path_stat,sep='\t')
                df_temp['hue_val'] = path[1]
                df_temp['dataset'] = path[2]
                df_temp['overlap'] = path[3]
                df_stat = df_stat.append(df_temp, ignore_index=True, sort=False)
            else:
                # df_temp = dummy
                print('no entrys for file %s'%(path_stat))

            counter += 1
            if counter == no_entys:
                df_stat = df_stat[df_stat.tool != 'test']
                # print(len(df_stat))
                # barplot
                #sns.set(style="whitegrid", font_scale=1)
                #fig, ax = plt.subplots()

                #sns.barplot(x="tool", y="accuracy", data=df_stat, hue="hue_val", ax=ax)

                #ax.set_xlabel('accuracy')
                #ax.set_title('experiment')
                #ax.set_ylabel('recovery')
                #ax.set_ylim([0,1])
                #ax.set_xticklabels(ax.get_xticklabels(), rotation=0)

                #fig.savefig(plot_dir + experiment + '_bar_accuracy.pdf', format='pdf', dpi=300, bbox_inches='tight')

                ##############################################lineplot
                # sns.set(style="whitegrid", font_scale=1)
                # fig, ax = plt.subplots()
                # #palette=[myblue,myblue,mygreen,mygreen],
                #
                # sns.pointplot(x="tool", y="accuracy", data=df_stat, hue="hue_val", ax=ax,
                #              linestyles=['-','--','-','--','-','--'],
                #              markers=['o','o','d','d','*','*'])
                #
                # ax.set_xlabel('accuracy')
                # #ax.legend(frameon=True, fancybox=True,  framealpha=1, shadow=True, title='')
                #
                # ax.set_title('experiment')
                # #ax.set_ylabel('recovery')
                # ax.set_ylim([0,1])
                # ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
                #
                # fig.savefig(plot_dir + experiment + '_line_accuracy.pdf', format='pdf', dpi=300, bbox_inches='tight')
                #
                sns.set(style="whitegrid", font_scale=1)
                fig, ax = plt.subplots()
                #palette=[myblue,myblue,mygreen,mygreen],

                sns.pointplot(x="hue_val", y="accuracy", data=df_stat, hue="tool", ax=ax, join=False, markers=['o','d','*','+'])

                ax.set_title('')
                #ax.legend(frameon=True, fancybox=True,  framealpha=1, shadow=True, title='')

                #ax.set_title('experiment')
                ax.set_ylim([-0.02,1.02])
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

                fig.savefig(plot_dir + experiment + '_line_tools_accuracy.pdf', format='pdf', dpi=300, bbox_inches='tight')

                #################################################################################
                # filterde data:
                # df_filtered_overlap_90 = df_stat[(df_stat['overlap'] == '0.9')]
                #
                # sns.set(style="whitegrid", font_scale=1)
                # fig, ax = plt.subplots()
                #
                # sns.barplot(x="dataset", y="accuracy", data=df_filtered_overlap_90, hue="tool", ax=ax)
                #
                # ax.set_xlabel('tools')
                # ax.set_ylabel('accuracy')
                # #ax.set_title(experiment + ' overlap 90%')
                # ax.set_ylim([0,1])
                # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
                #
                # fig.savefig(plot_dir + experiment + '_bar_filterd_accuracy.pdf', format='pdf', dpi=300, bbox_inches='tight')

                # filterde data:

                sns.set(style="whitegrid", font_scale=1)
                fig, ax = plt.subplots()

                sns.pointplot(x="hue_val", y="specificity", data=df_stat, hue="tool", ax=ax, join=False, markers=['o','d','*','+'])


                ax.set_title('')
                ax.set_ylim([-0.02,1.02])
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

                fig.savefig(plot_dir + experiment + '_line_tools_specificity.pdf', format='pdf', dpi=300, bbox_inches='tight')


                ##################################################################################

                sns.set(style="whitegrid", font_scale=1)
                fig, ax = plt.subplots()
                #palette=[myblue,myblue,mygreen,mygreen],

                sns.pointplot(x="hue_val", y="recall", data=df_stat, hue="tool", ax=ax, join=False, markers=['o','d','*','+'])

                ax.set_title('')
                #ax.legend(frameon=True, fancybox=True,  framealpha=1, shadow=True, title='')

                ax.set_ylim([-0.02,1.02])
                ax.set_xticklabels(ax.get_xticklabels(), rotation=45)

                fig.savefig(plot_dir + experiment + '_line_tools_recall.pdf', format='pdf', dpi=300, bbox_inches='tight', markers=['o','d','*','+'])

                # sns.set(style="whitegrid", font_scale=1)
                # fig, ax = plt.subplots()
                #
                # sns.barplot(x="dataset", y="accuracy", data=df_filtered_overlap_90, hue="tool", ax=ax)
                #
                # ax.set_xlabel('tools')
                # ax.set_ylabel('accuracy')
                # ax.set_title(experiment + ' overlap 90%')
                # #ax.set_ylabel('recovery')
                # ax.set_ylim([0,1])
                # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
                #
                # fig.savefig(plot_dir + experiment + '_bar_filterd_accuracy.pdf', format='pdf', dpi=300, bbox_inches='tight')
                #
                # sns.set(style="whitegrid", font_scale=1)
                # fig, ax = plt.subplots()
                # #palette=[myblue,myblue,mygreen,mygreen],
                #
                # sns.pointplot(x="hue_val", y="F1", data=df_stat, hue="tool", ax=ax)
                #
                # ax.set_xlabel('F1')
                # #ax.legend(frameon=True, fancybox=True,  framealpha=1, shadow=True, title='')
                #
                # ax.set_title(experiment)
                # #ax.set_ylabel('recovery')
                # ax.set_ylim([0,1])
                # ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
                #
                # fig.savefig(plot_dir + experiment + '_line_tools_F1.pdf', format='pdf', dpi=300, bbox_inches='tight')








if __name__ == '__main__':
    main()
