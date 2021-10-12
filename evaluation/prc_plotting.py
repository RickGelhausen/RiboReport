import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
# import seaborn as sns
import argparse
import os
import sys
from operator import itemgetter
import pickle
import glob
import numpy as np
from sklearn import metrics
from sklearn.metrics import precision_recall_curve


def compute_prc(labels, scores):
    precision, recall, thresholds = metrics.precision_recall_curve(
    labels, scores)
    auc_prc =metrics.auc(recall, precision)
    #auc_prc = 0
    return precision, recall, thresholds, auc_prc

def get_list(saved_dir, tool):
    dir_score_list = saved_dir + '/' + tool + '_score_list'
    if not os.path.isfile(dir_score_list):
        sys.exit('not existing file (PRC): %s' % (dir_score_list))

    with open(dir_score_list, 'rb') as handle:
        score_overlap_label_list = pickle.load(handle)

    return score_overlap_label_list

def get_ranks(score_overlap_label_list, tool):
    #first_indexes = len(score_list)
    rank_list = list(range(len(score_overlap_label_list)))
    #print('score list:\n',len(score_overlap_label_list))
    #print('rank list:\n',rank_list)
    list.reverse(rank_list)

    #diff = len(score_overlap_label_list)-(len(set(score_overlap_label_list)))
    sorted_list = sorted(score_overlap_label_list, key=itemgetter(0,1), reverse=True)


    last_score = 0
    last_overlap = 0
    last_rank = 0
    label_list = []
    count_pos = 0
    count_neg = 0
    counter_deepribo_neg_scores = 0

    for pos, trippel in enumerate(sorted_list):
        label_list.append(int(trippel[2]))
        #count number of labels
        if trippel[2] == '1':
            count_pos += 1
            #if tool == 'deepribo' and trippel[0]<0:
        elif trippel[2] == '0':
            count_neg += 1

        # check if the current score the same as the last score
        if trippel[0] == last_score:
            #check if overlaps are the same
            if trippel[1] == last_overlap:
                rank_list[pos] = last_rank
            else:
                last_rank = rank_list[pos]
                last_score = trippel[0]
                last_overlap = trippel[1]
        else:
            last_rank = rank_list[pos]
            last_score = trippel[0]
            last_overlap = trippel[1]

    #print ('#####Rank List 2:#############')
    #print ('Number of positive instances: %i'%(count_pos))
    #print ('Number of negative instances: %i'%(count_neg))
    #print ('#####END#############')
    base = count_pos/(count_pos + count_neg)

    return rank_list, label_list, base



def comput_roc_for_tool(saved_dir, tool):

    score_overlap_label_list = get_list(saved_dir, tool)
    #print('score list of %s'%(tool))
    #print(score_overlap_label_list)
    sort_list = sorted(score_overlap_label_list, key=itemgetter(0), reverse=True)
    #print(sort_list)

    score_list, label_list, base = get_ranks(score_overlap_label_list, tool)
    #print('score ranked list of %s'%(tool))
    #print(score_list)
    #print('lable ranked list of %s'%(tool))
    #print(label_list)

    precision, recall, thresholds, auc_prc = compute_prc(label_list, score_list)

    return precision, recall, base, auc_prc

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_path", action="store", dest="experiment_dict_path", required=True
                                           , help= "path list of trippels (score, overlap, label) for each tool")
    parser.add_argument("-s", "--species_lable", action="store", dest="species_lable", required=True
                                           , help= "name of the datasets backterial")
    parser.add_argument("-o", "--save_path", action="store", dest="save_path", required=True
                                      , help= "path to save data to.")
    parser.add_argument("-p", "--percent_overlap", action="store", dest="percent_overlap", required=True
                                      , help= "how much of the prediction should overlap with the Gene and vice verser")
    parser.add_argument("-e", "--experiment", action="store", dest="experiment", required=True
                                           , help= "name of subset")

    args = parser.parse_args()
    experiment_dict_path = args.experiment_dict_path
    species = args.species_lable
    plot_dir = args.save_path
    coverage_percent = args.percent_overlap
    experiment = args.experiment


    overlap = coverage_percent.replace(".", "")
    # 'deepribo', 'ribotish', 'reparation', 'irsom', 'spectre'

    #print('tool: deepribo\n')
    precision_deepribo, recall_deepribo, base, auc_prc_deepribo = comput_roc_for_tool(experiment_dict_path, 'deepribo')
    #print('tool: ribotish\n' )
    precision_ribotish, recall_ribotish, base, auc_prc_ribotish = comput_roc_for_tool(experiment_dict_path, 'ribotish',)
    #print('tool: reparation\n' )
    precision_reparation, recall_reparation, base, auc_prc_reparation = comput_roc_for_tool(experiment_dict_path, 'reparation')
    #print('tool: irsom\n' )
    precision_irsom, recall_irsom, base, auc_prc_irsom = comput_roc_for_tool(experiment_dict_path, 'irsom')
    #print('tool: spectre\n' )
    precision_spectre, recall_spectre, base, auc_prc_spectre = comput_roc_for_tool(experiment_dict_path, 'spectre')
    #print('tool: price\n' )
    precision_price, recall_price, base, auc_prc_price = comput_roc_for_tool(experiment_dict_path, 'price')
    #print('tool: ribotricer\n' )
    precision_ribotricer, recall_ribotricer, base, auc_prc_ribotricer = comput_roc_for_tool(experiment_dict_path, 'ribotricer')
    #print('tool: smorfer\n' )
    precision_smorfer, recall_smorfer, base, auc_prc_smorfer = comput_roc_for_tool(experiment_dict_path, 'smorfer')
    #print('deepribo AUC: %f' % (auc_prc_deepribo))
    #print('ribotish AUC: %f' %  (auc_prc_ribotish))
    #print('reparation AUC: %f' %(auc_prc_reparation))
    #print('irsom AUC: %f' %  (auc_prc_irsom))
    #print('spectre AUC: %f' %  (auc_prc_spectre))


    if species == 'EC':
        title = 'E. coli'
    elif species == 'LM':
        title = 'L. monocytogenes'
    elif species == 'PA':
        title = 'P. aeruginosa'
    elif species == 'ST':
        title = 'S. Typhimurium'
    elif species == 'HV':
        title = 'Haloferax volcanii'
    else:
        print('Error: unknown species label')
        title = 'unknown species label'


    # generat the legend information
    label_deepribo = 'DeepRibo AUC: %.2f' % (auc_prc_deepribo)
    label_ribotish = ('Ribo-TISH AUC: %.2f' %  (auc_prc_ribotish))
    label_reparation = ('Reparation AUC: %.2f' %(auc_prc_reparation))
    label_irsom = ('IRSOM AUC: %.2f' %  (auc_prc_irsom))
    label_spectre = ('SPECtre AUC: %.2f' %  (auc_prc_spectre))
    label_price = ('price AUC: %.2f' %  (auc_prc_price))
    label_ribotricer = ('ribotricer AUC: %.2f' %  (auc_prc_ribotricer))
    label_smorfer = ('smorfer AUC: %.2f' %  (auc_prc_smorfer))

    # choosing colure

    # plot the PRC into one plot
    plt.plot(recall_deepribo, precision_deepribo, color='#009E73', label=label_deepribo) # green
    plt.plot(recall_reparation, precision_reparation, color='#0072B2', label=label_reparation) #blue
    plt.plot(recall_ribotish, precision_ribotish, color='#F0E442', label=label_ribotish) #yellow
    plt.plot(recall_irsom, precision_irsom, color='#D55E00', label=label_irsom) #red
    plt.plot(recall_spectre, precision_spectre, color='#E69F00', label=label_spectre) #orange
    plt.plot(recall_price, precision_price, color='#000000', label=label_price) #black
    plt.plot(recall_ribotricer, precision_ribotricer, color='#56B4E9', label=label_ribotricer) #sky blue
    plt.plot(recall_smorfer, precision_smorfer, color='#CC79A7', label=label_smorfer) #reddish purple

    # plot the baseline
    plt.plot([0, 1], [base, base], color='#BBBBBB', linestyle='--') #grey #808080
    # set axix labels and title
    plt.xlabel('Recall', fontsize=14)
    plt.ylabel('Precision', fontsize=14)
    plt.xticks(fontsize=12, rotation=30)
    plt.yticks(fontsize=12)
    #plt.title('Precision Recall Characteristic (prc) Curve ' + species)
    plt.title(title)

    #fontsize=20
    plt.legend(prop={'size': 12})
    plt.ylim(-0.02, 1.02)
    #plt.show()

    # save plot
    save_prc_diag = plot_dir + species + '_prc_' + overlap + '_' + experiment + '_.pdf'
    plt.savefig(save_prc_diag, format='pdf', dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    main()
