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
    #dir_lable_list = saved_dir + '/' + tool + '_label_list'
    if not os.path.isfile(dir_score_list):
        sys.exit('not existing file (PRC): %s' % (dir_score_list))

    with open(dir_score_list, 'rb') as handle:
        score_overlap_label_list = pickle.load(handle)
    #with open(dir_lable_list, 'rb') as handle:
        #label_list = pickle.load(handle)
    #print(len(label_list))
    #print(len(score_overlap_label_list))

    return score_overlap_label_list

def get_ranks(score_overlap_label_list, tool):
    #first_indexes = len(score_list)
    rank_list = list(range(len(score_overlap_label_list)))
    list.reverse(rank_list)

    #diff = len(score_overlap_label_list)-(len(set(score_overlap_label_list)))

    sorted_list = sorted(score_overlap_label_list, key=itemgetter(0,1), reverse=True)

    #print ('#####Ranked List 1:#############')
    #print (count)
    #print ('#####END#############')
    #print ('#####sorted List#############')
    #print (sorted_list)
    #print ('#####END#############')

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

    score_list, label_list, base = get_ranks(score_overlap_label_list, tool)

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
    # 'deepribo', 'ribotish', 'reparation', 'irsom'


    precision_deepribo, recall_deepribo, base, auc_prc_deepribo = comput_roc_for_tool(experiment_dict_path, 'deepribo')
    precision_ribotish, recall_ribotish, base, auc_prc_ribotish = comput_roc_for_tool(experiment_dict_path, 'ribotish',)
    precision_reparation, recall_reparation, base, auc_prc_reparation = comput_roc_for_tool(experiment_dict_path, 'reparation')
    precision_irsom, recall_irsom, base, auc_prc_irsom = comput_roc_for_tool(experiment_dict_path, 'irsom')

    #print('deepribo AUC: %f' % (auc_deepribo))
    #print('ribotish AUC: %f' %  (auc_ribotish))
    #print('reparation AUC: %f' %(auc_reparation))
    #print('irsom AUC: %f' %  (auc_irsom))


#     label_deepribo = 'deepribo AUC: %f' % (auc_deepribo)
#     label_ribotish = ('ribotish AUC: %f' %  (auc_ribotish))
#     label_reparation = ('reparation AUC: %f' %(auc_reparation))
#     label_irsom = ('irsom AUC: %f' %  (auc_irsom))
#
#     print('++++\nfpr deepribo:')
#     print(fpr_deepribo)
#     print('+++++++++++++++++++++')
#     print('++++\ntpr deepribo:')
#     print(tpr_deepribo)
#     print('+++++++++++++++++++++')
#     plt.plot(fpr_deepribo, tpr_deepribo, color='orange', label=label_deepribo)
#     plt.plot(fpr_ribotish, tpr_ribotish, color='blue', label=label_ribotish)
#     plt.plot(fpr_reparation, tpr_reparation, color='red', label=label_reparation)
#     plt.plot(fpr_irsom, tpr_irsom, color='green', label=label_irsom)
#
#     plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
#     plt.xlabel('False Positive Rate')
#     plt.ylabel('True Positive Rate')
#     plt.title('Receiver Operating Characteristic (ROC) Curve ' + experiment)
#     plt.legend()
# #plt.show()
#     save_roc_diag = plot_dir + experiment + '_roc.pdf'
#     plt.savefig(save_roc_diag, format='pdf', dpi=300, bbox_inches='tight')

    #
    label_deepribo = 'DeepRibo AUC: %.2f' % (auc_prc_deepribo)
    label_ribotish = ('Ribo-TISH AUC: %.2f' %  (auc_prc_ribotish))
    label_reparation = ('Reparation AUC: %.2f' %(auc_prc_reparation))
    label_irsom = ('IRSOM AUC: %.2f' %  (auc_prc_irsom))

    plt.plot(recall_deepribo, precision_deepribo, color='green', label=label_deepribo)
    plt.plot(recall_reparation, precision_reparation, color='blue', label=label_reparation)
    plt.plot(recall_ribotish, precision_ribotish, color='yellow', label=label_ribotish)
    plt.plot(recall_irsom, precision_irsom, color='red', label=label_irsom)

    plt.plot([0, 1], [base, base], color='grey', linestyle='--')
    plt.xlabel('Recall', fontsize=14)
    plt.ylabel('Precision', fontsize=14)
    plt.xticks(fontsize=12, rotation=30)
    plt.yticks(fontsize=12)
    #plt.title('Precision Recall Characteristic (prc) Curve ' + species)
    plt.title(species)

    #fontsize=20
    plt.legend(prop={'size': 14})
    plt.ylim(-0.02, 1.02)
    #plt.show()
    save_prc_diag = plot_dir + species + '_prc_' + overlap + '_' + experiment+'_.pdf'
    plt.savefig(save_prc_diag, format='pdf', dpi=300, bbox_inches='tight')

if __name__ == '__main__':
    main()
