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


def compute_roc(labels, scores):
    #print(labels)
    #print(scores)
    fpr, tpr, thresholds = metrics.roc_curve(labels, scores, pos_label=1)
    auc =metrics.auc(fpr, tpr)
    return fpr, tpr, auc, thresholds


def get_list(saved_dir, tool):
    dir_score_list = saved_dir + '/' + tool + '_score_list'
    #dir_lable_list = saved_dir + '/' + tool + '_label_list'
    if not os.path.isfile(dir_score_list):
        sys.exit('not existing file (ROC): %s' % (dir_score_list))

    with open(dir_score_list, 'rb') as handle:
        score_overlap_label_list = pickle.load(handle)
    #with open(dir_lable_list, 'rb') as handle:
        #label_list = pickle.load(handle)
    #print(len(label_list))
    #print(len(score_overlap_label_list))

    return score_overlap_label_list

def get_ranks(score_overlap_label_list):
    #first_indexes = len(score_list)
    rank_list = list(range(len(score_overlap_label_list)))
    list.reverse(rank_list)

    #diff = len(score_overlap_label_list)-(len(set(score_overlap_label_list)))

    sorted_list = sorted(score_overlap_label_list, key=itemgetter(0,1), reverse=True)

    count = sorted_list.count((0,0,'0'))

    last_score = 0
    last_overlap = 0
    last_rank = 0
    label_list = []
    count_pos = 0
    count_neg = 0

    for pos,trippel in enumerate(sorted_list):
        label_list.append(int(trippel[2]))
        #count number of labels
        if trippel[2] == '1':
            count_pos += 1
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

    return rank_list, label_list



def comput_roc_for_tool(saved_dir, tool):

    score_overlap_label_list = get_list(saved_dir, tool)
    #print('original scores first entry: %i' %(score_overlap_label_list[0][0]))
    #print('length original scores: %i' %(len(score_overlap_label_list)))

    score_list, label_list = get_ranks(score_overlap_label_list)
    #print('length scores: %i' %(len(score_list)))
    #print('length labels: %i' %(len(label_list)))


    fpr, tpr, auc, thresholds = compute_roc(label_list, score_list)

    # get optimal treshold:
    optimal_idx = np.argmin(np.abs(tpr - fpr))
    #print(np.argmin(np.abs(tpr - fpr)))
    #print('optimal index: %f' %(optimal_idx))
    optimal_threshold = thresholds[optimal_idx]
    #print('optimal TH: %i' %(optimal_threshold))

    return fpr, tpr, auc, optimal_threshold

def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_path", action="store", dest="experiment_dict_path", required=True
                                           , help= "path list of trippels (score, overlap, label) for each tool")
    parser.add_argument("-e", "--experiment", action="store", dest="experiment", required=True
                                           , help= "experiment name")
    parser.add_argument("-o", "--save_path", action="store", dest="save_path", required=True
                                      , help= "path to save data to.")
    parser.add_argument("-p", "--percent_overlap", action="store", dest="percent_overlap", required=True
                                      , help= "how much of the prediction should overlap with the Gene and vice verser")

    args = parser.parse_args()
    experiment_dict_path = args.experiment_dict_path
    experiment = args.experiment
    plot_dir = args.save_path
    overlap = args.percent_overlap



    # 'deepribo', 'ribotish', 'reparation', 'irsom'


    fpr_deepribo, tpr_deepribo, auc_deepribo, opt_th_deepribo  = comput_roc_for_tool(experiment_dict_path, 'deepribo')
    fpr_ribotish, tpr_ribotish, auc_ribotish, opt_th_ribotish = comput_roc_for_tool(experiment_dict_path, 'ribotish',)
    fpr_reparation, tpr_reparation, auc_reparation, opt_th_reparation = comput_roc_for_tool(experiment_dict_path, 'reparation')
    fpr_irsom, tpr_irsom, auc_irsom, opt_th_irsom = comput_roc_for_tool(experiment_dict_path, 'irsom')

    #print('deepribo AUC: %f' % (auc_deepribo))
    #print('ribotish AUC: %f' %  (auc_ribotish))
    #print('reparation AUC: %f' %(auc_reparation))
    #print('irsom AUC: %f' %  (auc_irsom))


    #label_deepribo = ('DeepRibo AUC: %.2f and optTH: %i' % (auc_deepribo, opt_th_deepribo))
    #label_reparation = ('Reparation AUC: %.2f and optTH: %i' %(auc_reparation, opt_th_reparation))
    #label_ribotish = ('Ribo-TISH AUC: %.2f and optTH: %i' %  (auc_ribotish, opt_th_ribotish))
    #label_irsom = ('IRSOM AUC: %.2f and optTH: %i' %  (auc_irsom, opt_th_irsom))

    label_deepribo = ('DeepRibo AUC: %.2f' % (auc_deepribo))
    label_reparation = ('Reparation AUC: %.2f' %(auc_reparation))
    label_ribotish = ('Ribo-TISH AUC: %.2f' %  (auc_ribotish))
    label_irsom = ('IRSOM AUC: %.2f' %  (auc_irsom))

    #print('++++\nfpr deepribo:')
    #print(fpr_deepribo)
    #print('+++++++++++++++++++++')
    #print('++++\ntpr deepribo:')
    #print(tpr_deepribo)
    #print('+++++++++++++++++++++')
    plt.plot(fpr_deepribo, tpr_deepribo, color='green', label=label_deepribo)
    plt.plot(fpr_reparation, tpr_reparation, color='blue', label=label_reparation)
    plt.plot(fpr_ribotish, tpr_ribotish, color='yellow', label=label_ribotish)
    plt.plot(fpr_irsom, tpr_irsom, color='red', label=label_irsom)

    plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve ' + experiment)
    plt.legend()
#plt.show()
    save_roc_diag = plot_dir + experiment + '_roc_' + overlap + '_.pdf'
    plt.savefig(save_roc_diag, format='pdf', dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    main()
