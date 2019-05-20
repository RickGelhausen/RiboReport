#!/usr/bin/env python
import pandas as pd
#import mathe
import matplotlib as mpl
import seaborn
import os
import argparse
import csv
import collections
from operator import itemgetter
import sys


#filter on

score_cutoff = 1


def get_dict_with_empty_val(file, tool='annotation'):
    tool_dict = {}
    flag_header = 0
    f = open(file, 'r') # open the file for reading
    data = []
    for row_num, line in enumerate(f):
        values = line.strip().split('\t')
        #print values

        if flag_header == 1:
            header = values
            flag_header = 0
            #print header
            continue
        elif values[1] == tool:
            #print values
            key = values[3] + ':' + values[4] + ':' + values[6]
            #print key
            value = []
            #print value
            tool_dict[key] = value
        else:
            value = []
            #print('ceck if: ' + tool + ' is not ' + values[1] + '\n')
    return tool_dict

# assosiate genes:
def get_assosiated_genes(gen_dict, pred_predition_dict):
    assosiated_genes_dict = {}
    assosiated_predictions_dict = {}
    pred_predition_sub_dict = pred_predition_dict
    #print(gen_dict)
    #print(len(gen_dict))
    for k in gen_dict:
        #print(k)
        sorted_list = sorted(gen_dict[k], key=itemgetter(1))
        #print(sorted_list)
        i = 0
        # find gene:
        while len(sorted_list) != 0:
            id_pred = sorted_list.pop(0)[0]
            #print(id_pred)
            if id_pred in pred_predition_sub_dict.keys():
                #print(pred_predition_sub_dict.keys())
                try:
                    del pred_predition_sub_dict[id_pred]
                except KeyError:
                    pass
                assosiated_genes_dict[k] = id_pred
                assosiated_predictions_dict[id_pred]= k
                sorted_list = []
                #print(assosiated_genes_dict)
        #print(pred_predition_sub_dict)
    #print(assosiated_genes_dict)
    #print(pred_predition_sub_dict)
    return assosiated_genes_dict, pred_predition_sub_dict, assosiated_predictions_dict

def get_genes(gene_dict, df_overlap_score, tool):
    # get all assosiation of genes with a specific tool
    df_overlap = df_overlap_score[df_overlap_score['id_pred']==tool]

    #prepare the value of the dict in an extra colum of the dataframe
    df_overlap_test = df_overlap.astype({"start_pred": str, "stop_pred": str, "score_pred": str, "overlap_percent": str})
    df_overlap_test['tuple'] = df_overlap_test['start_pred']+ ':' + df_overlap_test['stop_pred'] +  ':' + df_overlap_test['strand_pred'] + '$' + df_overlap_test['score_pred'] + '$' + df_overlap_test['overlap_percent']

    # generate gene_prediction dataframe
    gene_predition_temp_dict = {k: list(v) for k, v in df_overlap_test.groupby(by=['start_ref', 'stop_ref', 'strand_ref'])['tuple']}

    gene_predition_dict ={str(k[0]) + ':' + str(k[1]) + ':' + k[2]: gene_predition_temp_dict[k] for k in gene_predition_temp_dict}

    for k in gene_predition_dict:
        gene_predition_dict[k] = [val.split('$') for val in gene_predition_dict[k]]

    # compute what genes are not predicted
    #not_pred_dict = gene_dict.keys() - gene_predition_dict.keys()

    return gene_predition_dict

def get_prediction(tool_dict, df_overlap_score, tool):
    # get all assosiation of genes with a specific tool
    df_overlap_tools = df_overlap_score[df_overlap_score['id_pred']==tool]

    #prepare the value of the dict in an extra colum of the dataframe
    df_overlap_test = df_overlap_tools.astype({"start_ref": str, "stop_ref": str, "overlap_percent": str})
    df_overlap_test['tuple'] = df_overlap_test['start_ref']+ ':' + df_overlap_test['stop_ref'] +  ':' + df_overlap_test['strand_ref'] + '-' + df_overlap_test['overlap_percent']

    # generate gene_prediction dataframe
    pred_predition_temp_dict = {k: list(v) for k, v in df_overlap_test.groupby(by=['start_pred', 'stop_pred', 'strand_pred'])['tuple']}

    pred_predition_dict ={str(k[0]) + ':' + str(k[1]) + ':' + k[2]: pred_predition_temp_dict[k] for k in pred_predition_temp_dict}

    for k in pred_predition_dict:
        pred_predition_dict[k] = [val.split('-') for val in pred_predition_dict[k]]

    # compute what genes are not predicted
    #not_pred_dict = tool_dict.keys() - pred_predition_dict.keys()

    return pred_predition_dict

def get_stat_for_tool(ref_file, tools_file, df_overlap_score, tool):

    # put all gegenes of the refernece inside a dictionary
    gene_dict = get_dict_with_empty_val(ref_file)
    #print(gene_dict)
    #print('\n')
    # put all predictions insed a dictionary
    tool_dict = get_dict_with_empty_val(tools_file, tool=tool)
    #print(tool_dict)
    #print('\n')
    # get all genes which are found by the tools prediction and the genes that are not found!
    gene_predicted_dict = get_genes(gene_dict, df_overlap_score, tool)
    #print(df_overlap_score)
    print(gene_predicted_dict)
    #print('\n')
    #get all predicitons of the tools that have a overlap with a gene of the reference
    predition_overlap_dict = get_prediction(tool_dict, df_overlap_score, tool)
    #print(predition_overlap_dict)
    #print('\n')
    # assosiatie each gene with just one prediciton. All othere genes
    assosiated_genes_dict, pred_predition_sub_list, assosiated_predictions_dict = get_assosiated_genes(gene_predicted_dict, predition_overlap_dict)
    #print(assosiated_genes_dict)
    #print('\n')
    #print(pred_predition_sub_list)
    #print('\n')


    #if not len(assosiated_genes_dict):
        #TP = 0
    #else:
    TP = len(assosiated_genes_dict)


    FP = len(tool_dict) - len(assosiated_predictions_dict) - len(pred_predition_sub_list)
    FP_with_suboptimals = len(tool_dict) - len(assosiated_predictions_dict)



    FN = len(gene_dict)  - len(assosiated_genes_dict)
    FP_prediction_list = set(tool_dict.keys()) - set(assosiated_predictions_dict.keys()) - set(pred_predition_sub_list.keys())
    FN_prediction_list = set(gene_dict.keys())  - set(assosiated_genes_dict.keys())


    if TP+FN != 0:
        recall = TP / (TP + FN)
    else:
        recall = 0 
	     
    if TP+FN != 0:        
        FNR = FN / (TP + FN)
    else:
        FNR = 0 
	   
    if TP+FP!= 0:
        precision = TP / (TP+FP)
    else:
        precision = 0
	    
    if TP+FP != 0:
        FDR = FP / (TP+FP)
    else:
        FDR = 0
	    
    if recall+precision != 0:
        F1 = (2* (recall*precision)) / (recall+precision)
    else:
        F1 = 0
        


    #if not pred_predition_sub_list:
        #subset = 0
    #else:
    suboptimals = len(pred_predition_sub_list)
    print("TOOL: %s, SUBOPTIMALS: %s" % (tool, suboptimals))
    #print('TP '+ str(TP) +'\n')
    #print('FP '+ str(FP) +'\n')
    #print('FP with suboptimals '+ str(FP_with_suboptimals) +'\n')
    #print('FN '+ str(FN) +'\n')
    #print('suboptimals '+ str(suboptimals) +'\n')

    gene_list = assosiated_genes_dict.keys()
    prediction_list = assosiated_predictions_dict.keys()
    stat_list = [TP, FP, FN, recall, FNR, precision, FDR, F1, len(pred_predition_sub_list), tool]
    return gene_list, prediction_list, stat_list, FP_prediction_list, FN_prediction_list


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-r", "--reference_data", action="store", dest="reference_path", required=True, default="./data/all_filtered_massspec.gtf"
                                           , help= "path to the reference data.")
    parser.add_argument("-t", "--tool_data", action="store", dest="tools_path", required=True, default="./data/combined.gtf"
                                      , help= "path to the tools data.")
    parser.add_argument("-o", "--save_path", action="store", dest="save_path", required=True, default="./result_dfs/"
                                      , help= "path to save data to.")
    parser.add_argument("-c", "--overlap_cutoff", action="store", dest="overlap_cutoff", required=True, default="./result_dfs/",
                                      type=float, help= "path to save data to.")



    args = parser.parse_args()
    overlap_cutoff= args.overlap_cutoff

    save_dir = args.save_path

    ref_file = args.reference_path
    tools_file = args.tools_path

    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)

    if not os.path.exists(args.tools_path):
        os.makedirs(args.tools_path)
    overlap_gtf = os.path.splitext(os.path.basename(args.reference_path))[0] +"_"+ os.path.splitext(os.path.basename(args.tools_path))[0] + "_overlap.gtf"
    cmd = 'bedtools intersect -a '+ ref_file +' -b ' + tools_file + ' -wo -s > ' + overlap_gtf
    os.system(cmd)
    
    #test if gtf file is emty:
    if os.stat(overlap_gtf).st_size == 0:
        sys.exit('no overlap was found for %s' %(overlap_gtf))
    #read gtf overlap file
    df_temp = pd.read_csv(overlap_gtf,header=None, sep="\t", comment="#")
    df_overlap = pd.DataFrame(df_temp.values, columns = ['chr','id_ref', 'feat_ref','start_ref','stop_ref',
                                                 'score_ref','strand_ref','i_ref','h_ref','chr_pred','id_pred',
                                                 'feat_pred','start_pred','stop_pred','score_pred',
                                                 'strand_pred','r_pred','s_pred', 'overlap'])

    df_overlap['ref_length'] = df_overlap['stop_ref'] - df_overlap['start_ref'] +1
    df_overlap['overlap_percent'] =  df_overlap['overlap']/df_overlap['ref_length']

    df_overlap = df_overlap.astype({"start_ref": int, "stop_ref": int, "start_pred": int, "stop_pred": int, "score_pred": float,  "overlap": int})

    df_overlap = df_overlap[df_overlap['overlap_percent']>= overlap_cutoff]

    df_overlap_score = df_overlap[df_overlap['score_pred']<= score_cutoff]

    if len(df_overlap_score[df_overlap_score['strand_ref'] != df_overlap_score['strand_pred']]) == 0:
        print('strands are equal')
    else:
        print('problem strands are not equal')

    #, FP_prediction_deepribo_list
    gene_deepribo_list, prediction_deepribo_list, stat_deepribo_list, FP_prediction_deepribo_list, FN_prediction_deepribo_list = get_stat_for_tool(ref_file, tools_file, df_overlap_score, 'deepribo')
    gene_ribotish_list, prediction_ribotish_list, stat_ribotish_list, FP_prediction_ribotish_list, FN_prediction_ribotish_list = get_stat_for_tool(ref_file, tools_file, df_overlap_score, 'ribotish')
    gene_reparation_list, prediction_reparation_list, stat_reparation_list, FP_prediction_reparation_list, FN_prediction_reparation_list = get_stat_for_tool(ref_file, tools_file, df_overlap_score, 'reparation')
    gene_irsom_list, prediction_irsom_list, stat_irsom_list, FP_prediction_irsom_list, FN_prediction_irsom_list = get_stat_for_tool(ref_file, tools_file, df_overlap_score, 'irsom')

    tims_deepribo_gen = max([len(gene_deepribo_list),len(gene_ribotish_list),len(gene_reparation_list),len(gene_irsom_list)]) - len(gene_deepribo_list)
    tims_ribotish_gen = max([len(gene_deepribo_list),len(gene_ribotish_list),len(gene_reparation_list),len(gene_irsom_list)]) - len(gene_ribotish_list)
    tims_reparation_gen = max([len(gene_deepribo_list),len(gene_ribotish_list),len(gene_reparation_list),len(gene_irsom_list)]) - len(gene_reparation_list)
    tims_irsom_gen = max([len(gene_deepribo_list),len(gene_ribotish_list),len(gene_reparation_list),len(gene_irsom_list)]) - len(gene_irsom_list)


    tims_deepribo_pred = max([len(prediction_deepribo_list),len(prediction_ribotish_list),len(prediction_reparation_list),len(prediction_irsom_list)]) - len(prediction_deepribo_list)
    tims_ribotish_pred = max([len(prediction_deepribo_list),len(prediction_ribotish_list),len(prediction_reparation_list),len(prediction_irsom_list)]) - len(prediction_ribotish_list)
    tims_reparation_pred = max([len(prediction_deepribo_list),len(prediction_ribotish_list),len(prediction_reparation_list),len(prediction_irsom_list)]) - len(prediction_reparation_list)
    tims_irsom_pred = max([len(prediction_deepribo_list),len(prediction_ribotish_list),len(prediction_reparation_list),len(prediction_irsom_list)]) - len(prediction_irsom_list)


    tims_deepribo_not_pred = max([len(FP_prediction_deepribo_list),len(FP_prediction_ribotish_list),len(FP_prediction_reparation_list),len(FP_prediction_irsom_list)]) - len(FP_prediction_deepribo_list)
    tims_ribotish_not_pred = max([len(FP_prediction_deepribo_list),len(FP_prediction_ribotish_list),len(FP_prediction_reparation_list),len(FP_prediction_irsom_list)]) - len(FP_prediction_ribotish_list)
    tims_reparation_not_pred = max([len(FP_prediction_deepribo_list),len(FP_prediction_ribotish_list),len(FP_prediction_reparation_list),len(FP_prediction_irsom_list)]) - len(FP_prediction_reparation_list)
    tims_irsom_not_pred = max([len(FP_prediction_deepribo_list),len(FP_prediction_ribotish_list),len(FP_prediction_reparation_list),len(FP_prediction_irsom_list)]) - len(FP_prediction_irsom_list)


    tims_deepribo_not_gene = max([len(FN_prediction_deepribo_list),len(FN_prediction_ribotish_list),len(FN_prediction_reparation_list),len(FN_prediction_irsom_list)]) - len(FN_prediction_deepribo_list)
    tims_ribotish_not_gene = max([len(FN_prediction_deepribo_list),len(FN_prediction_ribotish_list),len(FN_prediction_reparation_list),len(FN_prediction_irsom_list)]) - len(FN_prediction_ribotish_list)
    tims_reparation_not_gene = max([len(FN_prediction_deepribo_list),len(FN_prediction_ribotish_list),len(FN_prediction_reparation_list),len(FN_prediction_irsom_list)]) - len(FN_prediction_reparation_list)
    tims_irsom_not_gene = max([len(FN_prediction_deepribo_list),len(FN_prediction_ribotish_list),len(FN_prediction_reparation_list),len(FN_prediction_irsom_list)]) - len(FN_prediction_irsom_list)

    venn_genes_dict = {'deepribo': list(gene_deepribo_list) + ['NA'] * tims_deepribo_gen,
         'ribotish': list(gene_ribotish_list) + ['NA'] * tims_ribotish_gen,
         'reparation': list(gene_reparation_list) + ['NA'] * tims_reparation_gen,
         'irsom': list(gene_irsom_list) + ['NA'] * tims_irsom_gen}
    df_venn_genes = pd.DataFrame(venn_genes_dict)
    df_venn_genes.to_csv(path_or_buf=save_dir+'df_venn_genes.csv', index=False, sep='\t')

    venn_predictions_dict = {'deepribo': list(prediction_deepribo_list) + ['NA'] * tims_deepribo_pred,
         'ribotish': list(prediction_ribotish_list) + ['NA'] * tims_ribotish_pred,
         'reparation': list(prediction_reparation_list) + ['NA'] * tims_reparation_pred,
         'irsom': list(prediction_irsom_list) + ['NA'] * tims_irsom_pred}
    df_venn_predictions = pd.DataFrame(venn_predictions_dict)
    df_venn_predictions.to_csv(path_or_buf=save_dir+'df_venn_predictions.csv', index=False, sep='\t')



    venn_FP_predictions_dict = {'deepribo': list(FP_prediction_deepribo_list) + ['NA'] * tims_deepribo_not_pred,
         'ribotish': list(FP_prediction_ribotish_list) + ['NA'] * tims_ribotish_not_pred,
         'reparation': list(FP_prediction_reparation_list) + ['NA'] * tims_reparation_not_pred,
         'irsom': list(FP_prediction_irsom_list) + ['NA'] * tims_irsom_not_pred}
    df_venn_FP_predictions = pd.DataFrame(venn_FP_predictions_dict)
    df_venn_FP_predictions.to_csv(path_or_buf=save_dir+'df_venn_FP_predictions_dict.csv', index=False, sep='\t')

    venn_FN_gene_dict = {'deepribo': list(FN_prediction_deepribo_list) + ['NA'] * tims_deepribo_not_gene,
         'ribotish': list(FN_prediction_ribotish_list) + ['NA'] * tims_ribotish_not_gene,
         'reparation': list(FN_prediction_reparation_list) + ['NA'] * tims_reparation_not_gene,
         'irsom': list(FN_prediction_irsom_list) + ['NA'] * tims_irsom_not_gene}
    df_venn_FN_gene = pd.DataFrame(venn_FN_gene_dict)
    df_venn_FN_gene.to_csv(path_or_buf=save_dir+'df_venn_FN_gene_dict.csv', index=False, sep='\t')


    stat_list_header = ['TP', 'FP', 'FN', 'recall', 'FNR', 'precision', 'FDR', 'F1', 'subopt', 'tool']
    df_stat = pd.DataFrame([stat_reparation_list, stat_ribotish_list, stat_deepribo_list, stat_irsom_list], columns=stat_list_header)
    
    



    #pd.concat([df_stat, stat_reparation_list])
    df_stat
    df_stat.to_csv(path_or_buf=save_dir+'df_stat.csv', index=False, sep='\t')

if __name__ == '__main__':
    main()