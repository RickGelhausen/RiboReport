#!/usr/bin/env python
import pandas as pd
import matplotlib as mpl
import seaborn
import os
import argparse
import csv
import collections
from operator import itemgetter
import sys
import pickle
import glob


#filter on
# 1000000000 as cutoff is higher than all scores the tools would produce! Therefore the score_cutoff will be not used!
score_cutoff = 1000000000


def get_dict_with_empty_val(file, tool='RefSeq'):
    """Produces a dictionary with a key value pairs for every line of
    the given gtf file which contains the keyword tool.

    Parameters
    ----------
    tool : string
        annotation or either of the tool names

    Returns
    -------
    tool_dict
        dictionary with ID as key and value is the string parameter tool
        id of the start:stop:strand:assertionnumber
    """
    tool_dict = {}
    f = open(file, 'r')
    data = []
    for row_num, line in enumerate(f):
        values = line.strip().split('\t')
        #print values

        if values[1] == tool:
            #print values
            key = values[3] + ':' + values[4] + ':' + values[6] + ':' + values[0]
            #print key
            value = []
            #print value
            tool_dict[key] = value

    return tool_dict


def get_associated_genes(gene_dict, pred_prediction_dict, label):
    """function associates just one prediction to a gene.

        Parameters
        ----------
        gene_dict : dictionary
            geneIds with list of predictions
        pred_prediction_dict : dictionary
            predictionIds with lists of genes

        Returns
        -------
        associated_genes_dict
            associations of geneId -> predictionId (one to one relation)
        pred_prediction_sub_dict
            containing predictionsIds that are not associated to a gene
        associated_predictions_dict
            associations of predictionId -> geneId (one to one relation)
    """
    associated_genes_dict = {}
    associated_predictions_dict = {}
    pred_prediction_sub_dict = pred_prediction_dict
    score_list = []
    fn = 0

    for k in gene_dict:

        # sort list of predictions according to there score
        sorted_list = sorted(gene_dict[k], key=itemgetter(2), reverse=True)

        i = 0
        # find gene:
        while len(sorted_list) != 0:
            pred_info_list = sorted_list.pop(0)
            id_pred = pred_info_list[0]
            score = pred_info_list[1]
            overlap = pred_info_list[2]

            if id_pred in pred_prediction_sub_dict.keys():
                #print(id_pred)
                try:
                    del pred_prediction_sub_dict[id_pred]
                except KeyError:
                    pass
                associated_genes_dict[k] = id_pred
                associated_predictions_dict[id_pred]= k
                score_list.append((float(score), float(overlap), label))

                sorted_list = []
            elif sorted_list == []:
                fn += 1
                #score_list.append((0,0,label))
                #print(associated_genes_dict)
        #print(pred_prediction_sub_dict)


    return (associated_genes_dict, pred_prediction_sub_dict,
            associated_predictions_dict, fn, score_list)


def get_genes(gene_dict, df_overlap_score, tool):
    """function selects for a saves tool all predicted genes (key) and there
    associated predictions (value) in a dictionary

        Parameters
        ----------
        gene_dict : dictionary
            containing all gene IDs (not needed at the moment)
        df_overlap_score : datafame
            containing overlaps between the reference and prediction
        tool : string
            tells for what tool the overlap should be found

        Returns
        -------
        gene_prediction_dict
            A dictionary: {'geneId1': [['predictionId1', 'score', 'overlap'],
            ['predictionId2', 'score', 'overlap']...] 'geneId2':
            [['predictionId3', 'score', 'overlap'],...]...}
            Id: start:stop:strand
    """
    # get all association of genes with a specific tool
    df_overlap = df_overlap_score[df_overlap_score['id_pred'] == tool]

    #prepare the value of the dict in an extra column of the dataframe
    df_overlap_test = df_overlap.astype({"start_pred": str, "stop_pred": str,
                                        "score_pred": str,
                                        "overlap_percent": str,
                                        "chr":str})
    df_overlap_test['tuple'] = (df_overlap_test['start_pred'] + ':' +
                                df_overlap_test['stop_pred'] + ':' +
                                df_overlap_test['strand_pred'] + ':' +
                                df_overlap_test['chr']+ '$' +
                                df_overlap_test['score_pred'] + '$' +
                                df_overlap_test['overlap_percent'])
    #print(df_overlap_test)

    # dict key: gene_id; value: list of overlapping prediction (predID$score$overlap)
    gene_prediction_temp_dict = {k: list(v)
                                for k, v in df_overlap_test.groupby(
                                by=['start_ref', 'stop_ref', 'strand_ref', 'chr'])
                                ['tuple']}
    #print(gene_prediction_temp_dict)

    # making a sting out of the tuple saved as key
    gene_prediction_dict ={str(k[0]) + ':' + str(k[1]) + ':' +
                          k[2] + ':' + k[3]: gene_prediction_temp_dict[k]
                          for k in gene_prediction_temp_dict}

    for gene_id in gene_prediction_dict:
        gene_prediction_dict[gene_id] = [val.split('$')
                                        for val in gene_prediction_dict
                                        [gene_id]]

    # compute what genes are not predicted
    # not_pred_dict = gene_dict.keys() - gene_prediction_dict.keys()

    return gene_prediction_dict


def get_prediction(tool_dict, df_overlap_score, tool):
    """This fuction selects for a given tool all predictions (key) and there
    associated genes (value) and stores them in a dictionary

        Parameters
        ----------
        tool_dict : dictionary (not needed at the moment)
            keys containing all predictions ids but jet no values
        df_overlap_score : datafame
            containing overlaps between the reference and prediction
        tool : string
            tells for what tool the overlap should be found


        Returns
        -------
        pred_prediction_dict
            A dictionary: {'predId1': [['geneId1', 'overlap'],
            ['geneId2', ' 'overlap']...] 'predId2':
            [['geneId3', 'overlap'],...]...}
            Id: start:stop:strand
    """
    # get all association of genes with a specific tool
    df_overlap_tools = df_overlap_score[df_overlap_score['id_pred'] == tool]

    # prepare the value of the dict in an extra column of the dataframe
    df_overlap_test = df_overlap_tools.astype({
                                             "start_ref": str, "stop_ref": str,
                                             "overlap_percent": str})
    df_overlap_test['tuple'] = (df_overlap_test['start_ref'] + ':' +
                                df_overlap_test['stop_ref'] + ':' +
                                df_overlap_test['strand_ref'] + ':' +
                                df_overlap_test['chr'] + '-' +
                                df_overlap_test['overlap_percent'])

    # generate gene_prediction dataframe
    pred_prediction_temp_dict = {pred_Id: list(gene_ids) for pred_Id, gene_ids in
                                df_overlap_test.groupby(by=['start_pred',
                                'stop_pred', 'strand_pred', 'chr'])
                                ['tuple']}
    # convert pred_id tuple into one string
    pred_prediction_dict = {str(k[0]) + ':' + str(k[1]) + ':' +
                           k[2] + ':' + k[3]: pred_prediction_temp_dict[k]
                           for k in pred_prediction_temp_dict}


    for k in pred_prediction_dict:
        # split gene information stored in values into gene_id and overlap_percent
        pred_prediction_dict[k] = [val.split('-')
                                  for val in pred_prediction_dict[k]]

    # compute what genes are not predicted
    #not_pred_dict = tool_dict.keys() - pred_prediction_dict.keys()

    return pred_prediction_dict


def assort_genes_with_predictions(ref_file, tools_file, df_overlap_score, tool, label):
    """This function computes TP, FN and FP for a given tool.

        Parameters
        ----------
        ref_file : string
            file name of the reference file (.gtf)
        tools_file : string
            file name of the tools prediction (.gtf)
        df_overlap_score : dataframe
            bedtools produced overlaps between tool predictions and genes
        tool : string
            tool name
        label: string
            specifies if positive or negative data is analyzed


        Returns
        -------
        assorted_genes
            number of genes associated with one prediction
        not_assorted_genes
            number of genes without prediction tools_file
        suboptimals
            number of suboptimal predictions
        associated_gene_list
            list of all genes found by the currently investigated tool
        score_list
            list containing for the current tool (score, overlap, label) for each associated gene

    """
    # put all genes of the reference inside a dictionary
    gene_dict = get_dict_with_empty_val(ref_file)
    print('number of genes for %s: %i' % (label, len(gene_dict)))

    # put all predictions inside a dictionary
    tool_dict = get_dict_with_empty_val(tools_file, tool=tool)

    # get all genes which are found by the tools prediction and the genes that are not found!
    gene_predicted_dict = get_genes(gene_dict, df_overlap_score, tool)

    # get all predictions of the tools that have a overlap with a gene of the reference
    prediction_overlap_dict = get_prediction(tool_dict, df_overlap_score, tool)

    # associate each gene with just one prediction. All other genes
    associated_genes_dict, pred_prediction_sub_list, associated_predictions_dict, no_genes_lost_prediction, score_list = get_associated_genes(gene_predicted_dict, prediction_overlap_dict,label)


    assert (len(associated_genes_dict) == len(associated_predictions_dict)),\
            "list of overlapping genes and prediction is not equal"

    if not len(associated_genes_dict):
        assorted_genes = 0
        not_assorted_genes = len(gene_dict)
    else:
        assorted_genes = len(associated_predictions_dict)
        not_assorted_genes = len(gene_dict) - len(associated_predictions_dict)

    if not len(pred_prediction_sub_list):
        suboptimals = 0
    else:
        suboptimals = len(pred_prediction_sub_list)

    associated_gene_list =  list(associated_genes_dict.keys())

    return assorted_genes, not_assorted_genes, suboptimals, associated_gene_list, score_list


def get_stat_for_tool(ref_file, tools_file, df_overlap_score, ref_neg_file, df_neg_overlap, df_fp_overlap, flag_subopt, flag_no_gene, tool, save_dir):
    """This fuction computes TP, FN and FP for a given tool.

        Parameters
        ----------
        ref_file : string
            file name of the reference file (.gtf)
        tools_file : string
            file name of the tools prediction (.gtf)
        df_overlap_score : dataframe
            bedtools produced overlaps between tool predictions and genes
        tool : string
            tool name
        save_dir:
            directory where the list of scores overlaps and labels can be stores at


        Returns
        -------
        associated_TP_list
            list of genes that are found (TP)
        stat_list
            contains all statistical measurements:
            TP, FP, FN, recall, FNR, precision, FDR, F1, #subopt, tool]
    """
    TP, FN, suboptimals_TP, associated_TP_list, score_list_pos = assort_genes_with_predictions(ref_file, tools_file, df_overlap_score, tool, '1')
    FP, TN, suboptimals_FP, associated_FP_list, score_list_neg = assort_genes_with_predictions(ref_neg_file, tools_file, df_neg_overlap, tool, '0')
    #print(score_list_pos)

    # get list for ROC curve containing (score, overlap, label)
    save_score_list = save_dir + '/' + tool + '_score_list'
    #save_label_list = save_dir + '/' + tool + '_label_list'

    # generate (score, overlap, label) entrys for all genes which had no prediction
    if tool == 'deepribo':
        FN_pred_ROC = [(-999999,0,'1')]*FN
        TN_pred_ROC = [(-999999,0,'0')]*TN
    elif tool == 'price':
        FN_pred_ROC = [(1,0,'1')]*FN
        TN_pred_ROC = [(1,0,'0')]*TN
    else:
        FN_pred_ROC = [(0,0,'1')]*FN
        TN_pred_ROC = [(0,0,'0')]*TN

    #labels_list= label_pos_list + label_neg_list
    score_list = score_list_pos + score_list_neg + FN_pred_ROC + TN_pred_ROC

    # save the overlaps as input for a ROC curve
    with open(save_score_list, 'wb') as handle:
        pickle.dump(score_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

    df_fp_overlap_tool = df_fp_overlap[df_fp_overlap['id_pred'] == tool]
    #print(df_fp_overlap_tool.info)
    #print(len(df_fp_overlap_tool.index))

    # print(df_fp_overlap_tool)
    no_prediction_without_genes = len(df_fp_overlap_tool.index)

    if(flag_subopt == 1 and flag_no_gene == 1):
        FP = FP + int(suboptimals_TP) + int(no_prediction_without_genes)
    elif (flag_subopt == 0 and flag_no_gene == 1):
        FP = FP + int(no_prediction_without_genes)
    elif (flag_subopt == 1 and flag_no_gene == 0):
        FP = FP + int(suboptimals_TP)

    if TP+FN != 0:
        recall = TP / (TP + FN)
    else:
        recall = 0

    if TP+FN != 0:
        FNR = FN / (TP + FN)
    else:
        FNR = 0

    if TP+FP != 0:
        precision = TP / (TP+FP)
    else:
        precision = 0

    if TP+FP != 0:
        FDR = FP / (TP+FP)
    else:
        FDR = 0

    if recall+precision != 0:
        F1 = (2 * (recall*precision)) / (recall+precision)
    else:
        F1 = 0

    if (TP+FP+TN+FN) != 0:
        accuracy = (TP+TN) / (TP+FP+TN+FN)
    else:
        accuracy = 0

    if (TN+FP) != 0:
        specificity = TN/(TN+FP)
    else:
        specificity = 0


    # #if not pred_prediction_sub_list:
    #     #subset = 0
    # #else:
    # suboptimals = len(pred_prediction_sub_list)
    print("+++++++\nTOOL: %s, SUBOPTIMALS: %s\n" % (tool, suboptimals_TP))
    print('TP '+ str(TP) +'\n')
    print('FP '+ str(FP) +'\n')
    #print('FP with suboptimals '+ str(FP_with_suboptimals) +'\n')
    print('FN '+ str(FN) +'\n')
    print('TN '+ str(TN) +'\n')
    print('suboptimals FP '+ str(suboptimals_FP) +'\n')
    print('predictions not in genes '+ str(no_prediction_without_genes) +'\n+++++++++\n')
    no_prediction_without_genes
    #
    # gene_list = associated_genes_dict.keys()
    # prediction_list = associated_predictions_dict.keys()
    stat_list = [TP, TN, FP, FN, recall, specificity, FNR, precision, FDR, F1, accuracy,
                 suboptimals_TP, suboptimals_FP, no_prediction_without_genes, tool]
    return stat_list, associated_TP_list

def call_bedtools_overlap(reference_path, tools_path, overlap_cutoff, overlap_file, score_cutoff):
    """This function computes TP, FN and FP for a given tool.

        Parameters
        ----------
        reference_path : string
            file name of the reference file (.gtf)
        tools_path : string
            file name of the tools prediction (.gtf)
        overlap_file : string
            file name of the output overlaps (.gtf)
        score_cutoff : float
            gives


        Returns
        -------
        df_overlap_score
            dataframe containing overlapping regions
    """

    if os.stat(reference_path).st_size == 0:
        sys.exit('empty file: %s' % (reference_path))
    #print(which_set + '\n')
    cmd = ('bedtools intersect -a ' + reference_path + ' -b ' + tools_path + ' -wo -f ' + str(overlap_cutoff) + ' -r -s > ' + overlap_file)
    #print('##########')
    #print(cmd)
    #print('##########')

    os.system(cmd)

# read gtf overlap file

    if os.stat(overlap_file).st_size == 0:
        sys.exit('no overlap (th %f) was found for %s \n and tool %s' % (overlap_cutoff, overlap_file, tools_path))

    df_temp = pd.read_csv(overlap_file, header=None, sep="\t", comment="#")
    df_overlap = pd.DataFrame(df_temp.values, columns=['chr', 'id_ref',
                                                    'feat_ref', 'start_ref',
                                                    'stop_ref', 'score_ref',
                                                    'strand_ref', 'i_ref',
                                                    'h_ref', 'chr_pred',
                                                    'id_pred', 'feat_pred',
                                                    'start_pred', 'stop_pred',
                                                    'score_pred',
                                                    'strand_pred', 'r_pred',
                                                    's_pred', 'overlap'])

    df_overlap['ref_overlap_percent'] = (df_overlap['overlap'] /(df_overlap['stop_ref'] -
                                df_overlap['start_ref'] + 1))
    df_overlap['pred_overlap_percent'] = (df_overlap['overlap'] /(df_overlap['stop_pred'] -
                                df_overlap['start_pred'] + 1))
    #df_overlap['overlap_percent'] = min((df_overlap['overlap'] /df_overlap['ref_length']),(df_overlap['overlap'] /df_overlap['ref_length']))
    df_overlap['overlap_percent'] = df_overlap[['ref_overlap_percent','pred_overlap_percent']].min(axis =1)
    #print(df_overlap['score_pred'])
    df_overlap = df_overlap.astype({"start_ref": int, "stop_ref": int,
                                    "start_pred": int, "stop_pred": int,
                                    "score_pred": float,  "overlap": int})

    df_overlap_score = df_overlap[df_overlap['score_pred'] <= score_cutoff]

    # error handeling. But should not occure because of -s parameter of bedtools
    if (len(df_overlap_score[df_overlap_score['strand_ref'] !=
            df_overlap_score['strand_pred']]) == 0):
        print('strands are equal')
    else:
        print('problem strands are not equal')


    return df_overlap_score


def call_bedtools_not_overlap(reference_path, tools_path, overlap_cutoff, overlap_file, score_cutoff):
    """This function computes all genes or predictions (A) which are not overlapping with predictions or genes (B). Given a minimal overlap percentage.
    The result is computed using bedtools and the result is saved in a gtf format and a dictionary

        Parameters
        ----------
        reference_path : string
            file name of the reference file (.gtf)
        tools_path : string
            file name of the tools prediction (.gtf)
        overlap_file : string
            file name of the output overlaps (.gtf)
        score_cutoff : float
            gives


        Returns
        -------
        df_overlap_score
            dataframe containing not overlapping regions
    """

    #bedtools parameter
    # -f: Minimum overlap required as a fraction of A. Default is 1E-9 (i.e. 1bp).
    # -r: Require that the fraction of overlap be reciprocal for A and B. In other words, if -f is 0.90 and -r is used, this requires that B overlap at least 90% of A and that A also overlaps at least 90% of B.
    # -s: Force “strandedness”. That is, only report hits in B that overlap A on the same strand. By default, overlaps are reported without respect to strand.
    # -v: Only report those entries in A that have no overlap in B. Restricted by -f and -r.
    # -wo: Write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported. Restricted by -f and -r
    #
    cmd = ('bedtools intersect -a ' + reference_path + ' -b '  + tools_path + ' -wo -f ' + str(overlap_cutoff) + ' -r -s -v > ' + overlap_file)

    os.system(cmd)

    # test if gtf file is emty:
    if os.stat(overlap_file).st_size == 0:
        print('**********************\nno prediction exits that could not be associated with any gene\n**********************\n')
        df_overlap_score = score_cutoff
    else:
        df_temp = pd.read_csv(overlap_file, header=None, sep="\t", comment="#")
        df_overlap = pd.DataFrame(df_temp.values, columns=['chr_pred',
                                                        'id_pred', 'feat_pred',
                                                        'start_pred', 'stop_pred',
                                                        'score_pred',
                                                        'strand_pred', 'r_pred',
                                                        's_pred'])
        df_overlap = df_overlap.astype({"start_pred": int, "stop_pred": int,
                                        "score_pred": float})
        df_overlap_score = df_overlap

    return df_overlap_score



def get_times_gene_value(no_pred_tool_val):
    times_gene = (max([len(gene_deepribo_list),
                        len(gene_ribotish_list),
                        len(gene_reparation_list),
                        len(gene_irsom_list),
                        len(gene_spectre_list),
                        len(gene_price_list),
                        len(gene_ribotricer_list),
                        len(gene_smorfer_list)]) - no_pred_tool_val)
    return times_gene


def venn_prepataions(gene_deepribo_list, gene_ribotish_list, gene_reparation_list, gene_irsom_list, gene_spectre_list, gene_price_list, gene_ribotricer_list, gene_smorfer_list, save_dir):
    max_no_of_genes = max([len(gene_deepribo_list),
                        len(gene_ribotish_list),
                        len(gene_reparation_list),
                        len(gene_irsom_list),
                        len(gene_spectre_list),
                        len(gene_price_list),
                        len(gene_ribotricer_list),
                        len(gene_smorfer_list)])
    tims_deepribo_gen = max_no_of_genes - len(gene_deepribo_list)
    tims_ribotish_gen = max_no_of_genes - len(gene_ribotish_list)
    tims_reparation_gen = max_no_of_genes - len(gene_reparation_list)
    tims_irsom_gen =  max_no_of_genes - len(gene_irsom_list)
    tims_spectre_gen = max_no_of_genes - len(gene_spectre_list)
    tims_price_gen = max_no_of_genes - len(gene_price_list)
    tims_ribotricer_gen =max_no_of_genes - len(gene_ribotricer_list)
    tims_smorfer_gen = max_no_of_genes - len(gene_smorfer_list)

    venn_genes_dict = {'deepribo': list(gene_deepribo_list) + ['NA'] *
                       tims_deepribo_gen,
                       'spectre': list(gene_spectre_list) + ['NA'] *
                       tims_spectre_gen,
                       'ribotish': list(gene_ribotish_list) + ['NA'] *
                       tims_ribotish_gen,
                       'reparation': list(gene_reparation_list) + ['NA'] *
                       tims_reparation_gen,
                       'irsom': list(gene_irsom_list) + ['NA']*tims_irsom_gen,
                       'spectre': list(gene_spectre_list) + ['NA']*
                       tims_spectre_gen,
                       'price': list(gene_price_list) + ['NA']*
                       tims_price_gen,
                       'ribotricer': list(gene_ribotricer_list) + ['NA']*
                       tims_ribotricer_gen,
                       'smorfer': list(gene_smorfer_list) + ['NA']*
                       tims_smorfer_gen}
    df_venn_genes = pd.DataFrame(venn_genes_dict)
    df_venn_genes.to_csv(path_or_buf=save_dir+'df_venn_genes.csv',
                         index=False, sep='\t')



#def venn_prepataions_bm03(gene_deepribo_list, gene_ribotish_list,
#                     gene_reparation_list, gene_irsom_list, save_dir):
#    max_no_of_genes = max([len(gene_deepribo_list),
#                        len(gene_ribotish_list),
#                        len(gene_reparation_list),
#                        len(gene_irsom_list)])
#    tims_deepribo_gen = max_no_of_genes - len(gene_deepribo_list)
#    tims_ribotish_gen = max_no_of_genes - len(gene_ribotish_list)
#    tims_reparation_gen = max_no_of_genes - len(gene_reparation_list)
#    tims_irsom_gen =  max_no_of_genes - len(gene_irsom_list)
#
#    venn_genes_dict = {'deepribo': list(gene_deepribo_list) + ['NA'] *
#                       tims_deepribo_gen,
#                       'ribotish': list(gene_ribotish_list) + ['NA'] *
#                       tims_ribotish_gen,
#                       'reparation': list(gene_reparation_list) + ['NA'] *
#                       tims_reparation_gen,
#                       'irsom': list(gene_irsom_list) + ['NA'] * tims_irsom_gen}
#    df_venn_genes = pd.DataFrame(venn_genes_dict)
#    df_venn_genes.to_csv(path_or_buf=save_dir+'df_venn_genes.csv',
#                         index=False, sep='\t')


def venn_prepataions_bm03(gene_deepribo_list, gene_ribotish_list,
                     gene_reparation_list, gene_irsom_list,
                     gene_price_list, gene_ribotricer_list, save_dir):
    max_no_of_genes = max([len(gene_deepribo_list),
                        len(gene_ribotish_list),
                        len(gene_reparation_list),
                        len(gene_irsom_list),
                        len(gene_price_list),
                        len(gene_ribotricer_list)])
    tims_deepribo_gen = max_no_of_genes - len(gene_deepribo_list)
    tims_ribotish_gen = max_no_of_genes - len(gene_ribotish_list)
    tims_reparation_gen = max_no_of_genes - len(gene_reparation_list)
    tims_irsom_gen =  max_no_of_genes - len(gene_irsom_list)
    tims_price_gen = max_no_of_genes - len(gene_price_list)
    tims_ribotricer_gen =max_no_of_genes - len(gene_ribotricer_list)

    venn_genes_dict = {'deepribo': list(gene_deepribo_list) + ['NA'] *
                       tims_deepribo_gen,
                       'ribotish': list(gene_ribotish_list) + ['NA'] *
                       tims_ribotish_gen,
                       'reparation': list(gene_reparation_list) + ['NA'] *
                       tims_reparation_gen,
                       'irsom': list(gene_irsom_list) + ['NA'] * tims_irsom_gen,
                       'price': list(gene_price_list) + ['NA']*
                       tims_price_gen,
                       'ribotricer': list(gene_ribotricer_list) + ['NA']*
                       tims_ribotricer_gen}
    df_venn_genes = pd.DataFrame(venn_genes_dict)
    df_venn_genes.to_csv(path_or_buf=save_dir+'df_venn_genes.csv',
                         index=False, sep='\t')





def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-p", "--reference_pos_data", action="store",
                        dest="reference_pos_path", required=True,
                        default="./data/smallORFs_labels_pos.gff",
                        help="path to the positive reference data.")
    parser.add_argument("-n", "--reference_neg_data", action="store",
                        dest="reference_neg_path", required=True,
                        default="../data/smallORFs_labels_neg.gff",
                        help="path to the negative reference data.")
    parser.add_argument("-t", "--tool_data", action="store", dest="tools_path",
                        required=True, default="./data/predictions.gff",
                        help="path to the tools data.")
    parser.add_argument("-o", "--save_path", action="store", dest="save_path",
                        required=True, default="./result_dfs/",
                        help="path to save data to.")
    parser.add_argument("-c", "--overlap_cutoff", action="store",
                        dest="overlap_cutoff", required=True,
                        default="0", type=float, help="overlap used to accept if prediction is covering labeled ORF")
    parser.add_argument("-s", "--flag_subopt", action="store",
                        dest="flag_subopt", required=True,
                        default="0", type=int, help="set to 1 to enable flag. If not 0.")
    parser.add_argument("-g", "--flag_no_gene", action="store",
                        dest="flag_no_gene", required=True,
                        default="0", type=int, help="set to 1 to enable flag. If not 0.")


    args = parser.parse_args()
    overlap_cutoff = args.overlap_cutoff
    flag_subopt = args.flag_subopt
    flag_no_gene = args.flag_no_gene

    save_dir = args.save_path

    ref_pos_file = args.reference_pos_path
    ref_neg_file = args.reference_neg_path
    tools_file = args.tools_path

    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)

    if not os.path.exists(args.tools_path):
        os.makedirs(args.tools_path)
    pos_overlap_gtf = (os.path.splitext(os.path.basename(args.reference_pos_path))[0] +
                   "_" + os.path.splitext(os.path.basename(args.tools_path))[0] +
                   "_overlap_pos.gtf")
    neg_overlap_gtf = (os.path.splitext(os.path.basename(args.reference_pos_path))[0] +
                   "_" + os.path.splitext(os.path.basename(args.tools_path))[0] +
                   "_overlap_neg.gtf")

    temp_overlap_gtf = (os.path.splitext(os.path.basename(args.reference_pos_path))[0] +
               "_" + os.path.splitext(os.path.basename(args.tools_path))[0] +
               "_overlap_temp.gtf")
    fp_overlap_gtf = (os.path.splitext(os.path.basename(args.reference_pos_path))[0] +
               "_" + os.path.splitext(os.path.basename(args.tools_path))[0] +
               "_overlap_fp.gtf")


    # bed tools intersect parameter:
    # -wo Write the original A and B entries plus the number of base pairs of overlap between the two features.
    # Only A features with overlap are reported. Restricted by -f (min overlap) and -r (min overlap couts for A and B).
    # -s Force “strandedness”. That is, only report hits in B that overlap A on the same strand. Default, strandless.
    pos_overlap_gtf = (save_dir + pos_overlap_gtf)
    neg_overlap_gtf = (save_dir + neg_overlap_gtf)
    df_pos_overlap = call_bedtools_overlap(ref_pos_file, tools_file, overlap_cutoff, pos_overlap_gtf, score_cutoff)
    df_neg_overlap = call_bedtools_overlap(ref_neg_file, tools_file, overlap_cutoff, neg_overlap_gtf, score_cutoff)

    temp_overlap_gtf = (save_dir + temp_overlap_gtf)
    fp_overlap_gtf = (save_dir + fp_overlap_gtf)

    # get fp that are not overlapping with any gene
    df_temp_overlap = call_bedtools_not_overlap(tools_file, ref_pos_file, overlap_cutoff, temp_overlap_gtf, score_cutoff)
    # all predictions which are not overlapping with any gene for the given cutoff
    df_fp_overlap = call_bedtools_not_overlap(temp_overlap_gtf, ref_neg_file, overlap_cutoff, fp_overlap_gtf, score_cutoff)
    #print(df_fp_overlap.info)

    df_filterd = df_fp_overlap[df_fp_overlap.start_pred > df_fp_overlap.stop_pred].copy()
    print('#########filterd df')
    print(len(df_filterd))
    assert len(df_filterd) == 0, "error: some prediction have a greater start than end postion"
    #
    # find statistics values
    stat_list_ribotish, gene_ribotish_list = get_stat_for_tool(ref_pos_file, tools_file, df_pos_overlap, ref_neg_file, df_neg_overlap, df_fp_overlap, flag_subopt, flag_no_gene, 'ribotish', save_dir)
    stat_list_reparation, gene_reparation_list = get_stat_for_tool(ref_pos_file, tools_file, df_pos_overlap, ref_neg_file, df_neg_overlap, df_fp_overlap, flag_subopt, flag_no_gene, 'reparation', save_dir)
    stat_list_deepribo, gene_deepribo_list = get_stat_for_tool(ref_pos_file, tools_file, df_pos_overlap, ref_neg_file, df_neg_overlap, df_fp_overlap, flag_subopt, flag_no_gene, 'deepribo', save_dir)
    stat_list_irsom, gene_irsom_list = get_stat_for_tool(ref_pos_file, tools_file, df_pos_overlap, ref_neg_file, df_neg_overlap, df_fp_overlap, flag_subopt, flag_no_gene, 'irsom', save_dir)
    stat_list_spectre, gene_spectre_list = get_stat_for_tool(ref_pos_file, tools_file, df_pos_overlap, ref_neg_file, df_neg_overlap, df_fp_overlap, flag_subopt, flag_no_gene, 'spectre', save_dir)
    stat_list_price, gene_price_list = get_stat_for_tool(ref_pos_file, tools_file, df_pos_overlap, ref_neg_file, df_neg_overlap, df_fp_overlap, flag_subopt, flag_no_gene, 'price', save_dir)
    stat_list_ribotricer, gene_ribotricer_list = get_stat_for_tool(ref_pos_file, tools_file, df_pos_overlap, ref_neg_file, df_neg_overlap, df_fp_overlap, flag_subopt, flag_no_gene, 'ribotricer', save_dir)
    stat_list_smorfer, gene_smorfer_list = get_stat_for_tool(ref_pos_file, tools_file, df_pos_overlap, ref_neg_file, df_neg_overlap, df_fp_overlap, flag_subopt, flag_no_gene, 'smorfer', save_dir)

    # print('\n!!!!!!!!!\ngene list deepribo:\n', gene_deepribo_list)

    venn_prepataions(gene_deepribo_list, gene_ribotish_list, gene_reparation_list, gene_irsom_list, gene_spectre_list, gene_price_list, gene_ribotricer_list, gene_smorfer_list, save_dir)

    #venn_prepataions_bm03(gene_deepribo_list, gene_ribotish_list, gene_reparation_list, gene_irsom_list, gene_price_list, gene_ribotricer_list, save_dir)
    #venn_prepataions_bm03(gene_deepribo_list, gene_ribotish_list, gene_reparation_list, gene_irsom_list, save_dir)

    if(flag_subopt == 1):
        stat_file = 'df_stat_suboptimals.csv'
    else:
        stat_file = 'df_stat.csv'

    stat_list_header = ['TP', 'TN', 'FP', 'FN', 'recall', 'specificity', 'FNR', 'precision',
                        'FDR', 'F1', 'accuracy', 'subopt_tp', 'subopt_fp', 'no_genes', 'tool']

    df_stat = pd.DataFrame([stat_list_reparation, stat_list_ribotish,
                            stat_list_deepribo, stat_list_irsom, stat_list_spectre,
                            stat_list_price, stat_list_ribotricer, stat_list_smorfer],
                            columns=stat_list_header)
    #df_stat = pd.DataFrame([stat_list_reparation, stat_list_ribotish,
    #                        stat_list_deepribo, stat_list_irsom,
    #                        stat_list_price, stat_list_ribotricer],
    #                        columns=stat_list_header)
    df_stat.to_csv(path_or_buf=save_dir + stat_file, index=False, sep='\t')


    with open(save_dir + "/table_latex.txt", "a") as myfile:
        ref_pos_list = ref_pos_file.split('/')
        myfile.write('Results of species %s, subset %s and overlap %f\n'%(ref_pos_list[8],ref_pos_list[9], overlap_cutoff))
        df_stat_round = df_stat.round(2)
        myfile.write(df_stat_round.to_latex(index=False))



if __name__ == '__main__':
    main()
