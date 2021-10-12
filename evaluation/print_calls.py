#!/usr/bin/env python
import pandas as pd
#import mathe
import matplotlib as mpl
import os
import argparse
import csv
import collections
from operator import itemgetter
import sys
from datetime import datetime
import pickle
import glob


def main():
    # store commandline args
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-i", "--input_path", action="store",
                        dest="input_path", required=True,
                        default="./datasets/",
                        help="path to the input gffs.")
    parser.add_argument("-o", "--output_path", action="store",
                        dest="output_path", required=True,
                        default="./results/",
                        help="path to the negative reference data.")
    parser.add_argument("-s", "--python_script_dir", action="store",
                        dest="python_script_dir", required=True,
                        default="/home/teresa/Dokumente/benchmark_ribo_seq/RiboReport/evaluation/",
                        help="path to your '/RiboReport/evaluation/' forder where all evaluation scripts are stored.")




    args = parser.parse_args()

    input_dir = args.input_path
    output_dir = args.output_path
    python_script_dir = args.python_script_dir


    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    date = datetime.today().strftime('%Y_%m_%d')
    result_dir = output_dir + '/' + date + '_results'
    # os.mkdir(result_dir)
    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
    else:
        print('\n+++++++\noutput dir \'%s\' already exists\n+++++++++\n'%(result_dir))

    shell_script =  result_dir + '/evaluation_calls.sh'

    python_script_dir = '/home/teresa/Dokumente/benchmark_ribo_seq/RiboReport/evaluation/'


    f= open(shell_script,"w+")
    f.write("#!/usr/bin/env bash\n\ntrap ctrl_c INT\n\nfunction ctrl_c() {\necho \"** Trapped CTRL-C\"\nexit\n}\n\n")
    #f.write('test')
    # what to investigat test:
    #dataset_list = ['bm_01', 'bm_03', 'bm_06', 'bm_12']
    #dataset_list = ['bm_14']
    #experiment_list = ['smallORFs_labels']

    #experiment_list = ['CDS_labels', 'smallORFs_labels', 'operons_intersect_labels', 'operons_complement_labels']
    #overlap_list = ['0.01', '0.5']

    #bar_plot_call_dict = {}
    # full data:
    #dataset_list = ['bm_01', 'bm_03', 'bm_06', 'bm_12']
    dataset_list = ['bm_01']
    #dataset_list = ['bm_06', 'bm_12']
    #dataset_list = ['bm_03']

    experiment_list = ['CDS_labels', 'smallORFs_labels', 'operons_intersect_labels', 'operons_complement_labels']
    overlap_list =  ['0.01', '0.9',  '0.7']
    for dataset in dataset_list:
        bar_plot_call_dict = {}
        f.write('\n')
        save_dataset_dir = result_dir + '/' + dataset
        open_dataset_dir = input_dir + '/' + dataset
        #print(dataset)



        for experiment in experiment_list:
            #print(experiment)
            f.write('#%s with %s\n' % (dataset, experiment))
            save_experiment_dir = save_dataset_dir + '/' + experiment
            open_experiment_pos_file = open_dataset_dir + '/' + experiment + '_pos.gff'
            open_experiment_neg_file = open_dataset_dir + '/' + experiment + '_neg.gff'
            prediction_file = open_dataset_dir + '/predictions.gff'

            bar_key = save_experiment_dir + '/'
            bar_plot_call_dict[bar_key] = []
            for overlap in overlap_list:
                #print(overlap)
                bar_plot_call_dict[bar_key].append(overlap)
                save_overlap_dir = save_experiment_dir + '/' + overlap + '/'
                temp_dir = result_dir + '/' + dataset + '/' + experiment + '/' + overlap + '/'
                #os.mkdir(temp_dir)
                if not os.path.exists(temp_dir):
                    os.makedirs(temp_dir)
                if experiment == 'CDS_labels':
                    #call = ('python3 ' + python_script_dir + 'statistics_pos_neg.py -p ' + open_experiment_pos_file + ' -n '+
                            #open_experiment_neg_file + ' -t ' + prediction_file + ' -o ' +
                            #save_overlap_dir + ' -c ' + overlap + ' -s 0 -g 1')
                     call = ('python3 ' + python_script_dir + 'statistics_pos_neg.py -p ' + open_experiment_pos_file + ' -n '+
                            open_experiment_neg_file + ' -t ' + prediction_file + ' -o ' +
                            save_overlap_dir + ' -c ' + overlap + ' -s 0 -g 0')
                else:
                    call = ('python3 ' + python_script_dir + 'statistics_pos_neg.py -p ' + open_experiment_pos_file + ' -n '+
                            open_experiment_neg_file + ' -t ' + prediction_file + ' -o ' +
                            save_overlap_dir + ' -c ' + overlap + ' -s 0 -g 0')

                f.write('%s\n'%(call))
                call_subopt = ('python3 ' + python_script_dir + 'statistics_pos_neg.py -p ' + open_experiment_pos_file + ' -n '+
                        open_experiment_neg_file + ' -t ' + prediction_file + ' -o ' +
                        save_overlap_dir + ' -c ' + overlap + ' -s 1 -g 0')
                #f.write('%s\n'%(call_subopt))
                #print(call)
        #f.write('\n#barplos\n')



        for bar_path in bar_plot_call_dict:
            call_temp = ('python3  %s/plot_barplots_3_inputs.py '%(python_script_dir))
            for c, overlap in enumerate(bar_plot_call_dict[bar_path]):
                call_temp = ('%s-i%d=%s%s/df_stat.csv '%(call_temp, (c+1), bar_path, overlap))

            call_bar_plots = ('%s -o %s '%(call_temp,  bar_path))
            # print(call_bar_plots)
            #f.write('%s\n'%(call_bar_plots))

        # generate calls for plots:

    plot_dir = result_dir + '/plots/'
    #print(plot_dir)

    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    experiment_dict = {}
    for experiment in experiment_list:
        #print(experiment)
        experiment_dict[experiment]=[]
        for dataset in dataset_list:
            if dataset == 'bm_01':
                lable = 'EC'
            elif dataset == 'bm_03':
                lable = 'LM'
            elif dataset == 'bm_06':
                lable = 'PA'
            elif dataset == 'bm_12':
                lable = 'ST'
            elif dataset == 'bm_14':
                lable = 'HV'
            for overlap in overlap_list:
                #print(lable)
                stat_path = result_dir + '/' + dataset + '/' + experiment + '/' + overlap+ '/df_stat.csv'
                hue_symbel = lable + '_' + overlap
                experiment_dict[experiment].append((stat_path,hue_symbel, dataset, overlap))

                # Venn diagram call
                venn_path = result_dir + '/' + dataset + '/' + experiment + '/' + overlap+ '/'
                venn_input = venn_path + 'df_venn_genes.csv'


                call_script_venn = ('python3  %s/venn_diagram.py '%(python_script_dir))
                call_venn = ('%s -i %s -o %s -n %s -c %s'%(call_script_venn, venn_input, plot_dir, (lable + '_' + experiment) , overlap))
                #print(call_venn)
                #f.write('\n#####################\n%s\n'%(call_venn))

                # construct call for ROC cove
                if experiment == 'CDS_labels':
                    ###### Change!!!!!!
                    call_script_roc = ('python3  %s/roc_plotting.py '%(python_script_dir))
                    score_list_dir = result_dir + '/' + dataset + '/CDS_labels/' + overlap + '/'
                    call_roc = ('%s -i %s -e %s -o %s -p %s'%(call_script_roc, score_list_dir, lable, plot_dir, overlap))
                    #print(call_roc)
                    f.write('\n#######################\n%s\n'%(call_roc))
                    ######Change!!!!!!

                call_script_prc = ('python3  %s/prc_plotting.py '%(python_script_dir))
                score_list_dir = result_dir + '/' + dataset + '/' + experiment + '/' + overlap + '/'
                call_prc = ('%s -i %s -s %s -o %s -p %s -e %s'%(call_script_prc, score_list_dir, lable, plot_dir, overlap, experiment))
                #print(call_prc)
                f.write('\n#######################\n%s\n'%(call_prc))

    save_experiment_dict = plot_dir +'experiment_dict'

    with open(save_experiment_dict, 'wb') as handle:
        pickle.dump(experiment_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    call_script = ('python3  %s/experiment_barplots.py '%(python_script_dir))
    call_exp_barplots = ('%s -e %s -o %s '%(call_script, save_experiment_dict, plot_dir))
    #print(call_exp_barplots)
    # f.write('\n#######################\n%s\n'%(call_exp_barplots))




    f.close()
    print('to run the evaluation pipline call:\n%s'%(result_dir + '/evaluation_calls.sh'))
    os.chmod(shell_script, 0o777)

if __name__ == '__main__':
    main()
