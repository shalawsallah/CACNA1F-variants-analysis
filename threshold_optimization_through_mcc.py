#! / usr / bin / env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import auc, matthews_corrcoef
from scipy.stats import *

np.set_printoptions(precision=2)
#########################################################################
'''Identifying the threshold at which the best MCC is obtained using 80% of the data (stored in a file), and using the 20% to test the thresholds, compared to using all the data to find the best threshold elsewhere'''

genes = ['gjb1','avpr2','pdha1','btk','ocrl','ndp','hprt1']
data_set = 'training' # options = training or test
# genes = ['scn1a']
###########################################
'''plotting parameters'''
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rc('axes', labelsize=MEDIUM_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)
mydensity = 1
myalpha = 0.8
#############################################
'''if plotting all genes' MCC scores for each tool on one single plot'''
# fig, ax = plt.subplots()
##############################################
vest4_thresholds_identified_for_genes = []
revel_thresholds_identified_for_genes = []
vest4_mcc_identified_for_genes = []
revel_mcc_identified_for_genes = []
for i in genes:

    gene = i
    print('\n*************************   ' + gene.upper())
    pathogenicity_threshold_numbers = 51
    pathogenicity_threshold = [i for i in np.linspace(0,1, pathogenicity_threshold_numbers)]
    # data = pd.read_excel('/Users/mdefsss2/' + gene + '/' + gene + '_analysis.xlsx')
    data = pd.read_csv('/Users/mdefsss2/x_linked_genes/'+gene+'/'+data_set+'_set.csv')
    genes_transcripts = pd.read_excel('/Users/mdefsss2/x_linked_genes/genes_list.xlsx')

    var_class = data['class']
    classes = [data['class_binary'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
    ensembl_transcripts = [data['Ensembl_transcriptid'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
    coordinates = [data['pos(1-based)'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
    revel_scores = [data['REVEL_score'].tolist()[m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
    vest4_scores = [data['VEST4_score'].tolist()[m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
    # pph2_scores = [data['Polyphen2_HVAR_score'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
    transcript = ''
    transcript_found = True  # the predictions for the variant in transcript og interest not always available
    for i in range(len(genes_transcripts)): #obtaining the trascript used
        if genes_transcripts['gene'][i] == gene:
            transcript += genes_transcripts['transcript'][i].strip()
    print('transcript retrieved: ', transcript)

    vest4_predictions = []
    # # pph2_predictions = []
    for i in range(len(ensembl_transcripts)):
        strings = str(ensembl_transcripts[i])
        if ';' in strings:
            num_of_predictions = ensembl_transcripts[i].split(';')
            vest4_pred = vest4_scores[i].split(';')
            # pph2_pred = pph2_scores[i].split(';')
            for k in range(len(num_of_predictions)):
                if num_of_predictions[k] == transcript:
                    try:
                        vest4_predictions.append(float(vest4_pred[k]))
                    except ValueError as e:
                        transcript_found = False
                        # print(e) # the transcript of interest wasnt available & the first transcript prediction is therefore chosen instead
                        vest4_predictions.append(float(vest4_pred[k+1]))# the next transcript is chosen but not considered for analysis
                elif transcript not in num_of_predictions:
                    print('Transcript of inerest MISSING!! at: ', coordinates[k])
        elif strings.strip() == transcript.strip():
            vest4_predictions.append(float(vest4_scores[i]))
    ########################################
    '''improving MCC scores for VEST4 in cases where there are multiple predictions for multiplr transcripts, where sorting by scores in EXCEL doesn't do the job. 
    Here, the scores for transcript of interest are written to a new file along with the classes so they can be sorted to identify a better pathogenicity threshold'''

    # np.savetxt('/Users/mdefsss2/dbNSFP4/'+gene+'_vest_scores.csv', [p for p in zip(classes, vest4_predictions)], delimiter=',', fmt='%s')
    ########################################
    revel_predicted_pathogenic = []
    revel_predicted_nonpathogenic = []
    vest4_predicted_pathogenic = []
    vest4_predicted_nonpathogenic = []
    # pph2_predicted_pathogenic = []
    # pph2_predicted_nonpathogenic = []

    '''Mann Whitney U test'''

    for i in range(len(classes)):
        if classes[i] == 1:
            revel_predicted_pathogenic.append(revel_scores[i])
            vest4_predicted_pathogenic.append(vest4_predictions[i])
        elif classes[i] == 0:
            revel_predicted_nonpathogenic.append(revel_scores[i])
            vest4_predicted_nonpathogenic.append(vest4_predictions[i])
    print('number of pathogenic vars: ', len(revel_predicted_pathogenic))
    print('number of non-pathogenic vars: ', len(revel_predicted_nonpathogenic))

    print('MannU test for REVEL: ', mannwhitneyu(revel_predicted_pathogenic, revel_predicted_nonpathogenic))
    print('MannU test for VEST4: ', mannwhitneyu(vest4_predicted_pathogenic, vest4_predicted_nonpathogenic))
     #########################################################################
    '''for a number of pathogenicity thresholds, calculate the MCC scores for both tools'''
    binary_class = []
    for i in range(len(classes)):
        if classes[i] == 0 or classes[i] == 1:
            binary_class.append(classes[i])

    revel_MCC_at_diff_thresholds = []
    revel_max_mcc = 0
    revel_max_mcc_threshold = 0
    prec_revel = []
    rec_revel = []

    for k in range(len(pathogenicity_threshold)):
        revel_pred = []
        for i in range(len(revel_scores)):
                if float(revel_scores[i]) >= pathogenicity_threshold[k]:
                    i = 1
                    revel_pred.append(i)
                elif float(revel_scores[i]) < pathogenicity_threshold[k]:
                    i = 0
                    revel_pred.append(i)
        counter_tp = 0
        counter_fn = 0
        counter_fp = 0
        counter_tn = 0
        for m in range(len(revel_predicted_pathogenic)):
            if revel_predicted_pathogenic[m] >= pathogenicity_threshold[k]:
                counter_tp += 1
            elif revel_predicted_pathogenic[m] < pathogenicity_threshold[k]:
                counter_fn += 1
        for l in range(len(revel_predicted_nonpathogenic)):
            if revel_predicted_nonpathogenic[l] >= pathogenicity_threshold[k]:
                counter_fp += 1
                # print('counter_fp: ', counter_fp)
            elif revel_predicted_nonpathogenic[l] < pathogenicity_threshold[k]:
                counter_tn += 1
        prec_revel.append(counter_tp / (counter_tp + counter_fp) if counter_tp + counter_fp != 0 and counter_tp != 0 else 0)
        rec_revel.append(counter_tp / (counter_tp + counter_fn))
        accuracy = (counter_tp + counter_tn) / (counter_tp + counter_tn + counter_fp + counter_fn)
        revel_mcc = matthews_corrcoef(binary_class, revel_pred)
        if pathogenicity_threshold[k] == 0.5: # suggested threshold
            # print('REVEL Precision @ ', pathogenicity_threshold[k], '=', prec_revel[k])
            # print('REVEL Recall @ ', pathogenicity_threshold[k], '=', rec_revel[k])
            print(gene.upper(), ' REVEL MCC: ', revel_mcc, ' @ ', pathogenicity_threshold[k])
        # print('matthews_corr_coef (REVEL) @ '+str(pathogenicity_threshold[k])+': ', revel_mcc)
        # print('Precision', prec_revel[k])
        # print('Recall', rec_revel[k])
        revel_MCC_at_diff_thresholds.append(revel_mcc)
        for d in range(len(revel_MCC_at_diff_thresholds) - 1):  # the optimum MCC score
            if d == 0:
                revel_max_mcc += revel_MCC_at_diff_thresholds[d]
            else:
                if revel_MCC_at_diff_thresholds[d + 1] > revel_max_mcc:
                    revel_max_mcc = revel_MCC_at_diff_thresholds[d + 1]
                    revel_max_mcc_threshold = pathogenicity_threshold[k]
        # if pathogenicity_threshold[k] == revel_max_mcc_threshold:
        #     print('REVEL Precision @ ', revel_max_mcc_threshold, '=', prec_revel[k])
        #     print('REVEL Recall @ ', revel_max_mcc_threshold, '=', rec_revel[k])
    print(gene.upper(), ' REVEL MCC: ', revel_max_mcc, ' @ ', revel_max_mcc_threshold)
    revel_thresholds_identified_for_genes.append(revel_max_mcc_threshold)
    revel_mcc_identified_for_genes.append(revel_max_mcc)
    vest4_MCC_at_diff_thresholds = []
    vest4_max_mcc = 0
    vest4_max_mcc_threshold = 0
    prec_vest4 = []
    rec_vest4 = []

    for k in range(len(pathogenicity_threshold)):
        vest4_pred = []
        for i in range(len(revel_scores)):
                if float(vest4_predictions[i]) >= pathogenicity_threshold[k]:
                    i = 1
                    vest4_pred.append(i)
                elif float(vest4_predictions[i]) < pathogenicity_threshold[k]:
                    i = 0
                    vest4_pred.append(i)
        counter_tp = 0
        counter_fn = 0
        counter_fp = 0
        counter_tn = 0
        for m in range(len(vest4_predicted_pathogenic)):
            if vest4_predicted_pathogenic[m] >= pathogenicity_threshold[k]:
                counter_tp += 1
            elif vest4_predicted_pathogenic[m] < pathogenicity_threshold[k]:
                counter_fn += 1
        for l in range(len(vest4_predicted_nonpathogenic)):
            if vest4_predicted_nonpathogenic[l] >= pathogenicity_threshold[k]:
                counter_fp += 1
                # print('counter_fp: ', counter_fp)
            elif vest4_predicted_nonpathogenic[l] < pathogenicity_threshold[k]:
                counter_tn += 1
        prec_vest4.append(
            counter_tp / (counter_tp + counter_fp) if counter_tp + counter_fp != 0 and counter_tp != 0 else 0)
        rec_vest4.append(counter_tp / (counter_tp + counter_fn))
        accuracy = (counter_tp + counter_tn) / (counter_tp + counter_tn + counter_fp + counter_fn)
        vest4_mcc = matthews_corrcoef(binary_class, vest4_pred)
        # vest4_mcc_list.append(vest4_mcc)
        # print('matthews_corr_coef (VEST4) @ '+str(pathogenicity_threshold[k])+': ', vest4_mcc)
        if pathogenicity_threshold[k] == 0.5:
            # print('VEST4 Precision @ ', pathogenicity_threshold[k], '=', prec_vest4[k])
            # print('VEST4 Recall @ ', pathogenicity_threshold[k], '=', rec_vest4[k])
            print(gene.upper(), ' VEST4 MCC: ', vest4_mcc, ' @ ', pathogenicity_threshold[k])
        vest4_MCC_at_diff_thresholds.append(vest4_mcc)
        for d in range(len(vest4_MCC_at_diff_thresholds) - 1):  # the optimum MCC score
            if d == 0:
                vest4_max_mcc += vest4_MCC_at_diff_thresholds[d]
            else:
                if vest4_MCC_at_diff_thresholds[d + 1] > vest4_max_mcc:
                    vest4_max_mcc = vest4_MCC_at_diff_thresholds[d + 1]
                    vest4_max_mcc_threshold = pathogenicity_threshold[k]
        # if pathogenicity_threshold[k] == vest4_max_mcc_threshold:
            # print('VEST4 Precision @ ', vest4_max_mcc_threshold, '=', prec_vest4[k])
            # print('VEST4 Recall @ ', vest4_max_mcc_threshold, '=', rec_vest4[k])
    if transcript_found == False:
        print('VEST4 prediction NOT found for Transcript of interest!!')
    print(gene.upper(), ' VEST4 MCC: ', vest4_max_mcc, ' @ ', vest4_max_mcc_threshold)
    vest4_thresholds_identified_for_genes.append(vest4_max_mcc_threshold)
    if transcript_found == True:
        vest4_mcc_identified_for_genes.append(vest4_max_mcc)
    elif transcript_found == False:
        vest4_mcc_identified_for_genes.append('-')

    #######################################################
    '''comment out if plotting all genes' MCC on one graph for a tool'''
    # fig = plt.figure()

#######################################################
    '''plotting VEST4 & REVEL MCC scores for each gene for different thresholds on one plot'''
    # ax = fig.add_subplot(1,1,1)
    # ax.plot(pathogenicity_threshold, revel_MCC_at_diff_thresholds, label='REVEL')
    # # ax.plot(pathogenicity_threshold, vest4_MCC_at_diff_thresholds, label='VEST4')
    # ax.set_xlabel('Thresholds')
    # ax.set_ylabel('MCC')
    # ax.set_title(gene)
    # plt.legend()
    # plt.show()
#######################################################
    '''plotting each tool's MCC scores for all genes together on one plot'''
#     # plt.plot(pathogenicity_threshold, revel_MCC_at_diff_thresholds, label='')
#     ax.plot(pathogenicity_threshold, vest4_MCC_at_diff_thresholds, label='')
#     plt.xlabel('Thresholds')
#     plt.ylabel('MCC')
#     plt.title('VEST4')
#     plt.legend()
# plt.show()
#######################################################
    '''showing the best MCC scores for each gene in a scatter plot'''

plt.close('all')
print('####################################################')
print('REVEL Thresholds & MCCc for each gene:\n')
for i in range(len(revel_thresholds_identified_for_genes)):
    print(revel_thresholds_identified_for_genes[i], round(revel_mcc_identified_for_genes[i], 2))

print('####################################################')
print('VEST4 Thresholds & MCCs for each gene:\n')
for i in range(len(vest4_thresholds_identified_for_genes)):
    if vest4_mcc_identified_for_genes[i] == float or vest4_mcc_identified_for_genes[i] == int:
        print(vest4_thresholds_identified_for_genes[i], round(vest4_mcc_identified_for_genes[i], 2))
    else:
        print(vest4_thresholds_identified_for_genes[i], vest4_mcc_identified_for_genes[i])