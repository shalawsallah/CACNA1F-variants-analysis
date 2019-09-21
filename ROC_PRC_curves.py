#! /usr/bin/env python3
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics
# from sklearn.metrics import auc, precision_recall_curve, average_precision_score
from scipy.stats import *
from matplotlib import rc
# import seaborn as sns
import math
# from sklearn.metrics import auc, matthews_corrcoef

'''A precision of 0.5 (on Y axis) at a point (on X axis) on the PRC means 1 correct positive for 1 false positive.
So the expected PR curve for a random classifier is just a rectangle with side lengths "proportion of true positives" x 1. For example, if your dataset contains 10% positive cases and 90% negative cases, the expected auPR under chance is 0.1.

To obtain the AUC of all the folds in the cross validation as well as their average curve : https://scikit-learn.org/stable/auto_examples/model_selection/plot_roc_crossval.html#sphx-glr-auto-examples-model-selection-plot-roc-crossval-py'''

all_data = pd.read_csv('/Users/mdefsss2/cacna1f/WEKA_old_results/combined_train_test_data.csv')
model_data = pd.read_csv('/Users/mdefsss2/cacna1f/E_J_human_genetics/weka_results/prediction_AUCs.csv')
model_predictions = pd.read_csv('/Users/mdefsss2/cacna1f/E_J_human_genetics/weka_results/logistic_predictions.csv')
################################################

roc_slicing_numbers = 1001
print('number of variants: ', len(all_data))
# data_dropna = ml_data.dropna()
# print(len(data_dropna))
thresholds = model_data['Threshold']
model_tpr = model_data["'True Positive Rate'"]
model_fpr = model_data["'False Positive Rate'"]
model_prec = model_data["Precision"]
model_rec = model_data["Recall"]
pph2_predicted_pathogenic = []
pph2_predicted_nonpathogenic = []
cadd_predicted_pathogenic = []
cadd_predicted_nonpathogenic = []
condel_predicted_pathogenic = []
condel_predicted_nonpathogenic = []
ma_predicted_pathogenic = []
ma_predicted_nonpathogenic = []
fathmm_predicted_pathogenic = []
fathmm_predicted_nonpathogenic = []
provean_predicted_pathogenic = []
provean_predicted_nonpathogenic = []
sift_predicted_pathogenic = []
sift_predicted_nonpathogenic = []
variants = all_data['variants']
target_class = all_data['class'].tolist()
target_class_binary = all_data['binary_class'].tolist()
# provean_scores = all_data['Provean SCORE'].tolist()
pph2_scores = all_data[' pph2_prob'].tolist()
cadd_scores = all_data['PHRED'].tolist()
condel_scores = all_data['CONDEL'].tolist()
sift_scores = all_data['1-SIFT Score'].tolist()
# ma_scores = all_data['Provean SCORE'].tolist()
# fathmm_scores = all_data['Provean SCORE'].tolist()
# exit(2)
for i in range(len(variants)):
    if target_class[i] == 'pathogenic':
        # provean_predicted_pathogenic.append(provean_scores[i])
        pph2_predicted_pathogenic.append(pph2_scores[i])
        cadd_predicted_pathogenic.append(cadd_scores[i])
        condel_predicted_pathogenic.append(condel_scores[i])
        sift_predicted_pathogenic.append(sift_scores[i])
        # ma_predicted_pathogenic.append(ma_scores[i])
        # fathmm_predicted_pathogenic.append(fathmm_scores[i])
    elif target_class[i] == 'benign':
        # provean_predicted_nonpathogenic.append(provean_scores[i])
        pph2_predicted_nonpathogenic.append(pph2_scores[i])
        cadd_predicted_nonpathogenic.append(cadd_scores[i])
        condel_predicted_nonpathogenic.append(condel_scores[i])
        sift_predicted_nonpathogenic.append(sift_scores[i])
        # ma_predicted_nonpathogenic.append(ma_scores[i])
        # fathmm_predicted_nonpathogenic.append(fathmm_scores[i])

# model_scores = model_predictions['prediction_score']
# model_class = model_predictions['class']
# model_predicted_pathogenic = []
# model_predicted_nonpathogenic = []
# for i in range(len(model_scores)):
#     if model_class[i] == 'pathogenic':
#         model_predicted_pathogenic.append(model_scores[i])
#     if model_class[i] == 'benign':
#         model_predicted_nonpathogenic.append(model_scores[i])

mynormed = 1
myalpha = 0.3
mybins = 10
# font = {'family' : 'normal',
#         'weight' : 'bold',
#         'size'   : 14}
# rc('font', **font)
# plt.rcParams.update({'font.size': 14})
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

# plt.hist(model_predicted_pathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Pathogenic')
# plt.hist(model_predicted_nonpathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Non-pathogenic')
# plt.xlabel('Pathogenicity scores')
# plt.ylabel('Frequency')
# plt.title('Model')
# plt.legend(loc='best')
# plt.show()

plt.hist(pph2_predicted_pathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Pathogenic')
plt.hist(pph2_predicted_nonpathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Non-pathogenic')
plt.xlabel('Pathogenicity scores')
plt.ylabel('Frequency')
plt.title('PolyPhen2')
plt.legend(loc='best')
plt.show()

plt.hist(cadd_predicted_pathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Pathogenic')
plt.hist(cadd_predicted_nonpathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Non-pathogenic')
plt.xlabel('Pathogenicity scores')
plt.ylabel('Frequency')
plt.title('CADD')
plt.legend(loc='best')
plt.show()

plt.hist(condel_predicted_pathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Pathogenic')
plt.hist(condel_predicted_nonpathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Non-pathogenic')
plt.xlabel('Pathogenicity scores')
plt.ylabel('Frequency')
plt.title('CONDEL')
plt.legend(loc='best')
plt.show()

# plt.hist(ma_predicted_pathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Pathogenic')
# plt.hist(ma_predicted_nonpathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Non-pathogenic')
# plt.xlabel('Pathogenicity scores')
# plt.ylabel('Frequency')
# plt.title('MA')
# plt.legend(loc='best')
# plt.show()
#
# plt.hist(fathmm_predicted_pathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Pathogenic')
# plt.hist(fathmm_predicted_nonpathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Non-pathogenic')
# plt.xlabel('Pathogenicity scores')
# plt.ylabel('Frequency')
# plt.title('FATHMM')
# plt.legend(loc='best')
# plt.show()

plt.hist(sift_predicted_pathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Pathogenic')
plt.hist(sift_predicted_nonpathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Non-pathogenic')
plt.xlabel('Pathogenicity scores')
plt.ylabel('Frequency')
plt.title('SIFT')
plt.legend(loc='best')
plt.show()

# plt.hist(provean_predicted_pathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Pathogenic')
# plt.hist(provean_predicted_nonpathogenic, bins=mybins, alpha=myalpha, density=mynormed, label='Non-pathogenic')
# plt.xlabel('Pathogenicity scores')
# plt.ylabel('Frequency')
# plt.title('Provean')
# plt.legend(loc='best')
# plt.show()

##############################################
'''ROCs'''

'''WEKA model'''

model_roc_auc = metrics.auc(model_fpr, model_tpr)
print('model ROC: ', model_roc_auc)
####################################################
'''plot all ROCs using python-calculated fpr & tpr & ROC AUC'''


plt.plot(model_fpr, model_tpr, label='CACNA1F-vp (area = %0.2f)' % (model_roc_auc))

fpr, tpr, thresholds = metrics.roc_curve(target_class_binary, sift_scores, pos_label=1)
# print('SIFT FPR: ', fpr)
# print('SIFT TPR: ', tpr)
roc_auc = metrics.auc(fpr, tpr)
plt.plot(fpr, tpr, label='SIFT (area = %0.2f)' % (roc_auc))

# fpr, tpr, thresholds = metrics.roc_curve(target_class_binary, provean_scores, pos_label=0)
# roc_auc = metrics.auc(fpr, tpr)
# plt.plot(fpr, tpr, label='Provean (area = %0.2f)' % (roc_auc))

# fpr, tpr, thresholds = metrics.roc_curve(target_class_binary, fathmm_scores, pos_label=0)
# roc_auc = metrics.auc(fpr, tpr)
# plt.plot(fpr, tpr, label='FATHMM (area = %0.2f)' % (roc_auc))

# fpr, tpr, thresholds = metrics.roc_curve(target_class_binary, ma_scores, pos_label=1)
# roc_auc = metrics.auc(fpr, tpr)
# plt.plot(fpr, tpr, label='MA (area = %0.2f)' % (roc_auc))

fpr, tpr, thresholds = metrics.roc_curve(target_class_binary, pph2_scores, pos_label=1)
roc_auc = metrics.auc(fpr, tpr)
plt.plot(fpr, tpr, label='PolyPhen2 (area = %0.2f)' % (roc_auc))

fpr, tpr, thresholds = metrics.roc_curve(target_class_binary, cadd_scores, pos_label=1)
# print('CADD FPR: ', fpr)
# print('CADD TPR: ', tpr)
roc_auc = metrics.auc(fpr, tpr)
plt.plot(fpr, tpr, label='CADD (area = %0.2f)' % (roc_auc))

fpr, tpr, thresholds = metrics.roc_curve(target_class_binary, condel_scores, pos_label=1)
roc_auc = metrics.auc(fpr, tpr)
plt.plot(fpr, tpr, label='CONDEL (area = %0.2f)' % (roc_auc))

plt.legend(loc='lower right')
plt.plot([0, 1], [0, 1], 'r--')
plt.xlim([0.0, 1.02])
plt.ylim([0.0, 1.02])
plt.ylabel('True Positive Rate (sensitivity)')
plt.xlabel('False Positive Rate (1-specificity)')
plt.show()

sift_roc = metrics.roc_auc_score(target_class_binary, sift_scores)
print('sift_roc score: ', 1-sift_roc)
# provean_roc = metrics.roc_auc_score(target_class_binary, provean_scores)
# print('provean_roc score: ', 1-provean_roc)
# fathmm_roc = metrics.roc_auc_score(target_class_binary, fathmm_scores)
# print('fathmm_roc score: ', 1-fathmm_roc)
# ma_roc = metrics.roc_auc_score(target_class_binary, ma_scores)
# print('MA_roc score: ', ma_roc)
pph2_roc = metrics.roc_auc_score(target_class_binary, pph2_scores)
print('PPh2_roc score: ', pph2_roc)
cadd_roc = metrics.roc_auc_score(target_class_binary, cadd_scores)
print('CADD_roc score: ', cadd_roc)
condel_roc = metrics.roc_auc_score(target_class_binary, condel_scores)
print('CONDEL_roc score: ', condel_roc)

# # create the axis of thresholds (scores)
# plt.plot(fpr, tpr)  # add the ROC curve to it
# ax2 = plt.gca().twinx()
# ax2.plot(fpr, thresholds, markeredgecolor='r', linestyle='dashed', color='r')
# ax2.set_ylabel('Threshold', color='r')
# ax2.set_ylim([thresholds[0], thresholds[-1]])
# ax2.set_xlim([fpr[0], fpr[-1]])
# plt.show()
##############################
'''PRC of all tools using python-calculated precision & recall & PRC AUC'''

model_prc_auc = metrics.auc(model_rec, model_prec)
print('model PRC: ', model_prc_auc)
plt.plot(model_rec, model_prec, label='CACNA1F-vp (area = %0.2f)' % (model_prc_auc))

precision, recall, thresholds = metrics.precision_recall_curve(target_class_binary, sift_scores, pos_label=1)
prc_auc = metrics.auc(recall,precision)
plt.plot(recall, precision, label='SIFT (area = %0.2f)' % (prc_auc))

# precision, recall, thresholds = metrics.precision_recall_curve(target_class_binary, provean_scores, pos_label=1)
# prc_auc = metrics.auc(recall,precision)
# plt.plot(recall, precision, label='Provean (area = %0.2f)' % (prc_auc))

# precision, recall, thresholds = metrics.precision_recall_curve(target_class_binary, fathmm_scores, pos_label=1)
# prc_auc = metrics.auc(recall,precision)
# plt.plot(recall, precision, label='FATHMM (area = %0.2f)' % (prc_auc))

# precision, recall, thresholds = metrics.precision_recall_curve(target_class_binary, ma_scores, pos_label=1)
# prc_auc = metrics.auc(recall,precision)
# plt.plot(recall, precision, label='MA (area = %0.2f)' % (prc_auc))

precision, recall, thresholds = metrics.precision_recall_curve(target_class_binary, pph2_scores, pos_label=1)
prc_auc = metrics.auc(recall,precision)
plt.plot(recall, precision, label='PolyPhen2 (area = %0.2f)' % (prc_auc))

precision, recall, thresholds = metrics.precision_recall_curve(target_class_binary, cadd_scores, pos_label=1)
# print('CADD precision: ', precision)
# print('CADD recall: ', recall)
prc_auc = metrics.auc(recall,precision)
plt.plot(recall, precision, label='CADD (area = %0.2f)' % (prc_auc))

precision, recall, thresholds = metrics.precision_recall_curve(target_class_binary, condel_scores, pos_label=1)
prc_auc = metrics.auc(recall,precision)
plt.plot(recall, precision, label='CONDEL (area = %0.2f)' % (prc_auc))

plt.legend(loc='best')
# plt.plot([0, 1], [0, 1], 'r--')
plt.xlim([0.0, 1.02])
plt.ylim([0.0, 1.02])
plt.ylabel('Precision')
plt.xlabel('Recall')
plt.show()

sift_precision = metrics.average_precision_score(np.array(target_class_binary), np.array(sift_scores))
print('SIFT PRC score: ', sift_precision)
# provean_precision = metrics.average_precision_score(np.array(target_class_binary), np.array(provean_scores))
# print('Provean PRC score: ', provean_precision)
# fathmm_precision = metrics.average_precision_score(np.array(target_class_binary), np.array(fathmm_scores))
# print('FATHMM PRC score: ', fathmm_precision)
# ma_precision = metrics.average_precision_score(np.array(target_class_binary), np.array(ma_scores))
# print('ma PRC score: ', ma_precision)
pph2_precision = metrics.average_precision_score(np.array(target_class_binary), np.array(pph2_scores))
print('PPh2 PRC score: ', pph2_precision)
cadd_precision = metrics.average_precision_score(np.array(target_class_binary), np.array(cadd_scores))
print('CADD PRC score: ', cadd_precision)
condel_precision = metrics.average_precision_score(np.array(target_class_binary), np.array(condel_scores))
print('CONDEL PRC score: ', condel_precision)

##################################################################
##################################################################
'''Provean'''

# tpr_provean = []
# fpr_provean = []
# prec_provean = []
# rec_provean = []
# thresholds = [i for i in np.linspace(min(provean_scores), max(provean_scores), roc_slicing_numbers)]  # (0, 1.25, 41) when adding the weights & -137, 25, 41) when multiplying, & (-1.3, 1.1, 41) when using LMT
# # print('thresholds: ', thresholds)
# for i in range(len(thresholds)):
#     counter_tp = 0
#     counter_fn = 0
#     counter_fp = 0
#     counter_tn = 0
#     for k in range(len(provean_predicted_pathogenic)):
#         if provean_predicted_pathogenic[k] >= thresholds[i]:
#             counter_fn += 1
#         elif provean_predicted_pathogenic[k] < thresholds[i]:
#             counter_tp += 1
#     # print('threshold: ', thresholds[i])
#     # print('counter_fn: ', counter_fn)
#     # print('counter_tp: ', counter_tp)
#     # print('tpr_provean: ', counter_tp / (counter_tp + counter_fn))
#     tpr_provean.append(counter_tp / (counter_tp + counter_fn) if counter_tp + counter_fn != 0 and counter_tp != 0 else 0)
#     # if thresholds[i] == 0.5:
#     #     print('tpr_provean at 0.5: ', tpr1[i])
#
#     for l in range(len(provean_predicted_nonpathogenic)):
#         if provean_predicted_nonpathogenic[l] >= thresholds[i]:
#             counter_tn += 1
#         elif provean_predicted_nonpathogenic[l] < thresholds[i]:
#             counter_fp += 1
#     # print('threshold: ', thresholds[i])
#     # print('counter_fp: ', counter_fp)
#     # print('counter_tn: ', counter_tn)
#     # print('fpr_provean:', counter_fp / (counter_fp + counter_tn))
#     fpr_provean.append(counter_fp / (counter_fp + counter_tn) if counter_fp + counter_tn != 0 and counter_fp != 0 else 0)
#     # if thresholds[i] == 0.5:
#     #     print('fpr_provean at 0.5: ', fpr1[i])
#     prec_provean.append(counter_tp / (counter_tp + counter_fp) if counter_tp + counter_fp != 0 and counter_tp != 0 else 0)
#     rec_provean.append(counter_tp / (counter_tp + counter_fn))
# # print('tpr_provean: ', tpr_provean)
# # print('fpr_provean: ', fpr_provean)
# roc_auc_provean = metrics.auc(fpr_provean, tpr_provean)
# prc_auc_provean = metrics.auc(rec_provean, prec_provean)
# print('Provean ROC: ', roc_auc_provean)
# print('Provean PRC: ', prc_auc_provean)
#
# accuracy = (counter_tp + counter_tn) / (counter_tp + counter_tn + counter_fp + counter_fn)
# print('Provean accuracy: ', accuracy)
# print('counter_fp: ', counter_fp)
# print('counter_tn: ', counter_tn)
# print('counter_fn: ', counter_fn)
# print('counter_tp: ', counter_tp)
# mcc = ((counter_tp*counter_tn) - (counter_fp*counter_fn)) / math.sqrt( (counter_tp + counter_fp)*(counter_tp + counter_fn)*(counter_tn + counter_fp)*(counter_tn + counter_fn) )
# print('MCC Provean: ', mcc)

#############################################################
'''SIFT ROC'''
tpr_sift = []
fpr_sift = []
prec_sift = []
rec_sift = []
thresholds = [i for i in np.linspace(min(sift_scores), max(sift_scores), roc_slicing_numbers)]
for i in range(len(thresholds)):
    counter_tp = 0
    counter_fn = 0
    counter_fp = 0
    counter_tn = 0
    optimum_threshold_tp = []
    optimum_threshold_tn = []
    optimum_threshold_fp = []
    optimum_threshold_fn = []
    for k in range(len(sift_predicted_pathogenic)):
        if sift_predicted_pathogenic[k] >= thresholds[i]:
            counter_fn += 1
        elif sift_predicted_pathogenic[k] < thresholds[i]:
            counter_tp += 1
    if float(f'{thresholds[i]:.2f}') == 0.05:
        print('SIFT counter_tp', counter_tp)
        print('SIFT counter_fn', counter_fn)
        print('SIFT TPR: ', counter_tp / (counter_tp + counter_fn) if counter_tp + counter_fn != 0 and counter_tp != 0 else 0)
        optimum_threshold_tp.append(counter_tp)
        optimum_threshold_fn.append(counter_fn)
    # print('threshold: ', thresholds[i])
    # print('counter_fn: ', counter_fn)
    # print('counter_tp: ', counter_tp)
    # print('tpr SIFT: ', counter_tp / (counter_tp + counter_fn))
    tpr_sift.append(counter_tp / (counter_tp + counter_fn) if counter_tp + counter_fn != 0 and counter_tp != 0 else 0)
    for l in range(len(sift_predicted_nonpathogenic)):
        if sift_predicted_nonpathogenic[l] >= thresholds[i]:
            counter_tn += 1
            # print('counter_fp: ', counter_fp)
        elif sift_predicted_nonpathogenic[l] < thresholds[i]:
            counter_fp += 1
    if float(f'{thresholds[i]:.2f}') == 0.05:
        print('SIFT counter_tn', counter_tn)
        print('SIFT counter_fp', counter_fp)
        print('SIFT FPR: ', counter_fp / (counter_fp + counter_tn) if counter_fp + counter_tn != 0 and counter_fp != 0 else 0)
        optimum_threshold_tn.append(counter_tn)
        optimum_threshold_fp.append(counter_fp)
            # print('counter_tn: ', counter_tn)
    # print('threshold: ', thresholds[i])
    # print('counter_fp: ', counter_fp)
    # print('counter_tn: ', counter_tn)
    # print('fpr:', counter_fp / (counter_fp + counter_tn))
    fpr_sift.append(counter_fp / (counter_fp + counter_tn) if counter_fp + counter_tn != 0 and counter_fp != 0 else 0)
    prec_sift.append(counter_tp / (counter_tp + counter_fp) if counter_tp + counter_fp != 0 and counter_tp != 0 else 0)
    rec_sift.append(counter_tp / (counter_tp + counter_fn))
roc_auc_sift = metrics.auc(fpr_sift, tpr_sift)
prc_auc_sift = metrics.auc(rec_sift, prec_sift)

print('SIFT ROC: ', roc_auc_sift)
print('SIFT PRC: ', prc_auc_sift)

# accuracy = (counter_tp + counter_tn) / (counter_tp + counter_tn + counter_fp + counter_fn)
# mcc = (counter_tp*counter_tn - counter_fp*counter_fn) / math.sqrt( (counter_tp + counter_fp)*(counter_tp + counter_fn)*(counter_tn + counter_fp)*(counter_tn + counter_fn))

#############################################################
'''FATHMM ROC'''
# tpr_fathmm = []
# fpr_fathmm = []
# prec_fathmm = []
# rec_fathmm = []
# thresholds = [i for i in np.linspace(min(fathmm_scores), max(fathmm_scores), roc_slicing_numbers)]
# for i in range(len(thresholds)):
#     counter_tp = 0
#     counter_fn = 0
#     counter_fp = 0
#     counter_tn = 0
#     for k in range(len(fathmm_predicted_pathogenic)):
#         if fathmm_predicted_pathogenic[k] >= thresholds[i]:
#             counter_fn += 1
#         elif fathmm_predicted_pathogenic[k] < thresholds[i]:
#             counter_tp += 1
#     # print('threshold: ', thresholds[i])
#     # print('counter_fn: ', counter_fn)
#     # print('counter_tp: ', counter_tp)
#     # print('tpr: ', counter_tp / (counter_tp + counter_fn))
#     tpr_fathmm.append(counter_tp / (counter_tp + counter_fn) if counter_tp + counter_fn != 0 and counter_tp != 0 else 0)
#
#     for l in range(len(fathmm_predicted_nonpathogenic)):
#         if fathmm_predicted_nonpathogenic[l] >= thresholds[i]:
#             counter_tn += 1
#             # print('counter_fp: ', counter_fp)
#         elif fathmm_predicted_nonpathogenic[l] < thresholds[i]:
#             counter_fp += 1
#             # print('counter_tn: ', counter_tn)
#     # print('threshold: ', thresholds[i])
#     # print('counter_fp: ', counter_fp)
#     # print('counter_tn: ', counter_tn)
#     # print('fpr:', counter_fp / (counter_fp + counter_tn))
#     fpr_fathmm.append(counter_fp / (counter_fp + counter_tn) if counter_fp + counter_tn != 0 and counter_fp != 0 else 0)
#     prec_fathmm.append(counter_tp / (counter_tp + counter_fp) if counter_tp + counter_fp != 0 and counter_tp != 0 else 0)
#     rec_fathmm.append(counter_tp / (counter_tp + counter_fn))
# roc_auc_fathmm = metrics.auc(fpr_fathmm, tpr_fathmm)
# prc_auc_fathmm = metrics.auc(rec_fathmm, prec_fathmm)
#
# print('FATHMM ROC AUC1: ', roc_auc_fathmm)
# print('FATHMM PRC: ', prc_auc_fathmm)
#
# accuracy = (counter_tp + counter_tn) / (counter_tp + counter_tn + counter_fp + counter_fn)
# print('FATHMM accuracy: ', accuracy)
#
# mcc = (counter_tp*counter_tn - counter_fp*counter_fn) / math.sqrt( (counter_tp + counter_fp)*(counter_tp + counter_fn)*(counter_tn + counter_fp)*(counter_tn + counter_fn) )
# print('MCC FATHMM: ', mcc)
# #############################################################
# '''ma ROC'''
#
# tpr_ma = []
# fpr_ma = []
# prec_ma = []
# rec_ma = []
# thresholds = [i for i in np.linspace(min(ma_scores), max(ma_scores), roc_slicing_numbers)]
# for i in range(len(thresholds)):
#     counter_tp = 0
#     counter_fn = 0
#     counter_fp = 0
#     counter_tn = 0
#     for k in range(len(ma_predicted_pathogenic)):
#         if ma_predicted_pathogenic[k] >= thresholds[i]:
#             counter_tp += 1
#         elif ma_predicted_pathogenic[k] < thresholds[i]:
#             counter_fn += 1
#     # print('threshold: ', thresholds[i])
#     # print('counter_fn: ', counter_fn)
#     # print('counter_tp: ', counter_tp)
#     # print('tpr: ', counter_tp / (counter_tp + counter_fn))
#     tpr_ma.append(counter_tp / (counter_tp + counter_fn) if counter_tp + counter_fn != 0 and counter_tp != 0 else 0)
#
#     for l in range(len(ma_predicted_nonpathogenic)):
#         if ma_predicted_nonpathogenic[l] >= thresholds[i]:
#             counter_fp += 1
#         elif ma_predicted_nonpathogenic[l] < thresholds[i]:
#             counter_tn += 1
#     # print('threshold: ', thresholds[i])
#     # print('counter_fp: ', counter_fp)
#     # print('counter_tn: ', counter_tn)
#     # print('fpr:', counter_fp / (counter_fp + counter_tn))
#     fpr_ma.append(counter_fp / (counter_fp + counter_tn) if counter_fp + counter_tn != 0 and counter_fp != 0 else 0)
#     prec_ma.append(counter_tp / (counter_tp + counter_fp) if counter_tp + counter_fp != 0 and counter_tp != 0 else 0)
#     rec_ma.append(counter_tp / (counter_tp + counter_fn))
#     # print(counter_tp / (counter_tp + counter_fp))
#     # print(counter_tp / (counter_tp + counter_fn))
# # print('tpr: ', tpr_ma)
# # print('fpr: ', fpr_ma)
# roc_auc_ma = metrics.auc(fpr_ma, tpr_ma)
# prc_auc_ma = metrics.auc(rec_ma, prec_ma)
#
# print('MA ROC: ', roc_auc_ma)
# print('MA PRC: ', prc_auc_ma)
#
# accuracy = (counter_tp + counter_tn) / (counter_tp + counter_tn + counter_fp + counter_fn)
# print('MA accuracy: ', accuracy)
#
# mcc = (counter_tp*counter_tn - counter_fp*counter_fn) / math.sqrt( (counter_tp + counter_fp)*(counter_tp + counter_fn)*(counter_tn + counter_fp)*(counter_tn + counter_fn) )
# print('MCC MA: ', mcc)
#############################################################
'''PPh2 ROC'''

tpr_pph2 = []
fpr_pph2 = []
prec_pph2 = []
rec_pph2 = []
thresholds = [i for i in np.linspace(min(pph2_scores), max(pph2_scores), roc_slicing_numbers)]
for i in range(len(thresholds)):
    counter_tp = 0
    counter_fn = 0
    counter_fp = 0
    counter_tn = 0
    for k in range(len(pph2_predicted_pathogenic)):
        if pph2_predicted_pathogenic[k] >= thresholds[i]:
            counter_tp += 1
        elif pph2_predicted_pathogenic[k] < thresholds[i]:
            counter_fn += 1
    '''measuring performance at optimum threshold'''
    if float(f'{thresholds[i]:.2f}') == 0.85:
        print('PPh2 counter_tp', counter_tp)
        print('PPh2 counter_fn', counter_fn)
        print('PPh2 TPR: ', counter_tp / (counter_tp + counter_fn) if counter_tp + counter_fn != 0 and counter_tp != 0 else 0)
    # print('threshold: ', thresholds[i])
    # print('counter_fn: ', counter_fn)
    # print('counter_tp: ', counter_tp)
    # print('tpr PPh2: ', counter_tp / (counter_tp + counter_fn))
    tpr_pph2.append(counter_tp / (counter_tp + counter_fn) if counter_tp + counter_fn != 0 and counter_tp != 0 else 0)

    for l in range(len(pph2_predicted_nonpathogenic)):
        if pph2_predicted_nonpathogenic[l] >= thresholds[i]:
            counter_fp += 1
        elif pph2_predicted_nonpathogenic[l] < thresholds[i]:
            counter_tn += 1
    if float(f'{thresholds[i]:.2f}') == 0.85:
        print('PPh2 counter_tn', counter_tn)
        print('PPh2 counter_fp', counter_fp)
        print('PPh2 FPR: ', counter_fp / (counter_fp + counter_tn) if counter_fp + counter_tn != 0 and counter_fp != 0 else 0)
    # print('threshold: ', thresholds[i])
    # print('counter_fp: ', counter_fp)
    # print('counter_tn: ', counter_tn)
    # print('fpr:', counter_fp / (counter_fp + counter_tn))
    fpr_pph2.append(counter_fp / (counter_fp + counter_tn) if counter_fp + counter_tn != 0 and counter_fp != 0 else 0)
    prec_pph2.append(counter_tp / (counter_tp + counter_fp) if counter_tp + counter_fp != 0 and counter_tp != 0 else 0)
    rec_pph2.append(counter_tp / (counter_tp + counter_fn))
    # print(counter_tp / (counter_tp + counter_fp))
    # print(counter_tp / (counter_tp + counter_fn))
# print('tpr: ', tpr_pph2)
# print('fpr: ', fpr_pph2)
roc_auc_pph2 = metrics.auc(fpr_pph2, tpr_pph2)
prc_auc_pph2 = metrics.auc(rec_pph2, prec_pph2)

print('Polyphen2 ROC: ', roc_auc_pph2)
print('Polyphen2 PRC: ', prc_auc_pph2)

# accuracy = (counter_tp + counter_tn) / (counter_tp + counter_tn + counter_fp + counter_fn)
# print('Polyphen2 accuracy: ', accuracy)
# print('counter_fp: ', counter_fp)
# print('counter_tn: ', counter_tn)
# print('counter_fn: ', counter_fn)
# print('counter_tp: ', counter_tp)
# mcc = (counter_tp*counter_tn - counter_fp*counter_fn) / math.sqrt( (counter_tp + counter_fp)*(counter_tp + counter_fn)*(counter_tn + counter_fp)*(counter_tn + counter_fn) )
# print('MCC PPH2: ', mcc)
#############################################

'''CADD ROC'''
tpr_cadd = []
fpr_cadd = []
prec_cadd = []
rec_cadd = []
thresholds = [i for i in np.linspace(min(cadd_scores), max(cadd_scores), roc_slicing_numbers)]
for i in range(len(thresholds)):
    counter_tp = 0
    counter_fn = 0
    counter_fp = 0
    counter_tn = 0
    for k in range(len(cadd_predicted_pathogenic)):
        if cadd_predicted_pathogenic[k] >= thresholds[i]:
            counter_tp += 1
        elif cadd_predicted_pathogenic[k] < thresholds[i]:
            counter_fn += 1
    '''measuring performance at optimum threshold'''
    if float(f'{thresholds[i]:.0f}') == 15:
        print('CADD counter_tp', counter_tp)
        print('CADD counter_fn', counter_fn)
        print('CADD TPR: ', counter_tp / (counter_tp + counter_fn) if counter_tp + counter_fn != 0 and counter_tp != 0 else 0)
    # print('threshold: ', thresholds[i])
    # print('counter_fn: ', counter_fn) 'CADD (area = %0.2f)' % (prc_auc)
    # print('counter_tp: ', counter_tp)
    # print('tpr CADD: ', counter_tp / (counter_tp + counter_fn))
    tpr_cadd.append(counter_tp / (counter_tp + counter_fn) if counter_tp + counter_fn != 0 and counter_tp != 0 else 0)

    for l in range(len(cadd_predicted_nonpathogenic)):
        if cadd_predicted_nonpathogenic[l] >= thresholds[i]:
            counter_fp += 1
            # print('counter_fp: ', counter_fp)
        elif cadd_predicted_nonpathogenic[l] < thresholds[i]:
            counter_tn += 1
    if float(f'{thresholds[i]:.0f}') == 15:
        print('CADD counter_tn', counter_tn)
        print('CADD counter_fp', counter_fp)
        print('CADD FPR:', counter_fp / (counter_fp + counter_tn) if counter_fp + counter_tn != 0 and counter_fp != 0 else 0)
            # print('counter_tn: ', counter_tn)
    # print('threshold: ', thresholds[i])
    # print('counter_fp: ', counter_fp)
    # print('counter_tn: ', counter_tn)
    # print('fpr:', counter_fp / (counter_fp + counter_tn))
    fpr_cadd.append(counter_fp / (counter_fp + counter_tn) if counter_fp + counter_tn != 0 and counter_fp != 0 else 0)
    # if thresholds[i] == 0.5:
    #     print('FPR at 0.5: ', fpr1[i])
    prec_cadd.append(counter_tp / (counter_tp + counter_fp) if counter_tp + counter_fp != 0 and counter_tp != 0 else 0)
    rec_cadd.append(counter_tp / (counter_tp + counter_fn))
# print('tpr: ', tpr_pph2)
# print('fpr: ', fpr_pph2)
roc_auc_cadd = metrics.auc(fpr_cadd, tpr_cadd)
prc_auc_cadd = metrics.auc(rec_cadd, prec_cadd)

print('CADD ROC: ', roc_auc_cadd)
print('CADD PRC: ', prc_auc_cadd)
# accuracy = (counter_tp + counter_tn) / (counter_tp + counter_tn + counter_fp + counter_fn)
# print('CADD accuracy: ', accuracy)
# print('counter_fp: ', counter_fp)
# print('counter_tn: ', counter_tn)
# print('counter_fn: ', counter_fn)
# print('counter_tp: ', counter_tp)
# mcc = (counter_tp*counter_tn - counter_fp*counter_fn) / math.sqrt( (counter_tp + counter_fp)*(counter_tp + counter_fn)*(counter_tn + counter_fp)*(counter_tn + counter_fn) )
# print('MCC CADD: ', mcc)
#################################################

'''CONDEL ROC'''
tpr_condel = []
fpr_condel = []
prec_condel = []
rec_condel = []
thresholds = [i for i in np.linspace(min(condel_scores), max(condel_scores), roc_slicing_numbers)]  # (0, 1.25, 41) when adding the weights & -137, 25, 41) when multiplying, & (-1.3, 1.1, 41) when using LMT
# print('thresholds: ', thresholds)
for i in range(len(thresholds)):
    counter_tp_condel = 0
    counter_fn_condel = 0
    counter_fp_condel = 0
    counter_tn_condel = 0
    for k in range(len(condel_predicted_pathogenic)):
        if condel_predicted_pathogenic[k] >= thresholds[i]:
            counter_tp_condel += 1
        elif condel_predicted_pathogenic[k] < thresholds[i]:
            counter_fn_condel += 1
    '''measuring performance at optimum threshold'''
    if float(f'{thresholds[i]:.3f}') == 0.522:
        print('CONDEL counter_tp', counter_tp_condel)
        print('CONDEL counter_fn', counter_fn_condel)
        print('CONDEL TPR: ', counter_tp_condel / (counter_tp_condel + counter_fn_condel) if counter_tp_condel + counter_fn_condel != 0 and counter_tp_condel != 0 else 0)
    # print('threshold: ', thresholds[i])
    # print('counter_fn_condel: ', counter_fn_condel)
    # print('counter_tp_condel: ', counter_tp_condel)
    # print('tpr: ', counter_tp_condel / (counter_tp_condel + counter_fn_condel))
    tpr_condel.append(counter_tp_condel / (counter_tp_condel + counter_fn_condel) if counter_tp_condel + counter_fn_condel != 0 and counter_tp_condel != 0 else 0)

    for l in range(len(condel_predicted_nonpathogenic)):
        if condel_predicted_nonpathogenic[l] >= thresholds[i]:
            counter_fp_condel += 1
        elif condel_predicted_nonpathogenic[l] < thresholds[i]:
            counter_tn_condel += 1
    if float(f'{thresholds[i]:.3f}') == 0.522:
        print('CONDEL counter_tn', counter_tn_condel)
        print('CONDEL counter_fp', counter_fp_condel)
        print('CONDEL FPR: ', counter_fp_condel / (counter_fp_condel + counter_tn_condel) if counter_fp_condel + counter_tn_condel != 0 and counter_fp_condel != 0 else 0)
            # print('counter_tn_condel: ', counter_tn_condel)
    # print('threshold: ', thresholds[i])
    # print('counter_fp_condel: ', counter_fp_condel)
    # print('counter_tn_condel: ', counter_tn_condel)
    # print('fpr:', counter_fp_condel / (counter_fp_condel + counter_tn_condel))
    fpr_condel.append(counter_fp_condel / (counter_fp_condel + counter_tn_condel) if counter_fp_condel + counter_tn_condel != 0 and counter_fp_condel != 0 else 0)
    # if thresholds[i] == 0.5:
    #     print('FPR at 0.5: ', fpr1[i])
    prec_condel.append(counter_tp_condel / (counter_tp_condel + counter_fp_condel) if counter_tp_condel + counter_fp_condel != 0 and counter_tp_condel != 0 else 0)
    rec_condel.append(counter_tp_condel / (counter_tp_condel + counter_fn_condel))

roc_auc_condel = metrics.auc(fpr_condel, tpr_condel)
prc_auc_condel = metrics.auc(rec_condel, prec_condel)

print('CONDEL ROC: ', roc_auc_condel)
print('CONDEL PRC: ', prc_auc_condel)

# accuracy = (counter_tp_condel + counter_tn_condel) / (counter_tp_condel + counter_tn_condel + counter_fp_condel + counter_fn_condel)
# mcc_condel = ((counter_tp_condel*counter_tn_condel) - (counter_fp_condel*counter_fn_condel)) / math.sqrt( (counter_tp_condel + counter_fp_condel)*(counter_tp_condel + counter_fn_condel)*(counter_tn_condel + counter_fp_condel)*(counter_tn_condel + counter_fn_condel) )
# print('MCC CONDEL: ', mcc_condel)

####################################################
'''Re-plot all ROCs using manually calculated fpr & tpr added to it is the manually calculated ROC AUC: 

this is, different using a slicing threshold of 21, equal to the one calculated by python (metrics.auc(fpr, tpr)) when using a threshold of 2001'''

# plt.plot(fpr_provean, tpr_provean, label='Provean (area = %0.2f)' % (roc_auc_provean))
# plt.plot(fpr_sift, tpr_sift, label='SIFT (area = %0.2f)' % (roc_auc_sift))
# plt.plot(fpr_fathmm, tpr_fathmm, label='FATHMM (area = %0.2f)' % (roc_auc_fathmm))
# plt.plot(fpr_ma, tpr_ma, label='MA (area = %0.2f)' % (roc_auc_ma))
# plt.plot(fpr_pph2, tpr_pph2, label='PolyPhen2 (area = %0.2f)' % (roc_auc_pph2))
# plt.plot(fpr_cadd, tpr_cadd, label='CADD (area = %0.2f)' % (roc_auc_cadd))
# plt.plot(fpr_condel, tpr_condel, label='CONDEL (area = %0.2f)' % (roc_auc_condel))
# plt.legend(loc='lower right')
# # plt.plot([0, 1], [0, 1], 'r--')
# plt.xlim([0.0, 1.02])
# plt.ylim([0.0, 1.02])
# plt.ylabel('True Positive Rate (sensitivity)')
# plt.xlabel('False Positive Rate (1-specificity)')
# plt.show()

# # create the axis of thresholds (scores)
# plt.plot(fpr, tpr)  # add the ROC curve to it
# ax2 = plt.gca().twinx()
# ax2.plot(fpr, thresholds, markeredgecolor='r', linestyle='dashed', color='r')
# ax2.set_ylabel('Threshold', color='r')
# ax2.set_ylim([thresholds[0], thresholds[-1]])
# ax2.set_xlim([fpr[0], fpr[-1]])
# plt.show()
####################################
'''PRC of all tools using manually-calculated precision & recall added to it is the manually calculated PRC AUC

However, these results vary from those calculated by Python using its functions, above, possibly because of metrics.auc() not being accurate (refer to documentation)'''

# plt.plot(rec_provean, prec_provean, label='Provean (area = %0.2f)' % (prc_auc_provean))
# plt.plot(rec_sift, prec_sift, label='SIFT (area = %0.2f)' % (prc_auc_sift))
# plt.plot(rec_fathmm, prec_fathmm, label='FATHMM (area = %0.2f)' % (prc_auc_fathmm))
# plt.plot(rec_ma, prec_ma, label='MA (area = %0.2f)' % (prc_auc_ma))
# plt.plot(rec_pph2, prec_pph2, label='PolyPhen2 (area = %0.2f)' % (prc_auc_pph2))
# plt.plot(rec_cadd, prec_cadd, label='CADD (area = %0.2f)' % (prc_auc_cadd))
# plt.plot(rec_condel, prec_condel, label='CONDEL (area = %0.2f)' % (prc_auc_condel))
#
# # plt.plot([0, 1], [0, 1], 'r--')
# plt.legend(loc='best')
# plt.xlim([0.0, 1.02])
# plt.ylim([0.0, 1.02])
# plt.ylabel('Precision')
# plt.xlabel('Recall')
# plt.show()
######################################################


# sns.distplot(pathogenic_var_scores, kde=True, label='Pathogenic')
# sns.distplot(nonPathogenic_var_scores, kde=True, label='Non-pathogenic')
# plt.title('Classifier scores using WEKA-given weights')
# plt.xlabel('Pathogenicity score')
# plt.ylabel('Frequency')
# plt.show()

# sns.countplot(hgmd_var_scores, color='red', alpha=myalpha)
# sns.countplot(exac_var_scores, color='green', alpha=myalpha)
# plt.show()

# sns.boxplot(hgmd_var_scores, color='red', fliersize=0.3, width=0.3)
# sns.boxplot(exac_var_scores, color='green', fliersize=0.5, width=0.5)
# plt.show()
