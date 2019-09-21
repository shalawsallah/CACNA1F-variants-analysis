#! /usr/bin/env python3
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics
from scipy.stats import *
from matplotlib import rc
import math

all_data = pd.read_csv('/Users/mdefsss2/cacna1f/WEKA_old_results/combined_train_test_data.csv')
model_data = pd.read_csv('/Users/mdefsss2/cacna1f/E_J_human_genetics/weka_results/prediction_AUCs.csv')
model_predictions = pd.read_csv('/Users/mdefsss2/cacna1f/E_J_human_genetics/weka_results/logistic_predictions.csv')
################################################

roc_slicing_numbers = 1001
print('number of variants: ', len(all_data))
thresholds = model_data['Threshold']
model_tpr = model_data["'True Positive Rate'"]
model_fpr = model_data["'False Positive Rate'"]
model_prec = model_data["Precision"]
model_rec = model_data["Recall"]

variants = all_data['variants']
target_class = all_data['class'].tolist()
target_class_binary = all_data['binary_class'].tolist()
pph2_scores = all_data[' pph2_prob'].tolist()
cadd_scores = all_data['PHRED'].tolist()
condel_scores = all_data['CONDEL'].tolist()
sift_scores = all_data['1-SIFT Score'].tolist()

for i in range(len(variants)):
    if target_class[i] == 'pathogenic':
        pph2_predicted_pathogenic.append(pph2_scores[i])
        cadd_predicted_pathogenic.append(cadd_scores[i])
        condel_predicted_pathogenic.append(condel_scores[i])
        sift_predicted_pathogenic.append(sift_scores[i])
    elif target_class[i] == 'benign':
        pph2_predicted_nonpathogenic.append(pph2_scores[i])
        cadd_predicted_nonpathogenic.append(cadd_scores[i])
        condel_predicted_nonpathogenic.append(condel_scores[i])
        sift_predicted_nonpathogenic.append(sift_scores[i])

mynormed = 1
myalpha = 0.3
mybins = 10

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

##############################################
'''ROCs'''

'''WEKA model'''

model_roc_auc = metrics.auc(model_fpr, model_tpr)
print('model ROC: ', model_roc_auc)
####################################################
'''plot all ROCs using python-calculated fpr & tpr & ROC AUC'''


plt.plot(model_fpr, model_tpr, label='CACNA1F-vp (area = %0.2f)' % (model_roc_auc))

fpr, tpr, thresholds = metrics.roc_curve(target_class_binary, sift_scores, pos_label=1)
roc_auc = metrics.auc(fpr, tpr)
plt.plot(fpr, tpr, label='SIFT (area = %0.2f)' % (roc_auc))

fpr, tpr, thresholds = metrics.roc_curve(target_class_binary, pph2_scores, pos_label=1)
roc_auc = metrics.auc(fpr, tpr)
plt.plot(fpr, tpr, label='PolyPhen2 (area = %0.2f)' % (roc_auc))

fpr, tpr, thresholds = metrics.roc_curve(target_class_binary, cadd_scores, pos_label=1)
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
pph2_roc = metrics.roc_auc_score(target_class_binary, pph2_scores)
print('PPh2_roc score: ', pph2_roc)
cadd_roc = metrics.roc_auc_score(target_class_binary, cadd_scores)
print('CADD_roc score: ', cadd_roc)
condel_roc = metrics.roc_auc_score(target_class_binary, condel_scores)
print('CONDEL_roc score: ', condel_roc)

##############################
'''PRC of all tools using python-calculated precision & recall & PRC AUC'''

model_prc_auc = metrics.auc(model_rec, model_prec)
print('model PRC: ', model_prc_auc)
plt.plot(model_rec, model_prec, label='CACNA1F-vp (area = %0.2f)' % (model_prc_auc))

precision, recall, thresholds = metrics.precision_recall_curve(target_class_binary, sift_scores, pos_label=1)
prc_auc = metrics.auc(recall,precision)
plt.plot(recall, precision, label='SIFT (area = %0.2f)' % (prc_auc))

precision, recall, thresholds = metrics.precision_recall_curve(target_class_binary, pph2_scores, pos_label=1)
prc_auc = metrics.auc(recall,precision)
plt.plot(recall, precision, label='PolyPhen2 (area = %0.2f)' % (prc_auc))

precision, recall, thresholds = metrics.precision_recall_curve(target_class_binary, cadd_scores, pos_label=1)
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
pph2_precision = metrics.average_precision_score(np.array(target_class_binary), np.array(pph2_scores))
print('PPh2 PRC score: ', pph2_precision)
cadd_precision = metrics.average_precision_score(np.array(target_class_binary), np.array(cadd_scores))
print('CADD PRC score: ', cadd_precision)
condel_precision = metrics.average_precision_score(np.array(target_class_binary), np.array(condel_scores))
print('CONDEL PRC score: ', condel_precision)

##################################################################
