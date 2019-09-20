#! /usr/bin/env python3

# import matplotlib.pyplot as plt
# from matplotlib import style
import numpy as np
import pandas as pd
# from pandas import Series, DataFrame
# from sklearn.metrics import auc
from scipy.stats import *
# from scipy import stats
from mlxtend.evaluate import mcnemar_table, mcnemar
from sklearn.metrics import matthews_corrcoef, roc_auc_score
from mlxtend.plotting import checkerboard_plot
import matplotlib.pyplot as plt
''' 19/10/18
McNemar chi squared test for spec/sens analysis of the different tests'''
#######################
'''Analysis of Sensitivity & Specificity using Chi squared:

If you have two independent samples, then you can use Pearson's chi-square* to compare the two sensitivities. But if you have one sample of individuals who were given both diagnostic tests, then you have paired observations, and need to use the McNemar chi-square. Here are examples using SAS:
https://onlinecourses.science.psu.edu/stat509/node/152

If the expected counts are too low to justify using Pearson's chi-square, use the N-1 chi-square instead. See http://www.iancampbell.co.uk/twobytwo/twobytwo.htm for details and an online calculator.

If diagnostic tests were studied on two independent groups of patients, then two-sample tests for binomial proportions are appropriate (chi-square, Fisher's exact test). If both diagnostic tests were performed on each patient, then paired data result and methods that account for the correlated binary outcomes are necessary (McNemar's test)'''
#######################

# model_data = pd.read_csv('/Users/mdefsss2/cacna1f/E_J_human_genetics/weka_results/logistic_predictions.csv')
data = pd.read_csv('/Users/mdefsss2/cacna1f/E_J_human_genetics/weka_results/combined_train_test_data.csv')

# class_model = data['class']
model= data["'predicted class'"]
# predictions_model_binary= data["binary_predicted_class"]
# classifier_scores = model_data['score']
true_class_binary = data['binary_class']
true_class = data['class']
sift = data['SIFT SCORE']
pph2 = data[' pph2_prob']
cadd = data['PHRED']
condel = data['CONDEL']
# provean = predictions_others['Provean SCORE']
# ma = predictions_others['MA']
# fathmm = predictions_others['FATHMM']

# classifier_predictions = model_data['overall_prediction']

# target_index = predictions_others.columns.get_loc('class')
# sift_index = predictions_others.columns.get_loc('sift_binary')
# pph2_index = predictions_others.columns.get_loc('pph2_binary')
# cadd_index = predictions_others.columns.get_loc('cadd_binary')
# condel_index = predictions_others.columns.get_loc('condel_binary')
# model_prediction_index = model_data.columns.get_loc('binary_prediction')
# model_target_index = model_data.columns.get_loc('class')
# binary_model = []
# for i in true_class:
#     if i == 'pathogenic':
#         i = 1
#         binary_model.append(i)
#     elif i == 'benign':
#         i = 0
#         binary_model.append(i)

model_binary = []
counter_1 = 0
counter_0 = 0
for i in model:
    if i == 'pathogenic':
        i = 1
        counter_1 +=1
        model_binary.append(i)
    elif i == 'benign':
        # print(i)
        i = 0
        counter_0 +=1
        model_binary.append(i)
# print('number of predicted and true (model): ', len(model_binary), len(class_model))
print('matthews_corr_coef (model): ', matthews_corrcoef(true_class_binary, model_binary))
#####################################################################

cadd_binary = []
counter_1 = 0
counter_0 = 0
for i in cadd:
    if i < 15:
        i = 0
        counter_0 +=1
        cadd_binary.append(i)
    elif i >= 15:
        # print(i)
        i = 1
        counter_1 +=1
        cadd_binary.append(i)
print('number of predicted and true (CADD): ', len(cadd_binary), len(true_class_binary))
print('matthews_corr_coef (CADD): ', matthews_corrcoef(true_class_binary, cadd_binary))
#####################################################################

sift_binary = []
counter_1 = 0
counter_0 = 0
for i in sift:
    if i <= 0.05:
        i = 1
        counter_1 += 1
        sift_binary.append(i)
        # print(1)
    elif i > 0.05:
        i = 0
        counter_0 += 1
        sift_binary.append(i)
        # print(0)
print('matthews_corr_coef (SIFT): ', matthews_corrcoef(true_class_binary, sift_binary))
#####################################################################

pph2_binary = []
counter_1 = 0
counter_0 = 0
for i in pph2:
    if i >= 0.85:
        i = 1
        counter_1 += 1
        pph2_binary.append(i)
        # print(1)
    elif i < 0.85:
        i = 0
        counter_0 += 1
        pph2_binary.append(i)
        # print(0)
print('matthews_corr_coef (PPh2): ', matthews_corrcoef(true_class_binary, pph2_binary))
#####################################################################

condel_binary = []
counter_1 = 0
counter_0 = 0
for i in condel:
    if i >= 0.522:
        i = 1
        counter_1 += 1
        condel_binary.append(i)
        # print(1)
    elif i < 0.522:
        i = 0
        counter_0 += 1
        condel_binary.append(i)
        # print(0)
print('matthews_corr_coef (condel): ', matthews_corrcoef(true_class_binary, condel_binary))
#################################################################
'''forming a confusion matrix/contig table of the classifier & each of the other tests, and calculate Chi squared test, i.e. if the predictions from one classifier are significantly different from those of the other classifier (p<0.05) 
(IMPORTANT: use exact=True for sample sizes < 25 in boxes b+c, ie. upper right + lower left boxes in the confusion matrix)

[https://rasbt.github.io/mlxtend/user_guide/evaluate/mcnemar/]'''
#####################################################################
'''the target/real scores of the classifier'''
# target_scores = []
# for i in range(len(model_data)):
#     k = model_data.iloc[i, model_target_index]
#     target_scores.append(k)
# print(target_scores)
# target_scores_array = np.array(target_scores)
# print(target_scores_array)
#
# '''the classifier'''
# classifier_scores = []
# for i in range(len(model_data)):
#     k = model_prediction_index.iloc[i, model_prediction_index]
#     classifier_scores.append(k)
# # print(classifier_scores)
# classifier_scores_array = np.array(classifier_scores)
# print(classifier_scores_array)

'''FOR all other classifiers in a different file/format'''
# target_scores = []
# for i in range(len(predictions_others)):
#     k = predictions_others.iloc[i, target_index]
#     target_scores.append(k)
# # print(target_scores)
# target_scores_array = np.array(target_scores)
# print(target_scores_array)

'''SIFT'''
# sift_scores = []
# for i in range(len(predictions_others)):
#     k = predictions_others.iloc[i, sift_index]
#     sift_scores.append(k)
# sift_scores_array = np.array(sift_scores)

sift_and_model = mcnemar_table(y_target=np.array(true_class_binary), y_model1=np.array(model_binary), y_model2=np.array(sift_binary))
print('model & sift: ', '\n', sift_and_model)
chi2, p = mcnemar(ary=sift_and_model, corrected=True)
print(' chi_squared: ', chi2)
print(' p-value: ', p)

brd = checkerboard_plot(sift_and_model,
                        figsize=(2, 2),
                        fmt='%d',
                        col_labels=['model 2 wrong', 'model 2 right'],
                        row_labels=['model 1 wrong', 'model 1 right'])
plt.show()

'''PPH2'''

pph2_and_model = mcnemar_table(y_target=np.array(true_class_binary), y_model1=np.array(model_binary), y_model2=np.array(pph2_binary))
print('model & pph2: ', '\n', pph2_and_model)
chi2, p = mcnemar(ary=pph2_and_model, corrected=True)
print(' chi_squared: ', chi2)
print(' p-value: ', p)

brd = checkerboard_plot(pph2_and_model,
                        figsize=(2, 2),
                        fmt='%d',
                        col_labels=['model 2 wrong', 'model 2 right'],
                        row_labels=['model 1 wrong', 'model 1 right'])
plt.show()

'''CADD'''

cadd_and_model = mcnemar_table(y_target=np.array(true_class_binary), y_model1=np.array(model_binary), y_model2=np.array(cadd_binary))
print('model & cadd: ', '\n', cadd_and_model)
chi2, p = mcnemar(ary=cadd_and_model, corrected=True)
print(' chi_squared: ', chi2)
print(' p-value: ', p)

brd = checkerboard_plot(cadd_and_model,
                        figsize=(2, 2),
                        fmt='%d',
                        col_labels=['model 2 wrong', 'model 2 right'],
                        row_labels=['model 1 wrong', 'model 1 right'])
plt.show()

'''CONDEL'''

condel_and_model = mcnemar_table(y_target=np.array(true_class_binary), y_model1=np.array(model_binary), y_model2=np.array(condel_binary))
print('model & condel: ', '\n', condel_and_model)
chi2, p = mcnemar(ary=condel_and_model, corrected=True)
print(' chi_squared: ', chi2)
print(' p-value: ', p)

brd = checkerboard_plot(condel_and_model,
                        figsize=(2, 2),
                        fmt='%d',
                        col_labels=['model 2 wrong', 'model 2 right'],
                        row_labels=['model 1 wrong', 'model 1 right'])
plt.show()
