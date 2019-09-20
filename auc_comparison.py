#! /usr/bin/env python3

# import matplotlib.pyplot as plt
# from matplotlib import style
import numpy as np
import pandas as pd
from scipy.stats import *
from mlxtend.evaluate import mcnemar_table, mcnemar
from sklearn.metrics import matthews_corrcoef, roc_auc_score
from mlxtend.plotting import checkerboard_plot
import matplotlib.pyplot as plt
''' 19/10/18
McNemar chi square test for spec/sens analysis of the different tools'''

data = pd.read_csv('/Users/mdefsss2/cacna1f/E_J_human_genetics/weka_results/combined_train_test_data.csv')


model= data["'predicted class'"]
true_class_binary = data['binary_class']
true_class = data['class']
sift = data['SIFT SCORE']
pph2 = data[' pph2_prob']
cadd = data['PHRED']
condel = data['CONDEL']

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

'''SIFT'''


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
