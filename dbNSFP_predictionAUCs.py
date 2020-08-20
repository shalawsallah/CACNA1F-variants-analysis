#! /usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.metrics import matthews_corrcoef

'''predictions from https://sites.google.com/site/jpopgen/dbNSFP. info about predictors @ https://sites.google.com/site/jpopgen/dbNSFP'''
genes = ['cryba1', 'cryba2', 'cryba4', 'crybb1', 'crybb2', 'crybb3', 'crygc', 'crygd', 'crygs']# gene = 'scn1a'
# genes = ['scn1a','scn2a','scn4a','scn5a','scn8a','scn9a','scn10a','cacna1a','cacna1c','cacna1h', 'gjb1', 'gjb2', 'gjb3', 'gjb4', 'gjb6', 'gjc2', 'gja3', 'gja8', 'best1', 'pax6', 'rho']
group_1_predictors_only = False #the tools with one prediction for each variant

def roc(tool, predictions, scores_direction=int()):  # score _direction = 1 or 0
    '''plotting a ROC curve for the predictions of the tool'''
    try:
        fpr, tpr, thresholds = metrics.roc_curve(class_binary, predictions, pos_label=scores_direction)
        roc_auc = metrics.auc(fpr, tpr)
        # fig, ax = plt.subplots()
        plt.plot(fpr, tpr, label=str(tool) + ' (area = %0.2f)' % (roc_auc))
        # ax.legend()
        plt.legend(loc='lower right')
        plt.plot([0, 1], [0, 1], 'r--')
        plt.xlim([0.0, 1.02])
        plt.ylim([0.0, 1.02])
        plt.ylabel('True Positive Rate (sensitivity)')
        plt.xlabel('False Positive Rate (1-specificity)')
        plt.title('')
    except ValueError as e:
        print(e)
        pass

SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE)
plt.rc('axes', labelsize=BIGGER_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)
mydensity = 1
myalpha = 0.8

def prc(tool, predictions, scores_direction=int()):
    '''plotting a PR curve for the predictions of the tool'''
    try:
        precision, recall, thresholds = metrics.precision_recall_curve(class_binary, predictions,
                                                                       pos_label=scores_direction)
        prc_auc = metrics.auc(recall, precision)
        plt.plot(recall, precision, label=str(tool) + ' (area = %0.2f)' % (prc_auc))
        # ax.legend()
        plt.legend(loc='lower left')
        # plt.plot([0, 1], [0, 1], 'r--')
        plt.xlim([0.0, 1.02])
        plt.ylim([0.0, 1.02])
        plt.ylabel('Precision')
        plt.xlabel('Recall')
        plt.title('')
    except ValueError as e:
        print(e)
        pass

def tool_scores_histogram(coordinates, tool, predictions):
    '''plotting a histogram of the scores'''
    pathogenic_scores = []
    nonpathogenic_scores = []
    for i in range(len(coordinates)):
        if var_class[i] == 'pathogenic':
            pathogenic_scores.append(predictions[i])
        elif var_class[i] == 'benign':
            nonpathogenic_scores.append(predictions[i])
    fig, ax = plt.subplots()
    ax.hist(pathogenic_scores, density=mydensity, histtype='step', stacked=True, alpha=myalpha,
            label="Dataset P")
    ax.hist(nonpathogenic_scores, density=mydensity, histtype='step', stacked=True, alpha=myalpha,
            label="Dataset N")
    ax.legend()
    plt.xlabel(tool)
    plt.ylabel('Frequency')
    plt.legend(loc='upper left')
    plt.title('')
    # fig.tight_layout()

def mcc(true_binary, tool_scores, scores_direction=1, pathogenicity_threshold=0.5):
    '''MCC'''
    pred = []
    for i in range(len(tool_scores)):
        if scores_direction==1:
            if float(tool_scores[i]) >= pathogenicity_threshold:
                pred.append(1)
            elif float(tool_scores[i]) < pathogenicity_threshold:
                pred.append(0)
        elif scores_direction==0:
            if float(tool_scores[i]) >= pathogenicity_threshold:
                pred.append(0)
            elif float(tool_scores[i]) < pathogenicity_threshold:
                pred.append(1)
    return matthews_corrcoef(true_binary, pred)

#####################################################################
data = pd.read_table('/Users/mdefsss2/dbNSFP4/crystallins.out')
# data = pd.read_excel('/Users/mdefsss2/Desktop/test.xlsx')
genes_transcripts = pd.read_excel('/Users/mdefsss2/other_genes/genes_list.xlsx')
# variants = data['variants']
var_class = data['class']
print('All variants with TEST data: ', len(var_class))
class_binary = [data['class_binary'][m] for m in range(len(var_class)) if var_class[m] == 'pathogenic' or var_class[m] == 'benign']
print('Variants without TEST data: ', len(class_binary))
# coordinates = [data['pos(1-based)'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
ensembl_transcripts = [data['Ensembl_transcriptid'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
sift = [data['SIFT_score'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
pph2 = [data['Polyphen2_HVAR_score'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
# lrt = data['LRT_score']
# mut_taster = data['MutationTaster_score']
# mut_assessor = data['MutationAssessor_score']
fathmm = [data['FATHMM_score'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
provean = [data['PROVEAN_score'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
vest4 = [data['VEST4_score'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
metaSVM = [data['MetaSVM_score'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
# metaLR = data['MetaLR_score']
revel = [data['REVEL_score'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
cadd = [data['CADD_phred'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
dann = [data['DANN_score'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
gerp = [data['GERP++_RS'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
phylop = [data['phyloP100way_vertebrate'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
deogen2 = [data['DEOGEN2_score'][m] for m in range(len(var_class)) if var_class[m] =='pathogenic' or var_class[m] == 'benign']
# phastCons = data['phastCons100way_vertebrate']
# mvp = data['MVP_score']
# mpc = data['MPC_score']
# primateAI = data['PrimateAI_score']
sift_predictions = []
pph2_predictions = []
fathmm_predictions = []
provean_predictions = []
vest4_predictions = []
deogen2_predictions = []
for j in genes:
    transcript = ''
    for i in range(len(genes_transcripts['gene'])): #obtaining the trascript used
        if genes_transcripts['gene'][i] == j:
            transcript += genes_transcripts['transcript'][i].strip()
    print('transcript retrieved: ', transcript)

    if group_1_predictors_only == False:
        prediction_counted = False
        for i in range(len(ensembl_transcripts)):
            strings = str(ensembl_transcripts[i])
            if ';' in strings:
                # print(strings)
                prediction_counted = True
                transcripts_available = strings.split(';')
                # num_of_predictions = ensembl_transcripts[i].split(';')
                sift_pred = sift[i].split(';')
                pph2_pred = pph2[i].split(';')
                fathmm_pred = fathmm[i].split(';')
                provean_pred = provean[i].split(';')
                vest4_pred = vest4[i].split(';')
                deogen2_pred = deogen2[i].split(';')
                for k in range(len(transcripts_available)):
                    if str(transcripts_available[k]) == str(transcript):
                        try:
                            sift_predictions.append(float(sift_pred[k]))
                        except ValueError as e:
                            print(e)  # the transcript of interest wasnt available & the first transcript prediction is therefore chosen instead
                            # sift_predictions.append(float(sift_pred[k+1]))
                        try:
                            pph2_predictions.append(float(pph2_pred[k]))
                        except ValueError as e:
                            print(e)
                            # if type(pph2_pred[k+1]) == float:
                            #     pph2_predictions.append(float(pph2_pred[k+1]))
                        else:
                            pass
                        try:
                            fathmm_predictions.append(float(fathmm_pred[k]))
                        except ValueError as e:
                            print(e)
                            # fathmm_predictions.append(float(fathmm_pred[k+1]))
                        try:
                            provean_predictions.append(float(provean_pred[k]))
                        except ValueError as e:
                            print(e)
                            # provean_predictions.append(float(provean_pred[k+1]))
                        try:
                            vest4_predictions.append(float(vest4_pred[k]))
                        except ValueError as e:
                            print(e)
                            vest4_predictions.append(float(vest4_pred[k+1]))
                        try:
                            deogen2_predictions.append(float(deogen2_pred[k]))
                        except ValueError as e:
                            print(e)
                    # if transcript not in num_of_predictions:
                    #     print('Transcript of inerest MISSING!! at: ', i)
                    # else:
                    #     print('transcript not found among predictions for variant: ', i)
            elif ';' not in strings:
                if strings == transcript:
                    # prediction_counted = True
                    sift_predictions.append(float(sift[i]))
                    pph2_predictions.append(float(pph2[i]))
                    fathmm_predictions.append(float(fathmm[i]))
                    provean_predictions.append(float(provean[i]))
                    vest4_predictions.append(float(vest4[i]))
                    deogen2_predictions.append(float(deogen2[i]))
            else:
                print(ValueError)
            #     # print(i)
            #     sift_predictions=sift
            #     pph2_predictions=pph2
            #     fathmm_predictions=fathmm
            #     provean_predictions=provean
            #     vest4_predictions=vest4
            #     # vest4_predictions=[float(vest4[m]) for m in vest4]

print('number of predictions (set1): ', len(sift_predictions), len(fathmm_predictions), len(provean_predictions), len(pph2_predictions), len(vest4_predictions), len(deogen2_predictions))

print('number of predictions (set2): ', len(metaSVM), len(revel), len(cadd), len(dann), len(gerp), len(phylop))

# tool_scores_histogram(coordinates, 'FATHMM', fathmm_predictions)
# plt.show()
# tool_scores_histogram(coordinates, 'Provean', provean_predictions)
# plt.show()

roc('REVEL', revel, 1)
roc('MetaSVM', metaSVM, 1)
roc('CADD', cadd, 1)
roc('DANN', dann, 1)
if group_1_predictors_only == False:
    roc('VEST4', vest4_predictions, 1)
    roc('PolyPhen2', pph2_predictions, 1)
    roc('SIFT', sift_predictions, 0)
    roc('PROVEAN', provean_predictions, 0)
    roc('FATHMM', fathmm_predictions, 0)
    roc('Deogen2', deogen2_predictions, 1)
# roc('LRT', lrt, 0)
# roc('PrimateAI', primateAI, 1)
roc('GERP++', gerp, 1)
roc('Phylop100', phylop, 1)
# plt.savefig('/Users/mdefsss2/x_linked_genes/'+j+'/ROCs_', dpi=300, transparent=True, bbox_inches='tight')
plt.show()

prc('REVEL', revel, 1)
prc('MetaSVM', metaSVM, 1)
prc('CADD', cadd, 1)
prc('DANN', dann, 1)
if group_1_predictors_only == False:
    prc('VEST4', vest4_predictions, 1)
    prc('PolyPhen2', pph2_predictions, 1)
    # prc('SIFT', sift_predictions, 0)
    # prc('PROVEAN', provean_predictions, 0)
    # prc('FATHMM', fathmm_predictions, 0)
    prc('Deogen2', deogen2_predictions, 1)
# prc('LRT', lrt, 0)
# prc('PrimateAI', primateAI, 1)
prc('GERP++', gerp, 1)
prc('Phylop100', phylop, 1)
# plt.savefig('/Users/mdefsss2/x_linked_genes/' + j + '/PRCs_', dpi=300, transparent=True, bbox_inches='tight')
plt.show()

'''MCC'''
print('REVEL MCC: ', mcc(class_binary, revel))
print('MetaSVM MCC: ', mcc(class_binary, metaSVM, pathogenicity_threshold=0))
print('CADD MCC: ', mcc(class_binary, cadd, pathogenicity_threshold=15))
print('DANN MCC: ', mcc(class_binary, dann))
if group_1_predictors_only == False:
    print('VEST4 MCC: ', mcc(class_binary, vest4_predictions))
    print('PPh2 MCC: ', mcc(class_binary, pph2_predictions, pathogenicity_threshold=0.85))
    # print('SIFT MCC: ', mcc(class_binary, sift_predictions))
    # print('PROVEAN MCC: ', mcc(class_binary, provean_predictions))
    # print('FATHMM MCC: ', mcc(class_binary, fathmm_predictions))
    print('DEOGEN2 MCC: ', mcc(class_binary, deogen2_predictions))
plt.close('all')
'''END of running/functioning code'''
###########################################
# three_letters_to_one_letter = {'Phe': 'F', 'Tyr': 'Y', 'Leu': 'L', 'His': 'H', 'Gln': 'Q', 'Ile': 'I', 'Asn': 'N',
#                                'Met': 'M', 'Val': 'V', 'Asp': 'D', 'Glu': 'E', 'Ser': 'S', 'Pro': 'P', 'Arg': 'R',
#                                'Thr': 'T', 'Lys': 'K', 'Gly': 'G', 'Ala': 'A', 'Cys': 'C', 'Trp': 'W'}
# one_to_three_letter = {v:k for k,v in three_letters_to_one_letter.items()}
# aa_positions = []
# for i in data['aapos']:
#     if type(i) != int:
#         each = i.split(';')
#         pos = each[0]
#         aa_positions.append(pos)
#     else:
#         aa_positions.append(i)
'''to import the class of the corresponding variant predicted from the original file'''
# for k in range(len(coordinates)):
#     flag = False
#     variant = str(one_to_three_letter[data['aaref'][k]])+str(aa_positions[k])+str(one_to_three_letter[data['aaalt'][k]])
#     for i in range(len(variants)):
#         if variant == variants[i]:
#             flag = True
#             print(var_class[i])
#     if flag == False:
#         print('none')
'''the above code needs improvement... to account for the duplicate variants for which the class is counted/printed twice'''
