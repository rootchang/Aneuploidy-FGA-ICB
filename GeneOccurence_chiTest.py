from collections import Counter
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import multipletests
import pandas as pd
import numpy as np

# gene mutation info.
geneMut_fn = './tmb_mskcc_2018/data_mutations.txt'
sample_mutation_dict = {}
gene_dict = {}
for line in open(geneMut_fn, 'r').readlines()[1:]:
    words= line.strip('\n').split('\t')
    geneNA = words[0]
    gene_dict[geneNA] = 1
    SAMPLE_ID = words[16]
    if SAMPLE_ID not in sample_mutation_dict:
        sample_mutation_dict[SAMPLE_ID] = {geneNA:1}
    else:
        sample_mutation_dict[SAMPLE_ID][geneNA] = 1
print('Total sample number: ', len(sample_mutation_dict))
print('Total gene number: ', len(gene_dict))

fnOut = open('Differential_co_occurrence_genes.txt','w')
fnOut.write('Gene\tFGA/AS_H\tFGA/AS_L\tpval\tadj_pval\n')

CNA_fn = 'MSK_allInfo.csv'
CNAscore_info = pd.read_csv(CNA_fn)

################################################## FGA0.2 ##################################################
melanoma_FGA2H = CNAscore_info.loc[np.logical_and(CNAscore_info['CANCER_TYPE']=='Melanoma', CNAscore_info['FGA0.2_score']=='high'), 'SAMPLE_ID'].tolist()
melanoma_FGA2L = CNAscore_info.loc[np.logical_and(CNAscore_info['CANCER_TYPE']=='Melanoma', CNAscore_info['FGA0.2_score']=='low'), 'SAMPLE_ID'].tolist()
NSCLC_FGA2H = CNAscore_info.loc[np.logical_and(CNAscore_info['CANCER_TYPE']=='Non-Small Cell Lung Cancer', CNAscore_info['FGA0.2_score']=='high'), 'SAMPLE_ID'].tolist()
NSCLC_FGA2L = CNAscore_info.loc[np.logical_and(CNAscore_info['CANCER_TYPE']=='Non-Small Cell Lung Cancer', CNAscore_info['FGA0.2_score']=='low'), 'SAMPLE_ID'].tolist()

melanoma_FGA2H_genes = [list(sample_mutation_dict[sample].keys()) for sample in melanoma_FGA2H if sample in sample_mutation_dict]
melanoma_FGA2H_genes = sum(melanoma_FGA2H_genes,[])
melanoma_FGA2L_genes = [list(sample_mutation_dict[sample].keys()) for sample in melanoma_FGA2L if sample in sample_mutation_dict]
melanoma_FGA2L_genes = sum(melanoma_FGA2L_genes,[])
NSCLC_FGA2H_genes = [list(sample_mutation_dict[sample].keys()) for sample in NSCLC_FGA2H if sample in sample_mutation_dict]
NSCLC_FGA2H_genes = sum(NSCLC_FGA2H_genes,[])
NSCLC_FGA2L_genes = [list(sample_mutation_dict[sample].keys()) for sample in NSCLC_FGA2L if sample in sample_mutation_dict]
NSCLC_FGA2L_genes = sum(NSCLC_FGA2L_genes,[])

# Create a Counter object to count the occurrences of each element
counts = Counter(NSCLC_FGA2H_genes)
total_count = len(NSCLC_FGA2H)
NSCLC_FGA2H_genes_prob = {element: count / total_count * 100 for element, count in counts.items()}

counts = Counter(NSCLC_FGA2L_genes)
total_count = len(NSCLC_FGA2L)
NSCLC_FGA2L_genes_prob = {element: count / total_count * 100 for element, count in counts.items()}

# merge the two dictionaries
NSCLC_genes_prob = {}
for key in NSCLC_FGA2H_genes_prob.keys() | NSCLC_FGA2L_genes_prob.keys():
    NSCLC_genes_prob[key] = [0 if NSCLC_FGA2H_genes_prob.get(key) is None else NSCLC_FGA2H_genes_prob.get(key), 0 if NSCLC_FGA2L_genes_prob.get(key) is None else NSCLC_FGA2L_genes_prob.get(key)]
# chi-test for the significant difference of gene mutation presence
mutations = []
significances = []
for gene in NSCLC_genes_prob:
    occurrence_H = NSCLC_genes_prob[gene][0]
    occurrence_L = NSCLC_genes_prob[gene][1]
    mutation_table = [[occurrence_H, 100 - occurrence_H], [occurrence_L, 100 - occurrence_L]]
    chi2_stat, p_val, dof, exp_freq = chi2_contingency(mutation_table)
    mutations.append([gene, occurrence_H, occurrence_L])
    significances.append(p_val)
#adjust the p-values for multiple testing using the Bonferroni correction
adjusted_pvals = multipletests(significances, method='bonferroni')[1]
mutations_significances = list(zip(mutations, significances, adjusted_pvals))
mutations_significances = [c for c in mutations_significances if c[1]<0.05]
mutations_significances = sorted(mutations_significances, key=lambda x:abs(x[0][1]-x[0][2]), reverse=True)
fnOut.write('NSCLC FGA0.2:\n')
for item in mutations_significances:
    content = item[0] + [item[1],item[2]]
    content = '\t'.join([str(c) for c in content])
    fnOut.write(content+'\n')

# Create a Counter object to count the occurrences of each element
counts = Counter(melanoma_FGA2H_genes)
total_count = len(melanoma_FGA2H)
melanoma_FGA2H_genes_prob = {element: count / total_count * 100 for element, count in counts.items()}

counts = Counter(melanoma_FGA2L_genes)
total_count = len(melanoma_FGA2L)
melanoma_FGA2L_genes_prob = {element: count / total_count * 100 for element, count in counts.items()}

# merge the two dictionaries
melanoma_genes_prob = {}
for key in melanoma_FGA2H_genes_prob.keys() | melanoma_FGA2L_genes_prob.keys():
    melanoma_genes_prob[key] = [0 if melanoma_FGA2H_genes_prob.get(key) is None else melanoma_FGA2H_genes_prob.get(key), 0 if melanoma_FGA2L_genes_prob.get(key) is None else melanoma_FGA2L_genes_prob.get(key)]
# chi-test for the significant difference of gene mutation presence
mutations = []
significances = []
for gene in melanoma_genes_prob:
    occurrence_H = melanoma_genes_prob[gene][0]
    occurrence_L = melanoma_genes_prob[gene][1]
    mutation_table = [[occurrence_H, 100 - occurrence_H], [occurrence_L, 100 - occurrence_L]]
    chi2_stat, p_val, dof, exp_freq = chi2_contingency(mutation_table)
    mutations.append([gene, occurrence_H, occurrence_L])
    significances.append(p_val)
# adjust the p-values for multiple testing using the Bonferroni correction
adjusted_pvals = multipletests(significances, method='bonferroni')[1]
mutations_significances = list(zip(mutations, significances, adjusted_pvals))
mutations_significances = [c for c in mutations_significances if c[1]<0.05]
mutations_significances = sorted(mutations_significances, key=lambda x:abs(x[0][1]-x[0][2]), reverse=True)
fnOut.write('melanoma FGA0.2:\n')
for item in mutations_significances:
    content = item[0] + [item[1], item[2]]
    content = '\t'.join([str(c) for c in content])
    fnOut.write(content+'\n')


################################################## AS0.2 ##################################################
melanoma_AS2H = CNAscore_info.loc[np.logical_and(CNAscore_info['CANCER_TYPE']=='Melanoma', CNAscore_info['AS0.2_score']=='high'), 'SAMPLE_ID'].tolist()
melanoma_AS2L = CNAscore_info.loc[np.logical_and(CNAscore_info['CANCER_TYPE']=='Melanoma', CNAscore_info['AS0.2_score']=='low'), 'SAMPLE_ID'].tolist()
NSCLC_AS2H = CNAscore_info.loc[np.logical_and(CNAscore_info['CANCER_TYPE']=='Non-Small Cell Lung Cancer', CNAscore_info['AS0.2_score']=='high'), 'SAMPLE_ID'].tolist()
NSCLC_AS2L = CNAscore_info.loc[np.logical_and(CNAscore_info['CANCER_TYPE']=='Non-Small Cell Lung Cancer', CNAscore_info['AS0.2_score']=='low'), 'SAMPLE_ID'].tolist()

melanoma_AS2H_genes = [list(sample_mutation_dict[sample].keys()) for sample in melanoma_AS2H if sample in sample_mutation_dict]
melanoma_AS2H_genes = sum(melanoma_AS2H_genes,[])
melanoma_AS2L_genes = [list(sample_mutation_dict[sample].keys()) for sample in melanoma_AS2L if sample in sample_mutation_dict]
melanoma_AS2L_genes = sum(melanoma_AS2L_genes,[])
NSCLC_AS2H_genes = [list(sample_mutation_dict[sample].keys()) for sample in NSCLC_AS2H if sample in sample_mutation_dict]
NSCLC_AS2H_genes = sum(NSCLC_AS2H_genes,[])
NSCLC_AS2L_genes = [list(sample_mutation_dict[sample].keys()) for sample in NSCLC_AS2L if sample in sample_mutation_dict]
NSCLC_AS2L_genes = sum(NSCLC_AS2L_genes,[])

# Create a Counter object to count the occurrences of each element
counts = Counter(NSCLC_AS2H_genes)
total_count = len(NSCLC_AS2H)
NSCLC_AS2H_genes_prob = {element: count / total_count * 100 for element, count in counts.items()}

counts = Counter(NSCLC_AS2L_genes)
total_count = len(NSCLC_AS2L)
NSCLC_AS2L_genes_prob = {element: count / total_count * 100 for element, count in counts.items()}

# merge the two dictionaries
NSCLC_genes_prob = {}
for key in NSCLC_AS2H_genes_prob.keys() | NSCLC_AS2L_genes_prob.keys():
    NSCLC_genes_prob[key] = [0 if NSCLC_AS2H_genes_prob.get(key) is None else NSCLC_AS2H_genes_prob.get(key), 0 if NSCLC_AS2L_genes_prob.get(key) is None else NSCLC_AS2L_genes_prob.get(key)]
# chi-test for the significant difference of gene mutation presence
mutations = []
significances = []
for gene in NSCLC_genes_prob:
    occurrence_H = NSCLC_genes_prob[gene][0]
    occurrence_L = NSCLC_genes_prob[gene][1]
    mutation_table = [[occurrence_H, 100 - occurrence_H], [occurrence_L, 100 - occurrence_L]]
    chi2_stat, p_val, dof, exp_freq = chi2_contingency(mutation_table)
    mutations.append([gene, occurrence_H, occurrence_L])
    significances.append(p_val)
# adjust the p-values for multiple testing using the Bonferroni correction
adjusted_pvals = multipletests(significances, method='bonferroni')[1]
mutations_significances = list(zip(mutations, significances, adjusted_pvals))
mutations_significances = [c for c in mutations_significances if c[1]<0.05]
mutations_significances = sorted(mutations_significances, key=lambda x:abs(x[0][1]-x[0][2]), reverse=True)
fnOut.write('NSCLC AS0.2:\n')
for item in mutations_significances:
    content = item[0] + [item[1], item[2]]
    content = '\t'.join([str(c) for c in content])
    fnOut.write(content+'\n')

# Create a Counter object to count the occurrences of each element
counts = Counter(melanoma_AS2H_genes)
total_count = len(melanoma_AS2H)
melanoma_AS2H_genes_prob = {element: count / total_count * 100 for element, count in counts.items()}

counts = Counter(melanoma_AS2L_genes)
total_count = len(melanoma_AS2L)
melanoma_AS2L_genes_prob = {element: count / total_count * 100 for element, count in counts.items()}

# merge the two dictionaries
melanoma_genes_prob = {}
for key in melanoma_AS2H_genes_prob.keys() | melanoma_AS2L_genes_prob.keys():
    melanoma_genes_prob[key] = [0 if melanoma_AS2H_genes_prob.get(key) is None else melanoma_AS2H_genes_prob.get(key), 0 if melanoma_AS2L_genes_prob.get(key) is None else melanoma_AS2L_genes_prob.get(key)]
# chi-test for the significant difference of gene mutation presence
mutations = []
significances = []
for gene in melanoma_genes_prob:
    occurrence_H = melanoma_genes_prob[gene][0]
    occurrence_L = melanoma_genes_prob[gene][1]
    mutation_table = [[occurrence_H, 100 - occurrence_H], [occurrence_L, 100 - occurrence_L]]
    chi2_stat, p_val, dof, exp_freq = chi2_contingency(mutation_table)
    mutations.append([gene, occurrence_H, occurrence_L])
    significances.append(p_val)
# adjust the p-values for multiple testing using the Bonferroni correction
adjusted_pvals = multipletests(significances, method='bonferroni')[1]
mutations_significances = list(zip(mutations, significances, adjusted_pvals))
mutations_significances = [c for c in mutations_significances if c[1]<0.05]
mutations_significances = sorted(mutations_significances, key=lambda x:abs(x[0][1]-x[0][2]), reverse=True)
fnOut.write('melanoma AS0.2:\n')
for item in mutations_significances:
    content = item[0] + [item[1], item[2]]
    content = '\t'.join([str(c) for c in content])
    fnOut.write(content+'\n')




################################################## AS0.1 ##################################################
melanoma_AS1H = CNAscore_info.loc[np.logical_and(CNAscore_info['CANCER_TYPE']=='Melanoma', CNAscore_info['AS0.1_score']=='high'), 'SAMPLE_ID'].tolist()
melanoma_AS1L = CNAscore_info.loc[np.logical_and(CNAscore_info['CANCER_TYPE']=='Melanoma', CNAscore_info['AS0.1_score']=='low'), 'SAMPLE_ID'].tolist()
NSCLC_AS1H = CNAscore_info.loc[np.logical_and(CNAscore_info['CANCER_TYPE']=='Non-Small Cell Lung Cancer', CNAscore_info['AS0.1_score']=='high'), 'SAMPLE_ID'].tolist()
NSCLC_AS1L = CNAscore_info.loc[np.logical_and(CNAscore_info['CANCER_TYPE']=='Non-Small Cell Lung Cancer', CNAscore_info['AS0.1_score']=='low'), 'SAMPLE_ID'].tolist()

melanoma_AS1H_genes = [list(sample_mutation_dict[sample].keys()) for sample in melanoma_AS1H if sample in sample_mutation_dict]
melanoma_AS1H_genes = sum(melanoma_AS1H_genes,[])
melanoma_AS1L_genes = [list(sample_mutation_dict[sample].keys()) for sample in melanoma_AS1L if sample in sample_mutation_dict]
melanoma_AS1L_genes = sum(melanoma_AS1L_genes,[])
NSCLC_AS1H_genes = [list(sample_mutation_dict[sample].keys()) for sample in NSCLC_AS1H if sample in sample_mutation_dict]
NSCLC_AS1H_genes = sum(NSCLC_AS1H_genes,[])
NSCLC_AS1L_genes = [list(sample_mutation_dict[sample].keys()) for sample in NSCLC_AS1L if sample in sample_mutation_dict]
NSCLC_AS1L_genes = sum(NSCLC_AS1L_genes,[])

# Create a Counter object to count the occurrences of each element
counts = Counter(NSCLC_AS1H_genes)
total_count = len(NSCLC_AS1H)
NSCLC_AS1H_genes_prob = {element: count / total_count * 100 for element, count in counts.items()}

counts = Counter(NSCLC_AS1L_genes)
total_count = len(NSCLC_AS1L)
NSCLC_AS1L_genes_prob = {element: count / total_count * 100 for element, count in counts.items()}

# merge the two dictionaries
NSCLC_genes_prob = {}
for key in NSCLC_AS1H_genes_prob.keys() | NSCLC_AS1L_genes_prob.keys():
    NSCLC_genes_prob[key] = [0 if NSCLC_AS1H_genes_prob.get(key) is None else NSCLC_AS1H_genes_prob.get(key), 0 if NSCLC_AS1L_genes_prob.get(key) is None else NSCLC_AS1L_genes_prob.get(key)]
# chi-test for the significant difference of gene mutation presence
mutations = []
significances = []
for gene in NSCLC_genes_prob:
    occurrence_H = NSCLC_genes_prob[gene][0]
    occurrence_L = NSCLC_genes_prob[gene][1]
    mutation_table = [[occurrence_H, 100 - occurrence_H], [occurrence_L, 100 - occurrence_L]]
    chi2_stat, p_val, dof, exp_freq = chi2_contingency(mutation_table)
    mutations.append([gene, occurrence_H, occurrence_L])
    significances.append(p_val)
# adjust the p-values for multiple testing using the Bonferroni correction
adjusted_pvals = multipletests(significances, method='bonferroni')[1]
mutations_significances = list(zip(mutations, significances, adjusted_pvals))
mutations_significances = [c for c in mutations_significances if c[1]<0.05]
mutations_significances = sorted(mutations_significances, key=lambda x:abs(x[0][1]-x[0][2]), reverse=True)
fnOut.write('NSCLC AS0.1:\n')
for item in mutations_significances:
    content = item[0] + [item[1], item[2]]
    content = '\t'.join([str(c) for c in content])
    fnOut.write(content+'\n')

# Create a Counter object to count the occurrences of each element
counts = Counter(melanoma_AS1H_genes)
total_count = len(melanoma_AS1H)
melanoma_AS1H_genes_prob = {element: count / total_count * 100 for element, count in counts.items()}

counts = Counter(melanoma_AS1L_genes)
total_count = len(melanoma_AS1L)
melanoma_AS1L_genes_prob = {element: count / total_count * 100 for element, count in counts.items()}

# merge the two dictionaries
melanoma_genes_prob = {}
for key in melanoma_AS1H_genes_prob.keys() | melanoma_AS1L_genes_prob.keys():
    melanoma_genes_prob[key] = [0 if melanoma_AS1H_genes_prob.get(key) is None else melanoma_AS1H_genes_prob.get(key), 0 if melanoma_AS1L_genes_prob.get(key) is None else melanoma_AS1L_genes_prob.get(key)]
# chi-test for the significant difference of gene mutation presence
mutations = []
significances = []
for gene in melanoma_genes_prob:
    occurrence_H = melanoma_genes_prob[gene][0]
    occurrence_L = melanoma_genes_prob[gene][1]
    mutation_table = [[occurrence_H, 100 - occurrence_H], [occurrence_L, 100 - occurrence_L]]
    chi2_stat, p_val, dof, exp_freq = chi2_contingency(mutation_table)
    mutations.append([gene, occurrence_H, occurrence_L])
    significances.append(p_val)
# adjust the p-values for multiple testing using the Bonferroni correction
adjusted_pvals = multipletests(significances, method='bonferroni')[1]
mutations_significances = list(zip(mutations, significances, adjusted_pvals))
mutations_significances = [c for c in mutations_significances if c[1]<0.05]
mutations_significances = sorted(mutations_significances, key=lambda x:abs(x[0][1]-x[0][2]), reverse=True)
fnOut.write('melanoma AS0.1:\n')
for item in mutations_significances:
    content = item[0] + [item[1], item[2]]
    content = '\t'.join([str(c) for c in content])
    fnOut.write(content+'\n')

fnOut.close()